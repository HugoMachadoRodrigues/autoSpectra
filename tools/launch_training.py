#!/usr/bin/env python3
"""Launch the autoSpectra VAE training as a supervised daemon.

Two failure modes motivated this:

1. The run died whenever the parent session was torn down: `nohup ... &` leaves
   the child in the session's process group, so a group-wide kill reaches it.
   Fixed by the full daemonisation dance (fork -> setsid -> fork), which leaves
   R orphaned to init with no controlling terminal.

2. TensorFlow aborts hard between models on macOS:
     libc++abi: terminating due to uncaught exception of type
     std::__1::system_error: condition_variable wait failed: Invalid argument
   train_ossl.R now calls clear_session() between properties, but a C++ abort
   cannot be caught from R. So we supervise: on a non-zero exit, relaunch.
   train_ossl.R skips properties that already have both a model and metrics on
   disk, so each restart resumes instead of redoing finished work.

Usage:
    python3 tools/launch_training.py <logfile> [max_attempts] -- <train_ossl.R args>

Example:
    python3 tools/launch_training.py /tmp/run.log 12 -- \
        --only-vae --props oc,clay.tot,ph.h2o,n.tot

Prints the supervisor PID on stdout.
"""
import os
import subprocess
import sys
import time

# Repository root: this file lives in tools/, so the parent of its directory.
WORKTREE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ENV = {
    # tensorflow-metal 1.2.0 + our masked_mse / Sampling layer produce
    # NaN-collapsed predictions on Apple Silicon -> force CPU.
    "DISABLE_TENSORFLOW_METAL": "1",
    "CUDA_VISIBLE_DEVICES": "-1",
    "TF_CPP_MIN_LOG_LEVEL": "2",
    # Keep BLAS from oversubscribing cores against Keras' own thread pool;
    # thread churn is implicated in the condition_variable abort above.
    "OMP_NUM_THREADS": "4",
    "VECLIB_MAXIMUM_THREADS": "4",
    # Constrain TensorFlow's own pools. The abort is a race inside TF's
    # threading layer ("condition_variable wait failed"), it recurs on the FIRST
    # model of a fresh process, and it is not memory (24.8 GB free of 64, no
    # swap) -- so bound the pools rather than hope. Fewer, longer-lived threads
    # trade a little throughput for not aborting.
    "TF_NUM_INTRAOP_THREADS": "4",
    "TF_NUM_INTEROP_THREADS": "2",
    "TF_ENABLE_ONEDNN_OPTS": "0",
}


def log(msg):
    sys.stdout.write("\n=== supervisor: %s ===\n" % msg)
    sys.stdout.flush()


def main():
    argv = sys.argv[1:]
    train_args = []
    if "--" in argv:
        i = argv.index("--")
        train_args = argv[i + 1:]
        argv = argv[:i]
    if not argv:
        sys.exit(__doc__)

    log_path = argv[0]
    max_attempts = int(argv[1]) if len(argv) > 1 else 12
    pid_path = log_path + ".pid"
    cmd = ["Rscript", "train_ossl.R"] + train_args

    # --- first fork: parent reports the pid and returns ---
    if os.fork() > 0:
        for _ in range(100):
            time.sleep(0.05)
            if os.path.exists(pid_path):
                with open(pid_path) as fh:
                    print(fh.read().strip())
                return
        print("ERROR: daemon did not report a pid", file=sys.stderr)
        sys.exit(1)

    # --- detach from the controlling terminal and the session process group ---
    os.setsid()

    # --- second fork: the daemon can never reacquire a terminal ---
    if os.fork() > 0:
        os._exit(0)

    os.chdir(WORKTREE)
    os.umask(0)

    log_fd = os.open(log_path, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
    devnull = os.open(os.devnull, os.O_RDONLY)
    os.dup2(devnull, 0)
    os.dup2(log_fd, 1)
    os.dup2(log_fd, 2)
    os.close(devnull)
    if log_fd > 2:
        os.close(log_fd)

    with open(pid_path, "w") as fh:
        fh.write(str(os.getpid()) + "\n")

    env = dict(os.environ)
    env.update(ENV)

    # --- supervise ---
    for attempt in range(1, max_attempts + 1):
        log("attempt %d/%d starting" % (attempt, max_attempts))
        rc = subprocess.call(cmd, env=env)
        if rc == 0:
            log("training finished cleanly on attempt %d" % attempt)
            break
        log("exited with code %d - relaunching to resume" % rc)
        time.sleep(5)
    else:
        log("giving up after %d attempts" % max_attempts)

    try:
        os.unlink(pid_path)
    except OSError:
        pass


if __name__ == "__main__":
    main()
