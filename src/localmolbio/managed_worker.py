"""Opt-in foreground service companion. No HTTP endpoint can start it."""
import multiprocessing
import os
from pathlib import Path
import signal
import threading
import time
from uuid import uuid4

STATES = {0:'starting',1:'idle',2:'running',3:'stopping',4:'stopped',5:'failed'}


def _run(stop, ready, state, completed, owner_pid):
    # Isolate the worker and its aligner child from terminal signals. The owner stops both.
    os.setsid()
    from .worker import run_once
    from .fastq_alignment import ADAPTER
    from .input_validation import InputValidationError
    worker_id='managed-'+str(uuid4())
    def parent_alive():return os.getppid()==owner_pid
    def pulse():
        if not parent_alive():raise InputValidationError('worker_parent_stopped',retryable=True)
    def claimed():state.value=2
    state.value=1;ready.set()
    try:
        while not stop.is_set() and parent_alive():
            result=run_once(worker_id,adapter=ADAPTER,on_claim=claimed,on_pulse=pulse)
            if result['status'] not in ('idle','succeeded'):
                state.value=5;return  # No unattended retry loop after the first failed attempt.
            if result['status']=='succeeded':completed.value+=1
            state.value=1
            if result['status']=='idle':stop.wait(.5)
        state.value=4
    except BaseException:
        state.value=5
        raise


class ManagedAlignmentWorker:
    def __init__(self,directory,shutdown_timeout=70):
        self.directory=Path(directory)
        self.shutdown_timeout=shutdown_timeout
        self.process=None;self.lock=None;self.monitor=None
        self.reap_lock=threading.Lock();self.group_reaped=False
        ctx=multiprocessing.get_context('spawn')
        self.event=ctx.Event();self.ready=ctx.Event()
        self.state=ctx.Value('i',0);self.completed=ctx.Value('i',0)
        self.context=ctx

    def acquire(self):
        try:import fcntl
        except ImportError:raise RuntimeError('managed_worker_requires_posix') from None
        self.directory.mkdir(parents=True,exist_ok=True)
        self.lock=(self.directory/'.alignment-worker.lock').open('a')
        try:fcntl.flock(self.lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except OSError:
            self.lock.close();self.lock=None
            raise RuntimeError('managed_alignment_worker_already_active') from None

    def start(self):
        if self.lock is None:raise RuntimeError('managed_worker_lock_required')
        if self.process is not None:raise RuntimeError('managed_worker_already_started')
        self.process=self.context.Process(target=_run,args=(self.event,self.ready,self.state,self.completed,os.getpid()),name='molbio-alignment')
        self.process.start()
        if not self.ready.wait(10):
            self.stop();raise RuntimeError('managed_worker_start_failed')
        self.monitor=threading.Thread(target=self._watch_exit,daemon=True)
        self.monitor.start()

    def _reap_group(self):
        # Exactly once, promptly after child death: do not leave its aligner running.
        with self.reap_lock:
            if self.group_reaped or not self.ready.is_set():return
            self.group_reaped=True
            try:os.killpg(self.process.pid,signal.SIGKILL)
            except ProcessLookupError:pass

    def _watch_exit(self):
        while self.process.is_alive():time.sleep(.1)
        self._reap_group()

    def status(self):
        state=STATES[self.state.value]
        if self.process is not None and not self.process.is_alive() and state not in ('stopped','failed'):
            state='failed'
        if self.event.is_set() and self.process is not None and self.process.is_alive():state='stopping'
        return {'enabled':True,'state':state,'completed_tasks':self.completed.value,
                'adapter':'alignment','failure_policy':'stop_until_service_restart'}

    def stop(self):
        self.event.set()
        if self.process is None or self.process.pid is None:return
        self.process.join(self.shutdown_timeout)
        if self.process.is_alive():
            # Child creates its own session before signalling ready; kill its entire group.
            if self.ready.is_set():self._reap_group()
            else:self.process.kill()
            self.process.join(5)
            self.state.value=5
        self._reap_group()
        if self.monitor is not None:self.monitor.join(2)

    def release(self):
        self.stop()
        if self.lock is not None:
            self.lock.close();self.lock=None
