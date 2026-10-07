"""Explicit local setup diagnostics, loopback server and bounded queue execution."""
import argparse
import importlib
from importlib.metadata import version, PackageNotFoundError
import json
import os
from pathlib import Path
import shutil
import socket
from uuid import uuid4

DEPENDENCIES = {'biopython':'Bio', 'primer3-py':'primer3', 'fastapi':'fastapi',
                'uvicorn':'uvicorn', 'python-multipart':'multipart'}


def inspect_environment(directory, binary=None):
    """Read-only checks: never create directories, migrate databases or claim jobs."""
    checks=[]
    target=Path(directory).expanduser().resolve()
    ancestor=target
    while not ancestor.exists():ancestor=ancestor.parent
    writable=ancestor.is_dir() and os.access(ancestor,os.W_OK|os.X_OK)
    checks.append({'name':'data_directory','ready':writable,
                   'detail':'writable parent found' if writable else 'data path or parent is not a writable directory'})
    for package,module in DEPENDENCIES.items():
        try:
            importlib.import_module(module);installed=version(package)
            checks.append({'name':package,'ready':True,'version':installed})
        except (ImportError,OSError,PackageNotFoundError):
            checks.append({'name':package,'ready':False,'detail':'dependency missing or unable to load; reinstall the package in this environment'})
    tool={'ready':False,'detail':'Set --minimap2 or MOLBIO_MINIMAP2 to validated minimap2 2.31-r1302; optional for Sanger and library work.'}
    if binary:
        try:
            from .alignment_evidence import tool_identity
            _,identity=tool_identity(Path(binary).expanduser())
            tool={'ready':True,**identity}
        except (OSError,ValueError):
            tool={'ready':False,'detail':'Configured minimap2 is unavailable or not validated version 2.31-r1302.'}
    return {'ready':all(c['ready'] for c in checks),'data_directory':str(target),
            'database_exists':(target/'workbench.sqlite3').is_file(),
            'free_bytes':shutil.disk_usage(ancestor).free,'checks':checks,'alignment':tool,
            'execution':'foreground only; server does not start workers'}


def positive(value):
    parsed=int(value)
    if not 1 <= parsed <= 100:raise argparse.ArgumentTypeError('Choose 1–100 tasks')
    return parsed


def port_number(value):
    parsed=int(value)
    if not 1 <= parsed <= 65535:raise argparse.ArgumentTypeError('Choose a port from 1–65535')
    return parsed


def emit(value):
    print(json.dumps(value,ensure_ascii=False),flush=True)


def main(argv=None):
    parser=argparse.ArgumentParser(description='Local MolBio: explicit data storage, local serving and bounded queue execution.')
    parser.add_argument('--data-dir',default=os.environ.get('MOLBIO_DATA_DIR'),help='Data directory (or MOLBIO_DATA_DIR); required, never silently chosen from current directory')
    parser.add_argument('--minimap2',default=os.environ.get('MOLBIO_MINIMAP2'),help='Path to validated minimap2 2.31-r1302; not downloaded automatically')
    sub=parser.add_subparsers(dest='command',required=True)
    doctor=sub.add_parser('doctor',help='Check dependencies, storage and optional aligner without changing data')
    doctor.add_argument('--require-alignment',action='store_true',help='Return failure unless the pinned aligner is usable')
    serve=sub.add_parser('serve',help='Run the UI/API on loopback in the foreground; no workers start automatically')
    serve.add_argument('--port',type=port_number,default=8000)
    work=sub.add_parser('work',help='Execute up to a finite number of queued tasks, stopping at idle or first failure')
    work.add_argument('--adapter',choices=('ab1','fastq','alignment'),required=True,
                      help='ab1/fastq only check registered input identity; alignment runs local FASTQ alignment')
    work.add_argument('--max-jobs',type=positive,default=1)
    args=parser.parse_args(argv)
    if not args.data_dir or not args.data_dir.strip():parser.error('Specify --data-dir or MOLBIO_DATA_DIR; your storage location must be explicit')
    try:check=inspect_environment(args.data_dir,args.minimap2)
    except OSError:
        emit({'status':'rejected','error':'storage_unavailable','next':'Check the selected data path and its parent permissions.'});return 2
    needs_alignment=(args.command=='doctor' and args.require_alignment) or (args.command=='work' and args.adapter=='alignment')
    if args.command=='doctor':
        emit(check);return 0 if check['ready'] and (not needs_alignment or check['alignment']['ready']) else 2
    if not check['ready'] or (needs_alignment and not check['alignment']['ready']):
        emit({'status':'rejected','error':'environment_not_ready','diagnostics':check});return 2
    os.environ['MOLBIO_DATA_DIR']=check['data_directory']
    if args.minimap2:os.environ['MOLBIO_MINIMAP2']=str(Path(args.minimap2).expanduser().resolve())
    if args.command=='work':
        if not check['database_exists']:
            emit({'status':'rejected','error':'database_not_initialized','next':'Start serve, import data and queue a task first.'});return 2
        from . import worker
        from .database import initialise
        adapters={'ab1':worker.ADAPTER,'fastq':worker.fastq_validation.ADAPTER,'alignment':worker.fastq_alignment.ADAPTER}
        initialise();worker_id='foreground-'+str(uuid4())
        try:
            for _ in range(args.max_jobs):
                result=worker.run_once(worker_id,adapter=adapters[args.adapter]);emit(result)
                if result['status']=='idle':return 0
                if result['status']!='succeeded':return 1
        except KeyboardInterrupt:
            emit({'status':'interrupted','next':'Any active lease must expire before another worker can reclaim it.'});return 130
        return 0
    # Reserve the port before startup/migrations; an occupied port must not modify data.
    with socket.socket(socket.AF_INET,socket.SOCK_STREAM) as sock:
        try:sock.bind(('127.0.0.1',args.port));sock.listen(128)
        except OSError:
            emit({'status':'rejected','error':'port_unavailable','next':'Choose another --port or stop the existing service.'});return 2
        import uvicorn
        emit({'status':'starting','url':f'http://127.0.0.1:{args.port}',
              'data_directory':check['data_directory'],'alignment_ready':check['alignment']['ready'],
              'worker':'not_started; use the explicit work command for queued tasks'})
        server=uvicorn.Server(uvicorn.Config('localmolbio.api:app',host='127.0.0.1',port=args.port))
        server.run(sockets=[sock])
        return 0 if server.started else 1


if __name__=='__main__':
    raise SystemExit(main())
