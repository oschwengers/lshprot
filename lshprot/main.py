import argparse
import logging
import os
import pathlib
import platform as pf
import sys

from datetime import datetime
from pathlib import Path

import lshprot
import lshprot.db
import lshprot.search

def main():
    parser = argparse.ArgumentParser(
        prog=f'lshprot',
        description='Rapid detection of highly identical protein sequences',
        epilog=f'Version: {lshprot.__version__}',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=True
    )
    subparsers = parser.add_subparsers(dest='command', title='Subcommands')

    search_parser = subparsers.add_parser('search', description='Search sequences against a database', help='Search protein sequences against a database')
    arg_group_io = search_parser.add_argument_group('Input / Output')
    arg_group_io.add_argument('--query', '-q', action='store', type=lambda p: pathlib.Path(p).absolute(), default=None, help='Fasta file with query sequences')
    arg_group_io.add_argument('--db', '-d', action='store', type=lambda p: pathlib.Path(p).absolute(), default=None, help='Database path (default = CWD/lshprot.db)')
    arg_group_search = search_parser.add_argument_group('Search')
    arg_group_search.add_argument('--jaccard', '-j', action='store', type=float, default=0.90, help='Minimal sequence identity (default = 0.90)')
    arg_group_search.add_argument('--max-sequences', '-m', action='store', type=int, default=10, dest='max_sequences', help='Maximal number of returned sequences (default = 10)')
    arg_group_search.add_argument('--query-coverage', '-c', action='store', type=float, default=0.90, dest='query_coverage', help='Minimal query coverage (default = 0.90)')
    arg_group_search.add_argument('--subject-coverage', action='store', type=float, default=0.90, dest='subject_coverage', help='Minimal subject coverage (default = 0.90)')
    arg_group_search.add_argument('--identity', '-i', action='store', type=float, default=0.95, dest='identity', help='Minimal sequence identity (default = 0.95)')
    arg_group_general = search_parser.add_argument_group('General')
    arg_group_general.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    arg_group_general.add_argument('--debug', action='store_true', help='Run lshprot in debug mode. Temp data will not be removed.')
    arg_group_general.add_argument('--version', action='version', version=f'%(prog)s {lshprot.__version__}')
    
    db_parser = subparsers.add_parser('db', description='Build a database', help='Build a database')
    arg_group_io = db_parser.add_argument_group('Input / Output')
    arg_group_io.add_argument('--input', '-i', action='store', type=lambda p: pathlib.Path(p).absolute(), default=None, help='Fasta file with input sequences')
    arg_group_io.add_argument('--db', '-d', action='store', type=lambda p: pathlib.Path(p).absolute(), default=Path().joinpath('lshprot.db'), help='Database path (default = CWD/lshprot.db)')
    arg_group_io.add_argument('--tmp-dir', action='store', type=lambda p: pathlib.Path(p).absolute(), default=Path('.'), dest='tmp_dir', help='Temp directory (default = CWD)')
    # arg_group_io.add_argument('--keep', '-k', action='store_true', help='Keep existing db and add new sequences')
    arg_group_db = db_parser.add_argument_group('DB Parameters')
    arg_group_db.add_argument('--shingle', '-s', action='store', type=int, default=5, help='Shingle size in AA')
    arg_group_db.add_argument('--permutations', '-p', action='store', type=int, default=128, help='Number of hash permutations')
    arg_group_general = db_parser.add_argument_group('General')
    arg_group_general.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    arg_group_general.add_argument('--debug', action='store_true', help='Run lshprot in debug mode. Temp data will not be removed.')
    arg_group_general.add_argument('--version', action='version', version=f'%(prog)s {lshprot.__version__}')
    args = parser.parse_args()

    if args.command == None:
        parser.print_help()
    
    logging.basicConfig(
        format='%(asctime)s.%(msecs)03d - %(levelname)s - %(name)s - %(message)s',
        datefmt='%H:%M:%S',
        level=logging.DEBUG if args.debug else logging.INFO
    )
    log = logging.getLogger('MAIN')
    log.info('version=%s', lshprot.__version__)
    log.info('developer: Oliver Schwengers, github.com/oschwengers')
    log.info('command: %s', ' '.join(sys.argv))
    log.info('local time: %s', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    log.info('machine: type=%s, cores=%s', pf.processor(), os.cpu_count())
    log.info('system: type=%s, release=%s', pf.system(), pf.release())
    log.info('python: version=%s, implementation=%s', pf.python_version(), pf.python_implementation())

    if args.command == 'search':
        lshprot.search.main(args)
    elif args.command == 'db':
        lshprot.db.main(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()