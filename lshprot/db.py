# import argparse
import logging
import pickle
import sqlite3
import sys
import os

from pathlib import Path

from alive_progress import alive_bar
from Bio import SeqIO
from xopen import xopen

import lshprot.hashing as lshhash
import lshprot.util as util


log = logging.getLogger('DB')


def main(args):
    # parameter checks
    db_path = Path(args.db).resolve()
    log.info('parameter: db=%s', db_path)
    tmp_path = Path(args.tmp_dir).resolve()
    log.info('parameter: tmp-dir=%s', tmp_path)
    # ToDo check db path
    if args.shingle < 2:
        print(f"Error: shingle size must be >= 2!")
        sys.exit(-1)
    log.info('parameter: shingle=%i', args.shingle)
    if args.permutations <= 0:
        print(f"Error: number of permutations must be >= 0!")
        sys.exit(-1)
    log.info('parameter: permutations=%i', args.permutations)

    os.environ['SQLITE_TMPDIR'] = str(tmp_path)  # set SQLITE tmp dir env var

    # init and fill LSH
    with sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn, xopen(str(args.input), threads=0) as fh:
        conn.execute('PRAGMA page_size = 4096;')
        conn.execute('PRAGMA cache_size = 100000;')
        conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
        conn.execute(f'PRAGMA mmap_size = {20 * 1024 * 1024 * 1024};')
        conn.execute('PRAGMA synchronous = OFF;')
        conn.execute('PRAGMA journal_mode = OFF')
        conn.execute('PRAGMA threads = 2;')

        # if not args.keep:
        log.info('drop existing db')
        conn.execute('DROP TABLE IF EXISTS signatures;')
        log.debug('discarded table signatures')
        log.info('create table signatures')
        stmt = '''CREATE TABLE signatures (
            id TEXT PRIMARY KEY,
            sequence TEXT NOT NULL,
            hashes BLOB NOT NULL
            ) WITHOUT ROWID;'''
        stmt = ' '.join(stmt.replace('\n', '').split())
        conn.execute(stmt)
        log.debug('created table signatures')
        conn.execute('DROP TABLE IF EXISTS parameters;')
        log.debug('discarded table parameters')
        log.info('create table parameters')
        stmt = '''CREATE TABLE parameters (
            parameter TEXT PRIMARY KEY,
            value NUMBER NOT NULL
            ) WITHOUT ROWID;'''
        stmt = ' '.join(stmt.replace('\n', '').split())
        conn.execute(stmt)
        log.debug('created table signatures')
        
        log.debug('write parameters')
        conn.execute('INSERT INTO parameters (parameter, value) VALUES (?,?)', ('shingle', args.shingle))
        log.debug(f'INSERT INTO parameters (parameter, value) VALUES (shingle,{args.shingle})')
        conn.execute('INSERT INTO parameters (parameter, value) VALUES (?,?)', ('permutations', args.permutations))
        log.debug(f'INSERT INTO parameters (parameter, value) VALUES (permutations,{args.permutations})')
        log.info('stored db parameters: shingle=%i, permutations=%i', args.shingle, args.permutations)

        i = 0
        log.info('build and store sequence hash signatures')
        with alive_bar(scale='SI') as bar:
            inserted_ids = set()
            for record in SeqIO.parse(fh, 'fasta'):
                id = record.id
                if(id in inserted_ids):
                    raise ValueError(f'Duplicated ID in Fasta file! id={id}. Please, make sure that there are no cuplicates!')
                else:
                    inserted_ids.add(id)
                seq = str(record.seq).upper()
                if(seq[-1] == '*'):  # remove trailing stop asterik
                    seq = seq[:-1]
                if(util.FASTA_AA_SEQUENCE_PATTERN.fullmatch(seq) is None):
                    raise ValueError(f'Fasta sequence contains invalid AA characters! id={id}')

                minhash = lshhash.create_minhash(seq, shingle=args.shingle, permutations=args.permutations)
                # print(f'size of minhash: {sys.getsizeof(minhash)}')
                minhash_pickled = pickle.dumps(minhash, protocol=5)
                # print(f'size of pickled minhash: {sys.getsizeof(minhash_pickled)}')
                conn.execute('INSERT INTO signatures (id, sequence, hashes) VALUES (?,?,?)', (id, seq, minhash_pickled))
                log.debug('added sequence: id=%s, length=%i, sequence=%s..%s', id, len(seq), seq[:10], seq[-10:])
                i += 1
                bar(count=1)
                if i%10_000 == 0:
                    conn.commit()
                    # log.info('processed %i seqs', i)
        log.info('written sequences: %i', i)
        conn.commit()
        conn.execute('VACUUM')
        conn.commit()
    print(f'successfully wrote {i} sequences into db')