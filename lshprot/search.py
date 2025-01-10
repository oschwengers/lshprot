import logging
import pickle
import sqlite3
import sys

from pathlib import Path

from Bio import SeqIO
from datasketch import MinHashLSH
import pyopal
from xopen import xopen

import lshprot.hashing as lshhash
import lshprot.util as util


log = logging.getLogger('SEARCH')


def main(args):
    # parameter checks
    if args.db is None:
        print(f"Error: wrong database path! Type 'lshprot search --help'")
        sys.exit(-1)
    db_path = Path(args.db).resolve()
    # ToDo check db path
    if args.jaccard <= 0 or args.jaccard > 1:
        print(f"Error: Jaccard threshold must be between 0 and 1!")
        sys.exit(-1)
    if args.query_coverage <= 0 or args.query_coverage > 1:
        print(f"Error: query coverage must be between 0 and 1!")
        sys.exit(-1)
    if args.subject_coverage <= 0 or args.subject_coverage > 1:
        print(f"Error: subject coverage must be between 0 and 1!")
        sys.exit(-1)
    if args.identity <= 0 or args.identity > 1:
        print(f"Error: identity must be between 0 and 1!")
        sys.exit(-1)


    # fetch parameters from DB
    log.debug('fetch parameters from DB')
    with sqlite3.connect(f"file:{db_path}?mode=ro&nolock=1&cache=shared", uri=True, check_same_thread=False) as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT value FROM parameters WHERE parameter='shingle';")
        res = cursor.fetchone()
        args.shingle = res[0]
        log.info('parameters: shingle=%i', args.shingle)
        
        cursor = conn.cursor()
        cursor.execute("SELECT value FROM parameters WHERE parameter='permutations';")
        res = cursor.fetchone()
        args.permutations = res[0]
        log.info('parameters: permutations=%i', args.permutations)


    # init and fill LSH index
    log.debug('init MinHashLSH')
    lsh_index = MinHashLSH(threshold=args.jaccard, num_perm=args.permutations, weights=(0.3,0.7))
    with sqlite3.connect(f"file:{db_path}?mode=ro&nolock=1&cache=shared", uri=True, check_same_thread=False) as conn:
        conn.execute(f'PRAGMA mmap_size = {128 * 1024 * 1024};')  # 128 Mb
        conn.execute('PRAGMA omit_readlock;')
        cursor = conn.cursor()
        cursor.execute("SELECT id, hashes FROM signatures;")  # Initialize LSH index with existing MinHash signatures
        i = 0
        for id, minhash_bytes in cursor.fetchall():
            minhash = pickle.loads(minhash_bytes)  # Deserialize MinHash object
            lsh_index.insert(id, minhash, check_duplication=False)  # Populate LSH index
            i += 1
        log.info('read DB sequences into index: seqs=%i', i)


    # search query sequences
    with xopen(str(args.query), threads=0) as fh, sqlite3.connect(f"file:{db_path}?mode=ro&nolock=1&cache=shared", uri=True, check_same_thread=False) as conn:
        for record in SeqIO.parse(fh, 'fasta'):
            query_id = record.id
            query_seq = str(record.seq).upper()
            if(query_seq[-1] == '*'):  # remove trailing stop asterik
                query_seq = query_seq[:-1]
            if(util.FASTA_AA_SEQUENCE_PATTERN.fullmatch(query_seq) is None):
                raise ValueError(f'Fasta sequence contains invalid AA characters! id={record.id}')
            log.debug('compute query seq: id=%s, seq=%s..%s', query_id, query_seq[:10], query_seq[-10:])

            log.debug('compute query minhash')
            query_minhash = lshhash.create_minhash(query_seq, args.shingle, args.permutations)
            hits = []
            log.debug('query LSH index')
            hit_ids = lsh_index.query(query_minhash)  # Query the LSH index with the user's MinHash to find similar sequences
            cursor = conn.cursor()
            log.info('LSH hits: %i', len(hit_ids))
            for id in hit_ids:
                log.debug('retrieve hit sequence and minhashes')
                cursor.execute("SELECT sequence, hashes FROM signatures WHERE id = ?;", (id,))
                hit_seq, minhash_bytes = cursor.fetchone()
                minhash = pickle.loads(minhash_bytes)  # Deserialize MinHash object
                hit = {  # store hit
                    'id': id,
                    'sequence': hit_seq
                }

                # Calculate the Jaccard similarity percentage
                log.debug('compute jaccard')
                jaccard = query_minhash.jaccard(minhash)
                hit['jaccard'] = jaccard
                log.info('jaccard index: id=%s, jaccard=%.2f', id, jaccard)
                if jaccard >= args.jaccard:
                    log.debug('compute alignment')
                    for res in pyopal.align(query_seq, algorithm='nw', mode='full', threads=1, database=[hit_seq]):
                        # if res.identity() >= args.identity  and  res.coverage('query') >= args.query_coverage  and  res.coverage('target') >= args.subject_coverage:
                        hit['identity'] = res.identity()
                        hit['query-cov'] = res.coverage('query')
                        hit['subject-cov'] = res.coverage('target')
                        hit['score'] = res.score
                        hit['alignment'] = res.alignment
                        hit['cigar'] = res.cigar()
                        log.info(
                            'alignment: id=%s, idenity=%f, query-cov=%f, subject-cov=%f, score=%f, cigar=%s',
                            id, hit['identity'], hit['query-cov'], hit['subject-cov'], hit['score'], hit['cigar']
                        )
                        if hit['identity'] < args.identity:
                            log.debug('discard hit: id=%s, identity=%f', id, hit['identity'])
                        elif hit['query-cov'] < args.query_coverage:
                            log.debug('discard hit: id=%s, query-cov=%f', id, hit['query-cov'])
                        elif hit['subject-cov'] < args.subject_coverage:
                            log.debug('discard hit: id=%s, subject-cov=%f', id, hit['subject-cov'])
                        else:
                            hits.append(hit)
                else:
                    log.debug('discard hit: id=%s, jaccard=%f', id, jaccard)
            cursor.close()

            for hit in sorted(hits, key=lambda k: k['score'], reverse=True )[:args.max_sequences]:
                print(f"{query_id}\t{hit['id']}\t{hit['jaccard']}\t{hit['identity']}\t{hit['query-cov']}\t{hit['subject-cov']}\t{hit['score']}\t{hit['cigar']}\t{hit['sequence']}")
