import csv
import argparse

import numpy

import vienna
import RNA
import common
import multiprocessing
import fit_fns


class RnaRecord:
    def __init__(self, id, pri, db):
        self.id = id
        self.pri = pri
        self.db = db


def archiveii_tsv_generator(file_path):
    file = open(file_path, "r")
    for row in csv.reader(file, delimiter="\t"):
        name, primary, db = row[0], row[1], row[2]
        yield RnaRecord(name, primary, db)
    file.close()


def get_top_k_folds(primary, k):
    ctx = vienna.ViennaContext(primary)
    subs = []
    delta = 0.0
    mfes = []
    while len(subs) < k and delta < 10:
        subs = ctx.subopt(delta)
        if delta == 0.0:
            mfes = subs
        delta += 0.2
    return subs, mfes


def remove_invalid_pairs(primary, db):
    new_db = ['.'] * len(db)
    match = common.db_to_matching(db)
    for i in range(len(match)):
        if abs(match[i]-i) < 4:
            continue
        if match[i] != i and not common.valid_pair(primary[i], primary[match[i]]):
            continue
        if match[i] < i:
            new_db[i] = ')'
        elif match[i] > i:
            new_db[i] = '('
    return ''.join(new_db)


def relative_metric_rank(primary, db, metric, subs, mfes, distance_threshold):
    ctx = vienna.ViennaContext(primary)
    sorted = []
    for s in subs:
        sorted.append((metric(ctx, s, mfes), s))
    sorted.sort()
    for i in range(len(sorted)):
        if sorted[i][1] == db:
            return i / len(sorted)

    min_d_rank = None
    min_d = 1e10
    for i in range(len(sorted)):
        dist = RNA.bp_distance(sorted[i][1], db)
        if dist < min_d:
            min_d = dist
            min_d_rank = i
    return None if min_d/len(db) > distance_threshold else min_d_rank/len(sorted)


def inv_probability(ctx, db, mfes):
    return 1.0 - ctx.prob(db)


def ensemble_defect(ctx, db, mfes):
    return ctx.ensemble_defect(db)


def mfe_distance(ctx, db, mfes):
    return RNA.bp_distance(ctx.mfe(), db)


def avg_mfe_bpd(ctx, db, mfes):
    sum = 0
    for s in mfes:
        sum += RNA.bp_distance(s, db)
    return sum / len(mfes)


def avg_mfe_hd(ctx, db, mfes):
    sum = 0
    for s in mfes:
        sum += fit_fns.hamming_dist(s, db)
    return sum / len(mfes)


def process(packed_args):
    args, record, metrics = packed_args
    db = remove_invalid_pairs(record.pri, record.db)
    top_k, mfes = get_top_k_folds(record.pri, args.k)
    res = []
    for metric in metrics:
        res.append(relative_metric_rank(
            record.pri, db, metric, top_k, mfes, args.distance))
    return res


def main():
    ap = argparse.ArgumentParser(
        description="Benchmark metrics on real RNA")
    ap.add_argument("-t", "--tsv", type=str, default="data/archiveII.tsv",
                          help="TSV input file. Must contain tsv rows containing name, primary, db")
    ap.add_argument("-k", "--k", type=int, default=200000,
                    help="Number of suboptimal to use.")
    ap.add_argument("-d", "--distance", type=float, default=0.05,
                    help='Base pair distance threshold to use when a perfect match is not found.')
    ap.add_argument("--processes", type=int, default=1,
                    help="Number of processes to use")
    args = ap.parse_args()
    metrics = [inv_probability, ensemble_defect, avg_mfe_hd]
    with multiprocessing.Pool(args.processes) as pool:
        records = list(archiveii_tsv_generator(args.tsv))
        ranks = pool.map(process, [(args, record, metrics)
                         for record in records])
    res = [[] for _ in range(len(ranks[0]))]
    for i in range(len(ranks)):
        print(records[i].id, ranks[i])
        if None in ranks[i]:
            continue
        for j in range(len(ranks[i])):
            res[j].append(ranks[i][j])

    for i in range(len(res)):
        r = res[i]
        print(f'Result for {metrics[i].__name__}')
        print('average:', numpy.average(r))
        print('median: ', sorted(r)[len(r)//2])
        print('std: ', numpy.std(r))


if __name__ == '__main__':
    main()
