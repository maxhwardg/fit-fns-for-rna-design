from itertools import starmap
import common
import random


def valid_random_pri(match):
    seq = []
    for i in range(len(match)):
        if match[i] < i:
            seq.append(random.choice(common.valid_pairs(seq[match[i]])))
        else:
            seq.append(random.choice(common.RNA_ALPHA))
    return ''.join(seq)


def valid_mutations(match, pri):
    muts = []
    for i in range(len(match)):
        if match[i] == i:
            for b in common.RNA_ALPHA:
                if b == pri[i]:
                    continue
                muts.append([(i, b)])
        elif match[i] < i:
            for b in common.RNA_ALPHA:
                for b2 in common.valid_pairs(b):
                    if b == pri[i] and b2 == pri[match[i]]:
                        continue
                    muts.append([(i, b), (match[i], b2)])
    return muts


def apply_mutation(pri, mut):
    mutated = list(pri)
    for i, b in mut:
        mutated[i] = b
    return ''.join(mutated)


def walk(db, metric, guide, steps, init_pri=None):
    match = common.db_to_matching(db)
    pri = valid_random_pri(match) if init_pri is None else init_pri
    scores = {}
    scores[pri] = metric(pri, db)
    for _ in range(steps):
        mutated = apply_mutation(pri, guide(
            pri, match, valid_mutations(match, pri)))
        if mutated not in scores:
            scores[mutated] = metric(mutated, db)
        if scores[mutated] > scores[pri]:
            pri = mutated
    return pri


def parallel_walk(dbs, metric, guide, steps, init_pris=None, processes=None):
    matches = [common.db_to_matching(db) for db in dbs]
    pris = init_pris
    if pris is None:
        pris = [valid_random_pri(match) for match in matches]
    metric_cache = {}
    from multiprocessing import Pool
    with Pool(processes=processes) as pool:
        def add_to_metric_cache(pris):
            metric_args = []
            for i in range(len(dbs)):
                if (pris[i], dbs[i]) not in metric_cache:
                    metric_args.append((pris[i], dbs[i]))
            metrics = pool.starmap(metric, metric_args)
            for i, (pri, db) in enumerate(metric_args):
                metric_cache[(pri, db)] = metrics[i]

        add_to_metric_cache(pris)

        for _ in range(steps):
            muts = []
            for i in range(len(dbs)):
                muts.append(
                    (pris[i], matches[i], valid_mutations(matches[i], pris[i])))
            mutated = list(starmap(apply_mutation, zip(pris, guide(muts))))
            add_to_metric_cache(mutated)
            for i in range(len(dbs)):
                if metric_cache[(mutated[i], dbs[i])] > metric_cache[(pris[i], dbs[i])]:
                    pris[i] = mutated[i]
    return pris

def random_guide(pri, match, muts):
    return random.choice(muts)


def random_batch_guide(muts):
    return list(starmap(random_guide, muts))

