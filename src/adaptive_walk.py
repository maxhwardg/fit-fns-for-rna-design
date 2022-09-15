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


def walk(db, metric, guide, steps, seed=0, init_pri=None):
    random.seed(seed)
    match = common.db_to_matching(db)
    pri = valid_random_pri(match) if init_pri is None else init_pri
    scores = {}
    scores[pri] = metric(pri, db)
    for _ in range(steps):
        mutated = apply_mutation(pri, guide(
            pri, match, valid_mutations(match, pri)))
        if mutated not in scores:
            scores[mutated] = metric(mutated, db)
        if scores[mutated] >= scores[pri]:
            pri = mutated
    return pri


def random_guide(pri, match, muts):
    return random.choice(muts)
