import common
import vienna
import RNA


def prob_diff(ctx, ss):
    top_2 = ctx.top_2()
    pr = ctx.prob(ss)
    if ss == top_2[0]:
        return 1.0 if top_2[1] == None else pr-ctx.prob(top_2[1])
    else:
        return pr-ctx.prob(top_2[0])


def normalised_prob_diff(pri, db):
    ctx = vienna.ViennaContext(pri)
    top2 = ctx.top_2()
    dbp = ctx.prob(db)
    v = 0
    if db != top2[0]:
        p0 = ctx.prob(top2[0])
        v = (1-(p0-dbp)/p0)/2
    else:
        if top2[1] is None:
            v = 1
        else:
            v = (dbp-ctx.prob(top2[1])+1)/2
    return v


def normalised_ensemble_antidefect(pri, db):
    ctx = vienna.ViennaContext(pri)
    return 1.0-ctx.ensemble_defect(db)


def probability(pri, db):
    ctx = vienna.ViennaContext(pri)
    return ctx.prob(db)


def free_energy(pri, db):
    ctx = vienna.ViennaContext(pri)
    return -ctx.free_energy(db)


def hamming_dist(db1, db2):
    cnt = 0
    m1 = common.db_to_matching(db1)
    m2 = common.db_to_matching(db2)
    for i in range(len(db1)):
        if m1[i] != m2[i]:
            cnt += 1
    return cnt


def mfe_dist(pri, db, dist):
    ctx = vienna.ViennaContext(pri)
    return len(db)-dist(ctx.mfe(), db)


def avg_mfe_dist(pri, db, dist):
    ctx = vienna.ViennaContext(pri)
    sum = 0
    subs = ctx.subopt(0)
    for s in subs:
        sum += dist(s, db)
    return len(db)-(sum / len(subs))


def min_mfe_dist(pri, db, dist):
    ctx = vienna.ViennaContext(pri)
    mn = 1e18
    subs = ctx.subopt(0)
    for s in subs:
        mn = min(mn, dist(s, db))
    return len(db)-mn


def mfe_bpd(pri, db):
    return mfe_dist(pri, db, RNA.bp_distance)


def avg_mfe_bpd(pri, db):
    return avg_mfe_dist(pri, db, RNA.bp_distance)


def min_mfe_bpd(pri, db):
    return min_mfe_dist(pri, db, RNA.bp_distance)


def mfe_hd(pri, db):
    return mfe_dist(pri, db, hamming_dist)


def avg_mfe_hd(pri, db):
    return avg_mfe_dist(pri, db, hamming_dist)


def min_mfe_hd(pri, db):
    return min_mfe_dist(pri, db, hamming_dist)
