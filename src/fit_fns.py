import common
import vienna
import RNA
import math


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


def inv_interaction_network_fidelity(pred, true):
    m1 = common.db_to_matching(pred)
    m2 = common.db_to_matching(true)

    def pairs(m):
        pairs = set()
        for i in range(len(m)):
            if m[i] > i:
                pairs.add((i, m[i]))
        return pairs
    pairs1 = pairs(m1)
    pairs2 = pairs(m2)
    tp = 0.0
    fp = 0.0
    for p in pairs1:
        if p in pairs2:
            tp += 1.0
        else:
            fp += 1.0
    fn = 0.0
    for p in pairs2:
        if not p in pairs1:
            fn += 1.0

    if tp+fp == 0:
        ppv = 1
    else:
        ppv = tp / (tp+fp)

    if tp+fn == 0:
        sens = 1
    else:
        sens = tp / (tp+fn)

    return 1-math.sqrt(ppv*sens)


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


def mfe_inf(pri, db):
    return mfe_dist(pri, db, inv_interaction_network_fidelity)


def avg_inf(pri, db):
    return avg_mfe_dist(pri, db, inv_interaction_network_fidelity)


def min_mfe_inf(pri, db):
    return min_mfe_dist(pri, db, inv_interaction_network_fidelity)

def gc_content(pri):
    gc = 0
    for b in pri:
        if b == 'G' or b == 'C':
            gc += 1
    return gc/len(pri)

def gc_max_control(fit_fn, gc_pcnt, scale=10):
    def fn(pri, db):
        fit = fit_fn(pri, db)
        gc = gc_content(pri)
        return fit - scale*(math.e**max(0, gc-gc_pcnt)-1)/(math.e-1)*abs(fit)
    return fn

def gc_control(fit_fn, gc_pcnt, scale=1):
    def fn(pri, db):
        fit = fit_fn(pri, db)
        gc = gc_content(pri)
        return fit - scale*(math.e**abs(gc-gc_pcnt)-1)/(math.e-1)*abs(fit)
    return fn


def fe_gc_max_50pcnt(pri, db):
    return gc_max_control(free_energy, 0.5)(pri, db)
