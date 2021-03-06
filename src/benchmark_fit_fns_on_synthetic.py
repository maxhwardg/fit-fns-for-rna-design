import argparse
import fit_fns
import adaptive_walk as aw
import vienna
import common
import random


def make_random_db(len, window):
    pri = common.random_primary(len)
    ctx = vienna.ViennaContext(pri)
    subs = ctx.subopt(window)
    return random.choice(subs)


def gc_pcnt(pri):
    gc = 0
    for c in pri:
        if c == 'G' or c == 'C':
            gc += 1
    return gc / len(pri)


def main():
    ap = argparse.ArgumentParser(
        "Benchmarks using an adaptive walk on various fitness functions using random RNA structure puzzles")
    ap.add_argument("--steps", type=int, default=300,
                    help="Number of mutations to make in the walk")
    ap.add_argument("--processes", type=int, default=1,
                    help="Number of processes to use")
    ap.add_argument("--rnas_per_epoch", type=int, default=32,
                    help="Number of random RNA structures to generate per epoch")
    ap.add_argument("--epochs", type=int, default=50,
                    help="Number of epochs to run")
    ap.add_argument("--window", type=float, default=1.0,
                    help="Window size for suboptimal structures in RNA structure generation")
    ap.add_argument("--rna_length", type=int, default=40,
                    help="Length of RNAs to generate")
    all_fit_fns = dict([('probability', fit_fns.probability), ('ensemble_defect', fit_fns.normalised_ensemble_antidefect),
                        ('free_energy', fit_fns.free_energy), ('probability_diff',
                                                               fit_fns.normalised_prob_diff),
                        ('mfe_bpd', fit_fns.mfe_bpd), ('avg_mfe_bpd',
                                                       fit_fns.avg_mfe_bpd), ('min_mfe_bpd', fit_fns.min_mfe_bpd),
                        ('mfe_hd', fit_fns.mfe_hd), ('avg_mfe_hd', fit_fns.avg_mfe_hd), ('min_mfe_hd', fit_fns.min_mfe_hd)])
    ap.add_argument("--fitness_fns", type=str, default=[], action='append',
                    help=f'Fitness functions to use. Leave empty for all. Should be from {list(all_fit_fns.keys())}')
    ap.add_argument("--seed", type=int, default=None,
                    help="Random seed. Leave empty for no seed.")

    args = ap.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    if args.fitness_fns == []:
        args.fitness_fns = list(all_fit_fns.keys())

    fitness_fns = {}
    for fn in args.fitness_fns:
        fitness_fns[fn] = all_fit_fns[fn]

    gc_pcnt_sum = {}
    fit_fn_solves = {}
    unique_solves = {}
    for name in fitness_fns.keys():
        gc_pcnt_sum[name] = 0
        fit_fn_solves[name] = 0
        unique_solves[name] = 0
    total_seqs = 0
    epoch_dbs = [[make_random_db(args.rna_length, args.window) for _ in range(
        args.rnas_per_epoch)] for _ in range(args.epochs)]
    for epoch in range(args.epochs):
        dbs = epoch_dbs[epoch]
        db_solved = []
        for _ in dbs:
            db_solved.append([])
        for fn_name, fit_fn in fitness_fns.items():
            res = aw.parallel_walk(
                dbs, fit_fn, aw.random_batch_guide, args.steps, processes=args.processes)
            for i in range(len(res)):
                ctx = vienna.ViennaContext(res[i])
                subs = ctx.subopt(0)
                gc_pcnt_sum[fn_name] = gc_pcnt_sum.get(
                    fn_name, 0) + gc_pcnt(res[i])
                if len(subs) == 1 and subs[0] == dbs[i]:
                    db_solved[i].append(fn_name)
                    fit_fn_solves[fn_name] = fit_fn_solves.get(
                        fn_name, 0) + 1
        print("Epoch {}: {}".format(epoch, fit_fn_solves))
        for i, fns in enumerate(db_solved):
            total_seqs += 1
            if len(fns) == 1:
                unique_solves[fns[0]] = unique_solves.get(fns[0], 0) + 1
            print('\t', dbs[i], fns)

    print(f"Total Sequences: {total_seqs}")
    print("Format is \"Function name: #solves, solve ratio, GC percent, unique solves\"")
    for fn_name, solves in fit_fn_solves.items():
        print(
            f"{fn_name}: {solves}, {solves / total_seqs}, {gc_pcnt_sum[fn_name] / total_seqs}, {unique_solves[fn_name]}")


if __name__ == '__main__':
    main()
