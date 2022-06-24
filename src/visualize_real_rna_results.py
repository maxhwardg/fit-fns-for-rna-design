# Usage is to run python3 src/visualize_real_rna_results.py < results_from_benchmark_fit_fns_on_real_rna_dot_py
import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as font_manager

def main():
    pr = [0]*101
    ed = [0]*101
    bp = [0]*101
    hd = [0]*101
    for ln in sys.stdin:
        if ln.startswith("Result"):
            break
        name, lst = ln.strip().split('[')
        lst = "[" + lst
        lst = eval(lst)
        if lst[0] is None:
            continue
        pr[round(lst[0]*100)] += 1
        ed[round(lst[1]*100)] += 1
        bp[round(lst[2]*100)] += 1
        hd[round(lst[3]*100)] += 1
    x_axis = np.arange(len(pr))

    # Very hacky: just comment/uncomment what you need for the graph
    p = plt.bar(x_axis, pr, color='b')
    # e = plt.bar(x_axis, ed, color='r', alpha=0.5)
    s = plt.bar(x_axis, hd, color='g', alpha=0.5)
    font = {'fontname':'serif'}
    # plt.legend((p, e), ('Probability', 'Ensemble Defect'), prop=font_manager.FontProperties(family='serif'))
    plt.legend((p, s), ('Probability', 'Structure Distance'), prop=font_manager.FontProperties(family='serif'))
    # plt.legend((e,s), ('Ensemble Defect', 'Structure Distance'))

    plt.xlabel('Relative Rank Percentile', fontdict=font)
    plt.ylabel('Number of RNAs', fontdict=font)
    plt.show()

if __name__ == '__main__':
    main()