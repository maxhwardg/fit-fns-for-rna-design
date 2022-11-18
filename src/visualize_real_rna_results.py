# Usage is to run python3 src/visualize_real_rna_results.py < results_from_benchmark_fit_fns_on_real_rna_dot_py
import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as font_manager

def main():
    pr = np.zeros(101)
    ed = np.zeros(101)
    bp = np.zeros(101)
    hd = np.zeros(101)
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

    p = plt.plot(x_axis, np.cumsum(pr), color='b')
    e = plt.plot(x_axis, np.cumsum(ed), color='r')
    s = plt.plot(x_axis, np.cumsum(hd), color='g')
    font = {'fontname':'serif'}
    plt.legend((p, e, s), ('Probability', 'Ensemble Defect' 'Structure Distance'), prop=font_manager.FontProperties(family='serif'))

    plt.xlabel('Relative Rank Percentile', fontdict=font)
    plt.ylabel('Cumulative Number of RNAs', fontdict=font)
    plt.show()

if __name__ == '__main__':
    main()