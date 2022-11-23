# Usage is python3 src/visualize_real_rna_results.py < results/real_rna_bench
import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as font_manager

def main():
    buckets = 100
    pr = np.zeros(buckets+1)
    ed = np.zeros(buckets+1)
    bp = np.zeros(buckets+1)
    hd = np.zeros(buckets+1)
    for ln in sys.stdin:
        if ln.startswith("Result") or ln.find(":") != -1:
            break
        name, lst = ln.strip().split('[')
        lst = "[" + lst
        lst = eval(lst)
        if lst[0] is None:
            continue
        pr[round(lst[0]*buckets)] += 1
        ed[round(lst[1]*buckets)] += 1
        bp[round(lst[2]*buckets)] += 1
        hd[round(lst[3]*buckets)] += 1
    x_axis = np.arange(len(pr))

    p, = plt.plot(x_axis, np.cumsum(pr), color='b')
    e, = plt.plot(x_axis, np.cumsum(ed), color='r')
    s, = plt.plot(x_axis, np.cumsum(hd), color='g')
    font = {'fontname':'serif'}
    plt.legend((p, e, s), ('Probability', 'Ensemble Defect', 'Structure Distance'), prop=font_manager.FontProperties(family='serif'))

    plt.xlabel('Relative Rank Percentile', fontdict=font)
    plt.ylabel('Cumulative Number of RNAs', fontdict=font)
    plt.show()

if __name__ == '__main__':
    main()