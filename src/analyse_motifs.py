# Takes the output of benchmark_fit_fns_on_synthetic and counts the number of each kind of motif.
# These can be hairpin loops, bulge loops, internal loops, stacks, or multi loops.
# Usage is: python3 src/analyse_motifs.py < results/80nt

import sys


def is_db(s):
    return all(c in ".()" for c in s)


def count_motifs(db):
    stk = [-1]
    par = {}
    m = [i for i in range(len(db))]
    for i in range(len(db)):
        if db[i] == '(':
            stk.append(i)
        elif db[i] == ')':
            l = stk.pop()
            m[i] = l
            m[l] = i
            par[l] = stk[-1]
    children = [[] for _ in range(len(db)+1)]
    for k in par.keys():
        children[par[k]].append(k)
    res = {"pairs": 0, "hairpin": 0, "stack": 0,
           "bulge": 0, "internal": 0, "multi": 0}

    def dfs(l):
        for c in children[l]:
            dfs(c)
        if l != -1:
            res["pairs"] = res["pairs"] + 1
            if len(children[l]) == 0:
                res["hairpin"] = res.get("hairpin")+1
            elif len(children[l]) == 1:
                if children[l][0] == l+1 and m[children[l][0]] == m[l]-1:
                    res["stack"] = res.get("stack")+1
                elif children[l][0] == l+1 or m[children[l][0]] == m[l]-1:
                    res["bulge"] = res.get("bulge")+1
                else:
                    res["internal"] = res.get("internal")+1
            else:
                res["multi"] = res.get("multi")+1
    dfs(-1)
    res["helix"] = res["hairpin"] + res["bulge"] + res["internal"] + res["multi"]
    return res


def main():
    tots = count_motifs("")
    dbs = 0
    nucs = 0
    for ln in sys.stdin:
        toks = ln.strip().split(" ")
        for t in toks:
            if is_db(t):
                dbs += 1
                nucs = len(t)
                res = count_motifs(t)
                for k in res.keys():
                    tots[k] = tots[k] + res[k]
    print("dbs seen: ", dbs)
    print("Raw totals", tots)
    print("Paired percent", tots["pairs"]*2/dbs/nucs)
    for k in tots.keys():
        tots[k] = round(tots[k]/dbs, 2)
    print("Averages", tots)


if __name__ == "__main__":
    main()
