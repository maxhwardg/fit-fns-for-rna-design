"""Useful functions for interacting with Vienna RNA"""
import RNA
import numpy as np


class ViennaContext:
    def __init__(self, rna, temp=None, dangles=2, noLPs=False) -> None:
        md = RNA.md()
        md.uniq_ML = 1
        md.dangles = dangles
        md.noLP = noLPs
        if temp is not None:
            md.temperature = temp
        self.fc = RNA.fold_compound(rna, md)
        _, mfe = self.fc.mfe()
        self.fc.exp_params_rescale(mfe)
        self.fc.pf()

    def free_energy(self, ss):
        return self.fc.eval_structure(ss)

    def prob(self, ss):
        return self.fc.pr_structure(ss)

    def make_bppt(self):
        bpp = self.fc.bpp()
        sz = self.fc.length
        res = np.zeros((sz, sz))
        for i in range(sz):
            for j in range(sz):
                if j < i:
                    res[i, j] = bpp[j+1][i+1]
                elif i < j:
                    res[i, j] = bpp[i+1][j+1]
        for i in range(sz):
            res[i, i] = 1-sum(res[i])
        return res

    def subopt(self, energy_delta):
        sub = self.fc.subopt(int(energy_delta*100))
        return [s.structure for s in sub]

    def ensemble_defect(self, ss):
        return self.fc.ensemble_defect(ss)

    def mfe(self):
        return self.fc.mfe()[0]

    def top_2(self):
        samples = []
        delta = 50
        while len(samples) < 2:
            samples = self.fc.subopt(delta)
            delta += 50
            if delta > 1000:
                return (samples[0].structure, None)
        return (samples[0].structure, samples[1].structure)

    def psample(self, samples=1, redundant=True):
        return self.fc.pbacktrack(samples, RNA.PBACKTRACK_DEFAULT if redundant else RNA.PBACKTRACK_NON_REDUNDANT)

