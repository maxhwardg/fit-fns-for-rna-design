# This script processes a copy and paste of the results into latex table code
import sys


def main():
    names = {"probability": "Probability", "ensemble_defect": "Ensemble Defect",
             "free_energy": "Free Energy", "avg_inf": "Structure Distance (INF; average tie breaking)",
             "avg_mfe_bpd": "Structure Distance (BPD; average tie breaking)",
             "avg_mfe_hd": "Structure Distance (HD; average tie breaking)",
             "mfe_bpd": "Structure Distance (BPD; arbitrary tie breaking)",
             "mfe_hd": "Structure Distance (HD; arbitrary tie breaking)",
             "mfe_inf": "Structure Distance (INF; arbitrary tie breaking)",
             "min_mfe_bpd": "Structure Distance (BPD; minimum tie breaking)",
             "min_mfe_hd": "Structure Distance (HD; minimum tie breaking)",
             "min_mfe_inf": "Structure Distance (INF; minimum tie breaking)",
             }
    records = []
    for ln in sys.stdin:
        toks = ln.strip().split(" ")
        records.append((names[toks[0][:-1]], int(toks[1][:-1]), round(float(toks[2]
                       [:-1]), 2), round(float(toks[3][:-1]), 2), int(toks[4])))
    records.sort()
    for name, correct, rate, gc, unique in records:
        print(f"{name}   &  {correct}  &  {rate:.2f}  &  {gc:.2f}  &  {unique}  \\\\")
        print("\hline")


if __name__ == "__main__":
    main()
