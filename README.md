# Source code and data for "Fitness Functions for RNA Structure Design"
All scripts are written assuming Python 3, and were tested using Python 3.8.10.
They also assume that the Python bindings for [Vienna RNA](https://www.tbi.univie.ac.at/RNA/) 2.5.0 are installed.

There are three main programs.

To run the benchmark on synthetic RNA use:
```
python3 src/benchmark_fit_fns_on_synthetic.py
```

To run on real RNA use:
```
python3 src/benchmark_fit_fns_on_real_rna.py > results
```

Then to make the graphs for the result:
```
python3 src/visualize_real_rna_results.py < results
```

Note that these scripts can be very slow. Using `--processes=X` can help to use multiple cores.
Also, command line arguments are used to modify the programs. Please see the help messages.