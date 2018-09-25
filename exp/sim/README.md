# Getting Started

1. Modify the first few lines of `./2comp/main.m` to point to your local copy of MIRT.

2. Run `./2comp/main.m`, which should execute within a couple minutes on a typical desktop. This script first simulates MESE and DESS data using a two-compartment model. The script then runs NNLS and RNNLS estimators on the MESE data, followed by ML and PERK estimators on the DESS data. NNLS, RNNLS, and PERK results should be the same as in the paper, but ML results are poorer because a coarse grid search is used by default. 

3. To replicate ML results also, comment out lines 256, 260, 264, 268, 272 and uncomment lines 257, 261, 265, 269, 273 before rerunning. Note however that this will substantially increase ML run time (to about 5 hours in the paper's experiments). 

4. Repeat steps 1-3 for `./3comp/main.m` to investigate estimator robustness to a more realistic three-compartment model. Again, NNLS, RNNLS, and PERK results should be the same as in the paper; to replicate ML results also, comment out lines 245, 249, 253, 257, 261 and uncomment lines 246, 250, 254, 258, 262 before rerunning. 
