# Bounding Box-based Multi-objective Bayesian Optimization of Risk Measures under Input Uncertainty (AISTATS2024)

This repository contains the experimental codes written by R of the paper named "Bounding Box-based Multi-objective Bayesian Optimization of Risk Measures under Input Uncertainty."  

This script computes, the computation time of the acquisition function, simple Pareto hypervolume regret, $I^{(i)}_t$ and $I^{(ii)}_t$ for each trial and iteration, and outputs these results as matrix files (AF_computation_time_mat.txt, simple_PHV_regret_mat.txt, inference_e_accuracy_mat.txt and inference_e_accuracy2_mat.txt).

### Description of the Real-world Docking Simulation

In the docking simulation folder, each folder contains additional four documents, Y.csv, input_X.csv, input_W.csv and information.csv. 

The values included in Y.csv are the simulator's calculated docking scores with multiplying minus one. 

The values included in input_X.csv and input_W.csv are the 10-dimensional and 41-dimensioanl physicochemical explanatory variables calculated by the simulator, respectively. 

The information.csv contains identifying names of compound and isomer pairs.
