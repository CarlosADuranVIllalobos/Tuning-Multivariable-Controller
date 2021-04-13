v# Tuning-Multivariable-Controller
Comparison of different tuning approaches for multivariable data

The code is divided in 2 folders: MPC(Model Predictive Control) and B2B(Batch-to-Batch Optimisation). The objective of the B2B optimisation is to obtain an optimal feeds starting with suboptimal initial feeds. The objective of the MPC is to reduce variability by modifing an optimal feed. The subfolders, scripts and functions are described below:

- Subfolder Indpensim: Provides the functions to run the simulation.

- Subfolder Indepmsim_Fast_MPC: Provides the functions to run the simulation but calls to the function fast_controller.m() in the middle of the ODE loop. Hence only 1 simulation is necessary to run the MPC with multiple controlpoints.

- Subfolder Results: Stores results.

- Script create_training.m: Create training batches 

- Script create_nominal.m: Create nominal batches 

- Script plot_results.m: Plot results to compare approaches

- Function convfactor.m: Function to convert indpensim output data to a struct of batches.

- Function pls.m: PLS algorithm 

- Function plsmccv: Monte-Carlo cross-validation using PLS 


Scripts and functions that you may need to change to try different approaches:

- Script B2B_main.m: Runs B2B optimisations of n-replicates.

- Script MPC_main.m: Runs MPC and displays the outputs at each controlpoint (Runs multiple simulations at each batch).

- Script MPC_main_fast.m: Runs MPC and display the outputs at each batch (Runs only 2 simulations for each batch).

- Function mvtoptu1.m: Optimisation of the MVT using a PLS model and missing data algorithms.

- Function fast_controller.m: Optimisation of the MVT using a PLS model and missing data algorithms that is called inside the ODE loop.

- Function orgplsmpc: Organizes data PLS matrices based on the location of the control point to optimise the MVT


