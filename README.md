# Tuning-Multivariable-Controller
Comparison of different tuning approaches for multivariable data using and end-point Model Predictive Controller(MPC) using simulations two biopharmaceutical manufacture: Penicillin and Chinese Hamster Ovary cells production.
The objective of the end-point MPC is to reduce variability by modifing an optimal feed. This work was used to provide the results provided in the article in review "Tuning of latent variable model predictive controllers for bio-pharmaceutical manufacture.", Journal of Process Control. 


To run the controller use the file MPC_main.m in each folder. The folders and subfolders are described below:

- Folder IndPensim: Provides the functions to run the controller for the IndPenSim simulation.
    - Subfolder Results: Stores and plots results.
    - Subfolder Indpensim_Fast_MPC: Provides the functions to run the simulation but calls to the function fast_controller.m() in the middle of the ODE loop. Hence only 1     simulation is necessary to run the MPC with multiple controlpoints.
    - Subfoder Indpensim: Provides the function to run the simulation withou the end-point MPC.
    - Subfolder Tuning_methods: Tuning methods used for the comparison
   
- Folder ChoSim: Provides the functions to run the controller for the ChoSim simulation.
    - Subfolder Results: Stores and plots results.
    - Subfolder CHOcells_Fast_MPC: Provides the functions to run the simulation but calls to the function fast_controller.m() in the middle of the ODE loop. Hence only 1     simulation is necessary to run the MPC with multiple controlpoints.
    - Subfoder CHOcells: Provides the function to run the simulation withou the end-point MPC.
    - Subfolder Tuning_methods: Tuning methods used for the comparison


