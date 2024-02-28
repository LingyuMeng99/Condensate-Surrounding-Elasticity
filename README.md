# Simulation codes of condensates with heterogeneous surrounding elasticity

The first code nucleation_simulation.m includes all details of the numerical simulations for the condensate growth model, including the heterogeneous elasticity and active dissolution. All condensates' volumes are recorded in this code. 

The second code BurstCal_Rep.m includes details of bursting kinetics calculation and comparison with the experimental data. It runs nucleation_GrowthTime_pipesimulation.m, which only records the time of condensates' growth and dissolution to save computer memory. The experimental data are in CAST_burst_data.csv, which are from Larsson et al., Genomic encoding of transcriptional burst kinetics, Nature 565, 251–254 (2019). https://doi.org/10.1038/s41586-018-0836-1
