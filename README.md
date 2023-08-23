# Sharma_Nitrofurantoin_PBPKmodeling

Read me file for Nitrofurantoin (NFT) PBPK modeling in three different species.
Note that to run the following files you must first install GNU MCSim on your system version 6.1.0. 

There are first 4 model versions:
1.	NFTrabbit_V1.model.R
2.	NFTrabbit_V2.model.R
3.	NFTrabbit_V3.model.R
4.	NFTrabbit_V4.model.R
	
The corresponding 4 simulation files:
1.	NFTrabbit_V1.in.R
2.	NFTrabbit_V2. in.R
3.	NFTrabbit_V3.in.R
4.	NFTrabbit_V4.in.R

The corresponding R simulation files:
1.	NFTrabbit_V1.R
2.	NFTrabbit_V2.R
3.	NFTrabbit_V3.R
4.	NFTrabbit_V4.R

The corresponding R parallel simulation files, that allow to fit the experimental data:
1.	NFTrabbitV1parallel_fitting.R
2.	NFTrabbitV2parallel_fitting.R
3.	NFTrabbitV3parallel_fitting.R
4.	NFTrabbitV4parallel_fitting.R

NFTrabbit_ratV4.model.R, ExtrapolatedV4_ratPosteriorsim.R: these files were used to simulate the rat kinetics using posterior parameters estimated using the rabbit data as an input using 
ExtrapolatedV4_ratEHR.model.R , ExtrapolatedV4_ratEHR.in.R , ExtrapolatedV4_ratEHR.R,  ExtrapolatedV4_ratEHRfitting.R to fit the EHR and simulate the rat kinetics data. 
Extrapolated_ratEHR_rabbit.model.R, Extrapolated_ratEHR_rabbitposteriorsimulation.R:  with these model files one can simulate the rabbit kinetics where enterohepatic recirculation is included but the parameters of that were fitted from rat data.
ExtrapolatedRatehrRabbitRenal_mixed.model.R, ExtrapolatedRatehrRabbitRenal_mixed2.in.R, ExtrapolatedRatehrRabbitRenal_mixed.R and ExtrapolatedRatehrRabbitRenal_mixedfitting.R: with these files one can fit both the rabbit and rat data included both EHR and renal clearance mechanism (we call this in the manuscript ModelV5a).
Extrapolated_ratEHR_rabbit.model.R, Rabbitposteriorsimulation_EHR_CLrmixed.R: files to simulate the rabbit kinetic data using output from the fitting of hierarchical modeling. 

EHRKO_ModelV5a_rat.R and EHRKO_ModelV5a_rabbit.R -> these model files were used to generate the simulation for the enterohepatic recirculation (EHR) knockout scenario. These file use the posterior distribution of parameters generated from “ExtrapolatedRatehrRabbitRenal_mixedfitting.R” and create a simulation file with knockout EHR process (simply make Vehc zero, i.e. parameters for efflux from hepatic to bile)

Rabbit_V5asensitivity.model.R and RabbitNFTPBPK_sensitivityanalysis.R: files used to run the sensitivity analysis. 

Extrapolation to human
Agedynamics_NFThuman.model.R is the model file that includes age based scaling of all the physiological parameters. 
Agedynamic_montecarlo.in.R is the simulation file that includes all the parameters mean value estimated using ModelV4(Renal clearance and rest of the body Partition coefficient) and ModelV5a (EHR related parameters and gut absorption parameters). These paramters that are calculated based on the in-silico partition coefficient algorithm. 
agemonteforward.R is the file that compiles the model file and run the simulation with an input provided in the simulation file (.in file) and also plot the data. 
GFRdiseasedAgedynamic_montecarlo.in.R: this is the simulation file for the GFR compromised scenario and here we create the different GFR conditions mentioned in the paper. 
GFRcompromiseAgeSimulaiton.R: this file uses the compiled model file of “Agedynamics_NFThuman.model.R’  and simulates the various GFR conditions. 

For the manuscript all the generated data from corresponding model file can be found in Posterior_chains (for each model version you will find four chains of posterior distributions) and one posterior simulation as an output. 
The folder “Manuscript_figure” includes all the figures and the code to generate the figures that uses the data available in Posterior chains folder. 

