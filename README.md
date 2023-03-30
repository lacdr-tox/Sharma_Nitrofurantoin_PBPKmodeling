# Nitrofurantoin_PBPKmodeling


Read me file for NFT PBPK modeling in three different species:
Note that to run the following file you must first install GNU MCSim on your system version 6.1.0. 
There are first 4 version of models 
1.	NFTrabbit_V1.model.R
2.	NFTrabbit_V2.model.R
3.	NFTrabbit_V3.model.R
4.	NFTrabbit_V4.model.R
Corresponding 4 simulation file
1.	NFTrabbit_V1.in.R
2.	NFTrabbit_V2. in.R
3.	NFTrabbit_V3.in.R
4.	NFTrabbit_V4.in.R
Corresponding r simulation file
1.	NFTrabbit_V1.R
2.	NFTrabbit_V2.R
3.	NFTrabbit_V3.R
4.	NFTrabbit_V4.R
Corresponding r parallel simulation file, that allow to fit the experimental data
1.	NFTrabbitV1parallel_fitting.R
2.	NFTrabbitV2parallel_fitting.R
3.	NFTrabbitV3parallel_fitting.R
4.	NFTrabbitV4parallel_fitting.R

NFTrabbit_ratV4.model.R, ExtrapolatedV4_ratPosteriorsim.R , these files were used to simulate the rat kinetics using posterior parameters estimated using the rabbit data as an input using 
ExtrapolatedV4_ratEHR.model.R , ExtrapolatedV4_ratEHR.in.R , ExtrapolatedV4_ratEHR.R,  ExtrapolatedV4_ratEHRfitting.R to fit the EHR and simulate the rat kinetics data. 
Extrapolated_ratEHR_rabbit.model.R, Extrapolated_ratEHR_rabbitposteriorsimulation.R:  this model file is to simulate the rabbit kinetics where entero hepatic process is included but the parameters of that were fitted from rat’s data 
ExtrapolatedRatehrRabbitRenal_mixed.model.R, ExtrapolatedRatehrRabbitRenal_mixed2.in.R, ExtrapolatedRatehrRabbitRenal_mixed.R and ExtrapolatedRatehrRabbitRenal_mixedfitting.R and and file is to fit both the rabbit and rat data included both EHR and renal clearance mechanism (we call this in the manuscript as ModelV5a)
Extrapolated_ratEHR_rabbit.model.R, Rabbitposteriorsimulation_EHR_CLrmixed.R to simulate the rabbit kinetic data using output from the fitting of hierarchical modeling. 

EHRKO_ModelV5a_rat.R and EHRKO_ModelV5a_rabbit.R -> these model files were used to generate the simulation for the enterohepatic recirculation (EHR) knockout scenario. These file uses the posterior distribution of parameters generated from “ExtrapolatedRatehrRabbitRenal_mixedfitting.R” and create a simulation file with knockout EHR process (simply make Vehc zero, parameters for efflux from hepatic to bile)

Rabbit_V5asensitivity.model.R and RabbitNFTPBPK_sensitivityanalysis.R the later file is used to run the sensitivity analysis. 

Extrapolation to human
Agedynamics_NFThuman.model.R is the model file that includes age based scaling of all the physiological parameters. 
Agedynamic_montecarlo.in.R is the simulation file that includes all the parameters mean value estimated using ModelV4(Renal clearance and rest of the body Partition coefficient) and ModelV5a (EHR related parameters and gut absorption parameters). And the paramters that are calculated based on the in-silico partition coefficient algorithm. 
agemonteforward.R is the file that compiled the model file and run the simulation with an input provided in the simulation file (.in file) and also plot the data. 
GFRdiseasedAgedynamic_montecarlo.in.R -> this is the simulation file for the GFR compromised scenario and here we are creating different GFR conditions mentioned in the paper. 
GFRcompromiseAgeSimulaiton.R -> this file used the compiled model file of “Agedynamics_NFThuman.model.R’  and simulate the various GFR condition. 

For the manuscript all the generated data from corresponding model file can be find in Posterior_chains (for each version of model you will find four chains of posterior distribution file) and one posterior simulation as an output. 
The folder “Mansucript_figure” included all the figures and the code to generate the figures that uses the data available in Posterior chains folder. 
