getwd()
setwd("~/MCSim/mod")
library(ggmcmc)
library(tidyverse)
library(rstan)
library(StanHeaders)
library(bayesplot)
library("ggformula")
rm(list = ls())
# run different batches testing
parallel::detectCores()
library(rstudioapi)

Modeling_dir= "NFT_PBPK2"

set.seed(1234)
jobRunScript(workingDir = getwd(), path = "NFT_PBPK2/NFTrabbitV4.R", importEnv = T, exportEnv = "job_5")
set.seed(4234)
jobRunScript(workingDir = getwd(), path = "NFT_PBPK2/NFTrabbitV4.R", importEnv = T, exportEnv = "job_6")
set.seed(5234)
jobRunScript(workingDir = getwd(), path = "NFT_PBPK2/NFTrabbitV4.R", importEnv = T, exportEnv = "job_7")
set.seed(6234)
jobRunScript(workingDir = getwd(), path = "NFT_PBPK2/NFTrabbitV4.R", importEnv = T, exportEnv = "job_8")


makemcsim <- function(model,
                      deSolve = F,
                      dir = Modeling_dir) {
  exe_file <- paste0("mcsim_", model, " ")
  
  if (file.exists(Modeling_dir) == F) {
    dir.create(Modeling_dir)
  }
  
  if (file.exists(model)) {
    # Move model file from working directory to Modeling folder
    invisible(file.copy(
      from = paste0(getwd(), "/", model),
      to = paste0(getwd(), "/",Modeling_dir,"/", model)
    ))
    invisible(file.remove(model))
    message(paste0("* The '", model, "' has been moved to Modeling folder."))
  }
  
  if (deSolve == T) {
    system(paste("./MCSim/mod -R ", dir, "/", model, " ", model, ".c", sep = ""))
    system (paste0("R CMD SHLIB ", model, ".c")) # create *.dll files
    dyn.load(paste(model, .Platform$dynlib.ext, sep = "")) # load *.dll
    source(paste0(model, "_inits.R"))
  } else {
    system(paste("../mod/mod ", dir, "/", model, " ", model, ".c", sep = ""))
    system(
      paste(
        "gcc -O3 -I.. -I../sim -o mcsim_",
        model,
        " ",
        model,
        ".c ../sim/*.c -lm ",
        sep = ""
      )
    )
    invisible(file.remove(paste0(model, ".c")))
    if (file.exists(exe_file))
      message(paste0("* Created executable program '", exe_file, "'."))
  }
}


mcsim <- function(model,
                  input,
                  Setsim,  #extra take setpoint file to mod file
                  dir = Modeling_dir,
                  parallel = F) {
  if (file.exists(Modeling_dir) == F) {
    dir.create(Modeling_dir)
  }
  
  exc = paste0("mcsim_", model, "")
  if (file.exists(exc) == F) {
    makemcsim(model, dir = dir)
    if (file.exists(exc) == F) {
      stop("* '", exc, "' is not exist .")
    }
  }
  
  if (file.exists(input)) {
    invisible(file.copy(
      from = paste0(getwd(), "/", input),
      to = paste0(getwd(), "/Modeling/", input)
    ))
    invisible(file.remove(input))
  }
  
  
  tx  <- readLines(paste0(dir, "/", input))
  MCMC_line <- grep("MCMC \\(", x = tx)
  MonteCarlo_line <- grep("MonteCarlo \\(", x = tx)
  SetPoints_line <- grep("SetPoints \\(", x = tx)
  
  if (length(MCMC_line) != 0) {
    #file_defore <- list.files()
    RandomSeed <- exp(runif(1, min = 0, max = log(2147483646.0)))
    tx2 <-
      gsub(pattern = "10101010",
           replace = paste(RandomSeed),
           x = tx)
    checkfile <- "MCMC.check.out"
    
    if (file.exists(checkfile)) {
      file.remove(checkfile)
    }
    
    if (parallel == T) {
      i <- sample(1111:9999, 1)
      name <- gsub("\\..*", "", input)
      mcmc_input <- paste0(name, "_", i, ".in")
      mcmc_output <- paste0(name, "_", i, ".out")
      tx3 <-
        gsub(pattern = "MCMC.default.out",
             replace = mcmc_output,
             x = tx2)
      writeLines(tx3, con = mcmc_input)
      system(paste("./mcsim_", model, " ", mcmc_input, sep = ""))
      
    } else{
      tmp <- "tmp_mcmc.in.R"
      writeLines(tx, con = paste0(dir, "/", input))
      writeLines(tx2, con = paste0(dir, "/", tmp))
      system(paste("./mcsim_", model, " ", tmp, sep = ""))
      outfile <- "MCMC.default.out"
      tx2 <- gsub(pattern = ",0,",
                  replace = ",1,",
                  x = tx)
      tx3 <- gsub(
        pattern = paste0("\"", outfile, "\",\"\""),
        replace = paste0("\"", checkfile, "\",\"", outfile, "\""),
        x = tx2
      )
      writeLines(tx3, con = paste0(dir, "/", tmp))
      
      system(paste("./mcsim_", model, " ", dir, "/", tmp, sep = ""))
      file.remove(paste0(dir, "/", tmp))
    }
    
    if (file.exists(checkfile)) {
      message(paste0("* Create '", checkfile, "' from the last iteration."))
    }
    
    if (parallel == T) {
      df <- read.delim(mcmc_output)
    } else {
      df <- read.delim("MCMC.default.out")
    }
    
  } else if (length(MonteCarlo_line) != 0) {
    RandomSeed <- runif(1, 0, 2147483646)
    tx2 <-
      gsub(pattern = "10101010",
           replace = paste(RandomSeed),
           x = tx)
    writeLines(tx2, con = paste0(dir, "/", input))
    message(paste("Execute:", " ./mcsim_", model, " ", dir, "/", input, sep = ""))
    
    system(paste("./mcsim_", model, " ", dir, "/", input, sep = ""))
    writeLines(tx, con = paste0(dir, "/", input))
    df <- read.delim("simmc.out")
  } else if (length(SetPoints_line) != 0) {
    invisible(file.copy(
      from = paste0(getwd(), "/",Modeling_dir,"/", Setsim),
      to = paste0(getwd(), "/", Setsim)
    ))
    message(paste("Execute:", " ./mcsim_", model, " ", dir, "/", input, sep = ""))
    system(paste("./mcsim_", model, " ", dir, "/", input, sep = ""))
    # df <- read.delim("simmc.out")
    df <- read.delim(paste0(ouput))
  } else {
    message(paste("Execute:", " ./mcsim_", model, " ", dir, "/", input, sep = ""))
    system(paste("./mcsim_", model, " ", dir, "/", input, sep = ""))
    df <- read.delim(paste0(ouput), skip = 1)
  }
  return(df)
}

mcmc_array <- function(data, start_sampling = 0) {
  n_chains <- length(data)
  sample_number <- dim(data[[1]])[1] - start_sampling
  dim <- c(sample_number, n_chains, dim(data[[1]])[2])
  n_iter <- dim(data[[1]])[1]
  n_param <- dim(data[[1]])[2]
  x <- array(sample_number:(n_iter * n_chains * n_param), dim = dim)
  for (i in 1:n_chains) {
    x[, i,] <- as.matrix(data[[i]][(start_sampling + 1):n_iter,])
  }
  dimnames(x)[[3]] <- names(data[[1]])
  x
}
clear <- function() {
  files <- c(
    dir(pattern = c("*.out")),
    dir(pattern = c("sim.in")),
    dir(pattern = c("*.R")),
    dir(pattern = c("*.R.so")),
    dir(pattern = c("*.R.o")),
    dir(pattern = c("*.R.dll")),
    dir(pattern = c("*.R.c")),
    dir(pattern = c("*.R_inits.R")),
    dir(pattern = c("*.perks"))
  )
  invisible(file.remove(files))
}

Species = "Rabbit"
modelingversion = "V4"
m  = 1   #Job number starts
n = 5   #job number ends

jobs = c(noquote(paste0("job", "_", seq(m,n, by = 1),"$","out", collapse = ",")))
sims <- mcmc_array(list(job_5$out,job_6$out,job_7$out, job_8$out))
#sims <- mcmc_array(list(job_5$out,job_6$out, job_8$out))
#sims <- mcmc_array(list(job_7$out))
endparms = length(names(job_8$out))-3
parms_name_1 <- paste((names(job_8$out)[2:endparms]),sep = ",")

# Files <- list.files(path = paste0(getwd(), "/", "NFT_PBPK2/Posterior_chains","/"),pattern = "NFTrabbit_V4.*\\.out$")
# Files = list.files(path = paste0(getwd()),pattern = "^NFT_SB1.*\\.out$",full.names=TRUE)


#myfiles = lapply(Files, read.delim)
#sims <- mcmc_array(list(myfiles[[1]],myfiles[[2]],myfiles[[3]],myfiles[[4]]),start_sampling = 0)
# endparms = length(names(myfiles[[1]]))-3
# parms_name_1 <- paste((names(myfiles[[1]])[2:endparms]),sep = ",")

str <- ceiling(nrow(sims)/2) + 1
end <- nrow(sims)
j <- c(str:end)
length(j)
color_scheme_set("mix-blue-red")
mcmc_trace(sims[j,,], pars = parms_name_1) + ggplot2::scale_color_discrete()#, facet_args = list(ncol = 1, strip.position = "left"))
# monitor(sims[,,parms_name_1], digit=4)
# 
monitor_file = (monitor(sims[,,parms_name_1], digit=4))

mcmc_dens_overlay(x = sims[j,,], pars = parms_name_1)
#mcmc_dens_overlay(x = sims[j,,], pars = parms_name)

mcmc_dens_overlay(x = sims[j,,], pars = "LnData")

dataf = data.frame((monitor_file))
str(dataf)
str(monitor_file)
write.csv(dataf[ ,c("mean","sd","X2.5.","X97.5.","Rhat")], file = paste0(getwd(),"/",Modeling_dir, "/", "Posterior_chains","/", Species,modelingversion, "_NFTPBPK.csv"))

X <- sims[j,,] %>% matrix(nrow = length(j)*4)

parms_name_2 <- paste((names(job_8$out)),sep = ",")



filenames <- list.files(path = paste0(getwd(), "/"),pattern = "NFTrabbit_V4.*\\.out$")
file.copy(from = paste0(getwd(), "/", filenames),
          to = paste0(getwd(), "/",Modeling_dir, "/", "Posterior_chains","/", filenames))

colnames(X) = parms_name_2
head(X)
write.table(X, file = paste0(getwd(),"/",Modeling_dir, "/", "Posterior_chains","/",Species,modelingversion, "posterior_hierar1.out"), row.names = F, sep = "\t")

posteriordisfile <- read.delim(paste0(getwd(),"/",Modeling_dir, "/", "Posterior_chains","/",Species,modelingversion,"posterior_hierar1.out"), header = TRUE, sep = "")

head(posteriordisfile)
params = data.frame(colnames(posteriordisfile))
colnames(params) = "parametersname"

file = params %>% separate(col = parametersname, c("Parms", "Subj", "levels"),sep = "([.])",fill = "right")
drug = c("NFT")
drugs = drug[1]
simulation = c("general", "hierar")[1]
level = ifelse(drugs == "NFT",1,"NA") #here no role of level as this is non heirarchical methods 
a = file %>% filter(Subj ==1)
b = file %>% filter(Subj ==1 & levels == "")
c = file %>% filter(Subj ==1&levels ==level)
#c$levels = paste0(c$levels, ".")  #to add a dot as it is required
d = c(Parmshier = "iter",Parms = "iter", Subj = "", levels = " ")
e = c(Parmshier = "LnPrior",Parms = "LnPrior", Subj = "", levels = " ")
f = c(Parmshier = "LnData",Parms = "LnData", Subj = "", levels = " ")
g = c(Parmshier = "LnPosterior",Parms = "LnPosterior", Subj = "", levels = " ")
file0 = rbind(c,b)
fileg = b
file1 = if(simulation == "hierar"){file0} else {fileg}  # if simulaiton is heirachical then rbind otherwise just keep c (without levels)
file2 =distinct(file1,Parms,Subj, .keep_all= TRUE) %>% unite(.,"parmshier", Parms:levels, sep= ".", remove = FALSE) %>% rbind(d,.,e,f,g)
filename = file2$parmshier  #
randomdata = sample_n(posteriordisfile, 500)
hierSul_setpoint = randomdata %>% select(filename) %>% write.table(., file = paste0(getwd(),"/", Modeling_dir,"/",drugs, Modeling_dir,"setSim", ".out"), row.names = F, sep = "\t")
filesetout <- read.delim(paste0(getwd(),"/", Modeling_dir,"/",drugs, Modeling_dir,"setSim", ".out"), header = TRUE, sep = "")

#check if names are matching or not
colnames(filesetout)


model <- "NFTrabbit_V4.model.R"


Setsim = paste0(drugs, Modeling_dir,"setSim", ".out")
#chemical specific inputs and outputs
input = paste0(drugs, Modeling_dir,"set.in", ".R")  # this file will contain the chemical inputs. 
ouput = paste0(drugs, Modeling_dir,"set", ".out")  #check this output file name if it is the same or not
paramsSet = paste(file2$Parms[-c(1,length(file2$Parms)-2,(length(file2$Parms)-1), length(file2$Parms))], collapse = ",") 


### parameters in the output file need to be in simulation file with inbetween space
#gsub(",([A-Za-z])", ", \\1", .)  #to add a space after comma 

setsimparm = paste0( 'SetPoints ','(', '"', ouput , '"',',','"', Setsim ,'"', ',' , 0 , ',', paramsSet,')',';') %>% gsub(',([A-Za-z])', ', \\1', .)


##########
#experimental data reading
####################
#creat data section
####################

file_name = "Rabbit_PKdata.csv"
myfiles <- read.csv(paste0(getwd(), "/", Modeling_dir, "/", file_name))
head(myfiles)
#Timeid = paste0(unique(myfiles$Time, incomparables = FALSE),sep = ",",collapse = " ") %>% gsub(",$", "", .)  #gsub to remove the last comma


# Times = paste(0, Timeid,sep = ",")   # to include time zero which does not exist in experiment data
# Final_Rdose = unique(myfiles$dose_uM, incomparables = FALSE)
# Replicates = unique(myfiles$replID, incomparables = FALSE)



proteins = c("cplasma","P_excreted","Massbalance","Filteration","Reabsorption","Secretion","Total_renal_Cl","GFR_contribution","Secr_contribution","Reabs_contribution","Re_GFR_contribution","cgut", "cliver", "cfat","ckidney","cfilterate", "crestbody", "Agutlumen","Adelay","Aurine") 

exposure_route = c("IV", "oraldose")
exp_route = exposure_route[1]

#remove the input file if alread exist and create new one
file.remove(paste0(getwd(),"/",Modeling_dir, "/",input))

# Final_Rdose2 = c(0.5,1.25, 2.5, 5, 10, 15)
# Final_Rdose3 = c(1.5, 10)   #oral dose per kg in rabbit

if (exp_route =="IV"){
  
  Final_Rdose1 = c(0.5,1.25, 2.5, 5, 10, 15)
  
} else {
  
  Final_Rdose1 = c(1.25, 10)
  
}
Times = paste(0, 4, 0.01,sep = ",") 
simulationprint = function (dose,biomarker1,biomarker2,biomarker3, biomarker4, biomarker5, biomarker6, biomarker7, biomarker8,biomarker9, biomarker10,biomarker11,biomarker12,biomarker13,biomarker14,biomarker15,biomarker16,biomarker17,biomarker18,biomarker19,biomarker20,Times){
  
  variable1 = ifelse(biomarker1 == F, "#", paste0( "PrintStep", "(", biomarker1,",",Times, ")",";"))
  variable2 = ifelse(biomarker2 == F, "#", paste0( "PrintStep", "(", biomarker2,",",Times, ")",";"))
  variable3 = ifelse(biomarker3 == F, "#", paste0( "PrintStep", "(", biomarker3,",",Times, ")",";"))
  variable4 = ifelse(biomarker4 == F, "#", paste0( "PrintStep", "(", biomarker4,",",Times, ")",";"))
  variable5 = ifelse(biomarker5 =="", "#", paste0( "PrintStep", "(", biomarker5,",",Times, ")",";"))
  variable6 = ifelse(biomarker6 == "", "#", paste0( "PrintStep", "(", biomarker6,",",Times, ")",";"))
  variable7 = ifelse(biomarker7 == "", "#", paste0( "PrintStep", "(", biomarker7,",",Times, ")",";"))
  variable8 = ifelse(biomarker8 == "", "#", paste0( "PrintStep", "(", biomarker8,",",Times, ")",";"))
  variable9 = ifelse(biomarker9 == "", "#", paste0( "PrintStep", "(", biomarker9,",",Times, ")",";"))
  variable10 = ifelse(biomarker10 == "", "#", paste0( "PrintStep", "(", biomarker10,",",Times, ")",";"))
  variable11 = ifelse(biomarker11 == "", "#", paste0( "PrintStep", "(", biomarker11,",",Times, ")",";"))
  variable12 = ifelse(biomarker12 == "", "#", paste0( "PrintStep", "(", biomarker12,",",Times, ")",";"))
  variable13 = ifelse(biomarker13 == "", "#", paste0( "PrintStep", "(", biomarker13,",",Times, ")",";"))
  variable14 = ifelse(biomarker14 == "", "#", paste0( "PrintStep", "(", biomarker14,",",Times, ")",";"))
  variable15 = ifelse(biomarker15 == "", "#", paste0( "PrintStep", "(", biomarker15,",",Times, ")",";"))
  
  variable16 = ifelse(biomarker16 == "", "#", paste0( "PrintStep", "(", biomarker16,",",Times, ")",";"))
  variable17 = ifelse(biomarker17 == "", "#", paste0( "PrintStep", "(", biomarker17,",",Times, ")",";"))
  variable18 = ifelse(biomarker18 == "", "#", paste0( "PrintStep", "(", biomarker18,",",Times, ")",";"))
  variable19 = ifelse(biomarker19 == "", "#", paste0( "PrintStep", "(", biomarker19,",",Times, ")",";"))
  variable20 = ifelse(biomarker20 == "", "#", paste0( "PrintStep", "(", biomarker20,",",Times, ")",";"))
  
  Exposure = ifelse(exp_route == "IV", paste0( "IVdosing", "=", dose, ";"),paste0("oraldose", "=", dose,";", "\n", "oraldose1", "=","PerDose", "(", dose, ",",24,",", 0,",", 0.001,")", ";"))
  
  
  return(paste0("Simulation", "{" , "\n", Exposure ,"\n",
                variable1,"\n",variable2,"\n",variable3,"\n",variable4,"\n",variable5,"\n",variable6,"\n",variable7,
                "\n",variable8,
                "\n",variable9,
                "\n",variable10,
                "\n",variable11,
                "\n",variable12,
                "\n",variable13,
                "\n",variable14,
                "\n",variable15,
                "\n",variable16,
                "\n",variable17,
                "\n",variable18,
                "\n",variable19,
                "\n",variable20,
                "\n","}"))
  
}


f = file(paste0(getwd(),"/",Modeling_dir, "/",input), open = 'a')

write(setsimparm, f)


for (i in 1:length(Final_Rdose1)){
  
  write(simulationprint(dose = Final_Rdose1[i],biomarker1 = proteins[1],biomarker2 = proteins[2],biomarker3 = proteins[3],
                        biomarker4 = proteins[4],biomarker5 = proteins[5], biomarker6 = proteins[6],biomarker7 = proteins[7],
                        biomarker8 = proteins[8],biomarker9 = proteins[9],
                        biomarker10 = proteins[10],biomarker11 = proteins[11],
                        biomarker12 = proteins[12],biomarker13 = proteins[13],
                        biomarker14 = proteins[14],biomarker15 = proteins[15],
                        biomarker16 = proteins[16],biomarker17 = proteins[17],
                        biomarker18 = proteins[18],biomarker19 = proteins[19],biomarker20 = proteins[20],
                        Times), f)
  
}



# write("End.", f)
close(f)



model <- model
makemcsim(model = model, deSolve = F, dir = Modeling_dir)
X_setpts <- mcsim(model, input, Setsim = Setsim)
invisible(file.remove(Setsim))    # to remove the output file
X_setpts1 = X_setpts



proteins = c("cplasma","P_excreted","Massbalance","Filteration","Reabsorption","Secretion","Total_renal_Cl","GFR_contribution","Secr_contribution","Reabs_contribution","Re_GFR_contribution","cgut", "cliver", "cfat","ckidney","cfilterate", "crestbody", "Agutlumen","Adelay","Aurine") 

# Keap1modified = grep(pattern = paste0("^",proteins[4],collapse = ""),vars)
# names(X_setpts1)[Keap1modified] = "modifiedKeap1"
#updated the vars again
vars <- names(X_setpts1)

Simulationfile = data.frame()
#udpate the Keap1modified

for (i in 1:length(proteins)){
  finder =  grep(pattern = paste0("^",proteins[i],collapse = ""),vars)
  file <- apply(X_setpts1[finder], 2, quantile,  c(0.5, 0.025, 0.975)) %>% t()
  x = data.frame(file)
  X_mean <- apply(X_setpts1[finder],2, mean)   # can apply for all the doses
  x$mean = X_mean
  Simulationfile <- rbind(Simulationfile,x)
}
nrow(Simulationfile)

simul = ifelse(exp_route =="IV",6, 2) #if we include dose zero to check the steady state
#rename in captial as experimental data 
biomarker = length(proteins)
Times = seq(0,4, by = 0.01) #experimental data time point
Simulationfile$Simulation = (rep(paste0(rep("Simu",simul),1:simul),each = length(Times)*biomarker))  # as we have 3 variables
Simulationfile$dose_uM = as.numeric(rep(paste0(rep(Final_Rdose1)),each = length(Times)))
Simulationfile$Chemical = rep(paste0("NFT"),each = nrow(Simulationfile))
Simulationfile$Times <- Times
Simulationfile$variables = (rep(paste0(rep(proteins)),each = length(Times)*simul))
colnames(Simulationfile) <- c("median", "LCL", "UCL","mean","Simulations","dose_uM","Chemicals","Times","variables")

Expdata = myfiles 

write.table(Simulationfile, file = paste0(getwd(),"/",Modeling_dir, "/", "Posterior_chains","/",Species,modelingversion, "Simulation.out"), row.names = F, sep = "\t")

# Species = "Rabbit"
# modelingversion = "V4" 
# Simulationfile =read.delim(paste0(getwd(),"/",Modeling_dir, "/", "Posterior_chains","/",Species,modelingversion, "Simulation.out"), header = TRUE, sep = "")
# head(Simulation)
Ylab = c("Plasma conc.(mg/L)", "Percentage of Dose excreted","Cumulative urine (mg)", "Total Mass Balance check (mg)", "Kidney conc.(mg/L)","Liver conc.(mg/L)")
Title1 = c("NFT concentrations in Rabbit following a single intravenous dose","Dose normalized cumulative urine excretion in Rabbit following a single IV",
           "Cumulative urine excretion in Rabbit following a single IV","NFT concentrations in Rabbit following a single OD",
           "Cumulative urine excretion in Rabbit following a single OD")


plot_function1 = function (variable1,variable2, route, mainplot, labplot) {
  Simulation = Simulationfile %>% filter(.,variables == variable1)
  #Simulation = Simulationfile %>% filter(.,variables == variable & Simulations != "Simu1")
  exp_data = Expdata %>% filter(Variable == variable2) %>% filter(exposure_route == route)
  
  ggplot() +
    geom_line(data =Simulation,aes(x = Times, y = mean, color="Simulated_mean"),size=1) +
    geom_line(data =Simulation,aes(x = Times, y = LCL, color="Simulated_p2.5"),size=1) +
    geom_line(data =Simulation,aes(x = Times, y = UCL, color="Simulated_p97.5"),size=1) +
    geom_point(data =exp_data,aes(x = Time, y = cplasma, color="Exp_mean"), size = 1.2) +
    # geom_line(aes(x = Times, y = Nrf2max_observed, color="Exp_max"),linetype=2) +
    # geom_line(aes(x = Times, y = Nrf2min_observed, color="Exp_min"),linetype=2) +
    # geom_errorbar(aes(x = Times,ymin=Nrf2_observed-Nrf2_SD, ymax= Nrf2_observed +Nrf2_SD,color = "Experiments_mean(Sd)"), size=0.5,
    #               width=.25)+
    facet_wrap(~dose_uM) +
    labs(x = "Time (hr)", y = labplot, sec.x="Exposure (mg/kg)",
         title = mainplot) +
    scale_y_continuous(trans = 'log10')+
    # scale_y_continuous(limits=c(0,3)) +
    scale_color_manual(name = "type",
                       breaks = c("Simulated_mean","Simulated_p2.5","Simulated_p97.5", "Exp_mean"),
                       values = c("Simulated_mean" = "blue","Simulated_p2.5" = "green","Simulated_p97.5" = "red","Exp_mean" = "blue"))+
    theme(legend.title = element_blank()) + theme(legend.position="bottom") +
    theme(legend.text = element_text(size= 8,face="bold")) 
}

Title = ifelse(exp_route == "oraldose", Title1[4], Title1[1])

plot_function1(variable1 = "cplasma",variable2 = "cplasma",route = exp_route,mainplot =Title, labplot = Ylab[1])

plot_function2 = function (variable1,variable2, route, mainplot, labplot) {
  Simulation = Simulationfile %>% filter(.,variables == variable1)
  #Simulation = Simulationfile %>% filter(.,variables == variable & Simulations != "Simu1")
  exp_data = Expdata %>% filter(Variable == variable2) %>% filter(exposure_route == route)
  
  ggplot() +
    geom_line(data =Simulation,aes(x = Times, y = mean, color="Simulated_mean"),size=1) +
    geom_line(data =Simulation,aes(x = Times, y = LCL, color="Simulated_p2.5"),size=1) +
    geom_line(data =Simulation,aes(x = Times, y = UCL, color="Simulated_p97.5"),size=1) +
    geom_point(data =exp_data,aes(x = Time, y = cplasma, color="Exp_mean"), size = 1.2) +
    # geom_line(aes(x = Times, y = Nrf2max_observed, color="Exp_max"),linetype=2) +
    # geom_line(aes(x = Times, y = Nrf2min_observed, color="Exp_min"),linetype=2) +
    # geom_errorbar(aes(x = Times,ymin=Nrf2_observed-Nrf2_SD, ymax= Nrf2_observed +Nrf2_SD,color = "Experiments_mean(Sd)"), size=0.5,
    #               width=.25)+
    facet_wrap(~dose_uM) +
    labs(x = "Time (hr)", y = labplot, sec.x="Exposure (mg/kg)",
         title = mainplot) +
    scale_y_continuous(trans = 'log10')+
    # scale_y_continuous(limits=c(0,3)) +
    scale_color_manual(name = "type",
                       breaks = c("Simulated_mean","Simulated_p2.5","Simulated_p97.5", "Exp_mean"),
                       values = c("Simulated_mean" = "blue","Simulated_p2.5" = "green","Simulated_p97.5" = "red","Exp_mean" = "blue"))+
    theme(legend.title = element_blank()) + theme(legend.position="bottom") +
    theme(legend.text = element_text(size= 8,face="bold")) 
}

Title = ifelse(exp_route == "oraldose", Title1[4], Title1[1])

plot_function2(variable1 = "P_excreted",variable2 = "P_excreted",route = exp_route,mainplot = Title1[2], labplot = Ylab[2])



plot_function3 = function (variable1,dose,mainplot,labplot) {
  Simulation = Simulationfile %>% filter(.,variables %in% variable1) %>% filter(., dose_uM %in% dose) 
  #Simulation = Simulationfile %>% filter(.,variables == variable & Simulations != "Simu1")
  #exp_data = Expdata %>% filter(Variable == variable2)
  
  ggplot() +
    geom_line(data =Simulation,aes(x = Times, y = mean, color= as.factor(dose_uM))) +
    # geom_line(data =Simulation,aes(x = Times, y = LCL, color="Simulated_p2.5"),size=1) +
    # geom_line(data =Simulation,aes(x = Times, y = UCL, color="Simulated_p97.5"),size=1) +
    
    facet_wrap(~variables, scales = "free") +
    labs(x = "Time (hr)", y = labplot, sec.x="Exposure (mg/kg)",
         title = mainplot) +
    #scale_y_continuous(trans = 'log10')+
    # scale_y_continuous(limits=c(0,3)) +
    # scale_color_manual(name = "type",
    #                    breaks = c("Simulated_mean","Simulated_p2.5","Simulated_p97.5"),
    #                    values = c("Simulated_mean" = "blue","Simulated_p2.5" = "green","Simulated_p97.5" = "red"))+
    theme(legend.title = element_blank()) + theme(legend.position="bottom") +
    theme(legend.text = element_text(size= 8,face="bold")) 
}

plot_function3(variable1 = proteins[7:11],dose = Final_Rdose1, mainplot = Title1[5] ,labplot = "NFT conc.(mg/L)")


hysteresis = Simulationfile  %>% select(mean, dose_uM,Times, variables) 
files = pivot_wider(hysteresis, names_from = variables, values_from = mean)%>% data.frame(.)

library(ggplot2)
library(ggnewscale)


plot_function4 = function (variable1,variable2,time,dose) {
  files[,"dose_uM"] <- factor(files[,"dose_uM"])
  p = files %>% filter(., dose_uM %in% dose) %>% filter(., Times %in% time) %>% 
    ggplot(aes_string(x = variable1, y = variable2,color = "dose_uM"))+
    geom_path(size= 0.4)+
    # scale_color_manual(values = c("black","green","yellow","grey70","red","blue"))+
    # new_scale_color() +
    geom_point(aes(color = as.factor(Times)),size = 1)+ 
    #geom_label(label = files$Times) +  #+scale_color_gradient(low = "blue", high = "red")+
    labs(x = paste0(variable1), y = paste0(variable2), sec.x="Exposure (mg/kg)",
         title = paste0("NFT concentrations in Rabbit plasma following a single intravenous dose")) +
    
    theme(legend.title = element_blank()) + theme(legend.position="bottom") +
    theme(legend.text = element_text(size= 8,face="bold")) 
  
  p + geom_text(aes(label = Times), hjust = 1, vjust = 0, size = 2)
}
plot_function4(variable1 =proteins[12] ,variable2 = proteins[8], time = c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1, 0.2, 0.3,0.4,0.8,1.2),dose = Final_Rdose1[c(5,6)])
#1.1,1.2,1.3,1.4,1.5,1.6,.17,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,3,3.5,4
proteins
c(0,0.1,0.2,0.3,0.4)
,0.8, 1, 2,3