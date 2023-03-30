rm(list = ls())

library(tidyverse)
library(rstan)
library(StanHeaders)
library(bayesplot)
library(corrplot)
library(sensitivity)
library(pksensi)
library(foreach)     
library(doParallel) 
library(iterators)
library(parallel)

setwd("~/MCSim/mod")

# install.packages(c("iterators","parallel"))
theme_set(theme_light())
getwd()

Modeling_dir= "NFT_PBPK2"

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
    message(paste("Execute:", " ./mcsim_", model, " ", dir, "/", input, sep = ""))
    system(paste("./mcsim_", model, " ", dir, "/", input, sep = ""))
    df <- read.delim("simmc.out")
  } else {
    message(paste("Execute:", " ./mcsim_", model, " ", dir, "/", input, sep = ""))
    system(paste("./mcsim_", model, " ", dir, "/", input, sep = ""))
    #df <- read.delim("sim.out", skip = 1)
    df <- read.delim(paste0(ouput),skip = 1)
  }
  return(df)
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

report <- function() {
  cat("\n\n-----Report started line-----\n\n")
  cat(Sys.getenv("PATH"), "\n")
  print(Sys.which("gcc"))
  system('gcc -v')
}

readsims <- function(x, exp = 1) {
  ncols <- ncol(x)
  index <- which(x[, 1] == "Time")
  str <- ifelse(exp == 1, 1, index[exp - 1] + 1)
  end <- ifelse(exp == length(index) + 1, nrow(x), index[exp] - 2)
  X <- x[c(str:end), ]
  ncolX <- ncol(X)
  X <- as.data.frame(matrix(as.numeric(as.matrix(X)), ncol = ncolX))
  if (exp > 1)
    names(X) <-
    as.matrix(x[index[exp - 1], ])[1:ncols]
  else
    names(X) <- names(x)
  X <- X[, colSums(is.na(X)) != nrow(X)]
  return(X)
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




model <- "Agedynamics_NFThuman.model.R"
input <- "Agedynamic_montecarlo.in.R"

output <- "MonteHumanPBPK.out"

# Parmeters = read.csv(paste0(Modeling_dir,"/","NFTPBPK_parms.csv"))
# # Integrate (Lsodes, 1e-6, 1e-6, 1);
# # # OutputFile ("HumanPBPKforward.out");
# # 
# # MonteCarlo("MonteHumanPBPK.out", 5000, 6876875.9);
# 
# scaling = 0.25;
# 
# iterationnumber = 2000
# randomseed = 10101010
# 
# text1 = paste0('MonteCarlo','(', output, ',', iterationnumber,',',randomseed,")",";")
# 
# 
# ParmsFun = function (Para, value1, value2, value3, value4){
#   #return(paste0("Distrib","(" ,Para,",", "TruncNormal",",",value1,",",  value2,",", value3,",", value4,")",";"))
#   return(paste0("Distrib","(" ,Para,",", "LogNormal",",",value1,",",  value2,")",";"))
# }
# 
# 
# 
#   
# ParmsFun(Para = Parmeters$Parms,value1 = Parmeters$Value, value2 = 1.17)
# 
# file_name = "NFT_female_PK.csv"
# myfiles <- read.csv(paste0(getwd(), "/", Modeling_dir, "/", file_name))
# head(myfiles)
# #Timeid = paste0(unique(myfiles$Time, incomparables = FALSE),sep = ",",collapse = " ") %>% gsub(",$", "", .)  #gsub to remove the last comma
# Final_Rdose1 = c(50,100,200)
# file.remove(paste0(getwd(),"/",Modeling_dir, "/",input))
# exp_route = "oraldose"
# proteins = c("cplasma","P_excreted","B_excreted","Massbalance","Filteration","Reabsorption","Secretion","Total_renal_Cl","GFR_contribution","Secr_contribution","Reabs_contribution","Re_GFR_contribution","cgut", "cliver", "cfat","ckidney","cfilterate", "crestbody", "Agutlumen","Adelay","Aurine") 
# 
# 
# Times = paste(0, 24, 0.01,sep = ",") 
# 
# simulationprint = function (dose,biomarker1,biomarker2,biomarker3, biomarker4, biomarker5, biomarker6, biomarker7, biomarker8,biomarker9, biomarker10,biomarker11,biomarker12,biomarker13,biomarker14,biomarker15,biomarker16,biomarker17,biomarker18,biomarker19,biomarker20,biomarker21,Times){
#   
#   variable1 = ifelse(biomarker1 == F, "#", paste0( "PrintStep", "(", biomarker1,",",Times, ")",";"))
#   variable2 = ifelse(biomarker2 == F, "#", paste0( "PrintStep", "(", biomarker2,",",Times, ")",";"))
#   variable3 = ifelse(biomarker3 == F, "#", paste0( "PrintStep", "(", biomarker3,",",Times, ")",";"))
#   variable4 = ifelse(biomarker4 == F, "#", paste0( "PrintStep", "(", biomarker4,",",Times, ")",";"))
#   variable5 = ifelse(biomarker5 =="", "#", paste0( "PrintStep", "(", biomarker5,",",Times, ")",";"))
#   variable6 = ifelse(biomarker6 == "", "#", paste0( "PrintStep", "(", biomarker6,",",Times, ")",";"))
#   variable7 = ifelse(biomarker7 == "", "#", paste0( "PrintStep", "(", biomarker7,",",Times, ")",";"))
#   variable8 = ifelse(biomarker8 == "", "#", paste0( "PrintStep", "(", biomarker8,",",Times, ")",";"))
#   variable9 = ifelse(biomarker9 == "", "#", paste0( "PrintStep", "(", biomarker9,",",Times, ")",";"))
#   variable10 = ifelse(biomarker10 == "", "#", paste0( "PrintStep", "(", biomarker10,",",Times, ")",";"))
#   variable11 = ifelse(biomarker11 == "", "#", paste0( "PrintStep", "(", biomarker11,",",Times, ")",";"))
#   variable12 = ifelse(biomarker12 == "", "#", paste0( "PrintStep", "(", biomarker12,",",Times, ")",";"))
#   variable13 = ifelse(biomarker13 == "", "#", paste0( "PrintStep", "(", biomarker13,",",Times, ")",";"))
#   variable14 = ifelse(biomarker14 == "", "#", paste0( "PrintStep", "(", biomarker14,",",Times, ")",";"))
#   variable15 = ifelse(biomarker15 == "", "#", paste0( "PrintStep", "(", biomarker15,",",Times, ")",";"))
#   
#   variable16 = ifelse(biomarker16 == "", "#", paste0( "PrintStep", "(", biomarker16,",",Times, ")",";"))
#   variable17 = ifelse(biomarker17 == "", "#", paste0( "PrintStep", "(", biomarker17,",",Times, ")",";"))
#   variable18 = ifelse(biomarker18 == "", "#", paste0( "PrintStep", "(", biomarker18,",",Times, ")",";"))
#   variable19 = ifelse(biomarker19 == "", "#", paste0( "PrintStep", "(", biomarker19,",",Times, ")",";"))
#   variable20 = ifelse(biomarker20 == "", "#", paste0( "PrintStep", "(", biomarker20,",",Times, ")",";"))
#   variable21 = ifelse(biomarker21 == "", "#", paste0( "PrintStep", "(", biomarker21,",",Times, ")",";"))
#   
#   Exposure = ifelse(exp_route == "IV", paste0( "IVdosing", "=", dose, ";"),paste0("oraldose", "=", dose,";", "\n", "oraldose1", "=","PerDose", "(", dose, ",",24,",", 0,",", 0.001,")", ";"))
#   
#   
#   return(paste0("Simulation", "{" , "\n", Exposure ,"\n",
#                 variable1,"\n",variable2,"\n",variable3,"\n",variable4,"\n",variable5,"\n",variable6,"\n",variable7,
#                 "\n",variable8,
#                 "\n",variable9,
#                 "\n",variable10,
#                 "\n",variable11,
#                 "\n",variable12,
#                 "\n",variable13,
#                 "\n",variable14,
#                 "\n",variable15,
#                 "\n",variable16,
#                 "\n",variable17,
#                 "\n",variable18,
#                 "\n",variable19,
#                 "\n",variable20,
#                 "\n",variable21,
#                 "\n","}"))
#   
# }
# 
# 
# f = file(paste0(getwd(),"/",Modeling_dir, "/",input), open = 'a')
# 
# 
# write(text1,f)
# 
# write(paste0("Distrib","(" ,"age",",", "TruncNormal",",",35,",",  10,",", 25,",", 40,")",";"),f)
# 
# write(ParmsFun(Para = Parmeters$Parms,value1 = Parmeters$Value, value2 = 1.17),f)
# 
# 
# for (i in 1:length(Final_Rdose1)){
#   
#   write(simulationprint(dose = Final_Rdose1[i],biomarker1 = proteins[1],biomarker2 = proteins[2],biomarker3 = proteins[3],
#                         biomarker4 = proteins[4],biomarker5 = proteins[5], biomarker6 = proteins[6],biomarker7 = proteins[7],
#                         biomarker8 = proteins[8],biomarker9 = proteins[9],
#                         biomarker10 = proteins[10],biomarker11 = proteins[11],
#                         biomarker12 = proteins[12],biomarker13 = proteins[13],
#                         biomarker14 = proteins[14],biomarker15 = proteins[15],
#                         biomarker16 = proteins[16],biomarker17 = proteins[17],
#                         biomarker18 = proteins[18],biomarker19 = proteins[19],biomarker20 = proteins[20],
#                         biomarker21 = proteins[21],
#                         Times), f)
#   
# }
# 
# 
# 
# # write("End.", f)
# close(f)


makemcsim(model = model, dir = Modeling_dir)

out <-mcsim(
  model = model,
  input = input,
  dir = Modeling_dir,
  parallel = F)

out = read.delim(output, header = TRUE, sep = "")
paramsnumber = 20
vars <- names(out)
params = data.frame((vars))
colnames(params) = "parametersname"


proteins = c("cplasma","P_excreted","B_excreted","Massbalance","Filteration","Reabsorption","Secretion","Total_renal_Cl","GFR_contribution","Secr_contribution","Reabs_contribution","Re_GFR_contribution","cgut", "cliver", "cfat","ckidney","cfilterate", "crestbody", "Agutlumen","Adelay","Aurine") 


vars <- names(out)

Simulationfile = data.frame()


for (i in 1:length(proteins)){
  finder =  grep(pattern = paste0("^",proteins[i],collapse = ""),vars)
  file <- apply(out[finder], 2, quantile,  c(0.5, 0.025, 0.975)) %>% t()
  x = data.frame(file)
  X_mean <- apply(out[finder],2, mean)   # can apply for all the doses
  x$mean = X_mean
  Simulationfile <- rbind(Simulationfile,x)
}
nrow(Simulationfile)

Final_Rdose1 = c(50,100,200)
df <- data.frame(Final_Rdose1)
biomarker = length(proteins)
Times = seq(0,24, by = 0.1) #experimental data time point
simul = 3 #if we include dose zero to check the steady state
#rename in captial as experimental data 
Simulationfile$Simulation = (rep(paste0(rep("Simu",simul),1:simul),each = length(Times)*biomarker))  # as we have 3 variables
Simulationfile$dose_uM = as.numeric(rep(paste0(rep(Final_Rdose1)),each = length(Times)))
Simulationfile$Chemical = rep(paste0("NFT"),each = nrow(Simulationfile))
Simulationfile$Times <- Times
Simulationfile$variables = (rep(paste0(rep(proteins)),each = length(Times)*simul))
colnames(Simulationfile) <- c("median", "LCL", "UCL","mean","Simulations","dose_uM","Chemicals","Times","variables")

#write.table(Simulationfile, file = paste0(getwd(),"/",Modeling_dir, "/", "HumanMontecarloSimulation.out"), row.names = F, sep = "\t")

write.table(Simulationfile, file = paste0(getwd(),"/",Modeling_dir, "/", "Posterior_chains","/", "HumanMontecarloSimulation.out"), row.names = F, sep = "\t")



file_name = "NFT_female_PK.csv"
myfiles <- read.csv(paste0(getwd(), "/", Modeling_dir, "/", file_name))
Expdata = myfiles


plot_function = function (variable1,variable2) {
  Simulation = Simulationfile %>% filter(.,variables == variable1)
  #Simulation = Simulationfile %>% filter(.,variables == variable & Simulations != "Simu1")
  exp_data = Expdata %>% filter(Variable == variable2) %>% filter(., dose_uM %in% c(Final_Rdose1))
  yaxis = ifelse(variable1 == "P_excreted", "Fraction of Dose excreted (%)","Concentration(mg/L)") 
  ggplot() +
    geom_line(data =Simulation,aes(x = Times, y = mean, color="Simulated_mean"),size=1) +
    geom_line(data =Simulation,aes(x = Times, y = LCL, color="Simulated_p2.5"),size=1) +
    geom_line(data =Simulation,aes(x = Times, y = UCL, color="Simulated_p95"),size=1) +
    geom_point(data =exp_data,aes(x = Time, y = cplasma, color="Exp_mean"), size = 1.2) +
    geom_errorbar(data =exp_data,aes(x = Time,ymin=cplasma-cplasmaSD, ymax= cplasma +cplasmaSD,color = "Experiments(Sd)"), size=0.5,
                  width=.25)+
    facet_wrap(~dose_uM) +
    labs(x = "Time (hr)", y = "Plasma conc.(mg/L)", sec.x="Exposure (mg)",
         title = paste0("NFT concentrations in Human plasma following a single oral dose")) +
    #scale_y_continuous(trans = 'log10')+
    #scale_y_continuous(limits=c(0,3)) +
    scale_color_manual(name = "type",
                       breaks = c("Simulated_mean","Simulated_p2.5","Simulated_p95", "Exp_mean","Experiments(Sd)"),
                       values = c("Simulated_mean" = "blue","Simulated_p2.5" = "green","Simulated_p95" = "red","Exp_mean" = "black","Experiments(Sd)" = "black"))+

    # scale_color_manual(name = "type",
    #                    breaks = c("Simulated_mean","Exp_mean(Sd)"),
    #                    values = c("Simulated_mean" = "blue","Exp_mean(Sd)" = "black"))+
    theme(legend.title = element_blank()) + theme(legend.position="bottom") +
    theme(legend.text = element_text(size= 8,face="bold")) 
}

plot_function(variable1 = "cplasma",variable2 = "cplasma")

plot_function(variable1 = "P_excreted",variable2 = "P_excreted")


colnames(Expdata) = c("Time", "cplasma", "cplasmaSD","dose_uM","treatment", "variables")

log10_trans <- function(x) log10(x)
identity_trans <- function(x) x


plot_function = function (variable1) {
  Simulation = Simulationfile %>% filter(.,variables %in%  variable1)
  Simulation$variables <- factor(Simulation$variables, levels = c("cplasma", "cliver", "ckidney", "Aurine","B_excreted", "P_excreted"))
  
  exp_data = Expdata %>% filter(.,variables %in%  variable1)%>% filter(., dose_uM %in% c(Final_Rdose1))
  yaxis = ifelse(variable1 == "P_excreted", "Fraction of Dose excreted (%)","Concentration(mg/L)") 
  ggplot() +
   geom_line(data =Simulation,aes(x = Times, y = mean, color="Simulated_mean"),size=1) +
    geom_line(data =Simulation,aes(x = Times, y = LCL, color="Simulated_p2.5"),size=1) +
    geom_line(data =Simulation,aes(x = Times, y = UCL, color="Simulated_p95"),size=1) +
    facet_grid(variables~dose_uM, scales = "free") +
    labs(x = "Time (hr)", y = "Plasma conc.(mg/L)", sec.x="Exposure (mg)",
         title = paste0("NFT concentrations in Human plasma following a single oral dose")) +
   scale_y_continuous(trans = trans_func) +
  
    scale_color_manual(name = "type",
                       breaks = c("Simulated_mean","Simulated_p2.5","Simulated_p95", "Exp_mean","Experiments(Sd)"),
                       values = c("Simulated_mean" = "blue","Simulated_p2.5" = "green","Simulated_p95" = "red","Exp_mean" = "black","Experiments(Sd)" = "black"))+

  # log10_trans <- function(x) log10(x)
  # identity_trans <- function(x) x
  # if (!(variable1 %in% c("Aurine","B_excreted","P_excreted"))) {
  #   trans_func <- log10_trans
  # } else {
  #   trans_func <- identity_trans
  # }
  # 
  
    theme(legend.title = element_blank()) + theme(legend.position="bottom") +
    theme(legend.text = element_text(size= 8,face="bold")) 
}

plot_function(variable1 = c("cplasma", "cliver", "ckidney", "Aurine","B_excreted", "P_excreted"))
###############################
# TypeSim = c("Normal_GFR", "infants", "Kid", "Old", "Mild_GFR", "Moderate_GFR", "Severe_GFR")
# 
# 
# Final_Rdose1 = c(50,100)
# df <- data.frame(TypeSim)
# df1 = df
# df <- data.frame(TypeSim)
# df2 = data.frame(Final_Rdose1)
# colnames(df2) = "dose_uM"
# z = list()
# experiment = length(TypeSim)
# for (i in 1:length(TypeSim)){
#   z[[i]] = readsims(out, exp = i)
# }
# 
# z = list()
# experiment = 14
# for (i in 1:14){
#   z[[i]] = readsims(out, exp = i)
# }
# v1<- do.call(rbind,z)
# str(v1)
# biomarkers = colnames(v1)[-1]
# Times = seq(0,24, by = 0.1) #experimental data time point
# colnames(df1)= c("TypeSim")
# v1$Simulation = (rep(paste0(rep("Simu",length(TypeSim)),1:length(TypeSim)),each = length(Times)))  # as we have 3 variables
# v1$Chemical = rep(paste0("NFT"),each = nrow(v1))
# v1$TypeSim = rep(df1$TypeSim,each = length(Times))
# v1$dose_uM = rep(df2$dose_uM,each = length(Times)*length(TypeSim))
# 
# v2 <- v1 %>%
#   mutate(dose_GFR = paste0(TypeSim, "_", dose_uM))
# 
# # modified_files <- files %>% 
# #   mutate(modified_name = if_else(dose_GFR == "Filteration", "GFR", 
# #                                  if_else(dose_GFR == "Reabsorption", "Tubular Reabsorption", 
# #                                          if_else(name == "Secretion", "Tubular Secretion", name))))
# 
# v3 = v2 %>%  gather(variable, value, cplasma:Simulation)
# 
# v3$dose_GFR <- as.factor(v3$dose_GFR)
# str(v3)
# 
# # mutate(grouping_column = as.factor(grouping_column)) %>%
# #   group_by(grouping_column)
# 
# # Simulation = data_long %>% filter(variable == "cplasma") %>%  filter(TypeSim %in% c("Normal_GFR", "Mild_GFR", "Moderate_GFR", "Severe_GFR"))
# # Simulation$TypeSim <- factor(Simulation$TypeSim, levels = c("Normal_GFR", "Mild_GFR", "Moderate_GFR", "Severe_GFR"))
# # ggplot(Simulation, aes(x = Time, y = value, color = dose_GFR)) +
# #   geom_line() +
# #   xlab("X Variable") +
# #   ylab("Y Variable") +
# #   ggtitle("Multiple Line Plot")
# 
# 
# plot_function1 = function (variable1) {
#   Simulation = v3 %>% filter(variable %in% variable1) %>%  filter(TypeSim %in% c("Normal_GFR", "Mild_GFR", "Moderate_GFR", "Severe_GFR"))
#   Simulation$TypeSim <- factor(Simulation$TypeSim, levels = c("Normal_GFR", "Mild_GFR", "Moderate_GFR", "Severe_GFR"))
#   Simulation$value <- as.numeric(Simulation$value)
#   Simulation$variable <- factor(Simulation$variable, levels = c("cplasma", "cliver", "ckidney"))
#   ggplot() +
#     geom_line(data =Simulation,aes(x = Time, y = value, group= dose_GFR,color = as.factor(dose_uM)),size=1) +
#     #facet_wrap(~TypeSim) +
#     facet_grid(variable~TypeSim) +
#     #scale_y_continuous(breaks = 5) +
#     labs(x = "Time (hr)", y = "Plasma conc.(mg/L)", sec.x="Conditions",
#          title = paste0("NFT in human plasma following a repeated OD")) +
#     
#     theme(legend.title = element_blank()) + theme(legend.position="bottom") +
#     theme(legend.text = element_text(size= 8,face="bold")) 
# }
# 
# plot_function1("cplasma","cliver")
# 
# 
# ggplot(data = Simulation, aes(x = Time, y = value, group = group_var, color = group_var)) +
#   geom_line(size = 1) +
#   scale_color_manual(values = c("group1" = "blue", "group2" = "red")) +
#   labs(x = "Time (hr)", y = "Plasma conc.(mg/L)", title = paste0("NFT concentrations in Human plasma following a repeated oral dose of 50 mg")) +
#   theme(legend.title = element_blank()) + theme(legend.position = "bottom") +
#   theme(legend.text = element_text(size = 8, face = "bold"))
# 
# 
# plot_function1 = function (variable1) {
#   #Simulation = v2 %>% filter(variable == variable1) %>%  filter(TypeSim %in% c("Normal_GFR", "Mild_GFR", "Moderate_GFR", "Severe_GFR"))
#   Simulation = v2 %>% filter(variable %in% variable1)%>%  filter(TypeSim %in% c("Normal_GFR", "Mild_GFR", "Moderate_GFR", "Severe_GFR"))
#   Simulation$TypeSim <- factor(Simulation$TypeSim, levels = c("Normal_GFR", "Mild_GFR", "Moderate_GFR", "Severe_GFR"))
#   Simulation$variable <- factor(Simulation$variable, levels = c("cplasma", "cliver", "ckidney"))
#   ggplot() +
#     geom_line(data =Simulation,aes(x = Time, y = value, color="Simulated_mean"),size=1) +
#     
#     facet_grid(variable~TypeSim) +
#     labs(x = "Time (hr)", y = "NFT conc.(mg/L)", sec.x="Conditions",
#          title = paste0("NFT in human plasma following a OD of 100 mg qid")) +
#     #scale_y_continuous(trans = 'log10')+
#     # scale_y_continuous(limits=c(0,3)) +
#     scale_color_manual(name = "type",
#                        breaks = c("Simulated_mean"),
#                        values = c("Simulated_mean" = "blue"))+
#     theme(legend.title = element_blank()) + theme(legend.position="bottom") +
#     theme(legend.text = element_text(size= 8,face="bold")) 
# }
# 
# plot_function1(variable1 = c("cplasma","cliver"))
# 
# 


#Final_Rdose1 = c("1.infants", "2.kids", "3.Teenagers", "4.Adult","5.old", "6. Advanced old age")
Final_Rdose1 = c("1.NormalGFR", "2.mildGFR", "3.moderateGFR", "4.SeverGFR")
df <- data.frame(Final_Rdose1)
z = list()
experiment = length(Final_Rdose1)
for (i in 1:length(Final_Rdose1)){
  z[[i]] = readsims(out, exp = i)
}
v1<- do.call(rbind,z)


df1 = df
Times = seq(0,24, by = 1) #experimental data time point
colnames(df1)= c("dose_uM")
v1$Simulation = (rep(paste0(rep("Simu",length(Final_Rdose1)),1:length(Final_Rdose1)),each = length(Times)))  # as we have 3 variables
v1$Chemical = rep(paste0("NFT"),each = nrow(v1))
v1$dose_uM = rep(df1$dose_uM,each = length(Times))

v2 = v1 %>%  gather(variable, value, cplasma:Aurine)
# file_name = "NFT_female_PK.csv"
# myfiles <- read.csv(paste0(getwd(), "/", Modeling_dir, "/", file_name))
# head(myfiles)

plot_function1 = function (variable1) {
  Simulation = v2 %>% filter(variable == variable1)
  ggplot() +
    geom_line(data =Simulation,aes(x = Time, y = value, color="Simulated_mean"),size=1) +
    
    facet_wrap(~dose_uM) +
    labs(x = "Time (hr)", y = "Plasma conc.(mg/L)", sec.x="Conditions",
         title = paste0("NFT concentrations in Human plasma following a repeated oral dose")) +
    #scale_y_continuous(trans = 'log10')+
    # scale_y_continuous(limits=c(0,3)) +
    scale_color_manual(name = "type",
                       breaks = c("Simulated_mean"),
                       values = c("Simulated_mean" = "blue"))+
    theme(legend.title = element_blank()) + theme(legend.position="bottom") +
    theme(legend.text = element_text(size= 8,face="bold")) 
}


plot_function2 = function (variable1) {
  Simulation = v2 %>% filter(variable == variable1)
  ggplot() +
    geom_line(data =Simulation,aes(x = Time, y = value, color="Simulated_mean"),size=1) +
    
    facet_wrap(~dose_uM) +
    labs(x = "Time (hr)", y = "Urine amount (mg)", sec.x="Conditions",
         title = paste0("Cumulative NFT amount in Urine following a repeated oral dose")) +
    #scale_y_continuous(trans = 'log10')+
    # scale_y_continuous(limits=c(0,3)) +
    scale_color_manual(name = "type",
                       breaks = c("Simulated_mean"),
                       values = c("Simulated_mean" = "blue"))+
    theme(legend.title = element_blank()) + theme(legend.position="bottom") +
    theme(legend.text = element_text(size= 8,face="bold")) 
}

plot_function1(variable1 = "cplasma")

plot_function2(variable1 = "Aurine")

plot_function3 = function (variable1) {
  Simulation = v2 %>% filter(variable == variable1)
  ggplot() +
    geom_line(data =Simulation,aes(x = Time, y = value, color="Simulated_mean"),size=1) +
    
    facet_wrap(~dose_uM) +
    labs(x = "Time (hr)", y = paste0(variable1, "conc.(mg/L)"), sec.x="Conditions",
         title = paste0("NFT concentrations in Kidney following a repeated oral dose")) +
    #scale_y_continuous(trans = 'log10')+
    # scale_y_continuous(limits=c(0,3)) +
    scale_color_manual(name = "type",
                       breaks = c("Simulated_mean"),
                       values = c("Simulated_mean" = "blue"))+
    theme(legend.title = element_blank()) + theme(legend.position="bottom") +
    theme(legend.text = element_text(size= 8,face="bold")) 
}

plot_function3(variable1 = "ckidney")
plot_function3(variable1 = "cliver")

plot_function4 = function (variable1) {
  Simulation = v2 %>% filter(variable == variable1)
  ggplot() +
    geom_line(data =Simulation,aes(x = Time, y = value, color="Simulated_mean"),size=1) +
    
    facet_wrap(~dose_uM) +
    labs(x = "Time (hr)", y = "Tubular conc.(mg/L)", sec.x="Conditions",
         title = paste0("NFT concentrations in Kidney tubules following a repeated oral dose")) +
    #scale_y_continuous(trans = 'log10')+
    # scale_y_continuous(limits=c(0,3)) +
    scale_color_manual(name = "type",
                       breaks = c("Simulated_mean"),
                       values = c("Simulated_mean" = "blue"))+
    theme(legend.title = element_blank()) + theme(legend.position="bottom") +
    theme(legend.text = element_text(size= 8,face="bold")) 
}

plot_function4(variable1 = "cfilterate")


plot_function5 = function (variable1) {
  Simulation = v2 %>% filter(variable == variable1)
  ggplot() +
    geom_line(data =Simulation,aes(x = Time, y = value, color="Simulated_mean"),size=1) +
    
    facet_wrap(~dose_uM) +
    labs(x = "Time (hr)", y = "Liver conc.(mg/L)", sec.x="Conditions",
         title = paste0("NFT concentrations in liver following a repeated oral dose")) +
    #scale_y_continuous(trans = 'log10')+
    # scale_y_continuous(limits=c(0,3)) +
    scale_color_manual(name = "type",
                       breaks = c("Simulated_mean"),
                       values = c("Simulated_mean" = "blue"))+
    theme(legend.title = element_blank()) + theme(legend.position="bottom") +
    theme(legend.text = element_text(size= 8,face="bold")) 
}

plot_function5(variable1 = "cliver")



plot_function2 = function (variable1,variable2) {
  Simulation = v2 %>% filter(variable == variable1)
  exp_data = myfiles %>% filter(Variable == variable2)
  ggplot() +
    geom_line(data =Simulation,aes(x = Time, y = value, color="Simulated_mean"),size=1) +
    #geom_point(data =exp_data,aes(x = Time, y = cplasma, color="Exp_mean"), size = 1.2) +
    
    facet_wrap(~dose_uM) +
    labs(x = "Time (hr)", y = "% of dose excreted", sec.x="Exposure (mg/kg)",
         title = paste0("Dose normalized cumulative urine excretion percentage in Human following a single oral dose")) +
    #scale_y_continuous(trans = 'log10')+
    # scale_y_continuous(limits=c(0,3)) +
    scale_color_manual(name = "type",
                       breaks = c("Simulated_mean", "Exp_mean"),
                       values = c("Simulated_mean" = "blue","Exp_mean" = "blue"))+
    theme(legend.title = element_blank()) + theme(legend.position="bottom") +
    theme(legend.text = element_text(size= 8,face="bold")) 
}

plot_function1(variable1 = "cplasma",variable2 = "cplasma")
plot_function2(variable1 = "P_excreted",variable2 = "P_excreted")

plot_function(variable1 = "cplasma",variable2 = "cplasma")

plot_function(variable1 = "P_excreted",variable2 = "P_excreted")

plot_function1(variable1 = "ckidney",variable2 = "cplasma")

plot_function(variable1 = "cfilterate",variable2 = "cplasma")

plot_function(variable1 = "Aurine",variable2 = "cplasma")

plot_function(variable1 = "cliver",variable2 = "cplasma")
