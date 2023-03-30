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
input <- "GFRdiseasedAgedynamic_montecarlo.in.R"


output <- "GFRMonteHumanPBPK.out"



makemcsim(model = model, dir = Modeling_dir)

out <-mcsim(
  model = model,
  input = input,
  dir = Modeling_dir,
  parallel = F)

out = read.delim(output, header = TRUE, sep = "")
# paramsnumber = 20
vars <- names(out)
params = data.frame((vars))
colnames(params) = "parametersname"


# Verify the result
library(magrittr)
file = params %>% separate(col = parametersname, c("Parms", "Subj", "levels"),sep = "_|\\.",fill = "right") %>% data.frame() 
       
variable = data.frame(unique(file$Parms)) %>% set_colnames(c("Parms")) %>%  mutate(modified_parms = if_else(Parms == "P", "P_excreted",Parms))

proteins = variable$modified_parms[-(1:20)]
Simulationfile1 = data.frame()
for (i in 1:length(proteins)){
  finder =  grep(pattern = paste0("^",proteins[i],collapse = ""),vars)
  file <- apply(out[finder], 2, quantile,  c(0.5, 0.025, 0.975)) %>% t()
  x = data.frame(file)
  X_mean <- apply(out[finder],2, mean)   # can apply for all the doses
  x$mean = X_mean
  Simulationfile1 <- rbind(Simulationfile1,x)
}
nrow(Simulationfile1)


# first 20 rows are paramters

Simulationfile = Simulationfile1
Final_Rdose1 = c(50, 100)  #now i am simulating only 50 mg per two dose
Final_Rdose1 = c(50)  #now i am simulating only 50 mg per two dose
Conditions = c("1.NormalGFR", "2.mildGFR", "3.moderateGFR", "4.SevereGFR")
df1 <- data.frame(Conditions)
df2 <- data.frame(Final_Rdose1)

Times = seq(0,120, by = 0.1) #experimental data time point
biomarker = length(proteins)
#Times = seq(0,24, by = 0.1) #experimental data time point
simul = length(Final_Rdose1)*length(Conditions) #if we include dose zero to check the steady state
#rename in captial as experimental data 
Simulationfile$Simulation = (rep(paste0(rep("Simu",simul),1:simul),each = length(Times)*biomarker))  # as we have 3 variables
Simulationfile$dose_uM = as.numeric(rep(paste0(rep(Final_Rdose1)),each = length(Times)*length(Conditions)))
Simulationfile$Conditions = (rep(paste0(rep(Conditions)),each = length(Times)))
Simulationfile$Chemical = rep(paste0("NFT"),each = nrow(Simulationfile))
Simulationfile$Times <- Times
Simulationfile$variables = (rep(paste0(rep(proteins)),each = length(Times)*simul))
colnames(Simulationfile) <- c("median", "LCL", "UCL","mean","Simulations","dose_uM","Conditions","Chemicals","Times","variables")

Simulationfile = Simulationfile %>%  mutate(.,Dose_condition = paste(dose_uM ,Conditions, sep = "_"))


#write.table(Simulationfile, file = paste0(getwd(),"/",Modeling_dir, "/", "GFRcompromisedSimulation.out"), row.names = F, sep = "\t")

write.table(Simulationfile, file = paste0(getwd(),"/",Modeling_dir, "/", "Posterior_chains","/", "GFRcompromisedSimulation.out"), row.names = F, sep = "\t")




##################################
# For single dosing scenario

input <- "PexcretedGFRdiseasedAgedynamic_montecarlo.in.R" #This is to extract for single dose
output <- "SingledoseGFRMonteHumanPBPK.out"


out <-mcsim(
  model = model,
  input = input,
  dir = Modeling_dir,
  parallel = F)

out = read.delim(output, header = TRUE, sep = "")
# paramsnumber = 20
vars <- names(out)
params = data.frame((vars))
colnames(params) = "parametersname"


# Verify the result
library(magrittr)
file = params %>% separate(col = parametersname, c("Parms", "Subj", "levels"),sep = "_|\\.",fill = "right") %>% data.frame() 

variable = data.frame(unique(file$Parms)) %>% set_colnames(c("Parms")) %>%  mutate(modified_parms = if_else(Parms == "P", "P_excreted",Parms))

proteins = variable$modified_parms[-(1:20)]
Simulationfile1 = data.frame()
for (i in 1:length(proteins)){
  finder =  grep(pattern = paste0("^",proteins[i],collapse = ""),vars)
  file <- apply(out[finder], 2, quantile,  c(0.5, 0.025, 0.975)) %>% t()
  x = data.frame(file)
  X_mean <- apply(out[finder],2, mean)   # can apply for all the doses
  x$mean = X_mean
  Simulationfile1 <- rbind(Simulationfile1,x)
}
nrow(Simulationfile1)


# first 20 rows are paramters

Simulationfile = Simulationfile1
Final_Rdose1 = c(200)
Conditions = c("1.NormalGFR", "2.mildGFR", "3.moderateGFR", "4.SevereGFR")
df1 <- data.frame(Conditions)
df2 <- data.frame(Final_Rdose1)

#Times = seq(0,120, by = 0.1) #experimental data time point
biomarker = length(proteins)
Times = seq(0,24, by = 0.1) #experimental data time point
simul = length(Final_Rdose1)*length(Conditions) #if we include dose zero to check the steady state
#rename in captial as experimental data 
Simulationfile$Simulation = (rep(paste0(rep("Simu",simul),1:simul),each = length(Times)*biomarker))  # as we have 3 variables
Simulationfile$dose_uM = as.numeric(rep(paste0(rep(Final_Rdose1)),each = length(Times)*length(Conditions)))
Simulationfile$Conditions = (rep(paste0(rep(Conditions)),each = length(Times)))
Simulationfile$Chemical = rep(paste0("NFT"),each = nrow(Simulationfile))
Simulationfile$Times <- Times
Simulationfile$variables = (rep(paste0(rep(proteins)),each = length(Times)*simul))
colnames(Simulationfile) <- c("median", "LCL", "UCL","mean","Simulations","dose_uM","Conditions","Chemicals","Times","variables")

Simulationfile = Simulationfile %>%  mutate(.,Dose_condition = paste(dose_uM ,Conditions, sep = "_"))


#write.table(Simulationfile, file = paste0(getwd(),"/",Modeling_dir, "/", "GFRcompromisedSimulation.out"), row.names = F, sep = "\t")

write.table(Simulationfile, file = paste0(getwd(),"/",Modeling_dir, "/", "Posterior_chains","/", "200mg_SingledoseGFRcompromisedSimulation.out"), row.names = F, sep = "\t")




plot_function = function (variable1,dose) {
  Simulation = Simulationfile %>% filter(.,variables %in% variable1) %>% filter(., dose_uM %in% dose)
  yaxis = ifelse(variable1 == "P_excreted", "Fraction of Dose excreted (%)","Concentration(mg/L)") 
  ggplot() +
    geom_line(data =Simulation,aes(x = Times, y = mean, group = Dose_condition, color="Simulated_mean"),size=1) +
    geom_line(data =Simulation,aes(x = Times, y = LCL,  group = Dose_condition, color="Simulated_p2.5"),size=1) +
    geom_line(data =Simulation,aes(x = Times, y = UCL,  group = Dose_condition,color="Simulated_p95"),size=1) +
    facet_grid(variables~Conditions) +
    labs(x = "Time (hr)", y = "Plasma conc.(mg/L)", sec.x="Conditions",
         title = paste0("NFT concentrations in Human plasma following a repeated oral dose")) +
    #scale_y_continuous(trans = 'log10')+
    #scale_y_continuous(limits=c(0,3)) +
    scale_color_manual(name = "type",
                       breaks = c("Simulated_mean","Simulated_p2.5","Simulated_p95"),
                       values = c("Simulated_mean" = "blue","Simulated_p2.5" = "green","Simulated_p95" = "red"))+
    
    theme(legend.title = element_blank()) + theme(legend.position="bottom") +
    theme(legend.text = element_text(size= 8,face="bold")) 
}

plot_function(variable1 = c("cplasma", "cliver"), dose = 100)



#####################################################################
# data analysis
#####################################################################

Simulationfile =read.delim(paste0(getwd(),"/",Modeling_dir, "/", "GFRcompromisedSimulation.out"), header = TRUE, sep = "")


proteins = unique(Simulationfile$variables)
Final_Rdose1 = unique(Simulationfile$dose_uM)
Conditions = unique(Simulationfile$Conditions)
glimpse(Simulationfile)

proteins = c("cplasma", "cliver")

df_grouped = Simulationfile %>% filter(., variables %in% proteins) %>% 
  select(variables, Conditions,dose_uM, Times, median, UCL, LCL) %>%
  pivot_longer(-c(variables, Conditions,dose_uM, Times), names_to="value_type", values_to="value") %>%
  group_by(variables, Conditions,dose_uM, time=cut(Times, breaks=c(0, 24, 48, 72, 96, 120), labels=c(24, 48, 72, 96, 120)))

df_grouped_max = df_grouped %>% 
  filter(Times != 0) %>% 
  group_by(Conditions, variables, dose_uM, value_type, time) %>% 
  slice(which.max(value)) %>% 
  bind_rows() %>% rename(Cmax = value)

df_grouped_AUC = df_grouped %>% filter(Times != 0) %>%
  group_by(Conditions, variables, dose_uM, value_type,time) %>% 
  do(data.frame(AUC = auc(.$Times, .$value, type = "spline"))) %>% 
  bind_rows() 



df_grouped_max_time = df_grouped %>% 
  filter(Times != 0) %>% 
  group_by(Conditions, variables, dose_uM, value_type,time) %>% 
  slice(which.max(value)) %>% 
  slice(1) %>%
  bind_rows()%>% 
  rename(Tmax = Times)



df_grouped_trough = df_grouped %>% 
  filter(Times != 0) %>% 
  # filter(Times > df_grouped_max_time$Tmax[df_grouped_max_time$time.tmax == tp]) %>% 
  filter(Times > df_grouped_max_time$Tmax[df_grouped_max_time$Conditions %in% Conditions &
                                            df_grouped_max_time$variables %in% variables &
                                            df_grouped_max_time$dose_uM %in% dose_uM &
                                            df_grouped_max_time$value_type %in% value_type &
                                            df_grouped_max_time$time %in% time]) %>% 
  group_by(Conditions, variables, dose_uM, value_type,time) %>% 
  slice(which.min(value)) %>% 
  slice(1) %>%
  bind_rows()%>%
  rename(ctrough = value) %>% 
  rename(troughtime= Times)


df_grouped_AUC = df_grouped_AUC %>% 
  rename(variables.AUC = variables,
         Conditions.AUC = Conditions,
         dose_uM.AUC = dose_uM,
         value_type.AUC = value_type,
         time.AUC = time)

df_grouped_max_time = df_grouped_max_time %>% 
  rename(variables.tmax = variables,
         Conditions.tmax = Conditions,
         dose_uM.tmax = dose_uM,
         value_type.tmax = value_type,
         time.tmax = time)

df_grouped_trough = df_grouped_trough %>% 
  rename(variables.trough = variables,
         Conditions.trough = Conditions,
         dose_uM.trough = dose_uM,
         value_type.trough = value_type,
         time.trough = time)


df_combined = bind_cols(df_grouped_max, df_grouped_AUC,df_grouped_max_time,df_grouped_trough)
df2 = df_combined %>%  select(variables, Conditions,  dose_uM,value_type,Cmax, time,AUC,ctrough)
  

df_summarized1 <-df2 %>% 
  group_by(time, Conditions, variables,dose_uM) %>% 
  summarise(LCL_Cmax = min(Cmax[value_type == "LCL"]),
            median_Cmax = median(Cmax[value_type == "median"]),
            UCL_Cmax = max(Cmax[value_type == "UCL"]),
            LCL_AUC = min(AUC[value_type == "LCL"]),
            median_AUC = median(AUC[value_type == "median"]),
            UCL_AUC = max(AUC[value_type == "UCL"]),
            LCL_ctrough = min(ctrough[value_type == "LCL"]),
            median_ctrough = median(ctrough[value_type == "median"]),
            UCL_ctrough = max(ctrough[value_type == "UCL"]))



df_summarized = df_summarized1

#df_summarized = df_summarized1 %>%  filter(variables == "cliver" & time %in% c(24, 120))

df_summarized = df_summarized1 %>%  filter(variables == "cliver"& dose_uM == 50 & time %in% c(24, 120))

library(gridExtra)


colors <- c("median" = "green", "LCL" = "blue", "UCL" = "red")

p1 <- ggplot(df_summarized, aes(x = factor(time), y = median_Cmax, ymin = LCL_Cmax, ymax = UCL_Cmax)) + 
  geom_point(aes(color = "median"), size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = LCL_Cmax, color = "LCL"), shape = 17, size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = UCL_Cmax, color = "UCL"), shape = 14, size = 4) +
  scale_color_manual(values = colors) +
  guides(color = guide_legend(title = "Value Type")) +
  facet_grid(variables~Conditions) +
  ylim(0, 2)+
  xlab("Time") + ylab("Cmax") +
  theme_classic()

p2 <- ggplot(df_summarized, aes(x = factor(time), y = median_AUC, ymin = LCL_AUC, ymax = UCL_AUC)) + 
  geom_point(aes(color = "median"), size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = LCL_AUC, color = "LCL"), shape = 17, size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = UCL_AUC, color = "UCL"), shape = 14, size = 4) +
  scale_color_manual(values = colors) + 
  guides(color = guide_legend(title = "Value Type")) +
  facet_grid(variables~Conditions) +
  ylim(0, 40)+
  xlab("Time") + ylab("AUC") +
  theme_classic()

p3 <- ggplot(df_summarized, aes(x = factor(time), y = median_ctrough, ymin = LCL_ctrough, ymax = UCL_ctrough)) + 
  geom_point(aes(color = "median"), size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = LCL_ctrough, color = "LCL"), shape = 17, size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = UCL_ctrough, color = "UCL"), shape = 14, size = 4) +
  scale_color_manual(values = colors) + 
  guides(color = guide_legend(title = "Simulation type")) +
  facet_grid(variables~Conditions) +
  ylim(0, 1.2)+
  xlab("Time") + ylab("Ctrough") +
  theme_classic()


grid.arrange(p1+theme(legend.position='hidden'), p2+theme(legend.position='hidden'),
             p3+theme(legend.position='right'))


# P2 = plot_function1(variable = c("cplasma", "P_excreted","B_excreted"),mainplot ="ModelV5.NFT kinetics in rats following a single oral dose")
# 
# tiff(paste0(directory_i,"/", "OralRatAnnexFigure4",".png"), units="in",width=8, height=9, res=300, compression = 'lzw') 
# print(P2)
# dev.off()

colors <- c("median" = "green", "LCL" = "blue", "UCL" = "red")

df_summarized = df_summarized1 %>%  filter(variables == "cliver"& dose_uM == 50 & time %in% c(24,48, 120))

p1 <- ggplot(df_summarized, aes(x = factor(time), y = median_Cmax, ymin = LCL_Cmax, ymax = UCL_Cmax)) + 
  geom_point(aes(color = "median"), size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = LCL_Cmax, color = "LCL"), shape = 17, size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = UCL_Cmax, color = "UCL"), shape = 14, size = 4) +
  scale_color_manual(values = colors) +
  guides(color = guide_legend(title = "Value Type")) +
  facet_wrap(~Conditions, ncol = 4) +
  ylim(0, 1.5) +
  theme_classic()+
  xlab("") + ylab("Cmax") +
  theme(axis.line.x = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
 

p2 <- ggplot(df_summarized, aes(x = factor(time), y = median_AUC, ymin = LCL_AUC, ymax = UCL_AUC)) + 
  geom_point(aes(color = "median"), size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = LCL_AUC, color = "LCL"), shape = 17, size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = UCL_AUC, color = "UCL"), shape = 14, size = 4) +
  scale_color_manual(values = colors) + 
  guides(color = guide_legend(title = "Value Type")) +
  facet_wrap(~Conditions, ncol = 4) +
  ylim(0, 30)+
  xlab("") + ylab("AUC") +
  theme_classic()

p3 <- ggplot(df_summarized, aes(x = factor(time), y = median_ctrough, ymin = LCL_ctrough, ymax = UCL_ctrough)) + 
  geom_point(aes(color = "median"), size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = LCL_ctrough, color = "LCL"), shape = 17, size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = UCL_ctrough, color = "UCL"), shape = 14, size = 4) +
  scale_color_manual(values = colors) + 
  guides(color = guide_legend(title = "Simulation")) +
  facet_wrap(~Conditions, ncol = 4) +
  ylim(0, 1.1)+
  xlab("Time") + ylab("Ctrough") +
  theme_classic()+
  theme(axis.line.x = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
 



legend <- cowplot::get_legend(p3+theme(legend.position="bottom"))

# Create the title
title <- ggdraw() + 
  draw_label("GFR-related differences in pharmacokinetic parameters in liver", fontface = "bold", size = 14, x = 0, hjust = -0.1)

# Remove the legend from the plots
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")

# Combine the plots into a single plot with a common legend
combined_plot <- grid.arrange(title,p1, p3, p2, legend, nrow = 5, heights = c(0.1,1, 1, 1, 0.2))






library(cowplot)


# Add a common legend
plot_annotation(title = "Common legend", theme = theme(plot.title = element_text(hjust = 0.5))) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "Value Type"))

  ggplot(df_summarized, aes(x = factor(time), y = ctrough_ctrough, ymin = LCL_ctrough, ymax = UCL_ctrough)) + 
    # geom_errorbar() + 
    geom_point(aes(color = "median_ctrough"), size = 4) +
    geom_point(data = df_summarized, aes(x = factor(time), y = LCL, color = "LCL"), shape = 17, size = 4) +
    geom_point(data = df_summarized, aes(x = factor(time), y = UCL, color = "UCL"), shape = 14, size = 4) +
    scale_color_manual(values = c("median_ctrough" = "green", "LCL_ctrough" = "blue", "UCL_ctrough" = "red")) +
    guides(color = guide_legend(title = "Value Type")) +
    facet_wrap(~Conditions, ncol = 4) +
    xlab("Time") + ylab("ctrough") +
    theme_classic()
  
  ggplot(df_summarized, aes(x = factor(time), y = median_Cmax, ymin = LCL_Cmax, ymax = UCL_Cmax)) + 
    # geom_errorbar() + 
    geom_point(aes(color = "median_Cmax"), size = 4) +
    geom_point(data = df_summarized, aes(x = factor(time), y = LCL_Cmax, color = "LCL_Cmax"), shape = 17, size = 4) +
    geom_point(data = df_summarized, aes(x = factor(time), y = UCL_Cmax, color = "UCL_Cmax"), shape = 14, size = 4) +
    geom_point(aes(x = factor(time), y = median_AUC, color = "median_AUC"), size = 4) +
    geom_point(data = df_summarized, aes(x = factor(time), y = LCL_AUC, color = "LCL_AUC"), shape = 17, size = 4) +
    geom_point(data = df_summarized, aes(x = factor(time), y = UCL_AUC, color = "UCL_AUC"), shape = 14, size = 4) +
    scale_color_manual(values = c("median_Cmax" = "green", "LCL_Cmax" = "blue", "UCL_Cmax" = "red", "median_AUC" = "orange", "LCL_AUC" = "purple", "UCL_AUC" = "black")) +
    guides(color = guide_legend(title = "Value Type")) +
    facet_wrap(~Conditions, ncol = 4) +
    xlab("Time") + ylab("Cmax and AUC") +
    theme_classic()
  
  
  
  
  
 
 

    temp = df_grouped %>% 
      filter(Times != 0) %>% 
      # filter(Times > df_grouped_max_time$Tmax[df_grouped_max_time$time.tmax == tp]) %>% 
      filter(Times > df_grouped_max_time$Tmax[df_grouped_max_time$Conditions %in% Conditions &
                                                df_grouped_max_time$variables %in% variables &
                                                df_grouped_max_time$dose_uM %in% dose_uM &
                                                df_grouped_max_time$value_type %in% value_type]) %>% 
      slice(which.min(value)) %>% 
      slice(1) %>%
      bind_rows()%>%
      rename(ctrough = value) %>% 
      rename(troughtime= Times)


#################
    # first 20 rows are paramters
    
    Simulationfile = Simulationfile1
    Final_Rdose1 = c(100)
    Conditions = c("1.NormalGFR", "2.mildGFR", "3.moderateGFR", "4.SeverGFR")
    df1 <- data.frame(Conditions)
    df2 <- data.frame(Final_Rdose1)
    
    Times = seq(0,24, by = 0.1) #experimental data time point
    biomarker = length(proteins)
    #Times = seq(0,24, by = 0.1) #experimental data time point
    simul = length(Final_Rdose1)*length(Conditions) #if we include dose zero to check the steady state
    #rename in captial as experimental data 
    Simulationfile$Simulation = (rep(paste0(rep("Simu",simul),1:simul),each = length(Times)*biomarker))  # as we have 3 variables
    Simulationfile$dose_uM = as.numeric(rep(paste0(rep(Final_Rdose1)),each = length(Times)*length(Conditions)))
    Simulationfile$Conditions = (rep(paste0(rep(Conditions)),each = length(Times)))
    Simulationfile$Chemical = rep(paste0("NFT"),each = nrow(Simulationfile))
    Simulationfile$Times <- Times
    Simulationfile$variables = (rep(paste0(rep(proteins)),each = length(Times)*simul))
    colnames(Simulationfile) <- c("median", "LCL", "UCL","mean","Simulations","dose_uM","Conditions","Chemicals","Times","variables")
    
    Simulationfile = Simulationfile %>%  mutate(.,Dose_condition = paste(dose_uM ,Conditions, sep = "_"))
    
    
    
    write.table(Simulationfile, file = paste0(getwd(),"/",Modeling_dir, "/", "SingledoseGFRcompromisedSimulation.out"), row.names = F, sep = "\t")
    


    plot_function = function (variable1) {
      Simulation = Simulationfile %>% filter(.,variables %in% variable1) 
      yaxis = ifelse(variable1 == "P_excreted", "Fraction of Dose excreted (%)","Concentration(mg/L)") 
      ggplot() +
        geom_line(data =Simulation,aes(x = Times, y = mean, group = Dose_condition, color="Simulated_mean"),size=1) +
        geom_line(data =Simulation,aes(x = Times, y = LCL,  group = Dose_condition, color="Simulated_p2.5"),size=1) +
        geom_line(data =Simulation,aes(x = Times, y = UCL,  group = Dose_condition,color="Simulated_p95"),size=1) +
        facet_grid(variables~Conditions) +
        labs(x = "Time (hr)", y = "Plasma conc.(mg/L)", sec.x="Conditions",
             title = paste0("NFT concentrations in Human plasma following a repeated oral dose")) +
        #scale_y_continuous(trans = 'log10')+
        #scale_y_continuous(limits=c(0,3)) +
        scale_color_manual(name = "type",
                           breaks = c("Simulated_mean","Simulated_p2.5","Simulated_p95"),
                           values = c("Simulated_mean" = "blue","Simulated_p2.5" = "green","Simulated_p95" = "red"))+
        
        theme(legend.title = element_blank()) + theme(legend.position="bottom") +
        theme(legend.text = element_text(size= 8,face="bold")) 
    }
    
    plot_function(variable1 = "P_excreted")




