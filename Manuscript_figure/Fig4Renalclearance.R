library(tidyverse)
library("gridExtra")
library(ggsci)
library(ggplot2)
library(dplyr)
library(MESS)
setwd("~/MCSim/mod")
rm(list = ls())
Species = "Rabbit"
modeling = c("V2","V3", "V4")
Modeling_dir = "NFT_PBPK2/Posterior_chains"

Simulation <- modeling %>%
  map_df(~ read.delim(paste0(getwd(), "/",  Modeling_dir, "/", Species, .x, "Simulation.out"), header = TRUE, sep = "") %>%
           mutate(modeling_version = .x)
  )


##########################
# For model 3

Simulationfile = Simulation %>% filter(., modeling_version == modeling[3])
proteins = unique(Simulationfile$variables)
Final_Rdose1 = unique(Simulationfile$dose_uM)

hysteresis = Simulationfile  %>% select(median, dose_uM,Times, variables) %>% filter(variables %in% c("cplasma","ckidney","cfilterate","Filteration","Secretion","Reabsorption") & Times !=0.00)
files = pivot_wider(hysteresis, names_from = variables, values_from = median)%>% mutate(., TotalCl = Filteration + Secretion)%>% 
        mutate(., GFR_contribution = (Filteration/TotalCl)*100) %>% mutate(., Secr_contribution = (Secretion/TotalCl)*100)%>% data.frame(.)


cplasma_max = list()
secretion_max = list()
Filteration_max = list()
ckidney_max = list()
GFR_contribution_max = list()
Secr_contribution_max = list()
Secr_contribution_minatcmax= list()
TotalCl_max = list()
tmax = list()
tmax_half_abs = list()
tmax_half_elimi = list()

for (i in 1:length(Final_Rdose1)){
  cplasma_max[[i]] = files %>% filter(., dose_uM == Final_Rdose1[i]) %>%  select(Times, cplasma)%>% dplyr::filter(cplasma == max(cplasma))
  secretion_max[[i]] = files %>% filter(., dose_uM == Final_Rdose1[i]) %>% select(Times, Secretion)%>% dplyr::filter(Secretion == max(Secretion))
  Filteration_max[[i]] = files %>% filter(., dose_uM == Final_Rdose1[i]) %>% select(Times, Filteration)%>% dplyr::filter(Filteration == max(Filteration)) 
  ckidney_max[[i]] = files %>% filter(., dose_uM == Final_Rdose1[i]) %>%  select(Times, ckidney)%>% dplyr::filter(ckidney == max(ckidney))
  tmax[[i]] = files %>% filter(., dose_uM == Final_Rdose1[i]) %>% select(Times, cplasma)%>% filter(.,cplasma == max(cplasma)) %>% select(Times) 
  tmax_half_abs[[i]] = files %>% filter(., dose_uM == Final_Rdose1[i]) %>% select(Times,cplasma)%>% dplyr::filter(Times < tmax[[i]][[1]][1]) %>% filter(abs(cplasma - max(cplasma)/2) == min(abs(cplasma - max(cplasma)/2))) 
  tmax_half_elimi[[i]] = files %>% filter(., dose_uM == Final_Rdose1[i]) %>% select(Times,cplasma)%>% dplyr::filter(Times > tmax[[i]][[1]][1]) %>% filter(abs(cplasma - max(cplasma)/2) == min(abs(cplasma - max(cplasma)/2))) 
  
  GFR_contribution_max[[i]] = files %>% filter(., dose_uM == Final_Rdose1[i]) %>%  select(Times, GFR_contribution)%>% dplyr::filter(GFR_contribution == max(GFR_contribution))
  Secr_contribution_max1 = files %>% filter(., dose_uM == Final_Rdose1[i]) %>%  select(Times, Secr_contribution)%>% dplyr::filter(Secr_contribution == max(Secr_contribution))
  Secr_contribution_max[[i]] = Secr_contribution_max1[1,]
  #extract the minimum after zero before less than 0.02 hr for different doses
  Secr_contribution_minatcmax[[i]] = files %>% filter(., dose_uM == Final_Rdose1[i]) %>%  select(Times, Secr_contribution)%>% dplyr::filter(Times >0 & Times <0.1) %>% dplyr::filter(Secr_contribution == min(Secr_contribution))
  TotalCl_max[[i]] = files %>% filter(., dose_uM == Final_Rdose1[i]) %>%  select(Times, TotalCl)%>% dplyr::filter(TotalCl == max(TotalCl))
  
}

max_tubule =do.call(rbind.data.frame, cplasma_max)
max_secretion = do.call(rbind.data.frame, secretion_max)
max_Filteration = do.call(rbind.data.frame, Filteration_max)
max_ckidney = do.call(rbind.data.frame, ckidney_max)
max_GFR_contribution = do.call(rbind.data.frame, GFR_contribution_max)
max_Secr_contribution = do.call(rbind.data.frame, Secr_contribution_max)
min_Secr_contributionatcmax = do.call(rbind.data.frame, Secr_contribution_minatcmax) 
colnames(min_Secr_contributionatcmax)[colnames(min_Secr_contributionatcmax) == "Secr_contribution"] <- "Secr_contribution_atCmax"
max_TotalCl = do.call(rbind.data.frame, TotalCl_max)

finaldata = cbind(max_tubule,max_ckidney, max_secretion,max_Filteration,max_GFR_contribution,max_Secr_contribution,min_Secr_contributionatcmax,max_TotalCl)

df <- finaldata[, unique(names(finaldata))]  #Remove the time as duplicate column

files = pivot_longer(df, cols = (Secretion:TotalCl)) %>% data.frame(.)

variables = unique(files$name)


# Modify values in 'name' column using dplyr

modified_files <- files %>% 
  mutate(modified_name = if_else(name == "Filteration", "GFR", 
                                         if_else(name == "Secretion", "Tubular Secretion",
                                                 if_else(name == "TotalCl", "Total Renal clearence",
                                                         if_else(name =="GFR_contribution", "GFR contribution (%)",
                                                                 if_else(name =="Secr_contribution", "Secr_contribution (%)",
                                                                 if_else(name =="Secr_contribution_atCmax", "ATS contribution (%)",
                                                                         name)))))))


# Check the modified values in the 'modified_name' column
modified_files$modified_name = factor(modified_files$modified_name, levels = c("GFR","Tubular Secretion","Total Renal clearence","GFR contribution (%)","Secr_contribution (%)","ATS contribution (%)"))

########################
#Figure 3 

p1 = modified_files %>%
  mutate(rounded_cplasma = round(cplasma, 2)) %>%
  filter(., modified_name %in% c("Tubular Secretion","GFR")) %>%
  ggplot(aes(x = factor(rounded_cplasma), y = value)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(x = "Plasma concentration max (mg/L)", y = "Rate (mg/L/hr)",
       title = "A.", size = 14) +
  facet_wrap(~ modified_name, nrow = 1, scales = "free")


p2 = modified_files %>% 
  mutate(rounded_cplasma = round(cplasma, 2)) %>%
  filter(., modified_name %in% c("GFR contribution (%)", "ATS contribution (%)")) %>% 
  ggplot(aes(x = factor(rounded_cplasma), y = value)) +
  geom_bar(stat = "identity", width = 0.5) +
  facet_wrap(~ modified_name, nrow = 1) +
  labs(x = "Plasma concentration max (mg/L)", y = "Percentage of Total Renal clearance",
       title = "B.", size = 14)


directory_i = "~/MCSim/mod/NFT_PBPK2/Manuscript_figure"


tiff(paste0(directory_i,"/", "Figure4",".png"), units="in", width=7, height=6, res=300, compression = 'lzw') 
grid.arrange(p1, p2, ncol = 1)
dev.off()


# Time and dose dependent increase in GFR and decrease in active tubular secretion
library(ggplot2)
library(ggnewscale)



Simulationfile = Simulation %>% filter(., modeling_version == modeling[3])
proteins = unique(Simulationfile$variables)
Final_Rdose1 = unique(Simulationfile$dose_uM)

hysteresis = Simulationfile  %>% select(median, dose_uM,Times, variables) %>% filter(variables %in% c("cplasma","ckidney","cfilterate","Filteration","Secretion","Reabsorption") & Times !=0.00)
files = pivot_wider(hysteresis, names_from = variables, values_from = median)%>% mutate(., TotalCl = Filteration + Secretion)%>% 
  mutate(., GFR_contribution = (Filteration/TotalCl)*100) %>% mutate(., Secr_contribution = (Secretion/TotalCl)*100)%>% data.frame(.)


hysteresis1 = Simulationfile  %>% select(median, dose_uM,Times, variables) 
# files = pivot_wider(hysteresis, names_from = variables, values_from = median)%>% mutate(., netGFR = Filteration -Reabsorption) %>% data.frame(.)
files
# timeplot = function (variable1, variable2){

p3 = hysteresis1 %>% filter(., variables == "Secretion") %>% 
  ggplot()+
  geom_line(aes(Times, median, col = as.factor(dose_uM)),size=1.2) +
  labs(x = "Time [h]", y = "Secretion rate (mg/L/hr)",
       title = "A.", size = 14) +
  theme(
    legend.position = "right",
    axis.text = element_text(size=12),
    plot.title = element_text(size = 12, face = "bold")) +
  scale_color_discrete(name = "Dose")
# }



plot_function4 = function (variable1,variable2,time,dose) {
  files[,"dose_uM"] <- factor(files[,"dose_uM"])
  p = files %>% filter(., dose_uM %in% dose) %>% filter(., Times %in% time) %>% 
    ggplot(aes_string(x = variable1, y = variable2,color = "dose_uM"))+
    geom_path(size= 1.2)+
    geom_point(aes(color = as.factor(Times)),size = 0.6)+ 
    facet_wrap(~dose_uM, scales = "free") +
    scale_x_continuous()+
    scale_y_continuous(trans = 'log10')+
    labs(x = "Plasma concentration (mg/L)", y = "Active secretion contribution (%)", sec.x="Exposure (mg/kg)",
         title = paste0("B.")) +
    theme(legend.title = element_blank()) + theme(legend.position="") +
    theme(plot.title = element_text(size = 12, face = "bold")) +
    theme(legend.text = element_text(size= 8,face="bold")) 
  
  p + geom_text(aes(label = Times), hjust = 0.5, vjust = 0, size = 3)
}



p4 = plot_function4(variable1 ="cplasma" ,variable2 = "Secr_contribution", 
                    time = c(0.01, 0.03, 0.2, 0.4, 0.8, 1.6), dose = Final_Rdose1[c(1, 2, 6)])
p4

tiff(paste0(directory_i,"/", "FigS2",".png"), units="in", width=7, height=6, res=300, compression = 'lzw') 
grid.arrange(p3, p4, ncol = 1)
dev.off()
# plot_function4(variable1 ="cfilterate" ,variable2 = "Secretion", time = c(0.01,0.03,0.2,0.4,0.8,1.6,4),dose = Final_Rdose1[c(6)])
# 
# plot_function4(variable1 ="cfilterate" ,variable2 = "Secr_contribution", time = c(0.01,0.04,0.08,0.2, 0.4,0.8,1.6),dose = Final_Rdose1[c(1,6)])
# 
# plot_function4(variable1 ="cfilterate" ,variable2 = "Secr_contribution", time = c(0.01,0.04,0.08,0.2, 0.4,0.8,1.6),dose = Final_Rdose1[c(1,6)])
# 
# plot_function4(variable1 ="cfilterate" ,variable2 = "GFR_contribution", time = c(0.01,0.04,0.08,0.2, 0.4,0.8,1.6),dose = Final_Rdose1)
# 
# plot_function4(variable1 ="Filteration" ,variable2 = "Secretion", time = c(0.01,0.04,0.08,0.2, 0.4,0.8,1.6),dose = Final_Rdose1)
# 
# #1.1,1.2,1.3,1.4,1.5,1.6,.17,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,3,3.5,4
# proteins



Ylab = c("Plasma conc.(mg/L)", "Percentage of Dose excreted","Cumulative urine (mg)", "Total Mass Balance check (mg)", "Kidney conc.(mg/L)","Liver conc.(mg/L)")
Title1 = c("NFT concentrations in Rabbit following a single intravenous dose","Dose normalized cumulative urine excretion in Rabbit following a single IV",
           "Cumulative urine excretion in Rabbit following a single IV","NFT concentrations in Rabbit following a single OD",
           "Cumulative urine excretion in Rabbit following a single OD")




plot_function1 = function (variable1,variable2, route, mainplot, labplot) {
  Simulation = Simulationfile %>% filter(.,variables == variable1)
  #Simulation = Simulationfile %>% filter(.,variables == variable & Simulations != "Simu1")
  exp_data = Expdata %>% filter(Variable == variable2) %>% filter(exposure_route == route)
  
  ggplot() +
    geom_line(data =Simulation,aes(x = Times, y = median, color="Simulated_median"),size=1) +
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
                       breaks = c("Simulated_median","Simulated_p2.5","Simulated_p97.5", "Exp_mean"),
                       values = c("Simulated_median" = "blue","Simulated_p2.5" = "green","Simulated_p97.5" = "red","Exp_mean" = "blue"))+
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
    geom_line(data =Simulation,aes(x = Times, y = median, color="Simulated_median"),size=1) +
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
    #scale_y_continuous(trans = 'log10')+
    # scale_y_continuous(limits=c(0,3)) +
    scale_color_manual(name = "type",
                       breaks = c("Simulated_median","Simulated_p2.5","Simulated_p97.5", "Exp_mean"),
                       values = c("Simulated_median" = "blue","Simulated_p2.5" = "green","Simulated_p97.5" = "red","Exp_mean" = "blue"))+
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
    geom_line(data =Simulation,aes(x = Times, y = median, color= as.factor(dose_uM))) +
    # geom_line(data =Simulation,aes(x = Times, y = LCL, color="Simulated_p2.5"),size=1) +
    # geom_line(data =Simulation,aes(x = Times, y = UCL, color="Simulated_p97.5"),size=1) +
    
    facet_wrap(~variables, scales = "free") +
    labs(x = "Time (hr)", y = labplot, sec.x="Exposure (mg/kg)",
         title = mainplot) +
    #scale_y_continuous(trans = 'log10')+
    # scale_y_continuous(limits=c(0,3)) +
    # scale_color_manual(name = "type",
    #                    breaks = c("Simulated_median","Simulated_p2.5","Simulated_p97.5"),
    #                    values = c("Simulated_median" = "blue","Simulated_p2.5" = "green","Simulated_p97.5" = "red"))+
    theme(legend.title = element_blank()) + theme(legend.position="bottom") +
    theme(legend.text = element_text(size= 8,face="bold")) 
}

plot_function3(variable1 = proteins[7:11],dose = Final_Rdose1, mainplot = Title1[5] ,labplot = "NFT conc.(mg/L)")

proteins = unique(Simulationfile$variables)

