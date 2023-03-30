#Forward simulation of the extrapolation from rabbit to rats using model V4 using posterior paramters distribution 
#This code would generate the Figure S4. 

library(tidyverse)
library("gridExtra")
library(ggsci)
library(ggplot2)
setwd("~/MCSim/mod")
rm(list = ls())

# run different batches testing
parallel::detectCores()
library(rstudioapi)

Modeling_dir = "NFT_PBPK2/Posterior_chains"

Species = "rats"
modelingversion = "V4"

Simulationfile =read.delim(paste0(getwd(),"/",Modeling_dir, "/",Species,modelingversion, "ExtrapolatedV4Simulation.out"), header = TRUE, sep = "")

file_name = "Rat_PKdata.csv"
myfiles <- read.csv(paste0(getwd(), "/", Modeling_dir, "/", file_name))
head(myfiles)

proteins = c("cplasma","P_excreted","Massbalance","Filteration","Reabsorption","Secretion","Total_renal_Cl","GFR_contribution","Secr_contribution","Reabs_contribution","Re_GFR_contribution","cgut", "cliver", "cfat","ckidney","cfilterate", "crestbody", "Agutlumen","Adelay","Aurine") 

exposure_route = c("IV", "oraldose")



Expdata = myfiles %>%  rename(variables = Variable) %>% 
  mutate(modified_name = if_else(variables == "cplasma", "Plasma conc. (mg/L)", 
                                 if_else(variables == "P_excreted", "Perecentage dose excreted (%)", variables)))

Expdata$modified_name = factor(Expdata$modified_name, levels = c("Plasma conc. (mg/L)", "Perecentage dose excreted (%)"))



modified_files <- Simulationfile %>% 
  mutate(modified_name = if_else(variables == "cplasma", "Plasma conc. (mg/L)", 
                                 if_else(variables == "P_excreted", "Perecentage dose excreted (%)", variables)))
modified_files$modified_name = factor(modified_files$modified_name, levels = c("Plasma conc. (mg/L)","Perecentage dose excreted (%)","Adelay","Agutlumen","Aurine" ,"cfat",                         
                                                                               "cfilterate","cgut" , "ckidney" ,  "cliver",                       
                                                                               "crestbody","Filterate" ,"Massbalance"))

Ylab = c("Plasma conc.(mg/L)", "Percentage of Dose excreted","Cumulative urine (mg)", "Total Mass Balance check (mg)", "Kidney conc.(mg/L)","Liver conc.(mg/L)")
Title1 = c("NFT concentrations in Rat following a single intravenous dose","Dose normalized cumulative urine excretion in Rat following a single IV",
           "Cumulative urine excretion in Rat following a single IV","NFT concentrations in Rat following a single OD",
           "Cumulative urine excretion in Rat following a single OD")


plot_function1 = function (variable, mainplot) {
  Simulation = modified_files %>% filter(.,variables %in% variable)
  exp_data = Expdata %>% filter(variables  %in% variable) %>% filter(exposure_route == "Ivdose")
  ggplot() +
    geom_line(data =Simulation,aes(x = Times, y = median, color="Simulated_median"),size=1.2) +
    geom_line(data =Simulation,aes(x = Times, y = LCL, color="Simulated_p2.5"),size=1.2) +
    geom_line(data =Simulation,aes(x = Times, y = UCL, color="Simulated_p97.5"),size=1.2) +
    geom_point(data =exp_data,aes(x = Time, y = cplasma, color="Exp_mean"), size = 1.5) +
    
    facet_grid(modified_name~dose_uM, scales = "free") +
    #facet_grid(modified_name~dose_uM, labeller = label_both(cplasma = "Plasma", P_excreted = "Excreted"))+
    labs(x = "Time [h]", y = "Compartments", sec.x="Exposure (mg/kg)",
         title = mainplot, size = 12) +
    scale_y_continuous(trans = 'log10')+
    scale_color_manual(name = "type",
                       breaks = c("Simulated_median","Simulated_p2.5","Simulated_p97.5", "Exp_mean"),
                       values = c("Simulated_median" = "blue","Simulated_p2.5" = "green","Simulated_p97.5" = "red","Exp_mean" = "blue"))+
    theme(strip.text = element_text(size = 8)) +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size= 8,face="bold"),
          axis.text = element_text(size=12),
          plot.title = element_text(size = 12, face = "bold")) 
}

P1 = plot_function1(variable = c("cplasma", "P_excreted"),mainplot ="")

directory_i = "~/MCSim/mod/NFT_PBPK2/Manuscript_figure"


tiff(paste0(directory_i,"/", "FigS3",".png"), units="in",width=8, height=6, res=300, compression = 'lzw') 
print(P1)
dev.off()

