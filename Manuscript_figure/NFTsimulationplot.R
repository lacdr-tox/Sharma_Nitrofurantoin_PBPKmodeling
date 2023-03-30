
#This code would generate fig.2 and fig.3 for the main manuscript and annex.1 figure


library(tidyverse)
library("gridExtra")
library(ggsci)
library(ggplot2)
setwd("~/MCSim/mod")
rm(list = ls())
Species = "Rabbit"
modeling = c("V1","V2","V3", "V4")
Modeling_dir = "NFT_PBPK2/Posterior_chains"
directory_i = "~/MCSim/mod/NFT_PBPK2/Manuscript_figure"

Simulation <- modeling %>%
  map_df(~ read.delim(paste0(getwd(), "/", Modeling_dir, "/", Species, .x, "Simulation.out"), header = TRUE, sep = "") %>%
           mutate(modeling_version = .x)
  )


proteins = unique(Simulation$variables)
Final_Rdose1 = unique(Simulation$dose_uM)


Simulationfile = Simulation %>% filter(., modeling_version == modeling[1])

proteins = Simulationfile%>% select(variables) %>% unique() %>% pull()

file_name = "Rabbit_PKdata.csv"
myfiles <- read.csv(paste0(getwd(), "/", Modeling_dir, "/", file_name))
head(myfiles)

Expdata = myfiles %>%  rename(variables = Variable) %>% 
  mutate(modified_name = if_else(variables == "cplasma", "Plasma conc. (mg/L)", 
                                 if_else(variables == "P_excreted", "Percentage of dose excreted (%)", variables)))

Expdata$modified_name = factor(Expdata$modified_name, levels = c("Plasma conc. (mg/L)", "Percentage of dose excreted (%)"))



modified_files <- Simulationfile %>% 
  mutate(modified_name = if_else(variables == "cplasma", "Plasma conc. (mg/L)", 
                                 if_else(variables == "P_excreted", "Percentage of dose excreted (%)", variables)))
modified_files$modified_name = factor(modified_files$modified_name, levels = c("Plasma conc. (mg/L)","Percentage of dose excreted (%)","Adelay","Agutlumen","Aurine" ,"cfat",                         
                                                                                       "cfilterate","cgut" , "ckidney" ,  "cliver",                       
                                                                                       "crestbody","Filterate" ,"Massbalance"))

Ylab = c("Plasma conc.(mg/L)", "Percentage of Dose excreted","Cumulative urine (mg)", "Total Mass Balance check (mg)", "Kidney conc.(mg/L)","Liver conc.(mg/L)")
Title1 = c("NFT concentrations in Rabbit following a single intravenous dose","Dose normalized cumulative urine excretion in Rabbit following a single IV",
           "Cumulative urine excretion in Rabbit following a single IV","NFT concentrations in Rabbit following a single OD",
           "Cumulative urine excretion in Rabbit following a single OD")


plot_function1 = function (variable, mainplot) {
  Simulation = modified_files %>% filter(.,variables %in% variable)
  exp_data = Expdata %>% filter(variables  %in% variable) %>% filter(exposure_route == "IV")
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

P1 = plot_function1(variable = c("cplasma", "P_excreted"),mainplot ="A.")



tiff(paste0(directory_i,"/", "Figure2",".png"), units="in",width=8, height=6, res=300, compression = 'lzw') 
print(P1)
dev.off()


#plot for modeling version 2 and 3 will go to annex and combined together

Simulationfile2 = Simulation %>% filter(., modeling_version == modeling[2])

proteins = Simulationfile2%>% select(variables) %>% unique() %>% pull()

file_name = "Rabbit_PKdata.csv"
myfiles <- read.csv(paste0(getwd(), "/", Modeling_dir, "/", file_name))
head(myfiles)

Expdata = myfiles %>%  rename(variables = Variable) %>% 
  mutate(modified_name = if_else(variables == "cplasma", "Plasma conc. (mg/L)", 
                                 if_else(variables == "P_excreted", "Percentage of dose excreted (%)", variables)))

Expdata$modified_name = factor(Expdata$modified_name, levels = c("Plasma conc. (mg/L)", "Percentage of dose excreted (%)"))



modified_files <- Simulationfile2 %>% 
  mutate(modified_name = if_else(variables == "cplasma", "Plasma conc. (mg/L)", 
                                 if_else(variables == "P_excreted", "Percentage of dose excreted (%)", variables)))
modified_files$modified_name = factor(modified_files$modified_name, levels = c("Plasma conc. (mg/L)","Percentage of dose excreted (%)","Adelay","Agutlumen","Aurine" ,"cfat",                         
                                                                               "cfilterate","cgut" , "ckidney" ,  "cliver",                       
                                                                               "crestbody","Filterate" ,"Massbalance"))

Ylab = c("Plasma conc.(mg/L)", "Percentage of Dose excreted","Cumulative urine (mg)", "Total Mass Balance check (mg)", "Kidney conc.(mg/L)","Liver conc.(mg/L)")
Title1 = c("NFT concentrations in Rabbit following a single intravenous dose","Dose normalized cumulative urine excretion in Rabbit following a single IV",
           "Cumulative urine excretion in Rabbit following a single IV","NFT concentrations in Rabbit following a single OD",
           "Cumulative urine excretion in Rabbit following a single OD")


plot_function1 = function (variable, mainplot) {
  Simulation = modified_files %>% filter(.,variables %in% variable)
  exp_data = Expdata %>% filter(variables  %in% variable) %>% filter(exposure_route == "IV")
  ggplot() +
    geom_line(data =Simulation,aes(x = Times, y = median, color="Simulated_median"),size=1.2) +
    geom_line(data =Simulation,aes(x = Times, y = LCL, color="Simulated_p2.5"),size=1.2) +
    geom_line(data =Simulation,aes(x = Times, y = UCL, color="Simulated_p97.5"),size=1.2) +
    geom_point(data =exp_data,aes(x = Time, y = cplasma, color="Exp_mean"), size = 1.5) +
    
    facet_grid(modified_name~dose_uM, scales = "free") +
    #facet_grid(modified_name~dose_uM, labeller = label_both(cplasma = "Plasma", P_excreted = "Excreted"))+
    labs(x = "", y = "Compartments", sec.x="Exposure (mg/kg)",
         title = mainplot, size = 12) +
    scale_y_continuous(trans = 'log10')+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    scale_color_manual(name = "type",
                       breaks = c("Simulated_median","Simulated_p2.5","Simulated_p97.5", "Exp_mean"),
                       values = c("Simulated_median" = "blue","Simulated_p2.5" = "green","Simulated_p97.5" = "red","Exp_mean" = "blue"))+
    theme(strip.text = element_text(size = 7)) +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size= 8,face="bold"),
          axis.text = element_text(size=12),
          plot.title = element_text(size = 12, face = "bold")) 
}

P2 = plot_function1(variable = c("cplasma", "P_excreted"),mainplot ="A.")

Simulationfile = Simulation %>% filter(., modeling_version == modeling[3])

proteins = Simulationfile%>% select(variables) %>% unique() %>% pull()

file_name = "Rabbit_PKdata.csv"
myfiles <- read.csv(paste0(getwd(), "/", Modeling_dir, "/", file_name))
head(myfiles)

Expdata = myfiles %>%  rename(variables = Variable) %>% 
  mutate(modified_name = if_else(variables == "cplasma", "Plasma conc. (mg/L)", 
                                 if_else(variables == "P_excreted", "Percentage of dose excreted (%)", variables)))

Expdata$modified_name = factor(Expdata$modified_name, levels = c("Plasma conc. (mg/L)", "Percentage of dose excreted (%)"))



modified_files <- Simulationfile %>% 
  mutate(modified_name = if_else(variables == "cplasma", "Plasma conc. (mg/L)", 
                                 if_else(variables == "P_excreted", "Percentage of dose excreted (%)", variables)))
modified_files$modified_name = factor(modified_files$modified_name, levels = c("Plasma conc. (mg/L)","Percentage of dose excreted (%)","Adelay","Agutlumen","Aurine" ,"cfat",                         
                                                                               "cfilterate","cgut" , "ckidney" ,  "cliver",                       
                                                                               "crestbody","Filterate" ,"Massbalance"))

Ylab = c("Plasma conc.(mg/L)", "Percentage of Dose excreted","Cumulative urine (mg)", "Total Mass Balance check (mg)", "Kidney conc.(mg/L)","Liver conc.(mg/L)")
Title1 = c("NFT concentrations in Rabbit following a single intravenous dose","Dose normalized cumulative urine excretion in Rabbit following a single IV",
           "Cumulative urine excretion in Rabbit following a single IV","NFT concentrations in Rabbit following a single OD",
           "Cumulative urine excretion in Rabbit following a single OD")


plot_function2 = function (variable, mainplot) {
  Simulation = modified_files %>% filter(.,variables %in% variable)
  exp_data = Expdata %>% filter(variables  %in% variable) %>% filter(exposure_route == "IV")
  ggplot() +
    geom_line(data =Simulation,aes(x = Times, y = median, color="Simulated_median"),size=1.2) +
    geom_line(data =Simulation,aes(x = Times, y = LCL, color="Simulated_p2.5"),size=1.2) +
    geom_line(data =Simulation,aes(x = Times, y = UCL, color="Simulated_p97.5"),size=1.2) +
    geom_point(data =exp_data,aes(x = Time, y = cplasma, color="Exp_mean"), size = 1.5) +
    
    facet_grid(modified_name~dose_uM, scales = "free") +
    #facet_grid(modified_name~dose_uM, labeller = label_both(cplasma = "Plasma", P_excreted = "Excreted"))+
    labs(x = "Time [h]", y = "Compartments", sec.x="",
         title = mainplot, size = 12) +
    scale_y_continuous(trans = 'log10')+
    scale_color_manual(name = "type",
                       breaks = c("Simulated_median","Simulated_p2.5","Simulated_p97.5", "Exp_mean"),
                       values = c("Simulated_median" = "blue","Simulated_p2.5" = "green","Simulated_p97.5" = "red","Exp_mean" = "blue"))+
    theme(strip.text = element_text(size = 8)) +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size= 7,face="bold"),
          axis.text = element_text(size=12),
          plot.title = element_text(size = 12, face = "bold")) 
}

P3 = plot_function2(variable = c("cplasma", "P_excreted"),mainplot ="B.")


tiff(paste0(directory_i,"/", "Fig.S1",".png"), units="in", width=8, height=10, res=300, compression = 'lzw') 
grid.arrange(arrangeGrob(P2 + theme(legend.position="none"),
                         P3 + theme(legend.position="bottom"), nrow=2,heights = c(4.5, 5.5)))
dev.off()


#######################################
#Fig3

Simulationfile = Simulation %>% filter(., modeling_version == modeling[4])
proteins = Simulationfile%>% select(variables) %>% unique() %>% pull()

file_name = "Rabbit_PKdata.csv"
myfiles <- read.csv(paste0(getwd(), "/", Modeling_dir, "/", file_name))
head(myfiles)

Expdata = myfiles %>%  rename(variables = Variable) %>% 
  mutate(modified_name = if_else(variables == "cplasma", "Plasma conc. (mg/L)", 
                                 if_else(variables == "P_excreted", "Percentage of dose excreted (%)", variables)))

Expdata$modified_name = factor(Expdata$modified_name, levels = c("Plasma conc. (mg/L)", "Percentage of dose excreted (%)"))



modified_files <- Simulationfile %>% 
  mutate(modified_name = if_else(variables == "cplasma", "Plasma conc. (mg/L)", 
                                 if_else(variables == "P_excreted", "Percentage of dose excreted (%)", variables)))
modified_files$modified_name = factor(modified_files$modified_name, levels = c("Plasma conc. (mg/L)","Percentage of dose excreted (%)","Adelay","Agutlumen","Aurine" ,"cfat",                         
                                                                               "cfilterate","cgut" , "ckidney" ,  "cliver",                       
                                                                               "crestbody","Filterate" ,"Massbalance"))

Ylab = c("Plasma conc.(mg/L)", "Percentage of Dose excreted","Cumulative urine (mg)", "Total Mass Balance check (mg)", "Kidney conc.(mg/L)","Liver conc.(mg/L)")
Title1 = c("NFT concentrations in Rabbit following a single intravenous dose","Dose normalized cumulative urine excretion in Rabbit following a single IV",
           "Cumulative urine excretion in Rabbit following a single IV","NFT concentrations in Rabbit following a single OD",
           "Cumulative urine excretion in Rabbit following a single OD")


plot_function1 = function (variable, mainplot) {
  Simulation = modified_files %>% filter(.,variables %in% variable)
  exp_data = Expdata %>% filter(variables  %in% variable) %>% filter(exposure_route == "IV")
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

P5 = plot_function1(variable = c("cplasma", "P_excreted"),mainplot ="")


tiff(paste0(directory_i,"/", "Figure3",".png"), units="in", width=8, height=6, res=300, compression = 'lzw') 
print(P5)
dev.off()
