rm(list = ls())

library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(MESS)
# library(PKNCA)
# library(AUC)
# library(pROC)
setwd("~/MCSim/mod")
library("ggformula")
Modeling_dir = "NFT_PBPK2/Posterior_chains"
directory_i = "~/MCSim/mod/NFT_PBPK2/Manuscript_figure"



#####################################################################
# data analysis
#####################################################################

Simulationfile =read.delim(paste0(getwd(),"/",Modeling_dir, "/", "GFRcompromisedSimulation.out"), header = TRUE, sep = "")
Simulationfile =read.delim(paste0(getwd(),"/",Modeling_dir, "/", "newGFRcompromisedSimulation.out"), header = TRUE, sep = "")

Simulationfile <- Simulationfile %>% 
  mutate(modified_conditions = if_else(Conditions == "4.SeverGFR", "4.SevereGFR", Conditions))

proteins = unique(Simulationfile$variables)
Final_Rdose1 = unique(Simulationfile$dose_uM)
Conditions = unique(Simulationfile$Conditions)
modified_conditions = unique(Simulationfile$modified_conditions)
glimpse(Simulationfile)

proteins = c("cplasma", "cliver")

df_grouped = Simulationfile %>% filter(., variables %in% proteins) %>% 
  select(variables, modified_conditions,dose_uM, Times, median, UCL, LCL) %>%
  pivot_longer(-c(variables, modified_conditions,dose_uM, Times), names_to="value_type", values_to="value") %>%
  group_by(variables, modified_conditions,dose_uM, time=cut(Times, breaks=c(0, 24, 48, 72, 96, 120), labels=c(24, 48, 72, 96, 120)))

df_grouped_max = df_grouped %>% 
  filter(Times != 0) %>% 
  group_by(modified_conditions, variables, dose_uM, value_type, time) %>% 
  slice(which.max(value)) %>% 
  bind_rows() %>% rename(Cmax = value)

df_grouped_AUC = df_grouped %>% filter(Times != 0) %>%
  group_by(modified_conditions, variables, dose_uM, value_type,time) %>% 
  do(data.frame(AUC = auc(.$Times, .$value, type = "spline"))) %>% 
  bind_rows() 


df_grouped_max_time = df_grouped %>% 
  filter(Times != 0) %>% 
  group_by(modified_conditions, variables, dose_uM, value_type,time) %>% 
  slice(which.max(value)) %>% 
  slice(1) %>%
  bind_rows()%>% 
  rename(Tmax = Times)


df_grouped_trough = df_grouped %>% 
  filter(Times != 0) %>% 
  # filter(Times > df_grouped_max_time$Tmax[df_grouped_max_time$time.tmax == tp]) %>% 
  filter(Times > df_grouped_max_time$Tmax[df_grouped_max_time$modified_conditions %in% modified_conditions &
                                            df_grouped_max_time$variables %in% variables &
                                            df_grouped_max_time$dose_uM %in% dose_uM &
                                            df_grouped_max_time$value_type %in% value_type &
                                            df_grouped_max_time$time %in% time]) %>% 
  group_by(modified_conditions, variables, dose_uM, value_type,time) %>% 
  slice(which.min(value)) %>% 
  slice(1) %>%
  bind_rows()%>%
  rename(ctrough = value) %>% 
  rename(troughtime= Times)


df_grouped_AUC = df_grouped_AUC %>% 
  rename(variables.AUC = variables,
         modified_conditions.AUC = modified_conditions,
         dose_uM.AUC = dose_uM,
         value_type.AUC = value_type,
         time.AUC = time)

df_grouped_max_time = df_grouped_max_time %>% 
  rename(variables.tmax = variables,
         modified_conditions.tmax = modified_conditions,
         dose_uM.tmax = dose_uM,
         value_type.tmax = value_type,
         time.tmax = time)

df_grouped_trough = df_grouped_trough %>% 
  rename(variables.trough = variables,
         modified_conditions.trough = modified_conditions,
         dose_uM.trough = dose_uM,
         value_type.trough = value_type,
         time.trough = time)


df_combined = bind_cols(df_grouped_max, df_grouped_AUC,df_grouped_max_time,df_grouped_trough)
df2 = df_combined %>%  select(variables, modified_conditions,  dose_uM,value_type,Cmax, time,AUC,ctrough)


df_summarized1 <-df2 %>% 
  group_by(time, modified_conditions, variables,dose_uM) %>% 
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

colors <- c("median" = "green", "LCL" = "blue", "UCL" = "red")

df_summarized = df_summarized1 %>%  filter(variables == "cliver"& dose_uM == 50 & time %in% c(24,48, 120))

p1 <- ggplot(df_summarized, aes(x = factor(time), y = median_Cmax, ymin = LCL_Cmax, ymax = UCL_Cmax)) + 
  geom_point(aes(color = "median"),shape = 17, size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = LCL_Cmax, color = "LCL"), shape = 17, size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = UCL_Cmax, color = "UCL"), shape = 17, size = 4) +
  scale_color_manual(values = colors) +
  guides(color = guide_legend(title = "Value Type")) +
  facet_wrap(~modified_conditions, ncol = 4) +
  ylim(0, 0.8)+
  theme_classic()+
  xlab("") + ylab("Cmax") +
  theme(axis.line.x = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())



p2 <- ggplot(df_summarized, aes(x = factor(time), y = median_AUC, ymin = LCL_AUC, ymax = UCL_AUC)) + 
  geom_point(aes(color = "median"),shape = 17, size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = LCL_AUC, color = "LCL"), shape = 17, size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = UCL_AUC, color = "UCL"), shape = 17, size = 4) +
  scale_color_manual(values = colors) + 
  guides(color = guide_legend(title = "Value Type")) +
  facet_wrap(~modified_conditions, ncol = 4) +
  ylim(0, 15)+
  xlab("") + ylab("AUC") +
  theme_classic()

p3 <- ggplot(df_summarized, aes(x = factor(time), y = median_ctrough, ymin = LCL_ctrough, ymax = UCL_ctrough)) + 
  geom_point(aes(color = "median"), shape = 17,size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = LCL_ctrough, color = "LCL"), shape = 17, size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = UCL_ctrough, color = "UCL"), shape = 17, size = 4) +
  scale_color_manual(values = colors) + 
  guides(color = guide_legend(title = "Simulation")) +
  facet_wrap(~modified_conditions, ncol = 4) +
  ylim(0, 0.30)+
  xlab("Time") + ylab("Ctrough") +
  theme_classic()+
  theme(axis.line.x = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())




legend <- cowplot::get_legend(p3+theme(legend.position="bottom"))

# # Create the title
# title <- ggdraw() + 
#   draw_label("", fontface = "bold", size = 14, x = 0, hjust = -0.1)

# Remove the legend from the plots
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")

# # Combine the plots into a single plot with a common legend
# combined_plot <- grid.arrange(title,p1, p3, p2, legend, nrow = 5, heights = c(0.1,1, 1, 1, 0.2))
library(gridExtra)

tiff(paste0(directory_i,"/", "FigS10_liver",".png"), units="in",width=8, height=7, res=300, compression = 'lzw') 
grid.arrange(p1, p3, p2, legend, nrow = 4, heights = c(1, 1, 1, 0.2))
dev.off()




#############################

colors <- c("median" = "green", "LCL" = "blue", "UCL" = "red")

df_summarized = df_summarized1 %>%  filter(variables == "cplasma"& dose_uM == 50 & time %in% c(24,48, 120))


p1 <- ggplot(df_summarized, aes(x = factor(time), y = median_Cmax, ymin = LCL_Cmax, ymax = UCL_Cmax)) + 
  geom_point(aes(color = "median"),shape = 17, size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = LCL_Cmax, color = "LCL"), shape = 17, size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = UCL_Cmax, color = "UCL"), shape = 17, size = 4) +
  scale_color_manual(values = colors) +
  guides(color = guide_legend(title = "Value Type")) +
  facet_wrap(~modified_conditions, ncol = 4) +
  ylim(0, 1)+
  theme_classic()+
  xlab("") + ylab("Cmax") +
  theme(axis.line.x = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())



p2 <- ggplot(df_summarized, aes(x = factor(time), y = median_AUC, ymin = LCL_AUC, ymax = UCL_AUC)) + 
  geom_point(aes(color = "median"),shape = 17, size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = LCL_AUC, color = "LCL"), shape = 17, size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = UCL_AUC, color = "UCL"), shape = 17, size = 4) +
  scale_color_manual(values = colors) + 
  guides(color = guide_legend(title = "Value Type")) +
  facet_wrap(~modified_conditions, ncol = 4) +
  ylim(0, 25)+
  xlab("") + ylab("AUC") +
  theme_classic()

p3 <- ggplot(df_summarized, aes(x = factor(time), y = median_ctrough, ymin = LCL_ctrough, ymax = UCL_ctrough)) + 
  geom_point(aes(color = "median"), shape = 17,size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = LCL_ctrough, color = "LCL"), shape = 17, size = 4) +
  geom_point(data = df_summarized, aes(x = factor(time), y = UCL_ctrough, color = "UCL"), shape = 17, size = 4) +
  scale_color_manual(values = colors) + 
  guides(color = guide_legend(title = "Simulation")) +
  facet_wrap(~modified_conditions, ncol = 4) +
  ylim(0, 0.4)+
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
  draw_label("GFR-related differences in pharmacokinetic parameters in plasma", fontface = "bold", size = 14, x = 0, hjust = -0.1)

# Remove the legend from the plots
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")

# # Combine the plots into a single plot with a common legend
# combined_plot <- grid.arrange(title,p1, p3, p2, legend, nrow = 5, heights = c(0.1,1, 1, 1, 0.2))
library(gridExtra)

tiff(paste0(directory_i,"/", "Fig.7Plasma",".png"), units="in",width=8, height=7, res=300, compression = 'lzw') 
grid.arrange(p1, p3, p2, legend, nrow = 4, heights = c(1, 1, 1, 0.1))
dev.off()



#######################
Simulationfile =read.delim(paste0(getwd(),"/",Modeling_dir, "/", "SingledoseGFRcompromisedSimulation.out"), header = TRUE, sep = "")
#Simulationfile =read.delim(paste0(getwd(),"/",Modeling_dir, "/", "200mg_SingledoseGFRcompromisedSimulation.out"), header = TRUE, sep = "")

Simulationfile <- Simulationfile %>% 
  mutate(modified_conditions = if_else(Conditions == "4.SeverGFR", "4.SevereGFR", Conditions))

proteins = unique(Simulationfile$variables)
Final_Rdose1 = unique(Simulationfile$dose_uM)
Conditions = unique(Simulationfile$Conditions)
modified_conditions = unique(Simulationfile$modified_conditions)

proteins = c("cplasma","P_excreted","B_excreted","Massbalance","Filteration","Reabsorption","Secretion","Total_renal_Cl","GFR_contribution","Secr_contribution","Reabs_contribution","Re_GFR_contribution","cgut", "cliver", "cfat","ckidney","cfilterate", "crestbody", "Agutlumen","Adelay","Aurine") 



modified_files <- Simulationfile %>% 
  mutate(modified_name = if_else(variables == "cplasma", "Plasma conc. (mg/L)", 
                                 if_else(variables == "P_excreted", "Dose excreted in urine (%)", 
                                         if_else(variables == "B_excreted", "Dose excreted in bile (%)",
                                                 if_else(variables == "cliver", "Liver conc. (mg/L)",
                                                         if_else(variables == "cfat", "Fat",
                                                                 if_else(variables == "ckidney", "Kidney",
                                                                         if_else(variables == "cfilterate", "Tubules",
                                                                                 if_else(variables == "crestbody","Rest of the body",variables)))))))))


modified_files$modified_name = factor(modified_files$modified_name, levels = c("Plasma conc. (mg/L)","Liver conc. (mg/L)", "Dose excreted in urine (%)","Dose excreted in bile (%)",
                                                                               "Massbalance","Filteration","Reabsorption","Secretion","Total_renal_Cl","GFR_contribution","Secr_contribution",
                                                                               "Reabs_contribution","Re_GFR_contribution","cgut", "Kidney","Tubules","Fat", "Rest of the body", 
                                                                               "Agutlumen","Adelay","Aurine"))



##################
plot_function = function (variable) {
  Simulation = modified_files %>% filter(.,variables %in% variable)
  #yaxis = ifelse(variable1 == "P_excreted", "Fraction of Dose excreted (%)","Concentration(mg/L)") 
  ggplot() +
    geom_line(data =Simulation,aes(x = Times, y = median, group = Dose_condition, color="Simulated_median"),size=1) +
    geom_line(data =Simulation,aes(x = Times, y = LCL,  group = Dose_condition, color="Simulated_p2.5"),size=1) +
    geom_line(data =Simulation,aes(x = Times, y = UCL,  group = Dose_condition,color="Simulated_p95"),size=1) +
    facet_grid(modified_name~modified_conditions, scales = "free") +
    labs(x = "Time (hr)", y = "", sec.x="Conditions",
         title = paste0("")) +
    #scale_y_continuous(trans = 'log10')+
    #scale_y_continuous(limits=c(0,3)) +
    scale_color_manual(name = "type",
                       breaks = c("Simulated_median","Simulated_p2.5","Simulated_p95"),
                       values = c("Simulated_median" = "blue","Simulated_p2.5" = "green","Simulated_p95" = "red"))+
    
    theme(legend.title = element_blank()) + theme(legend.position="bottom") +
    theme(legend.text = element_text(size= 8,face="bold")) 
}

p1 = plot_function(variable = c("cplasma","cliver","P_excreted"))

directory_i = "~/MCSim/mod/NFT_forwardPBPK/Manuscript_figure"

tiff(paste0(directory_i,"/", "FigS11",".png"), units="in",width=8, height=10, res=300, compression = 'lzw') 
print(p1)
dev.off()
