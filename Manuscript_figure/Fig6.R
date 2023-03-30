library("ggformula")
rm(list = ls())
# run different batches testing
parallel::detectCores()
library(rstudioapi)

Modeling_dir = "NFT_PBPK2/Posterior_chains"
directory_i = "~/MCSim/mod/NFT_PBPK2/Manuscript_figure"

Simulationfile =read.delim(paste0(getwd(),"/",Modeling_dir, "/", "HumanMontecarloSimulation.out"), header = TRUE, sep = "")



file_name = "NFT_female_PK.csv"
myfiles <- read.csv(paste0(getwd(), "/", Modeling_dir, "/", file_name))
head(myfiles)

Expdata = myfiles %>%  rename(variables = Variable) %>% 
  mutate(modified_name = if_else(variables == "cplasma", "Plasma conc. (mg/L)", 
                                 variables))

#Timeid = paste0(unique(myfiles$Time, incomparables = FALSE),sep = ",",collapse = " ") %>% gsub(",$", "", .)  #gsub to remove the last comma
Final_Rdose1 = c(50,100,200)
exp_route = "oraldose"


proteins = c("cplasma","P_excreted","B_excreted","Massbalance","Filteration","Reabsorption","Secretion","Total_renal_Cl","GFR_contribution","Secr_contribution","Reabs_contribution","Re_GFR_contribution","cgut", "cliver", "cfat","ckidney","cfilterate", "crestbody", "Agutlumen","Adelay","Aurine") 


modified_files <- Simulationfile %>% 
  mutate(modified_name = if_else(variables == "cplasma", "Plasma conc. (mg/L)", 
                                 if_else(variables == "P_excreted", "Dose excreted in urine (%)", 
                                         if_else(variables == "B_excreted", "Dose excreted in bile (%)",variables))))
modified_files$modified_name = factor(modified_files$modified_name, levels = c("Plasma conc. (mg/L)","Dose excreted in urine (%)","Dose excreted in bile (%)",
                                                                               "Massbalance","Filteration","Reabsorption","Secretion","Total_renal_Cl","GFR_contribution","Secr_contribution",
                                                                               "Reabs_contribution","Re_GFR_contribution","cgut", "cliver", "cfat","ckidney","cfilterate", "crestbody", 
                                                                               "Agutlumen","Adelay","Aurine"))


plot_function1 = function (variable, mainplot) {
  Simulation = modified_files %>% filter(.,variables %in% variable)
  exp_data = Expdata %>% filter(variables  %in% variable)   #oraldose
  ggplot() +
    geom_line(data =Simulation,aes(x = Times, y = median, color="Simulated_median"),size=1.2) +
    geom_line(data =Simulation,aes(x = Times, y = LCL, color="Simulated_p2.5"),size=1.2) +
    geom_line(data =Simulation,aes(x = Times, y = UCL, color="Simulated_p97.5"),size=1.2) +
    geom_point(data =exp_data,aes(x = Time, y = cplasma, color="Exp_mean"), size = 1.5) +
    geom_errorbar(data =exp_data,aes(x = Time,ymin=cplasma-cplasmaSD, ymax= cplasma +cplasmaSD,color = "Experiments(Sd)"), size=0.5,
                  width=.25)+
    facet_grid(modified_name~dose_uM, scales = "free") +
    #facet_grid(modified_name~dose_uM, labeller = label_both(cplasma = "Plasma", P_excreted = "Excreted"))+
    labs(x = "Time [h]", y = "Compartments", sec.x="Exposure (mg/kg)",
         title = mainplot, size = 12) +
    scale_y_continuous(trans = 'log10')+
    scale_x_continuous(breaks = c(2,6,12)) +
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


P1 = plot_function1(variable = c("cplasma"),mainplot ="")


tiff(paste0(directory_i,"/", "Figure6",".png"), units="in",width=8, height=8, res=300, compression = 'lzw') 
print(P1)
dev.off()


##############################
