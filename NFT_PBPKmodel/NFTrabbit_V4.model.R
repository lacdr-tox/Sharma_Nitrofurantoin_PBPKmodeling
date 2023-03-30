# PBPK model  for RABBITS FOR NFT
#######################################################################
#Development of PBPK using MCsim model for Nitrofurantoin
#date: 16/05/22
#######################################################################

#%%%%%%%%%%specific notes for models%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#models for Rabbit Nitrofurantoin  
#Model version 4: generic model with only GFR based elimination and liver based metabolism in rabbits +  
#now we have added the both the active secretion and reabsorption from the filterate comartment

#-----------------------------------------------------------------------------------

States =  {
  Agutlumen,
  Agut,		       #amount of Flutamide gut
  Aliver,		     #amount of Flutamide liver
  ALiverdeg,      # liver degradation 
  Afeces,      # liver degradation 
  Akidney,        #amount of Flutamide in kidney
  Afilterate,
  Adelay, 
  Afat,		       #amount of Flutamide fat
  Arestbody,      #amount of Flutamide rest of the body
  Aplasma,		     #amount of Flutamide blood plasma
  Aurine };       #amount of NFT in urine

Outputs = {cgut, cliver, cfat,ckidney,cfilterate, crestbody, cplasma,P_excreted,Massbalance,Filteration,Reabsorption,Secretion,Total_renal_Cl,GFR_contribution,Secr_contribution,Reabs_contribution,Re_GFR_contribution};	

Inputs = {oraldose1}; 	


BW = 2.5; 
####
#Physiolgical data like cardiac output and GFR and hematocrit fraction are taken from httk using physiologica.data command
#constant Fraction of blood flows to organs (blood flow rate)
QCC = 15.96;                  #Total Cardiac blood output (L/h/kg^0.75) not it is to the power 0.75 
HCT = 0.41;  				         #hematocrit percentage
FQliver = 0.1245; 		         #Fraction cardiac output going to liver 
FQkidney = 0.151;  		       #Fraction cardiac output going to kidney  
FQfat =  0.0604;  		         #Fraction cardiac output going to fat  
FQgut = 0.2094;   			       #Fraction cardiac output going to gut


#constant organ volume as a fraction of total body weight  # For rabit from hTTK
Fliver = 0.04; 		           #Fraction liver volume
Fkidney = 0.006; 			       #Fraction kidney volume 
Ffilterate=0.0006;           # 10 percent of kidney volume
Ffat = 0.048;                 #fractional volume of fat  
Fgut = 0.048;   				       #fractional volume of gut
Fplasma = 0.08594;  	         #fractional volume of plasma

expT = 0.001;

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# partition coefficient parameter of NFT in Rabbits
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kgut_plasma = 0.622;		 	     #Liver/blood partition coefficient
Kliver_plasma = 0.651;		 	   #Liver/blood partition coefficient  
Kkidney_plasma = 0.671;		   #kidney/blood partition coefficient
Kfat_plasma = 0.159;	       #Fat/blood partition coefficient  (fitted)
Krestbody_plasma = 0.65;       #Rest of the body/blood partition coefficient	(fitted)
kgutabs = 0; 
kurine = 1.25;                 						
frac = 0.5;
fu = 0.4149; #0.4149
QurineC = 0.0025;
kdegliver = 0.1;
VmaxC = 1.8;  #0.043; #mg/hr/g liver  (scaled from rat using human microsomal data see doc file for more detail scaling)
Km = 63;       #km value is not given
IVdosing = 0;
oraldose = 0;
kbp = 0.76;  #0.76;  #blood to plasma partioning
kfeces = 0;
Tmc = 25;
Kt = 1.25;  #mg/L
Trc = 139.86;



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Compile parameters to be computed  in initialized
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
QCblood;
QCplasma;
Qliver;
Qgut;
Qkidney;
Qfilterate;
Qfat;
Qrestbody;
vliver;
vgut;
vkidney;
vfilterate;
vfat;
vplasma;
vrestbody;
Qurine;
Vmax;
Tr;
Tm;


Initialize {
  Agutlumen = 0;
  Aplasma = IVdosing*BW;		                    #amount of Flutamide blood plasma
  QCblood = QCC*pow(BW,0.75); 		 			         #Initial cardiac output for blood L/h  
  #QCplasma =QCblood*(1-HCT);   					   #Adjust initial cardiac output for plasma flow  
  QCplasma =QCblood;   					   #Adjust initial cardiac output for plasma flow  
  Qliver= FQliver* QCplasma;   	 				   #Plasma flow to liver 									  	
  Qgut= FQgut* QCplasma; 
  Qkidney= FQkidney*QCplasma; 							 #Plasma flow to kidney
  Qfilterate= 0.631;   				                 #Glomerulation filteration rate L/hr/BW (rabbit paper) https://doi.org/10.1258/la.2012.011065
  Qfat = FQfat*QCplasma;                    #plasma flow to fat  
  Qrestbody = QCplasma  -(Qliver + Qkidney + Qfat+ Qgut);  #plasma flow to rest of the body 
  vgut = Fgut*BW;
  vliver = Fliver*BW;  									 #Liver Volume  
  vkidney = Fkidney*BW;										 #volume of kidney
  vfilterate =Ffilterate * BW;
  vfat = Ffat*BW;                   				 #Volume of  fat 
  vplasma = Fplasma*BW; 
  vrestbody = 0.84*BW - vliver- vkidney- vfat- vplasma -vgut -vfilterate;
  Qurine = QurineC*pow(BW,0.75);
  Vmax = VmaxC*1000*vliver;    # multiplied 1000 for converting per g liver into kg liver then to whole liver by vliver
  Tr = Trc*BW;  #mg/L/kg to whole body weight (mcg/ml to mg/l is equivalent)
  Tm = Tmc*BW;  #mg/L/kg to whole body weight (mcg/ml to mg/l is equivalent)
}

Dynamics {
  
  cgut = Agut/vgut;
  cliver = Aliver/vliver;
  ckidney = Akidney/vkidney;
  cfilterate=Afilterate/vfilterate;
  cfat = Afat/vfat;
  crestbody = Arestbody/vrestbody;
  cplasma = Aplasma/vplasma;
  
  Filteration = Qfilterate*ckidney*((fu/kbp)/Kkidney_plasma);
  Reabsorption = (Tr*cfilterate);   #/(Kr+cfilterate)
  Secretion = (Tm*(ckidney*((fu/kbp)/Kkidney_plasma)))/(Kt+(ckidney*((fu/kbp)/Kkidney_plasma)));
  #############
  
  dt(Agutlumen) = ((oraldose1*BW))/expT -kgutabs*Agutlumen;   
  
  dt(Agut) = kgutabs*Agutlumen + Qgut*(cplasma*(fu/kbp)- cgut*((fu/kbp)/Kgut_plasma)); 
  
  dt(Aliver) = Qliver*cplasma*(fu/kbp) + Qgut*cgut*((fu/kbp)/Kgut_plasma) - (Qliver +Qgut)*cliver*((fu/kbp)/Kliver_plasma)-(Vmax*cliver*fu)/(Km + cliver*fu);

  dt(ALiverdeg) = (Vmax*cliver*fu)/(Km + cliver*fu);
  
  # dt(Aliver) = Qliver*cplasma*(fu/kbp) + Qgut*cgut*((fu/kbp)/Kgut_plasma) - (Qliver +Qgut)*cliver*((fu/kbp)/Kliver_plasma)-kdegliver*Aliver;
  # 
  # dt(ALiverdeg) = kdegliver*Aliver;
  
  dt(Afeces) = kfeces*Agut;
  
  dt(Akidney) = Qkidney*cplasma*(fu/kbp) - Qkidney*ckidney*((fu/kbp)/Kkidney_plasma) -Filteration + Reabsorption -Secretion;
  
  dt(Afilterate) = Filteration - Qfilterate*cfilterate -Reabsorption + Secretion  ;  # ((Tm*(pow(cfilterate, nt)))/(pow(Kt,nt)+ pow(cfilterate,nt))); #(Tm*cfilterate)/(Kt+cfilterate); 
  
  dt(Adelay)  = Qfilterate*cfilterate -Qurine *Adelay;
  
  dt(Afat) = Qfat*cplasma*(fu/kbp) - Qfat*cfat*((fu/kbp)/Kfat_plasma);  #Filtrate_amount
  dt(Arestbody) = Qrestbody*cplasma*(fu/kbp) - Qrestbody*crestbody*((fu/kbp)/Krestbody_plasma);
  
  dt(Aplasma) = (Qliver +Qgut)*cliver*((fu/kbp)/Kliver_plasma)+ Qkidney*ckidney*((fu/kbp)/Kkidney_plasma) + 
    Qfat*cfat*((fu/kbp)/Kfat_plasma)+ Qrestbody*crestbody*((fu/kbp)/Krestbody_plasma)- QCplasma*cplasma*(fu/kbp); 
  
  dt(Aurine)  = Qurine *Adelay;    #Qurine is L/hr
  
}

CalcOutputs { 
  
  Total_renal_Cl = Filteration + Secretion - Reabsorption;
  
  GFR_contribution = (t>0?(Filteration/Total_renal_Cl)*100 : 0);
  Re_GFR_contribution = (t>0?((Filteration-Reabsorption)/Total_renal_Cl)*100 : 0);
  Secr_contribution = (t>0?(Secretion/Total_renal_Cl)*100: 0);
  Reabs_contribution = (t>0?(Reabsorption/Total_renal_Cl)*100: 0);

  #percentage fraction output of chemical in urine
  P_excreted = (Aurine/((IVdosing + oraldose)*BW))*100;
  
  
  Massbalance =  Agutlumen + Agut + Aliver + ALiverdeg +  Afeces + Akidney + Afilterate + Adelay + Afat + Arestbody + Aplasma + Aurine;
  
} # End of output calculation

End.