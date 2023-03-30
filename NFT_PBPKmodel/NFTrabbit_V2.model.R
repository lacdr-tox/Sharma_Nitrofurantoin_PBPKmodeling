# PBPK model  for RABBITS FOR NFT
#######################################################################
#Development of PBPK using MCsim model for Nitrofurantoin
#date: 16/05/22
#######################################################################

#%%%%%%%%%%specific notes for models%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#models for Rabbit Nitrofurantoin  
#Model version 2: generic model with only GFR based elimination and liver based metabolism in rabbits +  
#now we have added the active secretion in the filterate comartment

#-----------------------------------------------------------------------------------

States =  {
  Agutlumen,
  Agut,		       #amount of Flutamide gut
  Aliver,		     #amount of Flutamide liver
  ALiverdeg,      # liver degradation 
  Akidney,        #amount of Flutamide in kidney
  Afilterate,
  Adelay, 
  Afat,		       #amount of Flutamide fat
  Arestbody,      #amount of Flutamide rest of the body
  Aplasma,		     #amount of Flutamide blood plasma
  Aurine };       #amount of NFT in urine

Outputs = {cgut, cliver, cfat,ckidney,cfilterate, crestbody, cplasma,P_excreted,Massbalance,Secretion,Filteration,Total_renal_Cl,GFR_contribution,Secr_contribution};	

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

FQfilterate = 4.21;             # GFR rate in rabbits ml/min/kg  #https://journals.sagepub.com/doi/10.1258/la.2012.011065?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed#bibr-LA-11-065C17


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
Kbp = 0.76;  #0.76;  #blood to plasma partioning
kfeces = 0;
Tmc = 25;
Kt = 1.25;  #mg/L
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
Tm;


Initialize {
  Agutlumen = 0;
  Aplasma = IVdosing*BW;		                    #amount of Flutamide blood plasma
  QCblood = QCC*pow(BW,0.75); 		 			         #Initial cardiac output for blood L/h  
 # QCplasma =QCblood*(1-HCT);   					   #Adjust initial cardiac output for plasma flow 
  QCplasma =QCblood;   					   #Adjust initial cardiac output for plasma flow 
  Qliver= FQliver* QCplasma;   	 				   #Plasma flow to liver 									  	
  Qgut= FQgut* QCplasma; 
  Qkidney= FQkidney*QCplasma; 							 #Plasma flow to kidney
  Qfilterate= 0.63; 				                 #Glomerulation filteration rate L/hr/BW (rabbit paper) https://doi.org/10.1258/la.2012.011065 
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
  Secretion = (Tm*(ckidney*((fu/Kbp)/Kkidney_plasma)))/(Kt+ (ckidney*((fu/Kbp)/Kkidney_plasma)));
  Filteration = Qfilterate*ckidney*((fu/Kbp)/Kkidney_plasma);
  #############
  
  dt(Agutlumen) = ((oraldose1*BW))/expT -kgutabs*Agutlumen;   
  
  dt(Agut) = kgutabs*Agutlumen + Qgut*(cplasma*(fu/Kbp)- cgut*((fu/Kbp)/Kgut_plasma)); 
  
  dt(Aliver) = Qliver*cplasma*(fu/Kbp) + Qgut*cgut*((fu/Kbp)/Kgut_plasma) - (Qliver +Qgut)*cliver*((fu/Kbp)/Kliver_plasma)-(Vmax*cliver*fu)/(Km + cliver*fu);

  dt(ALiverdeg) = (Vmax*cliver*fu)/(Km + cliver*fu);
  
  dt(Akidney) = Qkidney*cplasma*(fu/Kbp) - Qkidney*ckidney*((fu/Kbp)/Kkidney_plasma) -Filteration -Secretion;
  
  dt(Afilterate) = Filteration - Qfilterate*cfilterate + Secretion ;
  
  dt(Adelay)  = Qfilterate*cfilterate -Qurine *Adelay ;
  
  dt(Afat) = Qfat*cplasma*(fu/Kbp) - Qfat*cfat*((fu/Kbp)/Kfat_plasma);  #Filtrate_amount
  dt(Arestbody) = Qrestbody*cplasma*(fu/Kbp) - Qrestbody*crestbody*((fu/Kbp)/Krestbody_plasma);
  
  dt(Aplasma) = (Qliver +Qgut)*cliver*((fu/Kbp)/Kliver_plasma)+ Qkidney*ckidney*((fu/Kbp)/Kkidney_plasma) + 
    Qfat*cfat*((fu/Kbp)/Kfat_plasma)+ Qrestbody*crestbody*((fu/Kbp)/Krestbody_plasma)- QCplasma*cplasma*(fu/Kbp); 
  
  dt(Aurine)  = Qurine *Adelay;    #Qurine is L/hr
  
}

CalcOutputs {   
  
  Total_renal_Cl = Filteration + Secretion;

  GFR_contribution = (t>0?(Filteration/Total_renal_Cl)*100 : 0);
  Secr_contribution = (t>0?(Secretion/Total_renal_Cl)*100: 0);

  #percentage fraction output of chemical in urine
  P_excreted = (Aurine/((IVdosing*BW + oraldose*BW)))*100;
  
  Massbalance =  Agutlumen + Agut + Aliver + ALiverdeg +  Akidney + Afilterate + Adelay + Afat + Arestbody + Aplasma + Aurine;
  
} # End of output calculation

End.