# PBPK model  for rats
#This model uses posterior paramter distribution of rabbit v4 model and simulate the rat kinetics


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

BWrabbit = 2.5;
BW = 0.25; 
####
#Physiolgical data like cardiac output and GFR and hematocrit fraction are taken from httk using physiologica.data command
#constant Fraction of blood flows to organs (blood flow rate)
QCC = 15.5;                  #Total Cardiac blood output (L/h/kg^0.75) not it is to the power 0.75   # check this
HCT = 0.46;  				         #hematocrit percentage
FQliver = 0.174; 		         #Fraction cardiac output going to liver 
FQkidney = 0.141;  		       #Fraction cardiac output going to kidney  
FQfat =  0.07;  		         #Fraction cardiac output going to fat  
FQgut = 0.075;   			       #Fraction cardiac output going to gut

#constant organ volume as a fraction of total body weight  # For rabit from hTTK
Fliver = 0.036; 		           #Fraction liver volume
Flung    = 0.006;             #fraction lung volume  		         		                    
Fkidney = 0.0073; 			       #Fraction kidney volume 
Ffilterate=0.00073;           # 10 percent of kidney volume
Ffat = 0.07;                 #fractional volume of fat  
Fgut = 0.027;   				       #fractional volume of gut
Fplasma = 0.074;  	         #fractional volume of plasma

# microsomal proteins content, mg / gr  human 
# MSPL = 52  				           #liver microsomal protein content 
# MSPG = 3     			           #mg/g intestine  
#cytosol protein content, mg/g 
# CYTPL = 80.7  			         #mg/g liver
# CYTPG = 18    			         #mg/g gut
# Tintestine = 2250    	       #total intestine weight in gram/Whole BW
expT = 0.001;

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# partition coefficient parameter of NFT in Rabbits
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kgut_plasma = 0.622;		 	     #Liver/blood partition coefficient
Kliver_plasma = 0.651;		 	   #Liver/blood partition coefficient  
Kkidney_plasma = 0.671;		   #kidney/blood partition coefficient
Kfat_plasma = 0.159;	       #Fat/blood partition coefficient  (fitted)
Krestbody_plasma = 0.159;       #Rest of the body/blood partition coefficient	(fitted)
kgutabs = 0; 
kurine = 1.25;                 						
frac = 0.5;
fu = 0.4149; #0.4149
QurineC = 11.36;
VmaxC = 0.527;  #0.043; #mg/hr/g liver  (scaled from rat using human microsomal data see doc file for more detail scaling)
Km = 59.16;       #km value is not given
IVdosing = 0;
oraldose = 0;
hepaticdose = 0;
kbp = 0.76;  #blood to plasma partioning
kfeces = 0;
Trc = 95.431;
Tmc = 1.516;
Kr = 1.34;  #mg/L
Kt = 0.0228;  #mg/L
scaling =1;
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
vlung;
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
  Qfilterate= 0.129;  				              # https://doi.org/10.1016/j.vascn.2009.10.002 Glomerulation filteration rate L/hr/BW rats (Onodera & Furuhama, 1983)
  Qfat = FQfat*QCplasma;                    #plasma flow to fat  
  Qrestbody = QCplasma  -(Qliver + Qkidney + Qfat+ Qgut);  #plasma flow to rest of the body 
  vgut = Fgut*BW;
  vliver = Fliver*BW;  									 #Liver Volume  
  vkidney = Fkidney*BW;										 #volume of kidney
  vfilterate =Ffilterate * BW;
  vfat = Ffat*BW;                   				 #Volume of  fat 
  vlung = Flung*BW; 
  vplasma = Fplasma*BW; 
  vrestbody = 0.84*BW - vliver- vkidney- vfat- vplasma -vgut -vfilterate;
  
  Qurine = QurineC*BW;
  Tm = Tmc*BW;  #mg/L/kg to whole body weight (mcg/ml to mg/l is equivalent)
  Vmax = VmaxC*1000*vliver; #multiplied 1000 for converting per g liver into kg liver then to whole liver by vliver
  Tr = Trc*BW;  #mg/L/kg to whole body weight (mcg/ml to mg/l is equivalent)
  
  # Qurine = QurineC*pow((BWrabbit/BW),scaling);
  # Tm = Tmc*pow((BWrabbit/BW),scaling1);  #mg/L/kg to whole body weight (mcg/ml to mg/l is equivalent)
  # Vmax = VmaxC*pow((BWrabbit/BW),scaling)*1000*vliver; #multiplied 1000 for converting per g liver into kg liver then to whole liver by vliver
  # Tr = Trc*pow((BWrabbit/BW),scaling);  #mg/L/kg to whole body weight (mcg/ml to mg/l is equivalent)
  # 
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
  
} 
  

End.