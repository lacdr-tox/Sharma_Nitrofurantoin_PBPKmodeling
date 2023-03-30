###############################################################
# put specific age if you want to, model 1 is for boys, and model 2 is for girls.
#Model has two things included first the age less than 12 has dose per kg BW after 12 the standard dose
# Second the model has included the possiblity of changing the  
# note we have calculated total renal clearance using the output of the model simulation, where
#we calculated total renal clearance using only fiteration that is (GFR) and tubular secretion that is secretion. 
#cifilterate, we refer this as a tubules in the manuscript. 
States =  {
  Agutlumen,
  Agut,		       #amount of Nitrofurantoin gut
  Aliver,		     #amount of Nitrofurantoin liver
  ALiverdeg,      # liver degradation 
  ABile,
  Afeces,      # liver degradation 
  Akidney,        #amount of Nitrofurantoin in kidney
  Afilterate,
  Adelay, 
  Afat,		       #amount of Nitrofurantoin fat
  Arestbody,      #amount of Nitrofurantoin rest of the body
  Aplasma,		     #amount of Nitrofurantoin blood plasma
  Aurine };       #amount of NFT in urine

Outputs = {cgut, cliver, cfat,ckidney,cfilterate, crestbody, cplasma,P_excreted,Massbalance,Filteration,Reabsorption,Secretion,Total_renal_Cl,GFR_contribution,Secr_contribution,Reabs_contribution,Re_GFR_contribution,B_excreted,DOSE};	

Inputs = {oraldose1}; 	

##############################################################################
#Both adult and pregnancy model 
BWrat = 0.25;
BWrabbit = 2.5;
expT = 0.001;

#################################################
#Physiological parameters
#######################################################
sex=1;
age = 30;
GFR = 1;
CKD = 0;  #(for normal GFR)


Ffilterate=0.0004;
#weight of non perfused bone at birth (3.5 kg)
Fgut_f=0.013;
Fgut_m=0.013;

#Fraction cardiac output

FQliver_f = 0.141;
FQliver_m = 0.125;

FQbrain_f=0.0812;                                               
FQbrain_m=0.057;
FQlung_f=0.034;
FQlung_m=0.034; 

FQkidney_f=0.0905;
FQkidney_m=0.0906;
FQfilterate_f=0.031;
FQfilterate_m=0.037;                                          
FQfat_f=0.044;
FQfat_m=0.024;

FQgut_f=0.14;
FQgut_m=0.18;



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# partition coefficient parameter of NFT in Rabbits
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kgut_plasma = 0.622;		 	     #Liver/blood partition coefficient
Kliver_plasma = 0.651;		 	   #Liver/blood partition coefficient  
Kkidney_plasma = 0.671;		   #kidney/blood partition coefficient
Kfat_plasma = 0.159;	       #Fat/blood partition coefficient  (fitted)
Krestbody_plasma = 0.423;       #Rest of the body/blood partition coefficient	(fitted)
kgutabs = 2.11; 
fu = 0.4149; #0.4149
QurineC = 13.45;
VmaxC = 0.47;  #0.043; #mg/hr/g liver  (scaled from rat using human microsomal data see doc file for more detail scaling)
Km = 5.83;       #km value is not given
IVdosing = 0;
oraldose = 0;
Kbp = 0.76;  #0.76;  #blood to plasma partioning
kfeces = 0;
Tmc = 8.024;
Kt = 0.059;  #mg/L
Trc = 1.334;
Vehrc = 0.527;
Kehr =0.017;
kbile = 0.25;
kstc = 1;
scaling = 0.25;

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Compile parameters to be computed  in initialized
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Initialize {
  Agutlumen = 0;
  Aplasma = IVdosing;		                    #amount of Nitrofurantoin blood plasma
 
}

Dynamics {
  
  # age = t/(365*24);
  model=((sex<=0)?1:(sex<=1)?2:0);
  
  
  
  HT=(model<=1?((5.373e+01)+(1.296e+01)*age-(5.506e-01)*pow(age,2)+(1.113e-02)*pow(age,3)-(1.106e-04)*pow(age,4)+(4.697e-07)*pow(age,5)-(4.416e-10)*pow(age,6)):
        model<=2?((5.869e+01)+(1.265e+01)*age-(4.665e-01)*pow(age,2)+(7.198e-03)*pow(age,3)-(3.224e-05)*pow(age,4)-(2.512e-07)*pow(age,5)+(2.071e-09)*pow(age,6)):0);
  
  BW=(model<=1?((2.354e+00)+(4.050e+00)*age-(3.240e-02)*pow(age,2)-(3.057e-03)*pow(age,3)+(9.353e-05)*pow(age,4)-(1.022e-06)*pow(age,5)+(3.918e-09)*pow(age,6)):
        model<=2?((3.382e+00)+(2.866e+00)*age+(1.694e-01)*pow(age,2)-(1.169e-02)*pow(age,3)+(2.577e-04)*pow(age,4)-(2.484e-06)*pow(age,5)+(8.891e-09)*pow(age,6)):0);
  
  SA=(0.007184*pow(BW,0.425)*pow(HT,0.725));
  
  vliver=(model<=1?(0.0017717-0.0030113*age+0.0253455*BW):model<=2?(-0.0143744-0.0044728*age+0.0264591*BW):0);
  

  vkidney=(model<=1?(0.0458676-0.0003957*age+0.0035115*BW):model<=2?(5.668e-02-4.962e-04*age+3.501e-03*BW):0);
 
  vfilterate =(Ffilterate * BW);
  
  vfat=(model<=1?(6.132e-01+(8.475e-02*age)+(8.151e-05*pow(age,2)) +(1.341e-01*BW)+(2.297e-03*pow(BW,2))):
          model<=2?(1.3054356+(0.3622685*age)-(0.0025165*pow(age,2)) +(0.0906119*BW)+(0.0001731*pow(BW,2))):0);
  
  Fbm_m=(age<=10?0.015:age<=90?0.05:0);
  Fbm_f=(age<=10?0.015:age<=90?0.05:0);
  vbm=(model<=1?(Fbm_f*BW):model<=2?(Fbm_m*BW):0);
  
  vgut=(model<=1?(Fgut_f*BW):model<=2?(Fgut_m*BW):0);
  
  vplasma=(model<=1?((1423*SA-194)/1000):model<=2?((1587*SA-304)/1000):0);
  
  vskin=(model<=1?(2.995e-01-(1.680e-02*age)+(7.151e-05*pow(age,2))+ (5.456e-02*BW)-(3.793e-04*pow(BW,2))):
           model<=2?(3.796e-01-(4.055e-02*age)+(2.564e-04*pow(age,2))+ (4.671e-02*BW)-(8.207e-05*pow(BW,2))):0);
  
  vrestbody=(model<=1?(BW-(vliver+vkidney+vfilterate+vfat+vgut+vplasma)):
               model<=2?(BW-(vliver+vkidney+vfilterate+vfat+vgut+vplasma)):0);
  
 
  ###############################################################################################################################################################################
  #Blood flow
  
  QCC =(model<=1?(5.528076-(2.834486*age)+(0.012591*pow(age,2))+(204.262351*SA)+(19.274290*pow(SA,2))):
          model<=2?(6.48370-1.59948*age+214.68572*SA):0);              																		#Total cardiac output in L/hr
  
  vhct=((0.00000112815*pow(age,3))-(0.000172362*pow(age,2))+(0.00851264*age)+0.327363);
  
  
  QCplasma=QCC;
  
  Qliver=(model<=1?((2.590e-01-1.042e-03*age+4.265e-04*BW)*QCplasma):model<=2?((2.502e-01-1.062e-03*age+ 3.439e-04 *BW)*QCplasma):0);   
  
  Qkidney=(model<=1?((1.819e-01-1.160e-03*age+ 4.873e-04*BW)*QCplasma):model<=2?((1.975e-01 -1.182e-03*age+ 3.937e-04*BW)*QCplasma):0);
  
  ## multiply with 0.06 to convert the ml/min to L/hr
  # also created GFR based on CKD (chronic kidney disease, if there any then define GFR rate in the simulaltion file)
  
  Qfilterate = (CKD == 0?(((-6.616*pow(SA,2)) + (99.054*SA) - 17.74)*0.06):GFR*0.06);   
  
  
  Qfat=(model<=1?((8.298e-02+6.850e-04*age-2.804e-04*BW)*QCplasma):model<=2?((5.075e-02+4.327e-04*age -1.401e-04 *BW)*QCplasma):0); 
  
 
  Qgut=(model<=1?(FQgut_f*QCplasma):model<=2?(FQgut_m*QCplasma):0);
  
  
  Qrestbody=(model<=1?(QCplasma-(Qliver+Qkidney+Qfat+Qgut)):
               model<=2?(QCplasma-(Qliver+Qkidney+Qfat+Qgut)):0);
  
  
  Qurine = QurineC*pow((BWrabbit/BW),scaling);
  Tm = Tmc*pow((BWrabbit/BW),scaling);  #mg/L/kg to whole body weight (mcg/ml to mg/l is equivalent)
  Vmax = VmaxC*pow((BWrabbit/BW),scaling)*1000*vliver; #multiplied 1000 for converting per g liver into kg liver then to whole liver by vliver
  Tr = Trc*pow((BWrabbit/BW),scaling);  #mg/L/kg to whole body weight (mcg/ml to mg/l is equivalent)
  
  Vehr = Vehrc*pow((BWrat/BW),scaling);
  kbileS = kbile*pow((BWrat/BW),scaling);
  
  kgutabsS = kgutabs*pow((BWrat/BW),scaling);
  
  #############
  
  cgut = Agut/vgut;
  cliver = Aliver/vliver;
  ckidney = Akidney/vkidney;
  cfilterate=Afilterate/vfilterate;
  cfat = Afat/vfat;
  crestbody = Arestbody/vrestbody;
  cplasma = Aplasma/vplasma;
  
  Filteration = Qfilterate*ckidney*((fu/Kbp)/Kkidney_plasma);
  Reabsorption = (Tr*cfilterate);   #/(Kr+cfilterate)
  Secretion = (Tm*(ckidney*((fu/Kbp)/Kkidney_plasma)))/(Kt+(ckidney*((fu/Kbp)/Kkidney_plasma)));
  Absorption =  kgutabsS*Agutlumen;
  #############
  
  DOSE=((age<=12)?(oraldose1*BW):(oraldose1));
  
  #dt(Agutlumen) = ((oraldose1*BW))/expT- Absorption- kfeces*Agutlumen +kbileS*ABile;   
  
  dt(Agutlumen) = (DOSE/expT)- Absorption- kfeces*Agutlumen +kbileS*ABile;
  
  dt(Agut) = Absorption + Qgut*(cplasma*(fu/Kbp)- cgut*((fu/Kbp)/Kgut_plasma)); 
  
  dt(Aliver) = Qliver*cplasma*(fu/Kbp) + Qgut*cgut*((fu/Kbp)/Kgut_plasma) - (Qliver +Qgut)*cliver*((fu/Kbp)/Kliver_plasma)-(Vmax*cliver*fu)/(Km + cliver*fu) - (Vehr*cliver*fu)/(Kehr + cliver*fu);
  
  dt(ALiverdeg) = (Vmax*cliver*(fu))/(Km + cliver*(fu));
  
  dt(ABile) =  (Vehr*cliver*(fu))/(Kehr + cliver*(fu)) - kbileS*ABile;
  
  dt(Afeces) = kfeces*Agutlumen;
  
  dt(Akidney) = Qkidney*cplasma*(fu/Kbp) - Qkidney*ckidney*((fu/Kbp)/Kkidney_plasma) -Filteration + Reabsorption -Secretion;
  
  dt(Afilterate) = Filteration - Qfilterate*cfilterate -Reabsorption + Secretion ;  
  
  dt(Adelay)  = Qfilterate*cfilterate  -Qurine *Adelay;
  
  dt(Afat) = Qfat*cplasma*(fu/Kbp) - Qfat*cfat*((fu/Kbp)/Kfat_plasma);  #Filtrate_amount
  dt(Arestbody) = Qrestbody*cplasma*(fu/Kbp) - Qrestbody*crestbody*((fu/Kbp)/Krestbody_plasma);
  
  dt(Aplasma) = (Qliver +Qgut)*cliver*((fu/Kbp)/Kliver_plasma)+ Qkidney*ckidney*((fu/Kbp)/Kkidney_plasma) + 
    Qfat*cfat*((fu/Kbp)/Kfat_plasma)+ Qrestbody*crestbody*((fu/Kbp)/Krestbody_plasma)- QCplasma*cplasma*(fu/Kbp); 
  
  dt(Aurine)  = Qurine *Adelay;    #Qurine is L/hr
  
  # age below 12 dose are given per body weight
  #percentage fraction output of chemical in urine
  P_excreted = ((age<=12)?(Aurine/((IVdosing + oraldose*BW)))*100:(Aurine/((IVdosing + oraldose)))*100);
  B_excreted = ((age<=12)?(ABile/(IVdosing + ((oraldose*BW))))*100:(ABile/((IVdosing + oraldose)))*100) ;
  
}

CalcOutputs {  
  
  Total_renal_Cl = Filteration + Secretion - Reabsorption;
  
  GFR_contribution = (t>0?(Filteration/Total_renal_Cl)*100 : 0);
  Re_GFR_contribution = (t>0?((Filteration-Reabsorption)/Total_renal_Cl)*100 : 0);
  Secr_contribution = (t>0?(Secretion/Total_renal_Cl)*100: 0);
  Reabs_contribution = (t>0?(Reabsorption/Total_renal_Cl)*100: 0);
  
  
  Massbalance =  Agutlumen + Agut + Aliver + ABile + ALiverdeg +  Afeces + Akidney + Afilterate + Adelay + Afat + Arestbody + Aplasma + Aurine;
  
} # End of output calculationn

End.