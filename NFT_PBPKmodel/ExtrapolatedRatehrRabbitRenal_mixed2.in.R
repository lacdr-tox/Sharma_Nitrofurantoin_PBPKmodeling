
# thsi contain both the data of rats and rabbit to furhter optimize the renal and ehr mechanism. 

MCMC ("MCMC.default.out","", # name of output and restart file
      "",           # name of data file
      300000,0,       # iterations, print predictions flag,
      1,10000,       # printing frequency, iters to print
      10101010);    # random seed (default )   #acetylated form switched to nonaceylated


Level { # top   priors on parameters for all chemicals
  
  
  # Distrib(Kgut_plasma, TruncNormal,2, 1, 1e-3, 10);                # exact mean
  # Distrib(Kliver_plasma, TruncNormal, 2, 1, 1, 3);           # exact mean
  # Distrib(Kkidney_plasma, TruncNormal,  4.347, 4,  1e-3, 20);
  #Distrib(Kfat_plasma , TruncNormal,  0.77,  3,  1e-3, 10);             # exact mean
  #Distrib(Krestbody_plasma, TruncNormal,  0.010, 0.1,  1e-3, 10);
  
  # Distrib(Krestbody_plasma, Uniform,1e-4, 1);
  # Distrib(Kfat_plasma, Uniform,1e-4, 1);

  # Distrib(Vehrc, TruncNormal, 100, 70,  1e-4, 900);
  # Distrib(Kehr, TruncNormal, 500, 100,  1e-4, 900);
  # 
  # Distrib(Kehr, TruncNormal, 0.09, 0.01,  1e-4, 800);
  # 
  # Distrib(kfeces, TruncNormal, 0.05, 0.01,  1e-4, 8);

  #Distrib(VmaxC, TruncNormal, 4, 3,  0.01, 100);
  
  # Distrib(VmaxC, Uniform, 1e-4, 200);
  # Distrib(Km, TruncNormal, 200, 100,  1e-4, 800);
  
  # Distrib(fu, Uniform, 0.25, 0.90);
  
  Distrib(Vehrc, TruncNormal, 4, 1,  1e-3, 2);
 
  #Distrib(Kehr, TruncNormal, 0.5, 1,  1e-3, 10);
  
  Distrib(Kehr, Uniform,3e-4, 20);
  
  #Distrib(kbile, TruncNormal, 0.5, 1,  1e-2, 20);
  
  # #Distrib(Kehr, TruncNormal, 10, 10,  1e-4, 800);
  
  Distrib(kfeces, Uniform,1e-5, 8);
  
  Distrib(kbile, Uniform,1e-5, 10);   #kst
  
  Distrib(kgutabs, TruncNormal, 0.5, 1,  1e-4, 20);
  
  #Distrib(kgutabs, Uniform, 1e-3, 100);
  
  Level {   # priors for individual chemicals (on fraction acetylated form)

  Distrib(Vehrc, TruncNormal, Vehrc, 0.5,  1e-3, 0.7);
    #Distrib(Kehr, TruncNormal, 10, 10,  1e-4, 800);
    
  Distrib(Kehr, TruncNormal,Kehr, 2, 1e-3, 10);
    
  Distrib(kfeces, TruncNormal, kfeces, 3,  1e-4, 8);
    
  Distrib(kbile, TruncNormal, kbile,  3,  1e-4, 20);   #kst
    
  Distrib(kgutabs, TruncNormal,kgutabs, 0.5,  1e-3, 3);
    
  # Distrib(fu, TruncNormal,fu, 0.2, 0.25, 0.90);
    
    
  Likelihood (Data(cplasma), Normal, Prediction(cplasma), 1.2);
  Likelihood (Data(B_excreted), Normal, Prediction(B_excreted), 1.5);
  Likelihood (Data(P_excreted), Normal, Prediction(P_excreted), 1.5);
  
  
  
  Level {    # NFT in rats
    
    
    Simulation {    # Dose 3   #dose = 2.5
      
      IVdosing = 1.5;
      oraldose = 0;
      expT = 0.001;
      # oraldose1 = PerDose (oraldose,24,0, expT);
      oraldose1 = 0;
      Print (cplasma,0.04,0.1);
      Data  (cplasma, -1, -1);
      
      Print(B_excreted, 0.25,0.5, 0.75,1,1.25,1.5,1.75,2);
      Data  (B_excreted,8.6,12.63,15.2,17,18,19,19.53,20);   
      Print(P_excreted, 0.25);
      Data(P_excreted,-1);
    }
    
    Simulation {    # Dose 2   #dose = 1.25
      
      IVdosing = 2;
      oraldose = 0;
      expT = 0.001;
      # oraldose1 = PerDose (oraldose,24,0, expT);
      oraldose1 = 0;
      Print (cplasma,0.04,0.1,0.25,0.5,1,2);
      Data  (cplasma,8.18,6.12,2.88,0.7,0.19,0.05);
      
      Print(B_excreted, 0.5);
      Data  (B_excreted,-1,);
      
      Print(P_excreted, 0.5);
      Data(P_excreted,-1);
      
    }
    
    

    Simulation {    # Dose 4   #dose = 2.5  (from different study of Paul et al. 1959)

      IVdosing = 3;    #only data of percentage dose recovery in urine.
      oraldose = 0;
      expT = 0.001;
      #oraldose1 = PerDose (oraldose,24,0, expT);
      oraldose1 = 0;

      Print (cplasma,0.98);
      Data  (cplasma, -1);

      Print(B_excreted, 0.25);
      Data  (B_excreted,-1);

      Print(P_excreted, 4);
      Data(P_excreted,37);
    }
    
    
    Simulation {    # Dose 4   #dose = 2.5  (from different study of Paul et al. 1959)

      IVdosing = 5;
      oraldose = 0;
      expT = 0.001;
      oraldose1 = PerDose (oraldose,24,0, expT);

      Print (cplasma,0.98);
      Data  (cplasma, 1);

      Print(B_excreted, 0.25);
      Data  (B_excreted,-1);

      Print(P_excreted, 0.25);
      Data(P_excreted,-1);
    }
    
    Simulation {    # Dose 4   #dose = 2.5  (from different study of Paul et al. 1959)

      IVdosing = 7;
      oraldose = 0;
      expT = 0.001;
      oraldose1 = PerDose (oraldose,24,0, expT);

      Print (cplasma,0.83,0.9,0.98);
      Data  (cplasma, 3.7,3.7,3.4);

      Print(B_excreted, 0.25);
      Data  (B_excreted,-1);

      Print(P_excreted, 0.25);
      Data(P_excreted,-1);
    }
    
    Simulation {    # Dose 4   #dose = 2.5  (from different study of Paul et al. 1959)

      IVdosing = 10;
      oraldose = 0;
      expT = 0.001;
      oraldose1 = PerDose(oraldose,24,0, expT);

      Print (cplasma,0.98);
      Data  (cplasma, 7.3);

      Print(B_excreted, 0.25);
      Data  (B_excreted,-1);

      Print(P_excreted, 4);
      Data(P_excreted,35);
    }
    
    Simulation {    # Dose 4   #dose = 2.5  (from different study of Paul et al. 1959)

      IVdosing = 17;
      oraldose = 0;
      expT = 0.001;
      oraldose1 = PerDose (oraldose,24,0, expT);

      Print (cplasma,0.98);
      Data  (cplasma, 13.1);

      Print(B_excreted, 0.25);
      Data  (B_excreted,-1);

      Print(P_excreted, 4);
      Data(P_excreted,-1);
    }

    Simulation {    # Dose 4   #dose = 2.5  (from different study of Paul et al. 1959)

      IVdosing = 25;
      oraldose = 0;
      expT = 0.001;
      oraldose1 = PerDose (oraldose,24,0, expT);

      Print (cplasma,0.83,0.9,0.98);
      Data  (cplasma, 21,21.3,24.6);

      Print(B_excreted, 0.25);
      Data  (B_excreted,-1);

      Print(P_excreted, 4);
      Data(P_excreted,24);
    }
    
    #######################
    # oral dosing
    
      Simulation {    # Dose 4   #dose = 2.5  (from different study of Paul et al. 1959)

        IVdosing = 0;
        oraldose = 3.5;
        expT = 0.001;
        oraldose1 = PerDose (oraldose,24,0, expT);

        Print (cplasma,0.83);
        Data  (cplasma, -1);

        Print(B_excreted, 0.25);
        Data  (B_excreted,-1);

        Print(P_excreted, 4);
        Data(P_excreted,24);
      }

      Simulation {    # Dose 1 # dose 0.5

        IVdosing = 0;
        expT = 0.001;
        oraldose = 10;
        oraldose1 = PerDose (oraldose,24,0, expT);

        Print (cplasma,0.08,0.25,0.5,1,2,4,6,8,12);
        Data  (cplasma,0.25,0.288,0.433,0.584,0.998,0.573,0.353,0.231,0.14);

        Print(B_excreted,0.5);
        Data  (B_excreted,-1);
        #
        # Print(P_excreted,4);
        # Data(P_excreted,-1);     # this urine data from other study with the same dose of 10 mg/kg Paul et al. 1959

        Print(P_excreted,4);
        Data(P_excreted,40);     # this urine data from other study with the same dose of 10 mg/kg Paul et al. 1959


      }
  }
      Level {    # NFT in rats     
    #######################################################################
    #THis the rabbit IV DATA OF URINE EXCRETION
    ####################################################################3
    Simulation {    # Dose 1 # dose 0.5
      BW = 2.5; 
      QCC = 15.96;                  #Total Cardiac blood output (L/h/kg^0.75) not it is to the power 0.75 
      HCT = 0.41;  				         #hematocrit percentage
      FQliver = 0.1245; 		         #Fraction cardiac output going to liver 
      FQkidney = 0.151;  		       #Fraction cardiac output going to kidney  
      FQfat =  0.0604;  		         #Fraction cardiac output going to fat  
      FQgut = 0.2094;   			       #Fraction cardiac output going to gut
      
      Qfilterate = 0.631;      # Rabbit GFR rate
      
      
      #constant organ volume as a fraction of total body weight  # For rabit from hTTK
      Fliver = 0.04; 		           #Fraction liver volume
      Fkidney = 0.006; 			       #Fraction kidney volume 
      Ffilterate=0.0006;           # 10 percent of kidney volume
      Ffat = 0.048;                 #fractional volume of fat  
      Fgut = 0.048;   				       #fractional volume of gut
      Fplasma = 0.08594;  	         #fractional volume of plasma
      
      IVdosing = 0.5;
      oraldose1 = 0;
      
      Print (cplasma,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.43,0.5,0.55);
      Data  (cplasma,0.92,0.71,0.62,0.52,0.38,0.35,0.28,0.22,0.16,0.14,0.11);
      
      Print(P_excreted, 0.5,0.75,1,1.5,2,3,4);
      Data  (P_excreted,33,47,54,57,60,61,62);
      
      Print(B_excreted, 0.25);
      Data  (B_excreted,-1);
      
    }
    
    Simulation {    # Dose 2   #dose = 1.25
      
      BW = 2.5; 
      QCC = 15.96;                  #Total Cardiac blood output (L/h/kg^0.75) not it is to the power 0.75 
      HCT = 0.41;  				         #hematocrit percentage
      FQliver = 0.1245; 		         #Fraction cardiac output going to liver 
      FQkidney = 0.151;  		       #Fraction cardiac output going to kidney  
      FQfat =  0.0604;  		         #Fraction cardiac output going to fat  
      FQgut = 0.2094;   			       #Fraction cardiac output going to gut
      Qfilterate = 0.631;      # Rabbit GFR rate
      
      
      #constant organ volume as a fraction of total body weight  # For rabit from hTTK
      Fliver = 0.04; 		           #Fraction liver volume
                   #fraction lung volume  		         		                    
      Fkidney = 0.006; 			       #Fraction kidney volume 
      Ffilterate=0.0006;           # 10 percent of kidney volume
      Ffat = 0.048;                 #fractional volume of fat  
      Fgut = 0.048;   				       #fractional volume of gut
      Fplasma = 0.08594;  	         #fractional volume of pla
      
      IVdosing = 1.25;
      oraldose1 = 0;
      
      Print (cplasma,0.08,0.16,0.25,0.33,0.42,0.5,0.58, 0.66, 0.75,0.83,0.91);
      Data  (cplasma,2.27,1.47,1.21,1.00,0.77,0.55,0.42,0.30,0.24,0.19,0.15);
      
      Print(P_excreted, 0.5,0.75,1,1.5,2);
      Data  (P_excreted,-1,-1,-1,-1,-1);
      Print(B_excreted, 0.25);
      Data  (B_excreted,-1);
      
    }
    
    
    Simulation {    # Dose 3   #dose = 2.5
      
      BW = 2.5; 
      QCC = 15.96;                  #Total Cardiac blood output (L/h/kg^0.75) not it is to the power 0.75 
      HCT = 0.41;  				         #hematocrit percentage
      FQliver = 0.1245; 		         #Fraction cardiac output going to liver 
      FQkidney = 0.151;  		       #Fraction cardiac output going to kidney  
      FQfat =  0.0604;  		         #Fraction cardiac output going to fat  
      FQgut = 0.2094;   			       #Fraction cardiac output going to gut
      Qfilterate = 0.631;      # Rabbit GFR rate
      
      
      #constant organ volume as a fraction of total body weight  # For rabit from hTTK
      Fliver = 0.04; 		           #Fraction liver volume
                   #fraction lung volume  		         		                    
      Fkidney = 0.006; 			       #Fraction kidney volume 
      Ffilterate=0.0006;           # 10 percent of kidney volume
      Ffat = 0.048;                 #fractional volume of fat  
      Fgut = 0.048;   				       #fractional volume of gut
      Fplasma = 0.08594;  	         #fractional volume of pla
      
      IVdosing = 2.5;
      oraldose1 = 0;
      
      Print (cplasma,0.1,0.2,0.26,0.5,0.75,1,1.25,1.5);
      Data  (cplasma,4.45,3.39,2.69,2.23,1.40,0.61,0.26,0.13);
      
      Print(P_excreted, 0.5,0.75,1,1.5,2);
      Data  (P_excreted,28,37,42,44,45);
      Print(B_excreted, 0.25);
      Data  (B_excreted,-1);
    }
    
    Simulation {    # Dose 4   #dose = 5
      
      BW = 2.5; 
      QCC = 15.96;                  #Total Cardiac blood output (L/h/kg^0.75) not it is to the power 0.75 
      HCT = 0.41;  				         #hematocrit percentage
      FQliver = 0.1245; 		         #Fraction cardiac output going to liver 
      FQkidney = 0.151;  		       #Fraction cardiac output going to kidney  
      FQfat =  0.0604;  		         #Fraction cardiac output going to fat  
      FQgut = 0.2094;   			       #Fraction cardiac output going to gut
      Qfilterate = 0.631;      # Rabbit GFR rate
      
      #constant organ volume as a fraction of total body weight  # For rabit from hTTK
      Fliver = 0.04; 		           #Fraction liver volume
                   #fraction lung volume  		         		                    
      Fkidney = 0.006; 			       #Fraction kidney volume 
      Ffilterate=0.0006;           # 10 percent of kidney volume
      Ffat = 0.048;                 #fractional volume of fat  
      Fgut = 0.048;   				       #fractional volume of gut
      Fplasma = 0.08594;  	         #fractional volume of pla
      
      IVdosing = 5;
      oraldose1 = 0;
      Print (cplasma,0.1,0.2,0.26,0.5,0.75,1,1.25,1.5);
      Data  (cplasma,11.91,8.99,6.92,3.46,1.77,0.87,0.44,0.22);
      
      Print(P_excreted,0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4);
      Data(P_excreted,9.71,21,27,31,34,36,37, 38);
      Print(B_excreted, 0.25);
      Data  (B_excreted,-1);
      
    }
    Simulation {    # Dose 4   #dose = 5
      BW = 2.5; 
      QCC = 15.96;                  #Total Cardiac blood output (L/h/kg^0.75) not it is to the power 0.75 
      HCT = 0.41;  				         #hematocrit percentage
      FQliver = 0.1245; 		         #Fraction cardiac output going to liver 
      FQkidney = 0.151;  		       #Fraction cardiac output going to kidney  
      FQfat =  0.0604;  		         #Fraction cardiac output going to fat  
      FQgut = 0.2094;   			       #Fraction cardiac output going to gut
      Qfilterate = 0.631;      # Rabbit GFR rate
      
      #constant organ volume as a fraction of total body weight  # For rabit from hTTK
      Fliver = 0.04; 		           #Fraction liver volume
                   #fraction lung volume  		         		                    
      Fkidney = 0.006; 			       #Fraction kidney volume 
      Ffilterate=0.0006;           # 10 percent of kidney volume
      Ffat = 0.048;                 #fractional volume of fat  
      Fgut = 0.048;   				       #fractional volume of gut
      Fplasma = 0.08594;  	         #fractional volume of plasma
      
      IVdosing = 5;
      oraldose1 = 0;
      Print (cplasma,0.1,0.2,0.26,0.5,0.75,1,1.25,1.5);
      Data  (cplasma,11.91,8.99,6.92,3.46,1.77,0.87,0.44,0.22);
      
      Print(P_excreted,0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4);
      Data(P_excreted,9.71,21,27,31,34,36,37, 38);
      
      Print(B_excreted, 0.25);
      Data  (B_excreted,-1);
    }
  
    Simulation {    # Dose 5    #dose = 10
      
      BW = 2.5; 
      QCC = 15.96;                  #Total Cardiac blood output (L/h/kg^0.75) not it is to the power 0.75 
      HCT = 0.41;  				         #hematocrit percentage
      FQliver = 0.1245; 		         #Fraction cardiac output going to liver 
      FQkidney = 0.151;  		       #Fraction cardiac output going to kidney  
      FQfat =  0.0604;  		         #Fraction cardiac output going to fat  
      FQgut = 0.2094;   			       #Fraction cardiac output going to gut
      Qfilterate = 0.631;      # Rabbit GFR rate
      
      #constant organ volume as a fraction of total body weight  # For rabit from hTTK
      Fliver = 0.04; 		           #Fraction liver volume
                   #fraction lung volume  		         		                    
      Fkidney = 0.006; 			       #Fraction kidney volume 
      Ffilterate=0.0006;           # 10 percent of kidney volume
      Ffat = 0.048;                 #fractional volume of fat  
      Fgut = 0.048;   				       #fractional volume of gut
      Fplasma = 0.08594;  	         #fractional volume of plasma
      IVdosing = 10;
      oraldose1 = 0;
      Print (cplasma,0.1,0.2,0.26,0.5,0.75,1,1.25,1.5,1.75,2,2.25);
      Data (cplasma,23.34,19.09,15.46,10.14,6.08,3.53,2.23,1.39,0.59,0.31,0.15);
      
      Print(P_excreted,0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4);
      Data(P_excreted,7.11,7.5,12.5,15,18,19.65,20.37, 20.23);
      
      Print(B_excreted, 0.25);
      Data  (B_excreted,-1);
      
    }
    
    Simulation {    # Dose 6    #dose = 15
      
      BW = 2.5; 
      QCC = 15.96;                  #Total Cardiac blood output (L/h/kg^0.75) not it is to the power 0.75 
      HCT = 0.41;  				         #hematocrit percentage
      FQliver = 0.1245; 		         #Fraction cardiac output going to liver 
      FQkidney = 0.151;  		       #Fraction cardiac output going to kidney  
      FQfat =  0.0604;  		         #Fraction cardiac output going to fat  
      FQgut = 0.2094;   			       #Fraction cardiac output going to gut
      Qfilterate = 0.631;      # Rabbit GFR rate
      
      #constant organ volume as a fraction of total body weight  # For rabit from hTTK
      Fliver = 0.04; 		           #Fraction liver volume
                   #fraction lung volume  		         		                    
      Fkidney = 0.006; 			       #Fraction kidney volume 
      Ffilterate=0.0006;           # 10 percent of kidney volume
      Ffat = 0.048;                 #fractional volume of fat  
      Fgut = 0.048;   				       #fractional volume of gut
      Fplasma = 0.08594;  	         #fractional volume of plasma
      IVdosing = 15;
      oraldose1 = 0;
      Print (cplasma,0.1,0.2,0.26,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3);
      Data (cplasma,37.79,31.23,27.41,18.90,12.40,8.05,5.66,3.91,2.59,1.57,1.04,0.68,0.45,0.28);
      
      Print(P_excreted,0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4);
      Data(P_excreted,1.2,2.5,4.7,6.8,10.7,10.8,11.72,11.72);
      
      Print(B_excreted, 0.25);
      Data  (B_excreted,-1);
    }
    
    
    Simulation {    # Dose 1 # dose 0.5
      BW = 2.5; 
      QCC = 15.96;                  #Total Cardiac blood output (L/h/kg^0.75) not it is to the power 0.75 
      HCT = 0.41;  				         #hematocrit percentage
      FQliver = 0.1245; 		         #Fraction cardiac output going to liver 
      FQkidney = 0.151;  		       #Fraction cardiac output going to kidney  
      FQfat =  0.0604;  		         #Fraction cardiac output going to fat  
      FQgut = 0.2094;   			       #Fraction cardiac output going to gut
      Qfilterate = 0.631;      # Rabbit GFR rate
      
      #constant organ volume as a fraction of total body weight  # For rabit from hTTK
      Fliver = 0.04; 		           #Fraction liver volume
                   #fraction lung volume  		         		                    
      Fkidney = 0.006; 			       #Fraction kidney volume 
      Ffilterate=0.0006;           # 10 percent of kidney volume
      Ffat = 0.048;                 #fractional volume of fat  
      Fgut = 0.048;   				       #fractional volume of gut
      Fplasma = 0.08594;  	         #fractional volume of plasma
      IVdosing = 0;
      expT = 0.001;
      oraldose = 1.25;
      oraldose1 = PerDose (oraldose,24,0, expT);
      
      Print (cplasma,0.165,0.33,0.5,0.66,0.83,1,1.5,2,2.5, 3);
      Data  (cplasma,0.015,0.2,0.24,0.26,0.23,0.21,0.15,0.1,0.03,0.01);
      Print(B_excreted, 0.25);
      Data  (B_excreted,-1);
      Print(P_excreted, 0.5);
      Data  (P_excreted,-1);
      
    }
    
    Simulation {    # Dose 1 # dose 0.5
      
      BW = 2.5; 
      QCC = 15.96;                  #Total Cardiac blood output (L/h/kg^0.75) not it is to the power 0.75 
      HCT = 0.41;  				         #hematocrit percentage
      FQliver = 0.1245; 		         #Fraction cardiac output going to liver 
      FQkidney = 0.151;  		       #Fraction cardiac output going to kidney  
      FQfat =  0.0604;  		         #Fraction cardiac output going to fat  
      FQgut = 0.2094;   			       #Fraction cardiac output going to gut
      Qfilterate = 0.631;      # Rabbit GFR rate
      
      #constant organ volume as a fraction of total body weight  # For rabit from hTTK
      Fliver = 0.04; 		           #Fraction liver volume
                   #fraction lung volume  		         		                    
      Fkidney = 0.006; 			       #Fraction kidney volume 
      Ffilterate=0.0006;           # 10 percent of kidney volume
      Ffat = 0.048;                 #fractional volume of fat  
      Fgut = 0.048;   				       #fractional volume of gut
      Fplasma = 0.08594;  	         #fractional volume of plasma
      IVdosing = 0;
      expT = 0.001;
      oraldose = 10;
      oraldose1 = PerDose (oraldose,24,0, expT);
      
      Print (cplasma,0.083,0.165,0.33,0.5,0.66,0.83,1,1.5,2,2.5, 3);
      Data  (cplasma,0.17,0.77,1.41,1.78,1.27,1.15,1.06,0.87,0.754,0.528,0.27);
      Print(B_excreted, 0.25);
      Data  (B_excreted,-1);
      Print(P_excreted, 0.5);
      Data  (P_excreted,-1);
      
      
    }}
  }
}

End.
  ##############################################################
    #oral dose data
    
    
    # 
    # Simulation {    # Dose 4   #dose = 2.5  (from different study of Paul et al. 1959)
    # 
    #   IVdosing = 0;
    #   oraldose = 3.5;
    #   expT = 0.001;
    #   oraldose1 = PerDose (oraldose,24,0, expT);
    # 
    #   Print (cplasma,0.83);
    #   Data  (cplasma, -1);
    # 
    #   Print(B_excreted, 0.25);
    #   Data  (B_excreted,-1);
    # 
    #   Print(P_excreted, 4);
    #   Data(P_excreted,24);
    # }
    # 
    # Simulation {    # Dose 1 # dose 0.5
    # 
    #   IVdosing = 0;
    #   expT = 0.001;
    #   oraldose = 10;
    #   oraldose1 = PerDose (oraldose,24,0, expT);
    # 
    #   Print (cplasma,0.08,0.25,0.5,1,2,4,6,8,12);
    #   Data  (cplasma,0.25,0.288,0.433,0.584,0.998,0.573,0.353,0.231,0.14);
    # 
    #   Print(B_excreted,0.5);
    #   Data  (B_excreted,-1);
    #   #
    #   # Print(P_excreted,4);
    #   # Data(P_excreted,-1);     # this urine data from other study with the same dose of 10 mg/kg Paul et al. 1959
    # 
    #   Print(P_excreted,4);
    #   Data(P_excreted,40);     # this urine data from other study with the same dose of 10 mg/kg Paul et al. 1959
    # 
    # 
    # }
    # Simulation {    # Dose 4   #dose = 2.5  (from different study of Paul et al. 1959)
    # 
    #   IVdosing = 0;
    #   oraldose = 25;
    #   expT = 0.001;
    #   oraldose1 = PerDose (oraldose,24,0, expT);
    # 
    #   Print (cplasma,0.83);
    #   Data  (cplasma, -1);
    # 
    #   Print(B_excreted, 0.25);
    #   Data  (B_excreted,-1);
    # 
    #   Print(P_excreted, 4);
    #   Data(P_excreted,29);
    # }
    # 
    # 
    # ################################################################################################
    # #Rabbit oral dose
    # 
    # Simulation {    # Dose 1 # dose 0.5
    # 
    #   BW = 2.5; 
    #   QCC = 15.96;                  #Total Cardiac blood output (L/h/kg^0.75) not it is to the power 0.75 
    #   HCT = 0.41;  				         #hematocrit percentage
    #   FQliver = 0.1245; 		         #Fraction cardiac output going to liver 
    #   FQkidney = 0.151;  		       #Fraction cardiac output going to kidney  
    #   FQfat =  0.0604;  		         #Fraction cardiac output going to fat  
    #   FQgut = 0.2094;   			       #Fraction cardiac output going to gut
    #   #FQfilterate=0.1872;           #Fraction of glomeruler filteration L/hr/BW^0.75
    #   FQfilterate=0.0302;           #Fraction of glomeruler filteration L/hr/BW^0.75 (20 percent of kidney blood flow)
    #   
    #   
    #   #constant organ volume as a fraction of total body weight  # For rabit from hTTK
    #   Fliver = 0.04; 		           #Fraction liver volume
    #                #fraction lung volume  		         		                    
    #   Fkidney = 0.006; 			       #Fraction kidney volume 
    #   Ffilterate=0.0006;           # 10 percent of kidney volume
    #   Ffat = 0.048;                 #fractional volume of fat  
    #   Fgut = 0.048;   				       #fractional volume of gut
    #   Fplasma = 0.08594;  	         #fractional volume of plasma
    #   expT = 0.001;
    #   oraldose = 1.25;
    #   oraldose1 = PerDose (oraldose,24,0, expT);
    # 
    #   Print (cplasma,0.165,0.33,0.5,0.66,0.83,1,1.5,2,2.5, 3);
    #   Data  (cplasma,0.015,0.2,0.24,0.26,0.23,0.21,0.15,0.1,0.03,0.01);
    #   # Print (AUC,4);
    #   # Data  (AUC,-1);
    #   Print(P_excreted, 0.5);
    #   Data  (P_excreted,-1);
    #   Print(B_excreted, 0.25);
    #   Data  (B_excreted,-1);
    # 
    # }
    # 
    # Simulation {    # Dose 1 # dose 0.5
    #   
    #   BW = 2.5; 
    #   QCC = 15.96;                  #Total Cardiac blood output (L/h/kg^0.75) not it is to the power 0.75 
    #   HCT = 0.41;  				         #hematocrit percentage
    #   FQliver = 0.1245; 		         #Fraction cardiac output going to liver 
    #   FQkidney = 0.151;  		       #Fraction cardiac output going to kidney  
    #   FQfat =  0.0604;  		         #Fraction cardiac output going to fat  
    #   FQgut = 0.2094;   			       #Fraction cardiac output going to gut
    #   #FQfilterate=0.1872;           #Fraction of glomeruler filteration L/hr/BW^0.75
    #   FQfilterate=0.0302;           #Fraction of glomeruler filteration L/hr/BW^0.75 (20 percent of kidney blood flow)
    #   
    #   
    #   #constant organ volume as a fraction of total body weight  # For rabit from hTTK
    #   Fliver = 0.04; 		           #Fraction liver volume
    #                #fraction lung volume  		         		                    
    #   Fkidney = 0.006; 			       #Fraction kidney volume 
    #   Ffilterate=0.0006;           # 10 percent of kidney volume
    #   Ffat = 0.048;                 #fractional volume of fat  
    #   Fgut = 0.048;   				       #fractional volume of gut
    #   Fplasma = 0.08594;  	         #fractional volume of plasma
    # 
    #   expT = 0.001;
    #   oraldose = 10;
    #   oraldose1 = PerDose (oraldose,24,0, expT);
    # 
    #   Print (cplasma,0.083,0.165,0.33,0.5,0.66,0.83,1,1.5,2,2.5, 3);
    #   Data  (cplasma,0.17,0.77,1.41,1.78,1.27,1.15,1.06,0.87,0.754,0.528,0.27);
    #   # Print (AUC,4);
    #   # Data  (AUC,-1);
    #   Print(P_excreted, 0.5);
    #   Data  (P_excreted,-1);
    #   Print(B_excreted, 0.25);
    #   Data  (B_excreted,-1);
    # 
    # 
    # }
    # 
    
  }}

End. 

