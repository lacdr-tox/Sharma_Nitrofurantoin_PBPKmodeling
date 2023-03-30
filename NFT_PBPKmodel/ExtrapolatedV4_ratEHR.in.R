MCMC ("MCMC.default.out","", # name of output and restart file
      "",           # name of data file
      100000,0,       # iterations, print predictions flag,
      1,10000,       # printing frequency, iters to print
      10101010);    # random seed (default )   #acetylated form switched to nonaceylated


Level { # top   priors on parameters for all chemicals
  
  
  #Distrib(Kgut_plasma, TruncNormal,2, 1, 1e-3, 10);                # exact mean
  # Distrib(Kliver_plasma, TruncNormal, 2, 0.1, 1e-2, 3);           # exact mean
  #Distrib(Kkidney_plasma, TruncNormal,  4.347, 4,  1e-3, 20);
  # # # Distrib(Kfat_plasma , TruncNormal,  0.77,  3,  1e-3, 10);             # exact mean
  #Distrib(Krestbody_plasma, TruncNormal,  0.010, 0.1,  1e-3, 10);

  Distrib(Vehrc, TruncNormal, 0.6, 0.1,  1e-4, 200);
  #Distrib(Kehr, TruncNormal, 10, 10,  1e-4, 800);
  
  Distrib(Kehr, TruncNormal, 0.09, 0.01,  1e-4, 800);
  
  Distrib(kfeces, TruncNormal, 0.05, 0.01,  1e-4, 8);
 
  Distrib(kbile, TruncNormal, 0.02,  0.08,  1e-5, 3);   #kst
  
  
 Distrib(kgutabs, TruncNormal, 2,  1,  1e-3, 10);

  
  Likelihood (Data(cplasma), Normal, Prediction(cplasma), 1.1);
  Likelihood (Data(B_excreted), Normal, Prediction(B_excreted), 1.1);
  Likelihood (Data(P_excreted), Normal, Prediction(P_excreted), 1.1);
  
  
  
  Level {    # NFT in rats
    
#   
   
    
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
    
  #   Simulation {    # Dose 4   #dose = 2.5  (from different study of Paul et al. 1959)
  # 
  #     IVdosing = 0;
  #     oraldose = 3.5;
  #     expT = 0.001;
  #     oraldose1 = PerDose (oraldose,24,0, expT);
  # 
  #     Print (cplasma,0.83);
  #     Data  (cplasma, -1);
  # 
  #     Print(B_excreted, 0.25);
  #     Data  (B_excreted,-1);
  # 
  #     Print(P_excreted, 4);
  #     Data(P_excreted,24);
  #   }
  # 
  #   Simulation {    # Dose 1 # dose 0.5
  # 
  #     IVdosing = 0;
  #     expT = 0.001;
  #     oraldose = 10;
  #     oraldose1 = PerDose (oraldose,24,0, expT);
  # 
  #     Print (cplasma,0.08,0.25,0.5,1,2,4,6,8,12);
  #     Data  (cplasma,0.25,0.288,0.433,0.584,0.998,0.573,0.353,0.231,0.14);
  # 
  #     Print(B_excreted,0.5);
  #     Data  (B_excreted,-1);
  #     #
  #     # Print(P_excreted,4);
  #     # Data(P_excreted,-1);     # this urine data from other study with the same dose of 10 mg/kg Paul et al. 1959
  # 
  #     Print(P_excreted,4);
  #     Data(P_excreted,40);     # this urine data from other study with the same dose of 10 mg/kg Paul et al. 1959
  # 
  # 
  #   }
  #   Simulation {    # Dose 4   #dose = 2.5  (from different study of Paul et al. 1959)
  # 
  #     IVdosing = 0;
  #     oraldose = 25;
  #     expT = 0.001;
  #     oraldose1 = PerDose (oraldose,24,0, expT);
  # 
  #     Print (cplasma,0.83);
  #     Data  (cplasma, -1);
  # 
  #     Print(B_excreted, 0.25);
  #     Data  (B_excreted,-1);
  # 
  #     Print(P_excreted, 4);
  #     Data(P_excreted,29);
  #   }
  # 
  #   
  #   
  #   
   }}

End. 

