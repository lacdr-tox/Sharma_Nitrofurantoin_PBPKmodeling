MCMC ("MCMC.default.out","", # name of output and restart file
      "",           # name of data file
      100000,0,       # iterations, print predictions flag,
      1,10000,       # printing frequency, iters to print
      10101010);    # random seed (default )   #acetylated form switched to nonaceylated


Level { # top   priors on parameters for all chemicals
  
  
  # Distrib(Kgut_plasma, TruncNormal,1.961, 3, 1e-3, 10);                # exact mean
  # Distrib(Kliver_plasma, TruncNormal, 4.8, 3, 1, 5);           # exact mean
  # Distrib(Kkidney_plasma, TruncNormal,  4.347, 4,  1e-3, 20);
  # Distrib(Kfat_plasma , TruncNormal,  0.77,  3,  1e-3, 10);             # exact mean
   Distrib(Krestbody_plasma, TruncNormal,  0.4,  1,  1e-3, 10);

  Distrib(Tmc , TruncNormal,  2, 1,  1e-2, 500);             # exact mean

  # This is very sensitive and causing chains to diverge  so either fixed it to 0.0008 or see nftrabbitpbpkmy.in.r 
   Distrib(Kt , TruncLogNormal,  0.08, 1.2,  0.0008, 1);            
  

  Distrib(Trc , TruncNormal,  80, 30,  1e-5, 500);             # exact mean

  Distrib(VmaxC, TruncNormal, 0.043, 0.05,  0.0001, 200);
  Distrib(Km, TruncNormal, 60, 10,  1e-4, 800);


  Distrib(QurineC, TruncNormal, 5, 5,  1e-4, 70);
  # 
  Likelihood (Data(cplasma), Normal, Prediction(cplasma), 1.5);
  Likelihood (Data(P_excreted), Normal, Prediction(P_excreted), 1.5);
  
  
  Level {    # Sulforaphane chemicals
    
    
    Simulation {    # Dose 1 # dose 0.5
      
      IVdosing = 0.5;
      oraldose1 = 0;
      
      Print (cplasma,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.43,0.5,0.55);
      Data  (cplasma,0.92,0.71,0.62,0.52,0.38,0.35,0.28,0.22,0.16,0.14,0.11);
      
      Print(P_excreted, 0.5,0.75,1,1.5,2,3,4);
      Data  (P_excreted,33,47,54,57,60,61,62);
      
      
    }
    
    Simulation {    # Dose 2   #dose = 1.25
      
      IVdosing = 1.25;
      oraldose1 = 0;
      
      Print (cplasma,0.08,0.16,0.25,0.33,0.42,0.5,0.58, 0.66, 0.75,0.83,0.91);
      Data  (cplasma,2.27,1.47,1.21,1.00,0.77,0.55,0.42,0.30,0.24,0.19,0.15);
      
      Print(P_excreted, 0.5,0.75,1,1.5,2);
      Data  (P_excreted,-1,-1,-1,-1,-1);
      
    }
    
    
    Simulation {    # Dose 3   #dose = 2.5
      
      IVdosing = 2.5;
      oraldose1 = 0;
      
      Print (cplasma,0.1,0.2,0.26,0.5,0.75,1,1.25,1.5);
      Data  (cplasma,4.45,3.39,2.69,2.23,1.40,0.61,0.26,0.13);
      
      Print(P_excreted, 0.5,0.75,1,1.5,2);
      Data  (P_excreted,28,37,42,44,45);   
    }
    
    Simulation {    # Dose 4   #dose = 5
      
      IVdosing = 5;
      oraldose1 = 0;
      Print (cplasma,0.1,0.2,0.26,0.5,0.75,1,1.25,1.5);
      Data  (cplasma,11.91,8.99,6.92,3.46,1.77,0.87,0.44,0.22);
      
      Print(P_excreted,0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4);
      Data(P_excreted,9.71,21,27,31,34,36,37, 38);
      
    }
    Simulation {    # Dose 5    #dose = 10
      IVdosing = 10;
      oraldose1 = 0;
      Print (cplasma,0.1,0.2,0.26,0.5,0.75,1,1.25,1.5,1.75,2,2.25);
      Data (cplasma,23.34,19.09,15.46,10.14,6.08,3.53,2.23,1.39,0.59,0.31,0.15);
      
      Print(P_excreted,0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4);
      Data(P_excreted,7.11,7.5,12.5,15,18,19.65,20.37, 20.23);
      
    }
    
    Simulation {    # Dose 6    #dose = 15
      IVdosing = 15;
      oraldose1 = 0;
      Print (cplasma,0.1,0.2,0.26,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3);
      Data (cplasma,37.79,31.23,27.41,18.90,12.40,8.05,5.66,3.91,2.59,1.57,1.04,0.68,0.45,0.28);
      
      Print(P_excreted,0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4);
      Data(P_excreted,1.2,2.5,4.7,6.8,10.7,10.8,11.72,11.72);
    }
    
    
  }}

End.    
Simulation {    # Dose 1 # dose 0.5
  
  expT = 0.001;
  oraldose = 1.25;
  oraldose1 = PerDose (oraldose,24,0, expT);
  
  Print (cplasma,0.165,0.33,0.5,0.66,0.83,1,1.5,2,2.5, 3);
  Data  (cplasma,0.123,0.281,0.314,0.327,0.306,0.289,0.2229,0.193,0.134,0.062);
  Print(P_excreted, 0.5);
  Data  (P_excreted,-1);
  
}

Simulation {    # Dose 1 # dose 0.5
  
  expT = 0.001;
  oraldose = 1.25;
  oraldose1 = PerDose (oraldose,24,0, expT);
  
  Print (cplasma,0.083,0.165,0.33,0.5,0.66,0.83,1,1.5,2,2.5, 3);
  Data  (cplasma,0.2,0.865,1.86,2.24,1.816,1.65,1.32,0.9,0.754,0.528,0.32);
  Print(P_excreted, 0.5);
  Data  (P_excreted,-1);
  
  
}

}}

End. 


