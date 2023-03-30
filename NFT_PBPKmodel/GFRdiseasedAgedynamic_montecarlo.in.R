MonteCarlo("GFRMonteHumanPBPK.out",2000,10101010);

Distrib(age,Uniform,25,40);
Distrib(Kgut_plasma,LogNormal,0.622,1.17);
Distrib(Kliver_plasma,LogNormal,0.651,1.17);
Distrib(Kkidney_plasma,LogNormal,0.671,1.17);
Distrib(Kfat_plasma,LogNormal,0.159,1.17);
Distrib(Krestbody_plasma,LogNormal,0.42,1.17);
Distrib(fu,LogNormal,0.41,1.17);
Distrib(QurineC,LogNormal,13.45,1.17);
Distrib(Trc,LogNormal,1.334,1.17);
Distrib(Tmc,LogNormal,8.024,1.17);
Distrib(Kt,LogNormal,0.059,1.17);
Distrib(kfeces,LogNormal,0.0187,1.17);
Distrib(VmaxC,LogNormal,0.472,1.17);
Distrib(Km,LogNormal,5.83,1.17);
Distrib(Vehrc,LogNormal,0.527,1.17);
Distrib(Kehr,LogNormal,0.017,1.17);
Distrib(kbile,LogNormal,0.25,1.17);
Distrib(kgutabs,LogNormal,2.11,1.17);
Distrib(Kbp,LogNormal,0.76,1.17);


Simulation{   #normal
  sex=2;
  expT = 0.001;
  CKD = 0;   # There is a kidney disease (if zero then no kidney disease)
  GFR = 0;   # 60 ml/min for moderate kideny disease
  oraldose = 50;
  oraldose1 = PerDose (oraldose,6,0, expT);
  PrintStep(cplasma,0,120,0.1);
  PrintStep(P_excreted,0,120,0.1);
  PrintStep(cfat,0,120,0.1);
  PrintStep(cliver,0,120,0.1);
  PrintStep(cfilterate,0,120,0.1);
  PrintStep(ckidney,0,120,0.1);
  PrintStep(Aurine,0,120,0.1);
  PrintStep(crestbody,0,120,0.1);
  PrintStep(Filteration,0,120,0.1);
  PrintStep(Secretion,0,120,0.1);
  PrintStep(Reabsorption,0,120,0.1);
}

Simulation{   #mild
  sex=2;
  expT = 0.001;
  CKD = 1;   # There is a kidney disease (if zero then no kidney disease)
  GFR = 70;   # 60 ml/min for moderate kideny disease
  oraldose = 50;
  oraldose1 = PerDose (oraldose,6,0, expT);
  PrintStep(cplasma,0,120,0.1);
  PrintStep(P_excreted,0,120,0.1);
  PrintStep(cfat,0,120,0.1);
  PrintStep(cliver,0,120,0.1);
  PrintStep(cfilterate,0,120,0.1);
  PrintStep(ckidney,0,120,0.1);
  PrintStep(Aurine,0,120,0.1);
  PrintStep(crestbody,0,120,0.1);
  PrintStep(Filteration,0,120,0.1);
  PrintStep(Secretion,0,120,0.1);
  PrintStep(Reabsorption,0,120,0.1);
}

Simulation{   #moderate
  sex=2;
  expT = 0.001;
  CKD = 1;   # There is a kidney disease (if zero then no kidney disease)
  GFR = 45;   # 60 ml/min for moderate kideny disease
  oraldose = 50;
  oraldose1 = PerDose (oraldose,6,0, expT);
  PrintStep(cplasma,0,120,0.1);
  PrintStep(P_excreted,0,120,0.1);
  PrintStep(cfat,0,120,0.1);
  PrintStep(cliver,0,120,0.1);
  PrintStep(cfilterate,0,120,0.1);
  PrintStep(ckidney,0,120,0.1);
  PrintStep(Aurine,0,120,0.1);
  PrintStep(crestbody,0,120,0.1);
  PrintStep(Filteration,0,120,0.1);
  PrintStep(Secretion,0,120,0.1);
  PrintStep(Reabsorption,0,120,0.1);
}

Simulation{   #severe
  sex=2;
  expT = 0.001;
  CKD = 1;   # There is a kidney disease (if zero then no kidney disease)
  GFR = 20;   # 60 ml/min for moderate kideny disease
  oraldose = 50;
  oraldose1 = PerDose (oraldose,6,0, expT);
  PrintStep(cplasma,0,120,0.1);
  PrintStep(P_excreted,0,120,0.1);
  PrintStep(cfat,0,120,0.1);
  PrintStep(cliver,0,120,0.1);
  PrintStep(cfilterate,0,120,0.1);
  PrintStep(ckidney,0,120,0.1);
  PrintStep(Aurine,0,120,0.1);
  PrintStep(crestbody,0,120,0.1);
  PrintStep(Filteration,0,120,0.1);
  PrintStep(Secretion,0,120,0.1);
  PrintStep(Reabsorption,0,120,0.1);
}


Simulation{   #normal
  sex=2;
  expT = 0.001;
  CKD = 0;   # There is a kidney disease (if zero then no kidney disease)
  GFR = 70;   # 60 ml/min for moderate kideny disease
  oraldose = 100;
  oraldose1 = PerDose (oraldose,12,0, expT);
  PrintStep(cplasma,0,120,0.1);
  PrintStep(P_excreted,0,120,0.1);
  PrintStep(cfat,0,120,0.1);
  PrintStep(cliver,0,120,0.1);
  PrintStep(cfilterate,0,120,0.1);
  PrintStep(ckidney,0,120,0.1);
  PrintStep(Aurine,0,120,0.1);
  PrintStep(crestbody,0,120,0.1);
  PrintStep(Filteration,0,120,0.1);
  PrintStep(Secretion,0,120,0.1);
  PrintStep(Reabsorption,0,120,0.1);
}

Simulation{   #mild
  sex=2;
  expT = 0.001;
  CKD = 1;   # There is a kidney disease (if zero then no kidney disease)
  GFR = 70;   # 60 ml/min for moderate kideny disease
  oraldose = 100;
  oraldose1 = PerDose (oraldose,12,0, expT);
  PrintStep(cplasma,0,120,0.1);
  PrintStep(P_excreted,0,120,0.1);
  PrintStep(cfat,0,120,0.1);
  PrintStep(cliver,0,120,0.1);
  PrintStep(cfilterate,0,120,0.1);
  PrintStep(ckidney,0,120,0.1);
  PrintStep(Aurine,0,120,0.1);
  PrintStep(crestbody,0,120,0.1);
  PrintStep(Filteration,0,120,0.1);
  PrintStep(Secretion,0,120,0.1);
  PrintStep(Reabsorption,0,120,0.1);
}

Simulation{   #moderate
  sex=2;
  expT = 0.001;
  CKD = 1;   # There is a kidney disease (if zero then no kidney disease)
  GFR = 45;   # 60 ml/min for moderate kideny disease
  oraldose = 100;
  oraldose1 = PerDose (oraldose,12,0, expT);
  PrintStep(cplasma,0,120,0.1);
  PrintStep(P_excreted,0,120,0.1);
  PrintStep(cfat,0,120,0.1);
  PrintStep(cliver,0,120,0.1);
  PrintStep(cfilterate,0,120,0.1);
  PrintStep(ckidney,0,120,0.1);
  PrintStep(Aurine,0,120,0.1);
  PrintStep(crestbody,0,120,0.1);
  PrintStep(Filteration,0,120,0.1);
  PrintStep(Secretion,0,120,0.1);
  PrintStep(Reabsorption,0,120,0.1);
}

Simulation{   #severe
  sex=2;
  expT = 0.001;
  CKD = 1;   # There is a kidney disease (if zero then no kidney disease)
  GFR = 20;   # 60 ml/min for moderate kideny disease
  oraldose = 100;
  oraldose1 = PerDose (oraldose,12,0, expT);
  PrintStep(cplasma,0,120,0.1);
  PrintStep(P_excreted,0,120,0.1);
  PrintStep(cfat,0,120,0.1);
  PrintStep(cliver,0,120,0.1);
  PrintStep(cfilterate,0,120,0.1);
  PrintStep(ckidney,0,120,0.1);
  PrintStep(Aurine,0,120,0.1);
  PrintStep(crestbody,0,120,0.1);
  PrintStep(Filteration,0,120,0.1);
  PrintStep(Secretion,0,120,0.1);
  PrintStep(Reabsorption,0,120,0.1);
}
