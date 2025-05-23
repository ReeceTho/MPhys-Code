/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs  
================================*/
      
#define MASSES_INFO      
  /* Display information about mass spectrum  */
  
#define CONSTRAINTS 
//#define SMODELS
//#define MONOJET
//#define HIGGSBOUNDS 
//#define HIGGSSIGNALS
//#define LILITH
//#define SMODELS
  
#define OMEGA       /*  Calculate Freeze out relic density and display contribution of  individual channels */
//#define FREEZEIN  /*  Calculate relic density in Freeze-in scenario  */
   
#define INDIRECT_DETECTION  
  /* Compute spectra of gamma/positron/antiprotons/neutrinos for DM annihilation; 
     Calculate <sigma*v>;
     Integrate gamma signal over DM galactic squared density for given line 
     of sight; 
     Calculate galactic propagation of positrons and antiprotons.      
  */
      
//#define RESET_FORMFACTORS
  /* Modify default nucleus form factors, 
    DM velocity distribution,
    A-dependence of Fermi-dencity
  */     
#define CDM_NUCLEON     
  /* Calculate amplitudes and cross-sections for  CDM-mucleon collisions */  

#define CDM_NUCLEUS     
     // Calculate  exclusion rate for direct detection experiments Xenon1T, DarkSide50, CRESST, and PICO 
         
//#define NEUTRINO    
 /*  Neutrino signal of DM annihilation in Sun and Earth */

#define DECAYS

//#define CROSS_SECTIONS 
  
/*===== end of Modules  ======*/

/*===== Options ========*/
/*#define SHOWPLOTS*/
     /* Display  graphical plots on the screen */ 
#define CLEAN
/*===== End of DEFINE  settings ===== */


#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"

double rnd_slg(double xmin, double xmax) {
    
    double sgn=1.;
    if(xmin<0) { 
    xmin=1E-10;
    sgn= (rand() % 2 == 0) ? 1 : -1; // Randomly return +1 or -1 
    }
    
    double log_min = log(xmin);
    double log_max = log(xmax);
    
    // Generate a random value on the log scale
    double random_log = log_min + ((double)rand() / RAND_MAX) * (log_max - log_min);
   
    return sgn*exp(random_log);
      }

double rnd_lin(double xmin, double xmax) {
    return xmin + (xmax - xmin) * ((double)rand() / RAND_MAX);
}


int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   int spin2, charge3,cdim;

  ForceUG=0;  /* to Force Unitary Gauge assign 1 */
  //useSLHAwidth=0;
  VZdecay=0; VWdecay=0;  
 
/* ================the new ppiece of code starts here =================*/
   int a=0,n_points=10000;

   double MD1,ld345,DMP,DM3,DM2,ld,protonSI,Omh2,sigV,PvalDD,CMB_ID,MW,MZ,MH;
   char* expDD; 
   MW = 80.3692;
   MZ = 91.188;
   MH = 125.2;
                                
  FILE *fscan = fopen("scan.dat","w"); 
  fprintf(fscan,"%s %s %s %s %s %s %s %s %s %s %s %s \n", 
 "MD1          ","DMP          ","DM3          ","l345         ","ld          ",
 "Omegah2      ","sigV         ","protonSI     ","PvalDD       ","expName     ",
 "CMB_ID       ","brH_DMDM     ");

//	Set limit of scan here

double MD1_min  =1,  MD1_max  =1e+4;
double DMP_min  =1e-4, DMP_max  =1e+3; //DMP=Mh+-Mh1
double DM2_min  =1e-4,  DM2_max  =1e+3; //introducing DM2 which is Mh2-Mh1
double ld345_min=-10 ,  ld345_max=10;
double ld_min = 1e-4, ld_max = 10;

// Initialize the random seed using the current time
srand(time(NULL));

for (;;) { /* remember to close the loop in the end! */
 	MD1 = rnd_slg(MD1_min,MD1_max); DMP = rnd_slg(DMP_min,DMP_max);
	DM2 = rnd_slg(DM2_min,DM2_max); /* DM2 is just a MD2-MD1 mass difference, which I conver to DM3 value below */
 	ld345 = rnd_slg(ld345_min,ld345_max);
    ld = rnd_slg(ld_min, ld_max);
	DM3=DM2-DMP;  /* Note that Mh2-Mh+ =  DM2-DMP */

      printf("\n MD1=%.20E    \n",MD1);
      printf("DMP=%.20E    \n",DMP);
      printf("DM3=%.20E    \n",DM3);
      printf("ld345=%.20E    \n",ld345);
      printf("ld=%.20E    \n", ld);

	
	
	assignValW("MD1",MD1);
	assignValW("DMP",DMP);
	assignValW("DM3",DM3);
	assignValW("ld345",ld345);
	assignValW("ld",ld);
	
  	err=sortOddParticles(cdmName);
  	if(err) { printf("Can't calculate %s\n",cdmName); return 1;}
  
  	txtList HDecay;
	double HWidth=pWidth("H",&HDecay);
	double brH_DMDM = findBr(HDecay,"~h1,~h1");
	
/* ========= the new ppiece of code end  here =============*/



//  err=sortOddParticles(cdmName);
//  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}
  
  for(int k=1;k<=Ncdm;k++) 
  { 
     qNumbers(CDM[k], &spin2, &charge3, &cdim);
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM[k],  spin2,McdmN[k]); 
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }
  
#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF HIGGS AND ODD PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
}
#endif

#ifdef CONSTRAINTS
{ double csLim;
  if(Zinvisible()) printf("Excluded by Z->invizible\n");
  if(LspNlsp_LEP(&csLim)) printf("LEP excluded by e+,e- -> DM q q-\\bar  Cross Section= %.2E pb\n",csLim);
}
#endif

#ifdef SMODELS
{
  int combineSRs=0;
  char* combineAnas=NULL;//"ATLAS-SUSY-2018-05-ewk,ATLAS-SUSY-2019-08,ATLAS-SUSY-2019-09";
  
  int status=0, smodelsOK=0; 
  double Rvalue, Rexpected, SmoLsig, SmoLmax, SmoLSM;
  double CombRvalue, CombRexpected, CombSmoLsig, CombSmoLmax, CombSmoLSM;

  char analysis[50]={},topology[100]={},smodelsInfo[100];
  char CombAnalyses[200]={};
  int LHCrun=LHC8|LHC13;  //  LHC8 - 8TeV; LHC13 - 13TeV;   
//  int LHCrun=LHC13;  //  LHC13 - 13TeV only;   

  printf("\n\n=====  LHC constraints with SModelS  =====\n\n");

#include "../include/SMODELS.inc" // SLHA interface with SModelS

  printf("SModelS %s \n",smodelsInfo);
  if(smodelsOK) 
  { printf("\n highest r-value = %.2E",Rvalue); 

    if(Rvalue>0) 
    { printf(" from %s, topology: %s ",analysis,topology);
      if(Rexpected>0) 
      { printf("\n expected r = %.2E ",Rexpected);
        if(SmoLsig!=INFINITY)
        { printf("\n -2log (L_signal/L_max, L_SM/L_max) = %.2E %.2E",
                  2*(SmoLsig-SmoLmax), 2*(SmoLSM-SmoLmax) ); 
        }                  
      }
    }  
    if(status==1) { printf("\n excluded by SMS results"); }
    else if(status==0) printf("\n not excluded"); 
    else if(status==-1) printf("\n not not tested by results in SModelS database"); 
    printf("\n");

    // r-value and likelihoods from analysis cvombination
    if(CombRvalue>0) 
    { printf("\n Combination of %s",CombAnalyses);
      printf("\n r-value = %.2E (expected r = %.2E)",CombRvalue, CombRexpected); 
      if(CombRvalue>=1) printf("  --> excluded"); 
      else printf("  --> not excluded"); 
      printf("\n -2log (L_signal/L_max, L_SM/L_max) = %.2E %.2E \n\n", 
                    2*(CombSmoLsig-CombSmoLmax),2*(CombSmoLSM-CombSmoLmax));                     
    }

  } else system("cat smodels.err"); // problem: see smodels.err
}   

#endif 


#ifdef MONOJET
{ double CL=monoJet();
  printf(" Monojet signal exclusion CL is %.3e\n", CL);   
}  
#endif

#if defined(HIGGSBOUNDS) || defined(HIGGSSIGNALS)
{  int NH0=3, NHch=1; // number of neutral and charged Higgs particles.
   int HB_id[3]={0,0,0},HB_result[3];
   double  HB_obsratio[3],HS_observ=-1,HS_chi2, HS_pval;
   char HB_chan[3][100]={""}, HB_version[50], HS_version[50]; 
   NH0=hbBlocksMO("HB.in",&NHch); 
//    NH0= hbBlocksMDL("HB.in",&NHch); 
   system("echo 'BLOCK DMASS\n 25  2  '>> HB.in");
#include "../include/hBandS.inc"
#ifdef HIGGSBOUNDS
   printf("HiggsBounds(%s)\n", HB_version);
   for(int i=0;i<3;i++) if(HB_id[i]) printf("  id= %d  result = %d  obsratio=%.2E  channel= %s \n", HB_id[i],HB_result[i],HB_obsratio[i],HB_chan[i]);
#endif 
#ifdef HIGGSSIGNALS
   if(HS_observ>=0)
   {
     printf("HiggsSignals(%s)\n",HS_version); 
     printf("  Nobservables=%.0f chi^2 = %.2E pval= %.2E\n",HS_observ,HS_chi2, HS_pval);
   }
#endif   
}
#endif

#ifdef LILITH
{  double m2logL, m2logL_reference=0,pvalue;
   int exp_ndf,n_par=0,ndf;
   char Lilith_version[50];
   if(LilithMO("Lilith_in.xml"))
   {        
#include "../include/Lilith.inc"
      if(ndf)
      {
        printf("LILITH(DB%s):  -2*log(L): %.2f; -2*log(L_reference): %.2f; ndf: %d; p-value: %.2E \n",
        Lilith_version,m2logL,m2logL_reference,ndf,pvalue);
      }  
   } else printf("LILITH: there is no Higgs candidate\n");
}     
#endif


#ifdef SMODELS
{ int combineSRs=0;
  char* combineAnas=NULL;
  
  int status=0, smodelsOK=0; 
  double Rvalue, Rexpected, SmoLsig, SmoLmax, SmoLSM;
  double CombRvalue, CombRexpected, CombSmoLsig, CombSmoLmax, CombSmoLSM;

  char analysis[50]={},topology[100]={},smodelsInfo[100];
  char CombAnalyses[200]={};
  int LHCrun=LHC8|LHC13;  //  LHC8 - 8TeV; LHC13 - 13TeV;   
//  int LHCrun=LHC13;  //  LHC13 - 13TeV only;   

  printf("\n\n=====  LHC constraints with SModelS  =====\n\n");

#include "../include/SMODELS.inc" // SLHA interface with SModelS

  printf("SModelS %s \n",smodelsInfo);
  if(smodelsOK) 
  { printf("\n highest r-value = %.2E",Rvalue); 

    if(Rvalue>0) 
    { printf(" from %s, topology: %s ",analysis,topology);
      if(Rexpected>0) 
      { printf("\n expected r = %.2E ",Rexpected);
        if(SmoLsig!=INFINITY)
        { printf("\n -2log (L_signal/L_max, L_SM/L_max) = %.2E %.2E",
                  2*(SmoLsig-SmoLmax), 2*(SmoLSM-SmoLmax)); }
      }
    }  
    if(status==1) { printf("\n excluded by SMS results"); }
    else if(status==0) printf("\n not excluded"); 
    else if(status==-1) printf("\n not not tested by results in SModelS database"); 
    printf("\n");

    // r-value and likelihoods from analysis cvombination
    if(CombRvalue>0) 
    { printf("\n Combination of %s",CombAnalyses);
      printf("\n r-value = %.2E (expected r = %.2E)",CombRvalue, CombRexpected); 
      if(CombRvalue>=1) printf("  --> excluded"); 
      else printf("  --> not excluded"); 
      printf("\n -2log (L_signal/L_max, L_SM/L_max) = %.2E %.2E \n\n", 
                    -2*(CombSmoLsig-CombSmoLmax),-2*(CombSmoLSM-CombSmoLmax)); 
    }

  } else system("cat smodels.err"); // problem: see smodels.err
}   

#endif 


#ifdef OMEGA
{ int fast=1;
  double Beps=1.E-4, cut=0.01;
  double Omega;  
  int i,err; 
  printf("\n==== Calculation of relic density =====\n");   

  if(Ncdm==1) 
  {  double Xf;
     Omega=darkOmega(&Xf,fast,Beps,&err);
     printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);
     if(Omega>0)printChannels(Xf,cut,Beps,1,stdout);
  } else
  if(Ncdm==2)
  {
    Omega= darkOmega2(fast,Beps,&err);
    printf("Omega_1h^2=%.2E Omega_2h^2=%.2E err=%d \n", Omega*fracCDM[1], Omega*fracCDM[2],err);
  }else 
  {  
     Omega=darkOmegaN(fast,Beps,&err);
     printf("Omega=%.2E\n",Omega);
     for(int k=1;k<=Ncdm;k++) printf("   Omega_%d=%.2E\n",k,Omega*fracCDM[k]); 
  }
  
   if (isnan(Omega)) {
            // Skip the current iteration if the variable is NaN
            printf("Skipping NaN for Omega \n");
	    continue;
        }
   Omh2=Omega;//my addition   
}

#endif

#ifdef FREEZEIN
{
  double TR=1E6;
  double omegaFi;  
  toFeebleList(CDM[1]);
  VWdecay=0; VZdecay=0;
  
  omegaFi=darkOmegaFi(TR,CDM[1],&err);
  printf("omega freeze-in=%.3E\n", omegaFi);
  printChannelsFi(0,0,stdout);
}
#endif



#ifdef INDIRECT_DETECTION
{ 
  int err,i;
  double Emin=1,/* Energy cut  in GeV   */  sigmaV;
  double vcs_gz,vcs_gg;
  char txt[100];
  double SpA[NZ],SpE[NZ],SpP[NZ];
  double FluxA[NZ],FluxE[NZ],FluxP[NZ];
  double * SpNe=NULL,*SpNm=NULL,*SpNl=NULL;
  double Etest=Mcdm/2;
  
printf("\n==== Indirect detection =======\n");  

  sigmaV=calcSpectrum(1+2+4,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
  sigV=sigmaV ;//My addition 
  
    /* Returns sigma*v in cm^3/sec.     SpX - calculated spectra of annihilation.
       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.
       
       First parameter 1-includes W/Z polarization
                       2-includes gammas for 2->2+gamma
                       4-print cross sections             
    */

   CMB_ID=PlanckCMB(sigmaV , SpA, SpE); //My addition 
   printf("CMB ID/3.2x10^-28 (cm^3/s/GeV) = %.3e\n",CMB_ID); //My addition 


//  if(SpA)
  { 
     double fi=0.1,dfi=0.05; /* angle of sight and 1/2 of cone angle in [rad] */ 

     gammaFluxTab(fi,dfi, sigmaV, SpA,  FluxA);     
     printf("Photon flux  for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.2f[rad]\n",fi,2*dfi);
#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux for angle of sight %.2f[rad] and cone angle %.2f[rad]",fi,2*dfi);
     displayPlot(txt,"E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxA);
#endif
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, FluxA), Etest);       
  }

//  if(SpE)
  { 
    posiFluxTab(Emin, sigmaV, SpE,  FluxE);
#ifdef SHOWPLOTS     
    displayPlot("positron flux [cm^2 s sr GeV]^{-1}","E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxE);
#endif
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest);           
  }
  
//  if(SpP)
  { 
    pbarFluxTab(Emin, sigmaV, SpP,  FluxP  ); 
#ifdef SHOWPLOTS    
     displayPlot("antiproton flux [cm^2 s sr GeV]^{-1}","E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxP);
#endif
    printf("Antiproton flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxP),  Etest);             
  }
}  
#endif

#ifdef RESET_FORMFACTORS
{
/* 
   The user has approach to form factors  which specifies quark contents 
   of  proton and nucleon via global parametes like
      <Type>FF<Nucleon><q>
   where <Type> can be "Scalar", "pVector", and "Sigma"; 
         <Nucleon>     "P" or "N" for proton and neutron
         <q>            "d", "u","s"

   calcScalarQuarkFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigmaS[MeV])  
   calculates and rewrites Scalar form factors
*/
  printf("\n======== RESET_FORMFACTORS ======\n");
 
  printf("protonFF (default) d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(default) d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);
//                    To restore default form factors of  version 2  call 
     calcScalarQuarkFF(0.553,18.9,55.,243.5);


  printf("protonFF (new)     d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(new)     d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

//                    To restore default form factors  current version  call 
//  calcScalarQuarkFF(0.56,20.2,34,42);


}
#endif

#ifdef CDM_NUCLEON
{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;        
  double csSIp1,csSIn1,csSDp1,csSDn1, csSIp1_,csSIn1_,csSDp1_,csSDn1_;
  double csSIp2,csSIn2,csSDp2,csSDn2, csSIp2_,csSIn2_,csSDp2_,csSDn2_;
printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");   

  for(int k=1;k<=Ncdm;k++) 
  {  
    nucleonAmplitudes(CDM[k], pA0,pA5,nA0,nA5);
    printf("%s[%s]-nucleon micrOMEGAs amplitudes\n",CDM[k],antiParticle(CDM[k]));
    printf("proton:  SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",pA0[0], pA0[1],  pA5[0], pA5[1] );
    printf("neutron: SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",nA0[0], nA0[1],  nA5[0], nA5[1] ); 

    SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*McdmN[k]/(Nmass+ McdmN[k]),2.);
    csSIp1=  SCcoeff*pA0[0]*pA0[0];  csSIp1_=  SCcoeff*pA0[1]*pA0[1];
    csSDp1=3*SCcoeff*pA5[0]*pA5[0];  csSDp1_=3*SCcoeff*pA5[1]*pA5[1];
    csSIn1=  SCcoeff*nA0[0]*nA0[0];  csSIn1_=  SCcoeff*nA0[1]*nA0[1];
    csSDn1=3*SCcoeff*nA5[0]*nA5[0];  csSDn1_=3*SCcoeff*nA5[1]*nA5[1];
                    
    printf("%s[%s]-nucleon cross sections[pb]:\n",CDM[k],antiParticle(CDM[k]));
    printf(" proton  SI %.3E [%.3E] SD %.3E [%.3E]\n", csSIp1,csSIp1_,csSDp1,csSDp1_);
    printf(" neutron SI %.3E [%.3E] SD %.3E [%.3E]\n", csSIn1,csSIn1_,csSDn1,csSDn1_);     
  }
  
  protonSI=csSIp1; //My addition
  
}
#endif
  
#ifdef CDM_NUCLEUS
{ char* expName; 
  
  printf("\n===== Direct detection exclusion:======\n");
  double pval=DD_pval(AllDDexp, Maxwell, &expName);
       if(pval<0.1 ) printf("Excluded by %s  %.1f%%\n", expName, 100*(1-pval)); 
  else {printf("Not excluded by DD experiments  at 90%% level \n"); 
     expName="NotAppl";
     } 
  
  PvalDD=pval;//my addition
  expDD=expName;//my addition
  
}
#endif 

#ifdef NEUTRINO
if(!CDM[1] || !CDM[2])
{ double nu[NZ], nu_bar[NZ],mu[NZ];
  double Ntot;
  int forSun=1;
  double Emin=1;
  
 printf("\n===============Neutrino Telescope=======  for  "); 
 if(forSun) printf("Sun\n"); else printf("Earth\n");  

  err=neutrinoFlux(Maxwell,forSun, nu,nu_bar);
#ifdef SHOWPLOTS
  displayPlot("neutrino fluxes [1/Year/km^2/GeV]","E[GeV]",Emin,Mcdm,0, 2,"dnu/dE",0,SpectdNdE,nu,"dnu_bar/dE",0,SpectdNdE,nu_bar);
#endif
{ 
    printf(" E>%.1E GeV neutrino flux       %.2E [1/Year/km^2] \n",Emin,spectrInfo(Emin,nu,NULL));
    printf(" E>%.1E GeV anti-neutrino flux  %.2E [1/Year/km^2]\n",Emin,spectrInfo(Emin,nu_bar,NULL));  
} 
  
/* Upward events */
  
  muonUpward(nu,nu_bar, mu);
#ifdef SHOWPLOTS
  displayPlot("Upward muons[1/Year/km^2/GeV]","E",Emin,Mcdm/2, 0,1,"mu",0,SpectdNdE,mu);  
#endif
    printf(" E>%.1E GeV Upward muon flux    %.2E [1/Year/km^2]\n",Emin,spectrInfo(Emin,mu,NULL));
  
/* Contained events */
  muonContained(nu,nu_bar,1., mu);
#ifdef SHOWPLOTS 
  displayPlot("Contained  muons[1/Year/km^3/GeV]","E",Emin,Mcdm,0,1,"",0,SpectdNdE,mu); 
#endif
  printf(" E>%.1E GeV Contained muon flux %.2E [1/Year/km^3]\n",Emin,spectrInfo(Emin/Mcdm,mu,NULL));
}        
#endif 


#ifdef DECAYS
{ char*  pname = pdg2name(25);
  txtList L;
  double width; 
  if(pname)
  { 
    width=pWidth(pname,&L);  
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    printTxtList(L,stdout);
  } 
  
}            
#endif

#ifdef CROSS_SECTIONS
{
  char* next,next_;
  double nextM;
    
  next=nextOdd(1,&nextM); 
  if(next && nextM<1000)  
  { 
     double cs, Pcm=6500, Qren, Qfact, pTmin=0;
     int nf=3;
     char*next_=antiParticle(next);
     Qren=Qfact=nextM; 
 
     printf("\npp > nextOdd  at sqrt(s)=%.2E GeV\n",2*Pcm);  
  
     Qren=Qfact;
     cs=hCollider(Pcm,1,nf,Qren, Qfact, next,next_,pTmin,1);
     printf("Production of 'next' odd particle: cs(pp-> %s,%s)=%.2E[pb]\n",next,next_, cs);
  }  
}
 
#endif 

#ifdef CLEAN
  system("rm -f HB.* HB.* hb.* hs.*  debug_channels.txt debug_predratio.txt  Key.dat");
  system("rm -f Lilith_*   particles.py*");
//  system("rm -f   smodels.*");  
#endif 



  killPlots();
  
 /*  Adding last lines -- starts */ 

	fprintf(fscan,"%.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E   %.6E %s  %.6E  %.6E  \n",
	                 MD1, DMP, DM3, ld345, ld, Omh2,sigV, protonSI,  PvalDD,expDD,CMB_ID,brH_DMDM);
	
	protonSI,Omh2,sigV,PvalDD,CMB_ID;
	
	a++;
	printf("N_line= %d\n",a);
  }
  fclose(fscan);
 /*  Adding last lines -- ends */ 
  
  return 0;
}
