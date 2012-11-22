#include "TMath.h"
#include "TLorentzVector.h"
#include "ZZMatrixElement/MELA/interface/TUtil.hh"
#include <stdio.h>
#include <math.h>
#include <iostream>

using namespace std;

void SetEwkCoupligParameters(){
  
  ewinput_.Gf_inp=1.16639E-05;
  ewinput_.aemmz_inp=7.81751E-03;
  ewinput_.wmass_inp=79.956049884402844;
  ewinput_.zmass_inp=91.1876;
  ewinput_.xw_inp=0.23116864;

}


void My_choose(TVar::Process process){
 
//ZZ_4l
if(process==TVar::ZZ_2e2m ){ 
 
    //81 '  f(p1)+f(p2) --> Z^0(-->mu^-(p3)+mu^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))'
    //86 '  f(p1)+f(p2) --> Z^0(-->e^-(p5)+e^+(p6))+Z^0(-->mu^-(p3)+mu^+(p4)) (NO GAMMA*)'
    //    nproc_.nproc=81;  
    //    chooser_();
  
  // these settings are identical to use the chooser_() function
    npart_.npart=4;
    nqcdjets_.nqcdjets=0;

    vsymfact_.vsymfact=1.0;                                                                                                               
    interference_.interference=false;

    nwz_.nwz=0;
    bveg1_mcfm_.ndim=10;
    masses_mcfm_.mb=0;
    breit_.n2=1;
    breit_.n3=1;

    breit_.mass2=masses_mcfm_.zmass;
    breit_.width2=masses_mcfm_.zwidth;
    breit_.mass3=masses_mcfm_.zmass;
    breit_.width3=masses_mcfm_.zwidth;

    zcouple_.q1=-1.0;
    zcouple_.l1=zcouple_.le;
    zcouple_.r1=zcouple_.re;

    zcouple_.q2=-1.0;
    zcouple_.l2=zcouple_.le;
    zcouple_.r2=zcouple_.re;


 }  else if ( process == TVar::ZZ_4e) {

    // 90 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))' 'L'
    // nproc_.nproc=90;
    //  chooser_();

    // these settings are  from 
    // ProdHep/chooser.f
    npart_.npart=4;
    nqcdjets_.nqcdjets=0;

    vsymfact_.vsymfact=0.25;                                                                                                               
    interference_.interference=true;

    nwz_.nwz=0;
    bveg1_mcfm_.ndim=10;
    masses_mcfm_.mb=0;
    breit_.n2=1;
    breit_.n3=1;

    breit_.mass2=masses_mcfm_.zmass;
    breit_.width2=masses_mcfm_.zwidth;
    breit_.mass3=masses_mcfm_.zmass;
    breit_.width3=masses_mcfm_.zwidth;

    zcouple_.q1=-1.0;
    zcouple_.l1=zcouple_.le;
    zcouple_.r1=zcouple_.re;

    zcouple_.q2=-1.0;
    zcouple_.l2=zcouple_.le;
    zcouple_.r2=zcouple_.re;

    
 }  else if ( process == TVar::HZZ_4l) {

    // 114 '  f(p1)+f(p2) --> H(--> Z^0(mu^-(p3)+mu^+(p4)) + Z^0(e^-(p5)+e^+(p6))' 'N'
    // nproc_.nproc=114;
    // chooser_();
     
     npart_.npart=4;
     nqcdjets_.nqcdjets=0;

     bveg1_mcfm_.ndim=10;
     masses_mcfm_.mb=0;

     breit_.n2=1;
     breit_.n3=1;

     breit_.mass2 =masses_mcfm_.zmass;
     breit_.width2=masses_mcfm_.zwidth;
     breit_.mass3 =masses_mcfm_.zmass;
     breit_.width3=masses_mcfm_.zwidth;
     
     zcouple_.l1=zcouple_.le;
     zcouple_.r1=zcouple_.re;
     
     zcouple_.l2=zcouple_.le;
     zcouple_.r2=zcouple_.re;

 } 
 else{
     std::cerr <<"[My_choose]: Can't identify Process: " << process <<endl;
 } 
}

bool My_masscuts(double s[][12],TVar::Process process){

 double minZmassSqr=10*10;

 if(process==TVar::ZZ_2e2m){
   if(s[2][3]< minZmassSqr) return true;
   if(s[4][5]< minZmassSqr) return true;
 }
 return false;	 

}


bool My_smalls(double s[][12],int npart){

// Reject event if any s(i,j) is too small
// cutoff is defined in technical.Dat
	
      if ( 
       npart == 3 &&
       (
        (-s[5-1][1-1]< cutoff_.cutoff)  //gamma p1
     || (-s[5-1][2-1]< cutoff_.cutoff)  //gamma p2
     || (-s[4-1][1-1]< cutoff_.cutoff)  //e+    p1
     || (-s[4-1][2-1]< cutoff_.cutoff)  //e-    p2
     || (-s[3-1][1-1]< cutoff_.cutoff)  //nu    p1
     || (-s[3-1][2-1]< cutoff_.cutoff)  //nu    p2
     || (+s[5-1][4-1]< cutoff_.cutoff)  //gamma e+
     || (+s[5-1][3-1]< cutoff_.cutoff)  //gamma nu
     || (+s[4-1][3-1]< cutoff_.cutoff)  //e+    nu
	)	 
      ) 
        return true;
     
     else if (
       npart == 4 &&     
      (
        (-s[5-1][1-1]< cutoff_.cutoff)  //e-    p1
     || (-s[5-1][2-1]< cutoff_.cutoff)  //e-    p2
     || (-s[6-1][1-1]< cutoff_.cutoff)  //nb    p1
     || (-s[6-1][2-1]< cutoff_.cutoff)  //nb    p2
     || (+s[6-1][5-1]< cutoff_.cutoff)  //e-    nb
       )

     )
       
      return true;
     
     return false;
}




//Make sure
// 1. tot Energy Sum < 2EBEAM
// 2. PartonEnergy Fraction minimum<x0,x1<1
// 3. number of final state particle is defined
//
double SumMatrixElementPDF(TVar::Process process, mcfm_event_type* mcfm_event,double flavor_msq[nmsq][nmsq],double* flux, double EBEAM){

  int NPart=npart_.npart+2;
  double p4[4][12];
  double fx1[nmsq];
  double fx2[nmsq];
  double msq[nmsq][nmsq];
  
  
  //Parton Density Function is always evalualted at pT=0 frame
  //Make sure parton Level Energy fraction is [0,1]
  //phase space function already makes sure the parton energy fraction between [min,1]
  //  x0 EBeam =>   <= -x1 EBeam
  
  double sysPz=mcfm_event->p[0].Pz()    +mcfm_event->p[1].Pz();
  double sysE =mcfm_event->p[0].Energy()+mcfm_event->p[1].Energy();
  
  //Ignore the Pt doesn't make significant effect
  //double sysPt_sqr=sysPx*sysPx+sysPy*sysPy;
  //if(sysPt_sqr>=1.0E-10)  sysE=TMath::Sqrt(sysE*sysE-sysPt_sqr);
  
  double xx[2]={(sysE+sysPz)/EBEAM/2,(sysE-sysPz)/EBEAM/2};
  if(xx[0] > 1.0 || xx[0]<=xmin_.xmin) return 0.0;
  if(xx[1] > 1.0 || xx[1]<=xmin_.xmin) return 0.0;
  
  //Convert TLorentzVector into 4x12 Matrix
  //reverse sign of incident partons back
  for(int ipar=0;ipar<2;ipar++){    
    if(mcfm_event->p[ipar].Energy()>0){
      p4[0][ipar] = -mcfm_event->p[ipar].Px();
      p4[1][ipar] = -mcfm_event->p[ipar].Py();
      p4[2][ipar] = -mcfm_event->p[ipar].Pz();
      p4[3][ipar] = -mcfm_event->p[ipar].Energy();
    }
  }
  //initialize decayed particles
  for(int ipar=2;ipar<NPart;ipar++){
    
    p4[0][ipar] = mcfm_event->p[ipar].Px();
    p4[1][ipar] = mcfm_event->p[ipar].Py();
    p4[2][ipar] = mcfm_event->p[ipar].Pz();
    p4[3][ipar] = mcfm_event->p[ipar].Energy();
    
  }
  
  //calculate invariant masses between partons/final state particles
  double s[12][12];
  for(int jdx=0;jdx< NPart ;jdx++){
    s[jdx][jdx]=0;
    for(int kdx=jdx+1;kdx<NPart;kdx++){
      s[jdx][kdx]=2*(p4[3][jdx]*p4[3][kdx]-p4[2][jdx]*p4[2][kdx]-p4[1][jdx]*p4[1][kdx]-p4[0][jdx]*p4[0][kdx]);
      s[kdx][jdx]=s[jdx][kdx];
    }
  }
  
  
  //remove events has small invariant mass
  if(My_masscuts(s,process)) return 0.0;
  if(My_smalls(s,npart_.npart)) return 0.0;
  
  
  //Calculate Pdf
  //Always pass address through fortran function
  fdist_ (&density_.ih1, &xx[0], &scale_.scale, fx1); 
  fdist_ (&density_.ih2, &xx[1], &scale_.scale, fx2); 
  
  if( process==TVar::ZZ_2e2m || process==TVar::ZZ_4e )      qqb_zz_(p4[0],msq[0]);
  if( process==TVar::HZZ_4l)     qqb_hzz_(p4[0],msq[0]);
  
  double msqjk=0;
  for(int ii=0;ii<nmsq;ii++){
    for(int jj=0;jj<nmsq;jj++){
      
      //2-D matrix is reversed in fortran
      // msq[ parton2 ] [ parton1 ]
      //      flavor_msq[jj][ii] = fx1[ii]*fx2[jj]*msq[jj][ii];

      flavor_msq[jj][ii] = msq[jj][ii];
      //cout<<jj<<ii<<"="<<msq[jj][ii]<<"  ";
      msqjk+=flavor_msq[jj][ii];
    }//ii
    //    cout<<"\n";
  }//jj

  if( process==TVar::ZZ_2e2m) msqjk=msq[3][7]+msq[7][3];
  
  (*flux)=fbGeV2/(8*xx[0]*xx[1]*EBEAM*EBEAM);

  if(msqjk != msqjk || flux!=flux ){
    cout << "SumMatrixPDF: "<< TVar::ProcessName(process) << " msqjk="  << msqjk << " flux="<< flux <<endl;
    msqjk=0;
    flux=0;
  }
  return msqjk;
  
}

//
// Test code from Markus to calculate the HZZ cross-section
// 
double JHUGenMatEl(TVar::Process process, mcfm_event_type* mcfm_event, double MReso, double GaReso, double *xggcoupl, double *xvvcoupl) 
{
  // input unit = GeV/100 such that 125GeV is 1.25 in the code
  // this needs to be applied for all the p4
  MReso = MReso / 100.0;
  GaReso = GaReso /100.0;
  double p4[6][4];
  double MatElSq=0;
  int MYIDUP[4];

  int NPart = 6; 
  // p(i,0:3) = (E(i),px(i),py(i),pz(i))
  // i=0,1: glu1,glu2 (outgoing convention)
  // i=2,3: correspond to MY_IDUP(1),MY_IDUP(0)
  // i=4,5: correspond to MY_IDUP(3),MY_IDUP(2)
  for(int ipar=0;ipar<2;ipar++){   
    if(mcfm_event->p[ipar].Energy()>0){
      p4[ipar][0] = -mcfm_event->p[ipar].Energy()/100.;
      p4[ipar][1] = -mcfm_event->p[ipar].Px()/100.;
      p4[ipar][2] = -mcfm_event->p[ipar].Py()/100.;
      p4[ipar][3] = -mcfm_event->p[ipar].Pz()/100.;
    }
  }
  //initialize decayed particles
  for(int ipar=2;ipar<NPart;ipar++){
    p4[ipar][0] = mcfm_event->p[ipar].Energy()/100.;
    p4[ipar][1] = mcfm_event->p[ipar].Px()/100.;
    p4[ipar][2] = mcfm_event->p[ipar].Py()/100.;
    p4[ipar][3] = mcfm_event->p[ipar].Pz()/100.;
  }
  
  // particle ID: +7=e+,  -7=e-,  +8=mu+,  -8=mu-

  if ( TMath::Abs(mcfm_event->PdgCode[2]) == TMath::Abs(mcfm_event->PdgCode[3]) && 
       TMath::Abs(mcfm_event->PdgCode[3]) == TMath::Abs(mcfm_event->PdgCode[4]) && 
       TMath::Abs(mcfm_event->PdgCode[4]) == TMath::Abs(mcfm_event->PdgCode[5]) ) {
    if ( TMath::Abs(mcfm_event->PdgCode[2]) == 11  ) {
      MYIDUP[0]=+7;
      MYIDUP[1]=-7;
      MYIDUP[2]=+7;
      MYIDUP[3]=-7;
    } 
    if ( TMath::Abs(mcfm_event->PdgCode[2]) == 13  ) {
      MYIDUP[0]=+8;
      MYIDUP[1]=-8;
      MYIDUP[2]=+8;
      MYIDUP[3]=-8;
    } 
  } else {
      MYIDUP[0]=+7;
      MYIDUP[1]=-7;
      MYIDUP[2]=+8;
      MYIDUP[3]=-8;
  }
  if ( process == TVar::HZZ_4l || process == TVar::PSHZZ_4l ) {
    __modhiggs_MOD_evalamp_gg_h_vv(p4, &MReso,  &GaReso, xggcoupl, xvvcoupl, MYIDUP, &MatElSq);
  }
  if ( process == TVar::TZZ_4l ) {
    __modgraviton_MOD_evalamp_gg_g_vv(p4, &MReso,  &GaReso, xggcoupl, xvvcoupl, MYIDUP, &MatElSq);
  }
  if ( process == TVar::VZZ_4l ) {
    // -- YY: note that even if it is called xggcouplings, we are only testing xqq!
    __modzprime_MOD_evalamp_qqb_zprime_vv(p4, &MReso,  &GaReso, xggcoupl, xvvcoupl, MYIDUP, &MatElSq);
  }

  /*
  printf("\n ");
  std::cout << "resoance = " << MReso *100. << ", width = " << GaReso*100. << "\n";
  for ( int i=0; i<NPart;i++) {
    std::cout << "p["<<i<<"] (E, Px, Py, Pz) = (" << p4[i][0] << ", " << p4[i][1] << ", " << p4[i][2] << ", " << p4[i][3] << ")\n";
  }
  printf("Matr.el. squared: %20.17e \n ",MatElSq);
  */
  // 
  // This constant is needed to account for the different units used in 
  // JHUGen compared to the MCFM
  // 
  double constant = 1.45/pow(10, 8);
  return MatElSq*constant;

}

