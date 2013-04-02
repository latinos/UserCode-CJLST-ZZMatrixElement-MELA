#include <ZZMatrixElement/MELA/interface/SpinOneOddMELA.h>

#include "computeAngles.h"
#include "AngularPdfFactory.h"
#include "VectorPdfFactory.h"
#include <RooRealVar.h>

SpinOneOddMELA::SpinOneOddMELA(){

  z1mass_rrv = new RooRealVar("z1mass","m_{Z1}",0,180);
  z2mass_rrv = new RooRealVar("z2mass","m_{Z2}",0,120);
  costheta1_rrv = new RooRealVar("costheta1","cos#theta_{1}",-1,1);
  costheta2_rrv = new RooRealVar("costheta2","cos#theta_{2}",-1,1);
  phi_rrv= new RooRealVar("phi","#Phi",-3.1415,3.1415);
  costhetastar_rrv = new RooRealVar("costhetastar","cos#theta^{*}",-1,1);
  phistar1_rrv= new RooRealVar("phistar1","#Phi^{*}_{1}",-3.1415,3.1415);
  mzz_rrv= new RooRealVar("mzz","mZZ",80,1000);

  SMHiggs = new AngularPdfFactory(z1mass_rrv,z2mass_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,mzz_rrv);
  SMHiggs->makeSMHiggs();
  SMHiggs->makeParamsConst(true);

  sigAlt = new VectorPdfFactory(z1mass_rrv,z2mass_rrv,costhetastar_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,phistar1_rrv,mzz_rrv);

  sigAlt->makeZprime();
  sigAlt->makeParamsConst(true);

}

SpinOneOddMELA::~SpinOneOddMELA(){

  delete z1mass_rrv;
  delete z2mass_rrv;
  delete costheta1_rrv;
  delete costheta2_rrv;
  delete phi_rrv;
  delete costhetastar_rrv;
  delete phistar1_rrv;
  delete mzz_rrv;

  delete SMHiggs;
  delete sigAlt;

}

void SpinOneOddMELA::checkZorder(float& z1mass, float& z2mass,
			    float& costhetastar, float& costheta1,
			    float& costheta2, float& phi,
			    float& phistar1){
  
  float tempZ1mass=z1mass;
  float tempZ2mass=z2mass;
  float tempH1=costheta1;
  float tempH2=costheta2;
  float tempHs=costhetastar;
  float tempPhi1=phistar1;
  float tempPhi=phi;

  if(z2mass>z1mass){

    z1mass=tempZ2mass;
    z2mass=tempZ1mass;
    costhetastar=-tempHs;
    costheta1=tempH2;
    costheta2=tempH1;
    phi=tempPhi;
    phistar1=-tempPhi1-tempPhi;
    if(phistar1>3.1415)
      phistar1=phistar1-2*3.1415;
    if(phistar1<-3.1415)
      phistar1=phistar1+2*3.1415;

  }else

    return;

}

void SpinOneOddMELA::computeKD(TLorentzVector Z1_lept1, int Z1_lept1Id,
			      TLorentzVector Z1_lept2, int Z1_lept2Id,
			      TLorentzVector Z2_lept1, int Z2_lept1Id,
			      TLorentzVector Z2_lept2, int Z2_lept2Id,
			      float& kd, 
			      float& psig,
			      float& psigALT){

  //compute angles
  float m1=(Z1_lept1 + Z1_lept2).M();
  float m2=(Z2_lept1 + Z2_lept2).M();

  TLorentzVector ZZ = (Z1_lept1 + Z1_lept2 + Z2_lept1 + Z2_lept2);
  float mzz = ZZ.M();

  if(mzz<100.0 || mzz>179.9999){
    psig=0.0;
    psigALT=0.0;
    kd=0.0;
    return;
  }

  float costheta1,costheta2,costhetastar,phi,phistar1;

  mela::computeAngles(Z1_lept1, Z1_lept1Id, Z1_lept2, Z1_lept2Id, 
		      Z2_lept1, Z2_lept1Id, Z2_lept2, Z2_lept2Id,
		      costhetastar,costheta1,costheta2,phi,phistar1);

  //compute kd
  checkZorder(m1,m2,costhetastar,costheta1,costheta2,phi,phistar1);

  z1mass_rrv->setVal(m1);
  z2mass_rrv->setVal(m2);
  costhetastar_rrv->setVal(costhetastar);
  costheta1_rrv->setVal(costheta1);
  costheta2_rrv->setVal(costheta2);
  phi_rrv->setVal(phi);
  phistar1_rrv->setVal(phistar1);
  
  mzz_rrv->setVal(mzz);
  
  psig = SMHiggs->PDF->getVal();
  psigALT = sigAlt->PDF->getVal();
  kd = 1/(1+psigALT/psig);
  
}


void SpinOneOddMELA::computeKD(float zzmass, float z1mass, 
		       float z2mass, float costhetastar, 
		       float costheta1, float costheta2, 
		       float phi, float phistar1,
		       float& kd, float& psig, float& psigALT){

  if(zzmass<100.0 || zzmass>179.9999){
    psig=0.0;
    psigALT=0.0;
    kd=0.0;
    return;
  }

  checkZorder(z1mass,z2mass,costhetastar,costheta1,costheta2,phi,phistar1);

  z1mass_rrv->setVal(z1mass);
  z2mass_rrv->setVal(z2mass);
  costhetastar_rrv->setVal(costhetastar);
  costheta1_rrv->setVal(costheta1);
  costheta2_rrv->setVal(costheta2);
  phi_rrv->setVal(phi);
  phistar1_rrv->setVal(phistar1);

  mzz_rrv->setVal(zzmass);

  psig = SMHiggs->PDF->getVal();
  psigALT = sigAlt->PDF->getVal();
  kd = 1/(1+psigALT/psig);
    
}


