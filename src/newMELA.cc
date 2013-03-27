#include "ZZMatrixElement/MELA/interface/newMELA.h"
#include "ZZMatrixElement/MELA/src/computeAngles.h"

newMELA::newMELA(TVar::Process myModel_){
  
  mzz_rrv = new RooRealVar("mzz","m_{ZZ}",0.,1000.);
  z1mass_rrv = new RooRealVar("z1mass","m_{Z1}",0.,180.);
  z2mass_rrv = new RooRealVar("z2mass","m_{Z2}",0.,120.); 
  costhetastar_rrv = new RooRealVar("costhetastar","cos#theta^{*}",-1.,1.);  
  costheta1_rrv = new RooRealVar("costheta1","cos#theta_{1}",-1.,1.);  
  costheta2_rrv = new RooRealVar("costheta2","cos#theta_{2}",-1.,1.);
  phi_rrv= new RooRealVar("phi","#Phi",-3.1415,3.1415);
  phi1_rrv= new RooRealVar("phi1","#Phi_{1}",-3.1415,3.1415);

  spin0Model = new AngularPdfFactory(z1mass_rrv,z2mass_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,mzz_rrv);
  spin1Model = new VectorPdfFactory(z1mass_rrv,z2mass_rrv,costhetastar_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,phi1_rrv,mzz_rrv);
  spin2Model = new TensorPdfFactory(z1mass_rrv,z2mass_rrv,costhetastar_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,phi1_rrv,mzz_rrv);

  myModel = myModel_;
        
  if(!spin0Model->configure(myModel)) pdf = spin0Model->PDF;
  else{
    if(!spin1Model->configure(myModel)) pdf = spin1Model->PDF;
    else{
      if(!spin2Model->configure(myModel)) pdf = spin2Model->PDF;
      else{
	cout << "newMELA::newMELA ERROR model not found!!!" << endl; 
      }	
    }
  }
  
};

void newMELA::computeP(float mZZ, float mZ1, float mZ2, // input kinematics
		       float costhetastar,
		       float costheta1, 
		       float costheta2,
		       float phi,
		       float phi1,
		       float& prob){                   // output probability
    
  costhetastar_rrv->setVal(costhetastar);
  costheta1_rrv->setVal(costheta1);
  costheta2_rrv->setVal(costheta2);
  phi_rrv->setVal(phi);
  phi1_rrv->setVal(phi1);
    
  z1mass_rrv->setVal(mZ1);
  z2mass_rrv->setVal(mZ2);
  mzz_rrv->setVal(mZZ);

  if(mZZ>100.)
    prob = pdf->getVal();
  else 
    prob = -99.0;

}

void newMELA::computeP(TLorentzVector Z1_lept1, int Z1_lept1Id,  // input 4-vectors
		       TLorentzVector Z1_lept2, int Z1_lept2Id,  // 
		       TLorentzVector Z2_lept1, int Z2_lept1Id,
		       TLorentzVector Z2_lept2, int Z2_lept2Id,  
		       float& prob){                             // output probability
    
  //compute angles  
  float m1=(Z1_lept1 + Z1_lept2).M();
  float m2=(Z2_lept1 + Z2_lept2).M();
    
  TLorentzVector ZZ = (Z1_lept1 + Z1_lept2 + Z2_lept1 + Z2_lept2);
  float mzz = ZZ.M();
    
  // Skip candidates where KD is irrelevant.
  if (mzz<100.){
    prob = -99.0;
    return;
  }

  float costhetastar, costheta1, costheta2, phi, phi1;

  mela::computeAngles(Z1_lept1, Z1_lept1Id, Z1_lept2, Z1_lept2Id, 
		      Z2_lept1, Z2_lept1Id, Z2_lept2, Z2_lept2Id,
		      costhetastar,costheta1,costheta2,phi,phi1);

  costhetastar_rrv->setVal(costhetastar);
  costheta1_rrv->setVal(costheta1);
  costheta2_rrv->setVal(costheta2);
  phi_rrv->setVal(phi);
  phi1_rrv->setVal(phi1);
    
  z1mass_rrv->setVal(m1);
  z2mass_rrv->setVal(m2);
  mzz_rrv->setVal(mzz);

  if(mzz>100.)
    prob = pdf->getVal();
  else
    prob = -99.;

}

