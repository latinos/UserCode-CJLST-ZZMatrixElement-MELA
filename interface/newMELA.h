#ifndef NEWMELA
#define NEWMELA

#include "RooAbsPdf.h"
#include "TLorentzVector.h"
#include "ZZMatrixElement/MELA/interface/TVar.hh"
#include "ZZMatrixElement/MELA/src/AngularPdfFactory.h"
#include "ZZMatrixElement/MELA/src/VectorPdfFactory.h"
#include "ZZMatrixElement/MELA/src/TensorPdfFactory.h"

class newMELA{

public:
  
  TVar::Process myModel;
  RooAbsPdf* pdf;
  AngularPdfFactory* spin0Model;
  VectorPdfFactory* spin1Model;
  TensorPdfFactory* spin2Model;

  RooRealVar* mzz_rrv;
  RooRealVar* z1mass_rrv;
  RooRealVar* z2mass_rrv;
  RooRealVar* costhetastar_rrv; 
  RooRealVar* costheta1_rrv;
  RooRealVar* costheta2_rrv;
  RooRealVar* phi_rrv;
  RooRealVar* phi1_rrv;
  RooRealVar* pt_rrv;

  newMELA(){};

  newMELA(TVar::Process myModel_);

  void computeP(float mZZ, float mZ1, float mZ2, // input kinematics
		float costhetastar,
		float costheta1, 
		float costheta2,
		float phi,
		float phi1,
		float& prob);                   // output probability
    
  void computeP(TLorentzVector Z1_lept1, int Z1_lept1Id,  // input 4-vectors
		TLorentzVector Z1_lept2, int Z1_lept2Id,  // 
		TLorentzVector Z2_lept1, int Z2_lept1Id,
		TLorentzVector Z2_lept2, int Z2_lept2Id,  
		float& prob);                             // output probability
    


};

#endif

