#ifndef MELA_newMELA_h
#define MELA_newMELA_h

#include "TLorentzVector.h"
#include <vector>


class TFile; 
class TH1F; 
class TH2F;
class TH3F;
class RooRealVar;
class RooAbsPdf;
class RooArgSet;
class AngularPdfFactory;
class VectorPdfFactory;
class TensorPdfFactory;
class newZZMatrixElement;
class TGraph;


#include <ZZMatrixElement/MELA/interface/TVar.hh>
#include <ZZMatrixElement/MELA/interface/TEvtProb.hh>
#include <ZZMatrixElement/MELA/src/AngularPdfFactory.h>
#include <ZZMatrixElement/MELA/src/VectorPdfFactory.h>
#include <ZZMatrixElement/MELA/src/TensorPdfFactory.h>


class newMELA{

public:
  
  // newMELA(){};
  newMELA(int LHCsqrts=8, float mh=126); // higgs mass for supermela
  ~newMELA();
  
  void setProcess(TVar::Process myModel, TVar::MatrixElement myME, TVar::Production myProduction);
  
  void computeP(float mZZ, float mZ1, float mZ2, // input kinematics
		float costhetastar,
		float costheta1, 
		float costheta2,
		float phi,
		float phi1,
		int flavor,
		float& prob);                   // output probability
    
  void computeP(TLorentzVector Z1_lept1, int Z1_lept1Id,  // input 4-vectors
		TLorentzVector Z1_lept2, int Z1_lept2Id,  // 
		TLorentzVector Z2_lept1, int Z2_lept1Id,
		TLorentzVector Z2_lept2, int Z2_lept2Id,  
		float& prob);                             // output probability
    
  // Ordering of Z1/Z2 according to internal convention
  void checkZorder(float& z1mass, float& z2mass, float& costhetastar, float& costheta1, float& costheta2, float& phi, float& phistar1);
  
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


 private:

  // 
  // data memmbers 
  // 
  bool usePowhegTemplate_;
  int LHCsqrts;
  TVar::Process myModel_;
  TVar::MatrixElement myME_;
  TVar::Production myProduction_;
  newZZMatrixElement* ZZME;
  
  
  // 
  // functions 
  // 
  

};

#endif

