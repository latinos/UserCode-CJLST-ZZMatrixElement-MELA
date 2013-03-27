#ifndef VECTOR_PDF_FACTORY
#define VECTOR_PDF_FACTORY

#include "ZZMatrixElement/MELA/interface/TVar.hh"
#include "RooSpinOne_7D.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include <cmath>

#include "TF1.h"

class VectorPdfFactory{

public:
    
  RooRealVar* mZ;
  RooRealVar* gamZ;

  RooRealVar* R1Val;
  RooRealVar* R2Val;

  RooAbsPdf *PDF;

  int modelIndex;

  RooRealVar* g1Val;
  RooRealVar* g2Val;

  RooRealVar* aParam;
  
  TF1* fmZZNorm;

  VectorPdfFactory(){};

  VectorPdfFactory(RooRealVar* m1,RooRealVar* m2,RooRealVar* hs,RooRealVar* h1,RooRealVar* h2,RooRealVar* Phi,RooRealVar* Phi1,RooRealVar* mZZ){

    // Parameters
    mZ     = new RooRealVar("mZ","mZ",91.188);
    gamZ   = new RooRealVar("gamZ","gamZ",2.5);

    // related to tensor structure of V decays
    R1Val  = new RooRealVar("R1Val","R1Val",0.15);
    R2Val  = new RooRealVar("R2Val","R2Val",0.15);

    // dimensionless couplings
    g1Val = new RooRealVar("g1Val", "g1Val", 0.0);        
    g2Val = new RooRealVar("g2Val", "g2Val", 0.0);

    // random paramter (?)
    aParam = new RooRealVar("aParam","aParam",0.0);

    PDF = new RooSpinOne_7D("PDF","PDF", *mZZ, *m1, *m2, *h1, *h2,*hs, *Phi, *Phi1, 
			    *g1Val, *g2Val, *R1Val, *R2Val, *aParam, *mZ, *gamZ);

  };

  ~VectorPdfFactory(){

    delete g1Val;
    delete g2Val; 

    delete aParam;

    delete fmZZNorm;

    delete mZ;
    delete gamZ;

    delete R1Val;
    delete R2Val;

    delete PDF;

  };

  int configure(TVar::Process model_){

    switch (model_){
    case TVar::AVZZ_4l: makePseudoZprime(); return 0; break;
    case TVar::VZZ_4l : makeZprime(); return 0; break;
    default: makeZprime(); return 1; break;
    }

  };


  void makePseudoZprime(){  // NEED TO CALCULATE NORMALIZATIONS

    g1Val->setVal(0.0);
    g2Val->setVal(1.0); 

    modelIndex=0;

    fmZZNorm=new TF1("fmZZNorm","exp([0]+[1]*x+[2]*x*x+[3]*pow(x,3)+[4]*pow(x,4))",100,180);
    fmZZNorm->FixParameter(0,-48.5383);
    fmZZNorm->FixParameter(1,0.703812);
    fmZZNorm->FixParameter(2,0.000160962);
    fmZZNorm->FixParameter(3,-2.62737e-05);
    fmZZNorm->FixParameter(4,8.13202e-08);
  };

  void makeZprime(){  // NEED TO CALCULATE NORMALIZATIONS

    g1Val->setVal(1.0); 
    g2Val->setVal(0.0); 

    modelIndex=1;
    fmZZNorm=new TF1("fmZZNorm","exp([0]+[1]*x+[2]*x*x+[3]*pow(x,3)+[4]*pow(x,4))",100,180);
    fmZZNorm->FixParameter(0,-33.6399);
    fmZZNorm->FixParameter(1,0.474581);
    fmZZNorm->FixParameter(2,0.000350871);
    fmZZNorm->FixParameter(3,-1.67939e-05);
    fmZZNorm->FixParameter(4,4.86108e-08);
  };

  void makeParamsConst(bool yesNo=true){
    if(yesNo){

      g1Val->setConstant(kTRUE);
      g2Val->setConstant(kTRUE);

      gamZ->setConstant(kTRUE);
      mZ->setConstant(kTRUE);
      R1Val->setConstant(kTRUE);
      R2Val->setConstant(kTRUE);

    }else{

      g1Val->setConstant(kFALSE);
      g2Val->setConstant(kFALSE);

      gamZ->setConstant(kFALSE);
      mZ->setConstant(kFALSE);
      R1Val->setConstant(kFALSE);
      R2Val->setConstant(kFALSE);
    }
  };

  double getVal(double mZZ){

    double norm=-999;
    if(mZZ>180 || mZZ<100){
      //cout << "Normalization is not available for this value of mZZ: I'm extrapolating ..." << mZZ << endl;
    }
    norm = fmZZNorm->Eval(mZZ);
    return PDF->getVal()/norm;
  };

};

#endif


