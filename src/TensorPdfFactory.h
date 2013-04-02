#ifndef TENSOR_PDF_FACTORY
#define TENSOR_PDF_FACTORY

#include "ZZMatrixElement/MELA/interface/TVar.hh"
#include "RooSpinTwo_7D.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include <cmath>
#include "TF1.h"

class TensorPdfFactory{

public:
    
  RooRealVar* c1Val;
  RooRealVar* c2Val;
  RooRealVar* c3Val;
  RooRealVar* c4Val;
  RooRealVar* c5Val;
  RooRealVar* c6Val;
  RooRealVar* c7Val;
  RooRealVar* useGTerm;
  RooRealVar* g1Val;
  RooRealVar* g2Val;
  RooRealVar* g3Val;
  RooRealVar* g4Val;
  RooRealVar* g5Val;
  RooRealVar* g6Val;
  RooRealVar* g7Val;
  RooRealVar* g8Val;
  RooRealVar* g9Val;
  RooRealVar* g10Val;
  RooRealVar* fz1Val;
  RooRealVar* fz2Val;

  RooRealVar* mZ;
  RooRealVar* gamZ;

  RooRealVar* R1Val;
  RooRealVar* R2Val;

  RooAbsPdf *PDF;

  TensorPdfFactory(){};

  TensorPdfFactory(RooRealVar* m1,RooRealVar* m2,RooRealVar* hs,RooRealVar* h1,RooRealVar* h2,RooRealVar* Phi,RooRealVar* Phi1,RooRealVar* mZZ){

    // Parameters
    mZ     = new RooRealVar("mZ","mZ",91.188);
    gamZ   = new RooRealVar("gamZ","gamZ",2.5);

    // related to tensor structure of V decays
    R1Val  = new RooRealVar("R1Val","R1Val",0.15);
    R2Val  = new RooRealVar("R2Val","R2Val",0.15);

    // related to the gg/qq productions 
    fz1Val = new RooRealVar("fz1Val", "fz1Val", 0.);
    fz2Val = new RooRealVar("fz2Val", "fz2Val", 1.0);
           
    // minimal set of lorentz structures
    c1Val = new RooRealVar("c1Val", "c1Val", 0.0);
    c2Val = new RooRealVar("c2Val", "c2Val", 0.0);
    c3Val = new RooRealVar("c3Val", "c3Val", 0.0);
    c4Val = new RooRealVar("c4Val", "c4Val", 0.0);
    c5Val = new RooRealVar("c5Val", "c5Val", 0.0);
    c6Val = new RooRealVar("c6Val", "c6Val", 0.0);
    c7Val = new RooRealVar("c7Val", "c7Val", 0.0);

    useGTerm = new RooRealVar("useGTerm", "useGTerm",1.); // set to 1 if using g couplings
                                                          // set to -1 if using c couplings
    // dimensionless couplings
    g1Val = new RooRealVar("g1Val", "g1Val", 0.0);        
    g2Val = new RooRealVar("g2Val", "g2Val", 0.0);
    g3Val = new RooRealVar("g3Val", "g3Val", 0.0);
    g4Val = new RooRealVar("g4Val", "g4Val", 0.0);
    g5Val = new RooRealVar("g5Val", "g5Val", 0.0);
    g6Val = new RooRealVar("g6Val", "g6Val", 0.0);
    g7Val = new RooRealVar("g7Val", "g7Val", 0.0);
    g8Val = new RooRealVar("g8Val", "g8Val", 0.0);
    g9Val = new RooRealVar("g9Val", "g9Val", 0.0);
    g10Val = new RooRealVar("g10Val", "g10Val", 0.0);

    PDF = new RooSpinTwo_7D("PDF","PDF", *mZZ, *m1, *m2, *hs, *h1,*h2, *Phi, *Phi1, 
				  *c1Val, *c2Val, *c3Val, *c4Val, *c5Val, *c6Val, *c7Val, 
				  *useGTerm, *g1Val, *g2Val, *g3Val, *g4Val, *g5Val, *g6Val, *g7Val, *g8Val, *g9Val, *g10Val,
				  *fz1Val, *fz2Val, *R1Val, *R2Val, *mZ, *gamZ);

  };

  ~TensorPdfFactory(){

    delete fz1Val;
    delete fz2Val;

    delete c1Val; 
    delete c2Val; 
    delete c3Val; 
    delete c4Val; 
    delete c5Val; 
    delete c6Val; 
    delete c7Val; 

    delete useGTerm;
    delete g1Val;
    delete g2Val; 
    delete g3Val; 
    delete g4Val; 
    delete g5Val; 
    delete g6Val; 
    delete g7Val; 
    delete g8Val; 
    delete g9Val; 
    delete g10Val;

    delete mZ;
    delete gamZ;

    delete R1Val;
    delete R2Val;

    delete PDF;


  };

  int configure(TVar::Process model_,TVar::Production prod_){

    switch (model_){
    case TVar::TZZ_4l          : 
      if(prod_==TVar::GG || prod_==TVar::INDEPENDENT){
	makeMinGrav(); 
	return 0; break;
      }else if(prod_==TVar::QQB){
	makeqqMinGrav();
	return 0; break;
      }else{
	return 2; break;
      }
    case TVar::TZZ_2hplus_4l   : make2hPlus();  return 0; break;
    case TVar::PTZZ_2hminus_4l : make2hMinus(); return 0; break;
    case TVar::TZZ_2bplus_4l   : make2bPlus();  return 0; break;
    default: makeMinGrav(); return 1; break;
    }

  };

  void makeMinGrav(){      // Minimal coupling graviton produced through 
                           // gluon-gluon fusion
    fz1Val->setVal(0.0);
    fz2Val->setVal(1.0);

    g1Val->setVal(1.0); 
    g2Val->setVal(0.0); 
    g3Val->setVal(0.0); 
    g4Val->setVal(0.0); 
    g5Val->setVal(1.0); 
    g6Val->setVal(0.0); 
    g7Val->setVal(0.0); 
    g8Val->setVal(0.0); 

  };

  void makeqqMinGrav(){   // Minimal coupling graviton produced through 
                          // quark-anti-quark annihilation
    fz1Val->setVal(1.0);
    fz2Val->setVal(0.0);

    g1Val->setVal(1.0); 
    g2Val->setVal(0.0); 
    g3Val->setVal(0.0); 
    g4Val->setVal(0.0); 
    g5Val->setVal(1.0); 
    g6Val->setVal(0.0); 
    g7Val->setVal(0.0); 
    g8Val->setVal(0.0); 

  };

  void makeUnpolMinGrav(){  // unpolarized minimal coupling graviton

    fz1Val->setVal(0.4);
    fz2Val->setVal(0.4);

    g1Val->setVal(1.0); 
    g2Val->setVal(0.0); 
    g3Val->setVal(0.0); 
    g4Val->setVal(0.0); 
    g5Val->setVal(1.0); 
    g6Val->setVal(0.0); 
    g7Val->setVal(0.0); 
    g8Val->setVal(0.0); 

  };

  void make2hPlus(){      // exotic CP-even spin-2 resonance produced through 
                          // gluon-gluon fusion
    fz1Val->setVal(0.0);
    fz2Val->setVal(0.0);

    g1Val->setVal(0.0); 
    g2Val->setVal(0.0); 
    g3Val->setVal(0.0); 
    g4Val->setVal(1.0); 
    g5Val->setVal(0.0); 
    g6Val->setVal(0.0); 
    g7Val->setVal(0.0); 
    g8Val->setVal(0.0); 

  };

  void make2hMinus(){     // exotic CP-odd spin-2 resonance produced through 
                          // gluon-gluon fusion
    fz1Val->setVal(0.0);
    fz2Val->setVal(0.0);

    g1Val->setVal(0.0); 
    g2Val->setVal(0.0); 
    g3Val->setVal(0.0); 
    g4Val->setVal(0.0); 
    g5Val->setVal(0.0); 
    g6Val->setVal(0.0); 
    g7Val->setVal(0.0); 
    g8Val->setVal(1.0); 

  };

  void make2bPlus(){     // spin-2 with bulk
                          // gluon-gluon fusion
    fz1Val->setVal(0.0);
    fz2Val->setVal(1.0);

    g1Val->setVal(0.0); 
    g2Val->setVal(0.0); 
    g3Val->setVal(0.0); 
    g4Val->setVal(0.0); 
    g5Val->setVal(1.0); 
    g6Val->setVal(0.0); 
    g7Val->setVal(0.0); 
    g5Val->setVal(0.0); 

  };

  void makeParamsConst(bool yesNo=true){
    if(yesNo){
      c1Val->setConstant(kTRUE);
      c2Val->setConstant(kTRUE);
      c3Val->setConstant(kTRUE);
      c4Val->setConstant(kTRUE);
      c5Val->setConstant(kTRUE);
      c6Val->setConstant(kTRUE);
      c7Val->setConstant(kTRUE);

      useGTerm->setConstant(kTRUE);
      g1Val->setConstant(kTRUE);
      g2Val->setConstant(kTRUE);
      g3Val->setConstant(kTRUE);
      g4Val->setConstant(kTRUE);
      g5Val->setConstant(kTRUE);
      g6Val->setConstant(kTRUE);
      g7Val->setConstant(kTRUE);
      g8Val->setConstant(kTRUE);
      g9Val->setConstant(kTRUE);
      g10Val->setConstant(kTRUE);

      fz1Val->setConstant(kTRUE);
      fz2Val->setConstant(kTRUE);

      gamZ->setConstant(kTRUE);
      mZ->setConstant(kTRUE);
      R1Val->setConstant(kTRUE);
      R2Val->setConstant(kTRUE);

    }else{
      c1Val->setConstant(kFALSE);
      c2Val->setConstant(kFALSE);
      c3Val->setConstant(kFALSE);
      c4Val->setConstant(kFALSE);
      c5Val->setConstant(kFALSE);
      c6Val->setConstant(kFALSE);
      c7Val->setConstant(kFALSE);

      useGTerm->setConstant(kFALSE);
      g1Val->setConstant(kFALSE);
      g2Val->setConstant(kFALSE);
      g3Val->setConstant(kFALSE);
      g4Val->setConstant(kFALSE);
      g5Val->setConstant(kFALSE);
      g6Val->setConstant(kFALSE);
      g7Val->setConstant(kFALSE);
      g8Val->setConstant(kFALSE);
      g9Val->setConstant(kFALSE);
      g10Val->setConstant(kFALSE);

      fz1Val->setConstant(kFALSE);
      fz2Val->setConstant(kFALSE);

      gamZ->setConstant(kFALSE);
      mZ->setConstant(kFALSE);
      R1Val->setConstant(kFALSE);
      R2Val->setConstant(kFALSE);
    }
  };

};

#endif


