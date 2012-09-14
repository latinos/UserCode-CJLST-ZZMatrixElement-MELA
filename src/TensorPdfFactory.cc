#ifndef TENSOR_PDF_FACTORY
#define TENSOR_PDF_FACTORY

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

  TF1* fmZZNorm;

  RooRealVar* mZ;
  RooRealVar* gamZ;

  RooRealVar* R1Val;
  RooRealVar* R2Val;

  RooAbsPdf *PDF;

  int modelIndex;

  
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

    delete fmZZNorm;

    delete mZ;
    delete gamZ;

    delete R1Val;
    delete R2Val;

    delete PDF;


  };

  void makeMinGrav(){  // NEED TO CALCULATE NORMALIZATIONS

    fz1Val->setVal(0.0);
    fz2Val->setVal(1.0);

    g1Val->setVal(1.0);
    g5Val->setVal(1.0); 

    modelIndex=0;
 
    fmZZNorm=new TF1("fmZZNorm","exp([0]+[1]*x+[2]*x*x+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5) +[6]*pow(x,6))",100,180);
    fmZZNorm->FixParameter(0,-5.22866);
    fmZZNorm->FixParameter(1,0.170495);
    fmZZNorm->FixParameter(2,0.000877867);
    fmZZNorm->FixParameter(3,-4.01876e-07);
    fmZZNorm->FixParameter(4,-3.02114e-08);
    fmZZNorm->FixParameter(5,-1.43247e-10);
    fmZZNorm->FixParameter(6,1.04256e-12);
  };

  void make2hPlus(){  // NEED TO CALCULATE NORMALIZATIONS

    fz1Val->setVal(0.0);
    fz2Val->setVal(0.0);

    g4Val->setVal(1.0); 

    modelIndex=1;
    fmZZNorm=new TF1("fmZZNorm","exp([0]+[1]*x+[2]*x*x+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)+[6]*pow(x,6))",100,180);
    fmZZNorm->FixParameter(0,-8.59286);
    fmZZNorm->FixParameter(1,0.0337467);
    fmZZNorm->FixParameter(2,0.0004704);
    fmZZNorm->FixParameter(3,2.15492e-06);
    fmZZNorm->FixParameter(4,1.19595e-09);
    fmZZNorm->FixParameter(5,-4.19439e-11);
    fmZZNorm->FixParameter(6,-1.02067e-13);  
  };

  void make2hMinus(){  // NEED TO CALCULATE NORMALIZATIONS

    fz1Val->setVal(0.0);
    fz2Val->setVal(0.0);

    g8Val->setVal(1.0); 

    modelIndex=2;
    fmZZNorm=new TF1("fmZZNorm","exp([0]+[1]*x+[2]*x*x+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5))",100,180);
    
    fmZZNorm->FixParameter(0,-7.27873);
    fmZZNorm->FixParameter(1,0.0158159);
    fmZZNorm->FixParameter(2,0.000415876);
    fmZZNorm->FixParameter(3,2.86406e-06);
    fmZZNorm->FixParameter(4,6.2603e-09);
    fmZZNorm->FixParameter(5,-9.36201e-11);
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

  double getVal(double mZZ){
 
    double norm=-999;
    if(mZZ>180 || mZZ<100){
      cout << "Normalization is not available for this value of mZZ: I'm extrapolating ..." << mZZ << endl;
    }
    norm = fmZZNorm->Eval(mZZ);
    return PDF->getVal()/norm;
  };

};

#endif


