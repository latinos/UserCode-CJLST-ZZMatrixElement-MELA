#ifndef AngPdfFact_h
#define AngPdfFact_h

#include "RooRealVar.h"
#include "RooXZsZs_5D.h"
#include "ZZMatrixElement/MELA/interface/TVar.hh"
#include <cmath>
#include <string>

class AngularPdfFactory{

public:

  RooRealVar* mZ;     
  RooRealVar* mX;     
  RooRealVar* gamZ;   
    
  RooRealVar* a1Val;  
  RooRealVar* phi1Val;
  RooRealVar* a2Val;  
  RooRealVar* phi2Val;
  RooRealVar* a3Val;  
  RooRealVar* phi3Val;
    
  RooRealVar* R1Val;  
  RooRealVar* R2Val;  
  
  RooXZsZs_5D *PDF;

  int modelIndex;  //0 - SM Higgs, 1 - PS Higgs, 2 - Fully Longitudinal Scalar, -1 - Custom


  AngularPdfFactory(){};

  AngularPdfFactory(RooRealVar* m1,RooRealVar* m2,RooRealVar* h1,RooRealVar* h2,RooRealVar* Phi,RooRealVar* mZZ){

    // Parameters
    mZ     = new RooRealVar("mZ","mZ",91.188);
    gamZ   = new RooRealVar("gamZ","gamZ",2.5);
           
    a1Val  = new RooRealVar("a1Val","a1Val",1);
    phi1Val= new RooRealVar("phi1Val","phi1Val",0);
    a2Val  = new RooRealVar("a2Val","a2Val",0);
    phi2Val= new RooRealVar("phi2Val","phi2Val",0);
    a3Val  = new RooRealVar("a3Val","a3Val",0);
    phi3Val= new RooRealVar("phi3Val","phi3Val",0);
           
    R1Val  = new RooRealVar("R1Val","R1Val",0.15);
    R2Val  = new RooRealVar("R2Val","R2Val",0.15);

    PDF = new RooXZsZs_5D("PDF","PDF",*m1,*m2,*h1,*h2,*Phi,*a1Val,*phi1Val,*a2Val,*phi2Val,*a3Val,*phi3Val,*mZ,*gamZ,*mZZ,*R1Val,*R2Val);

  };

  ~AngularPdfFactory(){

    delete mZ;
    delete gamZ;
    delete a1Val;
    delete phi1Val;
    delete a2Val;
    delete phi2Val;
    delete a3Val;
    delete phi3Val;
    delete R1Val;
    delete R2Val;
    
    delete PDF;

  };

  int configure(TVar::Process model_){

    switch (model_){
    case TVar::HZZ_4l: makeSMHiggs(); return 0; break;
    case TVar::HDHZZ_4l : makeLGHiggs(); return 0; break;
    case TVar::PSHZZ_4l : makePSHiggs(); return 0; break;
    default: makeSMHiggs(); return 1; break;
    }


  };

  void makeSMHiggs(){
    a1Val->setVal(1.0);
    phi1Val->setVal(0.0);
    a2Val->setVal(0.0);
    phi2Val->setVal(0.0);
    a3Val->setVal(0.0);
    phi3Val->setVal(0.0);
    modelIndex=0;
  };

  void makeLGHiggs(){          
    a1Val->setVal(0.0);        // need to calculate the proper normalizations
    phi1Val->setVal(0.0);
    a2Val->setVal(1.0);
    phi2Val->setVal(0.0);
    a3Val->setVal(0.0);
    phi3Val->setVal(0.0);
    modelIndex=2;
  };


  void makePSHiggs(){
    a1Val->setVal(0.0);
    phi1Val->setVal(0.0);
    a2Val->setVal(0.0);
    phi2Val->setVal(0.0);
    a3Val->setVal(1.0);
    phi3Val->setVal(0.0);
    modelIndex=1;
  };

  void makeCustom(double a1, double a2, double a3,
		  double phi1, double phi2, double phi3){
    a1Val->setVal(a1);
    phi1Val->setVal(phi1);
    a2Val->setVal(a2);
    phi2Val->setVal(phi2);
    a3Val->setVal(a3);
    phi3Val->setVal(phi3);
    modelIndex=-1;
  };

  void makeParamsConst(bool yesNo=true){
    if(yesNo){
      a1Val->setConstant(kTRUE);
      phi1Val->setConstant(kTRUE);
      a2Val->setConstant(kTRUE);
      phi2Val->setConstant(kTRUE);
      a3Val->setConstant(kTRUE);
      phi3Val->setConstant(kTRUE);
      gamZ->setConstant(kTRUE);
      mZ->setConstant(kTRUE);
      R1Val->setConstant(kTRUE);
      R2Val->setConstant(kTRUE);
    }else{
      a1Val->setConstant(kFALSE);
      phi1Val->setConstant(kFALSE);
      a2Val->setConstant(kFALSE);
      phi2Val->setConstant(kFALSE);
      a3Val->setConstant(kFALSE);
      phi3Val->setConstant(kFALSE);
      gamZ->setConstant(kFALSE);
      mZ->setConstant(kFALSE);
      R1Val->setConstant(kFALSE);
      R2Val->setConstant(kFALSE);
    }
  };

  double getVal(double mZZ){

    double Norm[80];
    
    if(modelIndex==0){  // SMHiggs


      Norm[0 ]= 4*3.1415*0.715995;
      Norm[1 ]= 4*3.1415*0.807688;
      Norm[2 ]= 4*3.1415*0.914359;
      Norm[3 ]= 4*3.1415*1.04167;
      Norm[4 ]= 4*3.1415*1.19928;
      Norm[5 ]= 4*3.1415*1.39749;
      Norm[6 ]= 4*3.1415*1.64484;
      Norm[7 ]= 4*3.1415*1.95036;
      Norm[8 ]= 4*3.1415*2.32449;
      Norm[9 ]= 4*3.1415*2.77916;
      Norm[10]= 4*3.1415*3.9858;
      Norm[11]= 4*3.1415*3.32788;
      Norm[12]= 4*3.1415*4.76983;
      Norm[13]= 4*3.1415*5.69869;
      Norm[14]= 4*3.1415*6.79307;
      Norm[15]= 4*3.1415*8.07577;
      Norm[16]= 4*3.1415*9.57179;
      Norm[17]= 4*3.1415*11.3085;
      Norm[18]= 4*3.1415*13.3159;
      Norm[19]= 4*3.1415*15.6266;
      Norm[20]= 4*3.1415*18.2761;
      Norm[21]= 4*3.1415*21.3032;
      Norm[22]= 4*3.1415*24.7498;
      Norm[23]= 4*3.1415*28.6616;
      Norm[24]= 4*3.1415*33.0881;
      Norm[25]= 4*3.1415*38.0828;
      Norm[26]= 4*3.1415*43.704;
      Norm[27]= 4*3.1415*50.0145;
      Norm[28]= 4*3.1415*57.0824;
      Norm[29]= 4*3.1415*64.9817;
      Norm[30]= 4*3.1415*73.7921;
      Norm[31]= 4*3.1415*83.6001;
      Norm[32]= 4*3.1415*94.4997;
      Norm[33]= 4*3.1415*106.592;
      Norm[34]= 4*3.1415*119.988;
      Norm[35]= 4*3.1415*134.806;
      Norm[36]= 4*3.1415*151.177;
      Norm[37]= 4*3.1415*169.24;
      Norm[38]= 4*3.1415*189.15;
      Norm[39]= 4*3.1415*211.073;
      Norm[40]= 4*3.1415*235.19;
      Norm[41]= 4*3.1415*261.701;
      Norm[42]= 4*3.1415*290.82;
      Norm[43]= 4*3.1415*322.785;
      Norm[44]= 4*3.1415*357.855;
      Norm[45]= 4*3.1415*396.316;
      Norm[46]= 4*3.1415*438.479;
      Norm[47]= 4*3.1415*484.69;
      Norm[48]= 4*3.1415*535.329;
      Norm[49]= 4*3.1415*590.819;
      Norm[50]= 4*3.1415*651.625;
      Norm[51]= 4*3.1415*718.268;
      Norm[52]= 4*3.1415*791.328;
      Norm[53]= 4*3.1415*871.453;
      Norm[54]= 4*3.1415*959.373;
      Norm[55]= 4*3.1415*1055.91;
      Norm[56]= 4*3.1415*1161.98;
      Norm[57]= 4*3.1415*1278.65;
      Norm[58]= 4*3.1415*1407.12;
      Norm[59]= 4*3.1415*1548.77;
      Norm[60]= 4*3.1415*1705.19;
      Norm[61]= 4*3.1415*1878.21;
      Norm[62]= 4*3.1415*2069.98;
      Norm[63]= 4*3.1415*2283.03;
      Norm[64]= 4*3.1415*2520.35;
      Norm[65]= 4*3.1415*2785.48;
      Norm[66]= 4*3.1415*3082.73;
      Norm[67]= 4*3.1415*3417.3;
      Norm[68]= 4*3.1415*3795.62;
      Norm[69]= 4*3.1415*4225.66;
      Norm[70]= 4*3.1415*4717.54;
      Norm[71]= 4*3.1415*5284.21;
      Norm[72]= 4*3.1415*5942.62;
      Norm[73]= 4*3.1415*6715.39;
      Norm[74]= 4*3.1415*7633.45;
      Norm[75]= 4*3.1415*8740.24;
      Norm[76]= 4*3.1415*10098.8;
      Norm[77]= 4*3.1415*11803.7;
      Norm[78]= 4*3.1415*14002.8;
      Norm[79]= 4*3.1415*16935.5;
    }else{

      return PDF->getVal()/1.0;

    }

    if((int)floor(mZZ-100)>79){
      
      return PDF->getVal()/Norm[79];

    }if((int)floor(mZZ-100)<0){

      return PDF->getVal()/Norm[0];

    }

    return PDF->getVal()/Norm[(int)floor(mZZ-100)];

  };

  double getValIntegrOutAngles(RooRealVar* m1,RooRealVar* m2,RooRealVar* h1,RooRealVar* h2,RooRealVar* Phi,RooRealVar* mZZ){
    RooAbsPdf* PDFIntegratedOut =PDF->createProjection(RooArgSet(*h1,*h2,*Phi));
    double norm = PDFIntegratedOut->getNorm(RooArgSet(*m1, *m2, *mZZ));
    cout<<"norm "<<norm<<endl;
    double val = PDFIntegratedOut->getVal()/norm;
    cout<<"val "<<val<<endl;
   return val;
  }


};
#endif
