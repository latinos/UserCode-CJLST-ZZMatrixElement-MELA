#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZMatrixElement/MELA/interface/newZZMatrixElement.h>
#include <DataFormats/GeometryVector/interface/Pi.h>
#include <FWCore/ParameterSet/interface/FileInPath.h>

#include "computeAngles.h"
#include "AngularPdfFactory.h"
#include "VectorPdfFactory.h"
#include "TensorPdfFactory.h"
#include "RooqqZZ_JHU_ZgammaZZ_fast.h"
#include "RooqqZZ_JHU.h"
#include "ZZMatrixElement/MELA/interface/SuperMELA.h"

#include <RooMsgService.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TGraph.h>
#include <vector>

#include <string>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>


using namespace RooFit;

Mela::Mela(int LHCsqrts, float mh) 
{
  // Create symlinks to the required files, if these are not already present (do nothing otherwse)
  edm::FileInPath mcfmInput1("ZZMatrixElement/MELA/data/input.DAT");
  edm::FileInPath mcfmInput2("ZZMatrixElement/MELA/data/process.DAT");
  edm::FileInPath mcfmInput3("ZZMatrixElement/MELA/data/Pdfdata/cteq6l1.tbl");  
  symlink(mcfmInput1.fullPath().c_str(), "input.DAT");
  symlink(mcfmInput2.fullPath().c_str(), "process.DAT");
  mkdir("Pdfdata",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  symlink(mcfmInput3.fullPath().c_str(), "Pdfdata/cteq6l1.tbl");

  mzz_rrv = new RooRealVar("mzz","m_{ZZ}",0.,1000.);
  z1mass_rrv = new RooRealVar("z1mass","m_{Z1}",0.,180.);
  z2mass_rrv = new RooRealVar("z2mass","m_{Z2}",0.,120.); 
  costhetastar_rrv = new RooRealVar("costhetastar","cos#theta^{*}",-1.,1.);  
  costheta1_rrv = new RooRealVar("costheta1","cos#theta_{1}",-1.,1.);  
  costheta2_rrv = new RooRealVar("costheta2","cos#theta_{2}",-1.,1.);
  phi_rrv= new RooRealVar("phi","#Phi",-3.1415,3.1415);
  phi1_rrv= new RooRealVar("phi1","#Phi_{1}",-3.1415,3.1415);
  
  upFrac_rrv = new RooRealVar("upFrac","fraction up-quarks",.5,0.,1.);

  spin0Model = new AngularPdfFactory(z1mass_rrv,z2mass_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,mzz_rrv);
  spin1Model = new VectorPdfFactory(z1mass_rrv,z2mass_rrv,costhetastar_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,phi1_rrv,mzz_rrv);
  spin2Model = new TensorPdfFactory(z1mass_rrv,z2mass_rrv,costhetastar_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,phi1_rrv,mzz_rrv);
  qqZZmodel = new RooqqZZ_JHU_ZgammaZZ_fast("qqZZmodel","qqZZmodel",*z1mass_rrv,*z2mass_rrv,*costheta1_rrv,*costheta2_rrv,*phi_rrv,*costhetastar_rrv,*phi1_rrv,*mzz_rrv,*upFrac_rrv);

  edm::FileInPath HiggsWidthFile("Higgs/Higgs_CS_and_Width/txtFiles/8TeV-ggH.txt");
  std::string path = HiggsWidthFile.fullPath();
  //std::cout << path.substr(0,path.length()-12) << std::endl;
  ZZME = new  newZZMatrixElement(path.substr(0,path.length()-12 ).c_str(),1000.*LHCsqrts/2.,TVar::INFO);

  // 
  // configure the JHUGEn and MCFM calculations 
  // 
  

  // load TGraphs for VAMCFM scale factors
  edm::FileInPath ScaleFactorFile("ZZMatrixElement/MELA/data/scaleFactors.root");
  TFile* sf = TFile::Open(ScaleFactorFile.fullPath().c_str(),"r");
  vaScale_4e    = (TGraph*)sf->Get("scaleFactors_4e");
  vaScale_4mu   = (TGraph*)sf->Get("scaleFactors_4mu");
  vaScale_2e2mu = (TGraph*)sf->Get("scaleFactors_2e2mu");
  sf->Close(); 


  //
  // setup supermela
  //

  //deactivate generation messages
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  RooMsgService::instance().setStreamStatus(1,kFALSE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);// silence also the error messages, but should really be looked at.


  myR=new TRandom3(35797);
  //  std::cout << "before supermela" << std::endl;
  super = new SuperMELA(mh,"4mu",LHCsqrts); // preliminary intialization, we adjust the flavor later
  char cardpath[500];
  sprintf(cardpath,"HZZ4L_Combination/CombinationPy/CreateDatacards/SM_inputs_%dTeV/inputs_4mu.txt",LHCsqrts);
  //std::cout << "before supermela, pathToCards: " <<cardpath<< std::endl;
  edm::FileInPath cardfile(cardpath);
  std::string cpath=cardfile.fullPath();
  std::cout << cpath.substr(0,cpath.length()-14).c_str()  <<std::endl;
  super->SetPathToCards(cpath.substr(0,cpath.length()-14).c_str() );
  super->SetVerbosity(false);
  // std::cout << "starting superMELA initialization" << std::endl;
  super->init();
  //std::cout << "after supermela" << std::endl;


  
}


Mela::~Mela(){ 
  //std::cout << "begin destructor" << std::endl;  
  delete mzz_rrv;
  delete z1mass_rrv; 
  delete z2mass_rrv; 
  delete costhetastar_rrv;
  delete costheta1_rrv;
  delete costheta2_rrv;
  delete phi_rrv;
  delete phi1_rrv;
  delete upFrac_rrv;

  delete spin0Model;
  delete spin1Model;
  delete spin2Model;
  delete qqZZmodel;
  delete ZZME;
  delete super;
  delete myR;
}

void Mela::setProcess(TVar::Process myModel, TVar::MatrixElement myME, TVar::Production myProduction)
{
  myModel_ = myModel;
  myME_ = myME;
  myProduction_ = myProduction;

  // 
  // configure the analytical calculations 
  // 

  if(myModel_==TVar::ZZ_2e2m || myModel_==TVar::ZZ_4e)
    pdf = qqZZmodel;
  else if(!spin0Model->configure(myModel_)) pdf = spin0Model->PDF;
  else if(!spin1Model->configure(myModel_)) pdf = spin1Model->PDF;
  else if(!spin2Model->configure(myModel_,myProduction_)) pdf = spin2Model->PDF;
  else if(myME_ == TVar::ANALYTICAL)
    cout << "Mela::setProcess -> ERROR TVar::Process not found!!! " << myME_ << endl; 

}


// Re-order masses and angles as needed by likelihoodDiscriminant. 
// This follows a different convention than the usual Z1/Z2 definition!
void Mela::checkZorder(float& z1mass, float& z2mass,
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
    //cout<<"inverted"<<endl;
    z1mass=tempZ2mass;
    z2mass=tempZ1mass;
    costhetastar=-tempHs;
    costheta1=tempH2;
    costheta2=tempH1;
    phi=tempPhi;
    phistar1=-tempPhi1-tempPhi;
    if(phistar1>3.1415)
      phistar1=phistar1-2*Geom::pi();
    if(phistar1<-3.1415)
      phistar1=phistar1+2*Geom::pi();

  }else
    return;

}

void Mela::computeP(float mZZ, float mZ1, float mZ2, // input kinematics
		       float costhetastar,
		       float costheta1, 
		       float costheta2,
		       float phi,
		       float phi1,
		       int flavor, 
		       float& prob){                   // output probability
    
  costhetastar_rrv->setVal(costhetastar);
  costheta1_rrv->setVal(costheta1);
  costheta2_rrv->setVal(costheta2);
  phi_rrv->setVal(phi);
  phi1_rrv->setVal(phi1);
    
  z1mass_rrv->setVal(mZ1);
  z2mass_rrv->setVal(mZ2);
  mzz_rrv->setVal(mZZ);
  
  float constant = 1.;
  
  //
  // analytical calculations
  // 
  if ( myME_ == TVar::ANALYTICAL ) {
   
    if(mZZ>100.){

      if(myProduction_==TVar::INDEPENDENT){
	RooAbsPdf* integral = (RooAbsPdf*) pdf->createIntegral(RooArgSet(*costhetastar_rrv,*phi1_rrv));
	prob = integral->getVal();
	delete integral;
      }else{
	prob = pdf->getVal();
      }

    }else{
      prob = -99.0;
    }
  } 

  //
  // JHUGen or MCFM 
  //
  if ( myME_ == TVar::JHUGen || myME_ == TVar::MCFM ) {
    
    //initialize variables
    checkZorder(mZ1,mZ2,costhetastar,costheta1,costheta2,phi,phi1);
    ZZME->computeXS(mZZ,mZ1,mZ2,
		    costhetastar,costheta1,costheta2, 
		    phi, phi1, flavor,
		    myModel_, myME_,  myProduction_,  prob);

    // adding scale factors for MCMF calculation
    // -- taken from old code --
    
    if(flavor==1 && myME_ == TVar::MCFM){
      if(mZZ > 900)                   
	prob *= vaScale_4e->Eval(900.);
      else if (mZZ <  100 )
	prob *= vaScale_4e->Eval(100.);
      else
	prob *= vaScale_4e->Eval(mZZ);
    }

    if(flavor==2 && myME_ == TVar::MCFM){
      if(mZZ > 900)                   
	prob *= vaScale_4mu->Eval(900.);
      else if (mZZ <  100 )
	prob *= vaScale_4mu->Eval(100.);
      else
	prob *= vaScale_4mu->Eval(mZZ);
    }

    if(flavor==3 && myME_ == TVar::MCFM){
      if(mZZ > 900)                   
	prob *= vaScale_2e2mu->Eval(900.);
      else if (mZZ <  100 )
	prob *= vaScale_2e2mu->Eval(100.);
      else
	prob *= vaScale_2e2mu->Eval(mZZ);
    }

    // 
    // define the constants to be used on JHUGen
    // 
    
    // gg productions 
    if ( myME_ == TVar::JHUGen && myProduction_ == TVar::GG  ) {
      if ( flavor == 3 ) {
	if ( myModel_ == TVar::PSHZZ_4l )  constant = 6.0;
	if ( myModel_ == TVar::HDHZZ_4l )  constant = 2.1;
	if ( myModel_ == TVar::TZZ_4l )  constant = 0.6;
	if ( myModel_ == TVar::TZZ_2hplus_4l )  constant = 2.8;
	if ( myModel_ == TVar::PTZZ_2hminus_4l )  constant = 4.3;
	if ( myModel_ == TVar::TZZ_2bplus_4l )  constant = 1.0;
      }  else {
	if ( myModel_ == TVar::PSHZZ_4l )  constant = 7.0;
	if ( myModel_ == TVar::HDHZZ_4l )  constant = 2.3;
	if ( myModel_ == TVar::TZZ_4l )  constant = 1.4;
	if ( myModel_ == TVar::TZZ_2hplus_4l )  constant = 2.5;
	if ( myModel_ == TVar::PTZZ_2hminus_4l )  constant = 3.6;
	if ( myModel_ == TVar::TZZ_2bplus_4l )  constant = 1.2;
      }

    } 
    // qqb productions 
    if ( myME_ == TVar::JHUGen && myProduction_ == TVar::QQB  ) {
      if ( flavor == 3 ) {
	if ( myModel_ == TVar::VZZ_4l )  constant = 16;
	if ( myModel_ == TVar::AVZZ_4l )  constant = 13;
	if ( myModel_ == TVar::TZZ_4l )  constant = 13;
      } else {
	if ( myModel_ == TVar::VZZ_4l )  constant = 38;
	if ( myModel_ == TVar::AVZZ_4l )  constant = 28;
	if ( myModel_ == TVar::TZZ_4l )  constant = 30;
      }
    }
    // production independent calculations
    if ( myME_ == TVar::JHUGen && myProduction_ == TVar::INDEPENDENT  ) {
      if ( flavor == 3) {
	if ( myModel_ == TVar::VZZ_4l )  constant = 1.3e+10;
	if ( myModel_ == TVar::AVZZ_4l )  constant = 1.3e+10;
	if ( myModel_ == TVar::TZZ_4l )  constant = 1.6e+9;
      } else {
	if ( myModel_ == TVar::VZZ_4l )  constant = 1.6e+10;
	if ( myModel_ == TVar::AVZZ_4l )  constant = 1.4e+10;
	if ( myModel_ == TVar::TZZ_4l )  constant = 2.0e+9;
      }
    } 


    // 
    // define the constants to be used on analytical
    // 
    
    // gg productions 
    if ( myME_ == TVar::ANALYTICAL && myProduction_ == TVar::GG  ) {
      if ( flavor == 3 ) {
	if ( myModel_ == TVar::PSHZZ_4l )  constant = 6.4;
	if ( myModel_ == TVar::HDHZZ_4l )  constant = 2.2;
	if ( myModel_ == TVar::TZZ_4l )  constant = 9.5;
	if ( myModel_ == TVar::TZZ_2hplus_4l )  constant = 7.3e7;
	if ( myModel_ == TVar::PTZZ_2hminus_4l )  constant = 1.1e8;
	if ( myModel_ == TVar::TZZ_2bplus_4l )  constant = 16.3;
      }  else {
	if ( myModel_ == TVar::PSHZZ_4l )  constant = 6.5;
	if ( myModel_ == TVar::HDHZZ_4l )  constant = 2.2;
	if ( myModel_ == TVar::TZZ_4l )  constant = 9.3;
	if ( myModel_ == TVar::TZZ_2hplus_4l )  constant = 1.1e8;
	if ( myModel_ == TVar::PTZZ_2hminus_4l )  constant = 1.9e8;
	if ( myModel_ == TVar::TZZ_2bplus_4l )  constant = 15.6;
      }

    } 
    // qqb productions 
    if ( myME_ == TVar::ANALYTICAL && myProduction_ == TVar::QQB  ) {
      if ( flavor == 3 ) {
	if ( myModel_ == TVar::VZZ_4l )  constant = 4.6e5;
	if ( myModel_ == TVar::AVZZ_4l )  constant = 4.0e5;
	if ( myModel_ == TVar::TZZ_4l )  constant = 7.9;
      } else {
	if ( myModel_ == TVar::VZZ_4l )  constant = 4.6e5;
	if ( myModel_ == TVar::AVZZ_4l )  constant = 4.0e5;
	if ( myModel_ == TVar::TZZ_4l )  constant = 7.9;
      }
    }
    // production independent calculations
    if ( myME_ == TVar::ANALYTICAL && myProduction_ == TVar::INDEPENDENT  ) {
      if ( flavor == 3) {
	if ( myModel_ == TVar::VZZ_4l )  constant = 3.4e4;
	if ( myModel_ == TVar::AVZZ_4l )  constant = 3.4e4;
	if ( myModel_ == TVar::TZZ_4l )  constant = 0.66;
      } else {
	if ( myModel_ == TVar::VZZ_4l )  constant = 3.4e4;
	if ( myModel_ == TVar::AVZZ_4l )  constant = 3.4e4;
	if ( myModel_ == TVar::TZZ_4l )  constant = .66;
      }
    } 
    
    // ***
    // experimental for the ZZ decay 
    // ****
    
    if ( myME_ == TVar::MCFM 
	 && myProduction_ == TVar::INDEPENDENT 
	 && ( myModel_ == TVar::ZZ_2e2m || myModel_ == TVar::ZZ_4e )
	 )
      {
	prob = 0.;
	int gridsize_hs = 10; 
	double hs_min = -1.;
	double hs_max = 1.;
	double hs_step =( hs_max - hs_min ) / double (gridsize_hs); 
	
	int gridsize_phi1 = 10; 
	double phi1_min = -TMath::Pi();
	double phi1_max = TMath::Pi();
	double phi1_step =( phi1_max - phi1_min ) / double (gridsize_phi1); 
	
	for ( int i_hs = 0; i_hs < gridsize_hs + 1; i_hs ++ ) {
	  double hs_val = hs_min + i_hs * hs_step; 
	  for ( int i_phi1 = 0; i_phi1 < gridsize_phi1 +1 ; i_phi1 ++ ) {
	    double phi1_val = phi1_min + i_phi1 * phi1_step; 
	    float temp_prob(0.); 
	    // calculate the ZZ using MCFM
	    ZZME->computeXS(mZZ,mZ1,mZ2,
			    hs_val,costheta1,costheta2, 
			    phi, phi1_val, flavor,
			    myModel_, myME_,  myProduction_,  temp_prob);
	    prob += temp_prob;
	  }
	}
	prob =  prob / float ( (gridsize_hs + 1) * (gridsize_phi1 +1 )); 
      }
  }
  prob *= constant; 
}

void Mela::computeP(TLorentzVector Z1_lept1, int Z1_lept1Id,  // input 4-vectors
		       TLorentzVector Z1_lept2, int Z1_lept2Id,  // 
		       TLorentzVector Z2_lept1, int Z2_lept1Id,
		       TLorentzVector Z2_lept2, int Z2_lept2Id,  
		       float& prob){                             // output probability
    
  // get flavor type
  // NEED TO INCLUDE SOME PROTECTION SO THAT USER CANT                  
  // PASS FOUR-VECTORS IN WRONG ORDER.  FOR NOW ASSUMING                
  // THEY ARE PASSED AS e-,e+,mu-,mu+                                   
  // ------------------ channel ------------------------                
  int flavor;

  if(abs(Z1_lept1Id)==abs(Z1_lept2Id)&&
     abs(Z1_lept1Id)==abs(Z2_lept1Id)&&
     abs(Z1_lept1Id)==abs(Z2_lept2Id)){

    if(abs(Z1_lept1Id)==11) flavor=1;
    else flavor=2;

  }else flavor=3;

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

  computeP(mzz, m1, m2,
	   costhetastar,
	   costheta1,
	   costheta2,
	   phi,phi1,
	   flavor,
	   prob);

}

void Mela::computePM4l(TLorentzVector Z1_lept1, int Z1_lept1Id,  // input 4-vectors
		       TLorentzVector Z1_lept2, int Z1_lept2Id,  // 
		       TLorentzVector Z2_lept1, int Z2_lept1Id,
		       TLorentzVector Z2_lept2, int Z2_lept2Id,
		       TVar::SuperMelaSyst syst, 
		       float& prob){

  TLorentzVector ZZ = (Z1_lept1 + Z1_lept2 + Z2_lept1 + Z2_lept2);
  float mzz = ZZ.M();
  TVar::LeptonFlavor flavor = TVar::Flavor_Dummy;

  if( abs(Z1_lept1Id)==11 &&  abs(Z1_lept2Id)==11 &&
      abs(Z2_lept1Id)==11 &&  abs(Z1_lept2Id)==11 )
    flavor = TVar::Flavor_4e;
  
  if( abs(Z1_lept1Id)==13 &&  abs(Z1_lept2Id)==13 &&
      abs(Z2_lept1Id)==13 &&  abs(Z1_lept2Id)==13 )
    flavor = TVar::Flavor_4mu;

  if( abs(Z1_lept1Id)==11 &&  abs(Z1_lept2Id)==11 &&
      abs(Z2_lept1Id)==13 &&  abs(Z1_lept2Id)==13 )
    flavor = TVar::Flavor_2e2mu;

  if( abs(Z1_lept1Id)==13 &&  abs(Z1_lept2Id)==13 &&
      abs(Z2_lept1Id)==11 &&  abs(Z1_lept2Id)==11 )
    flavor = TVar::Flavor_2e2mu;
  
  
  computePM4l(mzz,flavor,syst,prob);
}

void Mela::computePM4l(float mZZ, TVar::LeptonFlavor flavor, TVar::SuperMelaSyst syst, float& prob){
  prob=-99;//default dummy.
  
  if(flavor == TVar::Flavor_Dummy) // only compute things if flavor determination succeded
    return;
  
  switch(flavor){
  case 1: super->SetDecayChannel("4e")   ;break;
  case 2: super->SetDecayChannel("4mu")  ;break;
  case 3: super->SetDecayChannel("2e2mu");break;
  default: std::cout << " unknown flavor: " << flavor << std::endl; exit(0);
  }


  if(syst == TVar::SMSyst_None){
    std::pair<double,double> m4lP = super->M4lProb(mZZ);
    if(myModel_ == TVar::HZZ_4l) // currently only supported signal is H(0+)
      prob = m4lP.first;
    if(myModel_ == TVar::SummedBackgrounds) // currently only supported background is summed paramterization
      prob = m4lP.second;
  }
  else{
    //systematics for p(m4l)
    float mZZtmp=mZZ;
    float meanErr=float(super->GetSigShapeSystematic("meanCB") );
    if( syst == TVar::SMSyst_ScaleUp ){
      mZZtmp = mZZ*(1.0+meanErr);
      if(mZZtmp>180.0 || mZZtmp<100)mZZtmp=mZZ;      
      std::pair<double,double> m4lPScaleUp = super->M4lProb(mZZtmp);
      if(myModel_ == TVar::HZZ_4l)
	prob = m4lPScaleUp.first; 
      if(myModel_ == TVar::SummedBackgrounds)
	prob = m4lPScaleUp.second;
    }

    if( syst == TVar::SMSyst_ScaleDown ){    
      mZZtmp = mZZ*(1.0-meanErr);
      if(mZZtmp>180.0 || mZZtmp<100)mZZtmp=mZZ;
      std::pair<double,double> m4lPScaleDown = super->M4lProb(mZZtmp);
      if(myModel_ == TVar::HZZ_4l)
	prob = m4lPScaleDown.first; 
      if(myModel_ == TVar::SummedBackgrounds)
	prob = m4lPScaleDown.second;
    }
    
    float sigmaErr=float(super->GetSigShapeSystematic("sigmaCB") );
    float sigmaCB=float(super->GetSigShapeParameter("sigmaCB"));    
    if( syst == TVar::SMSyst_ResUp || syst ==  TVar::SMSyst_ResDown ){
      mZZtmp= myR->Gaus(mZZ,sigmaErr*sigmaCB);
      if(mZZtmp>180.0 || mZZtmp<100) mZZtmp=mZZ;
      std::pair<double,double> m4lPResUp = super->M4lProb(mZZtmp);
      if(myModel_ == TVar::HZZ_4l)
	prob = m4lPResUp.first; 
      if(myModel_ == TVar::SummedBackgrounds)
	prob = m4lPResUp.second;      
    }
  }
}

void Mela::computeWeight(float mZZ, float mZ1, float mZ2, 
			 float costhetastar,
			 float costheta1, 
			 float costheta2,
			 float phi,
			 float phi1,
			 // return variables:
			  float& w
			 ){

  float dXsec_HZZ_JHU,dXsec_HZZ_JHU_interf; // temporary prob
  
  // calculate dXsec for 4e/4mu
  setProcess(TVar::HZZ_4l,TVar::JHUGen,TVar::GG);
  computeP(mZZ,mZ1,mZ2,
	   costhetastar,
	   costheta1,
	   costheta2,
	   phi,phi1,
	   1,dXsec_HZZ_JHU_interf);

  // calculate dXsec for 2e2mu
  setProcess(TVar::HZZ_4l,TVar::JHUGen,TVar::GG);
  computeP(mZZ,mZ1,mZ2,
	   costhetastar,
	   costheta1,
	   costheta2,
	   phi,phi1,
	   3,dXsec_HZZ_JHU);
  
  w = dXsec_HZZ_JHU_interf / dXsec_HZZ_JHU;

  // protect against anomalously large weights
  if (w>10.) w=0.;

}


