#include <ZZMatrixElement/MELA/interface/newMELA.h>
#include <ZZMatrixElement/MELA/interface/newZZMatrixElement.h>
#include <DataFormats/GeometryVector/interface/Pi.h>
#include <FWCore/ParameterSet/interface/FileInPath.h>

#include "computeAngles.h"
#include "AngularPdfFactory.h"
#include "VectorPdfFactory.h"
#include "TensorPdfFactory.h"
#include "RooqqZZ_JHU_ZgammaZZ_fast.h"
#include "RooqqZZ_JHU.h"

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

newMELA::newMELA(int LHCsqrts, float mh) 
{
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

  edm::FileInPath HiggsWidthFile("Higgs/Higgs_CS_and_Width/txtFiles/8TeV-ggH.txt");
  std::string path = HiggsWidthFile.fullPath();
  //std::cout << path.substr(0,path.length()-12) << std::endl;
  ZZME = new  newZZMatrixElement(path.substr(0,path.length()-12 ).c_str(),1000.*LHCsqrts/2.,TVar::INFO);

  // 
  // configure the JHUGEn and MCFM calculations 
  // 
  
  // Create symlinks to the required files, if these are not already present (do nothing otherwse)
  edm::FileInPath mcfmInput1("ZZMatrixElement/MELA/data/input.DAT");
  edm::FileInPath mcfmInput2("ZZMatrixElement/MELA/data/process.DAT");
  edm::FileInPath mcfmInput3("ZZMatrixElement/MELA/data/Pdfdata/cteq6l1.tbl");  
  symlink(mcfmInput1.fullPath().c_str(), "input.DAT");
  symlink(mcfmInput2.fullPath().c_str(), "process.DAT");
  mkdir("Pdfdata",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  symlink(mcfmInput3.fullPath().c_str(), "Pdfdata/cteq6l1.tbl");
  
}


newMELA::~newMELA(){ 
  //std::cout << "begin destructor" << std::endl;  
  delete mzz_rrv;
  delete z1mass_rrv; 
  delete z2mass_rrv; 
  delete costhetastar_rrv;
  delete costheta1_rrv;
  delete costheta2_rrv;
  delete phi_rrv;
  delete phi1_rrv;
  delete spin0Model;
  delete spin1Model;
  delete spin2Model;
  delete ZZME;
}

void newMELA::setProcess(TVar::Process myModel, TVar::MatrixElement myME, TVar::Production myProduction)
{
  myModel_ = myModel;
  myME_ = myME;
  myProduction_ = myProduction;

  // 
  // configure the analytical calculations 
  // 

  if(!spin0Model->configure(myModel_)) pdf = spin0Model->PDF;
  else{
    if(!spin1Model->configure(myModel_)) pdf = spin1Model->PDF;
    else{
      if(!spin2Model->configure(myModel_,myProduction_)) pdf = spin2Model->PDF;
      else{
	if(myME_ == TVar::ANALYTICAL)
	  cout << "newMELA::setProcess -> ERROR TVar::Process not found!!! " << myME_ << endl; 
      }	
    }
  }

}


// Re-order masses and angles as needed by likelihoodDiscriminant. 
// This follows a different convention than the usual Z1/Z2 definition!
void newMELA::checkZorder(float& z1mass, float& z2mass,
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

void newMELA::computeP(float mZZ, float mZ1, float mZ2, // input kinematics
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
	integral->getVal();
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
    
    // 
    // define the constants to be used on JHUGen
    // 
    
    // gg productions 
    if ( myME_ == TVar::JHUGen && myProduction_ == TVar::GG  ) {
      if ( flavor == 3 ) {
	if ( myModel_ == TVar::PSHZZ_4l )  constant = 6.0;
	if ( myModel_ == TVar::HDHZZ_4l )  constant = 2.1;
	if ( myModel_ == TVar::VZZ_4l )  constant = 16;
	if ( myModel_ == TVar::AVZZ_4l )  constant = 13;
	if ( myModel_ == TVar::TZZ_4l )  constant = 0.6;
      }  else {
	if ( myModel_ == TVar::PSHZZ_4l )  constant = 7.0;
	if ( myModel_ == TVar::HDHZZ_4l )  constant = 2.3;
	if ( myModel_ == TVar::VZZ_4l )  constant = 38;
	if ( myModel_ == TVar::AVZZ_4l )  constant = 28;
	if ( myModel_ == TVar::TZZ_4l )  constant = 1.4;
      }
      if ( myModel_ == TVar::PTZZ_2hminus_4l )  constant = 1e+10;
      if ( myModel_ == TVar::TZZ_2hplus_4l )  constant = 1e+10;
    } 
    // qqb productions 
    if ( myME_ == TVar::JHUGen && myProduction_ == TVar::QQB  ) {
      if ( flavor == 3 ) {
	if ( myModel_ == TVar::TZZ_4l )  constant = 13;
      } else {
	if ( myModel_ == TVar::TZZ_4l )  constant = 30;
      }
    }
    // production independent calculations
    if ( myME_ == TVar::JHUGen && myProduction_ == TVar::INDEPENDENT  ) {
      if ( myModel_ == TVar::VZZ_4l )  constant = 1e+10;
      if ( myModel_ == TVar::AVZZ_4l )  constant = 1e+10;
      if ( myModel_ == TVar::TZZ_4l )  constant = 1.6e+9;
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
  prob = prob * constant; 
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

