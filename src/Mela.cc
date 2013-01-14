/** \file
 *
 *  See header file for documentation.
 *
 *  $Date: 2013/01/14 15:30:44 $
 *  $Revision: 1.32 $
 */

#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZMatrixElement/MELA/interface/ZZMatrixElement.h>
#include <ZZMatrixElement/MELA/interface/SuperMELA.h>
#include <DataFormats/GeometryVector/interface/Pi.h>
#include <FWCore/ParameterSet/interface/FileInPath.h>

#include "computeAngles.h"
#include "AngularPdfFactory.h"
#include "TensorPdfFactory.h"
#include "RooqqZZ_JHU_ZgammaZZ_fast.h"
#include "RooqqZZ_JHU.h"
#include "RooTsallis.h"
//#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"  // replacement for RooTsallis
#include "RooTsallisExp.h"
#include "RooRapidityBkg.h"
#include "RooRapiditySig.h"

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

Mela::Mela(bool usePowhegTemplate, int LHCsqrts, float mh){ 

  usePowhegTemplate_=usePowhegTemplate;

  edm::FileInPath fip("ZZMatrixElement/MELA/data/my8DTemplateNotNorm.root");
  string fullPath = fip.fullPath();

  // Original code from KDProducer::beginJob
  //TFile* f =TFile::Open(fullPath.c_str(),"READ");
  f= new TFile(fullPath.c_str(),"READ");
  h_mzz= (TH1F*)(f->Get("h_mzz"));
  h_mzzm1m2= (TH3F*)(f->Get("h_mzzm1m2"));
  h_mzzcosthetastar= (TH2F*)(f->Get("h_mzzcosthetastar"));
  h_mzzcostheta1= (TH2F*)(f->Get("h_mzzcostheta1"));
  h_mzzcostheta2= (TH2F*)(f->Get("h_mzzcostheta2"));
  h_mzzphi1= (TH2F*)(f->Get("h_mzzphi1")); // This is phistar1
  h_mzzphi= (TH2F*)(f->Get("h_mzzphi"));
  // f->Close();

  z1mass_rrv = new RooRealVar("z1mass","m_{Z1}",0.,180.);
  z2mass_rrv = new RooRealVar("z2mass","m_{Z2}",0.,120.); 
  costhetastar_rrv = new RooRealVar("costhetastar","cos#theta^{*}",-1.,1.);  
  costheta1_rrv = new RooRealVar("costheta1","cos#theta_{1}",-1.,1.);  
  costheta2_rrv = new RooRealVar("costheta2","cos#theta_{2}",-1.,1.);
  phi_rrv= new RooRealVar("phi","#Phi",-3.1415,3.1415);
  phi1_rrv= new RooRealVar("phi1","#Phi_{1}",-3.1415,3.1415);
  pt_rrv= new RooRealVar("pt","p_{T}^{4l}",0.,1000);
  y_rrv= new RooRealVar("y","Y^{4l}",-4.0,4.0);
  sqrtS_rrv= new RooRealVar("sqrtS","#sqrt{s}",1000,14000);
  mzz_rrv= new RooRealVar("mzz","mZZ",80.,1000.);
  upFrac_rrv= new RooRealVar("upFrac","upFrac",.5);
  upFrac_rrv->setConstant(kTRUE);

  SMHiggs = new AngularPdfFactory(z1mass_rrv,z2mass_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,mzz_rrv);
  SMZgammaZZ = new RooqqZZ_JHU_ZgammaZZ_fast("SMZgammaZZ","SMZgammaZZ",*z1mass_rrv,*z2mass_rrv,*costheta1_rrv,*costheta2_rrv,*phi_rrv,*costhetastar_rrv,*phi1_rrv,*mzz_rrv,*upFrac_rrv);
  SMZZ = new RooqqZZ_JHU("SMZZ","SMZZ",*z1mass_rrv,*z2mass_rrv,*costheta1_rrv,*costheta2_rrv,*phi_rrv,*costhetastar_rrv,*phi1_rrv,*mzz_rrv);
  
  SMHiggs->makeSMHiggs();
  SMHiggs->makeParamsConst(true);
  

  PSHiggs = new AngularPdfFactory(z1mass_rrv,z2mass_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,mzz_rrv);

  PSHiggs->makePSHiggs();
  PSHiggs->makeParamsConst(true);


  minGrav = new TensorPdfFactory(z1mass_rrv,z2mass_rrv,costhetastar_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,phi1_rrv,mzz_rrv);

  minGrav->makeMinGrav();
  minGrav->makeParamsConst(true);

  // PDFs for Pt and Y

  static const int NptparamsS = 17;
  static const int NptparamsB = 11;

  string rrvnamesB[NptparamsB] = {"m","n0","n1","n2","ndue","bb0","bb1","bb2","T0","T1","T2"};
  allparamsB = new RooArgSet();
  for (int i = 0; i < NptparamsB; i++) {
    ptparamsB.push_back(new RooRealVar(rrvnamesB[i].c_str(),rrvnamesB[i].c_str(),-10000.,10000.));
    allparamsB->add(*ptparamsB[i]);
  }

  string rrvnamesS[NptparamsS] = {"ms","ns0","ns1","ns2","ndues","bbs0","bbs1","bbs2","Ts0","Ts1","Ts2","bbdues0","bbdues1","bbdues2","fexps0","fexps1","fexps2"};
  allparamsS = new RooArgSet();
  for (int i = 0; i < NptparamsS; i++) {
    ptparamsS.push_back(new RooRealVar(rrvnamesS[i].c_str(),rrvnamesS[i].c_str(),-10000.,10000.));
    allparamsS->add(*ptparamsS[i]);
  }
 
  char fileName[200];
  sprintf(fileName,"ZZMatrixElement/MELA/data/allParamsSig_%dTeV.txt",LHCsqrts);
  edm::FileInPath TsallisParams_Sig(fileName);
  fullPath = TsallisParams_Sig.fullPath();

  sigPt = new RooTsallisExp("sigPt","sigPt",*pt_rrv,*mzz_rrv,
					   *ptparamsS[0],*ptparamsS[1],*ptparamsS[2],
					   *ptparamsS[3],*ptparamsS[4],*ptparamsS[5],
					   *ptparamsS[6],*ptparamsS[7],*ptparamsS[8],
					   *ptparamsS[9],*ptparamsS[10],*ptparamsS[11],
					   *ptparamsS[12],*ptparamsS[13],*ptparamsS[14],
					   *ptparamsS[15],*ptparamsS[16]);

  allparamsS->readFromFile(fullPath.c_str(),0);

  sprintf(fileName,"ZZMatrixElement/MELA/data/allParamsBkg_%dTeV.txt",LHCsqrts);
  edm::FileInPath TsallisParam_Bkg(fileName);
  fullPath = TsallisParam_Bkg.fullPath();

  bkgPt = new RooTsallisSM("bkgPt","bkgPt",*pt_rrv,*mzz_rrv,
				     *ptparamsB[0],*ptparamsB[1],*ptparamsB[2],
				     *ptparamsB[3],*ptparamsB[4],*ptparamsB[5],
				     *ptparamsB[6],*ptparamsB[7],*ptparamsB[8],
				     *ptparamsB[9],*ptparamsB[10]);
  allparamsB->readFromFile(fullPath.c_str(),0);

  sigY = new RooRapiditySig("sigY", "sigY", *y_rrv, *mzz_rrv, *sqrtS_rrv);
  bkgY = new RooRapidityBkg("bkgY", "bkgY", *y_rrv, *mzz_rrv, *sqrtS_rrv);


  // ------

  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);

  edm::FileInPath HiggsWidthFile("Higgs/Higgs_CS_and_Width/txtFiles/8TeV-ggH.txt");
  std::string path = HiggsWidthFile.fullPath();
  //std::cout << path.substr(0,path.length()-12) << std::endl;
  
  //std::cout << "before supermela" << std::endl;
  super = new SuperMELA(mh,"4mu",LHCsqrts); // preliminary intialization, we adjust the flavor later
  char cardpath[500];
  sprintf(cardpath,"HZZ4L_Combination/CombinationPy/CreateDatacards/SM_inputs_%dTeV/inputs_4mu.txt",LHCsqrts);
  edm::FileInPath cardfile(cardpath);
  std::string cpath=cardfile.fullPath();
  std::cout << cpath.substr(0,cpath.length()-14).c_str()  <<std::endl;
  super->SetPathToCards(cpath.substr(0,cpath.length()-14).c_str() );
  super->SetVerbosity(false);
  super->init();
  //std::cout << "after supermela" << std::endl;

  // Create symlinks to the required files, if these are not already present (do nothing otherwse)
  edm::FileInPath mcfmInput1("ZZMatrixElement/MELA/data/input.DAT");
  edm::FileInPath mcfmInput2("ZZMatrixElement/MELA/data/process.DAT");
  edm::FileInPath mcfmInput3("ZZMatrixElement/MELA/data/Pdfdata/cteq6l1.tbl");  
  symlink(mcfmInput1.fullPath().c_str(), "input.DAT");
  symlink(mcfmInput2.fullPath().c_str(), "process.DAT");
  mkdir("Pdfdata",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  symlink(mcfmInput3.fullPath().c_str(), "Pdfdata/cteq6l1.tbl");

  ZZME = new  ZZMatrixElement(path.substr(0,path.length()-12 ).c_str(),1000.*LHCsqrts/2.,TVar::INFO);
  //std::cout << "after matrixele" << std::endl;
 
  edm::FileInPath ScaleFactorFile("ZZMatrixElement/MELA/data/scaleFactors.root");
  TFile* sf = TFile::Open(ScaleFactorFile.fullPath().c_str(),"r");
  vaScale_4e    = (TGraph*)sf->Get("scaleFactors_4e");
  vaScale_4mu   = (TGraph*)sf->Get("scaleFactors_4mu");
  vaScale_2e2mu = (TGraph*)sf->Get("scaleFactors_2e2mu");
  sf->Close(); 
  //vaScale->Print("v");
  //std::cout << "after scalefacts" << std::endl;

}

Mela::~Mela(){ 
  //std::cout << "begin destructor" << std::endl;  
  delete vaScale_4e;
  delete vaScale_4mu;
  delete vaScale_2e2mu;
  
  delete h_mzz;
  delete h_mzzm1m2;
  delete h_mzzcosthetastar;
  delete h_mzzcostheta1;
  delete h_mzzcostheta2;
  delete h_mzzphi1;
  delete h_mzzphi;
  delete f;

  delete z1mass_rrv; 
  delete z2mass_rrv; 
  delete costhetastar_rrv;
  delete costheta1_rrv;
  delete costheta2_rrv;
  delete phi_rrv;
  delete phi1_rrv;
  delete y_rrv;
  delete pt_rrv;
  delete sqrtS_rrv;
  delete mzz_rrv;
  delete upFrac_rrv;

  //std::cout << "destroy submes"<< std::endl;

  delete SMHiggs;
  delete SMZgammaZZ;
  delete SMZZ;
  delete PSHiggs;
  delete minGrav;
  
  delete ZZME;

  for(unsigned int i=0; i<ptparamsS.size(); i++){
    delete ptparamsS[i];
  }
  for(unsigned int i=0; i<ptparamsB.size(); i++){
    delete ptparamsB[i];
  }
  delete allparamsB;
  delete allparamsS;
  
  delete sigY;
  delete bkgY;

  delete bkgPt;
  delete sigPt;

  //std::cout << "end destructor"<< std::endl;
}

double Mela::sigPdfNorm(double mzz){

  double p0 =     -32.4817;
  double p1 =     0.567046;
  double p2 =  -0.00266482;
  double p3 =  7.04049e-06;
  double p4 = -1.13359e-08;
  double p5 =  1.15548e-11;
  double p6 = -7.49017e-15;
  double p7 =  2.99182e-18;
  double p8 = -6.71171e-22;
  double p9 =  6.46859e-26;

  double poly = p0+
    p1*mzz+
    p2*mzz*mzz+
    p3*mzz*mzz*mzz+
    p4*mzz*mzz*mzz*mzz+
    p5*mzz*mzz*mzz*mzz*mzz+
    p6*mzz*mzz*mzz*mzz*mzz*mzz+
    p7*mzz*mzz*mzz*mzz*mzz*mzz*mzz+
    p8*mzz*mzz*mzz*mzz*mzz*mzz*mzz*mzz+
    p9*mzz*mzz*mzz*mzz*mzz*mzz*mzz*mzz*mzz;

  return exp(poly);

}

double Mela::bkgPdfNorm(double mzz){

  mzz++;
    
  double p0;
  double p1;
  double p2;
  double p3;
  double p4;
  double p5;
  double p6;
  if(mzz<182){
 
    p0 =  1.36714e-09;
    p1 = -6.38091e-11;
    p2 =  1.22195e-12;
    p3 = -1.23122e-14;
    p4 =  6.89716e-17;
    p5 = -2.03987e-19;
    p6 =  2.49154e-22;

  }else if(mzz<210){

    p0 =  -4.83519e-08;
    p1 =   9.61952e-10;
    p2 =   -7.1795e-12;
    p3 =   2.38262e-14;
    p4 =  -2.96631e-17;
    p5 = 0.0; 
    p6 = 0.0;

  }else{

    p0 =  1.93511e-11;
    p1 = -8.65233e-14;
    p2 =  1.77038e-16;
    p3 = -1.95066e-19;
    p4 =  1.18991e-22;
    p5 = -3.78025e-26;
    p6 =  4.87841e-30;

  }
  return p0+
    p1*mzz+
    p2*mzz*mzz+
    p3*mzz*mzz*mzz+
    p4*mzz*mzz*mzz*mzz+
    p5*mzz*mzz*mzz*mzz*mzz+
    p6*mzz*mzz*mzz*mzz*mzz*mzz;

}

void Mela::computeWeight(float mZZ, float mZ1, float mZ2, 
			 float costhetastar,
			 float costheta1, 
			 float costheta2,
			 float phi,
			 float phistar1,
			 // return variables:
			 float& w
			 ){

  // vector algebra
  double dummy1,dummy2,dummy3,dummy4,dummy5,dummy6,dummy7,dummy8,dummy9,dummy10,dummy11;
  double dXsec_HZZ_JHU,dXsec_HZZ_JHU_interf; // temporary probabilities (FORTRAN functions will need double  not float)
  
  // calculate dXsec for 2e2mu
  ZZME->computeXS(mZZ,mZ1,mZ2,
		  costhetastar,costheta1,costheta2,
		  phi,phistar1,
		  3,
		  //output variables
		  dummy1, dummy2, dummy9,           // vars not used
		  dXsec_HZZ_JHU,
		  dummy3, dummy4, dummy5,    // vars not used
		  dummy6, dummy7, dummy8,dummy10,dummy11);   // vars not used

  // calculate dXsec for 4mu
  ZZME->computeXS(mZZ,mZ1,mZ2,
		  costhetastar,costheta1,costheta2,
		  phi,phistar1,
		  2,
		  //output variables
		  dummy1, dummy2, dummy9,      // vars not used
		  dXsec_HZZ_JHU_interf,
		  dummy3, dummy4, dummy5,    // vars not used
		  dummy6, dummy7, dummy8,dummy10,dummy11);   // vars not used
  
  w = dXsec_HZZ_JHU_interf / dXsec_HZZ_JHU;

  // protect against anomalously large weights
  if (w>10.) w=0.;

}


void Mela::computeKD(TLorentzVector Z1_lept1, int Z1_lept1Id,
		     TLorentzVector Z1_lept2, int Z1_lept2Id,
		     TLorentzVector Z2_lept1, int Z2_lept1Id,
		     TLorentzVector Z2_lept2, int Z2_lept2Id,
		     // return variables:
		     float& costhetastar,
		     float& costheta1, 
		     float& costheta2,
		     float& phi,
		     float& phistar1,
		     float& kd, 
		     float& psig,
		     float& pbkg,
		     bool withPt,
		     bool withY) {
  
  //compute angles  
  float m1=(Z1_lept1 + Z1_lept2).M();
  float m2=(Z2_lept1 + Z2_lept2).M();

  TLorentzVector ZZ = (Z1_lept1 + Z1_lept2 + Z2_lept1 + Z2_lept2);
  float mzz = ZZ.M();

  // Skip candidates where KD is irrelevant.
  if (mzz<100.){
    kd = 0;
    psig = 0;
    pbkg = 0;
    return;
  }

  float pt  = ZZ.Pt();
  float Y   = ZZ.Rapidity(); // Fixme: should probably protect against NaN?

  mela::computeAngles(Z1_lept1, Z1_lept1Id, Z1_lept2, Z1_lept2Id, 
		      Z2_lept1, Z2_lept1Id, Z2_lept2, Z2_lept2Id,
		      costhetastar,costheta1,costheta2,phi,phistar1);

  //compute kd
  pair<float,float> P = likelihoodDiscriminant(mzz,m1,m2,costhetastar,costheta1,costheta2,phi,phistar1,
					       withPt, pt, withY, Y);

  psig=P.first;
  pbkg=P.second;
  kd = psig/(psig+pbkg);
  
}



void Mela::computeKD(float mzz, float mZ1, float mZ2, 
		     float costhetastar,
		     float costheta1, 
		     float costheta2,
		     float phi,
		     float phistar1,
		     float& kd, 
		     float& psig,
		     float& pbkg,
		     bool withPt,
		     float pt4l,
		     bool withY,
		     float Y4l) {

  // Skip candidates where KD is irrelevant.
  if (mzz<100.){
    kd = 0;
    psig = 0;
    pbkg = 0;
    return;
  }

  pair<float,float> P = likelihoodDiscriminant(mzz, mZ1, mZ2, 
					       costhetastar, 
					       costheta1, 
					       costheta2, 
					       phi, 
					       phistar1,
					       withPt,pt4l,
					       withY, Y4l);


  psig = P.first;
  pbkg = P.second;
  kd = psig/(psig + pbkg);

}



vector<float> Mela::my8DTemplate(bool normalized,float mZZ, float m1, float m2, float costhetastar, float costheta1, float costheta2, float phi, float phistar1){

  //multiply the P values
  float n = h_mzz->GetBinContent(h_mzz->FindBin(mZZ));
  float Pmzzm1m2 = h_mzzm1m2->GetBinContent(h_mzzm1m2->FindBin(mZZ,m1,m2));

  // - - - - - - - - - - - - - - - whitbeck
  // if bin has no events: add 1
  // safety feature to prevent KD = 1 as a
  // result of low statistics

  if(Pmzzm1m2==0){
    Pmzzm1m2++;
    }
  // - - - - - - - - - - - - - - - 

  float Pmzzcosthetastar = h_mzzcosthetastar->GetBinContent(h_mzzcosthetastar->FindBin(mZZ,costhetastar));
  float Pmzzcostheta2 = h_mzzcostheta2->GetBinContent(h_mzzcostheta2->FindBin(mZZ,costheta2));
  float Pmzzcostheta1 = h_mzzcostheta1->GetBinContent(h_mzzcostheta1->FindBin(mZZ,costheta1));
  float Pmzzphi1 = h_mzzphi1->GetBinContent(h_mzzphi1->FindBin(mZZ,phistar1));
  float Pmzzphi = h_mzzphi->GetBinContent(h_mzzphi->FindBin(mZZ,phi));

  //normalization
  float binwidth_mzzm1m2 = h_mzzm1m2->GetYaxis()->GetBinWidth(1) * h_mzzm1m2->GetZaxis()->GetBinWidth(1);
  float binwidth_mzzcosthetastar = h_mzzcosthetastar->GetYaxis()->GetBinWidth(1);
  float binwidth_mzzcostheta1 = h_mzzcostheta1->GetYaxis()->GetBinWidth(1);
  float binwidth_mzzcostheta2 = h_mzzcostheta1->GetYaxis()->GetBinWidth(1);
  float binwidth_mzzphi1 = h_mzzphi1->GetYaxis()->GetBinWidth(1);
  float binwidth_mzzphi = h_mzzphi->GetYaxis()->GetBinWidth(1);

  float Pmzzm1m2_norm = Pmzzm1m2/(n*binwidth_mzzm1m2); 
  float Pmzzcosthetastar_norm = Pmzzcosthetastar/(n*binwidth_mzzcosthetastar);
  float Pmzzcostheta1_norm = Pmzzcostheta1/(n*binwidth_mzzcostheta1);
  float Pmzzcostheta2_norm = Pmzzcostheta2/(n*binwidth_mzzcostheta2);
  float Pmzzphi1_norm = Pmzzphi1/(n*binwidth_mzzphi1);
  float Pmzzphi_norm = Pmzzphi/(n*binwidth_mzzphi);

  vector <float> P;
  P.push_back(Pmzzm1m2);
  P.push_back(Pmzzcosthetastar);
  P.push_back(Pmzzcostheta1);
  P.push_back(Pmzzcostheta2);
  P.push_back(Pmzzphi);
  P.push_back(Pmzzphi1);

  vector <float> P_norm;
  P_norm.push_back(Pmzzm1m2_norm);
  P_norm.push_back(Pmzzcosthetastar_norm);
  P_norm.push_back(Pmzzcostheta1_norm);
  P_norm.push_back(Pmzzcostheta2_norm);
  P_norm.push_back(Pmzzphi_norm);
  P_norm.push_back(Pmzzphi1_norm);

  /*delete h_mzz;
  delete h_mzzm1m2;
  delete h_mzzcosthetastar;
  delete h_mzzcostheta1;
  delete h_mzzcostheta2;
  delete h_mzzphi1;
  delete h_mzzphi;*/
  
  if(normalized)
    return P_norm;
  else
    return P;
}

pair<float,float> Mela::likelihoodDiscriminant (float mZZ, float m1, float m2, float costhetastar, float costheta1, float costheta2, float phi, float phistar1,
						bool withPt, float pt, 
						bool withY, float y){


  checkZorder(m1,m2,costhetastar,costheta1,costheta2,phi,phistar1);
  
  z1mass_rrv->setVal(m1);  
  z2mass_rrv->setVal(m2);
  costhetastar_rrv->setVal(costhetastar);
  costheta1_rrv->setVal(costheta1);
  costheta2_rrv->setVal(costheta2);
  phi_rrv->setVal(phi);
  phi1_rrv->setVal(phistar1);
  mzz_rrv->setVal(mZZ);
  if(withPt)
    pt_rrv->setVal(pt);
  if(withY)
    y_rrv->setVal(y);

  vector <float> P=my8DTemplate(1, mZZ,  m1,  m2,  costhetastar,  costheta1,  costheta2,  phi,  phistar1);
  
  float Pbackg=-99;
  float Psig=-99; 

  RooAbsReal *tempIntegral_SMHiggs=0, *tempIntegral_SMZZ=0;
  RooAbsReal *tempIntegral_bkgPt=0, *tempIntegral_sigPt=0;
  RooAbsReal *tempIntegral_bkgY=0, *tempIntegral_sigY=0;

  if(usePowhegTemplate_){
    // using template background calculation
    if(mZZ>100 && mZZ<180){
      Pbackg = P[0]*P[1]*P[2]*P[3]*P[4]*P[5]*5.0;
      Psig=SMHiggs->getVal(mZZ);
    }
    if(mZZ>180&&mZZ<=2*91.188){
      z1mass_rrv->setVal(mZZ/2.-1e-9);
      z2mass_rrv->setVal(mZZ/2.-1e-9);
      tempIntegral_SMZZ = SMZZ->createIntegral(RooArgSet(*costhetastar_rrv,*costheta1_rrv,*costheta2_rrv,*phi_rrv,*phi1_rrv));
      Pbackg = SMZZ->getVal()/tempIntegral_SMZZ->getVal()*10.0;
      tempIntegral_SMHiggs=SMHiggs->PDF->createIntegral(RooArgSet(*costheta1_rrv,*costheta2_rrv,*phi_rrv));
      Psig = SMHiggs->PDF->getVal()/tempIntegral_SMHiggs->getVal();
    }
    if(mZZ>2*91.188){
      z1mass_rrv->setVal(91.188);
      z2mass_rrv->setVal(91.188);
      tempIntegral_SMZZ=SMZZ->createIntegral(RooArgSet(*costhetastar_rrv,*costheta1_rrv,*costheta2_rrv,*phi_rrv,*phi1_rrv));
      Pbackg = SMZZ->getVal()/tempIntegral_SMZZ->getVal()*10.0;
      tempIntegral_SMHiggs=SMHiggs->PDF->createIntegral(RooArgSet(*costheta1_rrv,*costheta2_rrv,*phi_rrv));
      Psig = SMHiggs->PDF->getVal()/tempIntegral_SMHiggs->getVal();
    }
  }else{
    
    // using analytic background calculation
    Pbackg = SMZgammaZZ->getVal()*2e3/bkgPdfNorm(mZZ); 
    Psig = SMHiggs->PDF->getVal()/sigPdfNorm(mZZ);
    
  }

  if (withPt) {
    //cout << "withPt ";
    tempIntegral_bkgPt=bkgPt->createIntegral(RooArgSet(*pt_rrv));
    Pbackg *= bkgPt->getVal()/tempIntegral_bkgPt->getVal();
    //cout << Pbackg << " " ;
    tempIntegral_sigPt=sigPt->createIntegral(RooArgSet(*pt_rrv));
    Psig *= sigPt->getVal()/tempIntegral_sigPt->getVal();
    //cout << Psig << " " << endl;
  }
  if(withY) {
    //cout << "withY ";
    tempIntegral_bkgY=bkgY->createIntegral(RooArgSet(*y_rrv));
    Pbackg *= bkgY->getVal()/tempIntegral_bkgY->getVal();
    //cout << Pbackg << " " ;
    tempIntegral_sigY=sigY->createIntegral(RooArgSet(*y_rrv));
    Psig *= sigY->getVal()/tempIntegral_sigY->getVal();
    //cout << Psig << endl;
  }

  delete tempIntegral_SMZZ;
  delete tempIntegral_SMHiggs;
  delete tempIntegral_sigPt;
  delete tempIntegral_bkgPt;
  delete tempIntegral_sigY;
  delete tempIntegral_bkgY;

  // - - - - - - - - - - - - - - - - - - - - - Whitbeck 
  // check whether P[i] is zero and print warning
  // message if so

  string varName[6]={"m1/m2","costhetastar","costheta1","coshteta2","phi","phistar1"};
  for(int iVar=0; iVar<6; iVar++){

    if(P[iVar]==0 && (m1+m2)<mZZ && m2>4 && mZZ>80 && mZZ<180)
	cout << " uh oh... Probability of " << varName[iVar] << " is zero." << endl;
  }
  // - - - - - - - - - - - - - - - - - - - - - 

  if(Psig<0 || Pbackg<0){
    cout<<"Mela::likelihoodDiscriminant() Error: KD not defined for this mzz (maybe mZZ<100 ?)"<<endl;
    cout << "=========================" << endl;
    cout << "psig: " << Psig << endl;
    cout << "pbkg: " << Pbackg << endl;
    cout << " - - - - - - - - - - - - " << endl;
    cout << "mzz: " << mZZ << endl;
    cout << "m1: " << m1 << endl;
    cout << "m2: " << m2 << endl;
    cout << "costheta1: " << costheta1 << endl;
    cout << "costheta2: " << costheta2 << endl;
    cout << "costhetastar: " << costhetastar << endl;
    cout << "phi: " << phi << endl;
    cout << "phi1: " << phistar1 << endl;
  }

  return make_pair(Psig,Pbackg);

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
void Mela::computeP(float mZZ, float mZ1, float mZ2, 
		    float costhetastar,
		    float costheta1, 
		    float costheta2,
		    float phi,
		    float phistar1,
		    //signal probabilities
		    float& p0plus_melaNorm,   // higgs, analytic distribution          
		    float& p0plus_mela,   // higgs, analytic distribution          
		    float& p0minus_mela,  // pseudoscalar, analytic distribution 
		    float& p0plus_VAJHU,  // higgs, vector algebra, JHUgen
		    float& p0minus_VAJHU, // pseudoscalar, vector algebra, JHUgen
		    float& p0plus_VAMCFM,// higgs, vector algebra, MCFM
		    float& p1_mela,  // zprime, analytic distribution 
		    float& p1_VAJHU, // zprime, vector algebra, JHUgen,
		    float& p1plus_VAJHU, // 1+ (axial vector), vector algebra, JHUgen,
		    float& p2_mela , // graviton, analytic distribution 
		    float& p2_VAJHU, // graviton, vector algebra, JHUgen,
		    float& p2qqb_VAJHU, // graviton produced by qqbar vector algebra, JHUgen,
		    //backgrounds
		    float& bkg_mela,  // background,  analytic distribution 
		    float& bkg_VAMCFM, // background, vector algebra, MCFM
		    float& ggzz_VAMCFM, // background, vector algebra, MCFM for ggzz
		    float& bkg_VAMCFMNorm, // background, vector algebra, MCFM Normalized
		    //pt/rapidity
		    float& p0_pt, // multiplicative probability for signal pt
		    float& p0_y, // multiplicative probability for signal y
		    float& bkg_pt, // multiplicative probability for bkg pt
		    float& bkg_y, // multiplicative probability for bkg y
		    // supermela
		    float& p0plus_m4l,  // signal m4l probability as in datacards
		    float& bkg_m4l,     // backgroun m4l probability as in datacards
		    //optional input parameters
		    float pt4l,
		    float Y4l,
		    int flavor // 1:4e, 2:4mu, 3:2e2mu (for interference effects)
		    ){
  

  //initialize variables
  checkZorder(mZ1,mZ2,costhetastar,costheta1,costheta2,phi,phistar1);


  //std mela variables (will initialize RooRealVars for other PDfs too)
  pair<float,float> P = likelihoodDiscriminant(mZZ, mZ1, mZ2, 
					       costhetastar, 
					       costheta1, 
					       costheta2, 
					       phi, 
					       phistar1,
					       false,0.0,
					       false,0.0);

  p0plus_melaNorm = P.first;
  bkg_mela    = P.second;
  p0plus_mela = SMHiggs->getVal(mZZ);

  // pseudo mela
  p0minus_mela = PSHiggs->getVal(mZZ);

  // Z'
  p1_mela =0.0 ; // not implemented yet.
  
  //graviMela
  p2_mela = minGrav->getVal(mZZ);


  //compute pt/Y probabilities:
  pt_rrv->setVal(pt4l);
  y_rrv->setVal(Y4l);

  RooAbsReal *tempIntegral_bkgPt=0, *tempIntegral_sigPt=0;
  RooAbsReal *tempIntegral_bkgY=0, *tempIntegral_sigY=0;

  tempIntegral_bkgPt=bkgPt->createIntegral(RooArgSet(*pt_rrv));
  bkg_pt = bkgPt->getVal()/tempIntegral_bkgPt->getVal();
  delete tempIntegral_bkgPt;
  tempIntegral_sigPt=sigPt->createIntegral(RooArgSet(*pt_rrv));
  p0_pt = sigPt->getVal()/tempIntegral_sigPt->getVal();
  delete tempIntegral_sigPt;

  tempIntegral_bkgY=bkgY->createIntegral(RooArgSet(*y_rrv));
  bkg_y = bkgY->getVal()/tempIntegral_bkgY->getVal();
  delete tempIntegral_bkgY;
  tempIntegral_sigY=sigY->createIntegral(RooArgSet(*y_rrv));
  p0_y = sigY->getVal()/tempIntegral_sigY->getVal();
  delete tempIntegral_sigY;


  // vector algebra
  double dummy1,dummy2,dummy3;
  double dXsec_ZZ_MCFM,dXsec_GGZZ_MCFM,dXsec_HZZ_MCFM,dXsec_HZZ_JHU,dXsec_PSHZZ_JHU,dXsec_VZZ_JHU,dXsec_AVZZ_JHU,dXsec_TZZ_JHU,dXsec_QQB_TZZ_JHU; // temporary probabilities (FORTRAN functions will need double, not float)
  ZZME->computeXS(mZZ,mZ1,mZ2,
		  costhetastar,costheta1,costheta2,
		  phi,phistar1,
		  flavor,
		  //output variables
		  dXsec_ZZ_MCFM,
		  dXsec_GGZZ_MCFM,
		  dXsec_HZZ_MCFM,
		  dXsec_HZZ_JHU,
		  dXsec_PSHZZ_JHU,
		  dXsec_VZZ_JHU,
		  dXsec_TZZ_JHU,
		  dXsec_AVZZ_JHU,
		  dXsec_QQB_TZZ_JHU,
		  dummy1,dummy2,dummy3);   // discriminants are not used 

  bkg_VAMCFM=dXsec_ZZ_MCFM;
  ggzz_VAMCFM=dXsec_GGZZ_MCFM;
  p0plus_VAMCFM=dXsec_HZZ_MCFM;
  p0plus_VAJHU=dXsec_HZZ_JHU;
  p0minus_VAJHU=dXsec_PSHZZ_JHU;
  p1_VAJHU=dXsec_VZZ_JHU;
  p1plus_VAJHU=dXsec_AVZZ_JHU;
  p2_VAJHU=dXsec_TZZ_JHU;
  p2qqb_VAJHU=dXsec_QQB_TZZ_JHU;

  if(flavor==1){
    if(mZZ > 900)                   
      bkg_VAMCFMNorm = bkg_VAMCFM * vaScale_4e->Eval(900.);
    else if (mZZ <  100 )
      bkg_VAMCFMNorm = bkg_VAMCFM * vaScale_4e->Eval(100.);
    else
      bkg_VAMCFMNorm = bkg_VAMCFM * vaScale_4e->Eval(mZZ);
  }

  if(flavor==2){
    if(mZZ > 900)                   
      bkg_VAMCFMNorm = bkg_VAMCFM * vaScale_4mu->Eval(900.);
    else if (mZZ <  100 )
      bkg_VAMCFMNorm = bkg_VAMCFM * vaScale_4mu->Eval(100.);
    else
      bkg_VAMCFMNorm = bkg_VAMCFM * vaScale_4mu->Eval(mZZ);
  }

  if(flavor==3){
    if(mZZ > 900)                   
      bkg_VAMCFMNorm = bkg_VAMCFM * vaScale_2e2mu->Eval(900.);
    else if (mZZ <  100 )
      bkg_VAMCFMNorm = bkg_VAMCFM * vaScale_2e2mu->Eval(100.);
    else
      bkg_VAMCFMNorm = bkg_VAMCFM * vaScale_2e2mu->Eval(mZZ);
  }

  // supermela
  switch(flavor){
  case 1: super->SetDecayChannel("4e")   ;break;
  case 2: super->SetDecayChannel("4mu")  ;break;
  case 3: super->SetDecayChannel("2e2mu");break;
  default: std::cout << " unknown flavor: " << flavor << std::endl; exit(0);
  }
  
  std::pair<double,double> m4lP = super->M4lProb(mZZ);
  p0plus_m4l = m4lP.first; 
  bkg_m4l = m4lP.second;


}
