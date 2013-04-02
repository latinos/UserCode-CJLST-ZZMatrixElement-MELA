///////////////////////////////
// Add Probabilities to tree //
///////////////////////////////

using namespace RooFit;
bool includePathIsSet = false;

void addProbtoTreeNew(char* inputFile,int flavor, int max=-1, int LHCsqrts=8){
  //flavor: 1:4e, 2:4mu, 3:2e2mu


  gSystem->Load("$CMSSW_BASE/src/ZZMatrixElement/MELA/data/$SCRAM_ARCH/libmcfm.so");
  gSystem->Load("$CMSSW_BASE/lib/slc5_amd64_gcc462/libZZMatrixElementMELA.so");

  // set up path for local cmssw includes
  // as well as roofit
  if (!includePathIsSet) {
    TString path = gSystem->GetIncludePath();
    path += "-I$CMSSW_BASE/src/ ";
    path += "-I$ROOFITSYS/include/ ";
    gSystem->SetIncludePath(path.Data());

    // this is awkward, but want to protect against 
    // multiple additions to the base include path
    // if this function is invoked multiple times per session
    includePathIsSet = true;
  }

  // load and create instance of myMela class
  gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/MELA/interface/newMELA.h+");

  newMELA myMELA(LHCsqrts,flavor);
  
  cout << "newMELA initialized" << endl;

  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);

  char inputFileName[500];
  char outputFileName[500];
  sprintf(inputFileName,"%s.root",inputFile);
  sprintf(outputFileName,"%s_withProbabilities_ABtest.root",inputFile);

  TFile* sigFile = new TFile(inputFileName);
  TTree* sigTree=0;
    if(sigFile)
        sigTree = (TTree*) sigFile->Get("SelectedTree");
    if(!sigTree){
      //2nd try with the name of data obs tree
      sigTree = (TTree*) sigFile->Get("data_obs");
      if(!sigTree){
	cout<<"ERROR could not find the tree!"<<endl;
	return;
      }
    }


  float m1,m2,mzz,h1,h2,hs,phi,phi1;                                    //angles
  float psig,pbkg,D;
  double oldD;                                                    //legacy probabilities
  float	p0plus_melaNorm,p0plus_mela,p0minus_mela,p0hplus_mela,p0plus_VAJHU,p0minus_VAJHU,p0plus_VAMCFM,p0hplus_VAJHU,p1_mela,p1plus_mela,p1_VAJHU,p1plus_VAJHU,p2_mela,p2qqb_mela,p2_VAJHU,p2qqb_VAJHU; // new signal probablities
  float bkg_mela, bkg_VAMCFM, ggzz_VAMCFM, bkg_VAMCFMNorm, bkg_decay_VAMCFM;    // background probabilities

  // probabilities for production independent calculations
  float p1_decay_VAJHU, p1plus_decay_VAJHU, p2_decay_VAJHU; 
  // probabilities for exotic spin-2 models
  float p2hminus_VAJHU, p2hplus_VAJHU, p2bplus_VAJHU;
  float p2hminus_mela, p2hplus_mela, p2bplus_mela,p2_decay_mela;

  // -------- CJLST TREES ---------------
  sigTree->SetBranchAddress("Z2Mass",&m2);
  sigTree->SetBranchAddress("Z1Mass",&m1);
  sigTree->SetBranchAddress("ZZMass",&mzz);
  sigTree->SetBranchAddress("costhetastar",&hs);
  sigTree->SetBranchAddress("helcosthetaZ1",&h1);
  sigTree->SetBranchAddress("helcosthetaZ2",&h2);
  sigTree->SetBranchAddress("helphi",&phi);
  sigTree->SetBranchAddress("phistarZ1",&phi1);
  sigTree->SetBranchAddress("ZZLD",&oldD);
  Y4l=0.0;
  pt4l=0.0;

  float weight;
  sigTree->SetBranchAddress("MC_weight_noxsec",&weight);
  //---------------------------------------*/

  TFile* newFile = new TFile(outputFileName,"RECREATE");
  TTree* newTree = sigTree->CloneTree(0);//new TTree("newTree","SelectedTree"); 

  newTree->Branch("p0plus_mela_NEW",&p0plus_mela,"p0plus_mela_NEW/F");  // higgs, vector algebra, JHUgen
  newTree->Branch("p0minus_mela_NEW",&p0minus_mela,"p0minus_mela_NEW/F");  // pseudoscalar, vector algebra, JHUgen
  newTree->Branch("p0hplus_mela_NEW",&p0hplus_mela,"p0hplus_mela_NEW/F");  // 0h+, analytic distribution 
  newTree->Branch("p0plus_VAJHU_NEW",&p0plus_VAJHU,"p0plus_VAJHU_NEW/F");  // higgs, vector algebra, JHUgen
  newTree->Branch("p0minus_VAJHU_NEW",&p0minus_VAJHU,"p0minus_VAJHU_NEW/F");  // pseudoscalar, vector algebra, JHUgen
  newTree->Branch("p0plus_VAMCFM",&p0plus_VAMCFM,"p0plus_VAMCFM/F");  // higgs, vector algebra, MCFM
  newTree->Branch("p0hplus_VAJHU_NEW",&p0hplus_VAJHU,"p0hplus_VAJHU_NEW/F");  // 0h+(high order dimension operator) , vector algebra, JHUgen
  newTree->Branch("p1_mela_NEW",&p1_mela,"p1_mela_NEW/F");  // zprime, analytic distribution 
  newTree->Branch("p1plus_mela_NEW",&p1plus_mela,"p1plus_mela_NEW/F");  // 1+, analytic distribution 
  newTree->Branch("p1_VAJHU_NEW",&p1_VAJHU,"p1_VAJHU_NEW/F");  // zprime, vector algebra, JHUgen,
  newTree->Branch("p1plus_VAJHU_NEW",&p1plus_VAJHU,"p1plus_VAJHU_NEW/F");  // 1+ (axial-vector), vector algebra, JHUgen,
  newTree->Branch("p2_mela_NEW",&p2_mela ,"p2_mela_NEW/F");  // graviton, analytic distribution 
  newTree->Branch("p2qqb_mela_NEW",&p2qqb_mela,"p2qqb_mela_NEW/F"); // graviton produced by qqbar vector algebra, analytical,
  newTree->Branch("p2_VAJHU_NEW",&p2_VAJHU,"p2_VAJHU_NEW/F");  // graviton produced by gg, vector algebra, JHUgen,
  newTree->Branch("p2qqb_VAJHU_NEW",&p2qqb_VAJHU,"p2qqb_VAJHU_NEW/F");  // graviton produced by qqbar, vector algebra, JHUgen,
  newTree->Branch("p1_decay_VAJHU_NEW",&p1_decay_VAJHU,"p1_decay_VAJHU_NEW/F");  // 1-, vector algebra, production indpendent JHUgen
  newTree->Branch("p1plus_decay_VAJHU_NEW",&p1plus_decay_VAJHU,"p1plus_decay_VAJHU_NEW/F");  // 1+, vector algebra, production indpendent JHUgen
  newTree->Branch("p2_decay_VAJHU_NEW",&p2_decay_VAJHU,"p2_decay_VAJHU_NEW/F");  // 2m+, vector algebra, production indpendent JHUgen
  newTree->Branch("p2_decay_mela_NEW",&p2_decay_mela,"p2_decay_mela_NEW/F");  // 2m+, analytical, production indpendent JHUgen
  newTree->Branch("p2hminus_VAJHU_NEW",&p2hminus_VAJHU,"p2hminus_VAJHU_NEW/F");  // 2h-, vector algebra, JHUgen,
  newTree->Branch("p2hplus_VAJHU_NEW",&p2hplus_VAJHU,"p2hplus_VAJHU_NEW/F");     // 2h+, vector algebra, JHUgen,
  newTree->Branch("p2hminus_mela_NEW",&p2hminus_mela,"p2hminus_mela_NEW/F");  // 2h-, vector algebra, JHUgen,
  newTree->Branch("p2hplus_mela_NEW",&p2hplus_mela,"p2hplus_mela_NEW/F");     // 2h+, vector algebra, analytical,
  newTree->Branch("p2bplus_VAJHU_NEW",&p2bplus_VAJHU,"p2bplus_VAJHU_NEW/F");     // 2b+, vector algebra, JHUgen,

  //backgrounds
  newTree->Branch("bkg_mela_NEW",&bkg_mela,"bkg_mela_NEW/F");  // background,  analytic distribution 
  newTree->Branch("bkg_VAMCFM",&bkg_VAMCFM,"bkg_VAMCFM/F");  // background, vector algebra, MCFM
  newTree->Branch("ggzz_VAMCFM",&ggzz_VAMCFM,"ggzz_VAMCFM/F");  // background, vector algebra, MCFM for ggzz
  newTree->Branch("bkg_VAMCFMNorm",&bkg_VAMCFMNorm,"bkg_VAMCFMNorm/F");  // background, vector algebra, MCFM Normalized
  newTree->Branch("bkg_decay_VAMCFM",&bkg_decay_VAMCFM,"bkg_decay_VAMCFM/F");  // background, vector algebra, MCFM integrating out the production angles
  
  //interference weight
  float interfWeight;
  newTree->Branch("interfWeight",&interfWeight,"interfWeight/F"); // weight to be used for interference reweighting

  cout << "branches set" << endl;
  
  for(int iEvt=0; iEvt<(max<0?sigTree->GetEntries():max); iEvt++){
    
    if(iEvt>=sigTree->GetEntries()) break;
    
    // if ( iEvt != 2 ) continue;
    if(iEvt%1000==1) {
      cout<<"---------\n event: "<<iEvt<<endl;
    }
    sigTree->GetEntry(iEvt);

    // 
    // ANALYTICAL calcualtions
    // 

    // 0+
    myMELA.setProcess(TVar::HZZ_4l, TVar::ANALYTICAL, TVar::GG);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, 
		    p0plus_mela);
    

    // 0-
    myMELA.setProcess(TVar::PSHZZ_4l, TVar::ANALYTICAL, TVar::GG);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, 
		    p0minus_mela);
    
    
    // 0h+
    myMELA.setProcess(TVar::HDHZZ_4l, TVar::ANALYTICAL, TVar::GG);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, 
		    p0hplus_mela);
    
    // 1-
    myMELA.setProcess(TVar::VZZ_4l, TVar::ANALYTICAL, TVar::QQB);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, 
		    p1_mela);
    
    // production independent
    myMELA.setProcess(TVar::VZZ_4l, TVar::ANALYTICAL, TVar::INDEPENDENT);
    //myMELA.computeP(mzz, m1, m2, 
    //hs,h1,h2,phi,phi1,flavor, 
    //p1_decay_mela);
    
    // 1+
    myMELA.setProcess(TVar::AVZZ_4l, TVar::ANALYTICAL, TVar::QQB);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, 
		    p1plus_mela);
    
    // production independent
    //myMELA.setProcess(TVar::AVZZ_4l, TVar::ANALYTICAL, TVar::INDEPENDENT);
    //myMELA.computeP(mzz, m1, m2, 
    //hs,h1,h2,phi,phi1,flavor, 
    //p1plus_decay_mela);
    
    // gg->2m+
    myMELA.setProcess(TVar::TZZ_4l, TVar::ANALYTICAL, TVar::GG);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, 
		    p2_mela);
    
    // qqb->2m+
    myMELA.setProcess(TVar::TZZ_4l, TVar::ANALYTICAL, TVar::QQB);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, 
		    p2qqb_mela);
    
    // qqb->2m+
    myMELA.setProcess(TVar::TZZ_4l, TVar::ANALYTICAL, TVar::INDEPENDENT);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, 
		    p2_decay_mela);
    
    // 2h-
    myMELA.setProcess(TVar::PTZZ_2hminus_4l, TVar::ANALYTICAL, TVar::GG);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, 
		    p2hminus_mela);
    
    // 2h+
    myMELA.setProcess(TVar::TZZ_2hplus_4l, TVar::ANALYTICAL, TVar::GG);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, 
		    p2hplus_mela);
    
    // 2b+
    myMELA.setProcess(TVar::TZZ_2bplus_4l, TVar::ANALYTICAL, TVar::GG);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, 
		    p2bplus_mela);

    
    // 
    // MCFM calcualtionsA
    // 

    // qqZZ 

    if ( flavor == 3 ) {
      myMELA.setProcess(TVar::ZZ_2e2m, TVar::MCFM, TVar::QQB);
      myMELA.computeP(mzz, m1, m2, 
		      hs,h1,h2,phi,phi1,flavor, bkg_VAMCFM);
    } else {
      myMELA.setProcess(TVar::ZZ_4e, TVar::MCFM, TVar::QQB);
      myMELA.computeP(mzz, m1, m2, 
		      hs,h1,h2,phi,phi1,flavor, bkg_VAMCFM);
    }

    // ggZZ
    myMELA.setProcess(TVar::GGZZ_4l, TVar::MCFM, TVar::GG);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, ggzz_VAMCFM);

    // ******  qqZZ-> production independent 
    
    if ( flavor == 3 ) {
      myMELA.setProcess(TVar::ZZ_2e2m, TVar::MCFM, TVar::INDEPENDENT);
      myMELA.computeP(mzz, m1, m2, 
		      hs,h1,h2,phi,phi1,flavor, bkg_decay_VAMCFM);
    } else {
      myMELA.setProcess(TVar::ZZ_4e, TVar::MCFM, TVar::INDEPENDENT);
      myMELA.computeP(mzz, m1, m2, 
		      hs,h1,h2,phi,phi1,flavor, bkg_decay_VAMCFM);
    }

    // 0+
    myMELA.setProcess(TVar::HZZ_4l, TVar::MCFM, TVar::GG);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, p0plus_VAMCFM);


    // 
    // JHU Gen based signal calculations 
    // 
    
    // 0+
    myMELA.setProcess(TVar::HZZ_4l, TVar::JHUGen, TVar::GG);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, p0plus_VAJHU);
    

    // 0-
    myMELA.setProcess(TVar::PSHZZ_4l, TVar::JHUGen, TVar::GG);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, p0minus_VAJHU);
    
    
    // 0h+
    myMELA.setProcess(TVar::HDHZZ_4l, TVar::JHUGen, TVar::GG);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, p0hplus_VAJHU);
    
    // 1-
    myMELA.setProcess(TVar::VZZ_4l, TVar::JHUGen, TVar::QQB);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, p1_VAJHU);
    
    // production independent
    myMELA.setProcess(TVar::VZZ_4l, TVar::JHUGen, TVar::INDEPENDENT);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, p1_decay_VAJHU);
    
    // 1+
    myMELA.setProcess(TVar::AVZZ_4l, TVar::JHUGen, TVar::QQB);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, p1plus_VAJHU);
    
    // production independent
    myMELA.setProcess(TVar::AVZZ_4l, TVar::JHUGen, TVar::INDEPENDENT);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, p1plus_decay_VAJHU);
    
    // gg->2m+
    myMELA.setProcess(TVar::TZZ_4l, TVar::JHUGen, TVar::GG);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, p2_VAJHU);
    
    // qqb->2m+
    myMELA.setProcess(TVar::TZZ_4l, TVar::JHUGen, TVar::QQB);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, p2qqb_VAJHU);
    
    // qqb->2m+
    myMELA.setProcess(TVar::TZZ_4l, TVar::JHUGen, TVar::INDEPENDENT);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, p2_decay_VAJHU);
    
    // 2h-
    myMELA.setProcess(TVar::PTZZ_2hminus_4l, TVar::JHUGen, TVar::GG);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, p2hminus_VAJHU);
    
    // 2h+
    myMELA.setProcess(TVar::TZZ_2hplus_4l, TVar::JHUGen, TVar::GG);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, p2hplus_VAJHU);
    
    // 2b+
    myMELA.setProcess(TVar::TZZ_2bplus_4l, TVar::JHUGen, TVar::GG);
    myMELA.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,flavor, p2bplus_VAJHU);

    newTree->Fill();
    
  }
  
  newFile->cd();
  newTree->Write("SelectedTree"); 
  newFile->Close();
  
}

