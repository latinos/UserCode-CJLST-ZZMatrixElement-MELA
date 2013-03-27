///////////////////////////////
// Add Probabilities to tree //
///////////////////////////////
// mela analytic:  p0plus_melaNorm /(p0plus_melaNorm +  bkg_mela)
// mela template:  p0plus_mela_postICHEP /(p0plus_mela_postICHEP +  bkg_mela_postICHEP)
// gravimela:      p0plus_mela/(p0plus_mela +  p2_mela)
// pseudomela:     p0plus_mela /(p0plus_mela +  p0minus_mela)
// mela with pt/Y: p0_pt*p0_y*p0plus_melaNorm /(p0_pt*p0_y*p0plus_melaNorm +  bkg_pt*bkg_y*bkg_mela)
// ME:             p0plus_VAJHU/(p0plus_VAJHU + 10.*bkg_VAMCFM) 
// pME:            p0plus_VAJHU/(p0plus_VAJHU + 6.*p0minus_VAJHU)
// graviME:        p0plus_VAJHU/(p0plus_VAJHU + 1.2*p2_VAJHU) 
// VAKD            p0plus_VAJHU/( bkg_VAMCFMNorm +  p0plus_VAJHU ); // =normalized ME
//
// 
//

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

  newMELA myMELA(false,LHCsqrts,flavor);
  
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
  float bkg_mela, bkg_VAMCFM, ggzz_VAMCFM, bkg_VAMCFMNorm;                                           // new background probabilities
  float pt4l, Y4l ,p0_pt,p0_y,p0_y,bkg_pt,bkg_y;                        // rapidity/pt
  float p0plus_m4l,bkg_m4l,smd;  //supermela
  float p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown,p0plus_m4l_ResUp,p0plus_m4l_ResDown;//alternative values of supermela used for systematic templates
  float bkg_m4l_ScaleUp, bkg_m4l_ScaleDown,bkg_m4l_ResUp,bkg_m4l_ResDown;//alternative values of supermela used for systematic templates
  float p1_decay_VAJHU, p1plus_decay_VAJHU, p2_decay_VAJHU; 
  float p2hminus_VAJHU, p2hplus_VAJHU, p2bplus_VAJHU;

  float VAKD=0;                         // MCFM/JHUGen kinimetatic discriminant

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

  newTree->Branch("p0_mela_postICHEP",&psig,"p0_mela_postICHEP/F");
  newTree->Branch("pbkg_mela_postICHEP",&pbkg,"p0_mela_postICHEP/F");

  newTree->Branch("p0plus_melaNorm",&p0plus_melaNorm,"p0plus_melaNorm/F");  // higgs, vector algebra, JHUgen
  newTree->Branch("p0plus_mela",&p0plus_mela,"p0plus_mela/F");  // higgs, vector algebra, JHUgen
  newTree->Branch("p0minus_mela",&p0minus_mela,"p0minus_mela/F");  // pseudoscalar, vector algebra, JHUgen
  newTree->Branch("p0hplus_mela",&p0hplus_mela,"p0hplus_mela/F");  // 0h+, analytic distribution 
  newTree->Branch("p0plus_VAJHU",&p0plus_VAJHU,"p0plus_VAJHU/F");  // higgs, vector algebra, JHUgen
  newTree->Branch("p0minus_VAJHU",&p0minus_VAJHU,"p0minus_VAJHU/F");  // pseudoscalar, vector algebra, JHUgen
  newTree->Branch("p0plus_VAMCFM",&p0plus_VAMCFM,"p0plus_VAMCFM/F");  // higgs, vector algebra, MCFM
  newTree->Branch("p0hplus_VAJHU",&p0hplus_VAJHU,"p0hplus_VAJHU/F");  // 0h+(high order dimension operator) , vector algebra, JHUgen
  newTree->Branch("p1_mela",&p1_mela,"p1_mela/F");  // zprime, analytic distribution 
  newTree->Branch("p1plus_mela",&p1plus_mela,"p1plus_mela/F");  // 1+, analytic distribution 
  newTree->Branch("p1_VAJHU",&p1_VAJHU,"p1_VAJHU/F");  // zprime, vector algebra, JHUgen,
  newTree->Branch("p1plus_VAJHU",&p1plus_VAJHU,"p1plus_VAJHU/F");  // 1+ (axial-vector), vector algebra, JHUgen,
  newTree->Branch("p2_mela",&p2_mela ,"p2_mela/F");  // graviton, analytic distribution 
  newTree->Branch("p2qqb_mela",&p2qqb_mela,"p2qqb_mela/F"); // graviton produced by qqbar vector algebra, analytical,
  newTree->Branch("p2_VAJHU",&p2_VAJHU,"p2_VAJHU/F");  // graviton produced by gg, vector algebra, JHUgen,
  newTree->Branch("p2qqb_VAJHU",&p2qqb_VAJHU,"p2qqb_VAJHU/F");  // graviton produced by qqbar, vector algebra, JHUgen,
  newTree->Branch("p1_decay_VAJHU",&p1_decay_VAJHU,"p1_decay_VAJHU/F");  // 1-, vector algebra, production indpendent JHUgen
  newTree->Branch("p1plus_decay_VAJHU",&p1plus_decay_VAJHU,"p1plus_decay_VAJHU/F");  // 1+, vector algebra, production indpendent JHUgen
  newTree->Branch("p2_decay_VAJHU",&p2_decay_VAJHU,"p2_decay_VAJHU/F");  // 2m+, vector algebra, production indpendent JHUgen
  newTree->Branch("p2hminus_VAJHU",&p2hminus_VAJHU,"p2hminus_VAJHU/F");  // 2h-, vector algebra, JHUgen,
  newTree->Branch("p2hplus_VAJHU",&p2hplus_VAJHU,"p2hplus_VAJHU/F");     // 2h+, vector algebra, JHUgen,
  newTree->Branch("p2bplus_VAJHU",&p2bplus_VAJHU,"p2bplus_VAJHU/F");     // 2b+, vector algebra, JHUgen,

  //backgrounds
  newTree->Branch("bkg_mela",&bkg_mela,"bkg_mela/F");  // background,  analytic distribution 
  newTree->Branch("bkg_VAMCFM",&bkg_VAMCFM,"bkg_VAMCFM/F");  // background, vector algebra, MCFM
  newTree->Branch("ggzz_VAMCFM",&ggzz_VAMCFM,"ggzz_VAMCFM/F");  // background, vector algebra, MCFM for ggzz
  newTree->Branch("bkg_VAMCFMNorm",&bkg_VAMCFMNorm,"bkg_VAMCFMNorm/F");  // background, vector algebra, MCFM Normalized
  //supermela
  newTree->Branch("p0plus_m4l",&p0plus_m4l,"p0plus_m4l/F"  );  
  newTree->Branch("bkg_m4l",   &bkg_m4l, "bkg_m4l/F");  
  newTree->Branch("superMELA",&smd,"superMELA/F"  ); 
   newTree->Branch("p0plus_m4l_ScaleUp",&p0plus_m4l_ScaleUp,"p0plus_m4l_ScaleUp/F"  );  
   newTree->Branch("p0plus_m4l_ScaleDown",&p0plus_m4l_ScaleDown,"p0plus_m4l_ScaleDown/F"  );  
   newTree->Branch("p0plus_m4l_ResUp",&p0plus_m4l_ResUp,"p0plus_m4l_ResUp/F"  );  
   newTree->Branch("p0plus_m4l_ResDown",&p0plus_m4l_ResDown,"p0plus_m4l_ResDown/F"  );  
   newTree->Branch("bkg_m4l_ScaleUp",&bkg_m4l_ScaleUp,"bkg_m4l_ScaleUp/F"  );  
   newTree->Branch("bkg_m4l_ScaleDown",&bkg_m4l_ScaleDown,"bkg_m4l_ScaleDown/F"  );  
   newTree->Branch("bkg_m4l_ResUp",&bkg_m4l_ResUp,"bkg_m4l_ResUp/F"  );  
   newTree->Branch("bkg_m4l_ResDown",&bkg_m4l_ResDown,"bkg_m4l_ResDown/F"  );   
  //pt/rapidity
  newTree->Branch("p0_pt",&p0_pt,"p0_pt/F");  // multiplicative probability for signal pt
  newTree->Branch("p0_y",&p0_y,"p0_y/F");  // multiplicative probability for signal y
  newTree->Branch("bkg_pt",&bkg_pt,"bkg_pt/F");  // multiplicative probability for bkg pt
  newTree->Branch("bkg_y",&bkg_y,"bkg_y/F");  // multiplicative probability for bkg y
  newTree->Branch("VAKD",&VAKD,"VAKD/F");  // discriminant

  //interference weight
  float interfWeight;
  newTree->Branch("interfWeight",&interfWeight,"interfWeight/F"); // weight to be used for interference reweighting
  
  for(int iEvt=0; iEvt<(max<0?sigTree->GetEntries():max); iEvt++){
    
    if(iEvt>=sigTree->GetEntries()) break;
    
    if(iEvt%1000==1) {
      cout<<"---------\n event: "<<iEvt<<endl;
    }
    sigTree->GetEntry(iEvt);
    
    // --------------------------------

    if(mzz>100.){

      // JHU Gen based calculations 
      myMELA.setProcess(TVar::HZZ_4l, TVar::JHUGen, TVar::GG);
      myMELA.computeP(mzz, m1, m2, 
		      hs,h1,h2,phi,phi1,flavor, p0plus_VAJHU);
      
      newTree->Fill();
      
    }
    
  }
  
  newFile->cd();
  newTree->Write("SelectedTree"); 
  newFile->Close();

}

