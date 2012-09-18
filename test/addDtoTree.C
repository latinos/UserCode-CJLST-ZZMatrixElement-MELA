
using namespace RooFit;


void addDtoTree(char* inputFile){


  gSystem->Load("$CMSSW_BASE/lib/slc5_amd64_gcc462/libZZMatrixElementMELA.so");
  gROOT->LoadMacro("../interface/Mela.h+");
  gROOT->LoadMacro("../interface/SpinTwoMinimalMELA.h+");
  gROOT->LoadMacro("../interface/PseudoMELA.h+");

  Mela myMELA;
  PseudoMELA myPseudoMELA;
  SpinTwoMinimalMELA mySpinTwoMinimalMELA;

  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);

  char inputFileName[100];
  char outputFileName[150];
  sprintf(inputFileName,"%s.root",inputFile);
  sprintf(outputFileName,"%s_withDiscriminants.root",inputFile);

  TFile* sigFile = new TFile(inputFileName);
  TTree* sigTree=0;
    if(sigFile)
        sigTree = (TTree*) sigFile->Get("passedEvents");
    if(!sigTree){
      cout<<"ERROR could not find the tree!"<<endl;
      return;
    }

  TFile* newFile = new TFile(outputFileName,"RECREATE");
  TTree* newTree = new TTree("newTree","angles"); 

  float m1,m2,mzz,h1,h2,hs,phi,phi1,psig,pbkg,D,pseudoD,graviD;

  double EL1,EL2,EL3,EL4;
  double pXL1,pXL2,pXL3,pXL4;
  double pYL1,pYL2,pYL3,pYL4;
  double pZL1,pZL2,pZL3,pZL4;

  sigTree->SetBranchAddress("EL1",&EL1);
  sigTree->SetBranchAddress("EL2",&EL2);
  sigTree->SetBranchAddress("EL3",&EL3);
  sigTree->SetBranchAddress("EL4",&EL4);

  sigTree->SetBranchAddress("pYL1",&pYL1);
  sigTree->SetBranchAddress("pYL2",&pYL2);
  sigTree->SetBranchAddress("pYL3",&pYL3);
  sigTree->SetBranchAddress("pYL4",&pYL4);

  sigTree->SetBranchAddress("pXL1",&pXL1);
  sigTree->SetBranchAddress("pXL2",&pXL2);
  sigTree->SetBranchAddress("pXL3",&pXL3);
  sigTree->SetBranchAddress("pXL4",&pXL4);

  sigTree->SetBranchAddress("pZL1",&pZL1);
  sigTree->SetBranchAddress("pZL2",&pZL2);
  sigTree->SetBranchAddress("pZL3",&pZL3);
  sigTree->SetBranchAddress("pZL4",&pZL4);

  newTree->Branch("z1mass",&m1,"z1mass/F");
  newTree->Branch("z2mass",&m2,"z2mass/F");
  newTree->Branch("zzmass",&mzz,"zzmass/F");
  newTree->Branch("costheta1",&h1,"costheta1/F"); 
  newTree->Branch("costheta2",&h2,"costheta2/F");
  newTree->Branch("costhetastar",&hs,"costhetastar/F");
  newTree->Branch("phi",&phi,"phi/F");  
  newTree->Branch("phistar1",&phi1,"phistar1/F");
  
  newTree->Branch("Psig",&psig,"Psig/F");  
  newTree->Branch("Pbkg",&pbkg,"Pbkg/F");  

  newTree->Branch("melaLD",&D,"melaLD/F");  
  newTree->Branch("pseudoMelaLD",&pseudoD,"pseudoMelaLD/F");  
  newTree->Branch("spinTwoMinimalMelaLD",&graviD,"spinTwoMinimalMelaLD/F");  

  for(int iEvt=0; iEvt<10000/*sigTree->GetEntries()*/; iEvt++){

    if(iEvt%1000==0) 
      cout << "event: " << iEvt << endl;

    sigTree->GetEntry(iEvt);

    // ---------------- calculate angles ================

    TLorentzVector l1_minus(pXL1,pYL1,pZL1,EL1);
    TLorentzVector l1_plus(pXL2,pYL2,pZL2,EL2);

    TLorentzVector p4Z1 = l1_minus+l1_plus;

    TLorentzVector l2_minus(pXL3,pYL3,pZL3,EL3);
    TLorentzVector l2_plus(pXL4,pYL4,pZL4,EL4);

    TLorentzVector p4Z2 = l2_minus+l2_plus;

    TLorentzVector p4H = p4Z1+p4Z2;

    if( p4Z1.M() > p4Z2.M() ){
      m1 = p4Z1.M();
      m2 = p4Z2.M();
    }else{
      m1 = p4Z2.M();
      m2 = p4Z1.M();
    }
    mzz = p4H.M();

    if(mzz>100. && mzz<1000.){

      //MELA LD
      myMELA.computeLD(l1_minus,l1_plus,l2_minus,l2_plus,hs,h1,h2,phi,phi1,D,psig,pbkg);
      myPseudoMELA.eval(l1_minus,l1_plus,l2_minus,l2_plus,pseudoD,psig,pbkg);
      mySpinTwoMinimalMELA.eval(l1_minus,l1_plus,l2_minus,l2_plus,graviD,psig,pbkg);

      newTree->Fill();
      
    }

   }

  newFile->cd();
  newTree->Write("angles"); 
  newFile->Close();

}

