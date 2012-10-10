
using namespace RooFit;


void addDtoTree(char* inputFile,int max=-1, int LHCsqrts=8){

  gSystem->Load("$CMSSW_BASE/lib/slc5_amd64_gcc462/libZZMatrixElementMELA.so");
  gROOT->LoadMacro("../interface/Mela.h+");
  gROOT->LoadMacro("../interface/SpinTwoMinimalMELA.h+");
  gROOT->LoadMacro("../interface/PseudoMELA.h+");

  Mela myMELA;
  Mela myMELA_ICHEP(true);
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
        sigTree = (TTree*) sigFile->Get("SelectedTree");
    if(!sigTree){
      cout<<"ERROR could not find the tree!"<<endl;
      return;
    }

  TFile* newFile = new TFile(outputFileName,"RECREATE");
  TTree* newTree = new TTree("newTree","SelectedTree"); 

  float m1,m2,mzz,h1,h2,hs,phi,phi1,psig,pbkg,D,pseudoD,graviD;
  float pt4l, Y4l;
  float oldD,D_PtY,D_postICHEP;

  // -------- CJLST TREES ---------------
  sigTree->SetBranchAddress("Z2Mass",&m2);
  sigTree->SetBranchAddress("Z1Mass",&m1);
  sigTree->SetBranchAddress("ZZMass",&mzz);
  sigTree->SetBranchAddress("costhetastar",&hs);
  sigTree->SetBranchAddress("helcosthetaZ1",&h1);
  sigTree->SetBranchAddress("helcosthetaZ2",&h2);
  sigTree->SetBranchAddress("helphi",&phi);
  sigTree->SetBranchAddress("phistarZ1",&phi1);
  sigTree->SetBranchAddress("ZZPt",&pt4l);
  sigTree->SetBranchAddress("ZZRapidity",&Y4l);
  sigTree->SetBranchAddress("ZZLD",&oldD);

  float weight;
  sigTree->SetBranchAddress("MC_weight_noxsec",&weight);
  //---------------------------------------*/


  newTree->Branch("Z1Mass",&m1,"Z1Mass/F");
  newTree->Branch("Z2Mass",&m2,"Z2Mass/F");
  newTree->Branch("ZZMass",&mzz,"ZZMass/F");
  newTree->Branch("helcosthetaZ1",&h1,"helcosthetaZ1/F"); 
  newTree->Branch("helcosthetaZ2",&h2,"helcosthetaZ2/F");
  newTree->Branch("costhetastar",&hs,"costhetastar/F");
  newTree->Branch("helphi",&phi,"helphi/F");  
  newTree->Branch("phistarZ1",&phi1,"phistarZ1/F");

  newTree->Branch("ZZPt",&pt4l,"ZZpt/F");
  newTree->Branch("ZZRapidity",&Y4l,"ZZRapidity/F");
  
  newTree->Branch("MC_weight_noxsec",&weight,"MC_weight_noxsec/F");

  newTree->Branch("ZZLD",&oldD,"ZZLD/F");

  newTree->Branch("Psig",&psig,"Psig/F");  
  newTree->Branch("Pbkg",&pbkg,"Pbkg/F");  

  newTree->Branch("ZZLD_analBkg",&D,"ZZLD_analBkg/F");  
  newTree->Branch("ZZLD_postICHEP",&D_postICHEP,"ZZLD_postICHEP/F");  
  newTree->Branch("ZZLD_PtY",&D_PtY,"ZZLD_PtY/F");  
  newTree->Branch("pseudoMelaLD",&pseudoD,"pseudoMelaLD/F");  
  newTree->Branch("spinTwoMinimalMelaLD",&graviD,"spinTwoMinimalMelaLD/F");  

  for(int iEvt=0; iEvt<(max<0?sigTree->GetEntries():max); iEvt++){

    if(iEvt>=sigTree->GetEntries()) break;

    if(iEvt%1000==0) 
      cout << "event: " << iEvt << endl;

    sigTree->GetEntry(iEvt);

    // --------------------------------

    if(mzz>100.){

      mySpinTwoMinimalMELA.computeKD(mzz,m1,m2,hs,h1,h2,phi,phi1,graviD,psig,pbkg);
      
      myPseudoMELA.computeKD(mzz,m1,m2,hs,h1,h2,phi,phi1,pseudoD,psig,pbkg);
      
      myMELA_ICHEP.computeKD(mzz,m1,m2,hs,h1,h2,phi,phi1,D_postICHEP,psig,pbkg,false,pt4l,false,Y4l,8);

      myMELA.computeKD(mzz,m1,m2,hs,h1,h2,phi,phi1,D_PtY,psig,pbkg,true,pt4l,true,Y4l,8);

      myMELA.computeKD(mzz,m1,m2,hs,h1,h2,phi,phi1,D,psig,pbkg,false,pt4l,false,Y4l,8);

      newTree->Fill();
      
    }

   }

  newFile->cd();
  newTree->Write("SelectedTree"); 
  newFile->Close();

}

