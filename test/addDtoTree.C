
using namespace RooFit;


void addDtoTree(char* inputFile,int max=-1, int LHCsqrts=8){

  gSystem->Load("/afs/cern.ch/user/y/yygao/public/libmcfm.so");
  gSystem->Load("$CMSSW_BASE/lib/slc5_amd64_gcc462/libZZMatrixElementMELA.so");
  gROOT->LoadMacro("../interface/SpinTwoMinimalMELA.h+");
  gROOT->LoadMacro("../interface/SpinTwoDecayMELA.h+");
  gROOT->LoadMacro("../interface/SpinTwoqqMinimalMELA.h+");  

  SpinTwoMinimalMELA mySpinTwoMinimalMELA;
  SpinTwoDecayMELA mySpinTwoDecayMELA;
  SpinTwoqqMinimalMELA mySpinTwoqqMinimalMELA;

  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);

  char inputFileName[100];
  char outputFileName[150];
  sprintf(inputFileName,"%s.root",inputFile);
  sprintf(outputFileName,"%s_withDiscriminants.root",inputFile);

  TChain* sigTree = new TChain("SelectedTree");
  sigTree->Add(inputFileName);
  
  if(!sigTree){
    cout<<"ERROR could not find the tree!"<<endl;
    return;
  }

  TFile* newFile = new TFile(outputFileName,"RECREATE");
  TTree* newTree = (TTree*) sigTree->CloneTree(0,"fast");

  float m1,m2,mzz,h1,h2,hs,phi,phi1;
  float pt4l, Y4l;
  float psig, pbkg, spinTwo, spinTwoDecay, spinTwoQQ;

  // -------- CJLST TREES ---------------
  sigTree->SetBranchAddress("Z2Mass",&m2);
  sigTree->SetBranchAddress("Z1Mass",&m1);
  sigTree->SetBranchAddress("ZZMass",&mzz);
  sigTree->SetBranchAddress("costhetastar",&hs);
  sigTree->SetBranchAddress("helcosthetaZ1",&h1);
  sigTree->SetBranchAddress("helcosthetaZ2",&h2);
  sigTree->SetBranchAddress("helphi",&phi);
  sigTree->SetBranchAddress("phistarZ1",&phi1);

  //---------------------------------------*/

  newTree->Branch("spinTwoDecay",&spinTwoDecay,"spinTwoDecay/F");  
  newTree->Branch("spinTwo",&spinTwo,"spinTwo/F");  
  newTree->Branch("spinTwoQQ",&spinTwoQQ,"spinTwoQQ/F");  

  for(int iEvt=0; iEvt<(max<0?sigTree->GetEntries():max); iEvt++){

    if(iEvt>=sigTree->GetEntries()) break;

    if(iEvt%1000==0) 
      cout << "event: " << iEvt << endl;

    sigTree->GetEntry(iEvt);

    // --------------------------------

    /*
    cout << mzz << endl;
    cout << m1 << endl;
dD    cout << m2 << endl;
    cout << h1 << endl;
    cout << h2 << endl;
    cout << hs << endl;
    cout << phi << endl;
    cout << phi1 << endl;
    */

    if(mzz>100.){

      

      mySpinTwoMinimalMELA.computeKD((float)mzz,(float)m1,(float)m2,(float)hs,(float)h1,(float)h2,(float)phi,(float)phi1,spinTwo,psig,pbkg);

      //std::cout << "Gravi "<< graviD << " " << psig << " " << pbkg <<std::endl;

      mySpinTwoqqMinimalMELA.computeKD((float)mzz,(float)m1,(float)m2,(float)hs,(float)h1,(float)h2,(float)phi,(float)phi1,spinTwoQQ,psig,pbkg);
      
      //std::cout << "Pseudo "<< pseudoD << " " << psig << " " << pbkg <<std::endl;
      
      mySpinTwoDecayMELA.computeKD((float)mzz,(float)m1,(float)m2,(float)hs,(float)h1,(float)h2,(float)phi,(float)phi1,spinTwoDecay,psig,pbkg);

      //std::cout << "ICHEP "<< D_postICHEP << " " << psig << " " << pbkg <<std::endl;

      newTree->Fill();
      
    }

   }

  newFile->cd();
  newTree->Write("SelectedTree"); 
  newFile->Close();

}

