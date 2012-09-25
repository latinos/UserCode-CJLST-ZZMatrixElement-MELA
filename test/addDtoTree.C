
using namespace RooFit;


void addDtoTree(char* inputFile,bool withPt=false, bool withY=false,int LHCsqrts=8){

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
        sigTree = (TTree*) sigFile->Get("SelectedTree");
    if(!sigTree){
      cout<<"ERROR could not find the tree!"<<endl;
      return;
    }

  TFile* newFile = new TFile(outputFileName,"RECREATE");
  TTree* newTree = new TTree("newTree","angles"); 

  float m1,m2,mzz,h1,h2,hs,phi,phi1,psig,pbkg,D,pseudoD,graviD;
  float pt4l, Y4l;
  float oldD;

  /*-------- UFL TREES --------
  double EL1,EL2,EL3,EL4;
  double pXL1,pXL2,pXL3,pXL4;
  double pYL1,pYL2,pYL3,pYL4;
  double pZL1,pZL2,pZL3,pZL4;
  double lomeD,oldD,mcfm_HZZ,mcfm_ZZ;

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

  sigTree->SetBranchAddress("lomeD",&lomeD);
  sigTree->SetBranchAddress("mcfm_HZZ",&mcfm_HZZ);
  sigTree->SetBranchAddress("mcfm_ZZ",&mcfm_ZZ);
  sigTree->SetBranchAddress("melaLD",&oldD);
  ------------------------------*/

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
  //---------------------------------------*/


  newTree->Branch("z1mass",&m1,"z1mass/F");
  newTree->Branch("z2mass",&m2,"z2mass/F");
  newTree->Branch("zzmass",&mzz,"zzmass/F");
  newTree->Branch("costheta1",&h1,"costheta1/F"); 
  newTree->Branch("costheta2",&h2,"costheta2/F");
  newTree->Branch("costhetastar",&hs,"costhetastar/F");
  newTree->Branch("phi",&phi,"phi/F");  
  newTree->Branch("phistar1",&phi1,"phistar1/F");

  newTree->Branch("ZZPt",&pt4l,"ZZpt/F");
  newTree->Branch("ZZRapidity",&Y4l,"ZZRapidity/F");
  
  newTree->Branch("ZZLD",&oldD,"ZZLD/F");
  //newTree->Branch("lomeD",&lomeD,"lomeD/D");
  //newTree->Branch("mcfm_HZZ",&mcfm_HZZ,"mcfm_HZZ/D");
  //newTree->Branch("mcfm_ZZ",&mcfm_ZZ,"mcfm_ZZ/D");

  newTree->Branch("Psig",&psig,"Psig/F");  
  newTree->Branch("Pbkg",&pbkg,"Pbkg/F");  

  newTree->Branch("melaLD",&D,"melaLD/F");  
  newTree->Branch("pseudoMelaLD",&pseudoD,"pseudoMelaLD/F");  
  newTree->Branch("spinTwoMinimalMelaLD",&graviD,"spinTwoMinimalMelaLD/F");  

  for(int iEvt=0; iEvt<sigTree->GetEntries(); iEvt++){

    if(iEvt%1000==0) 
      cout << "event: " << iEvt << endl;

    sigTree->GetEntry(iEvt);

    /*
    // ---------------- calculate angles ================
    // ---------------- from 4-vectors  
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
    */

    // --------------------------------

    if(mzz>350. && mzz<450.){

      // calculate discriminants from 4-vectors
      /*
      myPseudoMELA.eval(l1_minus, 11,
			l1_plus, -11,
			l2_minus, 13,
			l2_plus, -13,
			pseudoD,psig,pbkg);

      mySpinTwoMinimalMELA.eval(l1_minus, 11,
				l1_plus, -11,
				l2_minus, 13,
				l2_plus, -13,
				graviD,psig,pbkg);

      //MELA LD
      myMELA.computeKD(l1_minus, 11,
		       l1_plus, -11,
		       l2_minus, 13,
		       l2_plus, -13,
		       hs,h1,h2,phi,phi1,D,psig,pbkg,withPt,withY,LHCsqrts);

      // calculate discriminants from masses/angles
      */
      /*
      cout << "===================== " << endl;
      cout << "mzz: " << mzz << endl;
      cout << "m1: " << m1 << endl;
      cout << "m2: " << m2 << endl;
      cout << "h1: " << h1 << endl;
      cout << "h2: " << h2 << endl;
      cout << "phi: " << phi << endl;
      cout << "phi1: " << phi1 << endl;
      */

      //mySpinTwoMinimalMELA.eval(mzz,m1,m2,hs,h1,h2,phi,phi1,graviD,psig,pbkg);
      
      //myPseudoMELA.eval(mzz,m1,m2,hs,h1,h2,phi,phi1,pseudoD,psig,pbkg);

      myMELA.computeKD(mzz,m1,m2,hs,h1,h2,phi,phi1,D,psig,pbkg,withPt,pt4l,withY,Y4l,8);

      newTree->Fill();
      
    }

   }

  newFile->cd();
  newTree->Write("angles"); 
  newFile->Close();

}

