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


void addProbtoTree(char* inputFile,int flavor, int max=-1, int LHCsqrts=8){
  //flavor: 1:4e, 2:4mu, 3:2e2mu


  gSystem->Load("/afs/cern.ch/user/y/yygao/public/libmcfm.so");
  gSystem->Load("$CMSSW_BASE/lib/slc5_amd64_gcc462/libZZMatrixElementMELA.so");
  gROOT->LoadMacro("../interface/Mela.h+");

  Mela myMELA(false,LHCsqrts);
  Mela myMELA_ICHEP(true,LHCsqrts);

  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);

  char inputFileName[500];
  char outputFileName[500];
  sprintf(inputFileName,"%s.root",inputFile);
  sprintf(outputFileName,"%s_withProbabilities.root",inputFile);

  TFile* sigFile = new TFile(inputFileName);
  TTree* sigTree=0;
    if(sigFile)
        sigTree = (TTree*) sigFile->Get("SelectedTree");
    if(!sigTree){
      cout<<"ERROR could not find the tree!"<<endl;
      return;
    }


  float m1,m2,mzz,h1,h2,hs,phi,phi1;                                    //angles
  float psig,pbkg,D,oldD;                                                    //legacy probabilities
  float	p0plus_melaNorm,p0plus_mela,p0minus_mela,p0plus_VAJHU,p0minus_VAJHU,p0plus_VAMCFM,p1_mela,p1_VAJHU,p2_mela,p2_VAJHU; // new signal probablities
  float bkg_mela, bkg_VAMCFM,bkg_VAMCFMNorm;                                           // new background probabilities
  float pt4l, Y4l ,p0_pt,p0_y,p0_y,bkg_pt,bkg_y;                        // rapidity/pt
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
  sigTree->SetBranchAddress("ZZPt",&pt4l);
  sigTree->SetBranchAddress("ZZRapidity",&Y4l);
  sigTree->SetBranchAddress("ZZLD",&oldD);



  float weight;
  sigTree->SetBranchAddress("MC_weight_noxsec",&weight);
  //---------------------------------------*/

  TFile* newFile = new TFile(outputFileName,"RECREATE");
  TTree* newTree = sigTree->CloneTree(0);//new TTree("newTree","SelectedTree"); 

//   newTree->Branch("Z1Mass",&m1,"Z1Mass/F");
//   newTree->Branch("Z2Mass",&m2,"Z2Mass/F");
//   newTree->Branch("ZZMass",&mzz,"ZZMass/F");
//   newTree->Branch("helcosthetaZ1",&h1,"helcosthetaZ1/F"); 
//   newTree->Branch("helcosthetaZ2",&h2,"helcosthetaZ2/F");
//   newTree->Branch("costhetastar",&hs,"costhetastar/F");
//   newTree->Branch("helphi",&phi,"helphi/F");  
//   newTree->Branch("phistarZ1",&phi1,"phistarZ1/F");

//   newTree->Branch("ZZPt",&pt4l,"ZZpt/F");
//   newTree->Branch("ZZRapidity",&Y4l,"ZZRapidity/F");
  
//   newTree->Branch("MC_weight_noxsec",&weight,"MC_weight_noxsec/F");

  newTree->Branch("p0_mela_postICHEP",&psig,"p0_mela_postICHEP/F");
  newTree->Branch("pbkg_mela_postICHEP",&pbkg,"p0_mela_postICHEP/F");

  newTree->Branch("p0plus_melaNorm",&p0plus_mela,"p0plus_melaNorm/F");  // higgs, vector algebra, JHUgen
  newTree->Branch("p0plus_mela",&p0plus_mela,"p0plus_mela/F");  // higgs, vector algebra, JHUgen
  newTree->Branch("p0minus_mela",&p0minus_mela,"p0minus_mela/F");  // pseudoscalar, vector algebra, JHUgen
  newTree->Branch("p0plus_VAJHU",&p0plus_VAJHU,"p0plus_VAJHU/F");  // higgs, vector algebra, JHUgen
  newTree->Branch("p0minus_VAJHU",&p0minus_VAJHU,"p0minus_VAJHU/F");  // pseudoscalar, vector algebra, JHUgen
  newTree->Branch("p0plus_VAMCFM",&p0plus_VAMCFM,"p0plus_VAMCFM/F");  // higgs, vector algebra, MCFM
  newTree->Branch("p1_mela",&p1_mela,"p1_mela/F");  // zprime, analytic distribution 
  newTree->Branch("p1_VAJHU",&p1_VAJHU,"p1_VAJHU/F");  // zprime, vector algebra, JHUgen,
  newTree->Branch("p2_mela ",&p2_mela ,"p2_mela/F");  // graviton, analytic distribution 
  newTree->Branch("p2_VAJHU",&p2_VAJHU,"p2_VAJHU/F");  // graviton, vector algebra, JHUgen,
  //backgrounds
  newTree->Branch("bkg_mela",&bkg_mela,"bkg_mela/F");  // background,  analytic distribution 
  newTree->Branch("bkg_VAMCFM",&bkg_VAMCFM,"bkg_VAMCFM/F");  // background, vector algebra, MCFM
  newTree->Branch("bkg_VAMCFMNorm",&bkg_VAMCFMNorm,"bkg_VAMCFMNorm/F");  // background, vector algebra, MCFM Normalized
  //pt/rapidity
  newTree->Branch("p0_pt",&p0_pt,"p0_pt/F");  // multiplicative probability for signal pt
  newTree->Branch("p0_y",&p0_y,"p0_y/F");  // multiplicative probability for signal y
  newTree->Branch("bkg_pt",&bkg_pt,"bkg_pt/F");  // multiplicative probability for bkg pt
  newTree->Branch("bkg_y",&bkg_y,"bkg_y/F");  // multiplicative probability for bkg y
  newTree->Branch("VAKD",&VAKD,"VAKD/F");  // discriminant
  
  for(int iEvt=0; iEvt<(max<0?sigTree->GetEntries():max); iEvt++){
    
    if(iEvt>=sigTree->GetEntries()) break;
    
    if(iEvt%1000==0) 
      cout << "event: " << iEvt << endl;
    
    sigTree->GetEntry(iEvt);
    
    // --------------------------------

    if(mzz>100.){
      myMELA_ICHEP.computeKD(mzz,m1,m2,hs,h1,h2,phi,phi1,D,psig,pbkg,false,pt4l,false,Y4l); //keep legacy probabilities
 
      myMELA.computeP(mzz, m1, m2, 
		      hs,h1,h2,phi,phi1,
		      //signal probabilities
		      p0plus_melaNorm,   // higgs, analytic distribution, normalized       
		      p0plus_mela,   // higgs, analytic distribution          
		      p0minus_mela,  // pseudoscalar, analytic distribution 
		      p0plus_VAJHU,  // higgs, vector algebra, JHUgen
		      p0minus_VAJHU, // pseudoscalar, vector algebra, JHUgen
		      p0plus_VAMCFM,// higgs, vector algebra, MCFM
		      p1_mela,  // zprime, analytic distribution 
		      p1_VAJHU, // zprime, vector algebra, JHUgen,
		      p2_mela , // graviton, analytic distribution 
		      p2_VAJHU, // graviton, vector algebra, JHUgen,
		      //backgrounds
		      bkg_mela,  // background,  analytic distribution 
		      bkg_VAMCFM, // background, vector algebra, MCFM
		      bkg_VAMCFMNorm, // background, vector algebra, MCFM Norm
		      //pt/rapidity
		      p0_pt, // multiplicative probability for signal pt
		      p0_y, // multiplicative probability for signal y
		      bkg_pt, // multiplicative probability for bkg pt
		      bkg_y, // multiplicative probability for bkg y
		      //optional input parameters
		      pt4l,Y4l,flavor // 1:4e, 2:4mu, 3:2e2mu (for interference effects)
		      );

      VAKD = p0plus_VAJHU/( bkg_VAMCFMNorm +  p0plus_VAJHU );
      
// //       std::cout << "Gravi "  <<  p0plus_mela/(p0plus_mela +  p2_mela) << " " <<  p0plus_mela  << " " <<p2_mela  <<std::endl;
// //       std::cout << "Pseudo " << p0plus_mela /(p0plus_mela +  p0minus_mela) << " " <<p0plus_mela  << " " <<p0minus_mela  <<std::endl;
// //       std::cout << "ICHEP "  << psig /(psig +  pbkg) << " "<< p0plus_mela  << " " <<p0minus_mela  <<std::endl;
// //       std::cout << "PTY "    << p0_pt*p0_y*p0plus_melaNorm /(p0_pt*p0_y*p0plus_melaNorm +  bkg_pt*bkg_y*bkg_mela) << " " <<p0_pt*p0_y*p0plus_melaNorm  << " " <<  bkg_pt*bkg_y*bkg_mela <<std::endl;
// //       std::cout << "Nominal "<< p0plus_melaNorm /(p0plus_melaNorm +  bkg_mela) << " " <<p0plus_melaNorm  << " " <<  bkg_mela <<std::endl;

// //       std::cout << "ME "     << p0plus_VAJHU/(p0plus_VAJHU + 10.*bkg_VAMCFM) << " " << p0plus_VAJHU <<" " << bkg_VAMCFM<<std::endl;
// //       std::cout << "pME "     << p0plus_VAJHU/(p0plus_VAJHU + 6.*p0minus_VAJHU) << " " << p0plus_VAJHU <<" " <<p0minus_VAJHU <<std::endl;
// //       std::cout << "graviME "     << p0plus_VAJHU/(p0plus_VAJHU + 1.2*p2_VAJHU) << " " << p0plus_VAJHU <<" " << p2_VAJHU<<std::endl;

	
      newTree->Fill();
      
    }

   }

  newFile->cd();
  newTree->Write("SelectedTree"); 
  newFile->Close();

}

