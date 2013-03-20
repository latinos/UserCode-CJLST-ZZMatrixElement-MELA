#include "ZZMatrixElement/MELA/interface/ZZMatrixElement.h"
#include "TLorentzRotation.h"

ZZMatrixElement::ZZMatrixElement(const char* pathtoHiggsCSandWidth,double ebeam,TVar::VerbosityLevel verbosity):
  verb(verbosity),
  Xcal2(pathtoHiggsCSandWidth,ebeam),
  hzz4l_event(),
  EBEAM(ebeam)
{
  //std::cout << "TEST" << std::endl;
}

std::vector<TLorentzVector> ZZMatrixElement::Calculate4Momentum(double Mx,double M1,double M2,double theta,double theta1,double theta2,double Phi1,double Phi){
    double phi1,phi2;
    phi1=TMath::Pi()-Phi1;
    phi2=Phi1+Phi;

    double gamma1,gamma2,beta1,beta2;

    gamma1=(Mx*Mx+M1*M1-M2*M2)/(2*Mx*M1);
    gamma2=(Mx*Mx-M1*M1+M2*M2)/(2*Mx*M2);
    beta1=sqrt(1-1/(gamma1*gamma1));
    beta2=sqrt(1-1/(gamma2*gamma2));

    //gluon 4 vectors
    TLorentzVector p1CM(0,0,Mx/2,Mx/2);
    TLorentzVector p2CM(0,0,-Mx/2,Mx/2);

    //vector boson 4 vectors
    TLorentzVector kZ1(gamma1*M1*sin(theta)*beta1,0, gamma1*M1*cos(theta)*beta1,gamma1*M1*1);
    TLorentzVector kZ2(-gamma2*M2*sin(theta)*beta2,0, -gamma2*M2*cos(theta)*beta2,gamma2*M2*1);

    //Rotation and Boost matrices. Note gamma1*beta1*M1=gamma2*beta2*M2.

    TLorentzRotation Z1ToZ,Z2ToZ;

    Z1ToZ.Boost(0,0,beta1);
    Z2ToZ.Boost(0,0,beta2);
    Z1ToZ.RotateY(theta);
    Z2ToZ.RotateY(TMath::Pi()+theta);


    //fermion 4 vectors in vector boson rest frame

    TLorentzVector p3Z1((M1/2)*sin(theta1)*cos(phi1),(M1/2)*sin(theta1)*sin(phi1),(M1/2)*cos(theta1),(M1/2)*1);
    TLorentzVector p4Z1(-(M1/2)*sin(theta1)*cos(phi1),-(M1/2)*sin(theta1)*sin(phi1),-(M1/2)*cos(theta1),(M1/2)*1);
    TLorentzVector p5Z2((M2/2)*sin(theta2)*cos(phi2),(M2/2)*sin(theta2)*sin(phi2),(M2/2)*cos(theta2),(M2/2)*1);
    TLorentzVector p6Z2(-(M2/2)*sin(theta2)*cos(phi2),-(M2/2)*sin(theta2)*sin(phi2),-(M2/2)*cos(theta2),(M2/2)*1);

    // fermions 4 vectors in CM frame

    TLorentzVector p3CM,p4CM,p5CM,p6CM;

    p3CM=Z1ToZ*p3Z1;
    p4CM=Z1ToZ*p4Z1;
    p5CM=Z2ToZ*p5Z2;
    p6CM=Z2ToZ*p6Z2;

    vector<TLorentzVector> p;

    p.push_back(p3CM);
    p.push_back(p4CM);
    p.push_back(p5CM);
    p.push_back(p6CM);

    return p;
}

void ZZMatrixElement::computeXS(TLorentzVector Z1_lept1, int Z1_lept1Id,
	       TLorentzVector Z1_lept2, int Z1_lept2Id,
	       TLorentzVector Z2_lept1, int Z2_lept1Id,
	       TLorentzVector Z2_lept2, int Z2_lept2Id,
	       //output variables
	       double& dXsec_ZZ_MCFM,
	       double& dXsec_GGZZ_MCFM,
	       double& dXsec_HZZ_MCFM,
	       double& dXsec_HZZ_JHU,
	       double& dXsec_PSHZZ_JHU,
	       double& dXsec_VZZ_JHU,
	       double& dXsec_TZZ_JHU,
	       double& dXsec_AVZZ_JHU,
	       double& dXsec_QQB_TZZ_JHU,
	       double& dXsec_HDHZZ_JHU,
	       double& ME,
	       double& pseudoME,
	       double& graviME
	       ){
  
  //
  // initialise the differential cross-sections
  //
  dXsec_ZZ_MCFM = 0.;
  dXsec_GGZZ_MCFM = 0.;
  dXsec_HZZ_MCFM = 0.;
  dXsec_HZZ_JHU = 0.;
  dXsec_PSHZZ_JHU = 0.;
  dXsec_HDHZZ_JHU = 0.;
  dXsec_VZZ_JHU = 0.;
  dXsec_AVZZ_JHU = 0.;
  dXsec_TZZ_JHU = 0.;
  dXsec_QQB_TZZ_JHU = 0.;
  ME = 0.;
  pseudoME = 0.;
  graviME = 0.;
  

  hzz4l_event.p[0].SetXYZM(Z1_lept1.Px(), Z1_lept1.Py(), Z1_lept1.Pz(), 0.);
  hzz4l_event.p[1].SetXYZM(Z1_lept2.Px(), Z1_lept2.Py(), Z1_lept2.Pz(), 0.);
  hzz4l_event.p[2].SetXYZM(Z2_lept1.Px(), Z2_lept1.Py(), Z2_lept1.Pz(), 0.);
  hzz4l_event.p[3].SetXYZM(Z2_lept2.Px(), Z2_lept2.Py(), Z2_lept2.Pz(), 0.);

  double z1mass = (hzz4l_event.p[0]+hzz4l_event.p[1]).M();
  double z2mass = (hzz4l_event.p[2]+hzz4l_event.p[3]).M();
  double zzmass = (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).M();
  

  if (verb >= TVar::DEBUG) {
    cout << "four lepton p4 for ME calculations: ===========================" <<endl;
    printf("lep1 p3 = (%4.4f, %4.4f, %4.4f)  lep2 p3 = (%4.4f, %4.4f, %4.4f)\n",
	   Z1_lept1.Px(), Z1_lept1.Py(), Z1_lept1.Pz(), Z1_lept2.Px(), Z1_lept2.Py(), Z1_lept2.Pz());
    printf("lep3 p3 = (%4.4f, %4.4f, %4.4f)  lep4 p3 = (%4.4f, %4.4f, %4.4f)\n",
	   Z2_lept1.Px(), Z2_lept1.Py(), Z2_lept1.Pz(), Z2_lept2.Px(), Z2_lept2.Py(), Z2_lept2.Pz());
    std::cout << "ZZ system (pX, pY, pZ, E, mass) = ( "
	      << (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Px() << ", "
	      << (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Py() << ", "
	      << (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Pz() << ", "
	      << (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Energy()  << ", "
	      << zzmass << ")\n";
    std::cout << "Z1 mass = " << z1mass << "\tz2mass = " << z2mass << "\n";
    cout << "=========================================================\n";
  }

  hzz4l_event.PdgCode[0] =  Z1_lept1Id;
  hzz4l_event.PdgCode[1] =  Z1_lept2Id;
  hzz4l_event.PdgCode[2] =  Z2_lept1Id;
  hzz4l_event.PdgCode[3] =  Z2_lept2Id;

  if (verb >= TVar::DEBUG) {
    std::cout << "PDG code \n";
    std::cout << "hzz4l_event.PdgCode[0] = " << hzz4l_event.PdgCode[0] << "\n";
    std::cout << "hzz4l_event.PdgCode[1] = " << hzz4l_event.PdgCode[1] << "\n";
    std::cout << "hzz4l_event.PdgCode[2] = " << hzz4l_event.PdgCode[2] << "\n";
    std::cout << "hzz4l_event.PdgCode[3] = " << hzz4l_event.PdgCode[3] << "\n";
  }
  
  // ==== Begin the differential cross-section calculation
  Xcal2.SetHiggsMass(zzmass);
  
  //
  // MCFM based calculations
  //
  Xcal2.SetMatrixElement(TVar::MCFM);
  if(abs(Z1_lept1Id)==abs(Z2_lept1Id)){// 4e or 4mu
    dXsec_ZZ_MCFM = Xcal2.XsecCalc(TVar::ZZ_4e,hzz4l_event,verb);
  }
  else //2e2mu
    dXsec_ZZ_MCFM = Xcal2.XsecCalc(TVar::ZZ_2e2m,hzz4l_event,verb);
  // ggZZ
  dXsec_GGZZ_MCFM = Xcal2.XsecCalc(TVar::GGZZ_4l,hzz4l_event,verb);
  dXsec_HZZ_MCFM = Xcal2.XsecCalc(TVar::HZZ_4l,hzz4l_event,verb);

  //
  // JHUGen based calcualtions
  //
  // 0+
  Xcal2.SetMatrixElement(TVar::JHUGen);
  dXsec_HZZ_JHU = Xcal2.XsecCalc(TVar::HZZ_4l,hzz4l_event,verb);
  ME =  dXsec_HZZ_JHU / ( dXsec_HZZ_JHU + 10. * dXsec_ZZ_MCFM );
  
  // 0-
  dXsec_PSHZZ_JHU = Xcal2.XsecCalc(TVar::PSHZZ_4l,hzz4l_event,verb);
  pseudoME =  dXsec_HZZ_JHU / ( dXsec_HZZ_JHU + 6*dXsec_PSHZZ_JHU );

  // 0h+
  dXsec_HDHZZ_JHU = Xcal2.XsecCalc(TVar::HDHZZ_4l,hzz4l_event,verb);

  // 1-
  dXsec_VZZ_JHU = Xcal2.XsecCalc(TVar::VZZ_4l,hzz4l_event,verb);

  // 1+
  dXsec_AVZZ_JHU = Xcal2.XsecCalc(TVar::AVZZ_4l,hzz4l_event,verb);
  
  // 2m+
  dXsec_TZZ_JHU = Xcal2.XsecCalc(TVar::TZZ_4l,hzz4l_event,verb);
  graviME =  dXsec_HZZ_JHU / ( dXsec_HZZ_JHU + 1.2*dXsec_TZZ_JHU );

  // 2m+ by qqbar
  dXsec_QQB_TZZ_JHU = Xcal2.XsecCalc(TVar::QQB_TZZ_4l,hzz4l_event,verb);

  
  return;
}

void ZZMatrixElement::computeXS(float mZZ, float mZ1, float mZ2,
				float costhetastar,
				float costheta1,
				float costheta2,
				float phi,
				float phistar1,
				int flavor,
				//output variables
				double& dXsec_ZZ_MCFM,
				double& dXsec_GGZZ_MCFM,
				double& dXsec_HZZ_MCFM,
				double& dXsec_HZZ_JHU,
				double& dXsec_PSHZZ_JHU,
				double& dXsec_VZZ_JHU,
				double& dXsec_TZZ_JHU,
				double& dXsec_AVZZ_JHU,
				double& dXsec_QQB_TZZ_JHU,
				double& dXsec_HDHZZ_JHU,
				double& ME,
				double& pseudoME,
				double& graviME
				){
  std::vector<TLorentzVector> p;
  p=Calculate4Momentum(mZZ,mZ1,mZ2,acos(costhetastar),acos(costheta1),acos(costheta2),phistar1,phi);
  
  int pdg11,pdg12,pdg21,pdg22;
  assert(flavor==1 || flavor==2 || flavor==3 );
  if ( flavor == 1 ) {
    pdg11 = 11;
    pdg12 = -11;
    pdg21 = 11;
    pdg22 = -11;
  }
  if ( flavor == 2 ) {
    pdg11 = 13;
    pdg12 = -13;
    pdg21 = 13;
    pdg22 = -13;
  }
  if ( flavor == 3 ) {
    pdg11 = 11;
    pdg12 = -11;
    pdg21 = 13;
    pdg22 = -13;
  }
 


  computeXS(p[0], pdg11,
	    p[1], pdg12,
	    p[2], pdg21,
	    p[3], pdg22,
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
	    dXsec_HDHZZ_JHU,
	    ME,
	    pseudoME,
	    graviME);
  
  return;
}

void ZZMatrixElement::computeXSExtra(TLorentzVector Z1_lept1, int Z1_lept1Id,
	       TLorentzVector Z1_lept2, int Z1_lept2Id,
	       TLorentzVector Z2_lept1, int Z2_lept1Id,
	       TLorentzVector Z2_lept2, int Z2_lept2Id,
	       //output variables
   	       double& dXsec_VZZ_DECAY_JHU,
   	       double& dXsec_AVZZ_DECAY_JHU,
 	       double& dXsec_TZZ_DECAY_JHU,
   	       double& dXsec_PTZZ_2hminus_JHU,
   	       double& dXsec_TZZ_2hplus_JHU,
               double& dXsec_TZZ_2bplus_JHU
	       ){
  
  //
  // initialise the differential cross-sections
  //
  dXsec_VZZ_DECAY_JHU = 0.;
  dXsec_AVZZ_DECAY_JHU = 0.;
  dXsec_TZZ_DECAY_JHU = 0.;
  dXsec_PTZZ_2hminus_JHU = 0.;
  dXsec_TZZ_2hplus_JHU = 0.;
  dXsec_TZZ_2bplus_JHU = 0.;
  
  hzz4l_event.p[0].SetXYZM(Z1_lept1.Px(), Z1_lept1.Py(), Z1_lept1.Pz(), 0.);
  hzz4l_event.p[1].SetXYZM(Z1_lept2.Px(), Z1_lept2.Py(), Z1_lept2.Pz(), 0.);
  hzz4l_event.p[2].SetXYZM(Z2_lept1.Px(), Z2_lept1.Py(), Z2_lept1.Pz(), 0.);
  hzz4l_event.p[3].SetXYZM(Z2_lept2.Px(), Z2_lept2.Py(), Z2_lept2.Pz(), 0.);

  double z1mass = (hzz4l_event.p[0]+hzz4l_event.p[1]).M();
  double z2mass = (hzz4l_event.p[2]+hzz4l_event.p[3]).M();
  double zzmass = (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).M();
  

  if (verb >= TVar::DEBUG) {
    cout << "four lepton p4 for ME calculations: ===========================" <<endl;
    printf("lep1 p3 = (%4.4f, %4.4f, %4.4f)  lep2 p3 = (%4.4f, %4.4f, %4.4f)\n",
	   Z1_lept1.Px(), Z1_lept1.Py(), Z1_lept1.Pz(), Z1_lept2.Px(), Z1_lept2.Py(), Z1_lept2.Pz());
    printf("lep3 p3 = (%4.4f, %4.4f, %4.4f)  lep4 p3 = (%4.4f, %4.4f, %4.4f)\n",
	   Z2_lept1.Px(), Z2_lept1.Py(), Z2_lept1.Pz(), Z2_lept2.Px(), Z2_lept2.Py(), Z2_lept2.Pz());
    std::cout << "ZZ system (pX, pY, pZ, E, mass) = ( "
	      << (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Px() << ", "
	      << (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Py() << ", "
	      << (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Pz() << ", "
	      << (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Energy()  << ", "
	      << zzmass << ")\n";
    std::cout << "Z1 mass = " << z1mass << "\tz2mass = " << z2mass << "\n";
    cout << "=========================================================\n";
  }

  hzz4l_event.PdgCode[0] =  Z1_lept1Id;
  hzz4l_event.PdgCode[1] =  Z1_lept2Id;
  hzz4l_event.PdgCode[2] =  Z2_lept1Id;
  hzz4l_event.PdgCode[3] =  Z2_lept2Id;

  if (verb >= TVar::DEBUG) {
    std::cout << "PDG code \n";
    std::cout << "hzz4l_event.PdgCode[0] = " << hzz4l_event.PdgCode[0] << "\n";
    std::cout << "hzz4l_event.PdgCode[1] = " << hzz4l_event.PdgCode[1] << "\n";
    std::cout << "hzz4l_event.PdgCode[2] = " << hzz4l_event.PdgCode[2] << "\n";
    std::cout << "hzz4l_event.PdgCode[3] = " << hzz4l_event.PdgCode[3] << "\n";
  }
  
  // ==== Begin the differential cross-section calculation
  Xcal2.SetHiggsMass(zzmass);
  Xcal2.SetMatrixElement(TVar::JHUGen);
  
  // 1- decay 
  dXsec_VZZ_DECAY_JHU = Xcal2.XsecCalc(TVar::VZZ_DECAY_4l,hzz4l_event,verb);  

  // 1+ decay 
  dXsec_AVZZ_DECAY_JHU = Xcal2.XsecCalc(TVar::AVZZ_DECAY_4l,hzz4l_event,verb);  

  // 2m+ decay 
  dXsec_TZZ_DECAY_JHU = Xcal2.XsecCalc(TVar::TZZ_DECAY_4l,hzz4l_event,verb); 

  // 2h-
  dXsec_PTZZ_2hminus_JHU = Xcal2.XsecCalc(TVar::PTZZ_2hminus_4l,hzz4l_event,verb);  

  // 2h+
  dXsec_TZZ_2hplus_JHU = Xcal2.XsecCalc(TVar::TZZ_2hplus_4l,hzz4l_event,verb);  

  // 2b+ 
  dXsec_TZZ_2bplus_JHU = Xcal2.XsecCalc(TVar::TZZ_2bplus_4l,hzz4l_event,verb);   
  
  return;
}

void ZZMatrixElement::computeXSExtra(float mZZ, float mZ1, float mZ2,
				float costhetastar,
				float costheta1,
				float costheta2,
				float phi,
				float phistar1,
				int flavor,
				//output variables
   			        double& dXsec_VZZ_DECAY_JHU,
				double& dXsec_AVZZ_DECAY_JHU,
				double& dXsec_TZZ_DECAY_JHU,
  			        double& dXsec_PTZZ_2hminus_JHU,
				double& dXsec_TZZ_2hplus_JHU,
                                double& dXsec_TZZ_2bplus_JHU
				){
  std::vector<TLorentzVector> p;
  p=Calculate4Momentum(mZZ,mZ1,mZ2,acos(costhetastar),acos(costheta1),acos(costheta2),phistar1,phi);
  
  int pdg11,pdg12,pdg21,pdg22;
  assert(flavor==1 || flavor==2 || flavor==3 );
  if ( flavor == 1 ) {
    pdg11 = 11;
    pdg12 = -11;
    pdg21 = 11;
    pdg22 = -11;
  }
  if ( flavor == 2 ) {
    pdg11 = 13;
    pdg12 = -13;
    pdg21 = 13;
    pdg22 = -13;
  }
  if ( flavor == 3 ) {
    pdg11 = 11;
    pdg12 = -11;
    pdg21 = 13;
    pdg22 = -13;
  }
  
  computeXSExtra(p[0], pdg11,
	    p[1], pdg12,
	    p[2], pdg21,
	    p[3], pdg22,
	    //output variables
      	    dXsec_VZZ_DECAY_JHU,
            dXsec_AVZZ_DECAY_JHU,
            dXsec_TZZ_DECAY_JHU,
      	    dXsec_PTZZ_2hminus_JHU,
      	    dXsec_TZZ_2hplus_JHU,
	    dXsec_TZZ_2bplus_JHU
	    );
  
  return;
}

