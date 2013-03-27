#include "ZZMatrixElement/MELA/interface/newZZMatrixElement.h"
#include "TLorentzRotation.h"

newZZMatrixElement::newZZMatrixElement(const char* pathtoHiggsCSandWidth,
				       double ebeam,
				       TVar::VerbosityLevel verbosity):

  verb(verbosity),
  Xcal2(pathtoHiggsCSandWidth,ebeam),
  hzz4l_event(),
  EBEAM(ebeam)
				       {
  //std::cout << "TEST" << std::endl;
}

std::vector<TLorentzVector> newZZMatrixElement::Calculate4Momentum(double Mx,double M1,double M2,double theta,double theta1,double theta2,double Phi1,double Phi){
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


void newZZMatrixElement::computeXS(float mZZ, float mZ1, float mZ2,
				float costhetastar,
				float costheta1,
				float costheta2,
				float phi,
				float phistar1,
				int flavor,
  			        TVar::Process myModel,
  			        TVar::MatrixElement myME,
  			        TVar::Production myProduction,
  			        float &mevalue
				){
  std::vector<TLorentzVector> p;
  p=Calculate4Momentum(mZZ,mZ1,mZ2,acos(costhetastar),acos(costheta1),acos(costheta2),phistar1,phi);
  TLorentzVector Z1_lept1 = p[0];
  TLorentzVector Z1_lept2 = p[1];
  TLorentzVector Z2_lept1 = p[2];
  TLorentzVector Z2_lept2 = p[3];
  
  int Z1_lept1Id, Z1_lept2Id, Z2_lept1Id, Z2_lept2Id; 
  assert(flavor==1 || flavor==2 || flavor==3 );
  if ( flavor == 1 ) {
    Z1_lept1Id = 11;
    Z1_lept2Id = -11;
    Z2_lept1Id = 11;
    Z2_lept2Id = -11;
  }
  if ( flavor == 2 ) {
    Z1_lept1Id = 13;
    Z1_lept2Id = -13;
    Z2_lept1Id = 13;
    Z2_lept2Id = -13;
  }
  if ( flavor == 3 ) {
    Z1_lept1Id = 11;
    Z1_lept2Id = -11;
    Z2_lept1Id = 13;
    Z2_lept2Id = -13;
  }
  
  hzz4l_event.p[0].SetXYZM(Z1_lept1.Px(), Z1_lept1.Py(), Z1_lept1.Pz(), 0.);
  hzz4l_event.p[1].SetXYZM(Z1_lept2.Px(), Z1_lept2.Py(), Z1_lept2.Pz(), 0.);
  hzz4l_event.p[2].SetXYZM(Z2_lept1.Px(), Z2_lept1.Py(), Z2_lept1.Pz(), 0.);
  hzz4l_event.p[3].SetXYZM(Z2_lept2.Px(), Z2_lept2.Py(), Z2_lept2.Pz(), 0.);
  
  hzz4l_event.PdgCode[0] =  Z1_lept1Id;
  hzz4l_event.PdgCode[1] =  Z1_lept2Id;
  hzz4l_event.PdgCode[2] =  Z2_lept1Id;
  hzz4l_event.PdgCode[3] =  Z2_lept2Id;
  
  double z1mass = (hzz4l_event.p[0]+hzz4l_event.p[1]).M();
  double z2mass = (hzz4l_event.p[2]+hzz4l_event.p[3]).M();
  double zzmass = (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).M();
  
  if (verb >= TVar::DEBUG) {
    std::cout << "four lepton p4 for ME calculations: ===========================" <<endl;
    std::cout << "PDG code \n";
    std::cout << "hzz4l_event.PdgCode[0] = " << hzz4l_event.PdgCode[0] << "\n";
    std::cout << "hzz4l_event.PdgCode[1] = " << hzz4l_event.PdgCode[1] << "\n";
    std::cout << "hzz4l_event.PdgCode[2] = " << hzz4l_event.PdgCode[2] << "\n";
    std::cout << "hzz4l_event.PdgCode[3] = " << hzz4l_event.PdgCode[3] << "\n";
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
    std::cout << "=========================================================\n";
  }
  
  // ==== Begin the differential cross-section calculation
  Xcal2.SetHiggsMass(zzmass);
  Xcal2.SetMatrixElement(myME);
  Xcal2.SetProduction(myProduction);
  mevalue  = Xcal2.XsecCalc(myModel,myProduction,hzz4l_event,verb);
  
  return;
}
