#ifndef ZZMatrixElement_ZZMatrixElement_h
#define ZZMatrixElement_ZZMatrixElement_h

#include <vector>
#include <TLorentzVector.h>
#include "ZZMatrixElement/MELA/interface/TVar.hh"
#include "ZZMatrixElement/MELA/interface/TEvtProb.hh"

class  ZZMatrixElement{
public:
  //pathtoHiggsCSandWidth: path to the textfiles of the HiggsCSandWidth package
  ZZMatrixElement(const char* pathtoHiggsCSandWidth,double ebeam,TVar::VerbosityLevel verbosity);
  ~ZZMatrixElement(){};
  /// Compute angles and KD from the lepton four-vectors and pdgIds.
  /// FSR recollected photons must been added to the corresponding lepton's four-vector by the user.
  /// Theta1 is the angle corresponding to Z1.
  /// Z1_lept1 and  Z1_lept2 are supposed to come from the same Z.
  /// Leptons are re-ordered internally according to a standard convention:
  /// lept1 = negative-charged lepton (for OS pairs).
  void computeXS(TLorentzVector Z1_lept1, int Z1_lept1Id,
                 TLorentzVector Z1_lept2, int Z1_lept2Id,
                 TLorentzVector Z2_lept1, int Z2_lept1Id,
                 TLorentzVector Z2_lept2, int Z2_lept2Id,
                 // return variables:
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
		 double& graviME);
    
  /// Compute KD from masses and angles.
  /// The user must ensure that the order of m1/m2 matches the order of theta1/theta2.
  // flavor 1 for 4e, 2 for 4m, 3 for 2e2mu
  void computeXS(float mZZ, float mZ1, float mZ2,
                 float costhetastar,
                 float costheta1,
                 float costheta2,
                 float phistar1,
                 float phi,
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
		 double& graviME);
  //
  // The below two functions calculate the ME integrating out the production angles
  // 
  void computeXSDecay(TLorentzVector Z1_lept1, int Z1_lept1Id,
                 TLorentzVector Z1_lept2, int Z1_lept2Id,
                 TLorentzVector Z2_lept1, int Z2_lept1Id,
                 TLorentzVector Z2_lept2, int Z2_lept2Id,
                 // return variables:
    	         double& dXsec_VZZ_DECAY_JHU,
    	         double& dXsec_AVZZ_DECAY_JHU,
   	         double& dXsec_TZZ_DECAY_JHU);
    
  /// Compute KD from masses and angles.
  /// The user must ensure that the order of m1/m2 matches the order of theta1/theta2.
  // flavor 1 for 4e, 2 for 4m, 3 for 2e2mu
  void computeXSDecay(float mZZ, float mZ1, float mZ2,
                 float costhetastar,
                 float costheta1,
                 float costheta2,
                 float phistar1,
                 float phi,
		 int flavor,
		 //output variables
    	         double& dXsec_VZZ_DECAY_JHU,
    	         double& dXsec_AVZZ_DECAY_JHU,
                 double& dXsec_TZZ_DECAY_JHU);
  
  //compute four-momenta from angles only 
  // Nota bene: angles, not cos(theta)...
  std::vector<TLorentzVector> Calculate4Momentum(double Mx,double M1,double M2,double theta,double theta1,double theta2,double Phi1,double Phi);

private:
  
  TVar::VerbosityLevel verb;
  TEvtProb Xcal2;
  hzz4l_event_type hzz4l_event;
  double EBEAM;

};
#endif
