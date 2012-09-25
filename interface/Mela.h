#ifndef MELA_Mela_h
#define MELA_Mela_h

/** \class Mela
 *
 *  MELA discriminator
 *
 *
 *  $Date: 2012/09/24 12:15:07 $
 *  $Revision: 1.8 $
 *  \author JHU
 */

#include <vector>
#include <TLorentzVector.h>

class TFile; 
class TH1F; 
class TH2F;
class TH3F;


class Mela { 
public:
  Mela(bool usePowhegTemplate=false);

  ~Mela();

  /// Compute angles and KD from the lepton four-vectors and pdgIds.
  /// FSR recollected photons must been added to the corresponding lepton's four-vector by the user.
  /// Theta1 is the angle corresponding to Z1.
  /// Z1_lept1 and  Z1_lept2 are supposed to come from the same Z.
  /// Leptons are re-ordered internally according to a standard convention:
  /// lept1 = negative-charged lepton (for OS pairs).
  void computeKD(TLorentzVector Z1_lept1, int Z1_lept1Id,
		 TLorentzVector Z1_lept2, int Z1_lept2Id,
		 TLorentzVector Z2_lept1, int Z2_lept1Id,
		 TLorentzVector Z2_lept2, int Z2_lept2Id,
		 // return variables:
		 float& costhetastar,
		 float& costheta1, 
		 float& costheta2,
		 float& phi,
		 float& phistar1,
		 float& kd, 
		 float& psig,
		 float& pbkg,
		 bool withPt = false,
		 bool withY = false,
		 int LHCsqrts = 8);



  /// Compute KD from masses and angles. 
  /// The user must ensure that the order of m1/m2 matches the order of theta1/theta2.
  void computeKD(float mZZ, float mZ1, float mZ2, 
		 float costhetastar,
		 float costheta1, 
		 float costheta2,
		 float phi,
		 float phistar1,
		 float& kd, 
		 float& psig,
		 float& pbkg,
		 bool withPt = false,
		 float pt4l=0.0,
		 bool withY = false,
		 float Y4l=0.0,
		 int LHCsqrts=8);

private:

  double sigPdfNorm(double mzz);
  double bkgPdfNorm(double mzz);


  std::vector<float> my8DTemplate(bool normalized, float mZZ, float mZ1, float mZ2, 
				   float costhetastar, 
				   float costheta1, 
				   float costheta2, 
				   float phi, 
				   float phistar1);

  std::pair<float,float> likelihoodDiscriminant (float mZZ, float mZ1, float mZ2, 
						   float costhetastar, 
						   float costheta1, 
						   float costheta2, 
						   float phi, 
						   float phistar1,
						   int LHCsqrts=8, 
						   bool withPt = false, float pt = 0.0, 
						   bool withY = false, float y = 0.0);

  // Ordering of Z1/Z2 according to internal convention
  void checkZorder(float& z1mass, float& z2mass, float& costhetastar, float& costheta1, float& costheta2, float& phi, float& phistar1);

  TFile *f;
  TH1F *h_mzz;
  TH3F *h_mzzm1m2;
  TH2F *h_mzzcosthetastar;
  TH2F *h_mzzcostheta1;
  TH2F *h_mzzcostheta2;
  TH2F *h_mzzphi1;
  TH2F *h_mzzphi;
  
  bool usePowhegTemplate_;
 
};



#endif
