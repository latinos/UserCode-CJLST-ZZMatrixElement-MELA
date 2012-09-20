#ifndef MELA_Mela_h
#define MELA_Mela_h

/** \class Mela
 *
 *  MELA discriminator
 *
 *
 *  $Date: 2012/09/19 14:19:20 $
 *  $Revision: 1.5 $
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
  Mela();

  ~Mela();

  /// Compute angles and LD from the lepton four-vectors and pdgIds.
  /// Z1_lept1 and  Z1_lept2 are supposed to come from the same Z.
  /// Zs and leptons are re-ordered internally according to a defined convention:
  /// Z1 = closest to nominal mZ; lept1 = negative-charged lepton (for OS pairs).
  /// FSR recollected photons must been added to the corresponding lepton's four-vector by the user.
  void computeLD(TLorentzVector Z1_lept1, int Z1_lept1Id,
		 TLorentzVector Z1_lept2, int Z1_lept2Id,
		 TLorentzVector Z2_lept1, int Z2_lept1Id,
		 TLorentzVector Z2_lept2, int Z2_lept2Id,
		 // return variables:
		 float& costhetastar,
		 float& costheta1, 
		 float& costheta2,
		 float& phi,
		 float& phistar1,
		 float& ld, 
		 float& psig,
		 float& pbkg,
		 bool withPt = false,
		 bool withY = false,
		 int LHCsqrts = 8);



  /// Compute LD from masses and angles. It is assumed that the order of theta1/theta2 is the one defined 
  /// in the method above.
  void computeLD(float mZZ, float mZ1, float mZ2, 
		 float costhetastar,
		 float costheta1, 
		 float costheta2,
		 float phi,
		 float phistar1,
		 float& ld, 
		 float& psig,
		 float& pbkg,
		 bool withPt = false,
		 float pt4l=0.0,
		 bool withY = false,
		 float Y4l=0.0,
		 int LHCsqrts=8);

private:
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

  // Set internal ordering of Z1/Z2
  void checkZorder(float& z1mass, float& z2mass, float& costhetastar, float& costheta1, float& costheta2, float& phi, float& phistar1);

  TFile *f;
  TH1F *h_mzz;
  TH3F *h_mzzm1m2;
  TH2F *h_mzzcosthetastar;
  TH2F *h_mzzcostheta1;
  TH2F *h_mzzcostheta2;
  TH2F *h_mzzphi1;
  TH2F *h_mzzphi;
 
};



#endif
