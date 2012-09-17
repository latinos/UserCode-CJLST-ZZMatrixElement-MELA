#ifndef MELA_Mela_h
#define MELA_Mela_h

/** \class Mela
 *
 *  MELA 
 *
 *
 *  $Date: 2012/09/14 20:04:45 $
 *  $Revision: 1.2 $
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

  /// Compute angles and LD from the lepton four-vectors.
  void computeLD(TLorentzVector Z1_lept,
		 TLorentzVector Z1_antiLept,
		 TLorentzVector Z2_lept,
		 TLorentzVector Z2_antiLept,		 
		 // return variables:
		 float& costhetastar,
		 float& costheta1, 
		 float& costheta2,
		 float& phi,
		 float& phistar1,
		 float& ld, 
		 float& psig,
		 float& pbkg);



  /// Compute LD from masses and angles. 
  /// Note that the order of m1/m2 must be consistent with the order of theta1/theta2
  void computeLD(float mZZ, float m1, float m2, 
		 float costhetastar,
		 float costheta1, 
		 float costheta2,
		 float phi,
		 float phistar1,
		 float& ld, 
		 float& psig,
		 float& pbkg);

private:	
  std::vector<float> my8DTemplate(bool normalized, float mZZ, float m1, float m2, 
				   float costhetastar, 
				   float costheta1, 
				   float costheta2, 
				   float phi, 
				   float phistar1);

  std::pair<float,float> likelihoodDiscriminant (float mZZ, float m1, float m2, 
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
