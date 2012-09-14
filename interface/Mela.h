#ifndef MELA_Mela_h
#define MELA_Mela_h

/** \class Mela
 *
 *  MELA 
 *
 *
 *  $Date: 2012/08/19 20:46:50 $
 *  $Revision: 1.4 $
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
  void computeLD(TLorentzVector Z1L1,
		 TLorentzVector Z1L2,
		 TLorentzVector Z2L1,
		 TLorentzVector Z2L2,		 
		 // return variables:
		 double& costhetastar,
		 double& costheta1, 
		 double& costheta2,
		 double& phi,
		 double& phistar1,
		 float& ld, 
		 float& psig,
		 float& pbkg);



  /// Compute LD from masses and angles. 
  /// Note that the order of m1/m2 must be consistent with the order of theta1/theta2
  void computeLD(double mZZ, double m1, double m2, 
		 double costhetastar,
		 double costheta1, 
		 double costheta2,
		 double phi,
		 double phistar1,
		 float& ld, 
		 float& psig,
		 float& pbkg);

private:	
  std::vector<double> my8DTemplate(bool normalized, double mZZ, double m1, double m2, 
				   double costhetastar, 
				   double costheta1, 
				   double costheta2, 
				   double phi, 
				   double phistar1);

  std::pair<double,double> likelihoodDiscriminant (double mZZ, double m1, double m2, 
						   double costhetastar, 
						   double costheta1, 
						   double costheta2, 
						   double phi, 
						   double phistar1,
						   int LHCsqrts=8, 
						   bool withPt = false, double pt = 0.0, 
						   bool withY = false, double y = 0.0);

  // Set internal ordering of Z1/Z2
  void checkZorder(double& z1mass, double& z2mass, double& costhetastar, double& costheta1, double& costheta2, double& phi, double& phistar1);

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
