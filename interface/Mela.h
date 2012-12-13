#ifndef MELA_Mela_h
#define MELA_Mela_h

/** \class Mela
 *
 *  MELA discriminator
 *
 *
 *  $Date: 2012/12/13 13:19:02 $
 *  $Revision: 1.18 $
 *  \author JHU
 */

#include <vector>
#include <TLorentzVector.h>

class TFile; 
class TH1F; 
class TH2F;
class TH3F;
class RooRealVar;
class RooAbsPdf;
class RooArgSet;
class AngularPdfFactory;
class TensorPdfFactory;
class ZZMatrixElement;
class SuperMELA;
class TGraph;


class Mela { 
public:
  Mela(bool usePowhegTemplate=false, int LHCsqrts=8, float mh=126); // higgs mass for supermela

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
		 bool withY = false);

  void computeWeight(float mZZ, float mZ1, float mZ2, 
		     float costhetastar,
		     float costheta1, 
		     float costheta2,
		     float phi,
		     float phistar1,
		     // return variables:
		     float& w
		     );

  void computeP(float mZZ, float mZ1, float mZ2, 
		float costhetastar,
		float costheta1, 
		float costheta2,
		float phi,
		float phistar1,
		//signal probabilities
		float& p0plus_melaNorm,   // higgs, analytic distribution, normalized as for normal MELA distribution     
		float& p0plus_mela,   // higgs, analytic distribution          
		float& p0minus_mela,  // pseudoscalar, analytic distribution 
		float& p0plus_VAJHU,  // higgs, vector algebra, JHUgen
		float& p0minus_VAJHU, // pseudoscalar, vector algebra, JHUgen
		float& p0plus_VAMCFM,// higgs, vector algebra, MCFM
		float& p1_mela,  // zprime, analytic distribution 
		float& p1_VAJHU, // zprime, vector algebra, JHUgen,
		float& p2_mela , // graviton, analytic distribution 
		float& p2_VAJHU, // graviton, vector algebra, JHUgen,
		//backgrounds
		float& bkg_mela,  // background,  analytic distribution 
		float& bkg_VAMCFM, // background, vector algebra, MCFM
		float& ggzz_VAMCFM, // background, vector algebra, MCFM for ggZZ
		float& bkg_VAMCFMNorm, // background, vector algebra, MCFM, Normalized 
		//pt/rapidity
		float& p0_pt, // multiplicative probability for signal pt
		float& p0_y, // multiplicative probability for signal y
		float& bkg_pt, // multiplicative probability for bkg pt
		float& bkg_y, // multiplicative probability for bkg y
		//supermela
		float& p0plus_m4l,  // signal m4l probability as in datacards
		float& bkg_m4l,     // backgroun m4l probability as in datacards
		//optional input parameters
		float pt4l=0.0,
		float Y4l=0.0,
		int flavor=1 // 1:4e, 2:4mu, 3:2e2mu (for interference effects)
		);

  
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
		 float Y4l=0.0);

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
 
  RooRealVar* z1mass_rrv;
  RooRealVar* z2mass_rrv;
  RooRealVar* costhetastar_rrv; 
  RooRealVar* costheta1_rrv;
  RooRealVar* costheta2_rrv;
  RooRealVar* phi_rrv;
  RooRealVar* phi1_rrv;
  RooRealVar* pt_rrv;
  RooRealVar* y_rrv;
  RooRealVar* sqrtS_rrv;
  RooRealVar* mzz_rrv;
  RooRealVar* upFrac_rrv;
  
  AngularPdfFactory* SMHiggs;
  AngularPdfFactory* PSHiggs;
  TensorPdfFactory *minGrav;   
  RooAbsPdf* SMZgammaZZ;
  RooAbsPdf* SMZZ;
  
  ZZMatrixElement* ZZME;
  SuperMELA* super;

  std::vector<RooRealVar*> ptparamsS;
  std::vector<RooRealVar*> ptparamsB;
  RooArgSet* allparamsS;
  RooArgSet* allparamsB;

  RooAbsPdf* sigY; 
  RooAbsPdf* bkgY;
  
  RooAbsPdf* sigPt;
  RooAbsPdf* bkgPt;

  TGraph* vaScale_4e;
  TGraph* vaScale_4mu;
  TGraph* vaScale_2e2mu;

};



#endif
