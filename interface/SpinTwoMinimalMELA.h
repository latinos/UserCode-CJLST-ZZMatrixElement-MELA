#ifndef MELA_SpinTwoMinimalMela_h
#define MELA_SpinTwoMinimalMela_h

/** \class SpinTwoMinimalMela
 *
 *  SpinTwoMinimalMELA discriminator 
 *
 *
 *  $Date: 2012/09/18 15:14:25 $
 *  $Revision: 1.1 $
 *  \author JHU
 */

#include <TLorentzVector.h>

class AngularPdfFactory;
class TensorPdfFactory;
class RooRealVar;


class SpinTwoMinimalMELA{

public:

  SpinTwoMinimalMELA();

  ~SpinTwoMinimalMELA(){};

  void eval(TLorentzVector Z1_lept1, int Z1_lept1Id,
	    TLorentzVector Z1_lept2, int Z1_lept2Id,
	    TLorentzVector Z2_lept1, int Z2_lept1Id,
	    TLorentzVector Z2_lept2, int Z2_lept2Id,
	    float& ld, 
	    float& psig,
	    float& pbkg);
  
  void eval(float zzmass, float z1mass, float z2mass, 
	    float costhetstar, 
	    float costheta1, 
	    float costheta2, 
	    float phi, 
	    float phistar1,
	    float& ld, 
	    float& psig, 
	    float& psigALT);

private:
  void checkZorder(float& z1mass, float& z2mass,
				   float& costhetastar, float& costheta1,
				   float& costheta2, float& phi,
				   float& phistar1);

  RooRealVar* z1mass_rrv;
  RooRealVar* z2mass_rrv;
  RooRealVar* costheta1_rrv;
  RooRealVar* costheta2_rrv;
  RooRealVar* phi_rrv;
  RooRealVar* costhetastar_rrv;
  RooRealVar* phistar1_rrv;
  RooRealVar* mzz_rrv;

  AngularPdfFactory *SMHiggs;
  TensorPdfFactory *minGrav;   

};

#endif
