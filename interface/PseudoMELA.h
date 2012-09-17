#ifndef MELA_PseudoMela_h
#define MELA_PseudoMela_h

/** \class PseudoMela
 *
 *  PseudoMELA discriminator 
 *
 *
 *  $Date: 2012/09/14 16:29:02 $
 *  $Revision: 1.1 $
 *  \author JHU
 */

#include <TLorentzVector.h>

class AngularPdfFactory;
class RooRealVar;


class PseudoMELA{

public:

  PseudoMELA();

  ~PseudoMELA(){};

  void eval(TLorentzVector Z1_lept,
	    TLorentzVector Z1_antiLept,
	    TLorentzVector Z2_lept,
	    TLorentzVector Z2_antiLept,
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
  AngularPdfFactory *PSHiggs;


    
    

};

#endif
