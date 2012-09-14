#ifndef MELA_PseudoMela_h
#define MELA_PseudoMela_h

/** \class PseudoMela
 *
 *  PseudoMELA discriminator 
 *
 *
 *  $Date: 2012/08/19 20:46:50 $
 *  $Revision: 1.4 $
 *  \author JHU
 */


class AngularPdfFactory;
class RooRealVar;


class PseudoMELA{

public:

  PseudoMELA();

  ~PseudoMELA(){};


  float eval(float zzmass, float z1mass, float z2mass, 
	     float costhetstar, 
	     float costheta1, 
	     float costheta2, 
	     float phi, 
	     float phistar1);


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
