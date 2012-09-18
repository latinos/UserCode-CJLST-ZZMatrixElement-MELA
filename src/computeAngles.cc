/** \file
 *
 *  MELA - cf. http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/sbologne/MELAproject/
 *
 *  $Date: 2012/09/17 19:06:31 $
 *  $Revision: 1.2 $
 */

#include <ZZMatrixElement/MELA/src/computeAngles.h>
#include <TLorentzVector.h>
#include <algorithm>

const float PDGZmass = 91.1876;

using namespace std;

void mela::computeAngles(TLorentzVector p4M11, int Z1_lept1Id,
			 TLorentzVector p4M12, int Z1_lept2Id,
			 TLorentzVector p4M21, int Z2_lept1Id,
			 TLorentzVector p4M22, int Z2_lept2Id,
			 float& costhetastar, 
			 float& costheta1, 
			 float& costheta2, 
			 float& Phi, 
			 float& Phi1){

  //build Z 4-vectors
  TLorentzVector p4Z1 = p4M11 + p4M12;
  TLorentzVector p4Z2 = p4M21 + p4M22;

  // Sort Z1 and Z2 so that Z1 is the closest-to-mZ
  if (fabs(PDGZmass-p4Z1.M()) > fabs(PDGZmass-p4Z2.M()) ){
    swap(p4Z1, p4Z2);
    swap(p4M11, p4M21);
    swap(p4M12, p4M22);
  }

  // Sort Z1 leptons so that:
  if ( (Z1_lept1Id*Z1_lept2Id<0 && Z1_lept1Id<0) || // for OS pairs: lep1 must be the negative one
       (Z1_lept1Id*Z1_lept2Id>0 && p4M11.Phi()<=p4M12.Phi()) //for SS pairs: use random deterministic convention
       ) {
    swap(p4M11, p4M12);
  }
  
  // Same for Z2 leptons
  if ( (Z2_lept1Id*Z2_lept2Id<0 && Z2_lept1Id<0) ||
       (Z2_lept1Id*Z2_lept2Id>0 && p4M21.Phi()<=p4M22.Phi()) 
       ) {
    swap(p4M21, p4M22);
  }

  
  //build H 4-vectors
  TLorentzVector p4H = p4Z1 + p4Z2; 

  // -----------------------------------

  //// costhetastar
  TVector3 boostX = -(p4H.BoostVector());
  TLorentzVector thep4Z1inXFrame( p4Z1 );
  TLorentzVector thep4Z2inXFrame( p4Z2 );
  thep4Z1inXFrame.Boost( boostX );
  thep4Z2inXFrame.Boost( boostX );
  TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
  TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );    
  costhetastar = theZ1X_p3.CosTheta();

  //// --------------------------- costheta1
  TVector3 boostV1 = -(p4Z1.BoostVector());
  TLorentzVector p4M11_BV1( p4M11 );
  TLorentzVector p4M12_BV1( p4M12 );
  TLorentzVector p4M21_BV1( p4M21 );
  TLorentzVector p4M22_BV1( p4M22 );
  p4M11_BV1.Boost( boostV1 );
  p4M12_BV1.Boost( boostV1 );
  p4M21_BV1.Boost( boostV1 );
  p4M22_BV1.Boost( boostV1 );
    
  TLorentzVector p4V2_BV1 = p4M21_BV1 + p4M22_BV1;
  //// costheta1
  costheta1 = -p4V2_BV1.Vect().Dot( p4M11_BV1.Vect() )/p4V2_BV1.Vect().Mag()/p4M11_BV1.Vect().Mag();

  //// --------------------------- costheta2
  TVector3 boostV2 = -(p4Z2.BoostVector());
  TLorentzVector p4M11_BV2( p4M11 );
  TLorentzVector p4M12_BV2( p4M12 );
  TLorentzVector p4M21_BV2( p4M21 );
  TLorentzVector p4M22_BV2( p4M22 );
  p4M11_BV2.Boost( boostV2 );
  p4M12_BV2.Boost( boostV2 );
  p4M21_BV2.Boost( boostV2 );
  p4M22_BV2.Boost( boostV2 );
    
  TLorentzVector p4V1_BV2 = p4M11_BV2 + p4M12_BV2;
  //// costheta2
  costheta2 = -p4V1_BV2.Vect().Dot( p4M21_BV2.Vect() )/p4V1_BV2.Vect().Mag()/p4M21_BV2.Vect().Mag();
    
  //// --------------------------- Phi and Phi1 (old phistar1 - azimuthal production angle)
  //    TVector3 boostX = -(thep4H.BoostVector());
  TLorentzVector p4M11_BX( p4M11 );
  TLorentzVector p4M12_BX( p4M12 );
  TLorentzVector p4M21_BX( p4M21 );
  TLorentzVector p4M22_BX( p4M22 );
    
  p4M11_BX.Boost( boostX );
  p4M12_BX.Boost( boostX );
  p4M21_BX.Boost( boostX );
  p4M22_BX.Boost( boostX );
    
  TVector3 tmp1 = p4M11_BX.Vect().Cross( p4M12_BX.Vect() );
  TVector3 tmp2 = p4M21_BX.Vect().Cross( p4M22_BX.Vect() );    
    
  TVector3 normal1_BX( tmp1.X()/tmp1.Mag(), tmp1.Y()/tmp1.Mag(), tmp1.Z()/tmp1.Mag() ); 
  TVector3 normal2_BX( tmp2.X()/tmp2.Mag(), tmp2.Y()/tmp2.Mag(), tmp2.Z()/tmp2.Mag() ); 

  //// Phi
  TLorentzVector p4Z1_BX = p4M11_BX + p4M12_BX;    
  float tmpSgnPhi = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normal2_BX) );
  float sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
  Phi = sgnPhi * acos( -1.*normal1_BX.Dot( normal2_BX) );
    
    
  //////////////
    
  TVector3 beamAxis(0,0,1);
  TVector3 tmp3 = (p4M11_BX + p4M12_BX).Vect();
    
  TVector3 p3V1_BX( tmp3.X()/tmp3.Mag(), tmp3.Y()/tmp3.Mag(), tmp3.Z()/tmp3.Mag() );
  TVector3 tmp4 = beamAxis.Cross( p3V1_BX );
  TVector3 normalSC_BX( tmp4.X()/tmp4.Mag(), tmp4.Y()/tmp4.Mag(), tmp4.Z()/tmp4.Mag() );
        
  //// Phi1
  float tmpSgnPhi1 = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normalSC_BX) );
  float sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);    
  Phi1 = sgnPhi1 * acos( normal1_BX.Dot( normalSC_BX) );    
    
}

