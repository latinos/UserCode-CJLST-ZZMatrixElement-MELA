/** \file
 *
 *  MELA - cf. http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/sbologne/MELAproject/
 *  This class is adapted from:
 *  UserCode/scasasso/HZZ4lAnalysis/HZZ4lCommon/interface/HiggsCandidateFactory.h tag V00-00-00
 *
 * With the following modifications:
 * - input type changed to pat::CompositeCandidate
 * - leg1() -> daughter(0);
 * - leg2() -> daughter(1);
 * - Angles returned by reference.
 *
 *  $Date: 2012/08/31 11:47:41 $
 *  $Revision: 1.9 $
 */

#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <ZZMatrixElement/MELA/src/computeAngles.h>
#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <TLorentzVector.h>

#include <iostream>
using namespace std;

  //  Compute decay angles for a ZZ system. 
  //  Leptons are sorted so that M11/M21 are the negative-charged one, for OS pairs.
void mela::computeAngles(TLorentzVector p4M11,
			 TLorentzVector p4M12, 
			 TLorentzVector p4M21, 
			 TLorentzVector p4M22, 
			 double& costheta1, 
			 double& costheta2, 
			 double& Phi, 
			 double& costhetastar, 
			 double& Phi1){

  //build Z 4-vectors
  TLorentzVector p4Z1 = p4M11 + p4M12;
  TLorentzVector p4Z2 = p4M21 + p4M22;
  
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
  double tmpSgnPhi = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normal2_BX) );
  double sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
  Phi = sgnPhi * acos( -1.*normal1_BX.Dot( normal2_BX) );
    
    
  //////////////
    
  TVector3 beamAxis(0,0,1);
  TVector3 tmp3 = (p4M11_BX + p4M12_BX).Vect();
    
  TVector3 p3V1_BX( tmp3.X()/tmp3.Mag(), tmp3.Y()/tmp3.Mag(), tmp3.Z()/tmp3.Mag() );
  TVector3 tmp4 = beamAxis.Cross( p3V1_BX );
  TVector3 normalSC_BX( tmp4.X()/tmp4.Mag(), tmp4.Y()/tmp4.Mag(), tmp4.Z()/tmp4.Mag() );
        
  //// Phi1
  double tmpSgnPhi1 = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normalSC_BX) );
  double sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);    
  Phi1 = sgnPhi1 * acos( normal1_BX.Dot( normalSC_BX) );    
    
}

void mela::computeAngles(const pat::CompositeCandidate& obj,
			 double& costheta1, 
			 double& costheta2, 
			 double& phi, 
			 double& costhetastar, 
			 double& phistar1, 
			 double& phistar2
			 ){

  // These are irrelevant angles
  double phi1; 
  double phi2;

  //will convert to TLorentzVector
  math::XYZTLorentzVector leg1 = obj.daughter(0)->p4();
  math::XYZTLorentzVector leg2 = obj.daughter(1)->p4();
  math::XYZTLorentzVector leg11 = obj.daughter(0)->daughter(0)->p4();
  math::XYZTLorentzVector leg12 = obj.daughter(0)->daughter(1)->p4();
  math::XYZTLorentzVector leg21 = obj.daughter(1)->daughter(0)->p4();
  math::XYZTLorentzVector leg22 = obj.daughter(1)->daughter(1)->p4();


  // Check if there's FSR for either of the Zs
  int d0FSR = obj.userFloat("d0.dauWithFSR");
  if (d0FSR>=0) {
    if (obj.daughter(0)->numberOfDaughters()!=3) cout << "ERROR: mela::computeAngles: problem in FSR" << endl;
    //    leg1 = leg1 + obj.daughter(0)->daughter(2)->p4();
    if (d0FSR==0) {
      leg11 = leg11 + obj.daughter(0)->daughter(2)->p4();
    } else if (d0FSR==1){
      leg12 = leg12 + obj.daughter(0)->daughter(2)->p4();      
    }
  }
  
  int d1FSR = obj.userFloat("d1.dauWithFSR");
  if (d1FSR>=0) {
    if (obj.daughter(1)->numberOfDaughters()!=3) cout << "ERROR: mela::computeAngles: problem in FSR" << endl;
    //    leg2 = leg2 + obj.daughter(1)->daughter(2)->p4();
    if (d1FSR==0) {
      leg21 = leg21 + obj.daughter(1)->daughter(2)->p4();
    } else if (d1FSR==1){
      leg22 = leg22 + obj.daughter(1)->daughter(2)->p4();      
    }
  }


  TLorentzVector Z1(leg1.x(),leg1.y(),leg1.z(),leg1.t());
  TLorentzVector Z2(leg2.x(),leg2.y(),leg2.z(),leg2.t()) ;
  TLorentzVector thep4H = Z1+Z2;
  
  //define Z1 as the one nearest to nominal Z mass
  //  const double PDGZmass = 91.2; // Use consistent value with the rest of the code
  const float PDGZmass = 91.1876;

  TLorentzVector thep4Z1; TLorentzVector thep4M11; TLorentzVector thep4M12;
  TLorentzVector thep4Z2; TLorentzVector thep4M21; TLorentzVector thep4M22;
  
  if ( fabs(PDGZmass-Z1.M()) > fabs(PDGZmass-Z2.M()) ){	  			
    thep4Z1 = Z2; 
    //if the 2 leptons have not opposite charge (control-region) randomize the order using azimuthal angle in lab RF
    if((obj.daughter(1)->daughter(0)->charge())*(obj.daughter(1)->daughter(1)->charge()) <0){ //if opposite charge
      if(obj.daughter(1)->daughter(0)->charge()<0){
	thep4M11 = TLorentzVector(leg21.x(),leg21.y(),leg21.z(),leg21.t());// obj->leg2().leg1().p4(); 
	thep4M12 = TLorentzVector(leg22.x(),leg22.y(),leg22.z(),leg22.t());//obj->leg2().leg2().p4();
      }
      else{
	thep4M11 = TLorentzVector(leg22.x(),leg22.y(),leg22.z(),leg22.t());//obj->leg2().leg2().p4();
	thep4M12 = TLorentzVector(leg21.x(),leg21.y(),leg21.z(),leg21.t());// obj->leg2().leg1().p4(); 
      }
    }  //end if opposite charge
    else{
      if(obj.daughter(1)->daughter(0)->phi()>obj.daughter(1)->daughter(1)->phi()){  //if not opposite charge
	thep4M11 = TLorentzVector(leg21.x(),leg21.y(),leg21.z(),leg21.t());// obj->leg2().leg1().p4(); 
	thep4M12 = TLorentzVector(leg22.x(),leg22.y(),leg22.z(),leg22.t());//obj->leg2().leg2().p4();
      }
      else{
	thep4M11 = TLorentzVector(leg22.x(),leg22.y(),leg22.z(),leg22.t());//obj->leg2().leg2().p4();
	thep4M12 = TLorentzVector(leg21.x(),leg21.y(),leg21.z(),leg21.t());// obj->leg2().leg1().p4(); 
      }
    }  //end if not opposite charge
      
    thep4Z2 = Z1; 
    if((obj.daughter(0)->daughter(0)->charge())*(obj.daughter(0)->daughter(1)->charge()) <0){ //if opposite charge
      if(obj.daughter(0)->daughter(0)->charge()<0){
	thep4M21 = TLorentzVector(leg11.x(),leg11.y(),leg11.z(),leg11.t());//obj->leg1().leg1().p4(); 
	thep4M22 = TLorentzVector(leg12.x(),leg12.y(),leg12.z(),leg12.t());//obj->leg1().leg2().p4();
      }
      else{
	thep4M21 = TLorentzVector(leg12.x(),leg12.y(),leg12.z(),leg12.t());//obj->leg1().leg2().p4(); 
	thep4M22 = TLorentzVector(leg11.x(),leg11.y(),leg11.z(),leg11.t());//obj->leg1().leg1().p4();
      }
    } //end if opposite charge
    else{  
      if(obj.daughter(0)->daughter(0)->phi()>obj.daughter(0)->daughter(1)->phi()){
	thep4M21 = TLorentzVector(leg11.x(),leg11.y(),leg11.z(),leg11.t());//obj->leg1().leg1().p4(); 
	thep4M22 = TLorentzVector(leg12.x(),leg12.y(),leg12.z(),leg12.t());//obj->leg1().leg2().p4();
      }
      else{
	thep4M21 = TLorentzVector(leg12.x(),leg12.y(),leg12.z(),leg12.t());//obj->leg1().leg2().p4(); 
	thep4M22 = TLorentzVector(leg11.x(),leg11.y(),leg11.z(),leg11.t());//obj->leg1().leg1().p4();
      }
    } //end if not opposite charge
  } // end if fabs(PDGZmass-Z1.M()) > fabs(PDGZmass-Z2.M()) 
  else {
    thep4Z1 = Z1; 
    if((obj.daughter(0)->daughter(0)->charge())*(obj.daughter(0)->daughter(1)->charge()) <0){ //if opposite charge
      if(obj.daughter(0)->daughter(0)->charge()<0){
	thep4M11 =TLorentzVector(leg11.x(),leg11.y(),leg11.z(),leg11.t());// obj->leg1().leg1().p4(); 
	thep4M12 = TLorentzVector(leg12.x(),leg12.y(),leg12.z(),leg12.t());//obj->leg1().leg2().p4();
      }
      else{
	thep4M11 = TLorentzVector(leg12.x(),leg12.y(),leg12.z(),leg12.t());//obj->leg1().leg2().p4(); 
	thep4M12 = TLorentzVector(leg11.x(),leg11.y(),leg11.z(),leg11.t());//obj->leg1().leg1().p4();
      }
    } // end if opposite charge
    else{
     if(obj.daughter(0)->daughter(0)->phi()>obj.daughter(0)->daughter(1)->phi()){
	thep4M11 =TLorentzVector(leg11.x(),leg11.y(),leg11.z(),leg11.t());// obj->leg1().leg1().p4(); 
	thep4M12 = TLorentzVector(leg12.x(),leg12.y(),leg12.z(),leg12.t());//obj->leg1().leg2().p4();
      }
      else{
	thep4M11 = TLorentzVector(leg12.x(),leg12.y(),leg12.z(),leg12.t());//obj->leg1().leg2().p4(); 
	thep4M12 = TLorentzVector(leg11.x(),leg11.y(),leg11.z(),leg11.t());//obj->leg1().leg1().p4();
      }     
    } //end if not opposite charge
   thep4Z2 = Z2; 
   if((obj.daughter(1)->daughter(0)->charge())*(obj.daughter(1)->daughter(1)->charge()) <0){ //if opposite charge
     if(obj.daughter(1)->daughter(0)->charge()<0){
       thep4M21 = TLorentzVector(leg21.x(),leg21.y(),leg21.z(),leg21.t());// obj->leg2().leg1().p4(); 
       thep4M22 = TLorentzVector(leg22.x(),leg22.y(),leg22.z(),leg22.t());//obj->leg2().leg2().p4();
     }
     else{
       thep4M21 = TLorentzVector(leg22.x(),leg22.y(),leg22.z(),leg22.t());//obj->leg2().leg2().p4(); 
       thep4M22 = TLorentzVector(leg21.x(),leg21.y(),leg21.z(),leg21.t());//obj->leg2().leg1().p4();
     }
   } //end if opposite charge
   else{
     if(obj.daughter(1)->daughter(0)->phi()>obj.daughter(1)->daughter(1)->phi()){
       thep4M21 = TLorentzVector(leg21.x(),leg21.y(),leg21.z(),leg21.t());// obj->leg2().leg1().p4(); 
       thep4M22 = TLorentzVector(leg22.x(),leg22.y(),leg22.z(),leg22.t());//obj->leg2().leg2().p4();
     }
     else{
       thep4M21 = TLorentzVector(leg22.x(),leg22.y(),leg22.z(),leg22.t());//obj->leg2().leg2().p4(); 
       thep4M22 = TLorentzVector(leg21.x(),leg21.y(),leg21.z(),leg21.t());//obj->leg2().leg1().p4();
     }
   }//end if not opposite charge
  } // end if not fabs(PDGZmass-Z1.M()) > fabs(PDGZmass-Z2.M()) 
  double norm;
//   double costheta1=-99;
//   double costheta2=-99;
//   double phi=-99; 
//   double costhetastar=-99; 
//   double phistar1=-99 ;
//   double phistar2 =-99;
  double phistar12 = -99; 
//   double phi1=-99; 
//   double phi2=-99;
  
  TVector3 boostX = -(thep4H.BoostVector());
  TLorentzVector thep4Z1inXFrame( thep4Z1 );
  TLorentzVector thep4Z2inXFrame( thep4Z2 );	
  thep4Z1inXFrame.Boost( boostX );
  thep4Z2inXFrame.Boost( boostX );
  TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
  TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );
  
  // calculate phi1, phi2, costhetastar
  phi1 = theZ1X_p3.Phi();
  phi2 = theZ2X_p3.Phi();
  
  ///////////////////////////////////////////////
  // check for z1/z2 convention, redefine all 4 vectors with convention
  ///////////////////////////////////////////////	
  TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;
  p4H = thep4H;
  
  /* ORDER OF Z1 AND Z2 ALREADY CHOSEN IN MAIN FUNCTION!!!!!! - - - - - - 
     if ((phi1 < 0)&&(phi1 >= -TMath::Pi())){   // old convention based on phi
     p4Z1 = thep4Z2; p4M11 = thep4M21; p4M12 = thep4M22;
     p4Z2 = thep4Z1; p4M21 = thep4M11; p4M22 = thep4M12;		
     costhetastar = theZ2X_p3.CosTheta();
     }
     else{
     p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
     p4Z2 = thep4Z2; p4M21 = thep4M21; p4M22 = thep4M22;
     costhetastar = theZ1X_p3.CosTheta();
     }
     - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - -*/
  
  p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
  p4Z2 = thep4Z2; p4M21 = thep4M21; p4M22 = thep4M22;
  costhetastar = theZ1X_p3.CosTheta();
  
  // now helicity angles................................
  // ...................................................
  TVector3 boostZ1 = -(p4Z1.BoostVector());
  TLorentzVector p4Z2Z1(p4Z2);
  p4Z2Z1.Boost(boostZ1);
  //find the decay axis
  /////TVector3 unitx_1 = -Hep3Vector(p4Z2Z1);
  TVector3 unitx_1( -p4Z2Z1.X(), -p4Z2Z1.Y(), -p4Z2Z1.Z() );
  norm = 1/(unitx_1.Mag());
  unitx_1*=norm;
  //boost daughters of z2
  TLorentzVector p4M21Z1(p4M21);
  TLorentzVector p4M22Z1(p4M22);
  p4M21Z1.Boost(boostZ1);
  p4M22Z1.Boost(boostZ1);
  //create z and y axes
  /////TVector3 unitz_1 = Hep3Vector(p4M21Z1).cross(Hep3Vector(p4M22Z1));
  TVector3 p4M21Z1_p3( p4M21Z1.X(), p4M21Z1.Y(), p4M21Z1.Z() );
  TVector3 p4M22Z1_p3( p4M22Z1.X(), p4M22Z1.Y(), p4M22Z1.Z() );
  TVector3 unitz_1 = p4M21Z1_p3.Cross( p4M22Z1_p3 );
  norm = 1/(unitz_1.Mag());
  unitz_1 *= norm;
  TVector3 unity_1 = unitz_1.Cross(unitx_1);
	
  //caculate theta1
  TLorentzVector p4M11Z1(p4M11);
  p4M11Z1.Boost(boostZ1);
  TVector3 p3M11( p4M11Z1.X(), p4M11Z1.Y(), p4M11Z1.Z() );
  TVector3 unitM11 = p3M11.Unit();
  double x_m11 = unitM11.Dot(unitx_1); double y_m11 = unitM11.Dot(unity_1); double z_m11 = unitM11.Dot(unitz_1);
  TVector3 M11_Z1frame(y_m11, z_m11, x_m11);
  costheta1 = M11_Z1frame.CosTheta();
  //std::cout << "theta1: " << M11_Z1frame.Theta() << std::endl;
  //////-----------------------old way of calculating phi---------------/////////
  phi = M11_Z1frame.Phi();
	
  //set axes for other system
  TVector3 boostZ2 = -(p4Z2.BoostVector());
  TLorentzVector p4Z1Z2(p4Z1);
  p4Z1Z2.Boost(boostZ2);
  TVector3 unitx_2( -p4Z1Z2.X(), -p4Z1Z2.Y(), -p4Z1Z2.Z() );
  norm = 1/(unitx_2.Mag());
  unitx_2*=norm;
  //boost daughters of z2
  TLorentzVector p4M11Z2(p4M11);
  TLorentzVector p4M12Z2(p4M12);
  p4M11Z2.Boost(boostZ2);
  p4M12Z2.Boost(boostZ2);
  TVector3 p4M11Z2_p3( p4M11Z2.X(), p4M11Z2.Y(), p4M11Z2.Z() );
  TVector3 p4M12Z2_p3( p4M12Z2.X(), p4M12Z2.Y(), p4M12Z2.Z() );
  TVector3 unitz_2 = p4M11Z2_p3.Cross( p4M12Z2_p3 );
  norm = 1/(unitz_2.Mag());
  unitz_2*=norm;
  TVector3 unity_2 = unitz_2.Cross(unitx_2);
  //calcuate theta2
  TLorentzVector p4M21Z2(p4M21);
  p4M21Z2.Boost(boostZ2);
  TVector3 p3M21( p4M21Z2.X(), p4M21Z2.Y(), p4M21Z2.Z() );
  TVector3 unitM21 = p3M21.Unit();
  double x_m21 = unitM21.Dot(unitx_2); double y_m21 = unitM21.Dot(unity_2); double z_m21 = unitM21.Dot(unitz_2);
  TVector3 M21_Z2frame(y_m21, z_m21, x_m21);
  costheta2 = M21_Z2frame.CosTheta();
	
  // calculate phi
  //calculating phi_n
  TLorentzVector n_p4Z1inXFrame( p4Z1 );
  TLorentzVector n_p4M11inXFrame( p4M11 );
  n_p4Z1inXFrame.Boost( boostX );
  n_p4M11inXFrame.Boost( boostX );        
  TVector3 n_p4Z1inXFrame_unit = n_p4Z1inXFrame.Vect().Unit();
  TVector3 n_p4M11inXFrame_unit = n_p4M11inXFrame.Vect().Unit();  
  TVector3 n_unitz_1( n_p4Z1inXFrame_unit );
  //// y-axis is defined by neg lepton cross z-axis
  //// the subtle part is here...
  //////////TVector3 n_unity_1 = n_p4M11inXFrame_unit.Cross( n_unitz_1 );
  TVector3 n_unity_1 = n_unitz_1.Cross( n_p4M11inXFrame_unit );
  TVector3 n_unitx_1 = n_unity_1.Cross( n_unitz_1 );
	
  TLorentzVector n_p4M21inXFrame( p4M21 );
  n_p4M21inXFrame.Boost( boostX );
  TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();
  //rotate into other plane
  TVector3 n_p4M21inXFrame_unitprime( n_p4M21inXFrame_unit.Dot(n_unitx_1), n_p4M21inXFrame_unit.Dot(n_unity_1), n_p4M21inXFrame_unit.Dot(n_unitz_1) );
	
  ///////-----------------new way of calculating phi-----------------///////
  //double phi_n =  n_p4M21inXFrame_unitprime.Phi();
  /// and then calculate phistar1
  TVector3 n_p4PartoninXFrame_unit( 0.0, 0.0, 1.0 );
  TVector3 n_p4PartoninXFrame_unitprime( n_p4PartoninXFrame_unit.Dot(n_unitx_1), n_p4PartoninXFrame_unit.Dot(n_unity_1), n_p4PartoninXFrame_unit.Dot(n_unitz_1) );
  // negative sign is for arrow convention in paper
  phistar1 = (n_p4PartoninXFrame_unitprime.Phi());
	
  // and the calculate phistar2
  TLorentzVector n_p4Z2inXFrame( p4Z2 );
  n_p4Z2inXFrame.Boost( boostX );
  TVector3 n_p4Z2inXFrame_unit = n_p4Z2inXFrame.Vect().Unit();
  ///////TLorentzVector n_p4M21inXFrame( p4M21 );
  //////n_p4M21inXFrame.Boost( boostX );        
  ////TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();  
  TVector3 n_unitz_2( n_p4Z2inXFrame_unit );
  //// y-axis is defined by neg lepton cross z-axis
  //// the subtle part is here...
  //////TVector3 n_unity_2 = n_p4M21inXFrame_unit.Cross( n_unitz_2 );
  TVector3 n_unity_2 = n_unitz_2.Cross( n_p4M21inXFrame_unit );
  TVector3 n_unitx_2 = n_unity_2.Cross( n_unitz_2 );
  TVector3 n_p4PartoninZ2PlaneFrame_unitprime( n_p4PartoninXFrame_unit.Dot(n_unitx_2), n_p4PartoninXFrame_unit.Dot(n_unity_2), n_p4PartoninXFrame_unit.Dot(n_unitz_2) );
  phistar2 = (n_p4PartoninZ2PlaneFrame_unitprime.Phi());
	
  double phistar12_0 = phistar1 + phistar2;
  if (phistar12_0 > TMath::Pi()) phistar12 = phistar12_0 - 2*TMath::Pi();
  else if (phistar12_0 < (-1.)*TMath::Pi()) phistar12 = phistar12_0 + 2*TMath::Pi();
  else phistar12 = phistar12_0;
	
//   obj->costhetastar_ = costhetastar;
//   obj->helphi_ = phi;
//   obj->helphiZ1_ = phi1;
//   obj->helphiZ2_ = phi2;
//   obj->helcosthetaZ1_ = costheta1;
//   obj->helcosthetaZ2_ = costheta2;
//   obj->phistarZ1_ = phistar1;
//   obj->phistarZ2_ = phistar2; 

  // Get rid of -Wunused warnings
  if (phi1||phi2||phistar12) {}

}
