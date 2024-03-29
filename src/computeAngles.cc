/** \file
 *
 *  MELA - cf. http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/sbologne/MELAproject/
 *
 *  $Date: 2013/01/21 23:51:30 $
 *  $Revision: 1.6 $
 */

#include <ZZMatrixElement/MELA/src/computeAngles.h>
#include <TLorentzVector.h>
#include <algorithm>
#include <iostream>

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
  if (boostV2.Mag()>=1.) {
    cout << "Warning: Mela::computeAngles: Z2 boost with beta=1, scaling down" << endl;
    boostV2*=0.9999;
  }
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


  if (isnan(costhetastar) || isnan(costheta1) || isnan(costheta2) || isnan(Phi) || isnan(Phi1)){
    cout << "WARNING: NaN in computeAngles: " 
	 << costhetastar << " " 
	 << costheta1  << " " 
	 << costheta2  << " " 
	 << Phi  << " " 
	 << Phi1  << " " << endl;
    cout << "   boostV1: " <<boostV1.Pt() << " " << boostV1.Eta() << " " << boostV1.Phi() << " " << boostV1.Mag() << endl;
    cout << "   boostV2: " <<boostV2.Pt() << " " << boostV2.Eta() << " " << boostV2.Phi() << " " << boostV2.Mag() << endl;
  }
    
}

void mela::computeAnglesCS(TLorentzVector p4M11, int Z1_lept1Id,
						 TLorentzVector p4M12, int Z1_lept2Id,
						 TLorentzVector p4M21, int Z2_lept1Id,
						 TLorentzVector p4M22, int Z2_lept2Id,
						 float pbeam,  
						 float& costhetastar, 
						 float& costheta1, 
						 float& costheta2, 
						 float& Phi, 
						 float& Phi1){
	

	TVector3 LabXaxis( 1.0, 0.0, 0.0 );
	TVector3 LabYaxis( 0.0, 1.0, 0.0 );
	TVector3 LabZaxis( 0.0, 0.0, 1.0 );
	
	float Mprot = 0.938;
	float Ebeam = sqrt(pbeam*pbeam + Mprot*Mprot);
	
	TLorentzVector targ(0.,0.,-pbeam,Ebeam);
	TLorentzVector beam(0.,0., pbeam,Ebeam);
	
	//build Z 4-vectors
	TLorentzVector p4Z1 = p4M11 + p4M12;
	TLorentzVector p4Z2 = p4M21 + p4M22;
	
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
	TVector3 boostX = -(p4H.BoostVector());

	/////////////////////////////
	// Collin-Sopper calculation:
	// in the CS frame, the z-axis is along the bisectrice of one beam and the opposite of the other beam,
	// after their boost in X
	///////////////////////////////
	// Rotation for the CS Frame
	
	TRotation rotationCS;
	
	TLorentzVector beaminX(beam);
	TLorentzVector targinX(targ);
	targinX.Boost( boostX );
	beaminX.Boost( boostX );
	
	//Bisectrice: sum of unit vectors (remember: you need to invert one beam vector)
	TVector3 beam_targ_bisecinX((beaminX.Vect().Unit() - targinX.Vect().Unit()).Unit());
	
	// Define a rotationCS Matrix, with Z along the bisectric, 
	TVector3 newZaxisCS(beam_targ_bisecinX.Unit());
	TVector3 newYaxisCS(beaminX.Vect().Unit().Cross(newZaxisCS).Unit());
	TVector3 newXaxisCS(newYaxisCS.Unit().Cross(newZaxisCS).Unit());
	rotationCS.RotateAxes(newXaxisCS, newYaxisCS, newZaxisCS);
	rotationCS.Invert();
	
	//// costhetastar
	TLorentzVector thep4Z1inXFrame_rotCS( p4Z1 );
	TLorentzVector thep4Z2inXFrame_rotCS( p4Z2 );
	thep4Z1inXFrame_rotCS.Transform(rotationCS);
	thep4Z2inXFrame_rotCS.Transform(rotationCS);
	thep4Z1inXFrame_rotCS.Boost( boostX );
	thep4Z2inXFrame_rotCS.Boost( boostX );
	TVector3 theZ1XrotCS_p3 = TVector3( thep4Z1inXFrame_rotCS.X(), thep4Z1inXFrame_rotCS.Y(), thep4Z1inXFrame_rotCS.Z() );
	TVector3 theZ2XrotCS_p3 = TVector3( thep4Z2inXFrame_rotCS.X(), thep4Z2inXFrame_rotCS.Y(), thep4Z2inXFrame_rotCS.Z() );    
	costhetastar = theZ1XrotCS_p3.CosTheta();
	
	//// --------------------------- Phi and Phi1 (old phistar1 - azimuthal production angle)
	//    TVector3 boostX = -(thep4H.BoostVector());
	TLorentzVector p4M11_BX_rotCS( p4M11 );
	TLorentzVector p4M12_BX_rotCS( p4M12 );
	TLorentzVector p4M21_BX_rotCS( p4M21 );
	TLorentzVector p4M22_BX_rotCS( p4M22 );
    p4M11_BX_rotCS.Transform(rotationCS);
    p4M12_BX_rotCS.Transform(rotationCS);
    p4M21_BX_rotCS.Transform(rotationCS);
    p4M22_BX_rotCS.Transform(rotationCS);
	p4M11_BX_rotCS.Boost( boostX );
	p4M12_BX_rotCS.Boost( boostX );
	p4M21_BX_rotCS.Boost( boostX );
	p4M22_BX_rotCS.Boost( boostX );
    
	TVector3 tmp1 = p4M11_BX_rotCS.Vect().Cross( p4M12_BX_rotCS.Vect() );
	TVector3 tmp2 = p4M21_BX_rotCS.Vect().Cross( p4M22_BX_rotCS.Vect() );    
    
	TVector3 normal1_BX_rotCS( tmp1.X()/tmp1.Mag(), tmp1.Y()/tmp1.Mag(), tmp1.Z()/tmp1.Mag() ); 
	TVector3 normal2_BX_rotCS( tmp2.X()/tmp2.Mag(), tmp2.Y()/tmp2.Mag(), tmp2.Z()/tmp2.Mag() ); 
	
	//// Phi
	TLorentzVector p4Z1_BX_rotCS = p4M11_BX_rotCS + p4M12_BX_rotCS;    
	float tmpSgnPhi = p4Z1_BX_rotCS.Vect().Dot( normal1_BX_rotCS.Cross( normal2_BX_rotCS) );
	float sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
	Phi = sgnPhi * acos( -1.*normal1_BX_rotCS.Dot( normal2_BX_rotCS) );
    
    //////////////
    
	TVector3 beamAxis(0,0,1);
	TVector3 tmp3 = (p4M11_BX_rotCS + p4M12_BX_rotCS).Vect();
    
	TVector3 p3V1_BX_rotCS( tmp3.X()/tmp3.Mag(), tmp3.Y()/tmp3.Mag(), tmp3.Z()/tmp3.Mag() );
	TVector3 tmp4 = beamAxis.Cross( p3V1_BX_rotCS );
	TVector3 normalSC_BX_rotCS( tmp4.X()/tmp4.Mag(), tmp4.Y()/tmp4.Mag(), tmp4.Z()/tmp4.Mag() );
	
	//// Phi1
	float tmpSgnPhi1 = p4Z1_BX_rotCS.Vect().Dot( normal1_BX_rotCS.Cross( normalSC_BX_rotCS) );
	float sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);    
	Phi1 = sgnPhi1 * acos( normal1_BX_rotCS.Dot( normalSC_BX_rotCS) );  
	
	//// --------------------------- costheta1
	//define Z1 rotation 
	TRotation rotationZ1;
	TVector3 newZaxisZ1(thep4Z1inXFrame_rotCS.Vect().Unit());
	TVector3 newXaxisZ1(newYaxisCS.Cross(newZaxisZ1).Unit() );
	TVector3 newYaxisZ1(newZaxisZ1.Cross(newXaxisZ1).Unit() );
	rotationZ1.RotateAxes(newXaxisZ1, newYaxisZ1, newZaxisZ1);
	rotationZ1.Invert();
	
	TLorentzVector thep4Z1inXFrame_rotCS_rotZ1(thep4Z1inXFrame_rotCS);
	thep4Z1inXFrame_rotCS_rotZ1.Transform(rotationZ1); 
	TVector3 boostZ1inX_rotCS_rotZ1= -(thep4Z1inXFrame_rotCS_rotZ1.BoostVector());
	
	TLorentzVector p4M11_BX_rotCS_rotZ1(p4M11_BX_rotCS);
	TLorentzVector p4M12_BX_rotCS_rotZ1(p4M12_BX_rotCS);
	TLorentzVector p4M21_BX_rotCS_rotZ1(p4M21_BX_rotCS);
	TLorentzVector p4M22_BX_rotCS_rotZ1(p4M22_BX_rotCS);
	p4M11_BX_rotCS_rotZ1.Transform(rotationZ1); 
	p4M12_BX_rotCS_rotZ1.Transform(rotationZ1); 
	p4M21_BX_rotCS_rotZ1.Transform(rotationZ1); 
	p4M22_BX_rotCS_rotZ1.Transform(rotationZ1); 
	p4M11_BX_rotCS_rotZ1.Boost(boostZ1inX_rotCS_rotZ1);
	p4M12_BX_rotCS_rotZ1.Boost(boostZ1inX_rotCS_rotZ1);
	p4M21_BX_rotCS_rotZ1.Boost(boostZ1inX_rotCS_rotZ1);
	p4M22_BX_rotCS_rotZ1.Boost(boostZ1inX_rotCS_rotZ1);
    
	TLorentzVector p4V2_BX_rotCS_rotZ1 = p4M21_BX_rotCS_rotZ1 + p4M22_BX_rotCS_rotZ1;
	//// costheta1
	costheta1 = -p4V2_BX_rotCS_rotZ1.Vect().Dot( p4M11_BX_rotCS_rotZ1.Vect() )/p4V2_BX_rotCS_rotZ1.Vect().Mag()/p4M11_BX_rotCS_rotZ1.Vect().Mag();
	
	//// --------------------------- costheta2
	//define Z2 rotation 
	TRotation rotationZ2;
	TVector3 newZaxisZ2(thep4Z2inXFrame_rotCS.Vect().Unit());
	TVector3 newXaxisZ2(newYaxisCS.Cross(newZaxisZ2).Unit() );
	TVector3 newYaxisZ2(newZaxisZ2.Cross(newXaxisZ2).Unit() );
	rotationZ2.RotateAxes(newXaxisZ2, newYaxisZ2, newZaxisZ2);
	rotationZ2.Invert();
	
	TLorentzVector thep4Z2inXFrame_rotCS_rotZ2(thep4Z2inXFrame_rotCS);
	thep4Z2inXFrame_rotCS_rotZ2.Transform(rotationZ2); 
	TVector3 boostZ2inX_rotCS_rotZ2= -(thep4Z2inXFrame_rotCS_rotZ2.BoostVector());
	
	TLorentzVector p4M11_BX_rotCS_rotZ2(p4M11_BX_rotCS);
	TLorentzVector p4M12_BX_rotCS_rotZ2(p4M12_BX_rotCS);
	TLorentzVector p4M21_BX_rotCS_rotZ2(p4M21_BX_rotCS);
	TLorentzVector p4M22_BX_rotCS_rotZ2(p4M22_BX_rotCS);
	p4M11_BX_rotCS_rotZ2.Transform(rotationZ2); 
	p4M12_BX_rotCS_rotZ2.Transform(rotationZ2); 
	p4M21_BX_rotCS_rotZ2.Transform(rotationZ2); 
	p4M22_BX_rotCS_rotZ2.Transform(rotationZ2); 
	p4M11_BX_rotCS_rotZ2.Boost(boostZ2inX_rotCS_rotZ2);
	p4M12_BX_rotCS_rotZ2.Boost(boostZ2inX_rotCS_rotZ2);
	p4M21_BX_rotCS_rotZ2.Boost(boostZ2inX_rotCS_rotZ2);
	p4M22_BX_rotCS_rotZ2.Boost(boostZ2inX_rotCS_rotZ2);
	

	TLorentzVector p4V1_BX_rotCS_rotZ2= p4M11_BX_rotCS_rotZ2 + p4M12_BX_rotCS_rotZ2;
	//// costheta2
	costheta2 = -p4V1_BX_rotCS_rotZ2.Vect().Dot( p4M21_BX_rotCS_rotZ2.Vect() )/p4V1_BX_rotCS_rotZ2.Vect().Mag()/p4M21_BX_rotCS_rotZ2.Vect().Mag();


  if (isnan(costhetastar) || isnan(costheta1) || isnan(costheta2) || isnan(Phi) || isnan(Phi1)){
    cout << "WARNING: NaN in computeAngles: " 
	 << costhetastar << " " 
	 << costheta1  << " " 
	 << costheta2  << " " 
	 << Phi  << " " 
	 << Phi1  << " " << endl;
  }
}


