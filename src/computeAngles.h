#ifndef computeAngles_h
#define computeAngles_h

/*
 *  MELA - cf. http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/sbologne/MELAproject/
 *
 *  $Date: 2012/09/20 20:49:07 $
 *  $Revision: 1.4 $
 */

#include "TLorentzVector.h"

namespace mela {
  /// Compute decay angles from the lepton four-vectors and pdgIds.  
  /// Theta1 is the angle corresponding to Z1.
  /// Z1_lept1 and  Z1_lept2 are supposed to come from the same Z.
  /// Leptons are re-ordered internally according to a standard convention:
  /// lept1 = negative-charged lepton (for OS pairs).
  void computeAngles(TLorentzVector Z1_lept1, int Z1_lept1Id,
		     TLorentzVector Z1_lept2, int Z1_lept2Id,
		     TLorentzVector Z2_lept1, int Z2_lept1Id,
		     TLorentzVector Z2_lept2, int Z2_lept2Id,
		     float& costhetastar, 
		     float& costheta1, 
		     float& costheta2, 
		     float& Phi, 
		     float& Phi1);

  void computeAnglesCS(TLorentzVector Z1_lept1, int Z1_lept1Id,
			 TLorentzVector Z1_lept2, int Z1_lept2Id,
			 TLorentzVector Z2_lept1, int Z2_lept1Id,
			 TLorentzVector Z2_lept2, int Z2_lept2Id,
			 float pbeam,  
			 float& costhetastar, 
			 float& costheta1, 
			 float& costheta2, 
			 float& Phi, 
			 float& Phi1);	
}
#endif
