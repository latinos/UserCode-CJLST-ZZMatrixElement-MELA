#ifndef computeAngles_h
#define computeAngles_h

/*
 *  MELA - cf. http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/sbologne/MELAproject/
 *
 *  $Date: 2012/09/17 19:06:32 $
 *  $Revision: 1.2 $
 */

#include "TLorentzVector.h"

namespace mela {
  // Compute decay angles from the lepton four-vectors and pdgIds.  
  /// Z1_lept1 and  Z1_lept2 are supposed to come from the same Z.
  /// Zs and leptons are re-ordered internally according to a defined convention:
  /// Z1 = closest to nominal mZ; lept1 = negative-charged lepton (for OS pairs). 
  void computeAngles(TLorentzVector Z1_lept1, int Z1_lept1Id,
		     TLorentzVector Z1_lept2, int Z1_lept2Id,
		     TLorentzVector Z2_lept1, int Z2_lept1Id,
		     TLorentzVector Z2_lept2, int Z2_lept2Id,
		     float& costhetastar, 
		     float& costheta1, 
		     float& costheta2, 
		     float& Phi, 
		     float& Phi1);
}
#endif
