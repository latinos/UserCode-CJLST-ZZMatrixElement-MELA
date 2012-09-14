#ifndef computeAngles_h
#define computeAngles_h

/*
 *  MELA - cf. http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/sbologne/MELAproject/
 *
 *  This class is adapted from:
 *  UserCode/scasasso/HZZ4lAnalysis/HZZ4lCommon/interface/HiggsCandidateFactory.h tag V00-00-00
 *  $Date: 2012/08/31 11:42:44 $
 *  $Revision: 1.3 $
 */

#include "TLorentzVector.h"

namespace pat{
  class CompositeCandidate;
}

namespace mela {
  //  Compute decay angles for a ZZ system. Zs are taken from the CompositeCandidates, and are sorted
  //  so that Z1 = closest to nominal mZ. 
  //  Leptons are sorted so that 1 is the negative-charged one, for OS pairs.
  void computeAngles(const pat::CompositeCandidate& obj, 
		     double& costheta1, 
		     double& costheta2, 
		     double& phi, 
		     double& costhetastar, 
		     double& phistar1, 
		     double& phistar2
		     );

  //  Compute decay angles for a ZZ system. 
  //  Leptons are sorted so that M11/M21 are the negative-charged one, for OS pairs.
  void computeAngles(TLorentzVector p4M11,
		     TLorentzVector p4M12, 
		     TLorentzVector p4M21, 
		     TLorentzVector p4M22, 
		     double& costheta1, 
		     double& costheta2, 
		     double& Phi, 
		     double& costhetastar, 
		     double& Phi1);
}
#endif
