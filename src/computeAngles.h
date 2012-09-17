#ifndef computeAngles_h
#define computeAngles_h

/*
 *  MELA - cf. http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/sbologne/MELAproject/
 *
 *  This class is adapted from:
 *  UserCode/scasasso/HZZ4lAnalysis/HZZ4lCommon/interface/HiggsCandidateFactory.h tag V00-00-00
 *  $Date: 2012/09/14 20:04:45 $
 *  $Revision: 1.1 $
 */

#include "TLorentzVector.h"

namespace pat{
  class CompositeCandidate;
}

namespace mela {

  //  Compute decay angles for a ZZ system. 
  //  Leptons are sorted so that M11/M21 are the negative-charged one, for OS pairs.
  void computeAngles(TLorentzVector p4M11,
		     TLorentzVector p4M12, 
		     TLorentzVector p4M21, 
		     TLorentzVector p4M22, 
		     float& costheta1, 
		     float& costheta2, 
		     float& Phi, 
		     float& costhetastar, 
		     float& Phi1);
}
#endif
