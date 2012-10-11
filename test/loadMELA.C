//
// MELA root package loader - see testKD.C for instructions
//
{
  gSystem->Load("libZZMatrixElementMELA.so");
  gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/MELA/interface/Mela.h+");
 gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/MELA/interface/SuperMELA.h+");
 gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/MELA/interface/PseudoMELA.h+");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/ ");  
}
