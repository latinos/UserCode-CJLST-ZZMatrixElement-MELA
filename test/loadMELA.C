//
// MELA root package loader - see testKD.C for instructions
//
{
  gSystem->AddIncludePath("-I$ROOFITSYS/include/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/");
  gSystem->Load("libZZMatrixElementMELA.so");
  gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/MELA/interface/Mela.h+");
  gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/MELA/interface/TVar.hh+");
  gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/MELA/interface/newMELA.h+");

  gSystem->AddIncludePath("-I$CMSSW_BASE/src/ ");  
}
