// #include "TProof.h"
// #include "TProofDebug.h"


void
doit(const Long64_t nmbEntries = TTree::kMaxEntries)
{
	R__LOAD_LIBRARY(libDSelector.so)

	const size_t nmbProofThreads = 10;
	gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");

	const TString treeName        = "pippippimpimpmiss__B1_T1_U1_Effic_Tree";
	// const TString fileNamePattern = "./2017_01-ver04/batch02/tree_pippippimpimpmiss__B1_T1_U1_Effic/030730/tree_pippippimpimpmiss__B1_T1_U1_Effic_030730_000.root";
	// const TString fileNamePattern = "./2017_01-ver04/batch02/tree_pippippimpimpmiss__B1_T1_U1_Effic_030730.root";
	const TString fileNamePattern = "./tree_pippippimpimpmiss__B1_T1_U1_Effic_bggen_2017_01-ver03_batch01.root";
	// const TString treeName        = "pippimpmiss__B1_T1_U1_Effic_Tree";
	// const TString fileNamePattern = "./2017_01-ver04/batch02/tree_pippimpmiss__B1_T1_U1_Effic_030730.root";
	const TString selectorName    = "./DSelector_pippippimpimpmiss.C+";

	TChain* chain = new TChain(treeName);
	chain->Add(fileNamePattern);
	// run with proof
	gEnv->SetValue("ProofLite.Sandbox", "$PWD/.proof/");
	DPROOFLiteManager::Process_Chain(chain, selectorName.Data(), nmbProofThreads);
	// run interactively
	// chain->Process(selectorName, "", nmbEntries);
}
