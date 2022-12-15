// see https://halldweb.jlab.org/wiki/index.php/DSelector#Using_DSelector.27s_with_PROOF-Lite
void
doit(
	const bool     runPROOF   = false,
	const Long64_t nmbEntries = TTree::kMaxEntries)
{
	R__LOAD_LIBRARY(libDSelector.so)

	const size_t nmbProofThreads = 5;
	gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");

	const TString selectorName    = "./DSelector_pippippimpimpmiss.C+";
	// pi+pi+pi-pi-(p)
	const TString treeName        = "pippippimpimpmiss__B1_T1_U1_Effic_Tree";
	// const TString fileNamePattern = "./2017_01-ver04/batch02/tree_pippippimpimpmiss__B1_T1_U1_Effic/030730/tree_pippippimpimpmiss__B1_T1_U1_Effic_030730_000.root";
	// const TString fileNamePattern = "./2017_01-ver04/batch02/tree_pippippimpimpmiss__B1_T1_U1_Effic_030730.root";
	// const TString fileNamePattern = "./bggen_2017_01-ver03/batch01/tree_pippippimpimpmiss__B1_T1_U1_Effic/030285/tree_pippippimpimpmiss__B1_T1_U1_Effic_030285_000.root";
	const TString fileNamePattern = "./tree_pippippimpimpmiss__B1_T1_U1_Effic_bggen_2017_01-ver03_batch01.root";
	// pi+pi-(p)
	// const TString selectorName    = "./DSelector_pippimpmiss.C+";
	// const TString treeName        = "pippimpmiss__B1_T1_U1_Effic_Tree";
	// const TString fileNamePattern = "./2017_01-ver04/batch02/tree_pippimpmiss__B1_T1_U1_Effic_030730.root";

	TChain* chain = new TChain(treeName);
	chain->Add(fileNamePattern);
	cout << "processing tree '" << treeName << "' in file(s) '" << fileNamePattern << "' using selector '" << selectorName << "'" << endl;
	if (runPROOF) {
		// run with proof
		gEnv->SetValue("ProofLite.Sandbox", "$PWD/.proof/");
		DPROOFLiteManager::Process_Chain(chain, selectorName.Data(), nmbProofThreads);
	} else {
		// run interactively
		chain->Process(selectorName, "", nmbEntries);
	}
}
