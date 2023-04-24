// see https://halldweb.jlab.org/wiki/index.php/DSelector#Using_DSelector.27s_with_PROOF-Lite
//TODO convert to Python; rename output files
void
runSelector(
	const bool     runPROOF   = true,
	const Long64_t nmbEntries = TTree::kMaxEntries)
{
	R__LOAD_LIBRARY(libDSelector.so)

	const size_t nmbProofThreads = 20;
	gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");

	const TString selectorName    = "./DSelector_pippippimpimpmiss.C+";
	// pi+pi+pi-pi-(p)
	const TString treeName        = "pippippimpimpmiss__B1_T1_U1_Effic_Tree";
	// const TString fileNamePattern = "./data/MCbggen/2017_01-ver03/tree_pippippimpimpmiss__B1_T1_U1_Effic_MCbggen_2017_01-ver03_batch01.root";
	// const TString fileNamePattern = "./data/RD/2018_01-ver02/tree_pippippimpimpmiss__B1_T1_U1_Effic_RD_2018_01-ver02_041003.root";
	// const TString fileNamePattern = "./data/RD/2018_01-ver02/tree_pippippimpimpmiss__B1_T1_U1_Effic_RD_2018_01-ver02_042030.root";
	// const TString fileNamePattern = "./data/RD/2018_01-ver02/tree_pippippimpimpmiss__B1_T1_U1_Effic_RD_2018_01-ver02_042550.root";
	const TString fileNamePattern = "./data/MCbggen/2018_01-ver02/tree_pippippimpimpmiss__B1_T1_U1_Effic_MCbggen_2018_01-ver02.root";
	// pi+pi-(p)
	// const TString selectorName    = "./DSelector_pippimpmiss.C+";
	// const TString treeName        = "pippimpmiss__B1_T1_U1_Effic_Tree";
	// const TString fileNamePattern = "./data/2017_01-ver04/batch02/tree_pippimpmiss__B1_T1_U1_Effic_030730.root";

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
