void
doit(Long64_t nmbEntries = TTree::kMaxEntries)
{
	gROOT->ProcessLine(".x ./Load_DSelector.C");
	TFile* f = nullptr;
	// f = TFile::Open("./volatile/batch02/tree_pippippimpimpmiss__B1_T1_U1_Effic/030730/tree_pippippimpimpmiss__B1_T1_U1_Effic_030730_000.root");
	f = TFile::Open("./volatile/batch02/tree_pippippimpimpmiss__B1_T1_U1_Effic_030730.root");
	pippippimpimpmiss__B1_T1_U1_Effic_Tree->Process("./DSelector_pippippimpimpmiss.C+", "", nmbEntries);

	// f = TFile::Open("./volatile/batch02/tree_pippimpmiss__B1_T1_U1_Effic_030730.root");
	// pippimpmiss__B1_T1_U1_Effic_Tree->Process("./DSelector_pippimpmiss.C+", "", nmbEntries);
}
