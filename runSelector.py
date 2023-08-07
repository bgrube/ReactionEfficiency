#!/usr/bin/env python3


import ROOT


# see https://halldweb.jlab.org/wiki/index.php/DSelector#Using_DSelector.27s_with_PROOF-Lite
def runSelector(
  fileNamePattern:  str,
  treeName:         str,
  selectorFileName: str,
  nmbProofThreads:  int  = 20,
  runPROOF:         bool = True,
  nmbEntries:       int  = ROOT.TTree.kMaxEntries,  # type: ignore
) -> None:
  '''Runs DSelector over data set defined by given file name pattern'''
  # load required libraries
  ROOT.gSystem.Load("libDSelector.so")  # type: ignore
  ROOT.gROOT.ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C")  # type: ignore
  # create TChain with input data
  chain = ROOT.TChain(treeName)  # type: ignore
  chain.Add(fileNamePattern)
  print(f"processing tree '{treeName}' in file(s) '{fileNamePattern}' using selector '{selectorFileName}'")
  selector = selectorFileName if selectorFileName[-1] == "+" else selectorFileName + "+"  # ensure that selector is compiled
  # process TChain
  if runPROOF:
    # run with PROOF
    ROOT.gEnv.SetValue("ProofLite.Sandbox", "$PWD/.proof/")  # type: ignore
    ROOT.DPROOFLiteManager.Process_Chain(chain, selector, nmbProofThreads)  # type: ignore
  else:
    # run interactively
    chain.Process(selector, "", nmbEntries)


if __name__ == "__main__":
  #TODO add command-line interface
  ROOT.gROOT.SetBatch(True)  # type: ignore

  # pi+pi+pi-pi-(p)
  selectorFileName = "./DSelector_pippippimpimpmiss.C"
  # treeName        = "pippippimpimpmiss__B1_T1_U1_Effic_Tree"
  # fileNamePattern = "./data/MCbggen/2017_01-ver03/tree_pippippimpimpmiss__B1_T1_U1_Effic_MCbggen_2017_01-ver03_batch01.root"
  # fileNamePattern = "./data/RD/2018_01-ver02/tree_pippippimpimpmiss__B1_T1_U1_Effic_RD_2018_01-ver02_041003.root"
  # fileNamePattern = "./data/RD/2018_01-ver02/tree_pippippimpimpmiss__B1_T1_U1_Effic_RD_2018_01-ver02_042030.root"
  # fileNamePattern = "./data/RD/2018_01-ver02/tree_pippippimpimpmiss__B1_T1_U1_Effic_RD_2018_01-ver02_042550.root"
  # fileNamePattern = "./data/MCbggen/2018_01-ver02/tree_pippippimpimpmiss__B1_T1_U1_Effic_MCbggen_2018_01-ver02.root"
  treeName        = "pippippimpimmissprot__B1_T1_U1_Effic_Tree"
  fileNamePattern = "./data/RD/2019_11-ver01/tree_pippippimpimmissprot__B1_T1_U1_Effic_071592.root"
  # pi+pi-(p)
  # selectorFileName = "./DSelector_pippimpmiss.C"
  # treeName         = "pippimpmiss__B1_T1_U1_Effic_Tree"
  # fileNamePattern  = "./data/2017_01-ver04/batch02/tree_pippimpmiss__B1_T1_U1_Effic_030730.root"
  runSelector(fileNamePattern, treeName, selectorFileName)

  #TODO rename output files
