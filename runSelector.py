#!/usr/bin/env python3


import functools
import glob
import os

import ROOT


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


# load required libraries
ROOT.gSystem.Load("libDSelector.so")  # type: ignore
ROOT.gROOT.ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C")  # type: ignore


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
  # create TChain with input data
  chain = ROOT.TChain(treeName)  # type: ignore
  chain.Add(fileNamePattern)
  selector = selectorFileName.rstrip("+") + "+"  # ensure that selector is always compiled
                                                 # don't use "++", otherwise each proof job compiles the script and speed is much reduced
  print(f"processing tree '{treeName}' in file(s) '{fileNamePattern}' using selector '{selector}'")
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
  # treeName = "pippippimpimpmiss__B1_T1_U1_Effic_Tree"
  # # dataSets = [{"type" : "MCbggen", "period" : "2017_01-ver03"}]
  # dataSets = [{"type" : "RD",      "period" : "2018_01-ver02", "run" : "041003"}]
  # # dataSets = [{"type" : "RD",      "period" : "2018_01-ver02", "run" : "042030"}]
  # # dataSets = [{"type" : "RD",      "period" : "2018_01-ver02", "run" : "042550"}]
  # # dataSets = [{"type" : "MCbggen", "period" : "2018_01-ver02"}]
  # for dataSet in dataSets:
  #   dataSet.update({"fileName" : f"./data/{dataSet['type']}/{dataSet['period']}/tree_pippippimpimpmiss__B1_T1_U1_Effic_{dataSet['type']}_{dataSet['period']}"
  #   + (f"_{dataSet['run']}" if "run" in dataSet else "") + ".root"})
  treeName = "pippippimpimmissprot__B1_T1_U1_Effic_Tree"
  # inFiles = sorted(glob.glob("./data/RD/2019_11-ver01/*.root"))
  inFiles = ["./data/RD/2019_11-ver01/tree_pippippimpimmissprot__B1_T1_U1_Effic_071592.root"]
  dataSets = []
  for inFile in inFiles:
    runNumber = inFile.split(".")[-2].split("_")[-1]  # extract run number from file name of the form `tree_pippippimpimmissprot__B1_T1_U1_Effic_<run number>.root`
    dataSet = {"type" : "RD", "period" : "2019_11-ver01", "run" : runNumber}
    dataSet.update({"fileName" : f"./data/{dataSet['type']}/{dataSet['period']}/tree_pippippimpimmissprot__B1_T1_U1_Effic_{dataSet['run']}.root"})
    dataSets.append(dataSet)
  # pi+pi-(p)
  # selectorFileName = "./DSelector_pippimpmiss.C"
  # treeName = "pippimpmiss__B1_T1_U1_Effic_Tree"

  for dataSet in dataSets:
    runSelector(dataSet["fileName"], treeName, selectorFileName)
    # rename output files
    channel = "pippippimpimpmiss"
    histFileName = f"./{channel}.{dataSet['type']}_{dataSet['period']}" + (f"_{dataSet['run']}" if "run" in dataSet else "") + ".root"
    print(f"Writing histogram file to '{histFileName}'")
    os.replace(f"{channel}.root", histFileName)
    treeFileName = f"./{channel}_flatTree.{dataSet['type']}_{dataSet['period']}" + (f"_{dataSet['run']}" if "run" in dataSet else "") + ".root"
    print(f"Writing tree file to '{treeFileName}'")
    os.replace(f"{channel}_flatTree.root", treeFileName)
    print()

  #TODO convert to BruFit tree right inplace
