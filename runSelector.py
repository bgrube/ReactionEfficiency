#!/usr/bin/env python3


import functools
import glob
import os
from typing import List

import ROOT

import makeBruFitTree


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


# load required libraries
ROOT.gSystem.Load("libDSelector.so")
ROOT.gROOT.ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C")


# see https://halldweb.jlab.org/wiki/index.php/DSelector#Using_DSelector.27s_with_PROOF-Lite
def runSelector(
  fileNamePattern:  str,
  treeName:         str,
  selectorFileName: str,
  nmbProofThreads:  int  = 20,
  runPROOF:         bool = True,
  nmbEntries:       int  = ROOT.TTree.kMaxEntries,
) -> None:
  """Runs DSelector over data set defined by given file name pattern"""
  # create TChain with input data
  chain = ROOT.TChain(treeName)
  chain.Add(fileNamePattern)
  selector = selectorFileName.rstrip("+") + "+"  # ensure that selector is always compiled
                                                 # don't use "++", otherwise each proof job compiles the script and speed is much reduced
  print(f"processing tree '{treeName}' in file(s) '{fileNamePattern}' using selector '{selector}'")
  # process TChain
  if runPROOF:
    # run with PROOF
    ROOT.gEnv.SetValue("ProofLite.Sandbox", "$PWD/.proof/")
    ROOT.DPROOFLiteManager.Process_Chain(chain, selector, nmbProofThreads)
  else:
    # run interactively
    chain.Process(selector, "", nmbEntries)


if __name__ == "__main__":
  #TODO add command-line interface
  ROOT.gROOT.SetBatch(True)

  deleteFlatTreeFiles = False

  # pi+pi+pi-pi-(p)
  selectorFileName = "./DSelector_pippippimpimpmiss.C"
  # treeBaseName = "pippippimpimpmiss__B1_T1_U1_Effic"
  treeBaseName = "pippippimpimmissprot__B1_T1_U1_Effic"
  # dataSets = [{"type" : "MCbggen", "period" : "2017_01-ver03"}]
  # dataSets = [{"type" : "MCbggen", "period" : "2018_01-ver02"}]
  # dataSets = [{"type" : "RD",      "period" : "2018_01-ver02", "run" : "041003"}]
  # dataSets = [{"type" : "RD",      "period" : "2018_01-ver02", "run" : "042030"}]
  # dataSets = [{"type" : "RD",      "period" : "2018_01-ver02", "run" : "042550"}]
  # for dataSet in dataSets:
  #   dataSet.update({"fileName" : f"./data/{dataSet['type']}/{dataSet['period']}/tree_{treeBaseName}_{dataSet['type']}_{dataSet['period']}"
  #     + (f"_{dataSet['run']}" if "run" in dataSet else "") + ".root"})
  period = "2018_01-ver02"
  # period = "2019_11-ver01"
  inFileNames = sorted(glob.glob(f"./data/RD/{period}/tree_{treeBaseName}_RD_{period}_??????.root"))
  # inFiles = [f"./data/RD/2019_11-ver01/tree_{treeBaseName}_071592.root"]
  dataSets = []
  for inFile in inFileNames:
    runNumber = inFile.split(".")[-2].split("_")[-1]  # extract run number from file name of the form `tree_{treeBaseName}_<run number>.root`
    dataSet = {"type" : "RD", "period" : period, "run" : runNumber}
    dataSet.update({"fileName" : f"./data/{dataSet['type']}/{dataSet['period']}/tree_{treeBaseName}_{dataSet['type']}_{dataSet['period']}_{dataSet['run']}.root"})
    dataSets.append(dataSet)
  # pi+pi-(p)
  # selectorFileName = "./DSelector_pippimpmiss.C"
  # treeBaseName = "pippimpmiss__B1_T1_U1_Effic_Tree"

  # run selector and create flat-tree file
  channel = "pippippimpimpmiss"  # used for output tree and file names
  flatTreeFileNames: List[str] = []
  for dataSet in dataSets:
    runSelector(dataSet["fileName"], f"{treeBaseName}_Tree", selectorFileName)
    # rename output files
    histFileName = f"./data/{dataSet['type']}/{dataSet['period']}/{channel}.{dataSet['type']}_{dataSet['period']}" + (f"_{dataSet['run']}" if "run" in dataSet else "") + ".root"
    print(f"Moving histogram file to '{histFileName}'")
    os.replace(f"{channel}.root", histFileName)
    flatTreeFileName = f"./data/{dataSet['type']}/{dataSet['period']}/{channel}_flatTree.{dataSet['type']}_{dataSet['period']}" + (f"_{dataSet['run']}" if "run" in dataSet else "") + ".root"
    print(f"Moving flat-tree file to '{flatTreeFileName}'")
    os.replace(f"{channel}_flatTree.root", flatTreeFileName)
    flatTreeFileNames.append(flatTreeFileName)
    print()

  # convert flat tree to BruFit tree
  for flatTreeFileName in flatTreeFileNames:
    bruFitTreeFileName = f"{flatTreeFileName}.brufit"
    makeBruFitTree.makeBruFitTree(flatTreeFileName, outputFileName = bruFitTreeFileName)
    if deleteFlatTreeFiles:
      print(f"Removing flat-tree file '{flatTreeFileName}'")
      os.remove(flatTreeFileName)
