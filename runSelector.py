#!/usr/bin/env python3


from __future__ import annotations

import functools
import glob
import os

import ROOT

from makeBruFitTree import makeBruFitTree


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


# load required libraries
# ROOT.gSystem.Load("libDSelector.so")
#FIXME when DSelector is compiled inside Python environment, ACLiC also links some libraries from Python packages
#      this makes PROOF jobs fail complaining that they cannot find `libopenblas64_p-r0-0cf96a72.3.23.dev.so` (on ifarm)
#      this is a library from the user-installed numpy package
# print(f" GetLibraries() = {ROOT.gSystem.GetLibraries()}")
# adding the numpy directory to the dynamic library path has no effect
# ROOT.gSystem.AddDynamicPath(f"{os.path.expanduser('~')}/.local/lib/python3.9/site-packages/numpy.libs")
# print(f" GetDynamicPath() = {ROOT.gSystem.GetDynamicPath()}")
# workaround is to add the numpy directory to LD_LIBRARY_PATH
libraryPath = os.environ.get("LD_LIBRARY_PATH", "")
os.environ["LD_LIBRARY_PATH"] = f"{os.path.expanduser('~')}/.local/lib/python3.9/site-packages/numpy.libs" + (f":{libraryPath}" if libraryPath else "")
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
                                                 # don't use "++", otherwise each PROOF job compiles the script and speed is much reduced
  print(f"processing tree '{treeName}' in file(s) '{fileNamePattern}' using selector '{selector}'")
  # process TChain
  if runPROOF:
    # run with PROOF
    proofDirectory = f"{os.getcwd()}/.proof"
    # ROOT.DPROOFLiteManager.Set_SandBox(proofDirectory)  # has no effect
    ROOT.gEnv.SetValue("ProofLite.Sandbox", proofDirectory)
    ROOT.gEnv.SetValue("ProofServ.Sandbox", proofDirectory)
    ROOT.DPROOFLiteManager.Process_Chain(chain, selector, nmbProofThreads)
  else:
    # run interactively
    chain.Process(selector, "", nmbEntries)


if __name__ == "__main__":
  #TODO add command-line interface
  ROOT.gROOT.SetBatch(True)

  runPROOF = True
  nmbProofThreads = 10
  createBruFitTree = True
  deleteFlatTreeFiles = False
  writeBruFitFilesForAllInputFiles = False

  # define channel
  # # pi+pi+pi-pi-(p)
  # channel          = "pippippimpimpmiss"  # used for output tree and file names
  # treeName         = "pippippimpimmissprot__B1_T1_U1_Effic"
  # trueTopology     = "2#pi^{#plus}2#pi^{#minus}p"
  # cutVariableNames = ("NmbUnusedShowers",)
  # omega(p)
  channel          = "omegapmiss"
  treeName         = "omegamissprot__B1_T1_U1_Effic"
  trueTopology     = "2#gamma#pi^{#plus}#pi^{#minus}p[#pi^{0},#omega]"
  cutVariableNames = ("NmbUnusedShowers", "ThreePionMass")
  # # pi+pi-(p)
  # channel          = "pippimpmiss"
  # treeName         = "pippimpmiss__B1_T1_U1_Effic"
  # trueTopology     = "#pi^{#plus}#pi^{#minus}p"
  # cutVariableNames = ("NmbUnusedShowers",)
  selectorFileName = f"./DSelector_{channel}.C"

  # define input files
  dataPeriods = [
    "2017_01-ver03",
    # "2018_01-ver02",
    # "2018_08-ver02",
    # "2019_11-ver01",
    # "2017_01-ver03_goodToF",
    # "2018_01-ver02_goodToF",
    # "2018_08-ver02_goodToF",
  ]
  for dataType in ("MCbggen", "RD"):
    for dataPeriod in dataPeriods:
      dataDir = f"./data/{dataType}/{dataPeriod}"
      inFileNamePattern = f"{dataDir}/tree_{treeName}_{dataType}_{dataPeriod}*.root"
      print(f"Running DSelector over files '{inFileNamePattern}'")
      inFileNames = sorted(glob.glob(inFileNamePattern))
      assert inFileNames, f"Did not find any files matching the name pattern {inFileNamePattern}"

      # run selector and create flat-tree files
      flatTreeFileNames: list[str] = []
      for inFileName in inFileNames:
        runNumber = inFileName.split(".")[-2].split("_")[-1]  # extract run number from file name of the form `tree_{treeName}_<run number>.root`
        if not runNumber.isnumeric():
          runNumber = None
        runSelector(inFileName, f"{treeName}_Tree", selectorFileName, nmbProofThreads, runPROOF)
        # rename output files
        histOutFileName = f"{dataDir}/{channel}.{dataType}_{dataPeriod}" + ("" if runNumber is None else f"_{runNumber}") + ".root"
        print(f"Moving histogram file to '{histOutFileName}'")
        os.replace(f"{channel}.root", histOutFileName)
        flatTreeOutFileName = f"{dataDir}/{channel}_flatTree.{dataType}_{dataPeriod}" + ("" if runNumber is None else f"_{runNumber}") + ".root"
        print(f"Moving flat-tree file to '{flatTreeOutFileName}'")
        os.replace(f"{channel}_flatTree.root", flatTreeOutFileName)
        flatTreeFileNames.append(flatTreeOutFileName)
        print()

      if createBruFitTree:
        # merge all flat-tree files into single BruFit file
        makeBruFitTree(
          inputFileNames   = flatTreeFileNames,
          outputFileName   = f"{dataDir}/{channel}_flatTree.{dataType}_{dataPeriod}.root.brufit",
          treeName         = channel,
          trueTopology     = trueTopology,
          cutVariableNames = cutVariableNames,
        )
        if writeBruFitFilesForAllInputFiles and len(flatTreeFileNames) > 1:
          # write BruFit-tree file for each flat-tree file
          for flatTreeFileName in flatTreeFileNames:
            bruFitTreeFileName = f"{flatTreeFileName}.brufit"
            makeBruFitTree(
              inputFileNames   = flatTreeFileName,
              outputFileName   = bruFitTreeFileName,
              treeName         = channel,
              trueTopology     = trueTopology,
              cutVariableNames = cutVariableNames,
            )
        if deleteFlatTreeFiles:
          for flatTreeFileName in flatTreeFileNames:
            print(f"Removing flat-tree file '{flatTreeFileName}'")
            os.remove(flatTreeFileName)
