#!/usr/bin/env python3


import functools
import glob

import ROOT

import makePlots  # defines helper functions to generate histograms from data trees


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


def makeBruFitTree(
  inputFileName:  str,
  outputFileName: str,
  treeName:       str = "pippippimpimpmiss",
) -> None:
  """Converts tree in given file to BruFit format"""
  print(f"Converting tree '{treeName}' in '{inputFileName}' to BruFit format")
  print(f"!Note! BruFit trees contain a unique-ID column and hence must never be `hadd`ed")
  branchesToWrite = [
    "MissingMassSquared_Measured",  # fit variable
    "AccidWeightFactor",  # weight for removal of RF accidentals
    # cut variables
    "NmbUnusedShowers",
    "BestMissingMatchDistTOF",
    # binning variables
    "BeamEnergy",
    "MissingProtonP",
    "MissingProtonTheta",
    "MissingProtonPhi",
  ]
  # add columns with
  #   track-found flag
  #   combo ID (needed by BruFit; needs to be of type double)
  #   flag that indicates signal process in bggen MC
  # and write out new tree with only the needed branches
  # works currently only in single-threaded mode
  # see https://root-forum.cern.ch/t/accessing-entry-information-using-rdataframe/52378
  # and https://root.cern/doc/master/df007__snapshot_8C.html
  #TODO convert ThrownTopology only for bggen MC
  rdf = ROOT.RDataFrame(treeName, inputFileName) \
            .Define("TrackFound", makePlots.UNUSED_TRACK_FOUND_CONDITION) \
            .Define("ComboID", "(double)rdfentry_") \
            .Define("IsSignal", 'ThrownTopology.GetString() == "2#pi^{#plus}2#pi^{#minus}p"') \
            .Snapshot(treeName, outputFileName, ROOT.std.vector[ROOT.std.string](branchesToWrite + ["TrackFound", "ComboID", "IsSignal"])) \
            # .Range(0, 100000)
  print(f"Read {rdf.Count().GetValue()} entries from input tree")

  # print tree structure
  outputFile = ROOT.TFile(outputFileName, "READ")
  tree = outputFile.Get(treeName)
  print(f"Wrote BruFit file '{outputFileName}' with tree:")
  tree.Print()
  print()


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)

  # dataSets = []
  # dataSets = ["RD_2017_01-ver04_030730"]
  # dataSets = ["MCbggen_2017_01-ver03"]
  # dataSets = ["RD_2018_01-ver02_041003"]
  # dataSets = ["RD_2018_01-ver02_042030"]
  # dataSets = ["RD_2018_01-ver02_042550"]
  # dataSets = ["MCbggen_2018_01-ver02"]
  # if dataSets:
  #   inputFileNames = [f"./pippippimpimpmiss_flatTree.{dataSet}.root" for dataSet in dataSets]
  # else:
  #   inputFileNames = ["./pippippimpimpmiss_flatTree.root"]
  inputFileNames = sorted(glob.glob("./pippippimpimpmiss_flatTree.RD_2018_01-ver02*.root"))
  # inputFileNames = sorted(glob.glob("./pippippimpimpmiss_flatTree.RD_2019_11-ver01*.root"))

  for inputFileName in inputFileNames:
    makeBruFitTree(inputFileName, outputFileName = f"{inputFileName}.brufit")
