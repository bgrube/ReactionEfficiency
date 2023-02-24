#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import inspect

import ROOT

import makePlots  # defines helper functions to generate histograms from data trees


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)

  # dataset        = None
  dataset        = "bggen_2017_01-ver03"
  # dataset        = "030730"
  inputFileName  = f"../ReactionEfficiency/pippippimpimpmiss_flatTree.{dataset}.root" if dataset else "pippippimpimpmiss_flatTree.root"
  outputFileName = f"../ReactionEfficiency/{inputFileName}.brufit"
  treeName       = "pippippimpimpmiss"
  brufitVariable = "MissingMassSquared_Measured"  # branch name to be written


  # adds columns with "track found flag" and with event ID (needed by BruFit; needs to be of type double)
  # and writes out new tree with only the needed branches
  # works currently only in single-threaded mode
  # see https://root-forum.cern.ch/t/accessing-entry-information-using-rdataframe/52378
  # and https://root.cern/doc/master/df007__snapshot_8C.html
  branchesToWrite = ROOT.std.vector[ROOT.std.string]([brufitVariable, "TrackFound", "EventID"])
  ROOT.RDataFrame(treeName, inputFileName) \
      .Define("TrackFound", makePlots.UNUSED_TRACK_FOUND_CONDITION) \
      .Define("EventID", "(double)rdfentry_") \
      .Snapshot(treeName, outputFileName, branchesToWrite) \
      # .Range(0, 100000)

  outputFile = ROOT.TFile(outputFileName, "READ")
  tree = outputFile.Get(treeName)
  branches = [b.GetName() for b in tree.GetListOfBranches()]
  print(f"Wrote file '{outputFileName}' with tree '{treeName}' with branches {branches}")
