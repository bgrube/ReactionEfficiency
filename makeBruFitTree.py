#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import inspect

import ROOT

import makePlots  # defines helper functions to generate histograms from data trees


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)

  # dataset         = None
  dataset         = "bggen_2017_01-ver03"
  # dataset         = "030730"
  inputFileName   = f"../ReactionEfficiency/pippippimpimpmiss_flatTree.{dataset}.root" if dataset else "pippippimpimpmiss_flatTree.root"
  outputFileName  = f"../ReactionEfficiency/{inputFileName}.brufit"
  treeName        = "pippippimpimpmiss"
  branchesToWrite = [
    "MissingMassSquared_Measured",  # fit variable
    "AccidWeightFactor",  # weight for removal of RF accidentals
    # binning variables
    "BeamEnergy",
    "MissingProtonP",
    "MissingProtonTheta",
    "MissingProtonPhi"
  ]
  print(f"Converting tree '{treeName}' in '{inputFileName}' to BruFit format")

  # adds columns with "track-found flag" and with combo ID (needed by BruFit; needs to be of type double)
  # and writes out new tree with only the needed branches
  # works currently only in single-threaded mode
  # see https://root-forum.cern.ch/t/accessing-entry-information-using-rdataframe/52378
  # and https://root.cern/doc/master/df007__snapshot_8C.html
  rdf = ROOT.RDataFrame(treeName, inputFileName) \
            .Define("TrackFound", makePlots.UNUSED_TRACK_FOUND_CONDITION) \
            .Define("ComboID", "(double)rdfentry_") \
            .Snapshot(treeName, outputFileName, ROOT.std.vector[ROOT.std.string](branchesToWrite + ["TrackFound", "ComboID"])) \
            # .Range(0, 100000)
  print(f"Read {rdf.Count().GetValue()} entries from input tree")

  # print tree structure
  outputFile = ROOT.TFile(outputFileName, "READ")
  tree = outputFile.Get(treeName)
  print(f"Wrote BruFit file '{outputFileName}' with tree:")
  tree.Print()
