#!/usr/bin/env python3


from __future__ import annotations

from collections.abc import Sequence
import functools
import glob

import ROOT

import makePlots  # defines helper functions to generate histograms from data trees


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


def makeBruFitTree(
  inputFileNames: Sequence[str],
  outputFileName: str,
  treeName:       str = "pippippimpimpmiss",
  trueTopology:   str = "2#pi^{#plus}2#pi^{#minus}p",
) -> None:
  """Reads given files and writes out tree in BruFit format"""
  print(f"Converting tree '{treeName}' to BruFit format merging the data from {len(inputFileNames)} files: {inputFileNames}")
  print(f"!Note! BruFit trees contain a unique-ID column and hence must never be `hadd`ed")
  branchesToWrite = [
    "MissingMassSquared_Measured",  # fit variable
    "AccidWeightFactor",  # weight for removal of RF accidentals
    # cut variables
    "NmbUnusedShowers",
    # "BestMissingMatchDistTOF",
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
  rdf = ROOT.RDataFrame(treeName, inputFileNames) \
            .Define("TrackFound", makePlots.UNUSED_TRACK_FOUND_CONDITION) \
            .Define("ComboID", "(double)rdfentry_") \
            .Define("IsSignal", f'ThrownTopology.GetString() == "{trueTopology}"') \
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

  # treeName     = "pippippimpimpmiss"
  # trueTopology = "2#pi^{#plus}2#pi^{#minus}p"
  treeName     = "omegapmiss"
  trueTopology = "2#gamma#pi^{#plus}#pi^{#minus}p[#pi^{0},#omega]"
  dataPeriod   = "2017_01-ver03"
  # dataPeriod   = "2018_01-ver02"
  # dataPeriod   = "2018_08-ver02"

  inputFileNames = sorted(glob.glob(f"./data/RD/{dataPeriod}/{treeName}_flatTree.RD_{dataPeriod}*.root"))

  # for inputFileName in inputFileNames:
  #   makeBruFitTree(inputFileName, outputFileName = f"{inputFileName}.brufit")
  makeBruFitTree(
    inputFileNames = inputFileNames,
    outputFileName = f"./data/RD/{dataPeriod}/{treeName}_flatTree.RD_{dataPeriod}.root.brufit",
    treeName       = treeName,
    trueTopology   = trueTopology,
  )
