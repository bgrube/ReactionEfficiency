#!/usr/bin/env python3

"""
Create Test Dataset
===================

Create test dataset with random data.
"""

import sys

import ROOT

def write_dataset(i):
    df = ROOT.ROOT.RDataFrame(1000000)

    df1 = (
        df.Define("mu_pt", "500 * gRandom->Exp(0.1)")
          .Define("mu_eta", "gRandom->Gaus(0, 5)")
          .Define("mu_phi", "gRandom->Uniform(-3.14, 3.14)")
          .Define("jet_pt", "1000 * gRandom->Exp(0.1)")
          .Define("jet_eta", "gRandom->Gaus(0, 5)")
          .Define("jet_phi", "gRandom->Uniform(-3.14, 3.14)")
          .Define("weight", "gRandom->Gaus(1, 0.1)")
    )

    outfilename = f"dataset.{i+1:02}.root"
    print(f"Writing dataset to file '{outfilename}'")
    df1.Snapshot("dataset", outfilename)

def main():
    for i in range(10):
        write_dataset(i)

if __name__ == "__main__":
    ROOT.gROOT.SetBatch()
    sys.exit(main())
