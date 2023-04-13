#!/usr/bin/env python3

"""
RDataFrame Test 2
=================

Second test: define histograms and run RDataFrame event loop in child process. Call
ROOT.ROOT.EnableImplicitMT() in child process.
"""

import argparse
import multiprocessing
import sys

import ROOT


def parse_args():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--jobs", type=int, default=1, help="Number of jobs")
    parser.add_argument("-o", "--outfile", help="Path to output file")

    args = parser.parse_args()

    return args


@profile
def fill_hists(args):
    with multiprocessing.Manager() as manager:
        hists = manager.dict()
        proc = multiprocessing.Process(
            target=fill_hists_worker,
            args=(hists, args)
        )
        proc.start()
        proc.join()

    return hists


def fill_hists_worker(hists, args):
    if args.jobs > 1:
        ROOT.ROOT.EnableImplicitMT(args.jobs)
    
    # Read data created by `create_test_dataset.py`
    df = ROOT.ROOT.RDataFrame("dataset", "dataset.*.root")

    df_local = (
        df.Filter("mu_pt > 25", "Muon pT > 25 GeV")
          .Filter("abs(mu_eta) < 5", "Muon |eta| < 5")
          .Filter("jet_pt > 30", "Jet pT > 30 GeV")
          .Filter("abs(jet_eta) < 5", "Jet |eta| < 5")
          .Define("deltaR_lj", "ROOT::VecOps::DeltaR(mu_eta, jet_eta, mu_phi, jet_phi)")
          .Filter("deltaR_lj > 0.5", "deltaR_lj > 0.5")
    )

    rhists = []
    rhists.append(df_local.Histo1D(("mu_pt", "", 30, 25, 325), "mu_pt", "weight"))
    rhists.append(df_local.Histo1D(("mu_eta", "", 20, -5, 5), "mu_eta", "weight"))
    rhists.append(df_local.Histo1D(("mu_phi", "", 30, -3.14, 3.14), "mu_phi", "weight"))
    rhists.append(df_local.Histo1D(("jet_pt", "", 30, 30, 630), "jet_pt", "weight"))
    rhists.append(df_local.Histo1D(("jet_eta", "", 20, -5, 5), "jet_eta", "weight"))
    rhists.append(df_local.Histo1D(("jet_phi", "", 30, -3.14, 3.14), "jet_phi", "weight"))
    rhists.append(df_local.Histo1D(("deltaR_lj", "", 30, 0, 10), "deltaR_lj", "weight"))

    for rhist in rhists:
        hist = rhist.GetValue()
        hists[hist.GetName()] = hist


@profile
def write_hists(hists, outfilename):
    # Write histograms to file
    print(f"Writing histograms to file '{outfilename}'")
    outfile = ROOT.TFile.Open(outfilename, "RECREATE")
    try:
        for hist in hists.values():
            hist.Write()
    finally:
        outfile.Close()


def main():
    args = parse_args()

    # Fill hists many times to simulate many datasets
    for _ in range(1, 11):
        hists = fill_hists(args)

        if args.outfile:
            write_hists(hists, args.outfile)


if __name__ == "__main__":
    ROOT.gROOT.SetBatch()
    sys.exit(main())
