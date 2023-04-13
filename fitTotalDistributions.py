#!/usr/bin/env python3


import ROOT

import fitFunction

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetHistMinimumZero(True)


# file name : (hist name, output file base name)
inData = {
  "../ReactionEfficiency/pippimpmiss.30370_acc.root" : ("MissingMassSquared", "justin_Proton_2pi_mm2_30370_acc_fit"),
  "../pmatt/trackeff_Proton_2pi.30370_acc_Pval.root" : ("MissingMass/MissingMass", "paul_Proton_2pi_mm2_30370_acc_Pval_fit"),
  "../ReactionEfficiency/pippippimpimpmiss.30370_acc.root" : ("MissingMassSquared", "justin_Proton_4pi_mm2_30370_acc_fit"),
  "../pmatt/trackeff_Proton_4pi.30370_acc_Pval.root" : ("MissingMass/MissingMass", "paul_Proton_4pi_mm2_30370_acc_Pval_fit")}
rebinFactor = 20
fitRange = (-0.5, 4)

inFiles = [ROOT.TFile(inFileName) for inFileName in inData.keys()]
hists = [inFiles[index].Get(names[0]) for index, names in enumerate(inData.values())]
canvs = [ROOT.TCanvas(names[1], names[1], 1470, 891) for names in inData.values()]
for index, hist in enumerate(hists):
  hist.Rebin(rebinFactor)
  hist.GetXaxis().SetRangeUser(*fitRange)
  canv = canvs[index]
  canv.cd()
  fitFunction.fitDistribution(hist, particle = "Proton", fitRange = fitRange, forceCommonGaussianMean = False, fitMissingMassSquared = True)
  hist.Draw("E")
  # hist.DrawCopy("HIST E1 SAME")  # ensure points lie above curves
  canv.SaveAs(".pdf")
  print("\n\n")
