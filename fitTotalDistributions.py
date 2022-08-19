#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetHistMinimumZero(True);

import fitFunction


# file name : (hist name, output file base name)
inData = {
  "../ReactionEfficiency/volatile/batch02/hd_root_030730.root" : ("pippippimpimpmiss__B1_T1_U1_Effic/Custom_RecoilMass_pRecoil/HistRecoilMass", "justin_30730_mm_fit"),
  "../pmatt/trackeff_Proton_4pi.30370_mm_acc_Pval.root"        : ("MissingMass/MissingMass",                                                    "paul_30730_mm_acc_Pval_fit")}
rebinFactor = 17
mmRange = (0, 2)

inFiles = [ROOT.TFile(inFileName) for inFileName in inData.keys()]
hists = [inFiles[index].Get(names[0]) for index, names in enumerate(inData.values())]
print(hists)
canvs = [ROOT.TCanvas(names[1], names[1], 1470, 891) for names in inData.values()]
hists[1].Rebin(rebinFactor)
hists[1].GetXaxis().SetRangeUser(*mmRange)
for index, hist in enumerate(hists):
  canv = canvs[index]
  canv.cd()
  fitFunction.fitDistribution(hist, particle = "Proton", fitRange = mmRange, forceCommonGaussianMean = False, fitMissingMassSquared = False)
  hist.Draw("E")
  canv.SaveAs(".pdf")
  print("\n\n")
