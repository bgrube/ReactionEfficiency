#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TStyle.h"


void
plotTOF_MissingHits(
	// const int     periodIndex = 0,
	// const TString path        = "/cache/halld/gluex_simulations/REQUESTED_MC/2017_bggen_batch01_ver03_31_20220210010210pm/root/monitoring_hists"
	// const TString path        = "/cache/halld/gluex_simulations/REQUESTED_MC/2017_bggen_batch02_ver03_31_20220228053254pm/root/monitoring_hists"
	// const TString path        = "/cache/halld/gluex_simulations/REQUESTED_MC/2017_bggen_batch03_ver03_31_20220310021858pm/root/monitoring_hists"
	// const TString path        = "/cache/halld/gluex_simulations/REQUESTED_MC/2017_bggen_batch04_ver03_31_20220331023833pm/root/monitoring_hists"
	// const int     periodIndex = 1,
	// const TString path        = "/cache/halld/gluex_simulations/REQUESTED_MC/S2018_bggen_ver02_23_batch01_20220103074512pm/root/monitoring_hists"
	// const int     periodIndex = 2,
	// const TString path        = "/cache/halld/gluex_simulations/REQUESTED_MC/F2018_ver02_21_bggen_batch01_20211104064237pm/root/monitoring_hists"
	const int     periodIndex = 3,
	const TString path        = "/cache/halld/gluex_simulations/REQUESTED_MC/2019-11_bggen_batch01_3342/root/monitoring_hists"
) {
	gStyle->SetOptStat(0);

	const vector<TString> periodNames      = {"2017_01", "2018_01", "2018_08", "2019_11"};
	const vector<int    > minRunNmbPeriods = {30274,     40857,     50685,     71350};
	const vector<int    > maxRunNmbPeriods = {31057,     42559,     51768,     73266};

	const TString periodName = periodNames[periodIndex];
	const int minRun = minRunNmbPeriods[periodIndex];
	const int maxRun = maxRunNmbPeriods[periodIndex];
	TH1F* hTOFPointsVsRun = new TH1F("hTOFPointsVsRun", Form("%s;Run Number;Number of ToF Points", periodName.Data()), maxRun - minRun + 1, minRun - 0.5, maxRun + 0.5);
	TH1F* hTOFHitsVsRun   = new TH1F("hTOFHitsVsRun",   Form("%s;Run Number;Number of ToF Hits",   periodName.Data()), maxRun - minRun + 1, minRun - 0.5, maxRun + 0.5);

	for (int runNmb = minRun; runNmb <= maxRun; runNmb++) {
		const vector<int> runsToSkip = {51568, 51577, 51582, 51590, 51595, 51636};
		for (const int runToSkip : runsToSkip) {
			if (runNmb == runToSkip) {
				continue;
			}
		}

		const TString fileName = TString(Form("%s/hd_root_bggen_0%d", path.Data(), runNmb)) + ((periodName == "2019_11") ? ".root" : "_000.root");  // read only first chunk of each run for pre 2019_11 data
		TFile* f = TFile::Open(fileName);
		if (not f or f->IsZombie()) {
			continue;
		}
		TH1F* hTOFPoints = (TH1F*)f->Get("Independent/Hist_NumReconstructedObjects/NumTOFPoints");
		TH1F* hTOFHits   = (TH1F*)f->Get("Independent/Hist_NumReconstructedObjects/NumTOFHits");
		const double nmbTOFPoints = hTOFPoints->GetMean();
		const double nmbTOFHits   = hTOFHits->GetMean();
		hTOFPointsVsRun->SetBinContent(hTOFHitsVsRun->GetXaxis()->FindBin(runNmb), nmbTOFPoints);
		hTOFHitsVsRun->SetBinContent  (hTOFHitsVsRun->GetXaxis()->FindBin(runNmb), nmbTOFHits);
	}

	TCanvas* canvTofPoints = new TCanvas("tofPoints", "", 600, 400);
	hTOFPointsVsRun->SetMarkerStyle(20);
	hTOFPointsVsRun->Draw("P");

	TCanvas* canvTofHits = new TCanvas("tofHits", "", 600, 400);
	hTOFHitsVsRun->SetMarkerStyle(20);
	hTOFHitsVsRun->Draw("P");

	canvTofPoints->Print(Form("toffPointsVsRun_%s.pdf", periodName.Data()));
	canvTofHits->Print  (Form("toffHitsVsRun_%s.pdf",   periodName.Data()));
}
