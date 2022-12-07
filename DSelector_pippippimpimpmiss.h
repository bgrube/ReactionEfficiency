#ifndef DSelector_pippippimpimpmiss_h
#define DSelector_pippippimpimpmiss_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1D.h"
// #include "TH2D.h"

class DSelector_pippippimpimpmiss : public DSelector
{
	public:

		DSelector_pippippimpimpmiss(TTree* locTree = NULL)
		  : DSelector(locTree),
		    dSidebandSubtractAcc(true)
		{ }
		virtual ~DSelector_pippippimpimpmiss(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Get_ComboWrappers(void);
		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO

		bool dIsMC;

		bool dSidebandSubtractAcc;  // en/disables RF sideband subtraction

		// ANALYZE CUT ACTIONS
		// // Automatically makes mass histograms where one cut is missing
		DHistogramAction_AnalyzeCutActions* dAnalyzeCutActions;

		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS

		//Step 0
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DChargedTrackHypothesis* dPiPlus1Wrapper;
		DChargedTrackHypothesis* dPiPlus2Wrapper;
		DChargedTrackHypothesis* dPiMinus1Wrapper;
		DChargedTrackHypothesis* dPiMinus2Wrapper;
		DKinematicData* dMissingProtonWrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1D* dHist_BeamEnergy;
		TH1D* dHist_RFWeight;
		TH1D* dHist_MissingMass;
		TH1D* dHist_MissingMassSideband;
		TH2D* dHist_MissingParticle_MomVsTheta;
		TH2D* dHist_MissingParticle_PhiVsTheta;
		TH2D* dHist_MissingParticle_MomVsTheta_Measured;
		TH2D* dHist_MissingParticle_PhiVsTheta_Measured;

		TH1D* dHist_locMissingDeltaP;
		TH1D* dHist_locMissingDeltaPOverP;
		TH1D* dHist_locMissingDeltaTheta;
		TH1D* dHist_locMissingDeltaPhi;
		TH2D* dHist_locMissingProtonP_kinFitVsUnused;
		TH2D* dHist_locMissingProtonTheta_kinFitVsUnused;
		TH2D* dHist_locMissingProtonPhi_kinFitVsUnused;

		TH1D* dHist_MissingMassSquared;
		TH1D* dHist_MissingMassSquared_Found;
		TH1D* dHist_MissingMassSquared_Missing;
		TH1D* dHist_MissingMassSquaredSideband;
		TH2D* dHist_MissingMassSquaredVsBeamEnergy;
		TH2D* dHist_MissingMassSquaredVsBeamEnergy_Found;
		TH2D* dHist_MissingMassSquaredVsBeamEnergy_Missing;
		TH2D* dHist_MissingMassSquaredVsBeamEnergySideband;
		TH2D* dHist_MissingMassSquaredVsBeamEnergySideband_Found;
		TH2D* dHist_MissingMassSquaredVsBeamEnergySideband_Missing;

		// bggen MC histograms
		TH1D* dHist_ThrownTopologies;
		TH1D* dHist_ThrownTopologies_Found;
		TH1D* dHist_ThrownTopologies_Missing;
		map<TString, TH1D*> dHist_MissingMassSquared_ThrownTopology;
		map<TString, TH1D*> dHist_MissingMassSquared_ThrownTopology_Found;
		map<TString, TH1D*> dHist_MissingMassSquared_ThrownTopology_Missing;

	ClassDef(DSelector_pippippimpimpmiss, 0);
};

void DSelector_pippippimpimpmiss::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPiPlus1Wrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
	dPiPlus2Wrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
	dPiMinus1Wrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
	dPiMinus2Wrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(3));
	dMissingProtonWrapper = dStep0Wrapper->Get_FinalParticle(4);
}

#endif // DSelector_pippippimpimpmiss_h
