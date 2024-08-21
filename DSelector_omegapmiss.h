#ifndef DSelector_omegapmiss_h
#define DSelector_omegapmiss_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1F.h"
#include "TH3F.h"


class DSelector_omegapmiss : public DSelector
{
	public:

		DSelector_omegapmiss(TTree* locTree = NULL)
		  : DSelector(locTree),
		    dSidebandSubtractAcc(true)
		{ }
		virtual ~DSelector_omegapmiss(){}

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

		//Step 0: gamma p -> omega (p)
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DKinematicData* dMissingProtonWrapper;

		//Step 1: omega -> pi0 pi+ pi-
		DParticleComboStep* dStep1Wrapper;
		DChargedTrackHypothesis* dPiPlusWrapper;
		DChargedTrackHypothesis* dPiMinusWrapper;

		//Step 2: pi0 -> gamma gamma
		DParticleComboStep* dStep2Wrapper;
		DKinematicData* dDecayingPi0Wrapper;
		DNeutralParticleHypothesis* dPhoton1Wrapper;
		DNeutralParticleHypothesis* dPhoton2Wrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1F* dHist_BeamEnergy;
		TH1F* dHist_SignalTruthBeamEnergy;
		TH1F* dHist_SignalTruthThreePionMass;
		TH1F* dHist_SignalTruthProtonP;
		TH1F* dHist_SignalTruthProtonTheta;
		TH1F* dHist_SignalTruthProtonPhi;
		TH1F* dHist_SignalTruthBeamEnergy_BeamEnergyRange1;
		TH1F* dHist_SignalTruthThreePionMass_BeamEnergyRange1;
		TH1F* dHist_SignalTruthProtonP_BeamEnergyRange1;
		TH1F* dHist_SignalTruthProtonTheta_BeamEnergyRange1;
		TH1F* dHist_SignalTruthProtonPhi_BeamEnergyRange1;
		TH1F* dHist_SignalTruthBeamEnergy_BeamEnergyRange2;
		TH1F* dHist_SignalTruthThreePionMass_BeamEnergyRange2;
		TH1F* dHist_SignalTruthProtonP_BeamEnergyRange2;
		TH1F* dHist_SignalTruthProtonTheta_BeamEnergyRange2;
		TH1F* dHist_SignalTruthProtonPhi_BeamEnergyRange2;

		TH2F* dHist_MissingMassSquaredVsBeamEnergy;
		TH2F* dHist_MissingMassSquaredVsBeamEnergySideband;

		void setTreeKinematics(
			const TLorentzVector& P4,
			const string&         branchPrefix,
			const string&         branchSuffix = "",
			const int             arrayIndex = -1);

		void setTreeDelta(
			const TLorentzVector& missingProtonP4,
			const TLorentzVector& otherTrackP4,
			const int             arrayIndex,
			const string&         branchPrefix,
			const string&         branchSuffix = "");


	ClassDef(DSelector_omegapmiss, 0);
};

void DSelector_omegapmiss::Get_ComboWrappers(void)
{
	//Step 0: gamma p -> omega (p)
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dMissingProtonWrapper = dStep0Wrapper->Get_FinalParticle(1);  // index 0 is the omega

	//Step 1: omega -> pi0 pi+ pi-
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
	dPiPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));  // index 0 is the pi0
	dPiMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep1Wrapper->Get_FinalParticle(2));

	//Step 2: pi0 -> gamma gamma
	dStep2Wrapper = dComboWrapper->Get_ParticleComboStep(2);
	dDecayingPi0Wrapper = dStep2Wrapper->Get_InitialParticle();
	dPhoton1Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep2Wrapper->Get_FinalParticle(0));
	dPhoton2Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep2Wrapper->Get_FinalParticle(1));
}

#endif // DSelector_omegapmiss_h
