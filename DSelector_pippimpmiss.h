#ifndef DSelector_pippimpmiss_h
#define DSelector_pippimpmiss_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1D.h"
// #include "TH2D.h"

class DSelector_pippimpmiss : public DSelector
{
	public:

		DSelector_pippimpmiss(TTree* locTree = NULL)
		  : DSelector(locTree),
		    dSidebandSubtractAcc(true)
		{ }
		virtual ~DSelector_pippimpmiss(){}

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
		DChargedTrackHypothesis* dPiPlusWrapper;
		DChargedTrackHypothesis* dPiMinusWrapper;
		DKinematicData* dMissingProtonWrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1D* dHist_RFWeight;
		TH1D* dHist_MissingMass;
		TH1D* dHist_MissingMassSideband;
		TH1D* dHist_MissingMassSquared;
		TH1D* dHist_MissingMassSquaredSideband;
		TH1D* dHist_BeamEnergy;
		TH2D* dHist_MissingParticle_MomVsTheta;
		TH2D* dHist_MissingParticle_PhiVsTheta;

	ClassDef(DSelector_pippimpmiss, 0);
};

void DSelector_pippimpmiss::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPiPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
	dPiMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
	dMissingProtonWrapper = dStep0Wrapper->Get_FinalParticle(2);
}

#endif // DSelector_pippimpmiss_h
