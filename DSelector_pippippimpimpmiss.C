#include <cassert>

#include "TObjString.h"

#include "DSelector_pippippimpimpmiss.h"


// workaround: provide dummy conversion functions for unnecessary static cast to TLorentzVector in DSelector/DTreeInterface.h:662
namespace {

	struct myTObjString : TObjString
	{
		myTObjString(const TString& s = "")
			: TObjString(s)
		{ }

		operator TLorentzVector() const { return TLorentzVector(); }
	};

}


void DSelector_pippippimpimpmiss::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	const string baseName = "pippippimpimpmiss";
	dOutputFileName = baseName + ".root";  // "" for none
	dOutputTreeFileName = "";  // "" for none
	dFlatTreeFileName = baseName + "_flatTree.root";  // output flat tree (one combo per tree entry), "" for none
	dFlatTreeName = baseName;  // if blank, default name will be chosen
	dSaveDefaultFlatBranches = false;  // False: don't save default branches, reduce disk footprint.
	// dSaveTLorentzVectorsAsFundamentaFlatTree = false;  // Default (or false): save particles as TLorentzVector objects. True: save as four doubles instead.

	cout << "histogram output file name = '" << dOutputFileName << "'" << endl
	     << "flat tree output file name = '" << dFlatTreeFileName << "'" << endl
	     << "flat tree name = '"             << dFlatTreeName << "'" << endl
	     << "accidental subtraction = "      << dSidebandSubtractAcc << endl;

	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	Get_ComboWrappers();
	dPreviousRunNumber = 0;

	/*********************************** EXAMPLE USER INITIALIZATION: ANALYSIS ACTIONS **********************************/

	// EXAMPLE: Create deque for histogramming particle masses:
	// // For histogramming the phi mass in phi -> K+ K-
	// // Be sure to change this and dAnalyzeCutActions to match reaction
	// std::deque<Particle_t> MyPhi;
	// MyPhi.push_back(KPlus); MyPhi.push_back(KMinus);

	//ANALYSIS ACTIONS: //Executed in order if added to dAnalysisActions
	//false/true below: use measured/kinfit data

	//PID
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false));
	//below: value: +/- N ns, Unknown: All PIDs, SYS_NULL: all timing systems
	//dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, KPlus, SYS_BCAL));

	//PIDFOM (for charged tracks)
	dAnalysisActions.push_back(new DHistogramAction_PIDFOM(dComboWrapper));
	//dAnalysisActions.push_back(new DCutAction_PIDFOM(dComboWrapper, KPlus, 0.1));
	//dAnalysisActions.push_back(new DCutAction_EachPIDFOM(dComboWrapper, 0.1));

	//MASSES
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, 0, {PiPlus, PiMinus}, 1000, 0, 2, "TwoPi"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, 0, {PiPlus, PiPlus, PiMinus, PiMinus}, 1000, 0, 2, "FourPi"));
	dAnalysisActions.push_back(new DHistogramAction_MissingMass       (dComboWrapper, false, 5000, -0.5,  4.5));
	dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 5000, -0.5,  4.5));
	dAnalysisActions.push_back(new DHistogramAction_MissingP4(dComboWrapper, false, ""));

	//KINFIT RESULTS
	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

	//CUT MISSING MASS
	//dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.02));

	//CUT ON SHOWER QUALITY
	//dAnalysisActions.push_back(new DCutAction_ShowerQuality(dComboWrapper, SYS_FCAL, 0.5));

	//BEAM ENERGY
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
	//dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.2, 8.8));  // Coherent peak for runs in the range 30000-59999

	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

	// // ANALYZE CUT ACTIONS
	// // Change MyPhi to match reaction
	// dAnalyzeCutActions = new DHistogramAction_AnalyzeCutActions(dAnalysisActions, dComboWrapper, false, 0, MyPhi, 1000, 0.9, 2.4, "CutActionEffect");

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();
	// dAnalyzeCutActions->Initialize(); // manual action, must call Initialize()

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	//EXAMPLE MANUAL HISTOGRAMS:
	dHist_BeamEnergy = new TH1F("BeamEnergy", ";Beam Energy (GeV)", 1000,  2, 12);

	gDirectory->mkdir("MissingMassSquared", "MissingMassSquared");
	gDirectory->cd("MissingMassSquared");
	dHist_MissingMassSquaredVsBeamEnergy         = new TH2F("MissingMassSquaredVsBeamEnergy",         ";Beam Energy (GeV); Missing Mass Squared (GeV/c^{2})^{2}", 500, 2, 12, 5000, -0.5, 4.5);
	dHist_MissingMassSquaredVsBeamEnergySideband = new TH2F("MissingMassSquaredVsBeamEnergySideband", ";Beam Energy (GeV); Missing Mass Squared (GeV/c^{2})^{2}", 500, 2, 12, 5000, -0.5, 4.5);
	gDirectory->cd("..");

	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - MAIN TREE *************************/

	//EXAMPLE MAIN TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dTreeInterface->Create_Branch_Fundamental<Int_t>("my_int"); //fundamental = char, int, float, double, etc.
	dTreeInterface->Create_Branch_FundamentalArray<Int_t>("my_int_array", "my_int");
	dTreeInterface->Create_Branch_FundamentalArray<Float_t>("my_combo_array", "NumCombos");
	dTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("my_p4");
	dTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("my_p4_array");
	*/

	/************************** USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - FLAT TREE *************************/

	// RECOMMENDED: CREATE ACCIDENTAL WEIGHT BRANCH
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("AccidWeightFactor");

	//EXAMPLE FLAT TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>       ("BeamEnergy");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>       ("MissingMassSquared");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>       ("MissingProtonP");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>       ("MissingProtonTheta");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>       ("MissingProtonPhi");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>       ("MissingMassSquared_Measured");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>       ("MissingProtonP_Measured");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>       ("MissingProtonTheta_Measured");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>       ("MissingProtonPhi_Measured");
	dFlatTreeInterface->Create_Branch_Fundamental<Int_t>          ("NmbUnusedShowers");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>       ("EnergyUnusedShowers");
	dFlatTreeInterface->Create_Branch_NoSplitTObject<myTObjString>("ThrownTopology");
	//TODO write out only best match instead of all candidates?
	dFlatTreeInterface->Create_Branch_Fundamental<Int_t>          ("NmbTruthTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("TruthP",                    "NmbTruthTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("TruthTheta",                "NmbTruthTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("TruthPhi",                  "NmbTruthTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("TruthDeltaP",               "NmbTruthTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("TruthDeltaPOverP",          "NmbTruthTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("TruthDeltaTheta",           "NmbTruthTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("TruthDeltaPhi",             "NmbTruthTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("TruthDeltaP_Measured",      "NmbTruthTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("TruthDeltaPOverP_Measured", "NmbTruthTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("TruthDeltaTheta_Measured",  "NmbTruthTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("TruthDeltaPhi_Measured",    "NmbTruthTracks");
	dFlatTreeInterface->Create_Branch_Fundamental<Int_t>          ("NmbUnusedTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("UnusedP",                    "NmbUnusedTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("UnusedTheta",                "NmbUnusedTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("UnusedPhi",                  "NmbUnusedTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("UnusedDeltaP",               "NmbUnusedTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("UnusedDeltaPOverP",          "NmbUnusedTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("UnusedDeltaTheta",           "NmbUnusedTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("UnusedDeltaPhi",             "NmbUnusedTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("UnusedDeltaP_Measured",      "NmbUnusedTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("UnusedDeltaPOverP_Measured", "NmbUnusedTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("UnusedDeltaTheta_Measured",  "NmbUnusedTracks");
	dFlatTreeInterface->Create_Branch_FundamentalArray<Double_t>  ("UnusedDeltaPhi_Measured",    "NmbUnusedTracks");

	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want

	/************************************** DETERMINE IF ANALYZING SIMULATED DATA *************************************/

	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);
}


void
printTrack(const DKinematicData& track)
{
	cout << "ID = "  << track.Get_ID() << ", "
	     << "PID = " << track.Get_PID() << ", "
	     << "p_z = " << track.Get_P4().Pz() << endl;
}


void
DSelector_pippippimpimpmiss::fillTreeKinematics(
	const TLorentzVector& P4,
	const string&         branchPrefix,
	const string&         branchSuffix,
	const int             arrayIndex)
{
	// kinematic variables
	const TVector3 locP3    = P4.Vect();
	const double   locP     = locP3.Mag();
	const double   locTheta = locP3.Theta() * TMath::RadToDeg();
	const double   locPhi   = locP3.Phi()   * TMath::RadToDeg();
	if (arrayIndex < 0) {
		// fill scalar varable
		dFlatTreeInterface->Fill_Fundamental<Double_t>(branchPrefix + "P"     + branchSuffix, locP);
		dFlatTreeInterface->Fill_Fundamental<Double_t>(branchPrefix + "Theta" + branchSuffix, locTheta);
		dFlatTreeInterface->Fill_Fundamental<Double_t>(branchPrefix + "Phi"   + branchSuffix, locPhi);
	} else {
		// fill array variable
		dFlatTreeInterface->Fill_Fundamental<Double_t>(branchPrefix + "P"     + branchSuffix, locP,     arrayIndex);
		dFlatTreeInterface->Fill_Fundamental<Double_t>(branchPrefix + "Theta" + branchSuffix, locTheta, arrayIndex);
		dFlatTreeInterface->Fill_Fundamental<Double_t>(branchPrefix + "Phi"   + branchSuffix, locPhi,   arrayIndex);
	}
}


void
DSelector_pippippimpimpmiss::fillTreeDelta(
	const TLorentzVector& missingProtonP4,
	const TLorentzVector& otherTrackP4,
	const int             arrayIndex,
	const string&         branchPrefix,
	const string&         branchSuffix)
{
	// kinematic variables of missing proton
	const TVector3 locMissingProtonP3    = missingProtonP4.Vect();
	const double   locMissingProtonP     = locMissingProtonP3.Mag();
	const double   locMissingProtonTheta = locMissingProtonP3.Theta() * TMath::RadToDeg();
	const double   locMissingProtonPhi   = locMissingProtonP3.Phi()   * TMath::RadToDeg();
	// kinematic variables of other track
	const TVector3 locOtherTrackP3    = otherTrackP4.Vect();
	const double   locOtherTrackP     = locOtherTrackP3.Mag();
	const double   locOtherTrackTheta = locOtherTrackP3.Theta() * TMath::RadToDeg();
	const double   locOtherTrackPhi   = locOtherTrackP3.Phi()   * TMath::RadToDeg();
	// calculate deviation from missing-proton momentum in spherical coordinates
	const double locOtherDeltaP      = locOtherTrackP - locMissingProtonP;
	const double locOtherDeltaPOverP = locOtherDeltaP / locMissingProtonP;
	const double locOtherDeltaTheta  = locOtherTrackTheta - locMissingProtonTheta;
	double       locOtherDeltaPhi    = locOtherTrackPhi - locMissingProtonPhi;
	while (locOtherDeltaPhi > 180) {
		locOtherDeltaPhi -= 360;
	}
	while (locOtherDeltaPhi < -180) {
		locOtherDeltaPhi += 360;
	}
	dFlatTreeInterface->Fill_Fundamental<Double_t>(branchPrefix + "DeltaP"      + branchSuffix, locOtherDeltaP,      arrayIndex);
	dFlatTreeInterface->Fill_Fundamental<Double_t>(branchPrefix + "DeltaPOverP" + branchSuffix, locOtherDeltaPOverP, arrayIndex);
	dFlatTreeInterface->Fill_Fundamental<Double_t>(branchPrefix + "DeltaTheta"  + branchSuffix, locOtherDeltaTheta,  arrayIndex);
	dFlatTreeInterface->Fill_Fundamental<Double_t>(branchPrefix + "DeltaPhi"    + branchSuffix, locOtherDeltaPhi,    arrayIndex);
}


Bool_t DSelector_pippippimpimpmiss::Process(Long64_t locEntry)
{
	// The Process() function is called for each entry in the tree. The entry argument
	// specifies which entry in the currently loaded tree is to be processed.
	//
	// This function should contain the "body" of the analysis. It can contain
	// simple or elaborate selection criteria, run algorithms on the data
	// of the event and typically fill histograms.
	//
	// The processing can be stopped by calling Abort().
	// Use fStatus to set the return value of TTree::Process().
	// The return value is currently not used.

	if (locEntry % 100000 == 0)
		cout << "Event: " << locEntry << endl;

	//CALL THIS FIRST
	DSelector::Process(locEntry); //Gets the data from the tree for the entry
	//cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
	//TLorentzVector locProductionX4 = Get_X4_Production();

	/******************************************** GET POLARIZATION ORIENTATION ******************************************/

	//Only if the run number changes
	//RCDB environment must be setup in order for this to work! (Will return false otherwise)
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dPreviousRunNumber = locRunNumber;
	}

	/********************************************* SETUP UNIQUENESS TRACKING ********************************************/

	//ANALYSIS ACTIONS: Reset uniqueness tracking for each action
	//For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
	Reset_Actions_NewEvent();
	// dAnalyzeCutActions->Reset_NewEvent(); // manual action, must call Reset_NewEvent()

	//PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING
		//Sometimes, some content is the exact same between one combo and the next
			//e.g. maybe two combos have different beam particles, but the same data for the final-state
		//When histogramming, you don't want to double-count when this happens: artificially inflates your signal (or background)
		//So, for each quantity you histogram, keep track of what particles you used (for a given combo)
		//Then for each combo, just compare to what you used before, and make sure it's unique

	//EXAMPLE 1: Particle-specific info:
	set<Int_t> locUsedSoFar_BeamEnergy; //Int_t: Unique ID for beam particles. set: easy to use, fast to search

	//EXAMPLE 2: Combo-specific info:
		//In general: Could have multiple particles with the same PID: Use a set of Int_t's
		//In general: Multiple PIDs, so multiple sets: Contain within a map
		//Multiple combos: Contain maps within a set (easier, faster to search)
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_UnusedTrack;

	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE

	/****************************************** MC DATA: PARSE THROWN TOPOLOGY ******************************************/

	const TString locThrownTopology = Get_ThrownTopologyString();  // empty string if info is not available

	/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

	/*
	Int_t locMyInt = 7;
	dTreeInterface->Fill_Fundamental<Int_t>("my_int", locMyInt);

	TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);
	dTreeInterface->Fill_TObject<TLorentzVector>("my_p4", locMyP4);

	for(int loc_i = 0; loc_i < locMyInt; ++loc_i)
		dTreeInterface->Fill_Fundamental<Int_t>("my_int_array", 3*loc_i, loc_i); //2nd argument = value, 3rd = array index
	*/

	/************************************************* LOOP OVER COMBOS *************************************************/

	//Loop over combos
	for (UInt_t locComboIndex = 0; locComboIndex < Get_NumCombos(); ++locComboIndex) {
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(locComboIndex);

		// Is used to indicate when combos have been cut
		if (dComboWrapper->Get_IsComboCut()) // Is false when tree originally created
			continue; // Combo has been cut previously

		/********************************************** GET PARTICLE INDICES *********************************************/

		//Used for tracking uniqueness when filling histograms, and for determining unused particles

		//Step 0
		const Int_t locBeamID          = dComboBeamWrapper->Get_BeamID();
		const Int_t locPiPlus1TrackID  = dPiPlus1Wrapper->Get_TrackID();
		const Int_t locPiPlus2TrackID  = dPiPlus2Wrapper->Get_TrackID();
		const Int_t locPiMinus1TrackID = dPiMinus1Wrapper->Get_TrackID();
		const Int_t locPiMinus2TrackID = dPiMinus2Wrapper->Get_TrackID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: is kinfit if kinfit performed, else is measured
		// dTargetP4 is target p4
		// Step 0
		const TLorentzVector locBeamP4          = dComboBeamWrapper->Get_P4();
		const TLorentzVector locPiPlus1P4       = dPiPlus1Wrapper->Get_P4();
		const TLorentzVector locPiPlus2P4       = dPiPlus2Wrapper->Get_P4();
		const TLorentzVector locPiMinus1P4      = dPiMinus1Wrapper->Get_P4();
		const TLorentzVector locPiMinus2P4      = dPiMinus2Wrapper->Get_P4();
		const TLorentzVector locMissingProtonP4 = dMissingProtonWrapper->Get_P4();

		// Get Measured P4's:
		// Step 0
		const TLorentzVector locBeamP4_Measured     = dComboBeamWrapper->Get_P4_Measured();
		const TLorentzVector locPiPlus1P4_Measured  = dPiPlus1Wrapper->Get_P4_Measured();
		const TLorentzVector locPiPlus2P4_Measured  = dPiPlus2Wrapper->Get_P4_Measured();
		const TLorentzVector locPiMinus1P4_Measured = dPiMinus1Wrapper->Get_P4_Measured();
		const TLorentzVector locPiMinus2P4_Measured = dPiMinus2Wrapper->Get_P4_Measured();
		const TLorentzVector locMissingProtonP4_Measured
			= locBeamP4_Measured + dTargetP4 - (locPiPlus1P4_Measured + locPiPlus2P4_Measured + locPiMinus1P4_Measured + locPiMinus2P4_Measured);

		/********************************************* GET COMBO RF TIMING INFO *****************************************/

		const TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();
		// Double_t locBunchPeriod = dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
		// Double_t locDeltaT_RF = dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4_Measured, dComboWrapper);
		const Int_t locRelBeamBucket = dAnalysisUtilities.Get_RelativeBeamBucket(Get_RunNumber(), locBeamX4_Measured, dComboWrapper); // 0 for in-time events, non-zero integer for out-of-time photons
		const Int_t locNumOutOfTimeBunchesInTree = 1;  // YOU need to specify this number
			// Number of out-of-time beam bunches in tree (on a single side, so that total number out-of-time bunches accepted is 2 times this number for left + right bunches)
		const Bool_t locSkipNearestOutOfTimeBunch = false;  // True: skip events from nearest out-of-time bunch on either side (recommended).
		if (locSkipNearestOutOfTimeBunch and abs(locRelBeamBucket) == 1) {  // Skip nearest out-of-time bunch: tails of in-time distribution also leak in
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}

		const Int_t    locNumOutOfTimeBunchesToUse = (locSkipNearestOutOfTimeBunch) ? locNumOutOfTimeBunchesInTree - 1 : locNumOutOfTimeBunchesInTree;
		const Double_t locAccidentalScalingFactor  = dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(), locBeamP4.E(), dIsMC);  // Ideal value would be 1, but deviations require added factor, which is different for data and MC.
		// const Double_t locAccidentalScalingFactorError = dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber(), locBeamP4.E());  // Ideal value would be 1, but deviations observed, need added factor.
		const Double_t locHistAccidWeightFactor    = (dSidebandSubtractAcc) ? ((locRelBeamBucket == 0) ? 1 : -locAccidentalScalingFactor / (2 * locNumOutOfTimeBunchesToUse)) : 1;  // Weight by 1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// DO YOUR STUFF HERE
		// kinematic variables of missing proton
		// values from kinematic fit
		const double locMissingMassSquared = locMissingProtonP4.M2();
		// measured values
		const double locMissingMassSquared_Measured = locMissingProtonP4_Measured.M2();

		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		// dAnalyzeCutActions->Perform_Action(); // Must be executed before Execute_Actions()
		if (!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

		// ECAL info
		const size_t locNumUnusedShowers    = dComboWrapper->Get_NumUnusedShowers();
		const double locEnergyUnusedShowers = dComboWrapper->Get_Energy_UnusedShowers();

		//if you manually execute any actions, and it fails a cut, be sure to call:
			//dComboWrapper->Set_IsComboCut(true);

		/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

		/*
		TLorentzVector locMyComboP4(8.0, 7.0, 6.0, 5.0);
		//for arrays below: 2nd argument is value, 3rd is array index
		//NOTE: By filling here, AFTER the cuts above, some indices won't be updated (and will be whatever they were from the last event)
			//So, when you draw the branch, be sure to cut on "IsComboCut" to avoid these.
		dTreeInterface->Fill_Fundamental<Float_t>("my_combo_array", -2 * locComboIndex, locComboIndex);
		dTreeInterface->Fill_TObject<TLorentzVector>("my_p4_array", locMyComboP4, locComboIndex);
		*/

		/**************************************** EXAMPLE: HISTOGRAM BEAM ENERGY *****************************************/

		//Histogram beam energy (if haven't already)
		const double locBeamEnergy = locBeamP4.E();
		if (locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end()) {
			// dHist_BeamEnergy->Fill(locBeamEnergy); // Fills in-time and out-of-time beam photon combos
			dHist_BeamEnergy->Fill(locBeamEnergy, locHistAccidWeightFactor); // Alternate version with accidental subtraction

			locUsedSoFar_BeamEnergy.insert(locBeamID);
		}

		/****************************************** SET FLAT TREE BRANCHES ******************************************/

		// RECOMMENDED: FILL ACCIDENTAL WEIGHT
		dFlatTreeInterface->Fill_Fundamental<Double_t>("AccidWeightFactor", locHistAccidWeightFactor);

		// FILL ANY CUSTOM BRANCHES FIRST!!
		dFlatTreeInterface->Fill_Fundamental<Double_t>("BeamEnergy",                  locBeamEnergy);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("MissingMassSquared",          locMissingMassSquared);
		fillTreeKinematics(locMissingProtonP4, "MissingProton");
		dFlatTreeInterface->Fill_Fundamental<Double_t>("MissingMassSquared_Measured", locMissingMassSquared_Measured);
		fillTreeKinematics(locMissingProtonP4_Measured, "MissingProton", "_Measured");
		dFlatTreeInterface->Fill_Fundamental<Int_t>   ("NmbUnusedShowers",            locNumUnusedShowers);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("EnergyUnusedShowers",         locEnergyUnusedShowers);
		dFlatTreeInterface->Fill_TObject<myTObjString>("ThrownTopology",              myTObjString(locThrownTopology));

		/************************************ EXAMPLE: HISTOGRAM MISSING MASS SQUARED ************************************/

		// adapted from https://github.com/jrstevenjlab/wm_gluex/blob/master/analysis/omega_misspi/selector/DSelector_omega_misspi.C#L481
		// Loop over charged-track hypotheses and find unused track
		// NOTE! this assumes that there is maximum 1 unused track
		bool   locUnusedTrackExists = false;
		size_t nmbUnusedTracks      = 0;
		for (UInt_t locTrackIndex = 0; locTrackIndex < Get_NumChargedHypos(); ++locTrackIndex) {
			// Set branch array indices corresponding to this charged-track hypothesis
			dChargedHypoWrapper->Set_ArrayIndex(locTrackIndex);
			// Make sure it is the unused track
			if ((dChargedHypoWrapper->Get_PID() != Proton)
			    or (dChargedHypoWrapper->Get_TrackID() == locPiPlus1TrackID)  or (dChargedHypoWrapper->Get_TrackID() == locPiPlus2TrackID)
			    or (dChargedHypoWrapper->Get_TrackID() == locPiMinus1TrackID) or (dChargedHypoWrapper->Get_TrackID() == locPiMinus2TrackID)) {
				continue;
			}
			locUnusedTrackExists = true;
			++nmbUnusedTracks;

			//Uniqueness tracking: Build the map of particles used for the unused track
				//For beam: Don't want to group with final-state photons. Instead use "Unknown" PID (not ideal, but it's easy).
			map<Particle_t, set<Int_t> > locUsedThisCombo_UnusedTrack;
			locUsedThisCombo_UnusedTrack[Unknown].insert(locBeamID); //beam
			locUsedThisCombo_UnusedTrack[PiPlus ].insert(locPiPlus1TrackID);
			locUsedThisCombo_UnusedTrack[PiPlus ].insert(locPiPlus2TrackID);
			locUsedThisCombo_UnusedTrack[PiMinus].insert(locPiMinus1TrackID);
			locUsedThisCombo_UnusedTrack[PiMinus].insert(locPiMinus2TrackID);
			locUsedThisCombo_UnusedTrack[Proton ].insert(dChargedHypoWrapper->Get_TrackID());
			//compare to what's been used so far
			if (locUsedSoFar_UnusedTrack.find(locUsedThisCombo_UnusedTrack) == locUsedSoFar_UnusedTrack.end()) {
				// 4-vector of unused track
				const TLorentzVector locUnusedTrackP4 = dChargedHypoWrapper->Get_P4_Measured();
				// fill tree variables for unused track
				fillTreeKinematics(locUnusedTrackP4, "Unused", "", 0);
				fillTreeDelta(locMissingProtonP4,          locUnusedTrackP4, 0, "Unused");
				fillTreeDelta(locMissingProtonP4_Measured, locUnusedTrackP4, 0, "Unused", "_Measured");

				locUsedSoFar_UnusedTrack.insert(locUsedThisCombo_UnusedTrack);
			}
		}
		assert(nmbUnusedTracks <= 1);

		//Uniqueness tracking: Build the map of particles used for the missing mass
			//For beam: Don't want to group with final-state photons. Instead use "Unknown" PID (not ideal, but it's easy).
		map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass;
		locUsedThisCombo_MissingMass[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_MissingMass[PiPlus ].insert(locPiPlus1TrackID);
		locUsedThisCombo_MissingMass[PiPlus ].insert(locPiPlus2TrackID);
		locUsedThisCombo_MissingMass[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_MissingMass[PiMinus].insert(locPiMinus2TrackID);
		//compare to what's been used so far
		if (locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == locUsedSoFar_MissingMass.end()) {
			//unique missing mass combo: histogram it, and register this combo of particles
			dHist_MissingMassSquaredVsBeamEnergy->Fill        (locBeamEnergy, locMissingMassSquared_Measured, locHistAccidWeightFactor);
			dHist_MissingMassSquaredVsBeamEnergySideband->Fill(locBeamEnergy, locMissingMassSquared_Measured, 1 - locHistAccidWeightFactor);

			dFlatTreeInterface->Fill_Fundamental<Int_t>("NmbUnusedTracks", (locUnusedTrackExists) ? 1 : 0);  // indicate whether there was an unused track in the event or not
			// fill tree variables for truth track
			int locNmbTruthTracks = 0;
			for (UInt_t locTrackIndex = 0; locTrackIndex < Get_NumThrown(); ++locTrackIndex) {
				// Set branch array indices corresponding to this charged-track hypothesis
				dThrownWrapper->Set_ArrayIndex(locTrackIndex);
				// Make sure it is the unused track
				if (dThrownWrapper->Get_PID() != Proton) {
					continue;
				}
				const TLorentzVector locMissingProtonP4_Truth = dThrownWrapper->Get_P4_Measured();
				fillTreeKinematics(locMissingProtonP4_Truth, "Truth", "", locNmbTruthTracks);
				fillTreeDelta(locMissingProtonP4,          locMissingProtonP4_Truth, locNmbTruthTracks, "Truth");
				fillTreeDelta(locMissingProtonP4_Measured, locMissingProtonP4_Truth, locNmbTruthTracks, "Truth", "_Measured");

				++locNmbTruthTracks;
			}
			dFlatTreeInterface->Fill_Fundamental<Int_t>("NmbTruthTracks", locNmbTruthTracks);
			// FILL FLAT TREE
			Fill_FlatTree();

			locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);
		}

		//E.g. Cut
		//if((locMissingMassSquared < -0.04) || (locMissingMassSquared > 0.04))
		//{
		//	dComboWrapper->Set_IsComboCut(true);
		//	continue;
		//}

	} // end of combo loop

	//FILL HISTOGRAMS: Num combos / events surviving actions
	Fill_NumCombosSurvivedHists();

	/******************************************* LOOP OVER THROWN DATA (OPTIONAL) ***************************************/
/*
	//Thrown beam: just use directly
	if(dThrownBeam != NULL)
		double locEnergy = dThrownBeam->Get_P4().E();

	//Loop over throwns
	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dThrownWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
*/
	/****************************************** LOOP OVER OTHER ARRAYS (OPTIONAL) ***************************************/
/*
	//Loop over beam particles (note, only those appearing in combos are present)
	for(UInt_t loc_i = 0; loc_i < Get_NumBeam(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dBeamWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}

	//Loop over charged track hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dChargedHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}

	//Loop over neutral particle hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumNeutralHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dNeutralHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
*/

	/************************************ EXAMPLE: FILL CLONE OF TTREE HERE WITH CUTS APPLIED ************************************/
/*
	Bool_t locIsEventCut = true;
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);
		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut())
			continue;
		locIsEventCut = false; // At least one combo succeeded
		break;
	}
	if(!locIsEventCut && dOutputTreeFileName != "")
		Fill_OutputTree();
*/

	return kTRUE;
}

void DSelector_pippippimpimpmiss::Finalize(void)
{
	//Save anything to output here that you do not want to be in the default DSelector output ROOT file.

	//Otherwise, don't do anything else (especially if you are using PROOF).
		//If you are using PROOF, this function is called on each thread,
		//so anything you do will not have the combined information from the various threads.
		//Besides, it is best-practice to do post-processing (e.g. fitting) separately, in case there is a problem.

	//DO YOUR STUFF HERE

	//CALL THIS LAST
	DSelector::Finalize(); //Saves results to the output file
}
