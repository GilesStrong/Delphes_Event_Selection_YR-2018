/*
 * main.cc
 *
 *  Created on: 8 Apr 2016
 *      Author: giles giles.strong@outlook.com
 */
//Local
#include "main.hh"

bool debug = false;

TLorentzVector getHiggs2Taus(TClonesArray* mpt, TLorentzVector t_0, TLorentzVector t_1) {
	/*Returns 4-vector of Higgs->tau tau*/
	TLorentzVector higgs, mPT;
	mPT.SetPtEtaPhiM((MissingET*)(mpt->At(0))->MET, 0.0, (MissingET*)(mpt->At(0))->Phi, 0.0); //TODO Check this
	higgs = t_0 + t_1 + mPT;
	return higgs;
}

bool selectBJets(TClonesArray* jets, std::vector<int>* bJets, int* bJet_0, int* bJet_1) {
	/*Checks is a pair of b-jets exists, returning true if so and pointing bJet_0 and bJet_1 to
	selected jets. Selects pair of jets invariant mass closest to 125 GeV*/
	if (bJets->size() == 2) { //Only two b jets found
		*bJet_0 = (*bJets)[0];
		*bJet_1 = (*bJets)[1];
		if ((Jet*)(jets->at(bJet_0))->PT() < (Jet*)(jets->at(bJet_1))->PT()) {
			*bJet_1 = (*bJets)[0];
			*bJet_0 = (*bJets)[1];
		}
		return true;
	} else if (bJets->size() > 2) { //More than two b jets: select pair with invariant mass closest to 125 GeV
		double deltaMin = -1;
		double delta;
		TLorentzVector jet_i, jet_j, jet_combined;
		int iMin, jMin;
		for (int i : *bJets) {
			jet_i = jets->at(bJet_1)->P4;
			for (int j : *bJets) {
				if (i == j) continue;
				jet_j = jets->at(bJet_1)->P4;
				jet_combined = jet_i + jet_j;
				delta = std::abs(125-jet_combined.M());
				if (deltaMin > delta || deltaMin < 0) {
					deltaMin = delta;
					iMin = i;
					jMin = j;
				}
			}
		}
		*bJet_0 = iMin;
		*bJet_1 = jMin;
		if (jets->at(bJet_0)->PT() < jets->at(bJet_1)->PT()) {
			*bJet_1 = iMin;
			*bJet_0 = jMin;
		}
		return true;
	} else { //Less than two b jets found
		return false;
	}
}

void makeDirs(std::string outputName) {
	/*Makes directory structure for saving outputs*/
	std::vector<std::string> dirs;
	dirs.push_back("../outputs");
	dirs.push_back("../outputs/" + outputName);
	for (std::string dir : dirs) {
		system((char*)("mkdir -p " + dir).c_str());
	}
}

inline TLorentzVector getHiggs2Bs(TLorentzVector b_0, TLorentzVector b_1) {
	/*Returns 4-vector of Higgs->b b*/
	return b_0 + b_1;
}

inline TLorentzVector getDiHiggs(TLorentzVector higgs_0, TLorentzVector higgs_1) {
	/*Returns 4-vector of di_Higgs*/
	return higgs_0 + higgs_1;
}

inline double getMT(double pT, double mPT, double dphi) {
	return sqrt(2*pT*mPT*(1-cos(dphi)));
}

void showHelp() {
	/*Show help for input arguments*/
	std::cout << "-i : input mask\n";
	std::cout << "-o : output name\n";
	std::cout << "-t : use MC truth cut [0/1], default 0\n";
	std::cout << "-s : run event selection [0/1], default 1\n";
	std::cout << "-d : run in debug mode [0/1], default 0\n";
	std::cout << "-m : output information for MVA selection [0/1], default 0\n";
}

std::map<std::string, std::string> getOptions(int argc, char* argv[]) {
	/*Interpret input arguments*/
	std::map<std::string, std::string> options;
	options.insert(std::make_pair("-i", "")); //Input mask
	options.insert(std::make_pair("-o", "")); //Output name
	options.insert(std::make_pair("-t", "0")); //MC truth cut
	options.insert(std::make_pair("-s", "1")); //Event selection
	options.insert(std::make_pair("-d", "0")); //Debug mode
	options.insert(std::make_pair("-m", "0")); //Output information for MVA selection
	if (argc >= 2) {
		std::string option(argv[1]);
		if (option == "-h" || option == "--help") {
			showHelp();
			options.clear();
			return options;
		}
	}
	for (int i = 1; i < argc; i = i+2) {
		std::string option(argv[i]);
		std::string argument(argv[i+1]);
		if (option == "-h" || option == "--help" || argument == "-h" || argument == "--help") {
			showHelp();
			options.clear();
			return options;
		}
		options[option] = argument;
	}
	if (options["-i"] == "" || options["-o"] == "") {
		showHelp();
		options.clear();
		return options;
	}
	if (options["-d"] == "1") {
		debug = true;
		std::cout << "Running in debug mode\n";
	}
	return options;
}

TMatrixD decomposeVector(TLorentzVector in) {
	TMatrixD out(3, 3);
	out(0, 0) = in.Px()*in.Px();
	out(0, 1) = in.Px()*in.Py();
	out(0, 2) = in.Px()*in.Pz();
	out(1, 0) = in.Py()*in.Px();
	out(1, 1) = in.Py()*in.Py();
	out(1, 2) = in.Py()*in.Pz();
	out(2, 0) = in.Pz()*in.Px();
	out(2, 1) = in.Pz()*in.Py();
	out(2, 2) = in.Pz()*in.Pz();
	return out;
}

void appendSphericity(TMatrixD* mat, double* div, TLorentzVector mom) {
	TMatrixD decomp = decomposeVector(mom);
	*mat += decomp;
	*div += pow(mom.P(), 2);
}

void appendSpherocity(TMatrixD* mat, double* div, TLorentzVector mom) {
	TMatrixD decomp = decomposeVector(mom);
	decomp *= 1/std::abs(mom.P());
	*mat += decomp;
	*div += std::abs(mom.P());
}

std::vector<double> getEigenValues(TMatrixD in) {
	/*Return vector of sorted, nomalised eigenvalues of parssed matrix*/
	TMatrixD eigenMatrix = TMatrixDEigen(in).GetEigenValues();
	std::vector<double> eigenValues(3);
	eigenValues[0] = eigenMatrix(0, 0);
	eigenValues[1] = eigenMatrix(1, 1);
	eigenValues[2] = eigenMatrix(2, 2);
	std::sort(eigenValues.begin(), eigenValues.end(), std::greater<double>());
	double sum = 0;
	for (double n : eigenValues)
		sum += n;
	std::for_each(eigenValues.begin(), eigenValues.end(), [sum](double i) { return i/sum; });
	return eigenValues;
}

void getEventShapes(std::vector<double> sphericityV, std::vector<double> spherocityV,
		double* sphericity, double* spherocity,
		double* aplanarity, double* aplanority,
		double* upsilon, double* dShape) {
	*sphericity = (3/2)*(sphericityV[1]+sphericityV[2]);
	*spherocity = (3/2)*(spherocityV[1]+spherocityV[2]);
	*aplanarity = 3*sphericityV[2]/2;
	*aplanority = 3*spherocityV[2]/2;
	*upsilon = sqrt(3.0)*(sphericityV[1]-sphericityV[2])/2;
	*dShape = 27*spherocityV[0]*spherocityV[1]*spherocityV[2];
}

void getGlobalEventInfo(std::string input, Long64_t cEvent,
		double*  hT, double*  sT, double* centrality, double* eVis,
		int* nJets, int* nBJets, int* nTauJets,
		double* minJetPT, double* meanJetPT, double* maxJetPT,
		double* minJetMass, double* meanJetMass, double* maxJetMass,
		double* minJetEta, double* meanJetEta, double* maxJetEta,
		double* sphericity, double* spherocity,
		double* aplanarity, double* aplanority,
		double* upsilon, double* dShape) {
	/*Fills referenced variables with global event information*/
	if (debug) std::cout << "Getting global event info\n";
	//Load event info____________________________
	TChain *chain = new TChain("Delphes");
	chain->Add(input.c_str());
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchMuon = treeReader->UseBranch("MuonLoose");
	TClonesArray *branchJet = treeReader->UseBranch("Jet");
	TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
	treeReader->ReadEntry(cEvent);
	if (debug) std::cout << "Loaded info\n";
	//___________________________________________
	//Reset variables____________________________
	*hT = 0;
	*sT = 0;
	*centrality = 0;
	*eVis = 0;
	*nJets = 0;
	*nBJets = 0;
	*nTauJets = 0;
	*minJetPT = -1;
	*meanJetPT = 0;
	*maxJetPT = -1;
	*minJetMass = -1;
	*meanJetMass = 0;
	*maxJetMass = -1;
	*minJetEta = -1;
	*meanJetEta = 0;
	*maxJetEta = -1;
	*sphericity = 0;
	*spherocity = 0;
	*aplanarity = 0;
	*aplanority = 0;
	*upsilon = 0;
	*dShape = 0;
	//___________________________________________
	//Initialise holders_________________________
	Electron* electron;
	Muon* muon;
	Jet* jet;
	MissingET* mPT;
	TMatrixD sphericityT(3, 3), spherocityT(3, 3);
	double sphericityD = 0, spherocityD = 0;
	//___________________________________________
	//Loop through objects_____________________
	for (int i = 0; i < branchElectron->GetEntriesFast(); ++i) { //Loop over all electrons in event
		electron = (Electron*)branchElectron->At(i);
		*sT += electron->PT;
		*eVis += electron->P4().E();
		*centrality += electron->PT;
		appendSphericity(&sphericityT, &sphericityD, electron->P4());
		appendSpherocity(&spherocityT, &spherocityD, electron->P4());
	}
	for (int i = 0; i < branchMuon->GetEntriesFast(); ++i) { //Loop over all muons in event
		muon = (Muon*)branchMuon->At(i);
		*sT += muon->PT;
		*eVis += muon->P4().E();
		*centrality += muon->PT;
		appendSphericity(&sphericityT, &sphericityD, muon->P4());
		appendSpherocity(&spherocityT, &spherocityD, muon->P4());
	}
	mPT = (MissingET*)branchMissingET->At(0);
	*sT += mPT->MET;
	for (int i = 0; i < branchJet->GetEntriesFast(); ++i) { //Loop over all jets in event
		jet = (Jet*)branchJet->At(i);
		*hT += sqrt(pow(jet->Mass, 2)+pow(jet->PT, 2));
		*sT += sqrt(pow(jet->Mass, 2)+pow(jet->PT, 2));
		*centrality += jet->PT;
		*eVis += jet->P4().E();
		*nJets += 1;
		if (jet->TauTag) *nTauJets += 1;
		if (jet->BTag) *nBJets += 1;
		if (*minJetPT == -1 | jet->PT < *minJetPT) *minJetPT = jet->PT;
		*meanJetPT += jet->PT;
		if (jet->PT > *maxJetPT) *maxJetPT = jet->PT;
		if (*minJetMass == -1 | jet->Mass < *minJetMass) *minJetMass = jet->Mass;
		*meanJetMass += jet->Mass;
		if (jet->Mass > *maxJetMass) *maxJetMass = jet->Mass;
		if (*minJetEta == -1 | std::abs(jet->Eta) < *minJetEta) *minJetEta = std::abs(jet->Eta);
		*meanJetEta += jet->Eta;
		if (std::abs(jet->Eta) > *maxJetEta) *maxJetEta = std::abs(jet->Eta);
		appendSphericity(&sphericityT, &sphericityD, jet->P4());
		appendSpherocity(&spherocityT, &spherocityD, jet->P4());
	}
	//___________________________________________
	//Finalise variabales________________________
	*centrality /= *eVis;
	*meanJetPT /= *nJets;
	*meanJetMass /= *nJets;
	*meanJetEta /= *nJets;
	sphericityT *= 1/sphericityD;
	spherocityT *= 1/spherocityD;
	//___________________________________________
	//Calculate event shapes_____________________
	if (debug) std::cout << "Calculating global event shapes\n";
	std::vector<double> sphericityV = getEigenValues(sphericityT);
	std::vector<double> spherocityV = getEigenValues(spherocityT);
	getEventShapes(sphericityV, spherocityV,
			sphericity, spherocity,
			aplanarity, aplanority,
			upsilon, dShape);
	//___________________________________________
	chain->Delete();
	delete treeReader;
}

void getPrimaryEventShapes(TLorentzVector v_tau_0, TLorentzVector v_tau_1, TLorentzVector v_bJet_0, TLorentzVector v_bJet_1,
		double* sphericity, double* spherocity,
		double* aplanarity, double* aplanority,
		double* upsilon, double* dShape) {
	/*Sets values of referenced event-shape variables for final-states*/
	if (debug) std::cout << "Getting primary event shapes\n";
	TMatrixD sphericityT(3, 3), spherocityT(3, 3);
	double sphericityD = 0, spherocityD = 0;
	//Populate tensors___________________________
	appendSphericity(&sphericityT, &sphericityD, v_tau_0);
	appendSpherocity(&spherocityT, &spherocityD, v_tau_0);
	appendSphericity(&sphericityT, &sphericityD, v_tau_1);
	appendSpherocity(&spherocityT, &spherocityD, v_tau_1);
	appendSphericity(&sphericityT, &sphericityD, v_bJet_0);
	appendSpherocity(&spherocityT, &spherocityD, v_bJet_0);
	appendSphericity(&sphericityT, &sphericityD, v_bJet_1);
	appendSpherocity(&spherocityT, &spherocityD, v_bJet_1);
	sphericityT *= 1/sphericityD;
	spherocityT *= 1/spherocityD;
	//___________________________________________
	//Calculate event shapes_____________________
	if (debug) std::cout << "Calculating primary event shapes\n";
	std::vector<double> sphericityV = getEigenValues(sphericityT);
	std::vector<double> spherocityV = getEigenValues(spherocityT);
	getEventShapes(sphericityV, spherocityV,
			sphericity, spherocity,
			aplanarity, aplanority,
			upsilon, dShape);
	//___________________________________________
}

int main(int argc, char *argv[]) { //input, output, N events, truth
	std::map<std::string, std::string> options = getOptions(argc, argv);
	if (options.size() == 0) {
		return 1;
	}
	std::string outputName(options["-o"]);
	makeDirs(outputName);
	//ROOT settings______________________________
	gSystem->Load("libDelphes.so");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetPadGridX(kFALSE);
	gStyle->SetPadGridY(kFALSE);
	//___________________________________________
	//Initialise variables_______________________
	std::cout << "Initialising variables\n";
	int lepton_0, lepton_1, tau_0, tau_1, bJet_0, bJet_1;
	//Low-level variables________________________
	double t_0_pT, t_0_eta, t_0_phi, t_0_mass; //Tau 0 variables
	double t_1_pT, t_1_eta, t_1_phi, t_1_mass; //Tau 1 variables
	double b_0_pT, b_0_eta, b_0_phi, b_0_mass; //b-jet 0 variables
	double b_1_pT, b_1_eta, b_1_phi, b_1_mass; //b-jet 1 variables
	double mPT_pT, mPT_phi; //Missing ET variables
	//___________________________________________
	//Reconstructed variables____________________
	double h_tt_pT, h_tt_eta, h_tt_phi, h_tt_mass; //Higgs 0 variables
	double h_bb_pT, h_bb_eta, h_bb_phi, h_bb_mass; //Higgs 1 variables
	double diH_pT, diH_eta, diH_phi, diH_mass; //di-Higgs variables
	//___________________________________________
	//Global event variables_____________________
	double hT, sT, centrality, eVis; //Global kinematics
	int nJets, nBJets, nTauJets; //Jet multiplicities
	double minJetPT, meanJetPT, maxJetPT; //Global jet pTs
	double minJetMass, meanJetMass, maxJetMass; //Global jet masses
	double minJetEta, meanJetEta, maxJetEta; //Global jet etas
	double sphericityA, spherocityA, aplanarityA, aplanorityA, upsilonA, dShapeA; //Event shapes for all objects
	double sphericityP, spherocityP, aplanarityP, aplanorityP, upsilonP, dShapeP; //Event shapes for primary objects
	//___________________________________________
	//Generator-level variables for regression and cuts
	double gen_t_0_pT, gen_t_0_eta, gen_t_0_phi, gen_t_0_E; //Tau 0 variables
	double gen_t_1_pT, gen_t_1_eta, gen_t_1_phi, gen_t_1_E; //Tau 1 variables
	double gen_b_0_pT, gen_b_0_eta, gen_b_0_phi, gen_b_0_E; //b-jet 0 variables
	double gen_b_1_pT, gen_b_1_eta, gen_b_1_phi, gen_b_1_E; //b-jet 1 variables
	double gen_diH_pT, gen_diH_eta, gen_diH_phi, gen_diH_E, gen_diH_mass; //diHiggs variables
	double gen_h_bb_pT, gen_h_bb_eta, gen_h_bb_phi, gen_h_bb_E; //Higgs->bb variables
	double gen_h_tt_pT, gen_h_tt_eta, gen_h_tt_phi, gen_h_tt_E; //Higgs->tau tau variables
	bool gen_mctMatch; //MC truth match
	//___________________________________________
	double weight; //Event weight
	int nElectrons = 0, nMuons = 0;
	bool eventAccepted = false;
	TTree* e_tau_b_b = new TTree("e_tau_b_b", "e #tau b #bar{b}");
	e_tau_b_b->Branch("t_0_pT", &t_0_pT);
	e_tau_b_b->Branch("t_0_eta", &t_0_eta);
	e_tau_b_b->Branch("t_0_phi", &t_0_phi);
	e_tau_b_b->Branch("t_0_mass", &t_0_mass);
	e_tau_b_b->Branch("t_1_pT", &t_1_pT);
	e_tau_b_b->Branch("t_1_eta", &t_1_eta);
	e_tau_b_b->Branch("t_1_phi", &t_1_phi);
	e_tau_b_b->Branch("t_1_mass", &t_1_mass);
	e_tau_b_b->Branch("b_0_pT", &b_0_pT);
	e_tau_b_b->Branch("b_0_eta", &b_0_eta);
	e_tau_b_b->Branch("b_0_phi", &b_0_phi);
	e_tau_b_b->Branch("b_0_mass", &b_0_mass);
	e_tau_b_b->Branch("b_1_pT", &b_1_pT);
	e_tau_b_b->Branch("b_1_eta", &b_1_eta);
	e_tau_b_b->Branch("b_1_phi", &b_1_phi);
	e_tau_b_b->Branch("b_1_mass", &b_1_mass);
	e_tau_b_b->Branch("mPT_pT", &mPT_pT);
	e_tau_b_b->Branch("mPT_phi", &mPT_phi);
	e_tau_b_b->Branch("h_tt_pT", &h_tt_pT);
	e_tau_b_b->Branch("h_tt_eta", &h_tt_eta);
	e_tau_b_b->Branch("h_tt_phi", &h_tt_phi);
	e_tau_b_b->Branch("h_tt_mass", &h_tt_mass);
	e_tau_b_b->Branch("h_bb_pT", &h_bb_pT);
	e_tau_b_b->Branch("h_bb_eta", &h_bb_eta);
	e_tau_b_b->Branch("h_bb_phi", &h_bb_phi);
	e_tau_b_b->Branch("h_bb_mass", &h_bb_mass);
	e_tau_b_b->Branch("diH_pT", &diH_pT);
	e_tau_b_b->Branch("diH_eta", &diH_eta);
	e_tau_b_b->Branch("diH_phi", &diH_phi);
	e_tau_b_b->Branch("diH_mass", &diH_mass);
	e_tau_b_b->Branch("hT", &hT);
	e_tau_b_b->Branch("sT", &sT);
	e_tau_b_b->Branch("centrality", &centrality);
	e_tau_b_b->Branch("eVis", &eVis);
	e_tau_b_b->Branch("nJets", &nJets);
	e_tau_b_b->Branch("nBJets", &nBJets);
	e_tau_b_b->Branch("nTauJets", &nTauJets);
	e_tau_b_b->Branch("minJetPT", &minJetPT);
	e_tau_b_b->Branch("meanJetPT", &meanJetPT);
	e_tau_b_b->Branch("maxJetPT", &maxJetPT);
	e_tau_b_b->Branch("minJetMass", &minJetMass);
	e_tau_b_b->Branch("meanJetMass", &meanJetMass);
	e_tau_b_b->Branch("maxJetMass", &maxJetMass);
	e_tau_b_b->Branch("minJetEta", &minJetEta);
	e_tau_b_b->Branch("meanJetEta", &meanJetEta);
	e_tau_b_b->Branch("maxJetEta", &maxJetEta);
	e_tau_b_b->Branch("sphericityA", &sphericityA);
	e_tau_b_b->Branch("spherocityA", &spherocityA);
	e_tau_b_b->Branch("aplanarityA", &aplanarityA);
	e_tau_b_b->Branch("aplanorityA", &aplanorityA);
	e_tau_b_b->Branch("upsilonA", &upsilonA);
	e_tau_b_b->Branch("dShapeA", &dShapeA);
	e_tau_b_b->Branch("sphericityP", &sphericityP);
	e_tau_b_b->Branch("spherocityP", &spherocityP);
	e_tau_b_b->Branch("aplanarityP", &aplanarityP);
	e_tau_b_b->Branch("aplanorityP", &aplanorityP);
	e_tau_b_b->Branch("upsilonP", &upsilonP);
	e_tau_b_b->Branch("dShapeP", &dShapeP);
	e_tau_b_b->Branch("gen_t_0_pT", &gen_t_0_pT);
	e_tau_b_b->Branch("gen_t_0_eta", &gen_t_0_eta);
	e_tau_b_b->Branch("gen_t_0_phi", &gen_t_0_phi);
	e_tau_b_b->Branch("gen_t_0_E", &gen_t_0_E);
	e_tau_b_b->Branch("gen_t_1_pT", &gen_t_1_pT);
	e_tau_b_b->Branch("gen_t_1_eta", &gen_t_1_eta);
	e_tau_b_b->Branch("gen_t_1_phi", &gen_t_1_phi);
	e_tau_b_b->Branch("gen_t_1_E", &gen_t_1_E);
	e_tau_b_b->Branch("gen_b_0_pT", &gen_b_0_pT);
	e_tau_b_b->Branch("gen_b_0_eta", &gen_b_0_eta);
	e_tau_b_b->Branch("gen_b_0_phi", &gen_b_0_phi);
	e_tau_b_b->Branch("gen_b_0_E", &gen_b_0_E);
	e_tau_b_b->Branch("gen_b_1_pT", &gen_b_1_pT);
	e_tau_b_b->Branch("gen_b_1_eta", &gen_b_1_eta);
	e_tau_b_b->Branch("gen_b_1_phi", &gen_b_1_phi);
	e_tau_b_b->Branch("gen_b_1_E", &gen_b_1_E);
	e_tau_b_b->Branch("gen_diH_pT", &gen_diH_pT);
	e_tau_b_b->Branch("gen_diH_eta", &gen_diH_eta);
	e_tau_b_b->Branch("gen_diH_phi", &gen_diH_phi);
	e_tau_b_b->Branch("gen_diH_E", &gen_diH_E);
	e_tau_b_b->Branch("gen_diH_mass", &gen_diH_mass);
	e_tau_b_b->Branch("gen_h_bb_pT", &gen_h_bb_pT);
	e_tau_b_b->Branch("gen_h_bb_eta", &gen_h_bb_eta);
	e_tau_b_b->Branch("gen_h_bb_phi", &gen_h_bb_phi);
	e_tau_b_b->Branch("gen_h_bb_E", &gen_h_bb_E);
	e_tau_b_b->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	e_tau_b_b->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	e_tau_b_b->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	e_tau_b_b->Branch("gen_h_tt_E", &gen_h_tt_E);
	e_tau_b_b->Branch("gen_mctMatch", &gen_mctMatch);
	e_tau_b_b->Branch("gen_weight", &weight);
	TTree* mu_tau_b_b = new TTree("mu_tau_b_b", "#mu #tau_{h} b #bar{b}");
	mu_tau_b_b->Branch("t_0_pT", &t_0_pT);
	mu_tau_b_b->Branch("t_0_eta", &t_0_eta);
	mu_tau_b_b->Branch("t_0_phi", &t_0_phi);
	mu_tau_b_b->Branch("t_0_mass", &t_0_mass);
	mu_tau_b_b->Branch("t_1_pT", &t_1_pT);
	mu_tau_b_b->Branch("t_1_eta", &t_1_eta);
	mu_tau_b_b->Branch("t_1_phi", &t_1_phi);
	mu_tau_b_b->Branch("t_1_mass", &t_1_mass);
	mu_tau_b_b->Branch("b_0_pT", &b_0_pT);
	mu_tau_b_b->Branch("b_0_eta", &b_0_eta);
	mu_tau_b_b->Branch("b_0_phi", &b_0_phi);
	mu_tau_b_b->Branch("b_0_mass", &b_0_mass);
	mu_tau_b_b->Branch("b_1_pT", &b_1_pT);
	mu_tau_b_b->Branch("b_1_eta", &b_1_eta);
	mu_tau_b_b->Branch("b_1_phi", &b_1_phi);
	mu_tau_b_b->Branch("b_1_mass", &b_1_mass);
	mu_tau_b_b->Branch("mPT_pT", &mPT_pT);
	mu_tau_b_b->Branch("mPT_phi", &mPT_phi);
	mu_tau_b_b->Branch("h_tt_pT", &h_tt_pT);
	mu_tau_b_b->Branch("h_tt_eta", &h_tt_eta);
	mu_tau_b_b->Branch("h_tt_phi", &h_tt_phi);
	mu_tau_b_b->Branch("h_tt_mass", &h_tt_mass);
	mu_tau_b_b->Branch("h_bb_pT", &h_bb_pT);
	mu_tau_b_b->Branch("h_bb_eta", &h_bb_eta);
	mu_tau_b_b->Branch("h_bb_phi", &h_bb_phi);
	mu_tau_b_b->Branch("h_bb_mass", &h_bb_mass);
	mu_tau_b_b->Branch("diH_pT", &diH_pT);
	mu_tau_b_b->Branch("diH_eta", &diH_eta);
	mu_tau_b_b->Branch("diH_phi", &diH_phi);
	mu_tau_b_b->Branch("diH_mass", &diH_mass);
	mu_tau_b_b->Branch("hT", &hT);
	mu_tau_b_b->Branch("sT", &sT);
	mu_tau_b_b->Branch("centrality", &centrality);
	mu_tau_b_b->Branch("eVis", &eVis);
	mu_tau_b_b->Branch("nJets", &nJets);
	mu_tau_b_b->Branch("nBJets", &nBJets);
	mu_tau_b_b->Branch("nTauJets", &nTauJets);
	mu_tau_b_b->Branch("minJetPT", &minJetPT);
	mu_tau_b_b->Branch("meanJetPT", &meanJetPT);
	mu_tau_b_b->Branch("maxJetPT", &maxJetPT);
	mu_tau_b_b->Branch("minJetMass", &minJetMass);
	mu_tau_b_b->Branch("meanJetMass", &meanJetMass);
	mu_tau_b_b->Branch("maxJetMass", &maxJetMass);
	mu_tau_b_b->Branch("minJetEta", &minJetEta);
	mu_tau_b_b->Branch("meanJetEta", &meanJetEta);
	mu_tau_b_b->Branch("maxJetEta", &maxJetEta);
	mu_tau_b_b->Branch("sphericityA", &sphericityA);
	mu_tau_b_b->Branch("spherocityA", &spherocityA);
	mu_tau_b_b->Branch("aplanarityA", &aplanarityA);
	mu_tau_b_b->Branch("aplanorityA", &aplanorityA);
	mu_tau_b_b->Branch("upsilonA", &upsilonA);
	mu_tau_b_b->Branch("dShapeA", &dShapeA);
	mu_tau_b_b->Branch("sphericityP", &sphericityP);
	mu_tau_b_b->Branch("spherocityP", &spherocityP);
	mu_tau_b_b->Branch("aplanarityP", &aplanarityP);
	mu_tau_b_b->Branch("aplanorityP", &aplanorityP);
	mu_tau_b_b->Branch("upsilonP", &upsilonP);
	mu_tau_b_b->Branch("dShapeP", &dShapeP);
	mu_tau_b_b->Branch("gen_t_0_pT", &gen_t_0_pT);
	mu_tau_b_b->Branch("gen_t_0_eta", &gen_t_0_eta);
	mu_tau_b_b->Branch("gen_t_0_phi", &gen_t_0_phi);
	mu_tau_b_b->Branch("gen_t_0_E", &gen_t_0_E);
	mu_tau_b_b->Branch("gen_t_1_pT", &gen_t_1_pT);
	mu_tau_b_b->Branch("gen_t_1_eta", &gen_t_1_eta);
	mu_tau_b_b->Branch("gen_t_1_phi", &gen_t_1_phi);
	mu_tau_b_b->Branch("gen_t_1_E", &gen_t_1_E);
	mu_tau_b_b->Branch("gen_b_0_pT", &gen_b_0_pT);
	mu_tau_b_b->Branch("gen_b_0_eta", &gen_b_0_eta);
	mu_tau_b_b->Branch("gen_b_0_phi", &gen_b_0_phi);
	mu_tau_b_b->Branch("gen_b_0_E", &gen_b_0_E);
	mu_tau_b_b->Branch("gen_b_1_pT", &gen_b_1_pT);
	mu_tau_b_b->Branch("gen_b_1_eta", &gen_b_1_eta);
	mu_tau_b_b->Branch("gen_b_1_phi", &gen_b_1_phi);
	mu_tau_b_b->Branch("gen_b_1_E", &gen_b_1_E);
	mu_tau_b_b->Branch("gen_diH_pT", &gen_diH_pT);
	mu_tau_b_b->Branch("gen_diH_eta", &gen_diH_eta);
	mu_tau_b_b->Branch("gen_diH_phi", &gen_diH_phi);
	mu_tau_b_b->Branch("gen_diH_E", &gen_diH_E);
	mu_tau_b_b->Branch("gen_diH_mass", &gen_diH_mass);
	mu_tau_b_b->Branch("gen_h_bb_pT", &gen_h_bb_pT);
	mu_tau_b_b->Branch("gen_h_bb_eta", &gen_h_bb_eta);
	mu_tau_b_b->Branch("gen_h_bb_phi", &gen_h_bb_phi);
	mu_tau_b_b->Branch("gen_h_bb_E", &gen_h_bb_E);
	mu_tau_b_b->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	mu_tau_b_b->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	mu_tau_b_b->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	mu_tau_b_b->Branch("gen_h_tt_E", &gen_h_tt_E);
	mu_tau_b_b->Branch("gen_mctMatch", &gen_mctMatch);
	mu_tau_b_b->Branch("gen_weight", &weight);
	TTree* tau_tau_b_b = new TTree("tau_tau_b_b", "#tau_{h} #tau_{h} b #bar{b}");
	tau_tau_b_b->Branch("t_0_pT", &t_0_pT);
	tau_tau_b_b->Branch("t_0_eta", &t_0_eta);
	tau_tau_b_b->Branch("t_0_phi", &t_0_phi);
	tau_tau_b_b->Branch("t_0_mass", &t_0_mass);
	tau_tau_b_b->Branch("t_1_pT", &t_1_pT);
	tau_tau_b_b->Branch("t_1_eta", &t_1_eta);
	tau_tau_b_b->Branch("t_1_phi", &t_1_phi);
	tau_tau_b_b->Branch("t_1_mass", &t_1_mass);
	tau_tau_b_b->Branch("b_0_pT", &b_0_pT);
	tau_tau_b_b->Branch("b_0_eta", &b_0_eta);
	tau_tau_b_b->Branch("b_0_phi", &b_0_phi);
	tau_tau_b_b->Branch("b_0_mass", &b_0_mass);
	tau_tau_b_b->Branch("b_1_pT", &b_1_pT);
	tau_tau_b_b->Branch("b_1_eta", &b_1_eta);
	tau_tau_b_b->Branch("b_1_phi", &b_1_phi);
	tau_tau_b_b->Branch("b_1_mass", &b_1_mass);
	tau_tau_b_b->Branch("mPT_pT", &mPT_pT);
	tau_tau_b_b->Branch("mPT_phi", &mPT_phi);
	tau_tau_b_b->Branch("h_tt_pT", &h_tt_pT);
	tau_tau_b_b->Branch("h_tt_eta", &h_tt_eta);
	tau_tau_b_b->Branch("h_tt_phi", &h_tt_phi);
	tau_tau_b_b->Branch("h_tt_mass", &h_tt_mass);
	tau_tau_b_b->Branch("h_bb_pT", &h_bb_pT);
	tau_tau_b_b->Branch("h_bb_eta", &h_bb_eta);
	tau_tau_b_b->Branch("h_bb_phi", &h_bb_phi);
	tau_tau_b_b->Branch("h_bb_mass", &h_bb_mass);
	tau_tau_b_b->Branch("diH_pT", &diH_pT);
	tau_tau_b_b->Branch("diH_eta", &diH_eta);
	tau_tau_b_b->Branch("diH_phi", &diH_phi);
	tau_tau_b_b->Branch("diH_mass", &diH_mass);
	tau_tau_b_b->Branch("hT", &hT);
	tau_tau_b_b->Branch("sT", &sT);
	tau_tau_b_b->Branch("centrality", &centrality);
	tau_tau_b_b->Branch("eVis", &eVis);
	tau_tau_b_b->Branch("nJets", &nJets);
	tau_tau_b_b->Branch("nBJets", &nBJets);
	tau_tau_b_b->Branch("nTauJets", &nTauJets);
	tau_tau_b_b->Branch("minJetPT", &minJetPT);
	tau_tau_b_b->Branch("meanJetPT", &meanJetPT);
	tau_tau_b_b->Branch("maxJetPT", &maxJetPT);
	tau_tau_b_b->Branch("minJetMass", &minJetMass);
	tau_tau_b_b->Branch("meanJetMass", &meanJetMass);
	tau_tau_b_b->Branch("maxJetMass", &maxJetMass);
	tau_tau_b_b->Branch("minJetEta", &minJetEta);
	tau_tau_b_b->Branch("meanJetEta", &meanJetEta);
	tau_tau_b_b->Branch("maxJetEta", &maxJetEta);
	tau_tau_b_b->Branch("sphericityA", &sphericityA);
	tau_tau_b_b->Branch("spherocityA", &spherocityA);
	tau_tau_b_b->Branch("aplanarityA", &aplanarityA);
	tau_tau_b_b->Branch("aplanorityA", &aplanorityA);
	tau_tau_b_b->Branch("upsilonA", &upsilonA);
	tau_tau_b_b->Branch("dShapeA", &dShapeA);
	tau_tau_b_b->Branch("sphericityP", &sphericityP);
	tau_tau_b_b->Branch("spherocityP", &spherocityP);
	tau_tau_b_b->Branch("aplanarityP", &aplanarityP);
	tau_tau_b_b->Branch("aplanorityP", &aplanorityP);
	tau_tau_b_b->Branch("upsilonP", &upsilonP);
	tau_tau_b_b->Branch("dShapeP", &dShapeP);
	tau_tau_b_b->Branch("gen_t_0_pT", &gen_t_0_pT);
	tau_tau_b_b->Branch("gen_t_0_eta", &gen_t_0_eta);
	tau_tau_b_b->Branch("gen_t_0_phi", &gen_t_0_phi);
	tau_tau_b_b->Branch("gen_t_0_E", &gen_t_0_E);
	tau_tau_b_b->Branch("gen_t_1_pT", &gen_t_1_pT);
	tau_tau_b_b->Branch("gen_t_1_eta", &gen_t_1_eta);
	tau_tau_b_b->Branch("gen_t_1_phi", &gen_t_1_phi);
	tau_tau_b_b->Branch("gen_t_1_E", &gen_t_1_E);
	tau_tau_b_b->Branch("gen_b_0_pT", &gen_b_0_pT);
	tau_tau_b_b->Branch("gen_b_0_eta", &gen_b_0_eta);
	tau_tau_b_b->Branch("gen_b_0_phi", &gen_b_0_phi);
	tau_tau_b_b->Branch("gen_b_0_E", &gen_b_0_E);
	tau_tau_b_b->Branch("gen_b_1_pT", &gen_b_1_pT);
	tau_tau_b_b->Branch("gen_b_1_eta", &gen_b_1_eta);
	tau_tau_b_b->Branch("gen_b_1_phi", &gen_b_1_phi);
	tau_tau_b_b->Branch("gen_b_1_E", &gen_b_1_E);
	tau_tau_b_b->Branch("gen_diH_pT", &gen_diH_pT);
	tau_tau_b_b->Branch("gen_diH_eta", &gen_diH_eta);
	tau_tau_b_b->Branch("gen_diH_phi", &gen_diH_phi);
	tau_tau_b_b->Branch("gen_diH_E", &gen_diH_E);
	tau_tau_b_b->Branch("gen_diH_mass", &gen_diH_mass);
	tau_tau_b_b->Branch("gen_h_bb_pT", &gen_h_bb_pT);
	tau_tau_b_b->Branch("gen_h_bb_eta", &gen_h_bb_eta);
	tau_tau_b_b->Branch("gen_h_bb_phi", &gen_h_bb_phi);
	tau_tau_b_b->Branch("gen_h_bb_E", &gen_h_bb_E);
	tau_tau_b_b->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	tau_tau_b_b->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	tau_tau_b_b->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	tau_tau_b_b->Branch("gen_h_tt_E", &gen_h_tt_E);
	tau_tau_b_b->Branch("gen_mctMatch", &gen_mctMatch);
	tau_tau_b_b->Branch("gen_weight", &weight);
	TTree* e_e_b_b = new TTree("e_e_b_b", "e e b #bar{b}");
	e_e_b_b->Branch("t_0_pT", &t_0_pT);
	e_e_b_b->Branch("t_0_eta", &t_0_eta);
	e_e_b_b->Branch("t_0_phi", &t_0_phi);
	e_e_b_b->Branch("t_0_mass", &t_0_mass);
	e_e_b_b->Branch("t_1_pT", &t_1_pT);
	e_e_b_b->Branch("t_1_eta", &t_1_eta);
	e_e_b_b->Branch("t_1_phi", &t_1_phi);
	e_e_b_b->Branch("t_1_mass", &t_1_mass);
	e_e_b_b->Branch("b_0_pT", &b_0_pT);
	e_e_b_b->Branch("b_0_eta", &b_0_eta);
	e_e_b_b->Branch("b_0_phi", &b_0_phi);
	e_e_b_b->Branch("b_0_mass", &b_0_mass);
	e_e_b_b->Branch("b_1_pT", &b_1_pT);
	e_e_b_b->Branch("b_1_eta", &b_1_eta);
	e_e_b_b->Branch("b_1_phi", &b_1_phi);
	e_e_b_b->Branch("b_1_mass", &b_1_mass);
	e_e_b_b->Branch("mPT_pT", &mPT_pT);
	e_e_b_b->Branch("mPT_phi", &mPT_phi);
	e_e_b_b->Branch("h_tt_pT", &h_tt_pT);
	e_e_b_b->Branch("h_tt_eta", &h_tt_eta);
	e_e_b_b->Branch("h_tt_phi", &h_tt_phi);
	e_e_b_b->Branch("h_tt_mass", &h_tt_mass);
	e_e_b_b->Branch("h_bb_pT", &h_bb_pT);
	e_e_b_b->Branch("h_bb_eta", &h_bb_eta);
	e_e_b_b->Branch("h_bb_phi", &h_bb_phi);
	e_e_b_b->Branch("h_bb_mass", &h_bb_mass);
	e_e_b_b->Branch("diH_pT", &diH_pT);
	e_e_b_b->Branch("diH_eta", &diH_eta);
	e_e_b_b->Branch("diH_phi", &diH_phi);
	e_e_b_b->Branch("diH_mass", &diH_mass);
	e_e_b_b->Branch("hT", &hT);
	e_e_b_b->Branch("sT", &sT);
	e_e_b_b->Branch("centrality", &centrality);
	e_e_b_b->Branch("eVis", &eVis);
	e_e_b_b->Branch("nJets", &nJets);
	e_e_b_b->Branch("nBJets", &nBJets);
	e_e_b_b->Branch("nTauJets", &nTauJets);
	e_e_b_b->Branch("minJetPT", &minJetPT);
	e_e_b_b->Branch("meanJetPT", &meanJetPT);
	e_e_b_b->Branch("maxJetPT", &maxJetPT);
	e_e_b_b->Branch("minJetMass", &minJetMass);
	e_e_b_b->Branch("meanJetMass", &meanJetMass);
	e_e_b_b->Branch("maxJetMass", &maxJetMass);
	e_e_b_b->Branch("minJetEta", &minJetEta);
	e_e_b_b->Branch("meanJetEta", &meanJetEta);
	e_e_b_b->Branch("maxJetEta", &maxJetEta);
	e_e_b_b->Branch("sphericityA", &sphericityA);
	e_e_b_b->Branch("spherocityA", &spherocityA);
	e_e_b_b->Branch("aplanarityA", &aplanarityA);
	e_e_b_b->Branch("aplanorityA", &aplanorityA);
	e_e_b_b->Branch("upsilonA", &upsilonA);
	e_e_b_b->Branch("dShapeA", &dShapeA);
	e_e_b_b->Branch("sphericityP", &sphericityP);
	e_e_b_b->Branch("spherocityP", &spherocityP);
	e_e_b_b->Branch("aplanarityP", &aplanarityP);
	e_e_b_b->Branch("aplanorityP", &aplanorityP);
	e_e_b_b->Branch("upsilonP", &upsilonP);
	e_e_b_b->Branch("dShapeP", &dShapeP);
	e_e_b_b->Branch("gen_t_0_pT", &gen_t_0_pT);
	e_e_b_b->Branch("gen_t_0_eta", &gen_t_0_eta);
	e_e_b_b->Branch("gen_t_0_phi", &gen_t_0_phi);
	e_e_b_b->Branch("gen_t_0_E", &gen_t_0_E);
	e_e_b_b->Branch("gen_t_1_pT", &gen_t_1_pT);
	e_e_b_b->Branch("gen_t_1_eta", &gen_t_1_eta);
	e_e_b_b->Branch("gen_t_1_phi", &gen_t_1_phi);
	e_e_b_b->Branch("gen_t_1_E", &gen_t_1_E);
	e_e_b_b->Branch("gen_b_0_pT", &gen_b_0_pT);
	e_e_b_b->Branch("gen_b_0_eta", &gen_b_0_eta);
	e_e_b_b->Branch("gen_b_0_phi", &gen_b_0_phi);
	e_e_b_b->Branch("gen_b_0_E", &gen_b_0_E);
	e_e_b_b->Branch("gen_b_1_pT", &gen_b_1_pT);
	e_e_b_b->Branch("gen_b_1_eta", &gen_b_1_eta);
	e_e_b_b->Branch("gen_b_1_phi", &gen_b_1_phi);
	e_e_b_b->Branch("gen_b_1_E", &gen_b_1_E);
	e_e_b_b->Branch("gen_diH_pT", &gen_diH_pT);
	e_e_b_b->Branch("gen_diH_eta", &gen_diH_eta);
	e_e_b_b->Branch("gen_diH_phi", &gen_diH_phi);
	e_e_b_b->Branch("gen_diH_E", &gen_diH_E);
	e_e_b_b->Branch("gen_diH_mass", &gen_diH_mass);
	e_e_b_b->Branch("gen_h_bb_pT", &gen_h_bb_pT);
	e_e_b_b->Branch("gen_h_bb_eta", &gen_h_bb_eta);
	e_e_b_b->Branch("gen_h_bb_phi", &gen_h_bb_phi);
	e_e_b_b->Branch("gen_h_bb_E", &gen_h_bb_E);
	e_e_b_b->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	e_e_b_b->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	e_e_b_b->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	e_e_b_b->Branch("gen_h_tt_E", &gen_h_tt_E);
	e_e_b_b->Branch("gen_mctMatch", &gen_mctMatch);
	e_e_b_b->Branch("gen_weight", &weight);
	TTree* mu_mu_b_b = new TTree("mu_mu_b_b", "#mu #mu b #bar{b}");
	mu_mu_b_b->Branch("t_0_pT", &t_0_pT);
	mu_mu_b_b->Branch("t_0_eta", &t_0_eta);
	mu_mu_b_b->Branch("t_0_phi", &t_0_phi);
	mu_mu_b_b->Branch("t_0_mass", &t_0_mass);
	mu_mu_b_b->Branch("t_1_pT", &t_1_pT);
	mu_mu_b_b->Branch("t_1_eta", &t_1_eta);
	mu_mu_b_b->Branch("t_1_phi", &t_1_phi);
	mu_mu_b_b->Branch("t_1_mass", &t_1_mass);
	mu_mu_b_b->Branch("b_0_pT", &b_0_pT);
	mu_mu_b_b->Branch("b_0_eta", &b_0_eta);
	mu_mu_b_b->Branch("b_0_phi", &b_0_phi);
	mu_mu_b_b->Branch("b_0_mass", &b_0_mass);
	mu_mu_b_b->Branch("b_1_pT", &b_1_pT);
	mu_mu_b_b->Branch("b_1_eta", &b_1_eta);
	mu_mu_b_b->Branch("b_1_phi", &b_1_phi);
	mu_mu_b_b->Branch("b_1_mass", &b_1_mass);
	mu_mu_b_b->Branch("mPT_pT", &mPT_pT);
	mu_mu_b_b->Branch("mPT_phi", &mPT_phi);
	mu_mu_b_b->Branch("h_tt_pT", &h_tt_pT);
	mu_mu_b_b->Branch("h_tt_eta", &h_tt_eta);
	mu_mu_b_b->Branch("h_tt_phi", &h_tt_phi);
	mu_mu_b_b->Branch("h_tt_mass", &h_tt_mass);
	mu_mu_b_b->Branch("h_bb_pT", &h_bb_pT);
	mu_mu_b_b->Branch("h_bb_eta", &h_bb_eta);
	mu_mu_b_b->Branch("h_bb_phi", &h_bb_phi);
	mu_mu_b_b->Branch("h_bb_mass", &h_bb_mass);
	mu_mu_b_b->Branch("diH_pT", &diH_pT);
	mu_mu_b_b->Branch("diH_eta", &diH_eta);
	mu_mu_b_b->Branch("diH_phi", &diH_phi);
	mu_mu_b_b->Branch("diH_mass", &diH_mass);
	mu_mu_b_b->Branch("hT", &hT);
	mu_mu_b_b->Branch("sT", &sT);
	mu_mu_b_b->Branch("centrality", &centrality);
	mu_mu_b_b->Branch("eVis", &eVis);
	mu_mu_b_b->Branch("nJets", &nJets);
	mu_mu_b_b->Branch("nBJets", &nBJets);
	mu_mu_b_b->Branch("nTauJets", &nTauJets);
	mu_mu_b_b->Branch("minJetPT", &minJetPT);
	mu_mu_b_b->Branch("meanJetPT", &meanJetPT);
	mu_mu_b_b->Branch("maxJetPT", &maxJetPT);
	mu_mu_b_b->Branch("minJetMass", &minJetMass);
	mu_mu_b_b->Branch("meanJetMass", &meanJetMass);
	mu_mu_b_b->Branch("maxJetMass", &maxJetMass);
	mu_mu_b_b->Branch("minJetEta", &minJetEta);
	mu_mu_b_b->Branch("meanJetEta", &meanJetEta);
	mu_mu_b_b->Branch("maxJetEta", &maxJetEta);
	mu_mu_b_b->Branch("sphericityA", &sphericityA);
	mu_mu_b_b->Branch("spherocityA", &spherocityA);
	mu_mu_b_b->Branch("aplanarityA", &aplanarityA);
	mu_mu_b_b->Branch("aplanorityA", &aplanorityA);
	mu_mu_b_b->Branch("upsilonA", &upsilonA);
	mu_mu_b_b->Branch("dShapeA", &dShapeA);
	mu_mu_b_b->Branch("sphericityP", &sphericityP);
	mu_mu_b_b->Branch("spherocityP", &spherocityP);
	mu_mu_b_b->Branch("aplanarityP", &aplanarityP);
	mu_mu_b_b->Branch("aplanorityP", &aplanorityP);
	mu_mu_b_b->Branch("upsilonP", &upsilonP);
	mu_mu_b_b->Branch("dShapeP", &dShapeP);
	mu_mu_b_b->Branch("gen_t_0_pT", &gen_t_0_pT);
	mu_mu_b_b->Branch("gen_t_0_eta", &gen_t_0_eta);
	mu_mu_b_b->Branch("gen_t_0_phi", &gen_t_0_phi);
	mu_mu_b_b->Branch("gen_t_0_E", &gen_t_0_E);
	mu_mu_b_b->Branch("gen_t_1_pT", &gen_t_1_pT);
	mu_mu_b_b->Branch("gen_t_1_eta", &gen_t_1_eta);
	mu_mu_b_b->Branch("gen_t_1_phi", &gen_t_1_phi);
	mu_mu_b_b->Branch("gen_t_1_E", &gen_t_1_E);
	mu_mu_b_b->Branch("gen_b_0_pT", &gen_b_0_pT);
	mu_mu_b_b->Branch("gen_b_0_eta", &gen_b_0_eta);
	mu_mu_b_b->Branch("gen_b_0_phi", &gen_b_0_phi);
	mu_mu_b_b->Branch("gen_b_0_E", &gen_b_0_E);
	mu_mu_b_b->Branch("gen_b_1_pT", &gen_b_1_pT);
	mu_mu_b_b->Branch("gen_b_1_eta", &gen_b_1_eta);
	mu_mu_b_b->Branch("gen_b_1_phi", &gen_b_1_phi);
	mu_mu_b_b->Branch("gen_b_1_E", &gen_b_1_E);
	mu_mu_b_b->Branch("gen_diH_pT", &gen_diH_pT);
	mu_mu_b_b->Branch("gen_diH_eta", &gen_diH_eta);
	mu_mu_b_b->Branch("gen_diH_phi", &gen_diH_phi);
	mu_mu_b_b->Branch("gen_diH_E", &gen_diH_E);
	mu_mu_b_b->Branch("gen_diH_mass", &gen_diH_mass);
	mu_mu_b_b->Branch("gen_h_bb_pT", &gen_h_bb_pT);
	mu_mu_b_b->Branch("gen_h_bb_eta", &gen_h_bb_eta);
	mu_mu_b_b->Branch("gen_h_bb_phi", &gen_h_bb_phi);
	mu_mu_b_b->Branch("gen_h_bb_E", &gen_h_bb_E);
	mu_mu_b_b->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	mu_mu_b_b->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	mu_mu_b_b->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	mu_mu_b_b->Branch("gen_h_tt_E", &gen_h_tt_E);
	mu_mu_b_b->Branch("gen_mctMatch", &gen_mctMatch);
	mu_mu_b_b->Branch("gen_weight", &weight);
	TTree* e_mu_b_b = new TTree("e_mu_b_b", "e #mu b #bar{b}");
	e_mu_b_b->Branch("t_0_pT", &t_0_pT);
	e_mu_b_b->Branch("t_0_eta", &t_0_eta);
	e_mu_b_b->Branch("t_0_phi", &t_0_phi);
	e_mu_b_b->Branch("t_0_mass", &t_0_mass);
	e_mu_b_b->Branch("t_1_pT", &t_1_pT);
	e_mu_b_b->Branch("t_1_eta", &t_1_eta);
	e_mu_b_b->Branch("t_1_phi", &t_1_phi);
	e_mu_b_b->Branch("t_1_mass", &t_1_mass);
	e_mu_b_b->Branch("b_0_pT", &b_0_pT);
	e_mu_b_b->Branch("b_0_eta", &b_0_eta);
	e_mu_b_b->Branch("b_0_phi", &b_0_phi);
	e_mu_b_b->Branch("b_0_mass", &b_0_mass);
	e_mu_b_b->Branch("b_1_pT", &b_1_pT);
	e_mu_b_b->Branch("b_1_eta", &b_1_eta);
	e_mu_b_b->Branch("b_1_phi", &b_1_phi);
	e_mu_b_b->Branch("b_1_mass", &b_1_mass);
	e_mu_b_b->Branch("mPT_pT", &mPT_pT);
	e_mu_b_b->Branch("mPT_phi", &mPT_phi);
	e_mu_b_b->Branch("h_tt_pT", &h_tt_pT);
	e_mu_b_b->Branch("h_tt_eta", &h_tt_eta);
	e_mu_b_b->Branch("h_tt_phi", &h_tt_phi);
	e_mu_b_b->Branch("h_tt_mass", &h_tt_mass);
	e_mu_b_b->Branch("h_bb_pT", &h_bb_pT);
	e_mu_b_b->Branch("h_bb_eta", &h_bb_eta);
	e_mu_b_b->Branch("h_bb_phi", &h_bb_phi);
	e_mu_b_b->Branch("h_bb_mass", &h_bb_mass);
	e_mu_b_b->Branch("diH_pT", &diH_pT);
	e_mu_b_b->Branch("diH_eta", &diH_eta);
	e_mu_b_b->Branch("diH_phi", &diH_phi);
	e_mu_b_b->Branch("diH_mass", &diH_mass);
	e_mu_b_b->Branch("hT", &hT);
	e_mu_b_b->Branch("sT", &sT);
	e_mu_b_b->Branch("centrality", &centrality);
	e_mu_b_b->Branch("eVis", &eVis);
	e_mu_b_b->Branch("nJets", &nJets);
	e_mu_b_b->Branch("nBJets", &nBJets);
	e_mu_b_b->Branch("nTauJets", &nTauJets);
	e_mu_b_b->Branch("minJetPT", &minJetPT);
	e_mu_b_b->Branch("meanJetPT", &meanJetPT);
	e_mu_b_b->Branch("maxJetPT", &maxJetPT);
	e_mu_b_b->Branch("minJetMass", &minJetMass);
	e_mu_b_b->Branch("meanJetMass", &meanJetMass);
	e_mu_b_b->Branch("maxJetMass", &maxJetMass);
	e_mu_b_b->Branch("minJetEta", &minJetEta);
	e_mu_b_b->Branch("meanJetEta", &meanJetEta);
	e_mu_b_b->Branch("maxJetEta", &maxJetEta);
	e_mu_b_b->Branch("sphericityA", &sphericityA);
	e_mu_b_b->Branch("spherocityA", &spherocityA);
	e_mu_b_b->Branch("aplanarityA", &aplanarityA);
	e_mu_b_b->Branch("aplanorityA", &aplanorityA);
	e_mu_b_b->Branch("upsilonA", &upsilonA);
	e_mu_b_b->Branch("dShapeA", &dShapeA);
	e_mu_b_b->Branch("sphericityP", &sphericityP);
	e_mu_b_b->Branch("spherocityP", &spherocityP);
	e_mu_b_b->Branch("aplanarityP", &aplanarityP);
	e_mu_b_b->Branch("aplanorityP", &aplanorityP);
	e_mu_b_b->Branch("upsilonP", &upsilonP);
	e_mu_b_b->Branch("dShapeP", &dShapeP);
	e_mu_b_b->Branch("gen_t_0_pT", &gen_t_0_pT);
	e_mu_b_b->Branch("gen_t_0_eta", &gen_t_0_eta);
	e_mu_b_b->Branch("gen_t_0_phi", &gen_t_0_phi);
	e_mu_b_b->Branch("gen_t_0_E", &gen_t_0_E);
	e_mu_b_b->Branch("gen_t_1_pT", &gen_t_1_pT);
	e_mu_b_b->Branch("gen_t_1_eta", &gen_t_1_eta);
	e_mu_b_b->Branch("gen_t_1_phi", &gen_t_1_phi);
	e_mu_b_b->Branch("gen_t_1_E", &gen_t_1_E);
	e_mu_b_b->Branch("gen_b_0_pT", &gen_b_0_pT);
	e_mu_b_b->Branch("gen_b_0_eta", &gen_b_0_eta);
	e_mu_b_b->Branch("gen_b_0_phi", &gen_b_0_phi);
	e_mu_b_b->Branch("gen_b_0_E", &gen_b_0_E);
	e_mu_b_b->Branch("gen_b_1_pT", &gen_b_1_pT);
	e_mu_b_b->Branch("gen_b_1_eta", &gen_b_1_eta);
	e_mu_b_b->Branch("gen_b_1_phi", &gen_b_1_phi);
	e_mu_b_b->Branch("gen_b_1_E", &gen_b_1_E);
	e_mu_b_b->Branch("gen_diH_pT", &gen_diH_pT);
	e_mu_b_b->Branch("gen_diH_eta", &gen_diH_eta);
	e_mu_b_b->Branch("gen_diH_phi", &gen_diH_phi);
	e_mu_b_b->Branch("gen_diH_E", &gen_diH_E);
	e_mu_b_b->Branch("gen_diH_mass", &gen_diH_mass);
	e_mu_b_b->Branch("gen_h_bb_pT", &gen_h_bb_pT);
	e_mu_b_b->Branch("gen_h_bb_eta", &gen_h_bb_eta);
	e_mu_b_b->Branch("gen_h_bb_phi", &gen_h_bb_phi);
	e_mu_b_b->Branch("gen_h_bb_E", &gen_h_bb_E);
	e_mu_b_b->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	e_mu_b_b->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	e_mu_b_b->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	e_mu_b_b->Branch("gen_h_tt_E", &gen_h_tt_E);
	e_mu_b_b->Branch("gen_mctMatch", &gen_mctMatch);
	e_mu_b_b->Branch("gen_weight", &weight);
	std::cout << "Variables initialised\n";
	//___________________________________________
	//Initialise plots___________________________
	std::cout << "Initialising plot\n";
	TH1D* h_datasetSizes = new TH1D("Dataset_sizes", "Dataset sizes", 7, -0.7, 0.7);
	TH1D* h_e_tau_b_b_cutFlow;
	TH1D* h_mu_tau_b_b_cutFlow;
	TH1D* h_tau_tau_b_b_cutFlow;
	TH1D* h_e_e_b_b_cutFlow;
	TH1D* h_mu_mu_b_b_cutFlow;
	TH1D* h_e_mu_b_b_cutFlow;
	h_e_tau_b_b_cutFlow = new TH1D("e_tau_b_b_Cut_Flow", "e #tau_{h} b #bar{b} cut flow", 9, -0.9, 0.9);
	h_mu_tau_b_b_cutFlow = new TH1D("mu_tau_b_b_Cut_Flow", "#mu #tau_{h} b #bar{b} cut flow", 9, -0.9, 0.9);
	h_tau_tau_b_b_cutFlow = new TH1D("tau_tau_b_b_Cut_Flow", "#tau_{h} #tau_{h} b #bar{b} cut flow", 8, -0.8, 0.8);
	h_e_e_b_b_cutFlow = new TH1D("e_e_b_b_Cut_Flow", "e e b #bar{b} cut flow", 9, -0.9, 0.9);
	h_mu_mu_b_b_cutFlow = new TH1D("mu_mu_b_b_Cut_Flow", "#mu #mu b #bar{b} cut flow", 9, -0.9, 0.9);
	h_e_mu_b_b_cutFlow = new TH1D("e_mu_b_b_Cut_Flow", "e #mu b #bar{b} cut flow", 9, -0.9, 0.9);
	h_datasetSizes->GetXaxis()->SetBinLabel(1, "All");
	h_datasetSizes->GetXaxis()->SetBinLabel(2, "#mu #tau_{h} b #bar{b}");
	h_datasetSizes->GetXaxis()->SetBinLabel(3, "e #tau_{h} b #bar{b}");
	h_datasetSizes->GetXaxis()->SetBinLabel(4, "#tau_{h} #tau_{h} b #bar{b}");
	h_datasetSizes->GetXaxis()->SetBinLabel(5, "e e b #bar{b}");
	h_datasetSizes->GetXaxis()->SetBinLabel(6, "e #mu b #bar{b}");
	h_datasetSizes->GetXaxis()->SetBinLabel(7, "#mu #mu b #bar{b}");
	h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(1, "All");
	h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(2, "Quality #tau");
	h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(3, "Quality b#bar{b}");
	h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(4, "Quality e");
	h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(5, "1 e & 0 #mu");
	h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(6, "OS");
	h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(7, "m_{#tau#tau} Cut");
	h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(8, "m_{b#bar{b}} Cut");
	h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(1, "All");
	h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(2, "Quality #tau");
	h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(3, "Quality b#bar{b}");
	h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(4, "Quality #mu");
	h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(5, "1 #mu & 0 e");
	h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(6, "OS");
	h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(7, "m_{#tau#tau} Cut");
	h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(8, "m_{b#bar{b}} Cut");
	h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(1, "All");
	h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(2, "Quality #tau#tau");
	h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(3, "Quality b#bar{b}");
	h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(4, "0 e & 0 #mu");
	h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(5, "OS");
	h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(6, "m_{#tau#tau} Cut");
	h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(7, "m_{b#bar{b}} Cut");
	h_e_e_b_b_cutFlow->GetXaxis()->SetBinLabel(1, "All");
	h_e_e_b_b_cutFlow->GetXaxis()->SetBinLabel(2, "Quality di-e");
	h_e_e_b_b_cutFlow->GetXaxis()->SetBinLabel(3, "2 e & 0 #mu");
	h_e_e_b_b_cutFlow->GetXaxis()->SetBinLabel(4, "OS");
	h_e_e_b_b_cutFlow->GetXaxis()->SetBinLabel(5, "Quality b#bar{b}");
	h_e_e_b_b_cutFlow->GetXaxis()->SetBinLabel(6, "0 #tau_{h}");
	h_e_e_b_b_cutFlow->GetXaxis()->SetBinLabel(7, "m_{#tau#tau} Cut");
	h_e_e_b_b_cutFlow->GetXaxis()->SetBinLabel(8, "m_{b#bar{b}} Cut");
	h_e_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(1, "All");
	h_e_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(2, "Quality e and #mu");
	h_e_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(3, "1 e & 1 #mu");
	h_e_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(4, "OS");
	h_e_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(5, "Quality b#bar{b}");
	h_e_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(6, "0 #tau_{h}");
	h_e_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(7, "m_{#tau#tau} Cut");
	h_e_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(8, "m_{b#bar{b}} Cut");
	h_mu_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(1, "All");
	h_mu_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(2, "Quality di-#mu");
	h_mu_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(3, "2 #mu & 0 e");
	h_mu_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(4, "OS");
	h_mu_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(5, "Quality b#bar{b}");
	h_mu_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(6, "0 #tau_{h}");
	h_mu_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(7, "m_{#tau#tau} Cut");
	h_mu_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(8, "m_{b#bar{b}} Cut");
	std::cout << "Plots initialised\n";
	//___________________________________________
	//Load data__________________________________
	std::cout << "Running event selection\n";
	TChain *chain = new TChain("Delphes");
	chain->Add(options["-i"].c_str());
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	TClonesArray *branchMuon = treeReader->UseBranch("MuonLoose");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchJet = treeReader->UseBranch("Jet");
	TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
	std::cout << "Data loaded\n";
	//_______________________________________
	//Loop through events____________________
	Long64_t nEvents = treeReader->GetEntries();
	std::cout << "Total number of events in file: " << nEvents << "\n";
	std::vector<int> taus, bJets, electrons, muons;
	bool addMuon, addElectron;
	TLorentzVector v_tau_0, v_tau_1, v_bJet_0, v_bJet_1, v_higgs_tt, v_higgs_bb, v_diHiggs;
	Jet* tmpJet;
	Electron* tmpElectron;
	Muon* tmpMuon;
	std::cout << "Beginning event loop\n";
	for (Long64_t cEvent = 0; cEvent < nEvents; cEvent++) {
		if (debug) std::cout << "Loading event " << cEvent << "\n";
		treeReader->ReadEntry(entry); //Load next event
		if (debug) std::cout << "Event loaded, getting data\n";
		if (debug) std::cout << "Data loaded\n";
		if (cEvent % 1000 == 0) std::cout << "Loop: " << cEvent << "/" << nEvents << ", " <<
				100*cEvent/nEvents << "%\n";
		h_datasetSizes->Fill("All", 1);
		eventAccepted = false;
		//Check for mu tau b b finalstates___
		h_mu_tau_b_b_cutFlow->Fill("All", 1);
		finalstateSet("mu_tau_b_b");
		electrons.clear();
		muons.clear();
		taus.clear();
		bJets.clear();
		addMuon = false;
		addElectron = false;
		if (debug) std::cout << "Running mu tau b b\n";
		for (int i = 0; i < branchMuon->GetEntries(); i++) { //Loop through muons
			tmpMuon = (Muon*) branchMuon->At(i);
			if (tmpMuon->PT > muPTMin && std::abs(tmpMuon->Eta) < muEtaMax
					&& tmpMuon->IsolationVar < muIsoMax) { //Quality muon
				muons.push_back(i);
			} else if (tmpMuon->PT > muPTMinAdd && std::abs(tmpMuon->Eta) < muEtaMaxAdd
					&& tmpMuon->IsolationVar < muIsoMaxAdd) { //Additional muon
				addMuon = true;
				break;
			}
		}
		if (muons.size() == 1 && !addMuon) { //One quality muon found and no additional muons
			h_mu_tau_b_b_cutFlow->Fill("Quality #mu", 1);
			for (int i = 0; i < branchElectron->GetEntries(); i++) { //Loop through electrons
				tmpElectron = (Electron*) branchElectron->At(i);
				if (tmpElectron->PT > ePTMinAdd && std::abs(tmpElectron->Eta) < eEtaMaxAdd
						&& tmpElectron->IsolationVar < eIsoMaxAdd) { //Additional electron
					addElectron = true
					break;
				}
			}
			if (!addElectron) { //No additional electrons found
				h_mu_tau_b_b_cutFlow->Fill("1 #mu & 0 e", 1);
				for (int i = 0; i < branchJet->GetEntries(); i++) { //Loop through jets
					tmpJet = (Jet*) branchJet->At(i);
					if (tmpJet->TauTag == 1 && tmpJet->BTag == 0 && tmpJet->PT > tauPTMin
							&& std::abs(tmpJet->Eta) < tauEtaMax
							&& tmpJet->Charge != branchMuon->At(muons[0])->Charge) { //Quality  OS tau
						taus.push_back(i);
					}
					if (tmpJet->TauTag == 0 && tmpJet->BTag == 1 && tmpJet->PT > bJetPTMin
							&& std::abs(tmpJet->Eta) < bJetEtaMax) { //Quality b jet
						bJets.push_back(i);
					}
				}
				if (taus.size() >= 1) {//Quality tau
					h_mu_tau_b_b_cutFlow->Fill("Quality #tau", 1);
					if (bJets.size() >= 2) {//Quality b jets pairs found
						h_mu_tau_b_b_cutFlow->Fill("Quality b#bar{b}", 1);
						if (selectBJets(branchJet, &bJets, &bJet_0, &bJet_1) == true) { //Quality b-jet pair found
							v_tau_1 = branchMuon->At(muons[0])->P4;
							v_tau_0 = branchJet->At(taus[0])->P4;
							v_higgs_tt = getHiggs2Taus(branchMissingET, v_tau_0, v_tau_1);
							v_bJet_0 = getBJet(reader, bJet_0);
							v_bJet_1 = getBJet(reader, bJet_1);
							v_higgs_bb = getHiggs2Bs(v_bJet_0, v_bJet_1);
							v_diHiggs = getDiHiggs(v_higgs_tt, v_higgs_bb);
							if (debug) std::cout << "Accepted mu_tau_b_b event\n";
							t_0_pT = v_tau_0.Pt();
							t_0_eta = v_tau_0.Eta();
							t_0_phi = v_tau_0.Phi();
							t_0_mass = v_tau_0.M();
							t_1_pT = v_tau_1.Pt();
							t_1_eta = v_tau_1.Eta();
							t_1_phi = v_tau_1.Phi();
							t_1_mass = muMass;
							b_0_pT = v_bJet_0.Pt();
							b_0_eta = v_bJet_0.Eta();
							b_0_phi = v_bJet_0.Phi();
							b_0_mass = v_bJet_0.M();
							b_1_pT = v_bJet_1.Pt();
							b_1_eta = v_bJet_1.Eta();
							b_1_phi = v_bJet_1.Phi();
							b_1_mass = v_bJet_1.M();
							mPT_pT = branchMissingET->At(0)->MET;
							mPT_phi = branchMissingET->At(0)->Phi;
							h_tt_pT = v_higgs_tt.Pt();
							h_tt_eta = v_higgs_tt.Eta();
							h_tt_phi = v_higgs_tt.Phi();
							h_tt_mass = v_higgs_tt.M();
							h_bb_pT = v_higgs_bb.Pt();
							h_bb_eta = v_higgs_bb.Eta();
							h_bb_phi = v_higgs_bb.Phi();
							h_bb_mass = v_higgs_bb.M();
							diH_pT = v_diHiggs.Pt();
							diH_eta = v_diHiggs.Eta();
							diH_phi = v_diHiggs.Phi();
							diH_mass = v_diHiggs.M();
							getGlobalEventInfo(options["-i"], cEvent,
									&hT, &sT, &centrality, &eVis,
									&nJets, &nBJets, &nTauJets,
									&minJetPT, &meanJetPT, &maxJetPT,
									&minJetMass, &meanJetMass, &maxJetMass,
									&minJetEta, &meanJetEta, &maxJetEta,
									&sphericityA, &spherocityA,
									&aplanarityA, &aplanorityA,
									&upsilonA, &dShapeA);
							getPrimaryEventShapes(v_tau_0, v_tau_1, v_bJet_0, v_bJet_1,
									&sphericityP, &spherocityP,
									&aplanarityP, &aplanorityP,
									&upsilonP, &dShapeP);
							weight = (double)*reader->Event_Weight;
							mu_tau_b_b->Fill();
							h_datasetSizes->Fill("#mu #tau_{h} b #bar{b}", 1);
							eventAccepted = true;
						}
					}
				}
			}
		//___________________________________
		if (eventAccepted) continue;
		}
		std::cout << "Event loop complete\n";
	}
	std::cout << "All files complete\n";
	//___________________________________________
	//Writing plots______________________________
	TFile* outputFile = new TFile(("../outputs/" + outputName + "/" + outputName + ".root").c_str(), "recreate");
	outputFile->cd();
	std::cout << "Creating plots\n";
	TCanvas* c_datasetSizes = new TCanvas();
	c_datasetSizes->SetLogy();
	h_datasetSizes->GetXaxis()->SetTitle("Dataset");
	h_datasetSizes->GetYaxis()->SetTitle("Events");
	h_datasetSizes->Draw();
	h_datasetSizes->Write();
	c_datasetSizes->Print(("../outputs/" + outputName + "/datasetSizes.pdf").c_str());
	delete c_datasetSizes;
	TCanvas* c_e_tau_b_b_cutFlow = new TCanvas();
	c_e_tau_b_b_cutFlow->SetLogy();
	h_e_tau_b_b_cutFlow->GetXaxis()->SetTitle("Cuts");
	h_e_tau_b_b_cutFlow->GetYaxis()->SetTitle("Events");
	h_e_tau_b_b_cutFlow->Draw();
	h_e_tau_b_b_cutFlow->Write();
	c_e_tau_b_b_cutFlow->Print(("../outputs/" + outputName + "/e_tau_b_b_cutFlow.pdf").c_str());
	delete c_e_tau_b_b_cutFlow;
	TCanvas* c_mu_tau_b_b_cutFlow = new TCanvas();
	c_mu_tau_b_b_cutFlow->SetLogy();
	h_mu_tau_b_b_cutFlow->GetXaxis()->SetTitle("Cuts");
	h_mu_tau_b_b_cutFlow->GetYaxis()->SetTitle("Events");
	h_mu_tau_b_b_cutFlow->Draw();
	h_mu_tau_b_b_cutFlow->Write();
	c_mu_tau_b_b_cutFlow->Print(("../outputs/" + outputName + "/mu_tau_b_b_cutFlow.pdf").c_str());
	delete c_mu_tau_b_b_cutFlow;
	TCanvas* c_tau_tau_b_b_cutFlow = new TCanvas();
	c_tau_tau_b_b_cutFlow->SetLogy();
	h_tau_tau_b_b_cutFlow->GetXaxis()->SetTitle("Cuts");
	h_tau_tau_b_b_cutFlow->GetYaxis()->SetTitle("Events");
	h_tau_tau_b_b_cutFlow->Draw();
	h_tau_tau_b_b_cutFlow->Write();
	c_tau_tau_b_b_cutFlow->Print(("../outputs/" + outputName + "/tau_tau_b_b_cutFlow.pdf").c_str());
	delete c_tau_tau_b_b_cutFlow;
	TCanvas* c_mu_mu_b_b_cutFlow = new TCanvas();
	c_mu_mu_b_b_cutFlow->SetLogy();
	h_mu_mu_b_b_cutFlow->GetXaxis()->SetTitle("Cuts");
	h_mu_mu_b_b_cutFlow->GetYaxis()->SetTitle("Events");
	h_mu_mu_b_b_cutFlow->Draw();
	h_mu_mu_b_b_cutFlow->Write();
	c_mu_mu_b_b_cutFlow->Print(("../outputs/" + outputName + "/mu_mu_b_b_cutFlow.pdf").c_str());
	delete c_mu_mu_b_b_cutFlow;
	TCanvas* c_e_mu_b_b_cutFlow = new TCanvas();
	c_e_mu_b_b_cutFlow->SetLogy();
	h_e_mu_b_b_cutFlow->GetXaxis()->SetTitle("Cuts");
	h_e_mu_b_b_cutFlow->GetYaxis()->SetTitle("Events");
	h_e_mu_b_b_cutFlow->Draw();
	h_e_mu_b_b_cutFlow->Write();
	c_e_mu_b_b_cutFlow->Print(("../outputs/" + outputName + "/e_mu_b_b_cutFlow.pdf").c_str());
	delete c_e_mu_b_b_cutFlow;
	TCanvas* c_e_e_b_b_cutFlow = new TCanvas();
	c_e_e_b_b_cutFlow->SetLogy();
	h_e_e_b_b_cutFlow->GetXaxis()->SetTitle("Cuts");
	h_e_e_b_b_cutFlow->GetYaxis()->SetTitle("Events");
	h_e_e_b_b_cutFlow->Draw();
	h_e_e_b_b_cutFlow->Write();
	c_e_e_b_b_cutFlow->Print(("../outputs/" + outputName + "/e_e_b_b_cutFlow.pdf").c_str());
	delete c_e_e_b_b_cutFlow;
	std::cout << "Plots created\n";
	//___________________________________________
	//Save datasets______________________________
	std::cout << "Saving data\n";
	e_tau_b_b->Write();
	mu_tau_b_b->Write();
	tau_tau_b_b->Write();
	mu_mu_b_b->Write();
	e_mu_b_b->Write();
	e_e_b_b->Write();
	std::cout << "Data saved\n";
	outputFile->Close();
	delete outputFile;
	delete e_tau_b_b;
	delete mu_tau_b_b;
	delete tau_tau_b_b;
	delete mu_mu_b_b;
	delete e_mu_b_b;
	delete e_e_b_b;
	delete h_datasetSizes;
	delete h_e_tau_b_b_cutFlow;
	delete h_mu_tau_b_b_cutFlow;
	delete h_tau_tau_b_b_cutFlow;
	delete h_e_e_b_b_cutFlow;
	delete h_mu_mu_b_b_cutFlow;
	delete h_e_mu_b_b_cutFlow;
	//___________________________________________
}
