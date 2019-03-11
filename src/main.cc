/*
 * main.cc
 *
 *  Created on: 8 Apr 2016
 *      Author: giles giles.strong@outlook.com
 */
//Local
#include "main.hh"

bool debug = false;

std::map<int, double> tauEfficiency = {{21, 0.31}, {32, 0.45}, {54, 0.61}};

std::map<int, double> tauMistagRate30 = {{21, 0.0025}, {32, 0.0049}, {54, 0.0118}};
std::map<int, double> tauMistagRate23 = {{21, 0.0023}, {32, 0.0045}, {54, 0.0110}};
std::map<int, double> tauMistagRate14 = {{21, 0.0022}, {32, 0.0041}, {54, 0.0099}};

std::map<int, double> tauFakeFactor30 = {{21, 1.00}, {32, 2.05}, {54, 5.10}};
std::map<int, double> tauFakeFactor23 = {{21, 0.88}, {32, 1.80}, {54, 4.56}};
std::map<int, double> tauFakeFactor14 = {{21, 0.52}, {32, 1.00}, {54, 2.41}};

std::map<int, double>  mvaFakeFactorReduction = {{21, 3.0}, {32, 2.0}, {54, 2.0}};

#define n_Ai_coeffs 15

std::array<double, n_Ai_coeffs> A_integralXS = {
        2.100318379,
        10.2,
        0.287259045,
        0.098882779,
        1.321736614,
        -8.42431259,
        -1.388017366,
        2.8,
        0.518124457,
        -2.163473227,
        -0.550668596,
        5.871490593,
        0.296671491,
        -1.172793054,
        0.653429812
    };

int klambda_min = -5;
int klambda_max = 10;
double klambda_res = 0.25;

double getFakeRate(double pt, double eta) {
    double fakerate = (-8.33753e-03)
                     +((1.48065e-03)*pt)
                     -((3.23176e-05)*pow(pt,2))
                     +((2.91151e-07)*pow(pt,3))
                     -((1.20285e-09)*pow(pt,4));
    if (pt <= 190) {
        fakerate += ((1.88459e-12)*pow(pt,5))*0.25;
    } else {
        fakerate += 0.00058;
    }

    if (std::abs(eta) < 1.4) {
        fakerate *= tauFakeFactor14[tauWP];
    } else if (std::abs(eta) < 2.3) {
        fakerate *= tauFakeFactor23[tauWP];
    } else {
        fakerate *= tauFakeFactor30[tauWP];
    }

    if (useMVATaus) {
        fakerate /= mvaFakeFactorReduction[tauWP];
    }

    if (debug) std::cout << "Light jet with pt eta " << pt << " " << eta << ", fake rate of " << fakerate << "\n";
    return fakerate;
}

std::vector<bool> tagTaus_old(TClonesArray* jets) {
    /*Apply new tau tagging*/
    std::vector<bool> pass;
    Jet* tmpJet;

    for (int i = 0; i < jets->GetEntries(); i++) { //Loop through jets
        tmpJet = (Jet*) jets->At(i);
        if (tmpJet->TauTag) {
            pass.push_back(true);
        } else {
            pass.push_back(false);
        }
    }

    return pass;
}

std::vector<bool> tagTaus(TClonesArray* jets) {
    /*Apply new tau tagging*/
    std::vector<bool> pass;
    Jet* tmpJet;
    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0.0, 1.0);

    for (int i = 0; i < jets->GetEntries(); i++) { //Loop through jets
        tmpJet = (Jet*) jets->At(i);
        double tag = dist(e2);
        if (tmpJet->TauWeight > 0.1) { //Real tau
            if (tag <= tauEfficiency[tauWP]) { //Accept
                pass.push_back(true);
            } else {
                pass.push_back(false);
            }

        } else { //Light jet
            if (tag <= getFakeRate(tmpJet->PT, tmpJet->Eta)) { //Accept
                pass.push_back(true);
            } else {
                pass.push_back(false);
            }
        }
    }

    return pass;
}

std::vector<bool> tag_bjets(TClonesArray* jets, TClonesArray* gen_particles, double dR=0.4) {
    /*Match b-quarks to jets*/
    if (debug) std::cout << "Checking b jets\n";
    std::vector<TLorentzVector> gen_bquarks;
    GenParticle* tmpParticle;
    for (int i=0; i < gen_particles->GetEntries(); i++) {
        tmpParticle = (GenParticle*)gen_particles->At(i);
        if (std::abs(tmpParticle->PID) == 5) gen_bquarks.push_back((TLorentzVector)tmpParticle->P4());
    }
    if (debug) std::cout << gen_bquarks.size() << "  bquarks found\n";    
    std::vector<bool> bjet_real;
    Jet* tmpJet;
    bool real = false;
    for (int i = 0; i < jets->GetEntries(); i++) { //Loop through jets
        tmpJet = (Jet*) jets->At(i);
        real = false;
        for (TLorentzVector bquark : gen_bquarks) {
            if (bquark.DeltaR(tmpJet->P4()) < dR) {
                real = true;
                if (debug) std::cout << "b jet real\n";
                break;
            }
        }
        if (debug) std::cout << "b jet real? " << real << "\n";
        bjet_real.push_back(real);
    }
    return bjet_real;
}

double getMT2(TLorentzVector lepton1_p4, TLorentzVector lepton2_p4,
              TLorentzVector bjet_1, TLorentzVector bjet_2,
              TLorentzVector met_p4) {
    asymm_mt2_lester_bisect::disableCopyrightMessage();
    const double mVisA = bjet_1.M();
    const double pxA = bjet_1.Px();
    const double pyA = bjet_1.Py();
    const double mVisB = bjet_2.M();
    const double pxB = bjet_2.Px();
    const double pyB = bjet_2.Py();
    const double pxMiss = lepton1_p4.Px() + lepton2_p4.Px() + met_p4.Px();
    const double pyMiss = lepton1_p4.Py() + lepton2_p4.Py() + met_p4.Py();
    double chiA = lepton1_p4.M(); // hypothesised mass of invisible on side A.  Must be >=0.
    double chiB = lepton2_p4.M(); // hypothesised mass of invisible on side B.  Must be >=0.
    double MT2 =  asymm_mt2_lester_bisect::get_mT2(mVisA, pxA, pyA,mVisB, pxB, pyB,pxMiss, pyMiss,chiA, chiB,0);
    return MT2;
}

bool getOSTauTauPair(TClonesArray* jets, std::vector<int>* taus, int* tau_0, int* tau_1) {
    /*Checks whether an OS tau-tau pair exists, and if so returns true and points tau and lepton to
    the selected particles*/
    Jet *tau0, *tau1;
    std::vector<std::pair<int, int> > pairs; //Initialise array for OS tau pairs
    if (debug) std::cout << taus->size() << " tau jets found\n";
    for (int t0 : *taus) { //Loop through taus
        for (int t1 : *taus) {
            if (t0 == t1) continue;
            tau0 = (Jet*)jets->At(t0);
            tau1 = (Jet*)jets->At(t1);
            if (tau0->Charge != tau1->Charge) {
                pairs.push_back(std::make_pair(t0, t1));
            }
        }
    }
    if (debug) std::cout << pairs.size() << " OS tau-jet pairs found\n";
    if (pairs.size() == 1) { //Only one OS pair found
        tau0 = (Jet*)jets->At(pairs[0].first);
        tau1 = (Jet*)jets->At(pairs[0].second);
        if (tau0->PT < tau1->PT) { //Order taus by pT
            *tau_1 = pairs[0].first;
            *tau_0 = pairs[0].second;
        }
        return true;
    } else if (pairs.size() > 1) { //Multiple OS pairs: select highest summed |pT| pair
        double pT, highestPT = 0;
        std::pair<int, int> best;
        for (std::pair<int, int> p : pairs) { //Loop through pairs
            tau0 = (Jet*)jets->At(p.first);
            tau1 = (Jet*)jets->At(p.second);
            pT = tau0->PT+tau1->PT;
            if (pT > highestPT) {
                best = p;
                highestPT = pT;
            }
        }
        *tau_0 = best.first;
        *tau_1 = best.second;
        tau0 = (Jet*)jets->At(best.first);
        tau1 = (Jet*)jets->At(best.second);
        if (tau0->PT < tau1->PT) { //Order taus by pT
            *tau_1 = best.first;
            *tau_0 = best.second;
        }
        return true;
    } else { //No OS pairs found
        return false;
    }
    return false;
}

TLorentzVector getHiggs2Taus(MissingET* mpt, TLorentzVector t_0, TLorentzVector t_1) {
    /*Returns 4-vector of Higgs->tau tau*/
    TLorentzVector higgs, mPT;
    mPT.SetPtEtaPhiM(mpt->MET, 0.0, mpt->Phi, 0.0); //TODO Check this
    higgs = t_0 + t_1 + mPT;
    return higgs;
}

bool selectBJets(TClonesArray* jets, std::vector<int>* bJets, int* bJet_0, int* bJet_1) {
    /*Checks is a pair of b-jets exists, returning true if so and pointing bJet_0 and bJet_1 to
    selected jets. Selects pair of jets invariant mass closest to 111 GeV*/
    Jet *jet0, *jet1;
    if (bJets->size() == 2) { //Only two b jets found
        *bJet_0 = (*bJets)[0];
        *bJet_1 = (*bJets)[1];
        jet0 = (Jet*)jets->At(0);
        jet1 = (Jet*)jets->At(1);
        if (jet0->PT < jet1->PT) {
            *bJet_1 = (*bJets)[0];
            *bJet_0 = (*bJets)[1];
        }
        return true;
    } else if (bJets->size() > 2) { //More than two b jets: select pair with invariant mass closest to 111 GeV
        double deltaMin = -1;
        double delta;
        TLorentzVector jet_combined;
        int iMin, jMin;
        for (int i : *bJets) {
            jet0 = (Jet*)jets->At(i);
            for (int j : *bJets) {
                if (i == j) continue;
                jet1 = (Jet*)jets->At(j);
                jet_combined = jet1->P4() + jet0->P4();
                delta = std::abs(111-jet_combined.M());
                if (deltaMin > delta || deltaMin < 0) {
                    deltaMin = delta;
                    iMin = i;
                    jMin = j;
                }
            }
        }
        *bJet_0 = iMin;
        *bJet_1 = jMin;
        jet0 = (Jet*)jets->At(*bJet_0);
        jet1 = (Jet*)jets->At(*bJet_1);
        if (jet0->PT < jet1->PT) {
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
    if (debug) std::cout << "Shapes calculated\n";
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
    TClonesArray *branchMuon, *branchJet, *branchMissingET;
    if (input.find("13Te") != std::string::npos) {
        branchMuon = treeReader->UseBranch("Muon");
        branchJet = treeReader->UseBranch("Jet");
        branchMissingET = treeReader->UseBranch("MissingET");
    } else {
        branchMuon = treeReader->UseBranch("MuonLoose");
        branchJet = treeReader->UseBranch("JetPUPPI");
        branchMissingET = treeReader->UseBranch("PuppiMissingET");
    }
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
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
        if (jet->TauTag) *nTauJets += 1; // TODO: update to new skimming
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
    if (debug) std::cout << "Primary event shapes calculated\n";
    //___________________________________________
}


int moveToEnd(int p, TClonesArray* particles) {
    /*Selects particle at end of 'decay' chain in event record*/
    return p;
    GenParticle* mother = (GenParticle*)particles->At(p);
    while (((GenParticle*)particles->At(mother->D1))->PID == mother->PID &&
            ((GenParticle*)particles->At(mother->D2))->PID == mother->PID) {
        p = mother->D1;
        mother = (GenParticle*)particles->At(p);
    }
    return p;
}

bool correctDecayChannel(std::string input, Long64_t cEvent,
        std::map<std::string, TH1D*>* plots=NULL, int* hBB=NULL, int* hTauTau=NULL) {
    /*Make sure event is hh->bbtautau, and point hbb and htautau to the Higgs*/
    TChain *chain = new TChain("Delphes");
    chain->Add(input.c_str());
    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    treeReader->ReadEntry(cEvent);
    bool hBBFound = false, hTauTauFound = false;
    int nHiggs = 0;
    if (plots != NULL) (*plots)["cuts"]->Fill("hh->bb#tau#tau check", 1);
    for (int p = 0; p < branchParticle->GetEntriesFast(); ++p) {
        if (((GenParticle*)branchParticle->At(p))->PID == 25) { // && ((GenParticle*)branchParticle->At(p))->Status == 22) { //Particle is Higgs
            if (((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D1))->PID != 25 &&
                ((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D2))->PID != 25) {
                nHiggs++;
                if (plots != NULL) (*plots)["higgsDecay"]->Fill(std::abs(((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D1))->PID));
                if (plots != NULL) (*plots)["higgsDecay"]->Fill(std::abs(((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D2))->PID));
                if (std::abs(((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D1))->PID) == 5
                        && std::abs(((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D2))->PID) == 5) { //Daughters are b quarks
                    hBBFound = true;
                    if (hBB != NULL) *hBB = p; //Point to Higgs
                    if (hBBFound && hTauTauFound) { //h->bb and h->tautau found, so accept event
                        if (plots != NULL) (*plots)["cuts"]->Fill("hh->bb#tau#tau pass", 1);
                        chain->Delete();
                        delete treeReader;
                        return true;
                    }
                }
                if (std::abs(((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D1))->PID) == 15
                        && std::abs(((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D2))->PID) == 15) { //Daughters are taus
                    hTauTauFound = true;
                    if (hTauTau != NULL) *hTauTau = p; //Point to Higgs
                    if (hBBFound && hTauTauFound) { //h->bb and h->tautau found, so accept event
                        if (plots != NULL) (*plots)["cuts"]->Fill("hh->bb#tau#tau pass", 1);
                        chain->Delete();
                        delete treeReader;
                        return true;
                    }
                }
            }
            if (nHiggs >= 2) break; //Both Higgs found
        }
    }
    chain->Delete();
    delete treeReader;
    return false; //Both h->bb and h->tautau not found
}

bool getGenSystem(TClonesArray* branchParticle,
    TLorentzVector* v_gen_higgs_bb, TLorentzVector* v_gen_higgs_tt,
    TLorentzVector* v_gen_tau_0, TLorentzVector* v_gen_tau_1,
    TLorentzVector* v_gen_bJet_0, TLorentzVector* v_gen_bJet_1) {
    /*Simplified version of correctDecayChannel + truthCut designed for filtered Delphes, where mother-daughter links are broken*/
    GenParticle* tmpParticle;
    bool match = false;

    std::vector<TLorentzVector> higgs, bquarks, taus, anti_bquarks, anti_taus;
    for (int p = 0; p < branchParticle->GetEntriesFast(); ++p) {
        tmpParticle = (GenParticle*)branchParticle->At(p);
        if (tmpParticle->IsPU == true) continue;
        if (tmpParticle->Status == 22 && tmpParticle->PID == 25) {
            higgs.push_back(tmpParticle->P4());
        } // else if (tmpParticle->Status == 23) {
        //     if (tmpParticle->PID == 5) {//23
        //         bquarks.push_back(tmpParticle->P4());
        //     } else if (tmpParticle->PID == -5) {
        //         anti_bquarks.push_back(tmpParticle->P4());
        //     } else if (tmpParticle->PID == 15) {
        //         taus.push_back(tmpParticle->P4());
        //     } else if (tmpParticle->PID == -15) {
        //         anti_taus.push_back(tmpParticle->P4());
        //     }
        // }
    }
    if (higgs.size() != 2) {
        throw std::runtime_error("Signal doesn't contain exactly 2 status 22 Higgs");
        return false;
    }

    *v_gen_higgs_tt = higgs[0];
    *v_gen_higgs_bb = higgs[1];
    return false;

    //Doesn't work, status 22 appears to be purely rest-frame of di-Higgs
    // if (taus.size() >= 1 && anti_taus.size() >= 1 && bquarks.size() >= 1 && anti_bquarks.size() >= 1) {
    //     match = true;
    //     double min_dR = 999, dR;
    //     TLorentzVector sum;
    //     int htt = 0;

    // //Match taus to Higgs
    //     for (TLorentzVector i : taus) {
    //         for (TLorentzVector j : anti_taus) {
    //             sum = i+j;

    //             for (int k = 0; k < 2; k++) {
    //                 dR = sum.DeltaR(higgs[k]);
    //                 if (dR < min_dR) {
    //                     min_dR = dR;
    //                     *v_gen_higgs_tt = higgs[k];
    //                     htt = k;
    //                     if (i.Pt() > j.Pt()) {
    //                         *v_gen_tau_0 = i;
    //                         *v_gen_tau_1 = j;
    //                     } else {
    //                         *v_gen_tau_0 = j;
    //                         *v_gen_tau_1 = i;
    //                     }
    //                 }
    //             }
    //         }
    //     }

    //     //Match bquarks to Higgs
    //     min_dR = 999;
    //     for (TLorentzVector i : bquarks) {
    //         for (TLorentzVector j : anti_bquarks) {
    //             sum = i+j;

    //             dR = sum.DeltaR(higgs[1-htt]);
    //             if (dR < min_dR) {
    //                 min_dR = dR;
    //                 *v_gen_higgs_bb = higgs[1-htt];
    //                 if (i.Pt() > j.Pt()) {
    //                     *v_gen_bJet_0 = i;
    //                     *v_gen_bJet_1 = j;
    //                 } else {
    //                     *v_gen_bJet_0 = j;
    //                     *v_gen_bJet_1 = i;
    //                 }
    //             }
    //         }
    //     }
    // } else { //Gen-level final-states not all found
    //     *v_gen_higgs_tt = higgs[0];
    //     *v_gen_higgs_bb = higgs[1];
    // }
    // return match;
}
    

bool checkDiJet(TClonesArray* jets, TClonesArray* particles, int j_0, int j_1, int mother, int pID,
        int* swap, TH1D* dRPlot, double R) {
    /*Checks whether the particles are within their nearest jet*/
    //Associate particles to closest found jet___
    int p_0 = -1, p_1 = -1;
    Jet* jet_0 = (Jet*)jets->At(j_0);
    Jet* jet_1 = (Jet*)jets->At(j_1);
    GenParticle* higgs = (GenParticle*)particles->At(mother);
    if (std::abs(((GenParticle*)particles->At(moveToEnd(higgs->D1, particles)))->PID) != pID) { //Make sure decays products are correct
        std::cout << "Something's gone wrong in h->" + doubleToString(pID) + " -" + doubleToString(pID) + "!\n";
        return false;
    }
    if (((GenParticle*)particles->At(higgs->D1))->P4().DeltaR(jet_0->P4()) <
        ((GenParticle*)particles->At(higgs->D1))->P4().DeltaR(jet_1->P4())) {
        p_0 = higgs->D1;
        p_1 = higgs->D2;
        if (swap != NULL) *swap = 0;
    } else {
        p_0 = higgs->D2;
        p_1 = higgs->D1;
        if (swap != NULL) *swap = 1;
    }
    //___________________________________________
    //Check jets_________________________________
    double dR_0 = ((GenParticle*)particles->At(p_0))->P4().DeltaR(jet_0->P4());
    double dR_1 = ((GenParticle*)particles->At(p_1))->P4().DeltaR(jet_1->P4());
    if (dR_0 > R || dR_1 > R) { //particle(s) outside jet
        return false;
    }
    //___________________________________________
    //Accept association and fill plot___________
    dRPlot->Fill(dR_0);
    dRPlot->Fill(dR_1);
    //___________________________________________
    return true;
}

int ancestrySearch(GenParticle* child, GenParticle* parent_0, GenParticle* parent_1, TClonesArray* particles) {
    /*Recursive search through child's ancestry for parent 0 or 1. If found returns 0 or 1. If not found returns -1*/
    int ancestor = -1;
    if (debug) std::cout << "size: " << particles->GetEntries() << ", mothers: " << child->M1  << " " << child->M2 << "\n";
    GenParticle* mother;
    if (child->M1 > 0) {// && child->M1 < particles->GetEntries()) { //Check first mother
        mother = (GenParticle*)particles->At(child->M1);
        if (mother->GetUniqueID() == parent_0->GetUniqueID()) {
            return 0;
        } else if (mother->GetUniqueID() == parent_1->GetUniqueID()) {
            return 1;
        } else {
            ancestor = ancestrySearch(mother, parent_0, parent_1, particles); //Recursive search
        }
    }
    if (ancestor == -1 && child->M2 > 0) {// && child->M2 < particles->GetEntries()) { //Check second mother
        mother = (GenParticle*)particles->At(child->M2);
        if (mother->GetUniqueID() == parent_0->GetUniqueID()) {
            return 0;
        } else if (mother->GetUniqueID() == parent_1->GetUniqueID()) {
            return 1;
        } else {
            ancestor = ancestrySearch(mother, parent_0, parent_1, particles); //Recursive search
        }
    }
    return ancestor;
}

std::string typeLookup(std::vector<std::string> options) {
    /*Lookup for histogram bin name*/
    if (options[0] == "electron" && options[1] == "electron") return "ee";
    if (options[0] == "muon" && options[1] == "muon") return "#mu#mu";
    if (options[0] == "tau" && options[1] == "tau") return "#tau_{h}#tau_{h}";
    if (options[0] == "electron" && options[1] == "muon") return "e#mu";
    if (options[0] == "tau" && options[1] == "electron") return "e#tau_{h}";
    if (options[0] == "tau" && options[1] == "muon") return "#mu#tau_{h}";
    return "error";
}

bool truthCut(std::string input, Long64_t cEvent, int b_0, int b_1, int l_0, int l_1, int hBB, int hTauTau,
        std::vector<std::string> options, std::map<std::string, TH1D*>* plots,
        TLorentzVector* v_gen_higgs_bb, TLorentzVector* v_gen_higgs_tt,
        TLorentzVector* v_gen_tau_0, TLorentzVector* v_gen_tau_1,
        TLorentzVector* v_gen_bJet_0, TLorentzVector* v_gen_bJet_1) {
    /*Checks whether selected final states are correct*/
    double jetRadius = 0.5;
    TChain *chain = new TChain("Delphes");
    chain->Add(input.c_str());
    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchMuon, *branchJet;
    if (input.find("13Te") != std::string::npos) {
        branchMuon = treeReader->UseBranch("Muon");
        branchJet = treeReader->UseBranch("Jet");
    } else {
        branchMuon = treeReader->UseBranch("MuonLoose");
        branchJet = treeReader->UseBranch("JetPUPPI");
    }
    treeReader->ReadEntry(cEvent);
    if (debug) std::cout << "Loading data for MC truth cut on event mode " << options[0] << ":" << options[1] << "\n";
    //Check if selected final states are correct_
    (*plots)["cuts"]->Fill("MC-truth check", 1);
    int swap;
    //Check b jets_______________________________
    GenParticle *bJet_0, *bJet_1;
    GenParticle* higgs = (GenParticle*)branchParticle->At(hBB);
    (*plots)["cuts"]->Fill("b-jets check", 1);
    if (debug) std::cout << "Checking b-jets\n";
    if (!checkDiJet(branchJet, branchParticle, b_0, b_1, hBB, 5, &swap, (*plots)["bMatch"], jetRadius)) {
        if (debug) std::cout << "MC check fails due to di-Jet on b-jets check\n";
        chain->Delete();
        return false; //b-jet selection incorrect
    }
    if (debug) std::cout << "Both b jets confirmed\n";
    (*plots)["cuts"]->Fill("b-jets pass", 1);
    if (swap) {
        bJet_0 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D2, branchParticle));
        bJet_1 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D1, branchParticle));
    } else {
        bJet_0 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D1, branchParticle));
        bJet_1 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D2, branchParticle));
    }
    //___________________________________________
    //Check taus_________________________________
    if (debug) std::cout << "Checking taus\n";
    GenParticle *tau_0, *tau_1;
    higgs = (GenParticle*)branchParticle->At(hTauTau);
    if (options[0] == "tau" && options[1] == "tau") {
        //h->tau_h tau_h_________________________
        if (!checkDiJet(branchJet, branchParticle, l_0, l_1, hTauTau, 15, &swap, (*plots)["tauMatch"], jetRadius)) {
            if (debug) std::cout << "MC check fails due to di-Jet on tau-jets check\n";
            chain->Delete();
            return false; //tau-jet selection incorrect
        }
        if (swap) {
            tau_0 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D2, branchParticle));
            tau_1 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D1, branchParticle));
        } else {
            tau_0 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D1, branchParticle));
            tau_1 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D2, branchParticle));
        }
        (*plots)["cuts"]->Fill(("h->#tau#tau->" + typeLookup(options) + " pass").c_str(), 1);
        //_______________________________________
    } else if ((options[0] == "tau" && options[1] != "tau") || (options[0] != "tau" && options[1] == "tau")) {
        //h->tau_h light-lepton__________________
        //Load objects___________________________
        tau_0 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D1, branchParticle));
        tau_1 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D2, branchParticle));
        GenParticle* lightLepton;
        Jet* tauJet;
        for (int i = 0; i < 2; i++) {
            int l = l_0;
            if (i == 1) {
                l = l_1;
            }
            if (options[i] == "tau") {
                tauJet = (Jet*)branchJet->At(l);
            } else if (options[i] == "muon") {
                lightLepton = (GenParticle*)((Muon*)branchMuon->At(l))->Particle.GetObject();
            } else if (options[i] == "electron") {
                lightLepton = (GenParticle*)((Electron*)branchElectron->At(l))->Particle.GetObject();
            }
        }
        if (lightLepton) {
            if (debug) std::cout << "MC check fails due to light lepton being PU\n";
            return false;
        }
        //_______________________________________
        //Check taus_____________________________
        int leptonMother = ancestrySearch(lightLepton, tau_0, tau_1, branchParticle);
        if (leptonMother == -1) {
            if (debug) std::cout << "MC check fails due to ancestry check\n";
            chain->Delete();
            return false; //Light lepton did not come from tau decay
        }
        GenParticle* tauh;
        if (leptonMother == 0) {
            tauh = tau_1;
            tau_1 = tau_0; //Reassociate 0 to tau and 1 to lepton
            tau_0 = tauh;
        } else {
            tauh = tau_0;
        }
        if (tauh->P4().DeltaR(tauJet->P4()) > jetRadius) {
            if (debug) std::cout << "MC check fails due to tau-jet check\n";
            chain->Delete();
            return false; //Tau outside selected jet
        }
        (*plots)["cuts"]->Fill(("h->#tau#tau->" + typeLookup(options) + " pass").c_str(), 1);
        //_______________________________________
        //_______________________________________
    } else {
        //h->light-lepton light-lepton___________
        //Load objects___________________________
        GenParticle* higgs = (GenParticle*)branchParticle->At(hTauTau);
        tau_0 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D1, branchParticle));
        tau_1 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D2, branchParticle));
        GenParticle *lightLepton_0, *lightLepton_1;
        if (options[0] == "muon") {
            lightLepton_0 = (GenParticle*)((Muon*)branchMuon->At(l_0))->Particle.GetObject();
        } else if (options[0] == "electron") {
            lightLepton_0 = (GenParticle*)((Electron*)branchElectron->At(l_0))->Particle.GetObject();
        }
        if (options[1] == "muon") {
            lightLepton_1 = (GenParticle*)((Muon*)branchMuon->At(l_1))->Particle.GetObject();
        } else if (options[1] == "electron") {
            lightLepton_1 = (GenParticle*)((Electron*)branchElectron->At(l_1))->Particle.GetObject();
        }
        if (lightLepton_0->IsPU || lightLepton_1->IsPU) {
            if (debug) std::cout << "MC check fails due to light lepton being PU\n";
            return false;
        }
        //_______________________________________
        //Check taus_____________________________
        int leptonMother_0 = ancestrySearch(lightLepton_0, tau_0, tau_1, branchParticle);
        if (leptonMother_0 == -1) {
            if (debug) std::cout << "MC check fails due to ancestry check\n";
            chain->Delete();
            return false; //Light lepton 0 did not come from tau decay
        }
        int leptonMother_1 = ancestrySearch(lightLepton_1, tau_0, tau_1, branchParticle);
        if (leptonMother_1 == -1) {
            if (debug) std::cout << "MC check fails due to ancestry check\n";
            chain->Delete();
            return false; //Light lepton 1 did not come from tau decay
        }
        if (leptonMother_0 == leptonMother_1) {
            if (debug) std::cout << "MC check fails due to both leptons coming from same tau\n";
            chain->Delete();
            return false; //Leptons both came from same mother (somehow)
        }
        (*plots)["cuts"]->Fill(("h->#tau#tau->" + typeLookup(options) + " pass").c_str(), 1);
        if ((lightLepton_0->PT > lightLepton_1->PT & leptonMother_0 == 1) |
                (lightLepton_0->PT < lightLepton_1->PT & leptonMother_0 == 0)) {
            tau_0 = tau_1;
            tau_1 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D1, branchParticle));
        }
        //_______________________________________
        //_______________________________________
    }
    if (debug) std::cout << "Both taus confirmed\n";
    //___________________________________________
    if (debug) std::cout << "Event accepted\n";
    (*plots)["cuts"]->Fill("#taus pass", 1);
    (*plots)["cuts"]->Fill("MC-truth pass", 1);
    //Get vectors for regression_________________
    *v_gen_higgs_bb = ((GenParticle*)branchParticle->At(hBB))->P4();
    *v_gen_higgs_tt = ((GenParticle*)branchParticle->At(hTauTau))->P4();
    *v_gen_tau_0 = tau_0->P4();
    *v_gen_tau_1 = tau_1->P4();
    *v_gen_bJet_0 = bJet_0->P4();
    *v_gen_bJet_1 = bJet_1->P4();
    //___________________________________________
    chain->Delete();
    delete treeReader;
    return true;
}

double functionGF(double kl, double kt, double c2, double cg, double c2g, std::array<double, n_Ai_coeffs> const &A)
{
    // this can be extended to 5D coefficients; currently c2, cg, c2g are unused
    // return ( A1*pow(kt,4) + A3*pow(kt,2)*pow(kl,2) + A7*kl*pow(kt,3) );
    return ( A[0]*pow(kt,4) + A[1]*pow(c2,2) + (A[2]*pow(kt,2) + A[3]*pow(cg,2))*pow(kl,2) + A[4]*pow(c2g,2) + ( A[5]*c2 + A[6]*kt*kl )*pow(kt,2) + (A[7]*kt*kl + A[8]*cg*kl )*c2 + A[9]*c2*c2g + (A[10]*cg*kl + A[11]*c2g)*pow(kt,2)+ (A[12]*kl*cg + A[13]*c2g )*kt*kl + A[14]*cg*c2g*kl );
}

std::pair<int,int> get_bin_idx(double x, double y, TH2* histo)
{
    int ix = histo->GetXaxis()->FindBin(x);
    int iy = histo->GetYaxis()->FindBin(y);
    if (ix < 1) ix = 1;
    if (iy < 1) iy = 1;
    if (ix > histo->GetNbinsX()) ix = histo->GetNbinsX();
    if (iy > histo->GetNbinsY()) iy = histo->GetNbinsY();
    return std::make_pair(ix,iy);
}

int main(int argc, char *argv[]) { //input, output, N events, truth
    std::map<std::string, std::string> options = getOptions(argc, argv);
    if (options.size() == 0) {
        return 1;
    }
    //ROOT settings______________________________
    gSystem->Load("libDelphes.so");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadGridX(kFALSE);
    gStyle->SetPadGridY(kFALSE);
    //___________________________________________
    // read out the Ai histogram coefficients____
    TFile* f_14TeV_coeffs = new TFile("Coefficients_14TeV.root");
    std::array<TH2*,n_Ai_coeffs> histos_A;
    for (uint idx = 0; idx < n_Ai_coeffs; ++idx)
    {
        int Ai = idx + 1; // there is a naming shift, first coefficient is A1, second is A2 etc...
        std::string hName = std::string("A") + std::to_string(Ai) + std::string("_14TeV");
        TH2D* h = (TH2D*) f_14TeV_coeffs->Get(hName.c_str());
        histos_A.at(idx) = h;
    }
    TFile* f_HH_14TeV_histo = new TFile("HH_SM_2D_histo.root"); // done from the same tree as the HH_ME_info.root, but kept here for an easy usage
    TH2* HH_14TeV_histo = (TH2*) f_HH_14TeV_histo->Get("h_events_SM");
    //Build output files_________________________
    std::string outputName(options["-o"]);
    makeDirs(outputName);
    TFile* outputFile = new TFile((outputName + ".root").c_str(), "recreate");
    outputFile->cd();
    //Initialise variables_______________________
    std::cout << "Initialising variables\n";
    int lepton_0, lepton_1, tau_0, tau_1, bJet_0, bJet_1;
    //Low-level variables________________________
    double t_0_pT, t_0_eta, t_0_phi, t_0_mass, t_0_mT; //Tau 0 variables
    double t_1_pT, t_1_eta, t_1_phi, t_1_mass, t_1_mT; //Tau 1 variables
    double b_0_pT, b_0_eta, b_0_phi, b_0_mass; //b-jet 0 variables
    double b_1_pT, b_1_eta, b_1_phi, b_1_mass; //b-jet 1 variables
    double mPT_pT, mPT_phi; //Missing ET variables
    //___________________________________________
    //Reconstructed variables____________________
    double h_tt_pT, h_tt_eta, h_tt_phi, h_tt_mass; //Higgs 0 variables
    double h_bb_pT, h_bb_eta, h_bb_phi, h_bb_mass; //Higgs 1 variables
    double diH_pT, diH_eta, diH_phi, diH_mass, diH_mT2; //di-Higgs variables
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
    bool gen_mctMatch; //MC truth matc
    bool gen_t_0_real, gen_t_1_real, gen_b_0_real, gen_b_1_real; //Real/fake jets
    double gen_cosThetaStar;
    //___________________________________________
    //klambda reweighting________________________
    std::vector<double> gen_weight_klambda;
    //___________________________________________
    double weight; //Event weight
    int hBB = -1, hTauTau = -1;
    int nElectrons = 0, nMuons = 0;
    bool eventAccepted = false;
    TTree* e_tau_b_b = new TTree("e_tau_b_b", "e #tau b #bar{b}");
    e_tau_b_b->Branch("t_0_pT", &t_0_pT);
    e_tau_b_b->Branch("t_0_eta", &t_0_eta);
    e_tau_b_b->Branch("t_0_phi", &t_0_phi);
    e_tau_b_b->Branch("t_0_mass", &t_0_mass);
    e_tau_b_b->Branch("t_0_mT", &t_0_mT);
    e_tau_b_b->Branch("t_1_pT", &t_1_pT);
    e_tau_b_b->Branch("t_1_eta", &t_1_eta);
    e_tau_b_b->Branch("t_1_phi", &t_1_phi);
    e_tau_b_b->Branch("t_1_mass", &t_1_mass);
    e_tau_b_b->Branch("t_1_mT", &t_1_mT);
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
    e_tau_b_b->Branch("diH_mT2", &diH_mT2);
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
    e_tau_b_b->Branch("gen_cosThetaStar", &gen_cosThetaStar);
    e_tau_b_b->Branch("gen_weight", &weight);
    e_tau_b_b->Branch("gen_weight_klambda", &gen_weight_klambda);
    e_tau_b_b->Branch("gen_t_0_real", &gen_t_0_real);
    e_tau_b_b->Branch("gen_t_1_real", &gen_t_1_real);
    e_tau_b_b->Branch("gen_b_0_real", &gen_b_0_real);
    e_tau_b_b->Branch("gen_b_1_real", &gen_b_1_real);
    TTree* mu_tau_b_b = new TTree("mu_tau_b_b", "#mu #tau_{h} b #bar{b}");
    mu_tau_b_b->Branch("t_0_pT", &t_0_pT);
    mu_tau_b_b->Branch("t_0_eta", &t_0_eta);
    mu_tau_b_b->Branch("t_0_phi", &t_0_phi);
    mu_tau_b_b->Branch("t_0_mass", &t_0_mass);
    mu_tau_b_b->Branch("t_0_mT", &t_0_mT);
    mu_tau_b_b->Branch("t_1_pT", &t_1_pT);
    mu_tau_b_b->Branch("t_1_eta", &t_1_eta);
    mu_tau_b_b->Branch("t_1_phi", &t_1_phi);
    mu_tau_b_b->Branch("t_1_mass", &t_1_mass);
    mu_tau_b_b->Branch("t_1_mT", &t_1_mT);
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
    mu_tau_b_b->Branch("diH_mT2", &diH_mT2);
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
    mu_tau_b_b->Branch("gen_cosThetaStar", &gen_cosThetaStar);
    mu_tau_b_b->Branch("gen_weight", &weight);
    mu_tau_b_b->Branch("gen_weight_klambda", &gen_weight_klambda);
    mu_tau_b_b->Branch("gen_t_0_real", &gen_t_0_real);
    mu_tau_b_b->Branch("gen_t_1_real", &gen_t_1_real);
    mu_tau_b_b->Branch("gen_b_0_real", &gen_b_0_real);
    mu_tau_b_b->Branch("gen_b_1_real", &gen_b_1_real);
    TTree* tau_tau_b_b = new TTree("tau_tau_b_b", "#tau_{h} #tau_{h} b #bar{b}");
    tau_tau_b_b->Branch("t_0_pT", &t_0_pT);
    tau_tau_b_b->Branch("t_0_eta", &t_0_eta);
    tau_tau_b_b->Branch("t_0_phi", &t_0_phi);
    tau_tau_b_b->Branch("t_0_mass", &t_0_mass);
    tau_tau_b_b->Branch("t_0_mT", &t_0_mT);
    tau_tau_b_b->Branch("t_1_pT", &t_1_pT);
    tau_tau_b_b->Branch("t_1_eta", &t_1_eta);
    tau_tau_b_b->Branch("t_1_phi", &t_1_phi);
    tau_tau_b_b->Branch("t_1_mass", &t_1_mass);
    tau_tau_b_b->Branch("t_1_mT", &t_1_mT);
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
    tau_tau_b_b->Branch("diH_mT2", &diH_mT2);
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
    tau_tau_b_b->Branch("gen_cosThetaStar", &gen_cosThetaStar);
    tau_tau_b_b->Branch("gen_weight", &weight);
    tau_tau_b_b->Branch("gen_weight_klambda", &gen_weight_klambda);
    tau_tau_b_b->Branch("gen_t_0_real", &gen_t_0_real);
    tau_tau_b_b->Branch("gen_t_1_real", &gen_t_1_real);
    tau_tau_b_b->Branch("gen_b_0_real", &gen_b_0_real);
    tau_tau_b_b->Branch("gen_b_1_real", &gen_b_1_real);
    std::cout << "Variables initialised\n";
    //___________________________________________
    //Initialise plots___________________________
    std::cout << "Initialising plot\n";
    TH1D* h_datasetSizes = new TH1D("Dataset_sizes", "Dataset sizes", 4, -0.4, 0.4);
    TH1D* h_sum_w = new TH1D("Sum_w", "Sum klambda weights / SM_Evts", 1+((klambda_max-klambda_min)/klambda_res), klambda_min, klambda_max+klambda_res);
    TH1D* h_e_tau_b_b_cutFlow;
    TH1D* h_mu_tau_b_b_cutFlow;
    TH1D* h_tau_tau_b_b_cutFlow;
    std::map<std::string, TH1D*> mcTruthPlots;
    mcTruthPlots.insert(std::make_pair("cuts", new TH1D("mcTruth_cutFlow", "MC Truth Cuts", 20, -2.0, 2.0)));
    mcTruthPlots.insert(std::make_pair("bMatch", new TH1D("mcTruth_bJetMatching", "#DeltaR(b, jet)", 50, 0.0, 0.5)));
    mcTruthPlots.insert(std::make_pair("tauMatch", new TH1D("mcTruth_tauJetMatching", "#DeltaR(#tau, jet)", 50, 0.0, 0.5)));
    mcTruthPlots.insert(std::make_pair("higgsDecay", new TH1D("mcTruth_higgsDecay", "Higgs product |PID|", 50, 0, 50)));
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(1, "hh->bb#tau#tau check");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(2, "hh->bb#tau#tau pass");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(3, "MC-truth check");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(4, "MC-truth pass");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(5, "b-jets check");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(6, "b-jets pass");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(7, "#taus check");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(8, "#taus pass");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(9, "h->#tau#tau->ee check");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(10, "h->#tau#tau->ee pass");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(11, "h->#tau#tau->e#mu check");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(12, "h->#tau#tau->e#mu pass");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(13, "h->#tau#tau->#mu#mu check");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(14, "h->#tau#tau->#mu#mu pass");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(15, "h->#tau#tau->e#tau_{h} check");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(16, "h->#tau#tau->e#tau_{h} pass");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(17, "h->#tau#tau->#mu#tau_{h} check");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(18, "h->#tau#tau->#mu#tau_{h} pass");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(19, "h->#tau#tau->#tau_{h}#tau_{h} check");
    mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(20, "h->#tau#tau->#tau_{h}#tau_{h} pass");
    h_e_tau_b_b_cutFlow = new TH1D("e_tau_b_b_Cut_Flow", "e #tau_{h} b #bar{b} cut flow", 5, -0.5, 0.5);
    h_mu_tau_b_b_cutFlow = new TH1D("mu_tau_b_b_Cut_Flow", "#mu #tau_{h} b #bar{b} cut flow", 5, -0.5, 0.5);
    h_tau_tau_b_b_cutFlow = new TH1D("tau_tau_b_b_Cut_Flow", "#tau_{h} #tau_{h} b #bar{b} cut flow", 5, -0.5, 0.5);
    h_datasetSizes->GetXaxis()->SetBinLabel(1, "All");
    h_datasetSizes->GetXaxis()->SetBinLabel(2, "#mu #tau_{h} b #bar{b}");
    h_datasetSizes->GetXaxis()->SetBinLabel(3, "e #tau_{h} b #bar{b}");
    h_datasetSizes->GetXaxis()->SetBinLabel(4, "#tau_{h} #tau_{h} b #bar{b}");
    h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(1, "All");
    h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(2, "Quality e");
    h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(3, "1 e & 0 #mu");
    h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(4, "Quality #tau");
    h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(5, "Quality b#bar{b}");
    h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(1, "All");
    h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(2, "Quality #mu");
    h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(3, "1 #mu & 0 e");
    h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(4, "Quality #tau");
    h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(5, "Quality b#bar{b}");
    h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(1, "All");
    h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(2, "0 e & 0 #mu");
    h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(3, "Quality #tau#tau");
    h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(4, "OS");
    h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(5, "Quality b#bar{b}");

    std::cout << "Plots initialised\n";
    //___________________________________________
    //Load data__________________________________
    std::cout << "Running event selection\n";
    TChain *chain = new TChain("Delphes");
    chain->Add(options["-i"].c_str());
    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    TClonesArray *branchMuon, *branchJet, *branchMissingET;
    if (options["-i"].find("13Te") != std::string::npos) {
        branchMuon = treeReader->UseBranch("Muon");
        branchJet = treeReader->UseBranch("Jet");
        branchMissingET = treeReader->UseBranch("MissingET");
    } else {
        branchMuon = treeReader->UseBranch("MuonLoose");
        branchJet = treeReader->UseBranch("JetPUPPI");
        branchMissingET = treeReader->UseBranch("PuppiMissingET");
    }
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");   
    TClonesArray *branchWeights = treeReader->UseBranch("Weight");
    TClonesArray *branchParticle  = treeReader->UseBranch("Particle");
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
    MissingET* tmpMPT;
    Weight* tmpWeight;
    std::vector<bool> tauTags, bjets_real;


    std::cout << "Beginning event loop\n";
    for (Long64_t cEvent = 0; cEvent < nEvents; cEvent++) {
        if (debug) std::cout << "Loading event " << cEvent << "\n";
        treeReader->ReadEntry(cEvent); //Load next event
        if (debug) std::cout << "Event loaded, getting data\n";
        if (debug) std::cout << "Data loaded\n";
        if (cEvent % 1000 == 0) std::cout << "Loop: " << cEvent << "/" << nEvents << ", " <<
                100*cEvent/nEvents << "%\n";
        h_datasetSizes->Fill("All", 1);
        eventAccepted = false;
        if (options["-i"].find("13Te") != std::string::npos) {
            tauTags = tagTaus_old(branchJet);
        } else {
            tauTags = tagTaus(branchJet); //get new tau tags
        }
        hBB = -1;
        hTauTau = -1;
        gen_mctMatch = false;
        TLorentzVector v_gen_higgs_bb, v_gen_higgs_tt, v_gen_diHiggs, v_gen_tau_0, v_gen_tau_1, v_gen_bJet_0, v_gen_bJet_1;

        //klambda weights
        if (options["-i"].find("GluGluToHHTo2B2Tau_node_SM_14TeV") != std::string::npos) { //14TeV Signal
            //Gen system
            gen_mctMatch = getGenSystem(branchParticle, &v_gen_higgs_bb, &v_gen_higgs_tt, &v_gen_tau_0, &v_gen_tau_1, &v_gen_bJet_0, &v_gen_bJet_1);
            v_gen_diHiggs = getDiHiggs(v_gen_higgs_tt, v_gen_higgs_bb);
            gen_cosThetaStar = std::abs(v_gen_higgs_bb.CosTheta());
            gen_diH_mass = v_gen_diHiggs.M();

            //Binning
            auto bins = get_bin_idx(gen_diH_mass, std::abs(gen_cosThetaStar), histos_A.at(0));
            std::array<double,n_Ai_coeffs> values_A;
            for (uint idx = 0; idx < n_Ai_coeffs; ++idx)
            {
                values_A.at(idx) = histos_A.at(idx)->GetBinContent(bins.first, bins.second); // read the 2D bin content for these mHH, costh values
                values_A.at(idx) = A_integralXS.at(idx)*values_A.at(idx);
            }

            double w;
            double SM_evts = HH_14TeV_histo->GetBinContent(bins.first, bins.second);
            gen_weight_klambda.clear();
            for (double kl = klambda_min; kl < klambda_max+klambda_res; kl += klambda_res) {
                w = functionGF(kl, 1, 0, 0, 0, values_A)/SM_evts;
                gen_weight_klambda.push_back(w);
                h_sum_w->Fill(kl, w);
            }       
        }

        bjets_real = tag_bjets(branchJet, branchParticle, 0.4);

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
            tmpMuon = (Muon*) branchMuon->At(muons[0]);
            for (int i = 0; i < branchElectron->GetEntries(); i++) { //Loop through electrons
                tmpElectron = (Electron*) branchElectron->At(i);
                if (tmpElectron->PT > ePTMinAdd && std::abs(tmpElectron->Eta) < eEtaMaxAdd
                        && tmpElectron->IsolationVar < eIsoMaxAdd) { //Additional electron
                    addElectron = true;
                    break;
                }
            }
            if (!addElectron) { //No additional electrons found
                h_mu_tau_b_b_cutFlow->Fill("1 #mu & 0 e", 1);
                for (int i = 0; i < branchJet->GetEntries(); i++) { //Loop through jets
                    tmpJet = (Jet*) branchJet->At(i);
                    if (tauTags[i] == true && tmpJet->BTag == 0 && tmpJet->PT > tauPTMin
                            && std::abs(tmpJet->Eta) < tauEtaMax
                            && tmpJet->Charge != tmpMuon->Charge) { //Quality  OS tau
                        taus.push_back(i);
                    }
                    if (options["-i"].find("13Te") != std::string::npos) {
                        if (tauTags[i] == false && tmpJet->BTag == 1 && tmpJet->PT > bJetPTMin
                                && std::abs(tmpJet->Eta) < bJetEtaMax) { //Quality b jet
                            bJets.push_back(i);
                        }
                    } else {
                        if (tauTags[i] == false && (tmpJet->BTag & (1 << 4) ) && tmpJet->PT > bJetPTMin
                                && std::abs(tmpJet->Eta) < bJetEtaMax) { //Quality b jet
                            bJets.push_back(i);
                        }
                    }
                }
                if (taus.size() >= 1) {//Quality tau
                    h_mu_tau_b_b_cutFlow->Fill("Quality #tau", 1);
                    if (bJets.size() >= 2) {//Quality b jets pairs found
                        h_mu_tau_b_b_cutFlow->Fill("Quality b#bar{b}", 1);
                        if (selectBJets(branchJet, &bJets, &bJet_0, &bJet_1) == true) { //Quality b-jet pair found
                            if (debug) std::cout << "Accepting event\n";
                            v_tau_1 = tmpMuon->P4();
                            gen_t_1_real = true;
                            tmpJet = (Jet*)branchJet->At(taus[0]);
                            v_tau_0 = tmpJet->P4();
                            gen_t_0_real = (tmpJet->TauWeight > 0.1) ? true : false;
                            tmpMPT = (MissingET*)branchMissingET->At(0);
                            v_higgs_tt = getHiggs2Taus(tmpMPT, v_tau_0, v_tau_1);
                            tmpJet = (Jet*)branchJet->At(bJet_0);
                            v_bJet_0 = tmpJet->P4();
                            gen_b_0_real = bjets_real[bJet_0];
                            tmpJet = (Jet*)branchJet->At(bJet_1);
                            v_bJet_1 = tmpJet->P4();
                            gen_b_1_real = bjets_real[bJet_1];
                            v_higgs_bb = getHiggs2Bs(v_bJet_0, v_bJet_1);
                            v_diHiggs = getDiHiggs(v_higgs_tt, v_higgs_bb);
                            if (debug) std::cout << "Accepted mu_tau_b_b event\n";
                            mPT_pT = tmpMPT->MET;
                            mPT_phi = tmpMPT->Phi;
                            if (options["-i"].find("MG5_pp_hh_13TeV_10M_py8_Forced") != std::string::npos) { //13TeV Signal
                                if (correctDecayChannel(options["-i"], cEvent, &mcTruthPlots, &hBB, &hTauTau)) {
                                    gen_mctMatch = truthCut(options["-i"], cEvent, bJet_0, bJet_1, //Checks final-state selection was correct
                                                        taus[0], muons[0], hBB, hTauTau, {"tau", "muon"},
                                                        &mcTruthPlots, &v_gen_higgs_bb, &v_gen_higgs_tt,
                                                        &v_gen_tau_0, &v_gen_tau_1, &v_gen_bJet_0, &v_gen_bJet_1);
                                }
                            }
                            gen_t_0_pT = v_gen_tau_0.Pt();
                            gen_t_0_eta = v_gen_tau_0.Eta();
                            gen_t_0_phi = v_gen_tau_0.Phi();
                            gen_t_0_E = v_gen_tau_0.E();
                            gen_t_1_pT = v_gen_tau_1.Pt();
                            gen_t_1_eta = v_gen_tau_1.Eta();
                            gen_t_1_phi = v_gen_tau_1.Phi();
                            gen_t_1_E = v_gen_tau_1.E();
                            gen_b_0_pT = v_gen_bJet_0.Pt();
                            gen_b_0_eta = v_gen_bJet_0.Eta();
                            gen_b_0_phi = v_gen_bJet_0.Phi();
                            gen_b_0_E = v_gen_bJet_0.E();
                            gen_b_1_pT = v_gen_bJet_1.Pt();
                            gen_b_1_eta = v_gen_bJet_1.Eta();
                            gen_b_1_phi = v_gen_bJet_1.Phi();
                            gen_b_1_E = v_gen_bJet_1.E();
                            gen_diH_pT = v_gen_diHiggs.Pt();
                            gen_diH_eta = v_gen_diHiggs.Eta();
                            gen_diH_phi = v_gen_diHiggs.Phi();
                            gen_diH_E = v_gen_diHiggs.E();
                            gen_h_bb_pT = v_gen_higgs_bb.Pt();
                            gen_h_bb_eta = v_gen_higgs_bb.Eta();
                            gen_h_bb_phi = v_gen_higgs_bb.Phi();
                            gen_h_bb_E = v_gen_higgs_bb.E();
                            gen_h_tt_pT = v_gen_higgs_tt.Pt();
                            gen_h_tt_eta = v_gen_higgs_tt.Eta();
                            gen_h_tt_phi = v_gen_higgs_tt.Phi();
                            gen_h_tt_E = v_gen_higgs_tt.E();
                            t_0_pT = v_tau_0.Pt();
                            t_0_eta = v_tau_0.Eta();
                            t_0_phi = v_tau_0.Phi();
                            t_0_mass = v_tau_0.M();
                            t_0_mT = getMT(t_0_pT, mPT_pT, ROOT::Math::VectorUtil::DeltaPhi(v_tau_0, tmpMPT->P4()));
                            t_1_pT = v_tau_1.Pt();
                            t_1_eta = v_tau_1.Eta();
                            t_1_phi = v_tau_1.Phi();
                            t_1_mass = muMass;
                            t_1_mT = getMT(t_1_pT, mPT_pT, ROOT::Math::VectorUtil::DeltaPhi(v_tau_1, tmpMPT->P4()));
                            b_0_pT = v_bJet_0.Pt();
                            b_0_eta = v_bJet_0.Eta();
                            b_0_phi = v_bJet_0.Phi();
                            b_0_mass = v_bJet_0.M();
                            b_1_pT = v_bJet_1.Pt();
                            b_1_eta = v_bJet_1.Eta();
                            b_1_phi = v_bJet_1.Phi();
                            b_1_mass = v_bJet_1.M();
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
                            diH_mT2 = getMT2(v_tau_0, v_tau_1, v_bJet_0, v_bJet_1, tmpMPT->P4());
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
                            if (debug) std::cout << "1\n";

                            //tmpWeight = (Weight*)branchWeights->At(0);
                            //weight = tmpWeight->Weight;
                            mu_tau_b_b->Fill();
                            h_datasetSizes->Fill("#mu #tau_{h} b #bar{b}", 1);
                            eventAccepted = true;
                            if (debug) std::cout << "Written mu_tau_b_b event\n";
                        }
                    }
                }
            }
        }
        //___________________________________
        if (eventAccepted) continue;
        //Check for e tau b b finalstates___
        h_e_tau_b_b_cutFlow->Fill("All", 1);
        finalstateSet("e_tau_b_b");
        electrons.clear();
        muons.clear();
        taus.clear();
        bJets.clear();
        addMuon = false;
        addElectron = false;
        if (debug) std::cout << "Running e tau b b\n";
        for (int i = 0; i < branchElectron->GetEntries(); i++) { //Loop through Electons
            tmpElectron = (Electron*) branchElectron->At(i);
            if (tmpElectron->PT > ePTMin && std::abs(tmpElectron->Eta) < eEtaMax
                    && tmpElectron->IsolationVar < eIsoMax) { //Quality Electron
                electrons.push_back(i);
            } else if (tmpElectron->PT > ePTMinAdd && std::abs(tmpElectron->Eta) < eEtaMaxAdd
                    && tmpElectron->IsolationVar < eIsoMaxAdd) { //Additional electon
                addElectron = true;
                break;
            }
        }
        if (electrons.size() == 1 && !addElectron) { //One quality electrons found and no additional electrons
            h_e_tau_b_b_cutFlow->Fill("Quality e", 1);
            tmpElectron = (Electron*) branchElectron->At(electrons[0]);
            for (int i = 0; i < branchMuon->GetEntries(); i++) { //Loop through muons
                tmpMuon = (Muon*) branchMuon->At(i);
                if (tmpMuon->PT > muPTMinAdd && std::abs(tmpMuon->Eta) < muEtaMaxAdd
                        && tmpMuon->IsolationVar < muIsoMaxAdd) { //Additional muon
                    addMuon = true;
                    break;
                }
            }
            if (!addMuon) { //No additional muons found
                h_e_tau_b_b_cutFlow->Fill("1 e & 0 #mu", 1);
                for (int i = 0; i < branchJet->GetEntries(); i++) { //Loop through jets
                    tmpJet = (Jet*) branchJet->At(i);
                    if (tauTags[i] == true && tmpJet->BTag == 0 && tmpJet->PT > tauPTMin
                            && std::abs(tmpJet->Eta) < tauEtaMax
                            && tmpJet->Charge != tmpElectron->Charge) { //Quality  OS tau
                        taus.push_back(i);
                    }
                    if (options["-i"].find("13Te") != std::string::npos) {
                        if (tauTags[i] == false && tmpJet->BTag == 1 && tmpJet->PT > bJetPTMin
                                && std::abs(tmpJet->Eta) < bJetEtaMax) { //Quality b jet
                            bJets.push_back(i);
                        }
                    } else {
                        if (tauTags[i] == false && (tmpJet->BTag & (1 << 4) ) && tmpJet->PT > bJetPTMin
                                && std::abs(tmpJet->Eta) < bJetEtaMax) { //Quality b jet
                            bJets.push_back(i);
                        }
                    }
                }
                if (taus.size() >= 1) {//Quality tau
                    h_e_tau_b_b_cutFlow->Fill("Quality #tau", 1);
                    if (bJets.size() >= 2) {//Quality b jets pairs found
                        h_e_tau_b_b_cutFlow->Fill("Quality b#bar{b}", 1);
                        if (selectBJets(branchJet, &bJets, &bJet_0, &bJet_1) == true) { //Quality b-jet pair found
                            if (debug) std::cout << "Accepting event\n";
                            v_tau_1 = tmpElectron->P4();
                            gen_t_1_real = true;
                            tmpJet = (Jet*)branchJet->At(taus[0]);
                            v_tau_0 = tmpJet->P4();
                            gen_t_0_real = (tmpJet->TauWeight > 0.1) ? true : false;
                            tmpMPT = (MissingET*)branchMissingET->At(0);
                            v_higgs_tt = getHiggs2Taus(tmpMPT, v_tau_0, v_tau_1);
                            tmpJet = (Jet*)branchJet->At(bJet_0);
                            v_bJet_0 = tmpJet->P4();
                            gen_b_0_real = bjets_real[bJet_0];
                            tmpJet = (Jet*)branchJet->At(bJet_1);
                            v_bJet_1 = tmpJet->P4();
                            gen_b_1_real = bjets_real[bJet_1];
                            v_higgs_bb = getHiggs2Bs(v_bJet_0, v_bJet_1);
                            v_diHiggs = getDiHiggs(v_higgs_tt, v_higgs_bb);
                            if (debug) std::cout << "Accepted e_tau_b_b event\n";
                            if (options["-i"].find("MG5_pp_hh_13TeV_10M_py8_Forced") != std::string::npos) { //13TeV Signal
                                if (correctDecayChannel(options["-i"], cEvent, &mcTruthPlots, &hBB, &hTauTau)) {
                                    gen_mctMatch = truthCut(options["-i"], cEvent, bJet_0, bJet_1, //Checks final-state selection was correct
                                                        taus[0], electrons[0], hBB, hTauTau, {"tau", "electon"},
                                                        &mcTruthPlots, &v_gen_higgs_bb, &v_gen_higgs_tt,
                                                        &v_gen_tau_0, &v_gen_tau_1, &v_gen_bJet_0, &v_gen_bJet_1);
                                }
                            }
                            gen_t_0_pT = v_gen_tau_0.Pt();
                            gen_t_0_eta = v_gen_tau_0.Eta();
                            gen_t_0_phi = v_gen_tau_0.Phi();
                            gen_t_0_E = v_gen_tau_0.E();
                            gen_t_1_pT = v_gen_tau_1.Pt();
                            gen_t_1_eta = v_gen_tau_1.Eta();
                            gen_t_1_phi = v_gen_tau_1.Phi();
                            gen_t_1_E = v_gen_tau_1.E();
                            gen_b_0_pT = v_gen_bJet_0.Pt();
                            gen_b_0_eta = v_gen_bJet_0.Eta();
                            gen_b_0_phi = v_gen_bJet_0.Phi();
                            gen_b_0_E = v_gen_bJet_0.E();
                            gen_b_1_pT = v_gen_bJet_1.Pt();
                            gen_b_1_eta = v_gen_bJet_1.Eta();
                            gen_b_1_phi = v_gen_bJet_1.Phi();
                            gen_b_1_E = v_gen_bJet_1.E();
                            gen_diH_pT = v_gen_diHiggs.Pt();
                            gen_diH_eta = v_gen_diHiggs.Eta();
                            gen_diH_phi = v_gen_diHiggs.Phi();
                            gen_diH_E = v_gen_diHiggs.E();
                            gen_h_bb_pT = v_gen_higgs_bb.Pt();
                            gen_h_bb_eta = v_gen_higgs_bb.Eta();
                            gen_h_bb_phi = v_gen_higgs_bb.Phi();
                            gen_h_bb_E = v_gen_higgs_bb.E();
                            gen_h_tt_pT = v_gen_higgs_tt.Pt();
                            gen_h_tt_eta = v_gen_higgs_tt.Eta();
                            gen_h_tt_phi = v_gen_higgs_tt.Phi();
                            gen_h_tt_E = v_gen_higgs_tt.E();
                            mPT_pT = tmpMPT->MET;
                            mPT_phi = tmpMPT->Phi;
                            t_0_pT = v_tau_0.Pt();
                            t_0_eta = v_tau_0.Eta();
                            t_0_phi = v_tau_0.Phi();
                            t_0_mass = v_tau_0.M();
                            t_0_mT = getMT(t_0_pT, mPT_pT, ROOT::Math::VectorUtil::DeltaPhi(v_tau_0, tmpMPT->P4()));
                            t_1_pT = v_tau_1.Pt();
                            t_1_eta = v_tau_1.Eta();
                            t_1_phi = v_tau_1.Phi();
                            t_1_mass = eMass;
                            t_1_mT = getMT(t_1_pT, mPT_pT, ROOT::Math::VectorUtil::DeltaPhi(v_tau_1, tmpMPT->P4()));
                            b_0_pT = v_bJet_0.Pt();
                            b_0_eta = v_bJet_0.Eta();
                            b_0_phi = v_bJet_0.Phi();
                            b_0_mass = v_bJet_0.M();
                            b_1_pT = v_bJet_1.Pt();
                            b_1_eta = v_bJet_1.Eta();
                            b_1_phi = v_bJet_1.Phi();
                            b_1_mass = v_bJet_1.M();
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
                            diH_mT2 = getMT2(v_tau_0, v_tau_1, v_bJet_0, v_bJet_1, tmpMPT->P4());
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
                            //tmpWeight = (Weight*)branchWeights->At(0);
                            //weight = tmpWeight->Weight;
                            e_tau_b_b->Fill();
                            h_datasetSizes->Fill("e #tau_{h} b #bar{b}", 1);
                            eventAccepted = true;
                        }
                    }
                }
            }
        }
        //___________________________________
        if (eventAccepted) continue;
        //Check for tau tau b b finalstates___
        h_tau_tau_b_b_cutFlow->Fill("All", 1);
        finalstateSet("tau_tau_b_b");
        electrons.clear();
        muons.clear();
        taus.clear();
        bJets.clear();
        addMuon = false;
        addElectron = false;
        if (debug) std::cout << "Running tau tau b b\n";
        for (int i = 0; i < branchElectron->GetEntries(); i++) { //Loop through Electons
            tmpElectron = (Electron*) branchElectron->At(i);
            if (tmpElectron->PT > ePTMinAdd && std::abs(tmpElectron->Eta) < eEtaMaxAdd
                    && tmpElectron->IsolationVar < eIsoMaxAdd) { //Additional electon
                addElectron = true;
                break;
            }
        }
        if (!addElectron) { //No additional electrons
            for (int i = 0; i < branchMuon->GetEntries(); i++) { //Loop through muons
                tmpMuon = (Muon*) branchMuon->At(i);
                if (tmpMuon->PT > muPTMinAdd && std::abs(tmpMuon->Eta) < muEtaMaxAdd
                        && tmpMuon->IsolationVar < muIsoMaxAdd) { //Additional muon
                    addMuon = true;
                    break;
                }
            }
            if (!addMuon) { //No additional muons found
                h_tau_tau_b_b_cutFlow->Fill("0 e & 0 #mu", 1);
                for (int i = 0; i < branchJet->GetEntries(); i++) { //Loop through jets
                    tmpJet = (Jet*) branchJet->At(i);
                    if (tauTags[i] == true && tmpJet->BTag == 0 && tmpJet->PT > tauPTMin
                            && std::abs(tmpJet->Eta) < tauEtaMax) { //Quality tau
                        taus.push_back(i);
                    }
                    if (options["-i"].find("13Te") != std::string::npos) {
                        if (tauTags[i] == false && tmpJet->BTag == 1 && tmpJet->PT > bJetPTMin
                                && std::abs(tmpJet->Eta) < bJetEtaMax) { //Quality b jet
                            bJets.push_back(i);
                        }
                    } else {
                        if (tauTags[i] == false && (tmpJet->BTag & (1 << 4) ) && tmpJet->PT > bJetPTMin
                                && std::abs(tmpJet->Eta) < bJetEtaMax) { //Quality b jet
                            bJets.push_back(i);
                        }
                    }
                }
                if (taus.size() >= 2) {//2 quality taus
                    h_tau_tau_b_b_cutFlow->Fill("Quality #tau#tau", 1);
                    if (getOSTauTauPair(branchJet, &taus, &tau_0, &tau_1)) { //OS Tau pair
                        h_tau_tau_b_b_cutFlow->Fill("OS", 1);
                        if (bJets.size() >= 2) {//Quality b jets pairs found
                            h_tau_tau_b_b_cutFlow->Fill("Quality b#bar{b}", 1);
                            if (selectBJets(branchJet, &bJets, &bJet_0, &bJet_1) == true) { //Quality b-jet pair found
                                if (debug) std::cout << "Accepting event\n";
                                tmpJet = (Jet*)branchJet->At(tau_0);
                                v_tau_0 = tmpJet->P4();
                                gen_t_0_real = (tmpJet->TauWeight > 0.1) ? true : false;
                                tmpJet = (Jet*)branchJet->At(tau_1);
                                v_tau_1 = tmpJet->P4();
                                gen_t_1_real = (tmpJet->TauWeight > 0.1) ? true : false;
                                tmpMPT = (MissingET*)branchMissingET->At(0);
                                v_higgs_tt = getHiggs2Taus(tmpMPT, v_tau_0, v_tau_1);
                                tmpJet = (Jet*)branchJet->At(bJet_0);
                                v_bJet_0 = tmpJet->P4();
                                gen_b_0_real = bjets_real[bJet_0];
                                tmpJet = (Jet*)branchJet->At(bJet_1);
                                v_bJet_1 = tmpJet->P4();
                                gen_b_1_real = bjets_real[bJet_1];
                                v_higgs_bb = getHiggs2Bs(v_bJet_0, v_bJet_1);
                                v_diHiggs = getDiHiggs(v_higgs_tt, v_higgs_bb);
                                if (debug) std::cout << "Accepted tau_tau_b_b event\n";
                                if (options["-i"].find("MG5_pp_hh_13TeV_10M_py8_Forced") != std::string::npos) { //13TeV Signal
                                    if (correctDecayChannel(options["-i"], cEvent, &mcTruthPlots, &hBB, &hTauTau)) {
                                    gen_mctMatch = truthCut(options["-i"], cEvent, bJet_0, bJet_1, //Checks final-state selection was correct
                                                        tau_0, tau_1, hBB, hTauTau, {"tau", "tau"},
                                                        &mcTruthPlots, &v_gen_higgs_bb, &v_gen_higgs_tt,
                                                        &v_gen_tau_0, &v_gen_tau_1, &v_gen_bJet_0, &v_gen_bJet_1);
                                    }
                                }
                                gen_t_0_pT = v_gen_tau_0.Pt();
                                gen_t_0_eta = v_gen_tau_0.Eta();
                                gen_t_0_phi = v_gen_tau_0.Phi();
                                gen_t_0_E = v_gen_tau_0.E();
                                gen_t_1_pT = v_gen_tau_1.Pt();
                                gen_t_1_eta = v_gen_tau_1.Eta();
                                gen_t_1_phi = v_gen_tau_1.Phi();
                                gen_t_1_E = v_gen_tau_1.E();
                                gen_b_0_pT = v_gen_bJet_0.Pt();
                                gen_b_0_eta = v_gen_bJet_0.Eta();
                                gen_b_0_phi = v_gen_bJet_0.Phi();
                                gen_b_0_E = v_gen_bJet_0.E();
                                gen_b_1_pT = v_gen_bJet_1.Pt();
                                gen_b_1_eta = v_gen_bJet_1.Eta();
                                gen_b_1_phi = v_gen_bJet_1.Phi();
                                gen_b_1_E = v_gen_bJet_1.E();
                                gen_diH_pT = v_gen_diHiggs.Pt();
                                gen_diH_eta = v_gen_diHiggs.Eta();
                                gen_diH_phi = v_gen_diHiggs.Phi();
                                gen_diH_E = v_gen_diHiggs.E();
                                gen_h_bb_pT = v_gen_higgs_bb.Pt();
                                gen_h_bb_eta = v_gen_higgs_bb.Eta();
                                gen_h_bb_phi = v_gen_higgs_bb.Phi();
                                gen_h_bb_E = v_gen_higgs_bb.E();
                                gen_h_tt_pT = v_gen_higgs_tt.Pt();
                                gen_h_tt_eta = v_gen_higgs_tt.Eta();
                                gen_h_tt_phi = v_gen_higgs_tt.Phi();
                                gen_h_tt_E = v_gen_higgs_tt.E();
                                mPT_pT = tmpMPT->MET;
                                mPT_phi = tmpMPT->Phi;
                                t_0_pT = v_tau_0.Pt();
                                t_0_eta = v_tau_0.Eta();
                                t_0_phi = v_tau_0.Phi();
                                t_0_mass = v_tau_0.M();
                                t_0_mT = getMT(t_0_pT, mPT_pT, ROOT::Math::VectorUtil::DeltaPhi(v_tau_0, tmpMPT->P4()));
                                t_1_pT = v_tau_1.Pt();
                                t_1_eta = v_tau_1.Eta();
                                t_1_phi = v_tau_1.Phi();
                                t_1_mass = v_tau_0.M();;
                                t_1_mT = getMT(t_1_pT, mPT_pT, ROOT::Math::VectorUtil::DeltaPhi(v_tau_1, tmpMPT->P4()));
                                b_0_pT = v_bJet_0.Pt();
                                b_0_eta = v_bJet_0.Eta();
                                b_0_phi = v_bJet_0.Phi();
                                b_0_mass = v_bJet_0.M();
                                b_1_pT = v_bJet_1.Pt();
                                b_1_eta = v_bJet_1.Eta();
                                b_1_phi = v_bJet_1.Phi();
                                b_1_mass = v_bJet_1.M();
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
                                diH_mT2 = getMT2(v_tau_0, v_tau_1, v_bJet_0, v_bJet_1, tmpMPT->P4());
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
                                //tmpWeight = (Weight*)branchWeights->At(0);
                                //weight = tmpWeight->Weight;
                                tau_tau_b_b->Fill();
                                h_datasetSizes->Fill("#tau_{h} #tau_{h} b #bar{b}", 1);
                                eventAccepted = true;
                            }
                        }
                    }
                }
            }
        }
        //___________________________________
    }
    std::cout << "Event loop complete\n";
    //___________________________________________
    //Writing plots______________________________
    std::cout << "Creating plots\n";
    TCanvas* c_datasetSizes = new TCanvas();
    c_datasetSizes->SetLogy();
    h_datasetSizes->GetXaxis()->SetTitle("Dataset");
    h_datasetSizes->GetYaxis()->SetTitle("Events");
    h_datasetSizes->Draw();
    h_datasetSizes->Write();
    delete c_datasetSizes;
    TCanvas* c_e_tau_b_b_cutFlow = new TCanvas();
    c_e_tau_b_b_cutFlow->SetLogy();
    h_e_tau_b_b_cutFlow->GetXaxis()->SetTitle("Cuts");
    h_e_tau_b_b_cutFlow->GetYaxis()->SetTitle("Events");
    h_e_tau_b_b_cutFlow->Draw();
    h_e_tau_b_b_cutFlow->Write();
    delete c_e_tau_b_b_cutFlow;
    TCanvas* c_mu_tau_b_b_cutFlow = new TCanvas();
    c_mu_tau_b_b_cutFlow->SetLogy();
    h_mu_tau_b_b_cutFlow->GetXaxis()->SetTitle("Cuts");
    h_mu_tau_b_b_cutFlow->GetYaxis()->SetTitle("Events");
    h_mu_tau_b_b_cutFlow->Draw();
    h_mu_tau_b_b_cutFlow->Write();
    delete c_mu_tau_b_b_cutFlow;
    TCanvas* c_tau_tau_b_b_cutFlow = new TCanvas();
    c_tau_tau_b_b_cutFlow->SetLogy();
    h_tau_tau_b_b_cutFlow->GetXaxis()->SetTitle("Cuts");
    h_tau_tau_b_b_cutFlow->GetYaxis()->SetTitle("Events");
    h_tau_tau_b_b_cutFlow->Draw();
    h_tau_tau_b_b_cutFlow->Write();
    delete c_tau_tau_b_b_cutFlow;
    TCanvas* c_mcTruth_cutFlow = new TCanvas();
    mcTruthPlots["cuts"]->GetXaxis()->SetTitle("Cuts");
    mcTruthPlots["cuts"]->GetYaxis()->SetTitle("Events");
    mcTruthPlots["cuts"]->Draw();
    mcTruthPlots["cuts"]->Write();
    delete c_mcTruth_cutFlow;
    TCanvas* c_mcTruth_bJetMatch = new TCanvas();
    mcTruthPlots["bMatch"]->GetXaxis()->SetTitle("#DeltaR(b, jet)");
    mcTruthPlots["bMatch"]->GetYaxis()->SetTitle("Events");
    mcTruthPlots["bMatch"]->Draw();
    mcTruthPlots["bMatch"]->Write();
    delete c_mcTruth_bJetMatch;
    TCanvas* c_mcTruth_tauJetMatch = new TCanvas();
    mcTruthPlots["tauMatch"]->GetXaxis()->SetTitle("#DeltaR(#tau, jet)");
    mcTruthPlots["tauMatch"]->GetYaxis()->SetTitle("Events");
    mcTruthPlots["tauMatch"]->Draw();
    mcTruthPlots["tauMatch"]->Write();
    delete c_mcTruth_tauJetMatch;
    TCanvas* c_mcTruth_higgsDecay = new TCanvas();
    mcTruthPlots["higgsDecay"]->GetXaxis()->SetTitle("Higgs product |PID|");
    mcTruthPlots["higgsDecay"]->GetYaxis()->SetTitle("Events");
    mcTruthPlots["higgsDecay"]->Draw();
    mcTruthPlots["higgsDecay"]->Write();
    delete c_mcTruth_higgsDecay;
    TCanvas* c_sum_w = new TCanvas();
    h_sum_w->GetXaxis()->SetTitle("#kapa_#lambda");
    h_sum_w->GetYaxis()->SetTitle("#Sigma(w/sm_evts)");
    h_sum_w->Draw();
    h_sum_w->Write();
    delete c_sum_w;
    std::cout << "Plots created\n";
    //___________________________________________
    //Save datasets______________________________
    std::cout << "Saving data\n";
    e_tau_b_b->Write();
    mu_tau_b_b->Write();
    tau_tau_b_b->Write();
    std::cout << "Data saved\n";
    outputFile->Close();
    f_14TeV_coeffs->Close();
    f_HH_14TeV_histo->Close();
    delete outputFile;
    delete e_tau_b_b;
    delete mu_tau_b_b;
    delete tau_tau_b_b;
    delete h_datasetSizes;
    delete h_e_tau_b_b_cutFlow;
    delete h_mu_tau_b_b_cutFlow;
    delete h_tau_tau_b_b_cutFlow;
    delete HH_14TeV_histo;
    //___________________________________________
}