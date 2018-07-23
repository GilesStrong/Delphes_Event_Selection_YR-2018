/*
 * finalstateCuts.hh
 *
 *  Created on: 8 Apr 2016
 *      Author: giles
 */

#ifndef FINALSTATECUTS_HH_
#define FINALSTATECUTS_HH_

double ePTMin = 24.0; //Gev
double eEtaMax = 2.1;
double eIsoMax = 0.1;
double muPTMin = 19.0; //Gev
double muEtaMax = 2.1;
double muIsoMax = 0.1;
double ePTMinAdd = 10.0; //Gev
double eEtaMaxAdd = 2.5;
double eIsoMaxAdd = 0.3;
double muPTMinAdd = 10.0; //Gev
double muEtaMaxAdd = 2.4;
double muIsoMaxAdd = 0.3;
double tauPTMin; //Gev
double tauEtaMax;
double bJetPTMin = 30.0; //Gev
double bJetEtaMax = 2.4;
double higgsMassMin = 80; //GeV
double higgsMassMax = 160; //GeV

bool massCut = false;

void finalstateSet(std::string set) {
	if (set == "e_tau_b_b") {
		tauPTMin = 20.0;
		tauEtaMax = 2.3;
	} else if (set == "mu_tau_b_b") {
		tauPTMin = 20.0;
		tauEtaMax = 2.3;
	} else if (set == "tau_tau_b_b") {
		tauPTMin = 45.0;
		tauEtaMax = 2.1;
	} else if (set == "e_e_b_b") {
		tauPTMin = 10.0;
		tauEtaMax = 2.5;
	} else if (set == "e_mu_b_b") {
		tauPTMin = 10.0;
		tauEtaMax = 2.5;
	} else if (set == "mu_mu_b_b") {
		tauPTMin = 10.0;
		tauEtaMax = 2.5;
	}
}
#endif /* FINALSTATECUTS_HH_ */
