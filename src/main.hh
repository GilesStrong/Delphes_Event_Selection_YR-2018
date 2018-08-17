/*
 * main.hh
 *
 *  Created on: 8 Apr 2016
 *      Author: giles giles.strong@outlook.com
 */

#ifndef MAIN_HH_
#define MAIN_HH_
//C++
#include <string>
#include <iostream>
#include <vector>
#include <utility>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <random>
#include <cmath>
//Local
#include "myMethods.hh"
#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"
#include "finalstateCuts.hh"
//Root
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TLatex.h"
#include "RtypesCore.h"
#include "TString.h"
#include "TH1D.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TClonesArray.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooDecay.h"
#include "RooCBShape.h"
#include "RooFormula.h"
#include "RooGaussModel.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"
#include "RooAddModel.h"
#include "RooAddPdf.h"
#include "RooUniform.h"
#include "RooMoment.h"
#include "Math/GenVector/VectorUtil.h"
//Boost
#include "boost/regex.hpp"
#include "boost/filesystem.hpp"
#include "boost/algorithm/string.hpp"

const double eMass = 0.0005109989; //GeV
const double muMass = 0.1056583715; //GeV
#endif /* MAIN_HH_ */
