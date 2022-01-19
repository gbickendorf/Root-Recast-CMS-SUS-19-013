#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>
#include <random>

#include "PassedEvent.h"
#include "RootIO.h"

#include "ExRootAnalysis/ExRootClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"

#include <TClonesArray.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include "classes/DelphesClasses.h"
#include <boost/progress.hpp>
#include <algorithm>
#include <boost/math/special_functions/chebyshev.hpp>



class Analysis {
    private:
        static double sq(double x);
        static double excludeSignalRegion(Double_t *x, Double_t *par);
        static TF1 *fitFunc;
    public:
        static double AnalyseEvents(ExRootTreeReader *treeReader, vector<PassedEvent> &events , int MinNumEvents, int status);
        static double AnalyseEventsNew(ExRootTreeReader *treeReader, vector<PassedEvent> &events , int MinNumEvents, int status);
        static void AnalyseEventsSingleLeptonSample(ExRootTreeReader *treeReader, vector<PassedEvent> &events , int MinNumEvents);
        static void AnalyseEventsSinglePhotonSample(ExRootTreeReader *treeReader, vector<PassedEvent> &events , int MinNumEvents);
        static void AnalyseEvents2(ExRootTreeReader *treeReader, vector<PassedEvent> &events , int MinNumEvents, int status);
        static double FindNormalisation(const char * filename);
        static double FindNormalisation(const char * filename, int nFiles);
        static double FindNormalisation(const char * filename, int nFiles, vector<double> &mets);
        static void SampleFromEvents(vector<PassedEvent> &events);
        static void FitCherbyshev(vector<PassedEvent> events ,double &bNorm, double &bNormErr,vector<double> &params);
        static void CalcTransferFactor(vector<PassedEvent> events,double bNorm,double bNormErr, vector<double> &NCRi, double &transferFactor, double &transferFactorErr);
        static void GetEventsInSignalRegion(vector<PassedEvent> events,vector<PassedEvent> &signalevents);
        static void GetSRBinContent(vector<PassedEvent> &signalevents,vector<double> &binContent);
        static double TheoryXSection(double m_Gluino);
};
