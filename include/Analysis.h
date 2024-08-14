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
#include <TFitResultPtr.h>
#include <TPaveText.h>
#include <TStyle.h>
#include "classes/DelphesClasses.h"
#include <boost/progress.hpp>
#include <algorithm>
#include <boost/math/special_functions/chebyshev.hpp>



class Analysis {
    private:
        static double sq(double x);
        static double excludeSignalRegion(Double_t *x, Double_t *par);
        static TF1 *fitFunc;
        static double SRmjmin;
        static double SRmjmax;
    public:
        static double AnalyseEvents(ExRootTreeReader *treeReader, vector<PassedEvent> &events , int MinNumEvents, int status);
        static double AnalyseEventsNew(ExRootTreeReader *treeReader, vector<PassedEvent> &events , int MinNumEvents, int status);
        static void AnalyseEventsSingleLeptonSample(ExRootTreeReader *treeReader, vector<PassedEvent> &events , int MinNumEvents);
        static void AnalyseEventsSinglePhotonSample(ExRootTreeReader *treeReader, vector<PassedEvent> &events , int MinNumEvents);
        static void AnalyseEvents2(ExRootTreeReader *treeReader, vector<PassedEvent> &events , int MinNumEvents, int status);
        static void SampleFromEvents(vector<PassedEvent> &events);
        static void SampleFromEvents(vector<PassedEvent> &events, double intLumi);
        static void SampleFromEventsNewPassedEventDef(vector<PassedEvent> &events);
        static void SampleFromEventsNewPassedEventDef(vector<PassedEvent> &events, double intLumi);
        static void SampleFromEventsNewPassedEventDefHadronicBaseline(vector<PassedEvent> &events, double intLumi);
        static void FitCherbyshev(vector<PassedEvent> events ,double &bNorm, double &bNormErr,vector<double> &params, double mjmin, double mjmax);
        static void CalcTransferFactor(vector<PassedEvent> events,double bNorm,double bNormErr, vector<double> &NCRi, double &transferFactor, double &transferFactorErr, double mjmin, double mjmax);
        static void GetEventsInSignalRegion(vector<PassedEvent> events,vector<PassedEvent> &signalevents, double mjmin, double mjmax);
        static void GetSRBinContent(vector<PassedEvent> &signalevents,vector<double> &binContent);
        static double TheoryXSection(double m_Gluino);
        static double TheoryXSectionUpper(double m_Gluino);
        static double TheoryXSectionLower(double m_Gluino);
};
