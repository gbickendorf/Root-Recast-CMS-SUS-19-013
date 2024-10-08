#pragma once
#include "PassedEvent.h"
#include <vector>
#include <TH1.h>
#include <THStack.h>
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH2D.h>
#include <TF1.h>
#include <TRatioPlot.h>
#include <TGraphErrors.h>

#include<RootIO.h>
#include<Analysis.h>

class Plots {
    public:
        static void PlotAll();
        static void PlotPTMiss(vector<PassedEvent> events);
        static void PlotMJ1(vector<PassedEvent> events,vector<double> params);
        static void PlotPTShape(vector<PassedEvent> events);
        static void PlotSignalRegion(vector<PassedEvent> events, vector<double> NCRi, double transferFactor, double transferFactorErr);
        static void PlotPhotonLeptonValidation(vector<PassedEvent> events);
};