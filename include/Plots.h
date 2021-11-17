#pragma once
#include "PassedEvent.h"
#include <vector>
#include <TH1.h>
#include <THStack.h>
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>

class Plots {
    public:
        static void Plot1(vector<PassedEvent> events);
        static void Plot2(vector<PassedEvent> events);
};