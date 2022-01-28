
#include <iostream>
#include "TROOT.h"
#include <unistd.h>
#include <TLegend.h>
#include <vector>

#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"
#include "RooRandom.h"

#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"
#include "RooStats/ModelConfig.h"

#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"
#include "RooStats/NumEventsTestStat.h"

#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"

void Stats(int nobs = 3,        // number of observed events
           double b = 1,        // number of background events
           double sigmab = 0.2) // relative uncertainty in b
{
    RooWorkspace w("w");

    // make Poisson model * Gaussian constraint
    //RooRealVar sr("s","signal",3,0,15);
    w.factory("s[3,0,15]");
    w.factory("sum:nexp(s,b[1,0,10])");
    // Poisson of (n | s+b)
    w.factory("Poisson:pdf(nobs[0,50],nexp)");
    w.factory("Gaussian:constraint(b0[0,10],b,sigmab[1])");
    w.factory("PROD:model(pdf,constraint,constraint)");

    w.var("b0")->setVal(b);
    w.var("b0")->setConstant(true); // needed for being treated as global observables
    w.var("sigmab")->setVal(sigmab * b);

    RooStats::ModelConfig mc("ModelConfig", &w);
    mc.SetPdf(*w.pdf("model"));
    mc.SetParametersOfInterest(*w.var("s"));
    mc.SetObservables(*w.var("nobs"));
    mc.SetNuisanceParameters(*w.var("b"));

    // these are needed for the hypothesis tests
    mc.SetSnapshot(*w.var("s"));
    mc.SetGlobalObservables("b0");

    mc.Print();
    // import model in the workspace
    w.import(mc);

    // make data set with the namber of observed events
    RooDataSet data("data", "", *w.var("nobs"));
    w.var("nobs")->setVal(3);
    data.add(*w.var("nobs"));
    // import data set in workspace and save it in a file
    w.import(data);

    //w.Print();

    TString fileName = "CountingModel.root";
    w.Print();
    // write workspace in the file (recreate file if already existing)
    w.writeToFile(fileName, true);
    gApplication->Terminate();
}