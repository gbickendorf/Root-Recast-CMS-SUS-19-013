#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>
#include <random>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include "TROOT.h"
#include "TCanvas.h"

#include "PassedEvent.h"
#include "RootIO.h"

using namespace std;

class Stats
{
public:
    struct LogLikelihoodParams
    {
        vector<double> data;
        vector<double> signal;
        vector<double> CRi;
        double transferFac;

        vector<double> CRi_expected;
        vector<double> CRiErr_expected;
        double transferFac_expected;
        double transferFacErr_expected;
        double signalStrength;
        int mode;
    };
    static double LnLogNormal(double x, double xhat, double kappa);
    static double LnPoisson(double n, double nexpected);
    static double LogLikelihood(LogLikelihoodParams *par);
    static double LogLikelihoodMax(LogLikelihoodParams &par, double mu, int verbose);
    static double LogLikelihoodMax(LogLikelihoodParams &par, double mu_min, double mu_max,int verbose);
    static void ReadParameters(LogLikelihoodParams &par);
    static void Test();
    static void TestMinimise();

private:
    static double minFunc(const gsl_vector *v, void *params);
};
