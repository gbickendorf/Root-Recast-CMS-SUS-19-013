
#include <iostream>
#include "TROOT.h"
#include <unistd.h>
#include <TLegend.h>
#include <vector>
#include <TFile.h>

#include <TH1F.h>

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
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit;

void rf706_histpdf()
{
    // C r e a t e   p d f   f o r   s a m p l i n g
    // ---------------------------------------------

    RooRealVar x("x", "x", 0, 20);
    RooPolynomial p("p", "p", x, RooArgList(RooConst(0.01), RooConst(-0.01), RooConst(0.0004)));

    // C r e a t e   l o w   s t a t s   h i s t o g r a m
    // ---------------------------------------------------

    // Sample 500 events from p
    x.setBins(20);
    RooDataSet *data1 = p.generate(x, 500);

    // Create a binned dataset with 20 bins and 500 events
    RooDataHist *hist1 = data1->binnedClone();

    // Represent data in dh as pdf in x
    RooHistPdf histpdf1("histpdf1", "histpdf1", x, *hist1, 0);

    // Plot unbinned data and histogram pdf overlaid
    RooPlot *frame1 = x.frame(Title("Low statistics histogram pdf"), Bins(100));
    data1->plotOn(frame1);
    histpdf1.plotOn(frame1);

    // C r e a t e   h i g h   s t a t s   h i s t o g r a m
    // -----------------------------------------------------

    // Sample 100000 events from p
    x.setBins(10);
    RooDataSet *data2 = p.generate(x, 100000);

    // Create a binned dataset with 10 bins and 100K events
    RooDataHist *hist2 = data2->binnedClone();

    // Represent data in dh as pdf in x, apply 2nd order interpolation
    RooHistPdf histpdf2("histpdf2", "histpdf2", x, *hist2, 2);

    // Plot unbinned data and histogram pdf overlaid
    RooPlot *frame2 = x.frame(Title("High stats histogram pdf with interpolation"), Bins(100));
    data2->plotOn(frame2);
    histpdf2.plotOn(frame2);

    TCanvas *c = new TCanvas("rf706_histpdf", "rf706_histpdf", 800, 400);
    c->Divide(2);
    c->cd(1);
    gPad->SetLeftMargin(0.15);
    frame1->GetYaxis()->SetTitleOffset(1.4);
    frame1->Draw();
    c->cd(2);
    gPad->SetLeftMargin(0.15);
    frame2->GetYaxis()->SetTitleOffset(1.8);
    frame2->Draw();
}

void HistPDFTest()
{
    TFile f2("FinalSamples.root", "read");
    TVectorD *vSignalBins = (TVectorD *)f2.Get("SigBins");
    TVectorD *vSigErr = (TVectorD *)f2.Get("SigErr");
    TVectorD *vNCR_i = (TVectorD *)f2.Get("NCR_i");
    TVectorD *vDataBins = (TVectorD *)f2.Get("DataBins");
    TVectorD *vTransferFac = (TVectorD *)f2.Get("TransferFac");
    vector<double> data, signal, CRi, CRiErr;
    for (size_t i = 0; i < 6; i++)
    {
        data.push_back((*vDataBins)[i]);
        signal.push_back((*vSignalBins)[i]);
        CRi.push_back((*vNCR_i)[i]);
        CRiErr.push_back(sqrt((*vNCR_i)[i]));
    }
    double transferFac, transferFacErr;
    transferFac = (*vTransferFac)[0];
    transferFacErr = (*vTransferFac)[1];
    f2.Close();

    TH1F *hist = new TH1F("Hist", "hist", 6, 0, 6);
    for (size_t i = 0; i < 6; i++)
    {
        hist->SetBinContent(i + 1, CRi[i]);
    }
    //gPad->SetLogy();
    auto mycanvas = new TCanvas();
    gPad->SetLogy();
    mycanvas->cd(0);
    //hist->Draw("");
    RooRealVar ptmiss("ptmiss", "ptmiss", 0, 6);
    ptmiss.setBins(6);
    RooDataHist *dh = new RooDataHist("DataHist", "DataHist", RooArgSet(ptmiss), hist);
    RooHistPdf histpdf1("histpdf1", "histpdf1", ptmiss, *dh, 0);

    RooPlot *frame1 = ptmiss.frame(Title("Low statistics histogram pdf"), Bins(6));
    RooPlot *frame2 = ptmiss.frame(Title("High stats histogram pdf with interpolation"), Bins(6));

    RooDataHist *gendat = histpdf1.generateBinned(RooArgSet(ptmiss), 100000);
    gendat->plotOn(frame2);
    dh->plotOn(frame1);
    histpdf1.plotOn(frame1);
    mycanvas->Divide(2);
    mycanvas->cd(1);
    gPad->SetLeftMargin(0.15);
    frame1->GetYaxis()->SetTitleOffset(1.4);
    mycanvas->GetPad(1)->SetLogy();
    mycanvas->GetPad(2)->SetLogy();
    frame1->Draw();
    mycanvas->cd(2);
    gPad->SetLeftMargin(0.15);
    frame2->GetYaxis()->SetTitleOffset(1.8);
    frame2->Draw();
}

void Stats2Small()
{
    vector<double> SigBins, bgBins, dataBins;
    SigBins.push_back(1.0);
    SigBins.push_back(2.0);
    bgBins.push_back(5.0);
    bgBins.push_back(8.0);
    dataBins.push_back(4.0);
    dataBins.push_back(8.0);
    RooWorkspace w("w");
    /*
    RooRealVar mu("mu","mu",1,0,100);
    RooRealVar NSigExp0("NSigExp0","NSigExp0",SigBins[0]);
    RooRealVar NBGExp0("NBGExp0","NBGExp0",bgBins[0]);
    RooRealVar NData0("NData0","NData0",dataBins[0]);
    w.import(mu);
    w.import(NSigExp0);
    w.import(NBGExp0);
    w.import(NData0);
    w.factory("prod:NExpSig(NSigExp0,mu)");
    w.factory("sum:Nexp0(NBGExp0,NExpSig)");
    w.factory("Poisson:model(NData0,Nexp0)");
    */
    w.factory("mu[1,0,10]");
    w.factory("NSigExp0[0,0,100]");
    w.factory("NBGExp0[0,0,100]");
    w.factory("NData0[0,0,100]");
    w.factory("prod:NExpSig(NSigExp0,mu)");
    w.factory("sum:Nexp0(NBGExp0,NExpSig)");
    w.factory("Poisson:model(NData0,Nexp0)");
    w.Print();

    RooStats::ModelConfig mc("ModelConfigs", &w);
    mc.SetPdf(*w.pdf("model"));
    mc.SetParametersOfInterest(*w.var("mu"));
    mc.SetObservables("NData0");
    //mc.SetNuisanceParameters("");
    mc.SetSnapshot(*w.var("mu"));
    mc.SetGlobalObservables("NSigExp0,NBGExp0");
    mc.Print();
    RooDataSet d("data", "data", RooArgSet(*w.var("NData0")));
    d.add(RooArgSet(*w.var("NData0")));
    w.import(d);
    w.writeToFile("model2.root", true);
    w.Print();
    RooAbsData *dataw = w.data("data");
    RooStats::ModelConfig *sbModel = (RooStats::ModelConfig *)w.obj("ModelConfigs");
    cout << sbModel << endl;
    return;
    RooStats::ModelConfig *bModel = (RooStats::ModelConfig *)sbModel->Clone("BonlyModel");
    return;
    RooRealVar *poi = (RooRealVar *)bModel->GetParametersOfInterest()->first();
    poi->setVal(0);
    bModel->SetSnapshot(*poi);

    RooStats::AsymptoticCalculator asympCalc(*dataw, *bModel, *sbModel);

    asympCalc.SetOneSided(true);
    RooStats::HypoTestInverter inverter(asympCalc);
    inverter.SetConfidenceLevel(0.90);
    inverter.UseCLs(true);
    inverter.SetVerbose(false);
    inverter.SetFixedScan(5, 0.0, 6.0);
    RooStats::HypoTestInverterResult *result = inverter.GetInterval();
    cout << 100 * inverter.ConfidenceLevel() << "%  upper limit : " << result->UpperLimit() << endl;
}

void Stats2()
{
    Stats2Small();
    gApplication->Terminate();
    return;
    TFile f("../FinalSamples.root", "read");
    f.Print();
    TVectorD *vSignalBins = (TVectorD *)f.Get("SigBins");
    TVectorD *vSigErr = (TVectorD *)f.Get("SigErr");
    TVectorD *vNCR_i = (TVectorD *)f.Get("NCR_i");
    TVectorD *vDataBins = (TVectorD *)f.Get("DataBins");
    TVectorD *vTransferFac = (TVectorD *)f.Get("TransferFac");
    vector<double> data, signal, CRi, CRiErr;
    for (size_t i = 0; i < 6; i++)
    {
        data.push_back((*vDataBins)[i]);
        signal.push_back((*vSignalBins)[i]);
        CRi.push_back((*vNCR_i)[i]);
        CRiErr.push_back(sqrt((*vNCR_i)[i]));
    }
    double transferFac, transferFacErr;
    transferFac = (*vTransferFac)[0];
    transferFacErr = (*vTransferFac)[1];

    RooWorkspace w("w");
    w.factory("mu[1,0,10]");
    w.factory(TString("transFac[") + transferFac + ",0,10]");
    w.factory(TString("transFac_0[") + transferFac + "]");
    w.var("transFac_0")->setConstant(true);
    w.factory(TString("transFacErr[") + transferFacErr + ",1e-10,1]");
    w.var("transFacErr")->setConstant(true);
    vector<TString> PDFNames, ConstraintNames;
    for (size_t i = 0; i < 6; i++)
    {
        w.factory((TString("NCR") + i + "[" + CRi[i] + ",0,1e6]").Data());
        w.factory((TString("NCR") + i + "_0[" + CRi[i] + "]").Data());
        w.var((TString("NCR") + i + "_0").Data())->setConstant(true);
        w.factory((TString("NSignal") + i + "[" + signal[i] + ",0,1e6]").Data());
        w.factory((TString("NData") + i + "[" + data[i] + ",0,1e6]").Data());
        w.factory(TString("sum:nexp") + i + "(NCR" + i + "*transFac,mu*NSignal" + i + ")");
        PDFNames.push_back(TString("poisson_") + i);
        w.factory((TString("Poisson:") + PDFNames[i] + "(NData" + i + ",nexp" + i + ")").Data());
        //cout << (TString("Gaussian:constraintNC_")+i+"(NCR"+i+"_0,NCR"+i+",sqrt(NCR"+i+"))").Data() << endl;
        ConstraintNames.push_back(TString("constraintNC_") + i);
        w.factory((TString("Poisson:") + ConstraintNames[i] + "(NCR" + i + "_0,NCR" + i + ")").Data());
    }
    w.factory("Gaussian:constraintTransFac(transFac_0,transFac,transFacErr)");
    TString pdfs = "PROD:model(";
    for (size_t i = 0; i < PDFNames.size(); i++)
    {
        pdfs += PDFNames[i] + ",";
    }
    for (size_t i = 0; i < ConstraintNames.size(); i++)
    {
        pdfs += ConstraintNames[i] + ",";
    }
    pdfs += "constraintTransFac)";

    w.factory(pdfs.Data());

    RooStats::ModelConfig mc("ModelConfig", &w);
    mc.SetPdf(*w.pdf("model"));
    mc.SetParametersOfInterest(*w.var("mu"));
    mc.SetObservables("NData0,NData1,NData2,NData3,NData4,NData5");
    mc.SetNuisanceParameters("NCR0,NCR1,NCR2,NCR3,NCR4,NCR5,transFac");
    mc.SetSnapshot(*w.var("mu"));
    mc.SetGlobalObservables("NCR0_0,NCR1_0,NCR2_0,NCR3_0,NCR4_0,NCR5_0,transFac_0");
    mc.Print();

    w.import(mc);

    RooDataSet d("data", "data", RooArgSet(*w.var("NData0"), *w.var("NData1"), *w.var("NData2"), *w.var("NData3"), *w.var("NData4"), *w.var("NData5")));
    d.add(RooArgSet(*w.var("NData0"), *w.var("NData1"), *w.var("NData2"), *w.var("NData3"), *w.var("NData4"), *w.var("NData5")));

    w.import(d);
    /*
    w.factory("sum:nexp(s[3,0,15],b[1,0,10])");
    // Poisson of (n | s+b)
    w.factory("Poisson:pdf(nobs[0,50],nexp)");
    w.factory("Gaussian:constraint(b0[0,10],b,sigmab[1])");
    w.factory("PROD:model(pdf,constraint,constraint)");

    w.var("b0")->setVal(b);
    w.var("b0")->setConstant(true); // needed for being treated as global observables
    w.var("sigmab")->setVal(sigmab * b);

*/
    w.Print();
    w.writeToFile("model.root", true);
    RooAbsData *dataw = w.data("data");

    RooStats::ModelConfig *sbModel = (RooStats::ModelConfig *)w.obj("ModelConfig");
    RooStats::ModelConfig *bModel = (RooStats::ModelConfig *)sbModel->Clone("BonlyModel");
    RooRealVar *poi = (RooRealVar *)bModel->GetParametersOfInterest()->first();
    poi->setVal(0);
    bModel->SetSnapshot(*poi);

    RooStats::AsymptoticCalculator asympCalc(*dataw, *bModel, *sbModel);

    asympCalc.SetOneSided(true);
    RooStats::HypoTestInverter inverter(asympCalc);
    inverter.SetConfidenceLevel(0.90);
    inverter.UseCLs(true);
    inverter.SetVerbose(false);
    inverter.SetFixedScan(5, 0.0, 6.0);
    RooStats::HypoTestInverterResult *result = inverter.GetInterval();
    cout << 100 * inverter.ConfidenceLevel() << "%  upper limit : " << result->UpperLimit() << endl;

    cout << "DONE" << endl;
    gApplication->Terminate();
}