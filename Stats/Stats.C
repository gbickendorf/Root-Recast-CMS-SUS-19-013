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
//#include "RooStats/HistFactory/"

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

#include <iostream>
#include "TROOT.h"
#include <unistd.h>
#include <TLegend.h>
#include <vector>
#include <TFile.h>
#include <TVectorD.h>

#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "RooStats/HistFactory/Systematics.h"

#include "StandardHypoTestInvDemo.C"

#include <vector>
using namespace RooFit;
using namespace RooStats;
using namespace std;
using namespace RooStats::HistFactory;

void HistFactoryModel(double m_Go)
{
   ostringstream signalOutFile;
   signalOutFile << "/home/gerrit/Documents/Projects/Root-Recast-CMS-SUS-19-013/RootFiles/Signal" << m_Go << "GeV.root";

   RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
   TFile f1(signalOutFile.str().c_str(), "read");
   TFile f2("/home/gerrit/Documents/Projects/Root-Recast-CMS-SUS-19-013/RootFiles/FinalSamples.root", "read");
   TVectorD *vSignalBins = (TVectorD *)f1.Get("SigBins");
   TVectorD *vSigErr = (TVectorD *)f1.Get("SigErr");
   TVectorD *vM_Go = (TVectorD *)f1.Get("M_Go");
   TVectorD *vCrossSecErr = (TVectorD *)f1.Get("CrossSecErr");
   TVectorD *vNCRi = (TVectorD *)f2.Get("NCR_i");
   TVectorD *vDataBins = (TVectorD *)f2.Get("DataBins");
   TVectorD *vTransferFac = (TVectorD *)f2.Get("TransferFac");

   cout << "+++++++++++++++++++++" << (*vM_Go)[0] << "GeV+++++++++++++++++++++++++" << endl;
   vector<double> vdata, vsignal, vsignalErr, vCRi, vCRiErr, vShapeSyst;

   for (size_t i = 0; i < 6; i++)
   {
      vdata.push_back((*vDataBins)[i]);
      vsignal.push_back((*vSignalBins)[i]);
      vCRi.push_back((*vNCRi)[i]);
      vCRiErr.push_back(sqrt((*vNCRi)[i]));
      vsignalErr.push_back((*vSigErr)[i]);
   }

   double transferFac, transferFacErr;
   transferFac = (*vTransferFac)[0];
   transferFacErr = (*vTransferFac)[1];
   f1.Close();
   f2.Close();
   /*
      vdata = {237, 67, 20, 3, 3, 1};
      vsignal = {3.5, 4.3, 6.6, 7.2, 7.2, 11.6};
      vCRi = {1191, 320, 112, 16, 2, 1};
      vsignalErr = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
      vShapeSyst = {0.06, 0.05, 0.08, 0.15, 0.27, 0.3};
      transferFac = 0.198;
      transferFacErr = 0.009;*/
   // transferFacErr=0.3*transferFac;

   TH1D *hobsSR = new TH1D("hobsSR", "hobsSR", 6, 0, 6);
   TH1D *hobsCR = new TH1D("hobsCR", "hobsCR", 6, 0, 6);
   TH1D *hSignal = new TH1D("hSignal", "hSignal", 6, 0, 6);
   TH1D *hShapeSyst = new TH1D("hShapeSyst", "hShapeSyst", 6, 0, 6);
   for (size_t i = 0; i < 6; i++)
   {
      // vsignal[i]*=0.1;
      hobsSR->SetBinContent(i + 1, vdata[i]);
      hobsCR->SetBinContent(i + 1, vCRi[i]);
      hSignal->SetBinContent(i + 1, vsignal[i]);
      // hShapeSyst->SetBinContent(i + 1, vShapeSyst[i]);
   }

   // Create histogram templates for the background
   TH1 *h1_bSR = (TH1 *)hobsCR->Clone("h1_bSR");
   TH1 *h1_bCR = (TH1 *)hobsCR->Clone("h1_bCR");

   HistFactory::Measurement meas("HistModel", "HistModel");
   meas.SetPOI("mu");

   meas.SetLumi(1.0);
   meas.SetLumiRelErr(0.001);
   meas.AddConstantParam("Lumi");

   HistFactory::Channel channelSR("Signalregion");
   channelSR.SetData(hobsSR);

   // cout << (*vCrossSecErr)[2]/(*vCrossSecErr)[0]<<endl;
   // Create signal sample
   HistFactory::Sample signal("signal");
   signal.AddNormFactor("mu", 1, 0, 30);
   signal.AddOverallSys("signalTheoryConstraint", (*vCrossSecErr)[1] / (*vCrossSecErr)[0], (*vCrossSecErr)[2] / (*vCrossSecErr)[0]);
   signal.SetHisto(hSignal);
   channelSR.AddSample(signal);

   // Backgroundestimate in SR is CR scaled by the transferfactor
   RooStats::HistFactory::ShapeSys shapeSys_S;
   HistFactory::Sample backgSR("bA");
   backgSR.AddNormFactor("TransferFac", transferFac, transferFac - 10 * transferFacErr, transferFac + 10 * transferFacErr);
   // meas.AddConstantParam("TransferFac");
   backgSR.AddOverallSys("transFacConstraint", 1 - transferFacErr / transferFac, 1 + transferFacErr / transferFac);
   backgSR.SetHisto(h1_bSR);
   // shapeSys_S.SetErrorHist(hShapeSyst);
   // backgSR.AddShapeSys(shapeSys_S);

   channelSR.AddSample(backgSR);

   // Create channel for CR
   HistFactory::Channel channelCR("Bregion");
   channelCR.SetData(hobsCR);

   // Sample in CR
   HistFactory::Sample backgB("bB");
   backgB.SetHisto(h1_bCR);
   backgB.ActivateStatError();
   channelCR.AddSample(backgB);

   meas.AddChannel(channelSR);
   meas.AddChannel(channelCR);

   // meas.Print();

   // make the model
   RooWorkspace *w = HistFactory::MakeModelAndMeasurementFast(meas);

   RooStats::ModelConfig *sbModel = (RooStats::ModelConfig *)w->obj("ModelConfig");
   RooRealVar *poi = (RooRealVar *)sbModel->GetParametersOfInterest()->first();
   ModelConfig *bModel = (ModelConfig *)sbModel->Clone();
   bModel->SetName(TString(sbModel->GetName()) + TString("_with_poi_0"));
   double oldval = poi->getVal();
   poi->setVal(0);
   bModel->SetSnapshot(*poi);
   poi->setVal(oldval);
   w->import(*bModel);
   w->import(*sbModel);

   w->writeToFile("HistFactory.root");
   // w->Print();
}

int RunStats(int npoints, vector<double> &limits)
{
   RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
   const char *infile = "HistFactory.root";
   const char *wsName = "combined";
   const char *modelSBName = "ModelConfig";
   const char *modelBName = "ModelConfig_with_poi_0";
   const char *dataName = "obsData";
   int calculatorType = 2;
   int testStatType = 3;
   bool useCLs = true;
   // int npoints = 5000;
   double poimin = 0;
   double poimax = 300;
   vector<double> testlim;
   if (npoints > 0)
   {
      RunStats(0, testlim);
      poimin = testlim[0] / 10;
      poimax = testlim[0] * 10;
   }

   int ntoys = 1000;
   bool useNumberCounting = false;
   const char *nuisPriorName = 0;

   int mPrintLevel = -1;
   bool mEnableDetOutput = false;
   double poihat = 0;
   bool doFit = true;
   bool mOptimize = true;

   TFile *file = TFile::Open(infile);
   RooWorkspace *w = dynamic_cast<RooWorkspace *>(file->Get(wsName));

   RooAbsData *data = w->data(dataName);
   RooAbsData::setDefaultStorageType(RooAbsData::Vector);
   data->convertToVectorStore();
   ModelConfig *bModel = (ModelConfig *)w->obj(modelBName);
   ModelConfig *sbModel = (ModelConfig *)w->obj(modelSBName);

   sbModel->SetSnapshot(*sbModel->GetParametersOfInterest());

   RooArgSet initialParameters;
   RooArgSet *allParams = sbModel->GetPdf()->getParameters(*data);
   allParams->snapshot(initialParameters);
   delete allParams;

   const RooArgSet *poiSet = sbModel->GetParametersOfInterest();
   RooRealVar *poi = (RooRealVar *)poiSet->first();
   string mMinimizerType = ROOT::Math::MinimizerOptions::DefaultMinimizerType();

   if (doFit)
   {

      // do the fit : By doing a fit the POI snapshot (for S+B)  is set to the fit value
      // and the nuisance parameters nominal values will be set to the fit value.
      // This is relevant when using LEP test statistics
      RooArgSet constrainParams;
      if (sbModel->GetNuisanceParameters())
         constrainParams.add(*sbModel->GetNuisanceParameters());
      RooStats::RemoveConstantParameters(&constrainParams);

      RooFitResult *fitres = sbModel->GetPdf()->fitTo(
          *data, InitialHesse(false), Hesse(false), Minimizer(mMinimizerType.c_str(), "Migrad"), Strategy(0),
          PrintLevel(mPrintLevel), Constrain(constrainParams), Save(true), Offset(RooStats::IsNLLOffset()));
      if (fitres->status() != 0)
      {
         Warning("StandardHypoTestInvDemo",
                 "Fit to the model failed - try with strategy 1 and perform first an Hesse computation");
         fitres = sbModel->GetPdf()->fitTo(
             *data, InitialHesse(true), Hesse(false), Minimizer(mMinimizerType.c_str(), "Migrad"), Strategy(1),
             PrintLevel(mPrintLevel + 1), Constrain(constrainParams), Save(true), Offset(RooStats::IsNLLOffset()));
      }
      if (fitres->status() != 0)
         Warning("StandardHypoTestInvDemo", " Fit still failed - continue anyway.....");

      poihat = poi->getVal();
      std::cout << "StandardHypoTestInvDemo - Best Fit value : " << poi->GetName() << " = " << poihat << " +/- "
                << poi->getError() << std::endl;

      // save best fit value in the poi snapshot
      sbModel->SetSnapshot(*sbModel->GetParametersOfInterest());
   }

   SimpleLikelihoodRatioTestStat slrts(*sbModel->GetPdf(), *bModel->GetPdf());

   // null parameters must includes snapshot of poi plus the nuisance values
   RooArgSet nullParams(*sbModel->GetSnapshot());
   if (sbModel->GetNuisanceParameters())
      nullParams.add(*sbModel->GetNuisanceParameters());
   if (sbModel->GetSnapshot())
      slrts.SetNullParameters(nullParams);
   RooArgSet altParams(*bModel->GetSnapshot());
   if (bModel->GetNuisanceParameters())
      altParams.add(*bModel->GetNuisanceParameters());
   if (bModel->GetSnapshot())
      slrts.SetAltParameters(altParams);
   if (mEnableDetOutput)
      slrts.EnableDetailedOutput();

   // ratio of profile likelihood - need to pass snapshot for the alt
   RatioOfProfiledLikelihoodsTestStat ropl(*sbModel->GetPdf(), *bModel->GetPdf(), bModel->GetSnapshot());
   ropl.SetSubtractMLE(false);
   ropl.SetPrintLevel(mPrintLevel);
   ropl.SetMinimizer(mMinimizerType.c_str());
   if (mEnableDetOutput)
      ropl.EnableDetailedOutput();

   ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
   profll.SetOneSided(true);
   profll.SetMinimizer(mMinimizerType.c_str());
   profll.SetPrintLevel(mPrintLevel);
   if (mEnableDetOutput)
      profll.EnableDetailedOutput();

   profll.SetReuseNLL(mOptimize);
   slrts.SetReuseNLL(mOptimize);
   ropl.SetReuseNLL(mOptimize);

   if (mOptimize)
   {
      profll.SetStrategy(0);
      ropl.SetStrategy(0);
      ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
   }

   MaxLikelihoodEstimateTestStat maxll(*sbModel->GetPdf(), *poi);
   NumEventsTestStat nevtts;

   AsymptoticCalculator::SetPrintLevel(mPrintLevel);

   // create the HypoTest calculator class
   HypoTestCalculatorGeneric *hc = new AsymptoticCalculator(*data, *bModel, *sbModel, false);

   TestStatistic *testStat = &profll;
   ((AsymptoticCalculator *)hc)->SetOneSided(true);
   RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);

   HypoTestInverter calc(*hc);
   calc.SetConfidenceLevel(optHTInv.confLevel);

   calc.UseCLs(useCLs);
   calc.SetVerbose(true);

   if (npoints > 0)
   {
      if (poimin > poimax)
      {
         // if no min/max given scan between MLE and +4 sigma
         poimin = int(poihat);
         poimax = int(poihat + 4 * poi->getError());
      }
      calc.SetVerbose(0);
      calc.SetFixedScan(npoints, poimin, poimax);
   }

   HypoTestInverterResult *r = calc.GetInterval();
   std::cout << "The computed upper limit is: " << r->UpperLimit() << " +/- " << r->UpperLimitEstimatedError() << std::endl;

   std::cout << "Expected upper limits, using the B (alternate) model : " << std::endl;
   std::cout << " expected limit (median) " << r->GetExpectedUpperLimit(0) << std::endl;
   std::cout << " expected limit (-1 sig) " << r->GetExpectedUpperLimit(-1) << std::endl;
   std::cout << " expected limit (+1 sig) " << r->GetExpectedUpperLimit(1) << std::endl;
   std::cout << " expected limit (-2 sig) " << r->GetExpectedUpperLimit(-2) << std::endl;
   std::cout << " expected limit (+2 sig) " << r->GetExpectedUpperLimit(2) << std::endl;

   limits.push_back(r->UpperLimit());
   limits.push_back(r->GetExpectedUpperLimit(0));
   limits.push_back(r->GetExpectedUpperLimit(+1));
   limits.push_back(r->GetExpectedUpperLimit(-1));
   limits.push_back(r->GetExpectedUpperLimit(+2));
   limits.push_back(r->GetExpectedUpperLimit(-2));

   return 1;
}

int RunStats(vector<double> &limits)
{
   return RunStats(100, limits);
}

double TheoryXSection(double m_Gluino)
{
   return pow(10.0, 1.8821e-7 * m_Gluino * m_Gluino - 3.0554e-3 * m_Gluino + 2.3574);
}

double TheoryXSectionUpper(double m_Gluino)
{
   return pow(10.0, 2.1124e-7 * m_Gluino * m_Gluino - 3.1114e-3 * m_Gluino + 2.4538);
}

double TheoryXSectionLower(double m_Gluino)
{
   return pow(10.0, 2.3056e-7 * m_Gluino * m_Gluino - 3.2642e-3 * m_Gluino + 2.5049);
}

double TFuncCross(double *x, double *par)
{
   return TheoryXSection(x[0]);
}
double TFuncCrossUpper(double *x, double *par)
{
   return TheoryXSectionUpper(x[0]);
}
double TFuncCrossLower(double *x, double *par)
{
   return TheoryXSectionLower(x[0]);
}

void BrazilianPlot(vector<double> m_Go, vector<vector<double>> limits)
{

   auto c41 = new TCanvas("c41", "c41", 200, 10, 5000, 3000);
   gPad->SetTickx();
    gPad->SetTicky();
    c41->SetRightMargin(0.04);
    c41->SetLeftMargin(0.18);
    c41->SetBottomMargin(0.18);

   // double m[] = {0, 1, 2, 3, 4};
   vector<double> limit, expected, expected_m1sigma, expected_m2sigma, expected_p1sigma, expected_p2sigma, ex;
   for (size_t i = 0; i < limits.size(); i++)
   {
      limit.push_back(limits[i][0] * TheoryXSection(m_Go[i]));
      expected.push_back(limits[i][1] * TheoryXSection(m_Go[i]));
      expected_p1sigma.push_back((limits[i][2] - limits[i][1]) * TheoryXSection(m_Go[i]));
      expected_m1sigma.push_back((limits[i][1] - limits[i][3]) * TheoryXSection(m_Go[i]));
      expected_p2sigma.push_back((limits[i][4] - limits[i][1]) * TheoryXSection(m_Go[i]));
      expected_m2sigma.push_back((limits[i][1] - limits[i][5]) * TheoryXSection(m_Go[i]));
      ex.push_back(0.0);
   }

   auto ge = new TGraphAsymmErrors(limits.size(), &m_Go[0], &expected[0], &ex[0], &ex[0], &expected_m2sigma[0], &expected_p2sigma[0]);
   ge->SetTitle("");
   ge->GetXaxis()->SetTitleSize(0.07F);
   ge->GetXaxis()->SetTitle("m(#tilde{g}) [GeV]");
   ge->GetYaxis()->SetTitleSize(0.07F);
   ge->GetYaxis()->SetTitle("Corss section [pb]");
   ge->GetYaxis()->SetRangeUser(1e-4, 1e-1);
   ge->GetXaxis()->SetRangeUser(1310, 2490);
   ge->SetFillColor(kOrange+1);
   ge->Draw("SAME a3L");

   auto ge2 = new TGraphAsymmErrors(limits.size(), &m_Go[0], &expected[0], &ex[0], &ex[0], &expected_m1sigma[0], &expected_p1sigma[0]);
   ge2->SetFillColor(kGreen);
   ge2->Draw("SAME 3l");

   auto tg = new TGraph(limit.size(), &m_Go[0], &limit[0]);
   tg->SetLineWidth(2);

   auto tg2 = new TGraph(limit.size(), &m_Go[0], &expected[0]);
   tg2->SetLineStyle(kDashed);
   tg2->SetLineWidth(2);

   tg->Draw("L SAME");
   tg2->Draw("L SAME");

   TF1 *theory = new TF1("Theory", TFuncCross, m_Go[0], m_Go.back());
   TF1 *theoryupper = new TF1("Theory", TFuncCrossUpper, m_Go[0], m_Go.back());
   TF1 *theorylower = new TF1("Theory", TFuncCrossLower, m_Go[0], m_Go.back());

   theory->SetLineColor(kBlue);
   theoryupper->SetLineColor(kBlue);
   theorylower->SetLineColor(kBlue);

   theoryupper->SetLineStyle(kDashed);
   theorylower->SetLineStyle(kDashed);

   theory->Draw("SAME");
   theoryupper->Draw("SAME");
   theorylower->Draw("SAME");

   auto legend = new TLegend(0.48,0.6,0.9,0.9);
   legend->SetHeader("300 fb^{-1}");
   legend->AddEntry(tg,"observed");
   legend->AddEntry(tg2,"expected");
   legend->AddEntry(theory,"Theory #pm s.d._{theory}");
   legend->AddEntry(ge,"#pm 1 s.d.");
   legend->AddEntry(ge2,"#pm 2 s.d.");
   legend->Draw();
   gPad->SetLogy();
   gPad->SetTickx();
   gPad->SetTicky();
   gPad->RedrawAxis();
   c41->SaveAs("../Plots/Crosssections.pdf");
}

void RunStatsMacro()
{
   vector<double> m_Go = {1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500};
   vector<vector<double>> limits;
   vector<double> limit;
   for (size_t i = 0; i < m_Go.size(); i++)
   {
      HistFactoryModel(m_Go[i]);
      limit.clear();
      RunStats(limit);
      limits.push_back(limit);
   }
   TFile *fout = TFile::Open("StatsOutput.root", "RECREATE");
   fout->WriteObject(&limits, "limits");
   fout->WriteObject(&m_Go, "m_Go");
   fout->Close();
}

void DoPlots()
{
   TFile *fin = TFile::Open("StatsOutput.root", "READ");
   vector<vector<double>> limits;
   vector<double> m_Go;
   vector<vector<double>> *plimits;
   vector<double> *pm_Go;
   fin->GetObject("limits", plimits);
   fin->GetObject("m_Go", pm_Go);
   BrazilianPlot(*pm_Go,*plimits);
}

double SimpleSignificance(double m_Go)
{
   ostringstream signalOutFile;
   signalOutFile << "/home/gerrit/Documents/Projects/Root-Recast-CMS-SUS-19-013/RootFiles/Signal" << m_Go << "GeV.root";

   RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
   TFile f1(signalOutFile.str().c_str(), "read");
   TFile f2("/home/gerrit/Documents/Projects/Root-Recast-CMS-SUS-19-013/RootFiles/FinalSamples.root", "read");
   TVectorD *vSignalBins = (TVectorD *)f1.Get("SigBins");
   TVectorD *vSigErr = (TVectorD *)f1.Get("SigErr");
   TVectorD *vM_Go = (TVectorD *)f1.Get("M_Go");
   TVectorD *vCrossSecErr = (TVectorD *)f1.Get("CrossSecErr");
   TVectorD *vNCRi = (TVectorD *)f2.Get("NCR_i");
   TVectorD *vDataBins = (TVectorD *)f2.Get("DataBins");
   TVectorD *vTransferFac = (TVectorD *)f2.Get("TransferFac");

   cout << "+++++++++++++++++++++" << (*vM_Go)[0] << "GeV+++++++++++++++++++++++++" << endl;
   vector<double> vdata, vsignal, vsignalErr, vCRi, vCRiErr, vShapeSyst;
   double totalSI = 0.0;
   double totalBG= 0.0;
   for (size_t i = 5; i < 6; i++)
   {
      totalBG += (*vDataBins)[i];
      totalSI += (*vSignalBins)[i];
      vdata.push_back((*vDataBins)[i]);
      vsignal.push_back((*vSignalBins)[i]);
      vCRi.push_back((*vNCRi)[i]);
      vCRiErr.push_back(sqrt((*vNCRi)[i]));
      vsignalErr.push_back((*vSigErr)[i]);
   }

   double transferFac, transferFacErr;
   transferFac = (*vTransferFac)[0];
   transferFacErr = (*vTransferFac)[1];
   f1.Close();
   f2.Close();
   cout << totalSI/sqrt(totalBG) << endl;
   return totalSI/sqrt(totalBG);
}

void RunSimpleSign()
{
  vector<double> m_Go = {1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500};
  auto c41 = new TCanvas("c41", "c41", 200, 10, 2500, 1500);
  c41->SetRightMargin(0.04);
  c41->SetLeftMargin(0.18);
  c41->SetBottomMargin(0.18);
  //vector<vector<double>> limits;
  vector<double> significance;
  for (size_t i = 0; i < m_Go.size(); i++)
  {
     significance.push_back(SimpleSignificance(m_Go[i]));
  }
  double mgo[13];
  double sig[13];
  std::copy(m_Go.begin(),m_Go.end(), mgo);
  std::copy(significance.begin(),significance.end(), sig);
  auto g = new TGraph(m_Go.size(), mgo, sig);
  gPad->SetLogy();
  gPad->SetGrid();
  g->SetTitle("Significance 300 fb^{-1}");
  g->GetXaxis()->SetTitleSize(0.07F);
  g->GetXaxis()->SetTitle("m(#tilde{g}) [GeV]");
  g->GetYaxis()->SetTitleSize(0.07F);
  g->GetYaxis()->SetTitle("S / #sqrt{B}");
  //g->GetYaxis()->SetRangeUser(10, 1e-1);
  g->GetXaxis()->SetRangeUser(1310, 2490);
  g->Draw("AC*");
  c41->SaveAs("../Plots/Recast-SignificanceMETmin1200.pdf");
}


void HistFactoryModelMassSplitting(double m_Go)
{
   ostringstream signalOutFile;
   signalOutFile << "/home/gerrit/Documents/Projects/Root-Recast-CMS-SUS-19-013/RootFiles/SignalNeu2Mass" << m_Go << "GeV.root";

   RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
   TFile f1(signalOutFile.str().c_str(), "read");
   TFile f2("/home/gerrit/Documents/Projects/Root-Recast-CMS-SUS-19-013/RootFiles/FinalSamplesMassSplitt.root", "read");
   TVectorD *vSignalBins = (TVectorD *)f1.Get("SigBins");
   TVectorD *vSigErr = (TVectorD *)f1.Get("SigErr");
   TVectorD *vM_Go = (TVectorD *)f1.Get("M_Go");
   TVectorD *vCrossSecErr = (TVectorD *)f1.Get("CrossSecErr");
   TVectorD *vNCRi = (TVectorD *)f2.Get("NCR_i");
   TVectorD *vDataBins = (TVectorD *)f2.Get("DataBins");
   TVectorD *vTransferFac = (TVectorD *)f2.Get("TransferFac");

   cout << "+++++++++++++++++++++" << (*vM_Go)[0] << "GeV+++++++++++++++++++++++++" << endl;
   vector<double> vdata, vsignal, vsignalErr, vCRi, vCRiErr, vShapeSyst;

   for (size_t i = 0; i < 6; i++)
   {
      vdata.push_back((*vDataBins)[i]);
      vsignal.push_back((*vSignalBins)[i]);
      vCRi.push_back((*vNCRi)[i]);
      vCRiErr.push_back(sqrt((*vNCRi)[i]));
      vsignalErr.push_back((*vSigErr)[i]);
   }

   double transferFac, transferFacErr;
   transferFac = (*vTransferFac)[0];
   transferFacErr = (*vTransferFac)[1];
   f1.Close();
   f2.Close();
   /*
      vdata = {237, 67, 20, 3, 3, 1};
      vsignal = {3.5, 4.3, 6.6, 7.2, 7.2, 11.6};
      vCRi = {1191, 320, 112, 16, 2, 1};
      vsignalErr = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
      vShapeSyst = {0.06, 0.05, 0.08, 0.15, 0.27, 0.3};
      transferFac = 0.198;
      transferFacErr = 0.009;*/
   // transferFacErr=0.3*transferFac;

   TH1D *hobsSR = new TH1D("hobsSR", "hobsSR", 6, 0, 6);
   TH1D *hobsCR = new TH1D("hobsCR", "hobsCR", 6, 0, 6);
   TH1D *hSignal = new TH1D("hSignal", "hSignal", 6, 0, 6);
   TH1D *hShapeSyst = new TH1D("hShapeSyst", "hShapeSyst", 6, 0, 6);
   for (size_t i = 0; i < 6; i++)
   {
      // vsignal[i]*=0.1;
      hobsSR->SetBinContent(i + 1, vdata[i]);
      hobsCR->SetBinContent(i + 1, vCRi[i]);
      hSignal->SetBinContent(i + 1, vsignal[i]);
      // hShapeSyst->SetBinContent(i + 1, vShapeSyst[i]);
   }

   // Create histogram templates for the background
   TH1 *h1_bSR = (TH1 *)hobsCR->Clone("h1_bSR");
   TH1 *h1_bCR = (TH1 *)hobsCR->Clone("h1_bCR");

   HistFactory::Measurement meas("HistModel", "HistModel");
   meas.SetPOI("mu");

   meas.SetLumi(1.0);
   meas.SetLumiRelErr(0.001);
   meas.AddConstantParam("Lumi");

   HistFactory::Channel channelSR("Signalregion");
   channelSR.SetData(hobsSR);

   // cout << (*vCrossSecErr)[2]/(*vCrossSecErr)[0]<<endl;
   // Create signal sample
   HistFactory::Sample signal("signal");
   signal.AddNormFactor("mu", 1, 0, 300);
   signal.AddOverallSys("signalTheoryConstraint", (*vCrossSecErr)[1] / (*vCrossSecErr)[0], (*vCrossSecErr)[2] / (*vCrossSecErr)[0]);
   signal.SetHisto(hSignal);
   channelSR.AddSample(signal);

   // Backgroundestimate in SR is CR scaled by the transferfactor
   RooStats::HistFactory::ShapeSys shapeSys_S;
   HistFactory::Sample backgSR("bA");
   backgSR.AddNormFactor("TransferFac", transferFac, transferFac - 10 * transferFacErr, transferFac + 10 * transferFacErr);
   // meas.AddConstantParam("TransferFac");
   backgSR.AddOverallSys("transFacConstraint", 1 - transferFacErr / transferFac, 1 + transferFacErr / transferFac);
   backgSR.SetHisto(h1_bSR);
   // shapeSys_S.SetErrorHist(hShapeSyst);
   // backgSR.AddShapeSys(shapeSys_S);

   channelSR.AddSample(backgSR);

   // Create channel for CR
   HistFactory::Channel channelCR("Bregion");
   channelCR.SetData(hobsCR);

   // Sample in CR
   HistFactory::Sample backgB("bB");
   backgB.SetHisto(h1_bCR);
   backgB.ActivateStatError();
   channelCR.AddSample(backgB);

   meas.AddChannel(channelSR);
   meas.AddChannel(channelCR);

   // meas.Print();

   // make the model
   RooWorkspace *w = HistFactory::MakeModelAndMeasurementFast(meas);

   RooStats::ModelConfig *sbModel = (RooStats::ModelConfig *)w->obj("ModelConfig");
   RooRealVar *poi = (RooRealVar *)sbModel->GetParametersOfInterest()->first();
   ModelConfig *bModel = (ModelConfig *)sbModel->Clone();
   bModel->SetName(TString(sbModel->GetName()) + TString("_with_poi_0"));
   double oldval = poi->getVal();
   poi->setVal(0);
   bModel->SetSnapshot(*poi);
   poi->setVal(oldval);
   w->import(*bModel);
   w->import(*sbModel);

   w->writeToFile("HistFactory.root");
   // w->Print();
}


void Stats()
{
  RunSimpleSign();
  return;


  vector<double> massSplitting = {1900, 1800, 1600, 1250, 1000, 500, 250, 185};

  vector<double> limits;
  for(int i = 0; i < massSplitting.size(); i++)
  {
    HistFactoryModelMassSplitting(massSplitting[i]);
    RunStats(100, limits);
  }
  for (size_t i = 0; i < massSplitting.size(); i++) {
    cout << "MassSplitt " << massSplitting[i] << endl;
    for (size_t j = 0; j < 6; j++) {
      cout << limits[i*6+j] << endl;
      /* code */
    }
    /* code */
  }
  //RunSimpleSign();
  return;
   RunStatsMacro();
   DoPlots();
   // return;

   // BrazilianPlot(m_Go, limits);
}
/*
The computed upper limit is: 0.0662303 +/- 0
Expected upper limits, using the B (alternate) model :
 expected limit (median) 0.0662301
 expected limit (-1 sig) 0.0390162
 expected limit (+1 sig) 0.111859
 expected limit (-2 sig) 0.0255471
 expected limit (+2 sig) 0.177847
*/
