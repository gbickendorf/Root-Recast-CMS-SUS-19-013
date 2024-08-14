#include "PassedEvent.h"
#include "RootIO.h"
#include "Plots.h"
#include "Analysis.h"
#include "Stats.h"

#include "ExRootAnalysis/ExRootClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"
#include <iostream>
#include "TROOT.h"
#include <omp.h>
#include <unistd.h>
#include <TLegend.h>
#include <vector>

//#include <RooHist.h>
//#include <RooDataSet.h>
//#include <RooPlot.h>

using namespace std;

void ReadFilenames(const char *filename, vector<string> &FilePaths)
{
  std::ifstream file(filename);
  if (file.is_open())
  {
    std::string line;
    while (std::getline(file, line))
    {
      FilePaths.push_back(line);
    }
    file.close();
  }
}
void test()
{
  vector<string> paths;
  RootIO::GetRootFilePath("/media/gerrit/Files/RootFiles/A/ROOTFILES.txt", paths, -1);
  vector<PassedEvent> events;
  int total = 0;
  int inFile;
  for (size_t i = 0; i < paths.size(); i++)
  {
    RootIO::ReadEvents(paths[i].c_str(), events, 0, inFile);
    total += inFile;
  }
  cout << total << endl;
  RootIO::SaveEvents("/media/gerrit/Files/RootFiles/A/root.root", events, 1);
  events.clear();
  RootIO::ReadEvents("/media/gerrit/Files/RootFiles/A/root.root", events, 1);
  events.clear();
  RootIO::ReadEvents("RootFiles/A-Total.root", events, 1);
}

void ReadGluino(const char *filename, double m_Go, const char *outFilename, double mjmin, double mjmax, double intLumi)
{
  cout << "Mass " << m_Go;
  TVectorD xSect(3);
  xSect[0] = Analysis::TheoryXSection(m_Go);
  xSect[1] = Analysis::TheoryXSectionLower(m_Go);
  xSect[2] = Analysis::TheoryXSectionUpper(m_Go);

  TChain chain("Delphes");
  chain.Add(filename);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  vector<PassedEvent> events;
  int nTotal;
  double acceptance = Analysis::AnalyseEvents(treeReader, events, -1, 100);
  const double bins[7] = {300.0, 450.0, 600.0, 800.0, 1000.0, 1200.0, 20000.0};
  TH1 *histSignal = new TH1F("histSignal", "ptmiss CR", 6, bins);
  double integral = 0.0;
  for (size_t i = 0; i < events.size(); i++)
  {
    if (events[i].ptmiss > 300 && events[i].ptmiss < 20000)
      integral += 1.0;
    if (events[i].mj1 > mjmin && events[i].mj1 < mjmax && events[i].mj2 > mjmin && events[i].mj2 < mjmax)
      histSignal->Fill(events[i].ptmiss);
    /* code */
  }
  cout << histSignal->Integral() << endl;
  cout << integral << endl;
  histSignal->Scale(xSect[0] * intLumi * acceptance / integral);
  histSignal->Draw();
  FILE *pFile = fopen("Signal.csv", "w");
  fprintf(pFile, "Signal,stat\n");
  TVectorD binC(6);
  TVectorD binErr(6);
  TVectorD M(1);
  M[0] = m_Go;

  for (size_t i = 1; i < 7; i++)
  {
    binC[i - 1] = histSignal->GetBinContent(i);
    binErr[i - 1] = histSignal->GetBinError(i);
    cout << histSignal->GetBinCenter(i) << "   " << histSignal->GetBinContent(i) << " +- " << histSignal->GetBinError(i) << endl;
    fprintf(pFile, "%E,%E\n", histSignal->GetBinContent(i), histSignal->GetBinError(i));
    /* code */
  }
  fclose(pFile);

  TFile f(outFilename, "RECREATE");
  binC.Write("SigBins");
  binErr.Write("SigErr");
  M.Write("M_Go");
  xSect.Write("CrossSecErr");
  f.Close();
  delete treeReader;
  delete histSignal;
}

void CalcSignif(const char *filenameGo, double m_Go)
{
  double Oversample = 100.0;
  double OversampleBG= 1.0/0.485507;
  double intLumi=300.0e3;
  TVectorD xSect(3);
  xSect[0] = Analysis::TheoryXSection(m_Go);
  xSect[1] = Analysis::TheoryXSectionLower(m_Go);
  xSect[2] = Analysis::TheoryXSectionUpper(m_Go);
  cout << m_Go<<"  GeV"<< endl;
  vector<PassedEvent> BGEvents, SIEvents, BGEventsSR, SIEventsSR;
  TChain chain("Delphes");
  Analysis::SampleFromEventsNewPassedEventDef(BGEvents,intLumi*OversampleBG);

  chain.Add(filenameGo);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  
  int nTotal;
  double acceptance = Analysis::AnalyseEvents(treeReader, SIEvents, int(xSect[0]*intLumi*Oversample), 100);

  const double bins[7] = {300.0, 450.0, 600.0, 800.0, 1000.0, 1200.0, 20000.0};
  TH1 *histSignal = new TH1F("histSignal", "ptmiss CR", 6, bins);
  double integral = 0.0;
  for (size_t i = 0; i < SIEvents.size(); i++)
  {
    PassedEvent event =SIEvents[i];
    if (event.mj1 > 70.0 && event.mj1 < 100.0 && event.mj2 > 70.0 && event.mj2 < 100.0)
      SIEventsSR.push_back(event);
    /* code */
  }
  for (size_t i = 0; i < BGEvents.size(); i++)
  {
    PassedEvent event =BGEvents[i];
    if (event.mj1 > 70.0 && event.mj1 < 100.0 && event.mj2 > 70.0 && event.mj2 < 100.0)
      BGEventsSR.push_back(event);
    /* code */
  }
  for (size_t i = 0; i < 18; i++)
  {
    double metMIn = 300.0+i*100;
    int S = 0;
    int B = 0;
    for (size_t i = 0; i < SIEventsSR.size(); i++)
    {
      if(SIEventsSR[i].ptmiss > metMIn)
        S++;
    }
    for (size_t i = 0; i < BGEventsSR.size(); i++)
    {
      if(BGEventsSR[i].ptmiss > metMIn)
        B++;
    }
    cout << "Met > " << metMIn <<" S: " << S/Oversample << " B: " << B/OversampleBG << " Signif: " << S/sqrt(B/OversampleBG)/Oversample << endl;
  }
  
  return;

}

void CompleteAnalysis(vector<PassedEvent> &events, double mjmin, double mjmax)
{
  vector<PassedEvent> rawEvents = events;
  int cutLoop;
  int passed = 0;
  events.clear();
  for (size_t i = 0; i < rawEvents.size(); i++)
  {
    if (rawEvents[i].jAK8PT.size() < 2)
      continue;
    if (rawEvents[i].jAK8PT[0] < 200)
      continue;
    if (rawEvents[i].jAK8PT[1] < 200)
      continue;

    if (rawEvents[i].mj1 < mjmin-30 || rawEvents[i].mj1 > mjmax+40)
      continue;
    if (rawEvents[i].mj2 < mjmin-30 || rawEvents[i].mj2 > mjmax+40)
      continue;
    cutLoop = 0;
    for (size_t j = 0; j < rawEvents[i].jAK8AngSepBTag.size() / 2; j++)
    {
      if (rawEvents[i].jAK8AngSepBTag[j * 2 + 1] < 0.8)
        cutLoop++;
    }
    cutLoop=0;
    if (cutLoop)
      continue;

    passed++;
    events.push_back(rawEvents[i]);
  }
  cout << "Passed : " << passed << endl;
}
void OldMain(int argc, char **argv)
{
  double mjmin=70;
  double mjmax=100;
  vector<string> filenames = {"Run02.root", "Run03.root", "Run04.root", "Run05.root", "Run06.root", "Run07.root", "Run08.root", "Run09.root", "Run10.root", "Run11.root", "Run12.root", "Run13.root", "Run14.root"};
  vector<double> m_Go = {1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500};
  for (size_t n = 13; n < 13; n++)
  {
    //cout << "+++++++" << filenames[n] << "++++++++++++" << endl;
    ostringstream signalOutFile,signalinFile;
    signalinFile << "/media/gerrit/Files/RootFiles/Gluino/Jump/"<< filenames[n];
    cout << signalinFile.str() << endl;
    signalOutFile << "RootFiles/Gluino" << m_Go[n] << "GeV.root";

    ReadGluino(signalinFile.str().c_str(), m_Go[n], signalOutFile.str().c_str(),mjmin,mjmax,300e3);
    //return;
  }

  vector<PassedEvent> events, SREvents;
  Analysis::SampleFromEventsNewPassedEventDef(events,137.0e3);
  CompleteAnalysis(events,mjmin, mjmax);
  /*14556
Done Reading
Weight: 0.46997  Nexpected: 1073.41  Navail: 8915
Weight: 0.485507  Nexpected: 683.108  Navail: 3785
Weight: 0.406745  Nexpected: 476.706  Navail: 1856
Total events : 6833
*/
  //return;
  Plots::PlotPTShape(events, mjmin, mjmax);
  Plots::PlotPTMiss(events);
  return;

  vector<double> params;
  double bNorm, transferFactor, bNormErr, transferFactorErr;
  vector<double> NCRi, SRBinContent;
  Analysis::FitCherbyshev(events, bNorm, bNormErr, params, mjmin, mjmax);
  Analysis::CalcTransferFactor(events, bNorm, bNormErr, NCRi, transferFactor, transferFactorErr, mjmin, mjmax);
  Plots::PlotMJ1(events, params, mjmin, mjmax);
  Plots::PlotSignalRegion(events, NCRi, transferFactor, transferFactorErr, mjmin, mjmax);
  //Plots::PlotPhotonLeptonValidation(events, mjmin, mjmax);
  FILE *pFile = fopen("BackgroundEstimate.csv", "w");
  fprintf(pFile, "BG,stat,sys\n");

  TVectorD vNCRi(6);
  TVectorD vSRBin(6);
  TVectorD vtransferFactor(2);

  vtransferFactor[0] = transferFactor;
  vtransferFactor[1] = transferFactorErr;

  for (size_t i = 0; i < NCRi.size(); i++)
  {
    vNCRi[i] = NCRi[i];
    fprintf(pFile, "%E,%E,%E\n", NCRi[i] * transferFactor, sqrt(NCRi[i]) * transferFactor, NCRi[i] * transferFactorErr);
    /* code */
  }
  fclose(pFile);

  Analysis::GetEventsInSignalRegion(events, SREvents, mjmin, mjmax);
  Analysis::GetSRBinContent(SREvents, SRBinContent);
  pFile = fopen("Data.csv", "w");
  fprintf(pFile, "Data,stat\n");
  for (size_t i = 0; i < SRBinContent.size(); i++)
  {
    vSRBin[i] = SRBinContent[i];
    fprintf(pFile, "%E,%E\n", SRBinContent[i], sqrt(SRBinContent[i]));
    /* code */
  }
  fclose(pFile);
  TFile f("RootFiles/FinalSamples.root", "RECREATE");
  vNCRi.Write("NCR_i");
  vSRBin.Write("DataBins");
  vtransferFactor.Write("TransferFac");
  f.Close();
}

void RunSignif()
{
  vector<string> filenames = {"Run02.root", "Run03.root", "Run04.root", "Run05.root", "Run06.root", "Run07.root", "Run08.root", "Run09.root", "Run10.root", "Run11.root", "Run12.root", "Run13.root", "Run14.root"};
  vector<double> m_Go = {1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500};
  for (size_t n = 6; n < 13; n++)
  {
    cout << "+++++++" << filenames[n] << "++++++++++++" << endl;
    ostringstream signalOutFile;
    signalOutFile << "RootFiles/Signal" << m_Go[n] << "GeV.root";
    CalcSignif((string("/media/gerrit/Files/RootFiles/Gluino/Jump/") + filenames[n]).c_str(), m_Go[n]);
    return;
  }
}

int main(int argc, char **argv)
{
  double mjmin=70;
  double mjmax=100;
  OldMain(argc, argv);
  //RunSignif();
  return 0;

  vector<string> filenames = {"EventsRun1.root", "EventsRun2.root", "EventsRun3.root", "EventsRun4.root", "EventsRun5.root", "EventsRun6.root", "EventsRun7.root", "EventsRun8.root"};
  vector<double> m_n2 = {1900, 185, 1800, 1600, 1250, 1000, 500, 250};
  for (size_t n = 1; n < 8; n++)
  {
    cout << "+++++++" << filenames[n] << "++++++++++++" << endl;
    cout << m_n2[n] << endl;
    ostringstream signalOutFile;
    signalOutFile << "RootFiles/SignalNeu2Mass" << m_n2[n] << "GeV.root";
    ReadGluino((string("~/Documents/PHD/SUSY4Cathode/MassSplittings/Events/") + filenames[n]).c_str(), 1900, signalOutFile.str().c_str(),70,100,300e3);
  }

  vector<PassedEvent> events, SREvents;
  Analysis::SampleFromEventsNewPassedEventDef(events);
  Plots::PlotPTShape(events, mjmin, mjmax);
  Plots::PlotPTMiss(events);

  vector<double> params;
  double bNorm, transferFactor, bNormErr, transferFactorErr;
  vector<double> NCRi, SRBinContent;
  Analysis::FitCherbyshev(events, bNorm, bNormErr, params, mjmin, mjmax);
  Analysis::CalcTransferFactor(events, bNorm, bNormErr, NCRi, transferFactor, transferFactorErr, mjmin, mjmax);
  Plots::PlotMJ1(events, params, mjmin, mjmax);
  Plots::PlotSignalRegion(events, NCRi, transferFactor, transferFactorErr, mjmin, mjmax);
  Plots::PlotPhotonLeptonValidation(events, mjmin, mjmax);
  FILE *pFile = fopen("BackgroundEstimateMassSplitt.csv", "w");
  fprintf(pFile, "BG,stat,sys\n");

  TVectorD vNCRi(6);
  TVectorD vSRBin(6);
  TVectorD vtransferFactor(2);

  vtransferFactor[0] = transferFactor;
  vtransferFactor[1] = transferFactorErr;

  for (size_t i = 0; i < NCRi.size(); i++)
  {
    vNCRi[i] = NCRi[i];
    fprintf(pFile, "%E,%E,%E\n", NCRi[i] * transferFactor, sqrt(NCRi[i]) * transferFactor, NCRi[i] * transferFactorErr);
    /* code */
  }
  fclose(pFile);
  /*

  Analysis::GetEventsInSignalRegion(events, SREvents, mjmin, mjmax);
  Analysis::GetSRBinContent(SREvents, SRBinContent);
  pFile = fopen("DataMassSplitt.csv", "w");
  fprintf(pFile, "Data,stat\n");
  for (size_t i = 0; i < SRBinContent.size(); i++)
  {
    vSRBin[i] = SRBinContent[i];
    fprintf(pFile, "%E,%E\n", SRBinContent[i], sqrt(SRBinContent[i]));
    /* code */
  /*}
  fclose(pFile);
  TFile f("RootFiles/FinalSamplesMassSplitt.root", "RECREATE");
  vNCRi.Write("NCR_i");
  vSRBin.Write("DataBins");
  vtransferFactor.Write("TransferFac");
  f.Close();
  */
  
}
