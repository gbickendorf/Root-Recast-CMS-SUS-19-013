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


void ReadGluino(const char *filename, double m_Go, const char *outFilename)
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
    if (events[i].mj1 > 70.0 && events[i].mj1 < 100.0 && events[i].mj2 > 70.0 && events[i].mj2 < 100.0)
      histSignal->Fill(events[i].ptmiss);
    /* code */
  }
  cout << histSignal->Integral() << endl;
  cout << integral << endl;
  histSignal->Scale(xSect[0] * 300e3 * acceptance / integral);
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
  RunSignif();
  return 0;

  vector<string> filenames = {"EventsRun1.root", "EventsRun2.root", "EventsRun3.root", "EventsRun4.root", "EventsRun5.root", "EventsRun6.root", "EventsRun7.root", "EventsRun8.root"};
  vector<double> m_n2 = {1900, 185, 1800, 1600, 1250, 1000, 500, 250};
  for (size_t n = 1; n < 8; n++)
  {
    cout << "+++++++" << filenames[n] << "++++++++++++" << endl;
    cout << m_n2[n] << endl;
    ostringstream signalOutFile;
    signalOutFile << "RootFiles/SignalNeu2Mass" << m_n2[n] << "GeV.root";
    ReadGluino((string("~/Documents/PHD/SUSY4Cathode/MassSplittings/Events/") + filenames[n]).c_str(), 1900, signalOutFile.str().c_str());
  }

  vector<PassedEvent> events, SREvents;
  Analysis::SampleFromEventsNewPassedEventDef(events);
  Plots::PlotPTShape(events);
  Plots::PlotPTMiss(events);

  vector<double> params;
  double bNorm, transferFactor, bNormErr, transferFactorErr;
  vector<double> NCRi, SRBinContent;
  Analysis::FitCherbyshev(events, bNorm, bNormErr, params);
  Analysis::CalcTransferFactor(events, bNorm, bNormErr, NCRi, transferFactor, transferFactorErr);
  Plots::PlotMJ1(events, params);
  Plots::PlotSignalRegion(events, NCRi, transferFactor, transferFactorErr);
  Plots::PlotPhotonLeptonValidation(events);
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

  Analysis::GetEventsInSignalRegion(events, SREvents);
  Analysis::GetSRBinContent(SREvents, SRBinContent);
  pFile = fopen("DataMassSplitt.csv", "w");
  fprintf(pFile, "Data,stat\n");
  for (size_t i = 0; i < SRBinContent.size(); i++)
  {
    vSRBin[i] = SRBinContent[i];
    fprintf(pFile, "%E,%E\n", SRBinContent[i], sqrt(SRBinContent[i]));
    /* code */
  }
  fclose(pFile);
  TFile f("RootFiles/FinalSamplesMassSplitt.root", "RECREATE");
  vNCRi.Write("NCR_i");
  vSRBin.Write("DataBins");
  vtransferFactor.Write("TransferFac");
  f.Close();
}
