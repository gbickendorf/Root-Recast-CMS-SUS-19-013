#include "PassedEvent.h"
#include "RootIO.h"
#include "Plots.h"
#include "Analysis.h"

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

void TestNorm()
{
  TCanvas *c1 = new TCanvas("c1Trans", "c1");
  c1->cd();
  TH1F *hist = new TH1F("ht", "ht", 10, 0, 1000);
  TH1F *histNorm = new TH1F("htNorm", "ht", 10, 0, 1000);
  vector<double> metZ, metZNorm;
  Analysis::FindNormalisation("/media/gerrit/Files/DelphesEvents/W/ROOTFILES.txt", 10, metZ);
  Analysis::FindNormalisation("/media/gerrit/Files/DelphesEvents/WNorm/ROOTFILES.txt", 100, metZNorm);
  for (size_t i = 0; i < metZ.size(); i++)
  {
    hist->Fill(metZ[i]);
  }
  for (size_t i = 0; i < metZNorm.size(); i++)
  {
    histNorm->Fill(metZNorm[i]);
    /* code */
  }
  hist->Scale(1.0 / hist->Integral());
  histNorm->Scale(1.0 / histNorm->Integral());
  gPad->SetLogy();
  hist->SetLineColor(kRed);
  hist->Draw();
  histNorm->SetLineColor(kGreen);
  histNorm->Draw("SAME");
  auto legend = new TLegend(0.33, 0.80, 0.48, 0.9);
  legend->AddEntry(hist, "W");
  legend->AddEntry(histNorm, "WNorm");
  legend->Draw();
  c1->SaveAs("METCompareWMetMin.eps");
  gPad->SetLogy(0);
  gPad->SetGrid();
  hist->Sumw2();
  histNorm->Sumw2();
  hist->Divide(histNorm);
  hist->Draw();
  c1->SaveAs("DifferencesWMetMin.eps");
}

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
  RootIO::GetRootFilePath("/media/gerrit/Files/RootFiles/A/ROOTFILES.txt",paths,-1);
  vector<PassedEvent> events;
  int total=0;
  int inFile;
  for (size_t i = 0; i < paths.size(); i++)
  {
      RootIO::ReadEvents(paths[i].c_str(),events,0,inFile);
      total+=inFile;
  }
  cout << total << endl;
  RootIO::SaveEvents("/media/gerrit/Files/RootFiles/A/root.root",events,1);
  events.clear();
  RootIO::ReadEvents("/media/gerrit/Files/RootFiles/A/root.root",events,1);
  events.clear();
  RootIO::ReadEvents("RootFiles/A-Total.root",events,1);
  
}

int main(int argc, char **argv)
{
  vector<PassedEvent> events;
  Analysis::SampleFromEvents(events);
  Plots::PlotPTShape(events);
  //return 0;
  Plots::PlotPTMiss(events);

  vector<double> params;
  double bNorm, transferFactor;
  vector<double> bi;
  Analysis::FitCherbyshev(events, bNorm, params);
  Analysis::CalcTransferFactor(events, bNorm, bi, transferFactor);
  Plots::PlotMJ1(events, params);
  Plots::PlotSignalRegion(events, bi, transferFactor);
  Plots::PlotPhotonLeptonValidation(events);
  //Plot1(events);
  return 0;
  //RootIO::SaveEvents("Alltt.root",events);
  // Plots::PlotAll();
  //Analysis::RunBigAnalysis();
  /*
  vector<vector<PassedEvent>> pevents(omp_get_max_threads());
  vector<TChain*> pChains;
  for (size_t i = 0; i < omp_get_max_threads(); i++)
  {
    pChains.push_back(new TChain("Delphes"));
  }
  vector<ExRootTreeReader*> ptreeReader;
  vector<PassedEvent> events2;
  TChain chainLepton("Delphes");
  int numttEvents = 523456*4;
  int eventsInParallelLoop=0;
  int breaker=10000;
  vector<string> filenames;
  int n = filenames.size();
  RootIO::GetRootFilePath("/media/gerrit/Files/DelphesEvents/tt-Cluster/ROOTFILES.txt",filenames,-1);
  TChain testchain("Delphes");
  int overshotIndex=0;
  for (int i = 0; i < filenames.size(); i++)
  {
    testchain.Add(filenames[i].c_str());
    ExRootTreeReader *testreader = new ExRootTreeReader(&testchain);
    pChains[i%omp_get_max_threads()]->Add(filenames[i].c_str());
    eventsInParallelLoop+=testreader->GetEntries();
    if(eventsInParallelLoop>numttEvents)
    {
      overshotIndex=i%omp_get_max_threads();
      break;
    }
    testchain.Reset();
    delete testreader;
  }
  
  int totalEvents=0;
  #pragma omp parallel for
  for (size_t i = 0; i < omp_get_max_threads(); i++)
  {
    ExRootTreeReader *reader = new ExRootTreeReader((pChains[i]));
    if(!reader->GetEntries())
      continue;
    if(overshotIndex==i)
    {
      Analysis::AnalyseEventsSingleLeptonSample(reader,pevents[i],reader->GetEntries()-eventsInParallelLoop+numttEvents);
      totalEvents+=reader->GetEntries()-eventsInParallelLoop+numttEvents;
    }
    else
    {
      Analysis::AnalyseEventsSingleLeptonSample(reader,pevents[i],-1);
      totalEvents+=reader->GetEntries();
    }
    delete reader;
  }
  int nPassed=0;
  for (size_t i = 0; i < omp_get_max_threads(); i++)
  {
    nPassed+=pevents[i].size();
  }
  cout << nPassed << endl;
  cout << totalEvents<<endl;
  return 0;
  RootIO::GetRootFiles("/media/gerrit/Files/DelphesEvents/tt-Cluster/ROOTFILES.txt",chainLepton);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chainLepton);
  Analysis::AnalyseEventsSingleLeptonSample(treeReader,events2,numttEvents);
  chainLepton.Reset(); 
  return 0;
  //return 0;
  //cout << Analysis::FindNormalisation("/media/gerrit/Files/DelphesEvents/tt-Cluster/ROOTFILES.txt") <<endl;
  //cout << 0.00438544 << endl;
  //return 0;
  

  */
}
