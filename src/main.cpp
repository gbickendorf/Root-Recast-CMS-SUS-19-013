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
#include <omp.h>
#include <unistd.h>

#include <vector>
using namespace std;

int main()
{
  
  
  cout<<Analysis::FindNormalisation("/mnt/5c451946-c244-49ab-9fc5-1aaca8739b2a/ZNorm-Cluster/ROOTFILES.txt")/Analysis::FindNormalisation("/media/gerrit/Files/DelphesEvents/Z/ROOTFILES.txt");
  //cout<<Analysis::FindNormalisation("/media/gerrit/Files/DelphesEvents/Z/ROOTFILES.txt");
  return 0;
  /*
  return 0;///Analysis::FindNormalisation("/media/gerrit/Files/DelphesEvents/PythiaDelphesCards/q70.txt")<<endl;
  cout<<Analysis::FindNormalisation("/media/gerrit/Files/DelphesEvents/PythiaDelphesCards/Norm25.txt")/Analysis::FindNormalisation("/media/gerrit/Files/DelphesEvents/PythiaDelphesCards/q90.txt")<<endl;
  cout<<Analysis::FindNormalisation("/media/gerrit/Files/DelphesEvents/PythiaDelphesCards/Norm25.txt")/Analysis::FindNormalisation("/media/gerrit/Files/DelphesEvents/PythiaDelphesCards/ZXQ50.txt")<<endl;
  return 0;*/
  //Analysis::FindNormalisation("/home/gerrit/Documents/tmp/DelphesEventsIso/ROOTFILES.txt");
  //cout<<Analysis::FindNormalisation("/mnt/5c451946-c244-49ab-9fc5-1aaca8739b2a/ZNorm-Cluster/ROOTFILES.txt")/Analysis::FindNormalisation("/mnt/5c451946-c244-49ab-9fc5-1aaca8739b2a/Z-Cluster/ROOTFILES.txt")<<endl;
  //cout<<"Ideal : 0.002"<<endl;
  //return 0;

  //cpu_timer timers;
  
  int i = 0;
  vector<PassedEvent> events;
  //RootIO::ReadEvents("Normtest.root", events);
  vector<string> smallFiles={"/home/gerrit/Documents/tmp/Z-ClusterSmall/ROOTFILES.txt","/home/gerrit/Documents/tmp/Z-ClusterNormSmall/ROOTFILES.txt"};
  vector<string> filename = {"/mnt/5c451946-c244-49ab-9fc5-1aaca8739b2a/Z-Cluster/ROOTFILES.txt", "/media/gerrit/Files/DelphesEvents/W-Cluster/ROOTFILES.txt", "/media/gerrit/Files/DelphesEvents/tt-Cluster/ROOTFILES.txt"};
  filename.push_back("/mnt/5c451946-c244-49ab-9fc5-1aaca8739b2a/ZNorm-Cluster/ROOTFILES.txt");
  filename.push_back("/mnt/5c451946-c244-49ab-9fc5-1aaca8739b2a/WNorm-Cluster/ROOTFILES.txt");
  filename.push_back("/mnt/5c451946-c244-49ab-9fc5-1aaca8739b2a/ttNorm-Cluster/ROOTFILES.txt");
  filename.push_back("/home/gerrit/Documents/tmp/DelphesEventsIso/ROOTFILES.txt");
  filename.push_back("/media/gerrit/Files/DelphesEvents/PythiaDelphesCards/HT200.txt");
  filename.push_back("/media/gerrit/Files/DelphesEvents/Z/ROOTFILES.txt");
  for (i = 8; i < 9; i++)
  {
    vector<string> rootfiles;
    rootfiles.clear();
    RootIO::GetRootFilePath(filename[i].c_str(), rootfiles, 40);
    TChain chain("Delphes");
    cout << "Analyse all Events in " << filename[i] << endl;
    for (size_t i = 0; i < rootfiles.size(); i++)
    {
      chain.Add(rootfiles[i].c_str());
    }

    ExRootTreeReader *reader = new ExRootTreeReader(&chain);
    Analysis::AnalyseEventsNew(reader, events, -1, i);
    cout << events.size() << endl;
    //RootIO::SaveEvents("Normtest2.root", events);
  }
  events.clear();
  return 0;
  RootIO::ReadEvents("Normtest.root", events);
  TCanvas *c1 = new TCanvas("cLepPhoVal", "c1");
  c1->cd();
  gPad->SetLogy();
  TH1D *hist = new TH1D("ptmiss", "pt", 20, 100, 300);
  for (i = 0; i < events.size(); i++)
  {
    hist->Fill(events[i].ptmiss);
  }
  
  TF1 *f1 = (TF1 *)gROOT->GetFunction("expo");
  hist->Fit(f1);
  cout << f1->Eval(300)<<endl;
  hist->Draw();
  c1->SaveAs("ptmisssssss.eps");
  //cout << timer.format()<<endl;
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
  
  vector<PassedEvent> events;
  RootIO::ReadEvents("Normtest.root",events);
  Plots::PlotPTShape(events);
  //return 0;
  Plots::PlotPTMiss(events);
  
  vector<double> params;
  double bNorm,transferFactor;
  vector<double> bi;
  Analysis::FitCherbyshev(events,bNorm, params);
  Analysis::CalcTransferFactor(events,bNorm,bi,transferFactor);
  Plots::PlotMJ1(events,params);
  Plots::PlotSignalRegion(events,bi);
  Plots::PlotPhotonLeptonValidation(events);
  //Plot1(events);  
  */
}
