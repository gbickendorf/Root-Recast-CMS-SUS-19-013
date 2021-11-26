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
 // Plots::PlotAll();
  vector<string> files;
  RootIO::GetRootFilePath("/mnt/5c451946-c244-49ab-9fc5-1aaca8739b2a/Z-Cluster/ROOTFILES.txt",files,5);
  TChain testchain("Delphes");
  for (size_t i = 0; i < files.size(); i++)
  {
    testchain.Add(files[i].c_str());
  }
  ExRootTreeReader *reader = new ExRootTreeReader(&testchain);
  vector<PassedEvent> events;
  //Analysis::AnalyseEventsSinglePhotonSample(reader,events,-1);
  Analysis::AnalyseEvents(reader,events, -1,0);
  cout << events.size()<<endl;
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
  RootIO::ReadEvents(events);
  Plots::PlotPTShape(events);
  return 0;
  Plots::PlotPTMiss(events);
  
  vector<double> params;
  double bNorm,transferFactor;
  vector<double> bi;
  Analysis::FitCherbyshev(events,bNorm, params);
  Analysis::CalcTransferFactor(events,bNorm,bi,transferFactor);
  Plots::PlotMJ1(events,params);
  Plots::PlotSignalRegion(events,bi);
  //Plot1(events);  
  */
}
