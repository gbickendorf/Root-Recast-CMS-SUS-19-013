#pragma once
#include "PassedEvent.h"
#include <vector>
#include "TVectorD.h"
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;

class RootIO {
  public:  
    static void SaveEvents(const char * filename, vector<PassedEvent> events,int verbose);
    static void ReadEvents(const char * filename,vector<PassedEvent> &events,int verbose);
    static void SaveEvents(const char * filename, vector<PassedEvent> events);
    static void ReadEvents(const char * filename,vector<PassedEvent> &events);
    static void ReadEvents(const char * filename, vector<PassedEvent> &events,int verbose, int &nTotal);
    static void GetRootFiles(const char * filename,TChain &chain);
    static void GetRootFiles(const char * filename,TChain &chain, int nFiles);
    static void GetRootFilePath(const char * filename,vector<string> &rootFilePaths, int nFiles);
};