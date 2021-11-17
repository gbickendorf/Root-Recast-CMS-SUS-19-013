#pragma once
#include "PassedEvent.h"
#include <vector>
#include "TVectorD.h"
#include <TROOT.h>
#include <TFile.h>
using namespace std;

class RootIO {
  public:  
    static void SaveEvents(vector<PassedEvent> events);
    static void ReadEvents(vector<PassedEvent> &events);
};