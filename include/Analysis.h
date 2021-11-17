#pragma once
#include <vector>
#include <iostream>
#include <fstream>

#include "PassedEvent.h"
#include "RootIO.h"

#include "ExRootAnalysis/ExRootClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"

#include <TClonesArray.h>
#include <TChain.h>
#include "classes/DelphesClasses.h"
#include <boost/progress.hpp>



class Analysis {
    private:
        static double sq(double x);
    public:
        static void AnalyseEvents(ExRootTreeReader *treeReader, vector<PassedEvent> &events , size_t MinNumEvents, int status);
        static double FindNormalisation(const char * filename);
        static void RunBigAnalysis();
};
