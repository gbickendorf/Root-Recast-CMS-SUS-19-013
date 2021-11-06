#include <TROOT.h>
#include <iostream>
#include <TCanvas.h>
#include <TChain.h>
#include <TH1.h>
#include "TClonesArray.h"
#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootClasses.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"
using namespace std;

void exampleMacro()
{
	TCanvas *c1 = new TCanvas("c1","c1");
 c1->cd();
  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add("/home/gerrit/Documents/PHD/Tools/Delphes-3.5.0/delphesout.root");

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("FatJet");
	TClonesArray *weights = treeReader->UseBranch("Weight");
	TClonesArray *met = treeReader->UseBranch("MissingET");
  // Book histograms
  TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 100, 0.0, 1000.0);
  const double bins[7]= {300.0,450.0,600.0,800.0,1000.0,1200.0,2000.0};
  TH1 *histMet = new TH1F("ptmiss", "ptmiss", 6,bins);
  //gPad->SetLogy();

  // Loop over all events
  int passesOneJetCut=0;
  int passesAllJetCuts=0;
  int passesWeight=0;
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
	Weight *weight = (Weight*) weights->At(0);
	if(weight->Weight==0)
	{
		continue;
	}
	//printf("%E\n",weight->Weight);
    // If event contains at least 1 jet
    double metval=((MissingET*)met->At(0))->MET;
    if(branchJet->GetEntries() > 1)
    {
      // Take first jet
      Jet *jet = (Jet*) branchJet->At(0);

      if(jet->PT > 200)
      {

		  passesOneJetCut++;
		  Jet *jet2 = (Jet*) branchJet->At(1);
		  if(jet2->PT > 200)
		  {
			  histMet->Fill(metval);
			histJetPT->Fill(jet->PT);
			passesAllJetCuts++;
			}
	  }

      // Plot jet transverse momentum


      // Print jet transverse momentum
      //cout << jet->PT << endl;
    }
    //if(passesAllJetCuts>10)
    	//break;

  }
  printf("%lld\n%d\n%d\n",numberOfEntries,passesOneJetCut,passesAllJetCuts);

  // Show resulting histograms
  histJetPT->Draw();
  c1->SaveAs("PTs4.eps");
  gPad->SetLogy();
  histMet->Draw();
  c1->SaveAs("Met.eps");
}


int  main () {
   exampleMacro();
}
