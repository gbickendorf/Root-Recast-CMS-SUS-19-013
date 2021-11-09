#include <TROOT.h>
#include <iostream>
#include <fstream>
#include <TCanvas.h>
#include <TChain.h>
#include <TH1.h>
#include <TVector3.h>
#include "TClonesArray.h"
#include "TSystemDirectory.h"
#include "classes/DelphesClasses.h"
//#include "fastjet/ClusterSequence.hh"
#include "ExRootAnalysis/ExRootClasses.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"
using namespace std;
void exampleMacro()
{
 // vector<fastjet::PseudoJet> input_particles;
  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->cd();
  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add("/media/gerrit/Files/DelphesEvents/DelphesTestRun/mlm2.root");

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  // Get pointers to branches used in this analysis

  TClonesArray *branch_track = treeReader->UseBranch("Track");
  // Book histograms
  TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 100, 0.0, 1000.0);
  const double bins[7] = {300.0, 450.0, 600.0, 800.0, 1000.0, 1200.0, 2000.0};
  TH1 *histMet = new TH1F("ptmiss", "ptmiss", 6, bins);
  //gPad->SetLogy();

  // Loop over all events
  int passesOneJetCut = 0;
  int passesAllJetCuts = 0;
 // int passesWeight = 0;
 //numberOfEntries=1000;
  for (Int_t entry = 0; entry < numberOfEntries; ++entry)
  {

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    if(branch_track->GetEntries()==1)
      passesAllJetCuts++;
      printf("%lld\n",branch_track->GetEntries());

  }
  printf("%lld\n%d\n%d\n", numberOfEntries, passesOneJetCut, passesAllJetCuts);

  // Show resulting histograms
  histJetPT->Draw();
  c1->SaveAs("PTs4.eps");
  gPad->SetLogy();
  histMet->Draw();
  c1->SaveAs("Met.eps");
}

struct TestPlots
{
  TH1 *fElectronDeltaPT;
  TH1 *fElectronDeltaEta;

  TH1 *fPhotonDeltaPT;
  TH1 *fPhotonDeltaEta;
  TH1 *fPhotonDeltaE;

  TH1 *fMuonDeltaPT;
  TH1 *fMuonDeltaEta;

  TH1 *fJetDeltaPT;
};

double sq(double x)
{return x*x;}

void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots)
{
  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->cd();
  TH1 *histHT = new TH1F("HT", "HT", 100, 0.0, 1000.0);
  const double bins[7] = {300.0, 450.0, 600.0, 800.0, 1000.0, 1200.0, 2000.0};
  TH1 *histMet = new TH1F("ptmiss", "ptmiss", 6, bins);
  //TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  //TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  //TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  //TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray *branchAK4Jet = treeReader->UseBranch("AK4Jets");
  TClonesArray *branchAK8Jet = treeReader->UseBranch("AK8Jets");
  TClonesArray *branchMET = treeReader->UseBranch("MissingET");
  //TClonesArray *branchWeight = treeReader->UseBranch("LHEFEvent");
  //TClonesArray *branchHT = treeReader->UseBranch("ScalarHT");

  Long64_t allEntries = treeReader->GetEntries();
  //allEntries=10000;
  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  Electron *electron;
  Photon *photon;
  Muon *muon;

  Track *track;
  Tower *tower;

  Jet *jet;
  TObject *object;

  TLorentzVector momentum;

  Float_t Eem, Ehad;
  Bool_t skip;

  Long64_t entry;
  double MET, HT, angSepJ_HT, softDropMass;
  Int_t i, j, pdgCode, NJetAK4,Ntracks,ToggleAllCuts,cutLoop;
  Int_t NCut[10]={0,0,0,0,0,0,0,0,0,0};
  Int_t ToggleCut[10]={0,0,0,0,0,0,0,0,0,0};
  ToggleCut[0]=1; //Filter Zero Weight Events
  ToggleCut[1]=1; //NJetAK4 >= 2
  ToggleCut[2]=1; //PT_miss > 300GeV
  ToggleCut[3]=1; //HT > 400 GeV
  ToggleCut[4]=1; //Phi(j,HTMiss)>0.5(0.3)
  ToggleCut[5]=1; //~isolated Photon, Electron, Muon PT > 10 GeV
  ToggleCut[6]=1; //isolated Tracks mt> 100GeV, pt > 10GeV
  ToggleCut[7]=1; //2 AK8 Jets PT > 200 GeV
  ToggleCut[8]=1; //mJet of 2 AK8 Jets in [10,140]GeV
  ToggleCut[9]=1; //AngularSep of AK8 and bTagged
  ToggleAllCuts=0;
  TVector3 HTmiss;
 // allEntries=1000;

  Long64_t survived=0;
  // Loop over all events
  //allEntries=1000;
  for(entry = 0; entry < allEntries; ++entry)
  {
    if(entry % 1000 == 0)
      printf("%lld\n",entry);
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    //ScalarHT *scHT= (ScalarHT*)branchHT->At(0);
    //histHT->Fill(scHT->HT);
    //Weight *weight = (Weight*)branchWeight->At(0);
    //if(weight->Weight == 0&&ToggleCut[0])
    //{
    //  NCut[0]++;
    //  if(!ToggleAllCuts)
    //    continue;
    //}
    Ntracks= branchTrack->GetEntries();
    NJetAK4=branchAK4Jet->GetEntriesFast();
    if(NJetAK4<2 && ToggleCut[1])
    {
      NCut[1]++;
      if(!ToggleAllCuts)
        continue;
    }
    MET=((MissingET *)branchMET->At(0))->MET;
    if(MET<300.0 && ToggleCut[2])
    {
      NCut[2]++;
      if(!ToggleAllCuts)
        continue;    
    }
    // Loop over all electrons in event
    for(i = 0; i < branchElectron->GetEntriesFast(); ++i)
    {
      electron = (Electron*) branchElectron->At(i);
      particle = (GenParticle*) electron->Particle.GetObject();
      if(particle->PT > 10 && ToggleCut[5]&&electron->IsolationVar>0.1)
      {
        NCut[5]++;
        if(!ToggleAllCuts)
          continue;
      }
      //cout << "    Electrons pt: " << particle->PT << ", eta: " << particle->Eta << ", phi: " << particle->Phi << endl;
      //plots->fElectronDeltaPT->Fill((particle->PT - electron->PT)/particle->PT);
      //plots->fElectronDeltaEta->Fill((particle->Eta - electron->Eta)/particle->Eta);
    }

    // Loop over all photons in event
    for(i = 0; i < branchPhoton->GetEntriesFast(); ++i)
    {
      photon = (Photon*) branchPhoton->At(i);

      // skip photons with references to multiple particles
      if(photon->Particles.GetEntriesFast() != 1) continue;

      //particle = (GenParticle*) photon->Particles.At(0);
      if(photon->PT > 10 && ToggleCut[5])
      {
        NCut[5]++;
        if(!ToggleAllCuts)
          continue;
      }
      //cout << "    Photons pt: " << particle->PT << ", eta: " << particle->Eta << ", phi: " << particle->Phi << endl;
      //plots->fPhotonDeltaPT->Fill((particle->PT - photon->PT)/particle->PT);
      //plots->fPhotonDeltaEta->Fill((particle->Eta - photon->Eta)/particle->Eta);
      //plots->fPhotonDeltaE->Fill((particle->E - photon->E)/particle->E);
          
    }

    // Loop over all muons in event
    for(i = 0; i < branchMuon->GetEntriesFast(); ++i)
    {
      muon = (Muon*) branchMuon->At(i);
      particle = (GenParticle*) muon->Particle.GetObject();

      //plots->fMuonDeltaPT->Fill((particle->PT - muon->PT)/particle->PT);
      //plots->fMuonDeltaEta->Fill((particle->Eta - muon->Eta)/particle->Eta);
      if(particle->PT > 10 && ToggleCut[5]&&muon->IsolationVar>0.2)
      {
        NCut[5]++;
        if(!ToggleAllCuts)
          continue;
      }
      //cout << "    Muons pt: " << particle->PT << ", eta: " << particle->Eta << ", phi: " << particle->Phi << endl;
    }

    // cout << "--  New event -- " << endl;

    // Loop over all jets in event
    if(ToggleCut[6])
    {
      cutLoop=0;
      for (size_t itrack = 0; itrack < Ntracks; itrack++)
      {
        Track *track = (Track*)branchTrack->At(itrack);
        double teta=track->Eta;
        double tphi=track->Phi;
        if(sqrt(sq(track->Mass)+sq(track->PT))>100 && track->Eta >2.4)
          continue;
        
        double conePT=0.0;
        for (size_t icandtrack = 0; icandtrack < Ntracks; icandtrack++)
        {
          if(itrack==icandtrack)
            continue;
          Track *candTrack= (Track*)branchTrack->At(icandtrack);
          if(sqrt(sq(teta-candTrack->Eta)+sq(tphi-candTrack->Phi))<0.3)
            conePT+=candTrack->PT;
        }
        if(abs(track->PID)<20) 
        {
          if(track->PT>5&&conePT/track->PT<0.2)
            cutLoop=1;
        }
        else
        {
          if(track->PT>10&&conePT/track->PT<0.1)
            cutLoop=1;
        }
      }
      if(cutLoop)
      {
        NCut[6]++;
        continue;
      }
      
    }

    if(branchAK8Jet->GetEntriesFast() < 2)
    {
      //NCut[7]++;
      continue;
    }
    softDropMass=((Jet*) branchAK8Jet->At(0))->SoftDroppedJet.Mag();
    if((softDropMass<40.0 || softDropMass >140.0) && ToggleCut[8])
    {
      NCut[8]++;
      if(!ToggleAllCuts)
        continue;
    }
    softDropMass=((Jet*) branchAK8Jet->At(1))->SoftDroppedJet.Mag();
    if((softDropMass<40.0 || softDropMass >140.0) && ToggleCut[8])
    {
      NCut[8]++;
      if(!ToggleAllCuts)
        continue;
    }
    HT=0.0;
    HTmiss=TVector3(0.0,0.0,0.0);
    cutLoop=0;
    for(i = 0; i < NJetAK4; ++i)
    {
      jet = (Jet*) branchAK4Jet->At(i);
      HT+=jet->PT;
      HTmiss+=TVector3(jet->PT*cos(jet->Phi),jet->PT*sin(jet->Phi),jet->PT*sinh(jet->Eta));

      if(i==1&&jet->PT<200.0 && ToggleCut[7])
      {
        NCut[7]++;
        if(!ToggleAllCuts)
          continue;
      }

      // cout<<"Looping over jet constituents. Jet pt: "<<jet->PT<<", eta: "<<jet->Eta<<", phi: "<<jet->Phi<<endl;

      // Loop over all jet's constituents
      for(j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
      {
        object = jet->Constituents.At(j);

        // Check if the constituent is accessible
        if(object == 0) continue;

        if(object->IsA() == GenParticle::Class())
        {
          particle = (GenParticle*) object;
          cout << "    GenPart pt: " << particle->PT << ", eta: " << particle->Eta << ", phi: " << particle->Phi << endl;
          momentum += particle->P4();
        }
        else if(object->IsA() == Track::Class())
        {
          track = (Track*) object;
          //cout << "    Track pt: " << track->PT << ", eta: " << track->Eta << ", phi: " << track->Phi << endl;
          momentum += track->P4();
        }
        else if(object->IsA() == Tower::Class())
        {
          tower = (Tower*) object;
          //cout << "    Tower pt: " << tower->ET << ", eta: " << tower->Eta << ", phi: " << tower->Phi << endl;
          momentum += tower->P4();
        }
      }
      if(jet->BTag&&sqrt(sq(((Jet*) branchAK8Jet->At(1))->Eta-jet->Eta)+sq(((Jet*) branchAK8Jet->At(1))->Phi-jet->Phi))<0.8)
      {
        cutLoop++;
      }
      //plots->fJetDeltaPT->Fill((jet->PT - momentum.Pt())/jet->PT);
    }
    if(cutLoop&&ToggleCut[9])
    {
      NCut[9]++;
      if(!ToggleAllCuts)
        continue;
    }

    if(NJetAK4>1&&ToggleCut[4])
    {
      cutLoop=0;
      if(HTmiss.Phi()-((Jet*)branchAK4Jet->At(0))->Phi<0.5)
        cutLoop++;
      if(HTmiss.Phi()-((Jet*)branchAK4Jet->At(1))->Phi<0.5)
        cutLoop++;
      if(NJetAK4>2 && HTmiss.Phi()-((Jet*)branchAK4Jet->At(2))->Phi<0.3)
        cutLoop++;
      if(NJetAK4>3 && HTmiss.Phi()-((Jet*)branchAK4Jet->At(3))->Phi<0.3)
        cutLoop++;
      if(cutLoop)
      {
        //printf("%d\n",cutLoop);
        NCut[4]++;
        if(!ToggleAllCuts)
          continue;
      }
    }
    histHT->Fill(HT);
    if(HT<400&&ToggleCut[3])
    {
      NCut[3]++;
      if(!ToggleAllCuts)
        continue;
    }
    survived++;
    histMet->Fill(MET);
  }
  printf("\n\nTotal %lld\nSurvived %lld\nEfficiency %f\n",allEntries,survived,((float)survived)/allEntries);
  histHT->Draw();
  c1->SaveAs("HT.eps");
  histMet->Draw();
  c1->SaveAs("MET.eps");
  for(int i = 0; i < 10; i++)
  {
    printf("%d %d %d\n",i,ToggleCut[i],NCut[i]);
  }
}

int main()
{  TChain chains("Delphes");
  chains.Add("/media/gerrit/Files/DelphesEvents/DelphesTestRun/Events.root");
  ExRootTreeReader *treeReaders = new ExRootTreeReader(&chains);
  TestPlots *plotss = new TestPlots;
  AnalyseEvents(treeReaders,plotss);
  //exampleMacro();
  return 0;

    // Create chain of root trees
  TChain chain("Delphes");
  std::ifstream file("ROOTFILES.txt");
  int i = 3;
  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      chain.Add(line.c_str());
      i--;
      if(i == 0)
        break;
    }
    file.close();
  }
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  TestPlots *plots = new TestPlots;
  AnalyseEvents(treeReader,plots);
  //exampleMacro(); /cephfs/user/s6gebick/ClusterRes/Z/DelphesEvents/run_231/Events.root
}
