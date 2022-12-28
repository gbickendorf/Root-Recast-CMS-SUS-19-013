#include "Analysis.h"

double Analysis::sq(double x)
{
  return x * x;
}

double Analysis::AnalyseEventsNew(ExRootTreeReader *treeReader, vector<PassedEvent> &events, int MinNumEvents, int status)
{
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchAK4Jet = treeReader->UseBranch("AK4Jets");
  TClonesArray *branchAK8Jet = treeReader->UseBranch("AK8Jets");
  TClonesArray *branchMET = treeReader->UseBranch("MissingET");

  TClonesArray *branchElectrons[7];
  TClonesArray *branchMuons[7];
  double minPT[] = {13000.0, 200.0, 133.0, 100.0, 80.0, 66.0, 57.0, 10.0};
  for (size_t i = 0; i < 7; i++)
  {
    string bname = "electrons" + to_string(i);
    branchElectrons[i] = treeReader->UseBranch(("electrons" + to_string(i + 1)).c_str());
    branchMuons[i] = treeReader->UseBranch(("muons" + to_string(i + 1)).c_str());
  }
  Long64_t allEntries = treeReader->GetEntries();
  if (MinNumEvents > 0)
  {
    if (MinNumEvents > allEntries)
    {
      cerr << "Not enough Events!" << endl;
      return 0.0;
    }
    allEntries = MinNumEvents;
  }
  cout << "** Chain contains " << allEntries << " events" << endl;

  Electron *electron;
  Photon *photon;
  Muon *muon;

  Track *track;
  Tower *tower;

  Jet *jet;

  TLorentzVector momentum;

  Long64_t entry;
  double MET, HT, softDropMass;
  Int_t Ntracks, i, j;
  Int_t ToggleAllCuts, cutLoop;
  Int_t NCut[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Int_t ToggleCut[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  ToggleCut[0] = 1; // Filter Zero Weight Events
  ToggleCut[1] = 1; // NJetAK4 >= 2
  ToggleCut[2] = 1; // PT_miss > 300GeV
  ToggleCut[3] = 1; // HT > 400 GeV
  ToggleCut[4] = 1; // Phi(j,HTMiss)>0.5(0.3)
  ToggleCut[5] = 1; //~isolated Photon, Electron, Muon PT > 10 GeV
  ToggleCut[6] = 1; // isolated Tracks mt> 100GeV, pt > 10GeV
  ToggleCut[7] = 1; // 2 AK8 Jets PT > 200 GeV
  ToggleCut[8] = 1; // mJet of 2 AK8 Jets in [10,140]GeV
  ToggleCut[9] = 1; // AngularSep of AK8 and bTagged
  ToggleAllCuts = 0;
  TVector3 HTmiss;
  // allEntries=1000;
  vector<Jet *> AK4Jets;
  vector<Jet *> AK8Jets;

  Long64_t survived = 0;
  boost::progress_display *show_progress = new boost::progress_display(allEntries);
  for (entry = 0; entry < allEntries; ++entry)
  {
    AK4Jets.clear();
    AK8Jets.clear();
    PassedEvent event;
    int nLep = 0;
    int nPho = 0;
    double PTPho = 0.0;
    event.status = status;
    ++(*show_progress);
    vector<int> JetAK4Indices;
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    MET = ((MissingET *)branchMET->At(0))->MET;
    Ntracks = branchTrack->GetEntries();
    for (size_t i = 0; i < branchAK4Jet->GetEntriesFast(); i++)
    {
      jet = (Jet *)branchAK4Jet->At(i);
      if (abs(jet->Eta) < 2.4 && jet->PT > 30)
        AK4Jets.push_back(jet);
    }
    for (size_t i = 0; i < branchAK8Jet->GetEntriesFast(); i++)
    {
      jet = (Jet *)branchAK8Jet->At(i);
      if (abs(jet->Eta) < 2.4 && jet->PT > 30)
        AK8Jets.push_back(jet);
    }

    if (MET == 0.0)
    {
      NCut[0]++;
      continue;
    }

    if (AK4Jets.size() < 2 && ToggleCut[1])
    {
      NCut[1]++;
      if (!ToggleAllCuts)
        continue;
    }

    if (MET < 300 && ToggleCut[2])
    {
      NCut[2]++;
      continue;
    }

    HT = 0.0;
    HTmiss = TVector3(0.0, 0.0, 0.0);
    cutLoop = 0;
    int cutBTag = 0;
    for (i = 0; i < AK4Jets.size(); ++i)
    {
      HT += AK4Jets[i]->PT;
    }
    if (HT < 400 && ToggleCut[3])
    {
      NCut[3]++;
      continue;
    }

    for (i = 0; i < AK4Jets.size(); ++i)
    {
      jet = AK4Jets[i];
      HTmiss -= TVector3(jet->PT * cos(jet->Phi), jet->PT * sin(jet->Phi), jet->PT * sinh(jet->Eta));
    }

    cutLoop = 0;
    if (abs(HTmiss.Phi() - AK4Jets[0]->Phi) < 0.5)
    {
      NCut[4]++;
      continue;
    }
    if (abs(HTmiss.Phi() - AK4Jets[1]->Phi) < 0.5)
    {
      NCut[4]++;
      continue;
    }
    if (AK4Jets.size() > 2 && abs(HTmiss.Phi() - AK4Jets[2]->Phi) < 0.3)
    {
      NCut[4]++;
      continue;
    }
    if (AK4Jets.size() > 3 && abs(HTmiss.Phi() - AK4Jets[3]->Phi) < 0.3)
    {
      NCut[4]++;
      continue;
    }

    for (i = 0; i < 7; i++)
    {
      // cout<<branchElectrons[i]->GetEntriesFast()<<"  "<<i<<endl;
      for (size_t lepIndex = 0; lepIndex < branchElectrons[i]->GetEntriesFast(); lepIndex++)
      {
        electron = (Electron *)branchElectrons[i]->At(lepIndex);
        if (minPT[i] > electron->PT && minPT[i + 1] < electron->PT && electron->IsolationVar < 0.1)
          nLep++;
      }
      for (size_t lepIndex = 0; lepIndex < branchMuons[i]->GetEntriesFast(); lepIndex++)
      {
        muon = (Muon *)branchMuons[i]->At(lepIndex);
        if (minPT[i] > muon->PT && minPT[i + 1] < muon->PT && muon->IsolationVar < 0.2)
          nLep++;
      }
    }

    for (i = 0; i < branchPhoton->GetEntriesFast(); ++i)
    {
      photon = (Photon *)branchPhoton->At(i);
      // if (photon->Particles.GetEntriesFast() != 1)
      //  continue;
      double looseWP = 1.3 / photon->PT + 0.005;
      if (photon->PT > 10 && ToggleCut[5] && photon->IsolationVar < 1.3 / photon->PT + 0.005)
      {
        nPho++;
        PTPho = photon->PT;
      }
    }
    if (nLep + nPho > 0 && ToggleCut[5])
    {
      NCut[5]++;
      continue;
    }

    cutLoop = 0;
    for (int itrack = 0; itrack < Ntracks; itrack++)
    {
      Track *track = (Track *)branchTrack->At(itrack);
      MissingET *missingET = ((MissingET *)branchMET->At(0));
      double teta = track->Eta;
      double tphi = track->Phi;
      double mT = sqrt(2.0 * track->PT * missingET->MET * (1.0 - cos(missingET->Phi - track->Phi)));
      if (mT > 100 || abs(track->Eta) > 2.4)
        continue;

      double conePT = 0.0;
      for (int icandtrack = 0; icandtrack < Ntracks; icandtrack++)
      {
        if (itrack == icandtrack)
          continue;
        Track *candTrack = (Track *)branchTrack->At(icandtrack);
        if (sqrt(sq(teta - candTrack->Eta) + sq(tphi - candTrack->Phi)) < 0.3)
          conePT += candTrack->PT;
      }
      if (abs(track->PID) < 20)
      {
        if (track->PT > 5 && conePT / track->PT < 0.2)
          cutLoop = 1;
      }
      else
      {
        if (track->PT > 10 && conePT / track->PT < 0.1)
          cutLoop = 1;
      }
    }
    if (cutLoop && ToggleCut[6])
    {
      NCut[6]++;
      continue;
    }

    if ((AK8Jets.size() < 2 || AK8Jets[1]->PT < 200) && ToggleCut[7])
    {
      NCut[7]++;
      continue;
    }

    softDropMass = AK8Jets[0]->SoftDroppedJet.Mag();
    event.mj1 = softDropMass;
    if ((softDropMass < 40.0 || softDropMass > 140.0) && ToggleCut[8])
    {
      NCut[8]++;
      if (!ToggleAllCuts)
        continue;
    }
    softDropMass = AK8Jets[1]->SoftDroppedJet.Mag();
    event.mj2 = softDropMass;
    if ((softDropMass < 40.0 || softDropMass > 140.0) && ToggleCut[8])
    {
      NCut[8]++;
      if (!ToggleAllCuts)
        continue;
    }
    cutLoop = 0;
    for (i = 0; i < AK4Jets.size(); ++i)
    {
      jet = AK4Jets[i];
      if (jet->BTag && sqrt(sq(AK8Jets[1]->Eta - jet->Eta) + sq(AK8Jets[1]->Phi - jet->Phi)) < 0.8)
      {
        cutBTag++;
      }
    }
    for (i = 0; i < AK8Jets.size(); ++i)
    {
      if (i == 1)
        continue;
      jet = AK8Jets[i];
      if (jet->BTag && sqrt(sq(AK8Jets[1]->Eta - jet->Eta) + sq(AK8Jets[1]->Phi - jet->Phi)) < 0.8)
      {
        cutBTag++;
      }
    }
    if (cutBTag && ToggleCut[9])
    {
      NCut[9]++;
      continue;
    }

    if (MET > 200.0 && nLep == 1 && cutBTag == 0 && nPho == 0)
    {
      event.ptmiss = MET;
      event.status = 10;
      events.push_back(event);
      NCut[5]++;
      continue;
    }
    if (nLep == 0 && nPho == 1 && PTPho > 200)
    {
      event.ptmiss = PTPho;
      event.status = 11;
      events.push_back(event);
      NCut[5]++;
      if (cutBTag > 0 && ToggleCut[9])
        NCut[9]++;
      continue;
    }
    if (MET < 300 && ToggleCut[2])
    {
      continue;
    }
    if (nLep + nPho > 0 && ToggleCut[5])
    {
      NCut[5]++;
      continue;
    }
    if (cutBTag > 0 && ToggleCut[9])
    {
      NCut[9]++;
      continue;
    }

    event.ptmiss = MET;
    survived++;
    events.push_back(event);
  }
  printf("\n\nTotal %lld\nSurvived %lld\nEfficiency %f\n", allEntries, survived, ((float)survived) / allEntries);
  for (int i = 0; i < 10; i++)
  {
    printf("%d %d %d\n", i, ToggleCut[i], NCut[i]);
  }
  return ((double)survived) / allEntries;
}

double Analysis::AnalyseEvents(ExRootTreeReader *treeReader, vector<PassedEvent> &events, int MinNumEvents, int status)
{
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchAK4Jet = treeReader->UseBranch("AK4Jets");
  TClonesArray *branchAK8Jet = treeReader->UseBranch("AK8Jets");
  TClonesArray *branchMET = treeReader->UseBranch("MissingET");

  TClonesArray *branchElectrons[7];
  TClonesArray *branchMuons[7];
  double minPT[] = {13000.0, 200.0, 133.0, 100.0, 80.0, 66.0, 57.0, 10.0};
  for (size_t i = 0; i < 7; i++)
  {
    string bname = "electrons" + to_string(i);
    branchElectrons[i] = treeReader->UseBranch(("electrons" + to_string(i + 1)).c_str());
    branchMuons[i] = treeReader->UseBranch(("muons" + to_string(i + 1)).c_str());
  }
  Long64_t allEntries = treeReader->GetEntries();
  if (MinNumEvents > 0)
  {
    if (MinNumEvents > allEntries)
    {
      cerr << "Not enough Events!" << endl;
      return 0.0;
    }
    allEntries = MinNumEvents;
  }
  cout << "** Chain contains " << allEntries << " events" << endl;

  Electron *electron;
  Photon *photon;
  Muon *muon;

  Track *track;
  Tower *tower;

  Jet *jet;

  TLorentzVector momentum;

  Long64_t entry;
  double MET, HT, softDropMass;
  Int_t Ntracks, i, j;
  Int_t ToggleAllCuts, cutLoop;
  Int_t NCut[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Int_t ToggleCut[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  ToggleCut[0] = 1; // Filter Zero Weight Events
  ToggleCut[1] = 1; // NJetAK4 >= 2
  ToggleCut[2] = 1; // PT_miss > 300GeV
  ToggleCut[3] = 1; // HT > 400 GeV
  ToggleCut[4] = 1; // Phi(j,HTMiss)>0.5(0.3)
  ToggleCut[5] = 1; //~isolated Photon, Electron, Muon PT > 10 GeV
  ToggleCut[6] = 1; // isolated Tracks mt> 100GeV, pt > 10GeV
  ToggleCut[7] = 1; // 2 AK8 Jets PT > 200 GeV
  ToggleCut[8] = 1; // mJet of 2 AK8 Jets in [10,140]GeV
  ToggleCut[9] = 1; // AngularSep of AK8 and bTagged
  ToggleAllCuts = 0;
  TVector3 HTmiss;
  // allEntries=1000;
  vector<Jet *> AK4Jets;
  vector<Jet *> AK8Jets;

  Long64_t survived = 0;
  boost::progress_display *show_progress = new boost::progress_display(allEntries);
  for (entry = 0; entry < allEntries; ++entry)
  {
    AK4Jets.clear();
    AK8Jets.clear();
    PassedEvent event;
    int nLep = 0;
    int nPho = 0;
    double PTPho = 0.0;
    event.status = status;
    ++(*show_progress);
    vector<int> JetAK4Indices;
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    MET = ((MissingET *)branchMET->At(0))->MET;
    Ntracks = branchTrack->GetEntries();
    for (size_t i = 0; i < branchAK4Jet->GetEntriesFast(); i++)
    {
      jet = (Jet *)branchAK4Jet->At(i);
      if (abs(jet->Eta) < 2.4 && jet->PT > 30)
        AK4Jets.push_back(jet);
    }
    for (size_t i = 0; i < branchAK8Jet->GetEntriesFast(); i++)
    {
      jet = (Jet *)branchAK8Jet->At(i);
      if (abs(jet->Eta) < 2.4 && jet->PT > 30)
        AK8Jets.push_back(jet);
    }

    if (MET == 0.0)
    {
      NCut[0]++;
      continue;
    }

    if (AK4Jets.size() < 2 && ToggleCut[1])
    {
      NCut[1]++;
      if (!ToggleAllCuts)
        continue;
    }

    if (MET < 300.0 && ToggleCut[2])
    {
      NCut[2]++;
      continue;
    }

    HT = 0.0;
    HTmiss = TVector3(0.0, 0.0, 0.0);
    cutLoop = 0;
    int cutBTag = 0;
    for (i = 0; i < AK4Jets.size(); ++i)
    {
      HT += AK4Jets[i]->PT;
    }
    if (HT < 400 && ToggleCut[3])
    {
      NCut[3]++;
      continue;
    }

    for (i = 0; i < AK4Jets.size(); ++i)
    {
      jet = AK4Jets[i];
      HTmiss -= TVector3(jet->PT * cos(jet->Phi), jet->PT * sin(jet->Phi), jet->PT * sinh(jet->Eta));
    }

    cutLoop = 0;
    if (abs(HTmiss.Phi() - AK4Jets[0]->Phi) < 0.5)
    {
      NCut[4]++;
      continue;
    }
    if (abs(HTmiss.Phi() - AK4Jets[1]->Phi) < 0.5)
    {
      NCut[4]++;
      continue;
    }
    if (AK4Jets.size() > 2 && abs(HTmiss.Phi() - AK4Jets[2]->Phi) < 0.3)
    {
      NCut[4]++;
      continue;
    }
    if (AK4Jets.size() > 3 && abs(HTmiss.Phi() - AK4Jets[3]->Phi) < 0.3)
    {
      NCut[4]++;
      continue;
    }

    for (i = 0; i < 7; i++)
    {
      // cout<<branchElectrons[i]->GetEntriesFast()<<"  "<<i<<endl;
      for (size_t lepIndex = 0; lepIndex < branchElectrons[i]->GetEntriesFast(); lepIndex++)
      {
        electron = (Electron *)branchElectrons[i]->At(lepIndex);
        if (minPT[i] > electron->PT && minPT[i + 1] < electron->PT && electron->IsolationVar < 0.1)
          nLep++;
      }
      for (size_t lepIndex = 0; lepIndex < branchMuons[i]->GetEntriesFast(); lepIndex++)
      {
        muon = (Muon *)branchMuons[i]->At(lepIndex);
        if (minPT[i] > muon->PT && minPT[i + 1] < muon->PT && muon->IsolationVar < 0.2)
          nLep++;
      }
    }

    for (i = 0; i < branchPhoton->GetEntriesFast(); ++i)
    {
      photon = (Photon *)branchPhoton->At(i);
      // if (photon->Particles.GetEntriesFast() != 1)
      //  continue;
      double looseWP = 1.3 / photon->PT + 0.005;
      if (photon->PT > 10 && ToggleCut[5] && photon->IsolationVar < looseWP)
      {
        nPho++;
        PTPho = photon->PT;
      }
    }
    if (nLep + nPho > 0 && ToggleCut[5])
    {
      NCut[5]++;
      continue;
    }

    cutLoop = 0;
    for (int itrack = 0; itrack < Ntracks; itrack++)
    {
      Track *track = (Track *)branchTrack->At(itrack);
      MissingET *missingET = ((MissingET *)branchMET->At(0));
      double teta = track->Eta;
      double tphi = track->Phi;
      double mT = sqrt(2.0 * track->PT * missingET->MET * (1.0 - cos(missingET->Phi - track->Phi)));
      if (mT > 100 || abs(track->Eta) > 2.4)
        continue;

      double conePT = 0.0;
      for (int icandtrack = 0; icandtrack < Ntracks; icandtrack++)
      {
        if (itrack == icandtrack)
          continue;
        Track *candTrack = (Track *)branchTrack->At(icandtrack);
        if (sqrt(sq(teta - candTrack->Eta) + sq(tphi - candTrack->Phi)) < 0.3)
          conePT += candTrack->PT;
      }
      if (abs(track->PID) < 20)
      {
        if (track->PT > 5 && conePT / track->PT < 0.2)
          cutLoop = 1;
      }
      else
      {
        if (track->PT > 10 && conePT / track->PT < 0.1)
          cutLoop = 1;
      }
    }
    if (cutLoop && ToggleCut[6])
    {
      NCut[6]++;
      continue;
    }

    if ((AK8Jets.size() < 2 || AK8Jets[1]->PT < 200) && ToggleCut[7])
    {
      NCut[7]++;
      continue;
    }

    softDropMass = AK8Jets[0]->SoftDroppedJet.Mag();
    event.mj1 = softDropMass;
    if ((softDropMass < 40.0 || softDropMass > 140.0) && ToggleCut[8])
    {
      NCut[8]++;
      continue;
    }
    softDropMass = AK8Jets[1]->SoftDroppedJet.Mag();
    event.mj2 = softDropMass;
    if ((softDropMass < 40.0 || softDropMass > 140.0) && ToggleCut[8])
    {
      NCut[8]++;
      continue;
    }
    // cutLoop = 0;
    for (i = 0; i < AK4Jets.size(); ++i)
    {
      jet = AK4Jets[i];
      if (jet->BTag && sqrt(sq(AK8Jets[1]->Eta - jet->Eta) + sq(AK8Jets[1]->Phi - jet->Phi)) < 0.8)
      {
        cutBTag++;
      }
      if (jet->BTag && sqrt(sq(AK8Jets[0]->Eta - jet->Eta) + sq(AK8Jets[0]->Phi - jet->Phi)) < 0.8)
      {
        cutBTag++;
      }
    }
    if (cutBTag && ToggleCut[9])
    {
      NCut[9]++;
      continue;
    }
    event.ptmiss = MET;
    survived++;
    events.push_back(event);
  }
  printf("\n\nTotal %lld\nSurvived %lld\nEfficiency %f\n", allEntries, survived, ((float)survived) / allEntries);
  for (int i = 0; i < 10; i++)
  {
    printf("%d %d %d\n", i, ToggleCut[i], NCut[i]);
  }
  return ((double)survived) / allEntries;
}

void Analysis::AnalyseEventsSinglePhotonSample(ExRootTreeReader *treeReader, vector<PassedEvent> &events, int MinNumEvents)
{
  int status = 11;
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchAK4Jet = treeReader->UseBranch("AK4Jets");
  TClonesArray *branchAK8Jet = treeReader->UseBranch("AK8Jets");

  Long64_t allEntries = treeReader->GetEntries();
  if (MinNumEvents > 0)
  {
    if (MinNumEvents > allEntries)
    {
      cerr << "Not enough Events!" << endl;
      return;
    }
    allEntries = MinNumEvents;
  }
  cout << "** Chain contains " << allEntries << " events" << endl;

  Electron *electron;
  Photon *photon;
  Muon *muon;

  Track *track;
  Tower *tower;

  Jet *jet;
  TObject *object;

  TLorentzVector momentum;

  Long64_t entry;
  double HT, softDropMass;
  Int_t NJetAK4, Ntracks, i, j;
  Int_t ToggleAllCuts, cutLoop;
  Int_t NCut[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Int_t ToggleCut[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  ToggleCut[0] = 1; // Filter Zero Weight Events
  ToggleCut[1] = 1; // NJetAK4 >= 2
  ToggleCut[2] = 1; // PT_photon > 200GeV
  ToggleCut[3] = 1; // HT > 400 GeV
  ToggleCut[4] = 1; // Phi(j,HTMiss)>0.5(0.3)
  ToggleCut[5] = 1; //~isolated Photon, Electron, Muon PT > 10 GeV
  ToggleCut[6] = 1; // isolated Tracks mt> 100GeV, pt > 10GeV
  ToggleCut[7] = 1; // 2 AK8 Jets PT > 200 GeV
  ToggleCut[8] = 1; // mJet of 2 AK8 Jets in [10,140]GeV
  ToggleCut[9] = 0; // AngularSep of AK8 and bTagged
  ToggleAllCuts = 0;
  TVector3 HTmiss;
  // allEntries=1000;

  Long64_t survived = 0;
  boost::progress_display *show_progress = new boost::progress_display(allEntries);
  for (entry = 0; entry < allEntries; ++entry)
  {
    int nPhoton = 0;
    PassedEvent event;
    event.status = status;
    ++(*show_progress);
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    Ntracks = branchTrack->GetEntries();
    NJetAK4 = branchAK4Jet->GetEntriesFast();
    if (NJetAK4 < 2 && ToggleCut[1])
    {
      NCut[1]++;
      if (!ToggleAllCuts)
        continue;
    }
    cutLoop = 0;
    // Loop over all electrons in event
    for (i = 0; i < branchElectron->GetEntriesFast(); ++i)
    {
      electron = (Electron *)branchElectron->At(i);

      if (electron->PT > 10 && ToggleCut[5] && electron->IsolationVar > 0.1)
      {
        NCut[5]++;
        cutLoop++;
        if (!ToggleAllCuts)
          break;
      }
    }

    // Loop over all photons in event
    for (i = 0; i < branchPhoton->GetEntriesFast(); ++i)
    {
      photon = (Photon *)branchPhoton->At(i);

      // skip photons with references to multiple particles
      if (photon->Particles.GetEntriesFast() != 1)
        continue;
      // particle = (GenParticle*) photon->Particles.At(0);
      if (photon->PT > 10 && ToggleCut[5])
      {
        if (nPhoton)
        {
          NCut[5]++;
          cutLoop++;
          if (!ToggleAllCuts)
            break;
        }
        event.ptmiss = (double)photon->PT;
        nPhoton++;
      }
    }
    if (nPhoton != 1)
      continue;
    if (ToggleCut[2] && event.ptmiss < 200)
    {
      NCut[2]++;
      continue;
    }

    // Loop over all muons in event
    for (i = 0; i < branchMuon->GetEntriesFast(); ++i)
    {
      muon = (Muon *)branchMuon->At(i);
      if (muon->PT > 10 && ToggleCut[5] && muon->IsolationVar > 0.2)
      {
        NCut[5]++;
        cutLoop++;
        if (!ToggleAllCuts)
          break;
      }
    }
    if (cutLoop)
      continue;
    // Loop over all jets in event
    if (ToggleCut[6])
    {
      cutLoop = 0;
      for (int itrack = 0; itrack < Ntracks; itrack++)
      {
        Track *track = (Track *)branchTrack->At(itrack);
        double teta = track->Eta;
        double tphi = track->Phi;
        if (sqrt(sq(track->Mass) + sq(track->PT)) > 100 && track->Eta > 2.4)
          continue;

        double conePT = 0.0;
        for (int icandtrack = 0; icandtrack < Ntracks; icandtrack++)
        {
          if (itrack == icandtrack)
            continue;
          Track *candTrack = (Track *)branchTrack->At(icandtrack);
          if (sqrt(sq(teta - candTrack->Eta) + sq(tphi - candTrack->Phi)) < 0.3)
            conePT += candTrack->PT;
        }
        if (abs(track->PID) < 20)
        {
          if (track->PT > 5 && conePT / track->PT < 0.2)
            cutLoop = 1;
        }
        else
        {
          if (track->PT > 10 && conePT / track->PT < 0.1)
            cutLoop = 1;
        }
      }
      if (cutLoop)
      {
        NCut[6]++;
        continue;
      }
    }

    if (branchAK8Jet->GetEntriesFast() < 2)
    {
      // NCut[7]++;
      continue;
    }
    softDropMass = ((Jet *)branchAK8Jet->At(0))->SoftDroppedJet.Mag();
    event.mj1 = softDropMass;
    if ((softDropMass < 40.0 || softDropMass > 140.0) && ToggleCut[8])
    {
      NCut[8]++;
      if (!ToggleAllCuts)
        continue;
    }
    softDropMass = ((Jet *)branchAK8Jet->At(1))->SoftDroppedJet.Mag();
    event.mj2 = softDropMass;
    if ((softDropMass < 40.0 || softDropMass > 140.0) && ToggleCut[8])
    {
      NCut[8]++;
      if (!ToggleAllCuts)
        continue;
    }
    HT = 0.0;
    HTmiss = TVector3(0.0, 0.0, 0.0);
    cutLoop = 0;
    for (i = 0; i < NJetAK4; ++i)
    {
      jet = (Jet *)branchAK4Jet->At(i);
      HT += jet->PT;
      HTmiss += TVector3(jet->PT * cos(jet->Phi), jet->PT * sin(jet->Phi), jet->PT * sinh(jet->Eta));

      if (i == 1 && jet->PT < 200.0 && ToggleCut[7])
      {
        NCut[7]++;
        if (!ToggleAllCuts)
          continue;
      }
      for (j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
      {
        object = jet->Constituents.At(j);

        if (object == 0)
          continue;

        if (object->IsA() == GenParticle::Class())
        {
          printf("FUCK\n");
        }
        else if (object->IsA() == Track::Class())
        {
          track = (Track *)object;
          momentum += track->P4();
        }
        else if (object->IsA() == Tower::Class())
        {
          tower = (Tower *)object;
          momentum += tower->P4();
        }
      }
      if (jet->BTag && sqrt(sq(((Jet *)branchAK8Jet->At(1))->Eta - jet->Eta) + sq(((Jet *)branchAK8Jet->At(1))->Phi - jet->Phi)) < 0.8)
      {
        cutLoop++;
      }
    }
    if (cutLoop && ToggleCut[9])
    {
      NCut[9]++;
      if (!ToggleAllCuts)
        continue;
    }

    if (NJetAK4 > 1 && ToggleCut[4])
    {
      cutLoop = 0;
      if (abs(HTmiss.Phi() - ((Jet *)branchAK4Jet->At(0))->Phi) < 0.5)
        cutLoop++;
      if (abs(HTmiss.Phi() - ((Jet *)branchAK4Jet->At(1))->Phi) < 0.5)
        cutLoop++;
      if (NJetAK4 > 2 && abs(HTmiss.Phi() - ((Jet *)branchAK4Jet->At(2))->Phi) < 0.3)
        cutLoop++;
      if (NJetAK4 > 3 && abs(HTmiss.Phi() - ((Jet *)branchAK4Jet->At(3))->Phi) < 0.3)
        cutLoop++;
      if (cutLoop)
      {
        NCut[4]++;
        if (!ToggleAllCuts)
          continue;
      }
    }
    if (HT < 400 && ToggleCut[3])
    {
      NCut[3]++;
      if (!ToggleAllCuts)
        continue;
    }
    survived++;
    events.push_back(event);
    // histMet->Fill(MET);
  }
  printf("\n\nTotal %lld\nSurvived %lld\nEfficiency %f\n", allEntries, survived, ((float)survived) / allEntries);
  for (int i = 0; i < 10; i++)
  {
    printf("%d %d %d\n", i, ToggleCut[i], NCut[i]);
  }
}

void Analysis::AnalyseEventsSingleLeptonSample(ExRootTreeReader *treeReader, vector<PassedEvent> &events, int MinNumEvents)
{
  int status = 10;
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchAK4Jet = treeReader->UseBranch("AK4Jets");
  TClonesArray *branchAK8Jet = treeReader->UseBranch("AK8Jets");
  TClonesArray *branchMET = treeReader->UseBranch("MissingET");

  Long64_t allEntries = treeReader->GetEntries();
  if (MinNumEvents > 0)
  {
    if (MinNumEvents > allEntries)
    {
      cerr << "Not enough Events!" << endl;
      return;
    }
    allEntries = MinNumEvents;
  }
  cout << "** Chain contains " << allEntries << " events" << endl;

  Electron *electron;
  Photon *photon;
  Muon *muon;

  Track *track;
  Tower *tower;

  Jet *jet;
  TObject *object;

  TLorentzVector momentum;

  Long64_t entry;
  double MET, HT, softDropMass;
  Int_t NJetAK4, Ntracks, i, j;
  Int_t ToggleAllCuts, cutLoop;
  Int_t NCut[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Int_t ToggleCut[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  ToggleCut[0] = 1; // Filter Zero Weight Events
  ToggleCut[1] = 1; // NJetAK4 >= 2
  ToggleCut[2] = 1; // PT_miss > 200GeV
  ToggleCut[3] = 1; // HT > 400 GeV
  ToggleCut[4] = 1; // Phi(j,HTMiss)>0.5(0.3)
  ToggleCut[5] = 1; //~isolated Photon, Electron, Muon PT > 10 GeV
  ToggleCut[6] = 1; // isolated Tracks mt> 100GeV, pt > 10GeV
  ToggleCut[7] = 1; // 2 AK8 Jets PT > 200 GeV
  ToggleCut[8] = 1; // mJet of 2 AK8 Jets in [10,140]GeV
  ToggleCut[9] = 1; // AngularSep of AK8 and bTagged
  ToggleAllCuts = 0;
  TVector3 HTmiss;
  // allEntries=1000;

  Long64_t survived = 0;
  boost::progress_display *show_progress = new boost::progress_display(allEntries);
  for (entry = 0; entry < allEntries; ++entry)
  {
    int nLepton = 0;
    PassedEvent event;
    event.status = status;
    ++(*show_progress);
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    Ntracks = branchTrack->GetEntries();
    NJetAK4 = branchAK4Jet->GetEntriesFast();
    if (NJetAK4 < 2 && ToggleCut[1])
    {
      NCut[1]++;
      if (!ToggleAllCuts)
        continue;
    }
    MET = ((MissingET *)branchMET->At(0))->MET;
    if (MET < 200.0 && ToggleCut[2])
    {
      NCut[2]++;
      if (!ToggleAllCuts)
        continue;
    }
    event.ptmiss = MET;
    // Loop over all electrons in event
    for (i = 0; i < branchElectron->GetEntriesFast(); ++i)
    {
      electron = (Electron *)branchElectron->At(i);

      if (electron->PT > 10 && ToggleCut[5] && electron->IsolationVar > 0.1)
      {
        nLepton++;
      }
    }

    // Loop over all photons in event
    cutLoop = 0;
    for (i = 0; i < branchPhoton->GetEntriesFast(); ++i)
    {
      photon = (Photon *)branchPhoton->At(i);

      // skip photons with references to multiple particles
      if (photon->Particles.GetEntriesFast() != 1)
        continue;

      // particle = (GenParticle*) photon->Particles.At(0);
      if (photon->PT > 10 && ToggleCut[5])
        cutLoop++;
    }
    if (cutLoop)
    {
      NCut[5]++;
      continue;
    }

    // Loop over all muons in event
    for (i = 0; i < branchMuon->GetEntriesFast(); ++i)
    {
      muon = (Muon *)branchMuon->At(i);
      if (muon->PT > 10 && ToggleCut[5] && muon->IsolationVar > 0.2)
      {
        nLepton++;
      }
    }
    if (nLepton != 1)
      continue;
    // Loop over all jets in event
    if (ToggleCut[6])
    {
      cutLoop = 0;
      for (int itrack = 0; itrack < Ntracks; itrack++)
      {
        Track *track = (Track *)branchTrack->At(itrack);
        double teta = track->Eta;
        double tphi = track->Phi;
        if (sqrt(sq(track->Mass) + sq(track->PT)) > 100 && abs(track->Eta) > 2.4)
          continue;

        double conePT = 0.0;
        for (int icandtrack = 0; icandtrack < Ntracks; icandtrack++)
        {
          if (itrack == icandtrack)
            continue;
          Track *candTrack = (Track *)branchTrack->At(icandtrack);
          if (sqrt(sq(teta - candTrack->Eta) + sq(tphi - candTrack->Phi)) < 0.3)
            conePT += candTrack->PT;
        }
        if (abs(track->PID) < 20)
        {
          if (track->PT > 5 && conePT / track->PT < 0.2)
            cutLoop = 1;
        }
        else
        {
          if (track->PT > 10 && conePT / track->PT < 0.1)
            cutLoop = 1;
        }
      }
      if (cutLoop)
      {
        NCut[6]++;
        continue;
      }
    }

    if (branchAK8Jet->GetEntriesFast() < 2)
    {
      // NCut[7]++;
      continue;
    }
    softDropMass = ((Jet *)branchAK8Jet->At(0))->SoftDroppedJet.Mag();
    event.mj1 = softDropMass;
    if ((softDropMass < 40.0 || softDropMass > 140.0) && ToggleCut[8])
    {
      NCut[8]++;
      if (!ToggleAllCuts)
        continue;
    }
    softDropMass = ((Jet *)branchAK8Jet->At(1))->SoftDroppedJet.Mag();
    event.mj2 = softDropMass;
    if ((softDropMass < 40.0 || softDropMass > 140.0) && ToggleCut[8])
    {
      NCut[8]++;
      if (!ToggleAllCuts)
        continue;
    }
    HT = 0.0;
    HTmiss = TVector3(0.0, 0.0, 0.0);
    cutLoop = 0;
    for (i = 0; i < NJetAK4; ++i)
    {
      jet = (Jet *)branchAK4Jet->At(i);
      HT += jet->PT;
      HTmiss += TVector3(jet->PT * cos(jet->Phi), jet->PT * sin(jet->Phi), jet->PT * sinh(jet->Eta));

      if (i == 1 && jet->PT < 200.0 && ToggleCut[7])
      {
        NCut[7]++;
        if (!ToggleAllCuts)
          continue;
      }
      for (j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
      {
        object = jet->Constituents.At(j);

        if (object == 0)
          continue;

        if (object->IsA() == GenParticle::Class())
        {
          printf("FUCK\n");
        }
        else if (object->IsA() == Track::Class())
        {
          track = (Track *)object;
          momentum += track->P4();
        }
        else if (object->IsA() == Tower::Class())
        {
          tower = (Tower *)object;
          momentum += tower->P4();
        }
      }
      if (jet->BTag && sqrt(sq(((Jet *)branchAK8Jet->At(1))->Eta - jet->Eta) + sq(((Jet *)branchAK8Jet->At(1))->Phi - jet->Phi)) < 0.8)
      {
        cutLoop++;
      }
    }
    if (cutLoop && ToggleCut[9])
    {
      NCut[9]++;
      if (!ToggleAllCuts)
        continue;
    }

    if (NJetAK4 > 1 && ToggleCut[4])
    {
      cutLoop = 0;
      if (abs(HTmiss.Phi() - ((Jet *)branchAK4Jet->At(0))->Phi) < 0.5)
        cutLoop++;
      if (abs(HTmiss.Phi() - ((Jet *)branchAK4Jet->At(1))->Phi) < 0.5)
        cutLoop++;
      if (NJetAK4 > 2 && abs(HTmiss.Phi() - ((Jet *)branchAK4Jet->At(2))->Phi) < 0.3)
        cutLoop++;
      if (NJetAK4 > 3 && abs(HTmiss.Phi() - ((Jet *)branchAK4Jet->At(3))->Phi) < 0.3)
        cutLoop++;
      if (cutLoop)
      {
        NCut[4]++;
        if (!ToggleAllCuts)
          continue;
      }
    }
    if (HT < 400 && ToggleCut[3])
    {
      NCut[3]++;
      if (!ToggleAllCuts)
        continue;
    }
    survived++;
    events.push_back(event);
    // histMet->Fill(MET);
  }
  printf("\n\nTotal %lld\nSurvived %lld\nEfficiency %f\n", allEntries, survived, ((float)survived) / allEntries);
  for (int i = 0; i < 10; i++)
  {
    printf("%d %d %d\n", i, ToggleCut[i], NCut[i]);
  }
}

void Analysis::AnalyseEvents2(ExRootTreeReader *treeReader, vector<PassedEvent> &events, int MinNumEvents, int status)
{
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchAK4Jet = treeReader->UseBranch("AK4Jets");
  TClonesArray *branchAK8Jet = treeReader->UseBranch("AK8Jets");
  TClonesArray *branchMET = treeReader->UseBranch("MissingET");
  /*
  TClonesArray *branchElectrons[7];
  TClonesArray *branchMuons[7];
  double minPT[] = {13000.0, 200.0, 133.0, 100.0, 80.0, 66.0, 57.0, 10.0};
  for (size_t i = 0; i < 7; i++)
  {
    string bname = "electrons" + to_string(i);
    branchElectrons[i] = treeReader->UseBranch(("electrons" + to_string(i + 1)).c_str());
    branchMuons[i] = treeReader->UseBranch(("muons" + to_string(i + 1)).c_str());
  }
*/
  Long64_t allEntries = treeReader->GetEntries();
  if (MinNumEvents > 0)
  {
    if (MinNumEvents > allEntries)
    {
      cerr << "Not enough Events!" << endl;
      return;
    }
    allEntries = MinNumEvents;
  }
  cout << "** Chain contains " << allEntries << " events" << endl;

  Electron *electron;
  Photon *photon;
  Muon *muon;

  Track *track;
  Tower *tower;

  Jet *jet;
  TObject *object;

  TLorentzVector momentum;

  Long64_t entry;
  double MET, HT, softDropMass;
  Int_t Ntracks, i, j;
  Int_t ToggleAllCuts, cutLoop;
  Int_t NCut[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Int_t ToggleCut[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  ToggleCut[0] = 1; // Filter Zero Weight Events
  ToggleCut[1] = 1; // NJetAK4 >= 2
  ToggleCut[2] = 1; // PT_miss > 300GeV
  ToggleCut[3] = 1; // HT > 400 GeV
  ToggleCut[4] = 1; // Phi(j,HTMiss)>0.5(0.3)
  ToggleCut[5] = 1; //~isolated Photon, Electron, Muon PT > 10 GeV
  ToggleCut[6] = 1; // isolated Tracks mt> 100GeV, pt > 10GeV
  ToggleCut[7] = 1; // 2 AK8 Jets PT > 200 GeV
  ToggleCut[8] = 0; // mJet of 2 AK8 Jets in [10,140]GeV
  ToggleCut[9] = 0; // AngularSep of AK8 and bTagged
  ToggleAllCuts = 0;
  TVector3 HTmiss;
  // allEntries=1000;
  vector<Jet *> AK4Jets;
  vector<Jet *> AK8Jets;

  Long64_t survived = 0;
  boost::progress_display *show_progress = new boost::progress_display(allEntries);
  for (entry = 0; entry < allEntries; ++entry)
  {
    AK4Jets.clear();
    AK8Jets.clear();
    PassedEvent event;
    int nLep = 0;
    int nPho = 0;
    double PTPho = 0.0;
    event.status = status;
    ++(*show_progress);
    vector<int> JetAK4Indices;
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    MET = ((MissingET *)branchMET->At(0))->MET;
    Ntracks = branchTrack->GetEntries();
    for (size_t i = 0; i < branchAK4Jet->GetEntriesFast(); i++)
    {
      jet = (Jet *)branchAK4Jet->At(i);
      if (abs(jet->Eta) < 2.4 && jet->PT > 30)
        AK4Jets.push_back(jet);
    }
    for (size_t i = 0; i < branchAK8Jet->GetEntriesFast(); i++)
    {
      jet = (Jet *)branchAK8Jet->At(i);
      if (abs(jet->Eta) < 2.4 && jet->PT > 30)
        AK8Jets.push_back(jet);
    }

    if (MET == 0.0)
    {
      NCut[0];
      continue;
    }

    if (AK4Jets.size() < 2 && ToggleCut[1])
    {
      NCut[1]++;
      if (!ToggleAllCuts)
        continue;
    }

    if (MET < 100 && ToggleCut[2])
    {
      NCut[2]++;
      continue;
    }

    HT = 0.0;
    HTmiss = TVector3(0.0, 0.0, 0.0);
    cutLoop = 0;
    int cutBTag = 0;
    for (i = 0; i < AK4Jets.size(); ++i)
    {
      HT += AK4Jets[i]->PT;
    }
    if (HT < 400 && ToggleCut[3])
    {
      NCut[3]++;
      continue;
    }

    for (i = 0; i < AK4Jets.size(); ++i)
    {
      jet = AK4Jets[i];
      HTmiss -= TVector3(jet->PT * cos(jet->Phi), jet->PT * sin(jet->Phi), jet->PT * sinh(jet->Eta));
    }

    cutLoop = 0;
    if (abs(HTmiss.Phi() - AK4Jets[0]->Phi) < 0.5)
    {
      NCut[4]++;
      continue;
    }
    if (abs(HTmiss.Phi() - AK4Jets[1]->Phi) < 0.5)
    {
      NCut[4]++;
      continue;
    }
    if (AK4Jets.size() > 2 && abs(HTmiss.Phi() - AK4Jets[2]->Phi) < 0.3)
    {
      NCut[4]++;
      continue;
    }
    if (AK4Jets.size() > 3 && abs(HTmiss.Phi() - AK4Jets[3]->Phi) < 0.3)
    {
      NCut[4]++;
      continue;
    }
    /*
    for (size_t i = 0; i < 7; i++)
    {
      for (size_t lepIndex = 0; lepIndex < branchElectrons[i]->GetEntriesFast(); lepIndex++)
      {
        electron = (Electron *)branchElectrons[i]->At(lepIndex);
        if (minPT[i] > electron->PT && minPT[i + 1] < electron->PT && electron->IsolationVar < 0.1)
          nLep++;
      }
      for (size_t lepIndex = 0; lepIndex < branchMuons[i]->GetEntriesFast(); lepIndex++)
      {
        muon = (Muon *)branchMuons[i]->At(lepIndex);
        if (minPT[i] > muon->PT && minPT[i + 1] < muon->PT && muon->IsolationVar < 0.2)
          nLep++;
      }
    }
*/
    for (i = 0; i < branchPhoton->GetEntriesFast(); ++i)
    {
      photon = (Photon *)branchPhoton->At(i);
      if (photon->Particles.GetEntriesFast() != 1)
        continue;

      if (photon->PT > 10 && ToggleCut[5] && photon->IsolationVar < 1.3 / photon->PT + 0.005)
      {
        nPho++;
        PTPho = photon->PT;
      }
    }
    if (nLep + nPho != 0 && ToggleCut[5])
    {
      NCut[5]++;
      continue;
    }

    cutLoop = 0;
    for (int itrack = 0; itrack < Ntracks; itrack++)
    {
      Track *track = (Track *)branchTrack->At(itrack);
      MissingET *missingET = ((MissingET *)branchMET->At(0));
      double teta = track->Eta;
      double tphi = track->Phi;
      double mT = sqrt(2.0 * track->PT * missingET->MET * (1.0 - cos(missingET->Phi - track->Phi)));
      if (mT > 100 || abs(track->Eta) > 2.4)
        continue;

      double conePT = 0.0;
      for (int icandtrack = 0; icandtrack < Ntracks; icandtrack++)
      {
        if (itrack == icandtrack)
          continue;
        Track *candTrack = (Track *)branchTrack->At(icandtrack);
        if (sqrt(sq(teta - candTrack->Eta) + sq(tphi - candTrack->Phi)) < 0.3)
          conePT += candTrack->PT;
      }
      if (abs(track->PID) < 20)
      {
        if (track->PT > 5 && conePT / track->PT < 0.2)
          cutLoop = 1;
      }
      else
      {
        if (track->PT > 10 && conePT / track->PT < 0.1)
          cutLoop = 1;
      }
    }
    if (cutLoop && ToggleCut[6])
    {
      NCut[6]++;
      continue;
    }

    if ((AK8Jets.size() < 2 || AK8Jets[1]->PT < 200) && ToggleCut[7])
    {
      NCut[7]++;
      continue;
    }

    softDropMass = AK8Jets[0]->SoftDroppedJet.Mag();
    event.mj1 = softDropMass;
    if ((softDropMass < 40.0 || softDropMass > 140.0) && ToggleCut[8])
    {
      NCut[8]++;
      if (!ToggleAllCuts)
        continue;
    }
    softDropMass = AK8Jets[1]->SoftDroppedJet.Mag();
    event.mj2 = softDropMass;
    if ((softDropMass < 40.0 || softDropMass > 140.0) && ToggleCut[8])
    {
      NCut[8]++;
      if (!ToggleAllCuts)
        continue;
    }
    cutLoop = 0;
    for (i = 0; i < AK4Jets.size(); ++i)
    {
      jet = AK4Jets[i];
      if (jet->BTag && sqrt(sq(AK8Jets[1]->Eta - jet->Eta) + sq(AK8Jets[1]->Phi - jet->Phi)) < 0.8)
      {
        cutBTag++;
      }
    }
    for (i = 0; i < AK8Jets.size(); ++i)
    {
      if (i == 1)
        continue;
      jet = AK8Jets[i];
      if (jet->BTag && sqrt(sq(AK8Jets[1]->Eta - jet->Eta) + sq(AK8Jets[1]->Phi - jet->Phi)) < 0.8)
      {
        cutBTag++;
      }
    }
    if (cutBTag && ToggleCut[9])
    {
      NCut[9]++;
      continue;
    }

    event.ptmiss = MET;
    survived++;
    events.push_back(event);
  }
  printf("\n\nTotal %lld\nSurvived %lld\nEfficiency %f\n", allEntries, survived, ((float)survived) / allEntries);
  for (int i = 0; i < 10; i++)
  {
    printf("%d %d %d\n", i, ToggleCut[i], NCut[i]);
  }
}

void Analysis::SampleFromEventsNewPassedEventDef(vector<PassedEvent> &events, double intLumi)
{
  double loX[] = {4.248, 1.375, 4.584, 1.964};
  double nloX[] = {5.410, 1.773, 6.741, 5.218};
  double X[] = {181.0, 511.0, 78.0, 2344.0};
  double ST[] = {2284, 1407, 1172, 10};
  double NT[] = {147143876, 407148516, 84600510, 0};
  vector<vector<PassedEvent>> rawEvents(3);
  vector<PassedEvent> importedEvents;
  RootIO::ReadEventsNewPassedEventDef("/home/gerrit/Documents/PHD/SUSY4Cathode/ClusterTestRun/root-on-vscode/AllCuts.root", importedEvents, 0);

  for (size_t i = 0; i < importedEvents.size(); i++)
  {
    importedEvents[i].status--;
    rawEvents[importedEvents[i].status].push_back(importedEvents[i]);
  }

  cout << "Done Reading"<< endl;
  random_device rd;
  mt19937_64 gen(0); // 0
  uniform_real_distribution<double> dist(0.0, 1.0);
  for (size_t i = 0; i < 3; i++)
  {
    double Nexpected = nloX[i] / loX[i] * intLumi* X[i] * ST[i] / NT[i];
    cout << "Weight: " << Nexpected / ST[i] << "  Nexpected: " << Nexpected << "  Navail: " << rawEvents[i].size() << endl;
    if (Nexpected / ST[i] > 1)
      cout << "Not enough Events" << endl;
    for (size_t iEvent = 0; iEvent < rawEvents[i].size(); iEvent++)
    {
      if (dist(gen) < Nexpected / ST[i])
      {
        events.push_back(rawEvents[i][iEvent]);
      }
    }
  }
  cout << "Total events : " << events.size() << endl;
}
void Analysis::SampleFromEventsNewPassedEventDef(vector<PassedEvent> &events)
{
  Analysis::SampleFromEventsNewPassedEventDef(events, 300.0e3);
}

void Analysis::SampleFromEvents(vector<PassedEvent> &events, double intLumi)
{
  double loX[] = {4.248, 1.375, 4.584, 1.964};
  double nloX[] = {5.410, 1.773, 6.741, 5.218};
  double X[] = {181.0, 511.0, 78.0, 2344.0};
  double ST[] = {2284, 1407, 1172, 10};
  double NT[] = {36753049, 147113599, 42475284, 424699605};
  vector<vector<PassedEvent>> rawEvents(4);
  RootIO::ReadEvents("RootFiles/Z-Total.root", rawEvents[0], 0);
  RootIO::ReadEvents("RootFiles/W-Total.root", rawEvents[1], 0);
  RootIO::ReadEvents("RootFiles/tt-Total.root", rawEvents[2], 0);
  RootIO::ReadEvents("RootFiles/A-Total.root", rawEvents[3]);
  random_device rd;
  mt19937_64 gen(0); // 0
  uniform_real_distribution<double> dist(0.0, 1.0);
  for (size_t i = 0; i < 4; i++)
  {
    double Nexpected = nloX[i] / loX[i] * intLumi * X[i] * ST[i] / NT[i];
    cout << "Weight: " << Nexpected / ST[i] << "  Nexpected: " << Nexpected << "  Navail: " << rawEvents[i].size() << endl;
    if (Nexpected / ST[i] > 1)
      cout << "Not enough Events" << endl;
    for (size_t iEvent = 0; iEvent < rawEvents[i].size(); iEvent++)
    {
      if (dist(gen) < Nexpected / ST[i])
      {
        events.push_back(rawEvents[i][iEvent]);
      }
    }
  }
}

void Analysis::SampleFromEvents(vector<PassedEvent> &events)
{
  Analysis::SampleFromEvents(events, 300.0e3);
}

double Analysis::TheoryXSection(double m_Gluino)
{
  return pow(10.0, 1.8821e-7 * sq(m_Gluino) - 3.0554e-3 * m_Gluino + 2.3574);
}

double Analysis::TheoryXSectionUpper(double m_Gluino)
{
  return pow(10.0, 2.1124e-7 * sq(m_Gluino) - 3.1114e-3 * m_Gluino + 2.4538);
}

double Analysis::TheoryXSectionLower(double m_Gluino)
{
  return pow(10.0, 2.3056e-7 * sq(m_Gluino) - 3.2642e-3 * m_Gluino + 2.5049);
}

TF1 *Analysis::fitFunc;
double Analysis::excludeSignalRegion(Double_t *x, Double_t *par)
{
  if (x[0] > 70 && x[0] < 100)
  {
    TF1::RejectedPoint();
    return 0;
  }
  fitFunc->SetParameters(par);
  return fitFunc->Eval(x[0]);
}

void Analysis::FitCherbyshev(vector<PassedEvent> events, double &bNorm, double &bNormErr, vector<double> &params)
{
  double sysError = 0.0;
  double statError = 0.0;
  TCanvas *c1 = new TCanvas("c1FitCherbyshev", "c1",200,180);
  c1->cd();
  TH1F *h1 = new TH1F("h1", "", 20, 40.0, 140.0);

  for (size_t i = 0; i < events.size(); i++)
  {
    if (events[i].mj2 < 70.0 || events[i].mj2 > 100.0 || events[i].status >= 10)
      continue;
    h1->Fill(events[i].mj1);
  }
  vector<double> yield;
  for (int n = 1; n < 5; n++)
  {
    string polyName = "chebyshev" + to_string(n);
    fitFunc = (TF1 *)gROOT->GetFunction(polyName.c_str());
    TF1 *backGroundFit = new TF1(("Fit" + polyName).c_str(), excludeSignalRegion, 40, 140, fitFunc->GetNpar());
    for (size_t iPar = 0; iPar < fitFunc->GetNpar(); iPar++)
    {
      backGroundFit->SetParameter(iPar, 1.0);
    }

    cout << "++++++++++++++" << endl;
    cout << polyName << endl;
    h1->Fit(backGroundFit, "IQ");
    yield.push_back(fitFunc->Integral(70, 100) / 5);
    if (n == 1)
    {
      h1->Draw("E");
      h1->GetYaxis()->SetRangeUser(0, 150);
      params.push_back(backGroundFit->GetParameter(0));
      params.push_back(backGroundFit->GetParameter(1));
      params.push_back(backGroundFit->GetParError(0));
      params.push_back(backGroundFit->GetParError(1));
      c1->SaveAs("Plots/Ch1.pdf");
      TH1F *pseudoyield = new TH1F("PseudoExperiments", "Pseudoexperiments BG normalisation", 100, yield[0] - 100, yield[0] + 100);
      int nPseudo = 10000;
      TF1 *backgroundModel = new TF1("backgroundmodel", "[1]*x+[0]", 40.0, 140.0);
      backgroundModel->SetParameter(0, params[0]);
      backgroundModel->SetParameter(1, params[1]);
      for (int i = 0; i < nPseudo; i++)
      {
        TH1F *pseudo = new TH1F("pseudo1", "", 20, 40.0, 140.0);
        pseudo->FillRandom("backgroundmodel", h1->Integral());
        pseudo->Fit(backGroundFit, "IQ");
        pseudoyield->Fill(fitFunc->Integral(70, 100) / 5);
        delete pseudo;
      }
      TFitResultPtr fitptr=pseudoyield->Fit("gaus", "IQS");
      statError = pseudoyield->GetFunction("gaus")->GetParameter(2);
      fitptr->Print();
      pseudoyield->GetXaxis()->SetTitle("B_{norm}");
      gStyle->SetStatFontSize(0.1);
      pseudoyield->Draw();
      gPad->Update();
      

      gStyle->SetOptFit(1);
      c1->SaveAs("Plots/PseudoExperiments.pdf");
    }
  }
  bNorm = yield[0];
  for (int i = 0; i < 4; i++)
  {
    sysError = max(sysError, abs(yield[0] - yield[i]));
  }
  bNormErr = sqrt(sq(sysError) + sq(statError));
  cout << "B_Norm " << bNorm << "\nStat : " << statError << "\nSys : " << sysError << "\nTotal : " << bNormErr << endl;
}

void Analysis::CalcTransferFactor(vector<PassedEvent> events, double bNorm, double bNormErr, vector<double> &NCRi, double &transferFactor, double &transferFactorErr)
{
  TCanvas *c1 = new TCanvas("c1Trans", "c1");
  c1->cd();
  const double bins[7] = {300.0, 450.0, 600.0, 800.0, 1000.0, 1200.0, 2000.0};
  TH1 *histCR = new TH1F("ptmiss_CR", "ptmiss CR", 6, bins);
  for (size_t i = 0; i < events.size(); i++)
  {
    PassedEvent event = events[i];
    if (event.status >= 10)
      continue;
    if (event.mj1 < 70 || event.mj1 > 100)
      if (event.mj2 < 70 || event.mj2 > 100)
        histCR->Fill(event.ptmiss);
  }
  cout << histCR->GetBinContent(1) << endl;
  histCR->Draw();
  transferFactor = bNorm / histCR->Integral();
  transferFactorErr = bNormErr / histCR->Integral();
  NCRi.push_back(histCR->GetBinContent(1));
  NCRi.push_back(histCR->GetBinContent(2));
  NCRi.push_back(histCR->GetBinContent(3));
  NCRi.push_back(histCR->GetBinContent(4));
  NCRi.push_back(histCR->GetBinContent(5));
  NCRi.push_back(histCR->GetBinContent(6));

  //  c1->SaveAs("Plots/CR.pdf");
  cout << "\n\n++++++++++++++++++++++++++ \n TransferFactor: "<< transferFactor << "  +-  " << transferFactorErr<<endl;
}

void Analysis::GetEventsInSignalRegion(vector<PassedEvent> events, vector<PassedEvent> &signalevents)
{
  for (int i = 0; i < events.size(); i++)
  {
    if (events[i].status < 10 && events[i].mj1 > 70 && events[i].mj1 < 100 && events[i].mj2 > 70 && events[i].mj2 < 100)
      signalevents.push_back(events[i]);
  }
}

void Analysis::GetSRBinContent(vector<PassedEvent> &signalevents, vector<double> &binContent)
{
  const double bins[7] = {300.0, 450.0, 600.0, 800.0, 1000.0, 1200.0, 2000.0};
  TH1 *histSR = new TH1F("SignalRegionTMP", "SR", 6, bins);
  for (size_t i = 0; i < signalevents.size(); i++)
  {
    histSR->Fill(signalevents[i].ptmiss);
  }
  for (size_t i = 1; i < 7; i++)
  {
    binContent.push_back(histSR->GetBinContent(i));
  }

  delete histSR;
}