#include "Analysis.h"

double Analysis::sq(double x)
{
  return x * x;
}

void Analysis::AnalyseEvents(ExRootTreeReader *treeReader, vector<PassedEvent> &events, int MinNumEvents, int status)
{
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
  ToggleCut[0] = 1; //Filter Zero Weight Events
  ToggleCut[1] = 1; //NJetAK4 >= 2
  ToggleCut[2] = 1; //PT_miss > 300GeV
  ToggleCut[3] = 1; //HT > 400 GeV
  ToggleCut[4] = 1; //Phi(j,HTMiss)>0.5(0.3)
  ToggleCut[5] = 1; //~isolated Photon, Electron, Muon PT > 10 GeV
  ToggleCut[6] = 1; //isolated Tracks mt> 100GeV, pt > 10GeV
  ToggleCut[7] = 1; //2 AK8 Jets PT > 200 GeV
  ToggleCut[8] = 1; //mJet of 2 AK8 Jets in [10,140]GeV
  ToggleCut[9] = 1; //AngularSep of AK8 and bTagged
  ToggleAllCuts = 0;
  TVector3 HTmiss;
  // allEntries=1000;

  Long64_t survived = 0;
  boost::progress_display *show_progress = new boost::progress_display(allEntries);
  for (entry = 0; entry < allEntries; ++entry)
  {
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
    if (MET < 300.0 && ToggleCut[2])
    {
      NCut[2]++;
      if (!ToggleAllCuts)
        continue;
    }
    event.ptmiss = MET;
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

      //particle = (GenParticle*) photon->Particles.At(0);
      if (photon->PT > 10 && ToggleCut[5])
      {
        NCut[5]++;
        cutLoop++;
        if (!ToggleAllCuts)
          break;
      }
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
      //NCut[7]++;
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
    //histMet->Fill(MET);
  }
  printf("\n\nTotal %lld\nSurvived %lld\nEfficiency %f\n", allEntries, survived, ((float)survived) / allEntries);
  for (int i = 0; i < 10; i++)
  {
    printf("%d %d %d\n", i, ToggleCut[i], NCut[i]);
  }
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
  ToggleCut[0] = 1; //Filter Zero Weight Events
  ToggleCut[1] = 1; //NJetAK4 >= 2
  ToggleCut[2] = 1; //PT_photon > 200GeV
  ToggleCut[3] = 1; //HT > 400 GeV
  ToggleCut[4] = 1; //Phi(j,HTMiss)>0.5(0.3)
  ToggleCut[5] = 1; //~isolated Photon, Electron, Muon PT > 10 GeV
  ToggleCut[6] = 1; //isolated Tracks mt> 100GeV, pt > 10GeV
  ToggleCut[7] = 1; //2 AK8 Jets PT > 200 GeV
  ToggleCut[8] = 1; //mJet of 2 AK8 Jets in [10,140]GeV
  ToggleCut[9] = 0; //AngularSep of AK8 and bTagged
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
      event.ptmiss = max(event.ptmiss, (double)photon->PT);
      //particle = (GenParticle*) photon->Particles.At(0);
      if (photon->PT > 10 && ToggleCut[5])
      {
        if (nPhoton)
        {
          NCut[5]++;
          cutLoop++;
          if (!ToggleAllCuts)
            break;
        }
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
      //NCut[7]++;
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
    //histMet->Fill(MET);
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
  ToggleCut[0] = 1; //Filter Zero Weight Events
  ToggleCut[1] = 1; //NJetAK4 >= 2
  ToggleCut[2] = 1; //PT_miss > 200GeV
  ToggleCut[3] = 1; //HT > 400 GeV
  ToggleCut[4] = 1; //Phi(j,HTMiss)>0.5(0.3)
  ToggleCut[5] = 1; //~isolated Photon, Electron, Muon PT > 10 GeV
  ToggleCut[6] = 1; //isolated Tracks mt> 100GeV, pt > 10GeV
  ToggleCut[7] = 1; //2 AK8 Jets PT > 200 GeV
  ToggleCut[8] = 1; //mJet of 2 AK8 Jets in [10,140]GeV
  ToggleCut[9] = 1; //AngularSep of AK8 and bTagged
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
    cutLoop = 0;
    for (i = 0; i < branchElectron->GetEntriesFast(); ++i)
    {
      electron = (Electron *)branchElectron->At(i);

      if (electron->PT > 10 && ToggleCut[5] && electron->IsolationVar > 0.1)
      {
        if (nLepton)
        {
          NCut[5]++;
          if (!ToggleAllCuts)
            cutLoop++;
        }
        nLepton++;
      }
    }

    // Loop over all photons in event
    for (i = 0; i < branchPhoton->GetEntriesFast(); ++i)
    {
      photon = (Photon *)branchPhoton->At(i);

      // skip photons with references to multiple particles
      if (photon->Particles.GetEntriesFast() != 1)
        continue;

      //particle = (GenParticle*) photon->Particles.At(0);
      if (photon->PT > 10 && ToggleCut[5])
      {
        NCut[5]++;
        if (!ToggleAllCuts)
          continue;
      }
    }

    // Loop over all muons in event
    for (i = 0; i < branchMuon->GetEntriesFast(); ++i)
    {
      muon = (Muon *)branchMuon->At(i);
      if (muon->PT > 10 && ToggleCut[5] && muon->IsolationVar > 0.2)
      {
        if (nLepton)
        {
          NCut[5]++;
          if (!ToggleAllCuts)
            cutLoop++;
        }
        nLepton++;
      }
    }
    if (cutLoop && nLepton != 1)
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
      //NCut[7]++;
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
    //histMet->Fill(MET);
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
  TClonesArray *branchIsoEl[7];

  TClonesArray *branchIsoMu[7];
  string brnachname;
  double maxPT[8] = {1e6, 200, 133, 100, 80, 66, 57, 10};
  for (int n = 0; n < 7; n++)
  {
    brnachname = "electrons";
    brnachname += to_string(n + 1);
    branchIsoEl[n] = treeReader->UseBranch(brnachname.c_str());
    brnachname = "muons";
    brnachname += to_string(n + 1);
    branchIsoMu[n] = treeReader->UseBranch(brnachname.c_str());
  }

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
  ToggleCut[0] = 1; //Filter Zero Weight Events
  ToggleCut[1] = 1; //NJetAK4 >= 2
  ToggleCut[2] = 1; //PT_miss > 300GeV
  ToggleCut[3] = 1; //HT > 400 GeV
  ToggleCut[4] = 1; //Phi(j,HTMiss)>0.5(0.3)
  ToggleCut[5] = 1; //~isolated Photon, Electron, Muon PT > 10 GeV
  ToggleCut[6] = 1; //isolated Tracks mt> 100GeV, pt > 10GeV
  ToggleCut[7] = 1; //2 AK8 Jets PT > 200 GeV
  ToggleCut[8] = 1; //mJet of 2 AK8 Jets in [10,140]GeV
  ToggleCut[9] = 1; //AngularSep of AK8 and bTagged
  ToggleAllCuts = 0;
  TVector3 HTmiss;
  // allEntries=1000;

  Long64_t survived = 0;
  boost::progress_display *show_progress = new boost::progress_display(allEntries);
  for (entry = 0; entry < allEntries; ++entry)
  {

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
    if (false)
    {
      cutLoop = 0;
      // Loop over all electrons in event
      for (i = 0; i < branchElectron->GetEntriesFast(); ++i)
      {
        electron = (Electron *)branchElectron->At(i);

        if (electron->PT > 10 && ToggleCut[5] && electron->IsolationVar > 0.1)
        {
          cutLoop++;
        }
      }
      // Loop over all muons in event
      for (i = 0; i < branchMuon->GetEntriesFast(); ++i)
      {
        muon = (Muon *)branchMuon->At(i);
        if (muon->PT > 10 && ToggleCut[5] && muon->IsolationVar > 0.2)
        {
          cutLoop++;
        }
      }
      if (cutLoop)
      {
        NCut[5]++;
        if (!ToggleAllCuts)
          continue;
      }
    } //38
    else
    {
      cutLoop = 0;
      for (int n = 0; n < 7; n++)
      {
        for (i = 0; i < branchIsoEl[n]->GetEntriesFast(); i++)
        {
          electron = (Electron *)branchIsoEl[n]->At(i);
          if (electron->PT < maxPT[n] && electron->PT > maxPT[n + 1] && electron->IsolationVar < 0.1)
            cutLoop++;
        }
        for (i = 0; i < branchIsoMu[n]->GetEntriesFast(); i++)
        {
          muon = (Muon *)branchIsoMu[n]->At(i);
          if (muon->PT < maxPT[n] && muon->PT > maxPT[n + 1] && muon->IsolationVar < 0.2)
            cutLoop++;
        }
      }
      if (cutLoop)
        NCut[5]++;
    }
    // Loop over all photons in event
    for (i = 0; i < branchPhoton->GetEntriesFast(); ++i)
    {
      photon = (Photon *)branchPhoton->At(i);

      // skip photons with references to multiple particles
      if (photon->Particles.GetEntriesFast() != 1)
        continue;

      //particle = (GenParticle*) photon->Particles.At(0);
      if (photon->PT > 10 && ToggleCut[5])
      {
        NCut[5]++;
        if (!ToggleAllCuts)
          continue;
      }
    }

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
      //NCut[7]++;
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
    //histMet->Fill(MET);
  }
  printf("\n\nTotal %lld\nSurvived %lld\nEfficiency %f\n", allEntries, survived, ((float)survived) / allEntries);
  for (int i = 0; i < 10; i++)
  {
    printf("%d %d %d\n", i, ToggleCut[i], NCut[i]);
  }
}

double Analysis::FindNormalisation(const char *filename)
{

  TChain chain("Delphes");
  int maxFiles = 100;
  RootIO::GetRootFiles(filename, chain, maxFiles);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

  TClonesArray *branchAK4Jet = treeReader->UseBranch("AK4Jets");
  TClonesArray *branchMET = treeReader->UseBranch("MissingET");

  double HT = 0.0;
  Int_t NJetAK4;
  Long64_t allEntries = treeReader->GetEntries();
  //allEntries=10000;
  cout << "** Chain contains " << allEntries << " events" << endl;
  int PassesHTCut = 0;
  boost::progress_display *show_progress = new boost::progress_display(allEntries);
  for (Long64_t entry = 0; entry < allEntries; ++entry)
  {
    ++(*show_progress);
    HT = 0.0;
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    NJetAK4 = branchAK4Jet->GetEntriesFast();
    for (int i = 0; i < NJetAK4; ++i)
    {
      Jet *jet = (Jet *)branchAK4Jet->At(i);
      HT += jet->PT;
    }
    MissingET *met = (MissingET *)branchMET->At(0);
    if (HT > 400.0 && met->MET > 300)
      PassesHTCut++;
  }
  printf("HT:\n%E\n%d rel. uncert. %E\n%lld\n\n", 1.0 * PassesHTCut / allEntries, PassesHTCut, 1.0 / sqrt(PassesHTCut), allEntries);
  return 1.0 * PassesHTCut / allEntries;
}

void Analysis::RunBigAnalysis()
{
  ExRootTreeReader *treeReader;
  double integratedLuminosity = 137000.0; //in pb^-1
  double totalZcrosssection = 5.4100e4;
  double totalWcrosssection = 1.77300e5;
  double totalttcrosssection = 674.1;
  cout << "FIND NORMS" << endl;
  size_t numZEvents = 0;// integratedLuminosity * totalZcrosssection * FindNormalisation("/mnt/5c451946-c244-49ab-9fc5-1aaca8739b2a/ZNorm-Cluster/ROOTFILES.txt") / FindNormalisation("/mnt/5c451946-c244-49ab-9fc5-1aaca8739b2a/Z-Cluster/ROOTFILES.txt");
  size_t numWEvents = 0;// integratedLuminosity * totalWcrosssection * FindNormalisation("/mnt/5c451946-c244-49ab-9fc5-1aaca8739b2a/WNorm-Cluster/ROOTFILES.txt") / FindNormalisation("/media/gerrit/Files/DelphesEvents/W-Cluster/ROOTFILES.txt");
  size_t numttEvents = integratedLuminosity * totalttcrosssection * FindNormalisation("/mnt/5c451946-c244-49ab-9fc5-1aaca8739b2a/ttNorm-Cluster/ROOTFILES.txt") / FindNormalisation("/media/gerrit/Files/DelphesEvents/tt-Cluster/ROOTFILES.txt");

  vector<PassedEvent> events;


  cout << numZEvents << endl
       << numWEvents << endl;
  cout << "ANALYSE Z" << endl;
  TChain chainZ("Delphes");
  RootIO::GetRootFiles("/mnt/5c451946-c244-49ab-9fc5-1aaca8739b2a/Z-Cluster/ROOTFILES.txt", chainZ);
  treeReader = new ExRootTreeReader(&chainZ);
  cout << "Need : " << numZEvents << " Events" << endl;
  cout << "Have : " << treeReader->GetEntries() << " Events" << endl;
  Analysis::AnalyseEvents(treeReader, events, numZEvents, 0);
  RootIO::SaveEvents(events);
  chainZ.Reset();

  cout << "ANALYSE W" << endl;
  TChain chainW("Delphes");
  RootIO::GetRootFiles("/media/gerrit/Files/DelphesEvents/W-Cluster/ROOTFILES.txt", chainW);
  treeReader = new ExRootTreeReader(&chainW);
  cout << "Need : " << numWEvents << " Events" << endl;
  cout << "Have : " << treeReader->GetEntries() << " Events" << endl;
  Analysis::AnalyseEvents(treeReader, events, numWEvents, 1);
  RootIO::SaveEvents(events);
  chainW.Reset();

  cout << "ANALYSE tt" << endl;
  TChain chaintt("Delphes");
  RootIO::GetRootFiles("/media/gerrit/Files/DelphesEvents/tt-Cluster/ROOTFILES.txt", chaintt);
  treeReader = new ExRootTreeReader(&chaintt);
  cout << "Need : " << numttEvents << " Events" << endl;
  cout << "Have : " << treeReader->GetEntries() << " Events" << endl;
  Analysis::AnalyseEvents(treeReader, events, numttEvents, 2);
  RootIO::SaveEvents(events);
  chaintt.Reset();

  cout << "Single Lepton Z" << endl;
  TChain chainLepton("Delphes");
  RootIO::GetRootFiles("/mnt/5c451946-c244-49ab-9fc5-1aaca8739b2a/Z-Cluster/ROOTFILES.txt", chainLepton);
  treeReader = new ExRootTreeReader(&chainLepton);
  Analysis::AnalyseEventsSingleLeptonSample(treeReader, events, numZEvents);
  RootIO::SaveEvents(events);
  chainLepton.Reset();

  cout << "Single Lepton W" << endl;
  RootIO::GetRootFiles("/media/gerrit/Files/DelphesEvents/W-Cluster/ROOTFILES.txt", chainLepton);
  treeReader = new ExRootTreeReader(&chainLepton);
  Analysis::AnalyseEventsSingleLeptonSample(treeReader, events, numWEvents);
  RootIO::SaveEvents(events);
  chainLepton.Reset();

  cout << "Single Lepton tt" << endl;
  RootIO::GetRootFiles("/media/gerrit/Files/DelphesEvents/tt-Cluster/ROOTFILES.txt", chainLepton);
  treeReader = new ExRootTreeReader(&chainLepton);
  Analysis::AnalyseEventsSingleLeptonSample(treeReader, events, numttEvents);
  RootIO::SaveEvents(events);
  chainLepton.Reset();

  cout << "Single Photon Z" << endl;
  TChain chainPhoton("Delphes");
  RootIO::GetRootFiles("/mnt/5c451946-c244-49ab-9fc5-1aaca8739b2a/Z-Cluster/ROOTFILES.txt", chainPhoton);
  treeReader = new ExRootTreeReader(&chainPhoton);
  Analysis::AnalyseEventsSinglePhotonSample(treeReader, events, numZEvents);
  RootIO::SaveEvents(events);
  chainPhoton.Reset();

  cout << "Single Photon W" << endl;
  RootIO::GetRootFiles("/media/gerrit/Files/DelphesEvents/W-Cluster/ROOTFILES.txt", chainPhoton);
  treeReader = new ExRootTreeReader(&chainPhoton);
  Analysis::AnalyseEventsSinglePhotonSample(treeReader, events, numWEvents);
  RootIO::SaveEvents(events);
  chainPhoton.Reset();

  cout << "Single Photon tt" << endl;
  RootIO::GetRootFiles("/media/gerrit/Files/DelphesEvents/tt-Cluster/ROOTFILES.txt", chainPhoton);
  treeReader = new ExRootTreeReader(&chainPhoton);
  Analysis::AnalyseEventsSinglePhotonSample(treeReader, events, numttEvents);
  RootIO::SaveEvents(events);
  chainPhoton.Reset();

  RootIO::SaveEvents(events);
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

void Analysis::FitCherbyshev(vector<PassedEvent> &events, double &bNorm, vector<double> &params)
{
  double sysError = 0.0;
  double statError = 0.0;
  TCanvas *c1 = new TCanvas("c1FitCherbyshev", "c1");
  c1->cd();
  TH1F *h1 = new TH1F("h1", "", 20, 40.0, 140.0);
  for (size_t i = 0; i < events.size(); i++)
  {
    if (events[i].mj2 < 70.0 || events[i].mj2 > 100.0 || events[i].status >=10)
      continue;
    h1->Fill(events[i].mj1);
  }
  vector<double> yield;
  for (int n = 1; n < 5; n++)
  {
    string polyName = "chebyshev" + to_string(n);
    fitFunc = (TF1 *)gROOT->GetFunction(polyName.c_str());
    TF1 *backGroundFit = new TF1(("Fit" + polyName).c_str(), excludeSignalRegion, 40, 140, fitFunc->GetNpar());
    h1->Fit(backGroundFit, "IQ");
    yield.push_back(fitFunc->Integral(70, 100) / 5);
    if (n == 1)
    {
      h1->Draw();
      params.push_back(backGroundFit->GetParameter(0));
      params.push_back(backGroundFit->GetParameter(1));
      c1->SaveAs("Plots/Ch1.eps");
      TH1F *pseudoyield = new TH1F("pseudoExperiments", "", 100, yield[0] - 100, yield[0] + 100);
      int nPseudo = 1000;
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
      pseudoyield->Fit("gaus", "IQ");
      statError = pseudoyield->GetFunction("gaus")->GetParameter(2);
      pseudoyield->Draw();
      c1->SaveAs("Plots/PseudoExperiments.eps");
    }
  }
  bNorm = yield[0];
  for (int i = 1; i < 4; i++)
  {
    sysError = max(sysError, abs(yield[0] - yield[i]));
  }
  cout << "Stat : " << statError << "\nSys : " << sysError << "\nTotal : " << sqrt(sq(sysError) + sq(statError)) << endl;
}

void Analysis::CalcTransferFactor(vector<PassedEvent> &events, double &bNorm, vector<double> &bi, double &transferFactor)
{
  TCanvas *c1 = new TCanvas("c1Trans", "c1");
  c1->cd();
  const double bins[7] = {300.0, 450.0, 600.0, 800.0, 1000.0, 1200.0, 2000.0};
  TH1 *histCR = new TH1F("ptmiss_CR", "ptmiss CR", 6, bins);
  for (size_t i = 0; i < events.size(); i++)
  {
    PassedEvent event = events[i];
    if(event.status>=10)
      continue;
    if (event.mj1 < 70 || event.mj1 > 100)
      if (event.mj2 < 70 || event.mj2 > 100)
        histCR->Fill(event.ptmiss);
  }
  cout << histCR->GetBinContent(1) << endl;
  histCR->Draw();
  transferFactor = bNorm / histCR->Integral();
  bi.push_back(transferFactor * histCR->GetBinContent(1));
  bi.push_back(transferFactor * histCR->GetBinContent(2));
  bi.push_back(transferFactor * histCR->GetBinContent(3));
  bi.push_back(transferFactor * histCR->GetBinContent(4));
  bi.push_back(transferFactor * histCR->GetBinContent(5));
  bi.push_back(transferFactor * histCR->GetBinContent(6));

  c1->SaveAs("Plots/CR.eps");
}