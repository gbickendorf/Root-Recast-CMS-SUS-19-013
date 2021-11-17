#include "Analysis.h"

double Analysis::sq(double x)
{return x*x;}

void Analysis::AnalyseEvents(ExRootTreeReader *treeReader, vector<PassedEvent> &events , size_t MinNumEvents, int status)
{
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchAK4Jet = treeReader->UseBranch("AK4Jets");
  TClonesArray *branchAK8Jet = treeReader->UseBranch("AK8Jets");
  TClonesArray *branchMET = treeReader->UseBranch("MissingET");

  Long64_t allEntries = treeReader->GetEntries();
  cout << "Need : "<<MinNumEvents <<" Events"<<endl;
  cout << "Have : "<< allEntries << " Events"<<endl;
  if(MinNumEvents >0)
  {
    if(MinNumEvents>allEntries)
    {
      cerr<<"Not enough Events!"<<endl;
      return;
    }
    allEntries=MinNumEvents;
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
  size_t NJetAK4,Ntracks, i, j;
  Int_t ToggleAllCuts,cutLoop;
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
  boost::progress_display * show_progress = new boost::progress_display(allEntries);
  for(entry = 0; entry < allEntries; ++entry)  
  {
    PassedEvent event;
    event.status=status;
    ++(*show_progress);
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
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
    event.ptmiss=MET;
    // Loop over all electrons in event
    for(i = 0; i < branchElectron->GetEntriesFast(); ++i)
    {
      electron = (Electron*) branchElectron->At(i);
      
      if(electron->PT > 10 && ToggleCut[5]&&electron->IsolationVar>0.1)
      {
        NCut[5]++;
        if(!ToggleAllCuts)
          continue;
      }
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
    }

    // Loop over all muons in event
    for(i = 0; i < branchMuon->GetEntriesFast(); ++i)
    {
      muon = (Muon*) branchMuon->At(i);
      if(muon->PT > 10 && ToggleCut[5]&&muon->IsolationVar>0.2)
      {
        NCut[5]++;
        if(!ToggleAllCuts)
          continue;
      }
    }
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
    event.mj1=softDropMass;
    if((softDropMass<40.0 || softDropMass >140.0) && ToggleCut[8])
    {
      NCut[8]++;
      if(!ToggleAllCuts)
        continue;
    }
    softDropMass=((Jet*) branchAK8Jet->At(1))->SoftDroppedJet.Mag();
    event.mj2=softDropMass;
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
      for(j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
      {
        object = jet->Constituents.At(j);

        if(object == 0) continue;

        if(object->IsA() == GenParticle::Class())
        {
          printf("FUCK\n");
        }
        else if(object->IsA() == Track::Class())
        {
          track = (Track*) object;
          momentum += track->P4();
        }
        else if(object->IsA() == Tower::Class())
        {
          tower = (Tower*) object;
          momentum += tower->P4();
        }
      }
      if(jet->BTag&&sqrt(sq(((Jet*) branchAK8Jet->At(1))->Eta-jet->Eta)+sq(((Jet*) branchAK8Jet->At(1))->Phi-jet->Phi))<0.8)
      {
        cutLoop++;
      }
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
      if(abs(HTmiss.Phi()-((Jet*)branchAK4Jet->At(0))->Phi)<0.5)
        cutLoop++;
      if(abs(HTmiss.Phi()-((Jet*)branchAK4Jet->At(1))->Phi)<0.5)
        cutLoop++;
      if(NJetAK4>2 &&abs( HTmiss.Phi()-((Jet*)branchAK4Jet->At(2))->Phi)<0.3)
        cutLoop++;
      if(NJetAK4>3 &&abs( HTmiss.Phi()-((Jet*)branchAK4Jet->At(3))->Phi)<0.3)
        cutLoop++;
      if(cutLoop)
      {
        NCut[4]++;
        if(!ToggleAllCuts)
          continue;
      }
    }
    if(HT<400&&ToggleCut[3])
    {
      NCut[3]++;
      if(!ToggleAllCuts)
        continue;
    }
    survived++;
    events.push_back(event);
    //histMet->Fill(MET);
  }
  printf("\n\nTotal %lld\nSurvived %lld\nEfficiency %f\n",allEntries,survived,((float)survived)/allEntries);
  for(int i = 0; i < 10; i++)
  {
    printf("%d %d %d\n",i,ToggleCut[i],NCut[i]);
  }
}

double Analysis::FindNormalisation(const char * filename)
{
  
  TChain chain("Delphes");
  int maxFiles =10;
  std::ifstream file(filename);
  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      chain.Add(line.c_str());
      if(!--maxFiles)
        break;

    }
    file.close();
  }
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

  TClonesArray *branchAK4Jet = treeReader->UseBranch("AK4Jets");

  double HT=0.0;
  Int_t NJetAK4;
  Long64_t allEntries = treeReader->GetEntries();
  //allEntries=10000;
  cout << "** Chain contains " << allEntries << " events" << endl;
  int PassesHTCut=0;
  boost::progress_display * show_progress = new boost::progress_display(allEntries);
  for(Long64_t entry = 0; entry < allEntries; ++entry)  
  {
    ++(*show_progress);
    HT=0.0;
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    NJetAK4=branchAK4Jet->GetEntriesFast();
    for(int i = 0; i < NJetAK4; ++i)
    {
      Jet* jet = (Jet*) branchAK4Jet->At(i);
      HT+=jet->PT;
    }
    if(HT>400.0)
      PassesHTCut++;
  }
  printf("HT:\n%E\n%d rel. uncert. %E\n%lld\n\n",1.0*PassesHTCut/allEntries,PassesHTCut,1.0/sqrt(PassesHTCut),allEntries);
  return 1.0*PassesHTCut/allEntries;
}

void Analysis::RunBigAnalysis()
{
  ExRootTreeReader *treeReader;
  double integratedLuminosity = 137000.0; //in pb^-1
  double totalZcrosssection=5.4100e4;
  double totalWcrosssection=1.77300e5;
  size_t numZEvents=integratedLuminosity*totalZcrosssection*FindNormalisation("/media/gerrit/Files/DelphesEvents/ZNorm-Cluster/ROOTFILES.txt")/FindNormalisation("/media/gerrit/Files/DelphesEvents/Z-Cluster/ROOTFILES.txt");
  size_t numWEvents=integratedLuminosity*totalWcrosssection*FindNormalisation("/media/gerrit/Files/DelphesEvents/WNorm-Cluster/ROOTFILES.txt")/FindNormalisation("/media/gerrit/Files/DelphesEvents/W-Cluster/ROOTFILES.txt");
  vector<PassedEvent> events;
  cout<< numZEvents<<endl << numWEvents<<endl;
  
  TChain chain("Delphes");
  std::ifstream fileZ("/media/gerrit/Files/DelphesEvents/Z-Cluster/ROOTFILES.txt");
  if (fileZ.is_open()) {
    std::string line;
    while (std::getline(fileZ, line)) {
      chain.Add(line.c_str());
    }
    fileZ.close();
  }
  treeReader = new ExRootTreeReader(&chain);
  Analysis::AnalyseEvents(treeReader,events,numZEvents,0);

  chain.Reset();

  std::ifstream fileW("/media/gerrit/Files/DelphesEvents/W-Cluster/ROOTFILES.txt");
  if (fileW.is_open()) {
    std::string line;
    while (std::getline(fileW, line)) {
      chain.Add(line.c_str());
    }
    fileW.close();
  }
  treeReader = new ExRootTreeReader(&chain);
  Analysis::AnalyseEvents(treeReader,events,numWEvents,1);
  RootIO::SaveEvents(events);
}
