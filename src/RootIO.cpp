#include "RootIO.h"

void RootIO::SaveEvents(const char *filename, vector<PassedEvent> events, int verbose)
{
  if(verbose)
    cout << "Writing " << events.size() << " Events " << endl;
  TFile f(filename, "RECREATE");
  TVectorD mj1(events.size());
  TVectorD mj2(events.size());
  TVectorD ptmiss(events.size());
  TVectorD status(events.size());
  for (size_t i = 0; i < events.size(); i++)
  {
    mj1[i] = events[i].mj1;
    mj2[i] = events[i].mj2;
    ptmiss[i] = events[i].ptmiss;
    status[i] = 1.0 * events[i].status;
  }
  mj1.Write("mj1", TObject::kOverwrite);
  mj2.Write("mj2", TObject::kOverwrite);
  ptmiss.Write("ptmiss", TObject::kOverwrite);
  status.Write("status", TObject::kOverwrite);
  f.Close();
}

void RootIO::SaveEvents(const char *filename, vector<PassedEvent> events)
{
  RootIO::SaveEvents(filename, events, 0);
}

void RootIO::ReadEventsNewPassedEventDef(const char * filename,vector<PassedEvent> &events,int verbose)
{
  TChain *myChain = new TChain("Tree");
  myChain->Add(filename);
  TTreeReader myReader(myChain);
  // TFile f(filename.c_str(), "read");
  // TTreeReader myReader("Tree", &f);
  TTreeReaderValue<int> status(myReader, "Status");
  TTreeReaderValue<double> MET(myReader, "MET");
  TTreeReaderValue<double> Mj1(myReader, "Mj1");
  TTreeReaderValue<double> Mj2(myReader, "Mj2");
  int i = 0;
  while (myReader.Next())
  {
    PassedEvent event;
    // cout <<*MET<<endl;
    event.status = *status;
    event.ptmiss = *MET;
    event.mj1 = *Mj1;
    event.mj2 = *Mj2;

    events.push_back(event);
    // break;

    // cout << (*myVector).size()<<endl;
  }
  myReader.Clear();
  // f.Close();
  cout << events.size() << endl;
}

void RootIO::ReadEvents(const char * filename, vector<PassedEvent> &events, int verbose, int &nTotal)
{
  TFile f(filename);
  TVectorD *mj1 = (TVectorD *)f.Get("mj1");
  TVectorD *mj2 = (TVectorD *)f.Get("mj2");
  TVectorD *ptmiss = (TVectorD *)f.Get("ptmiss");
  TVectorD *status = (TVectorD *)f.Get("status");
  TVectorD *rawEvents=(TVectorD *)f.Get("nRawEvents");
  vector<int> diff;
  for (size_t i = 0; i < mj1->GetNoElements(); i++)
  {
    //mj1[0][0].Print();
    events.push_back(PassedEvent(mj1[0][i], mj2[0][i], ptmiss[0][i], (Int_t)status[0][i]));
    while (diff.size() < (Int_t)status[0][i] + 1)
    {
      diff.push_back(0);
    }
    diff[(Int_t)status[0][i]]++;
  }
  if (verbose)
  {
    cout << "Read " << mj1->GetNoElements() << " events" << endl;
    for (size_t i = 0; i < diff.size(); i++)
    {
      cout << "Status " << i << " " << diff[i] << endl;
    }
  }

  f.Close();
  nTotal=(int)rawEvents[0][0];
}

void RootIO::ReadEvents(const char *filename, vector<PassedEvent> &events, int verbose)
{
  TFile f(filename);
  TVectorD *mj1 = (TVectorD *)f.Get("mj1");
  TVectorD *mj2 = (TVectorD *)f.Get("mj2");
  TVectorD *ptmiss = (TVectorD *)f.Get("ptmiss");
  TVectorD *status = (TVectorD *)f.Get("status");
  vector<int> diff;
  for (size_t i = 0; i < mj1->GetNoElements(); i++)
  {
    //mj1[0][0].Print();
    events.push_back(PassedEvent(mj1[0][i], mj2[0][i], ptmiss[0][i], (Int_t)status[0][i]));
    while (diff.size() < (Int_t)status[0][i] + 1)
    {
      diff.push_back(0);
    }
    diff[(Int_t)status[0][i]]++;
  }
  if (verbose)
  {
    cout << "Read " << mj1->GetNoElements() << " events" << endl;
    for (size_t i = 0; i < diff.size(); i++)
    {
      cout << "Status " << i << " " << diff[i] << endl;
    }
  }

  f.Close();
}

void RootIO::ReadEvents(const char *filename, vector<PassedEvent> &events)
{
  RootIO::ReadEvents(filename, events, 0);
}

void RootIO::GetRootFilePath(const char *filename, vector<string> &rootFilePaths, int nFiles)
{
  std::ifstream file(filename);
  if (file.is_open())
  {
    std::string line;
    while (std::getline(file, line))
    {
      rootFilePaths.push_back(line);
      if (--nFiles == 0)
        break;
    }
    file.close();
  }
}

void RootIO::GetRootFiles(const char *filename, TChain &chain, int nFiles)
{
  vector<string> rootFilePaths;
  RootIO::GetRootFilePath(filename, rootFilePaths, nFiles);
  std::reverse(rootFilePaths.begin(), rootFilePaths.end());
  for (size_t i = 0; i < rootFilePaths.size(); i++)
  {
    chain.Add(rootFilePaths[i].c_str());
  }
}

void RootIO::GetRootFiles(const char *filename, TChain &chain)
{
  GetRootFiles(filename, chain, -1);
}
