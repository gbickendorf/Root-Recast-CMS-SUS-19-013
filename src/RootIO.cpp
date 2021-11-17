#include "RootIO.h"

void RootIO::SaveEvents(vector<PassedEvent> events)
{
  cout<<"Writing "<< events.size()<<" Events "<<endl;
  TFile f("Analysis.root","RECREATE");
  TVectorD mj1(events.size());
  TVectorD mj2(events.size());
  TVectorD ptmiss(events.size());
  TVectorD status(events.size());
  for (size_t i = 0; i < events.size(); i++)
  {
    mj1[i]=events[i].mj1;
    mj2[i]=events[i].mj2;
    ptmiss[i]=events[i].ptmiss;
    status[i]=1.0*events[i].status;
  }
  mj1.Write("mj1",TObject::kOverwrite);
  mj2.Write("mj2",TObject::kOverwrite);
  ptmiss.Write("ptmiss",TObject::kOverwrite);
  status.Write("status",TObject::kOverwrite);
  f.Close();
}

void RootIO::ReadEvents(vector<PassedEvent> &events)
{
  TFile f("Analysis2.root");
  TVectorD *mj1 = (TVectorD*)f.Get("mj1");
  TVectorD *mj2 = (TVectorD*)f.Get("mj2");
  TVectorD *ptmiss = (TVectorD*)f.Get("ptmiss");
  TVectorD *status = (TVectorD*)f.Get("status");
  for (size_t i = 0; i < mj1->GetNoElements(); i++)
  {
    //mj1[0][0].Print();
    events.push_back(PassedEvent(mj1[0][i],mj2[0][i],ptmiss[0][i],(Int_t)status[0][i]));
  }
  cout << "Read "<< mj1->GetNoElements()<<" events"<<endl;
  f.Close();
  
}
