#include "PassedEvent.h"
//ClassImp(PassedEvent);

PassedEvent::PassedEvent()
{
    //cout << " Events"<<endl;
}
PassedEvent::PassedEvent(Double_t _mj1, Double_t _mj2, Double_t _ptmiss, Int_t _status)
{
    mj1=_mj1;
    mj2=_mj2;
    ptmiss=_ptmiss;
    status=_status;
}
void PassedEvent::Print()
{
    cout << "mj1    : "<< mj1 << "\nmj2    : " <<mj2<<"\nptmiss : " <<ptmiss<<"\nstatus : "<< status<<endl;
}