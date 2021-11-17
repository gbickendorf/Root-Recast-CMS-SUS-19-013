#pragma once
#include <TROOT.h>
#include "TObject.h"
#include <iostream>
using namespace std;

class PassedEvent : public TObject {
  public:  
    PassedEvent();
    PassedEvent(Double_t _mj1, Double_t _mj2, Double_t _ptmiss, Int_t _status);
//    ~PassedEvent();
    Double_t mj1;
    Double_t mj2;
    Double_t ptmiss;
    Int_t status;//0:Z 1:W 2:tt 3:t 
    void Print();
    //ClassDef(PassedEvent,1)
};