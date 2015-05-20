#ifndef __TotalInclusive_h__
#define __TotalInclusive_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "PPhysics.h"
#include "TH2.h"
#include "TH3.h"
#include "TRandom3.h"

class	TotalInclusive  : public PPhysics
{
private:

    GH1*	time;
    GH1*        time_all;
    GH1*	time_cut;

    GH1*        FPD_TI_hel0;
    GH1*        FPD_TI_hel1;
    GH1*        FPD_TI_CB_hel0;
    GH1*        FPD_TI_CB_hel1;
    GH1*        FPD_TI_TAPS_hel0;
    GH1*        FPD_TI_TAPS_hel1;

    GH1*        CBtrigg;
    GH1*        TAPStrigg;
    GH1*        MIXtrigg;

    TH1*	TaggerAccScal;

    string  	config;

    Int_t       evtNum;

protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void    ProcessScalerRead();
    virtual Bool_t  Write();

public:
    TotalInclusive();
    virtual ~TotalInclusive();
    virtual Bool_t  Init();

    void        FillFPD_TI(const GTreeTrack& tree, GH1* gHist0, GH1* gHist1, 
			   GH1* gCBHist0, GH1* gCBHist1, GH1* gTAPSHist0, GH1* gTAPSHist1);

    void        FillEsum(const GTreeTrack& tree, GH1* gHist0, GH1* gHist1, GH1* gHist2);
};
#endif
