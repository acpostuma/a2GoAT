#ifndef __PTaggBack_h__
#define __PTaggBack_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "PPhysics.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"

class	PTaggBack  : public PPhysics
{
private:
    Int_t nTaggerChannels = 0;

    TH1*	TaggerAccScal;
    TH1*	TaggerPreScal;
    TH1*	TaggerCurScal;
    TH1*	TaggerFirScal;
    TH1*	LiveTimeScal;
    Bool_t  FreeScalers;

    TGraphErrors TaggBackground;
    Double_t clock_pre = 0;
    Int_t nScalIter = 0;
    Int_t nScalRead = 1;
    Int_t back_bin_lo = 1;
    Int_t back_bin_hi = 328;
    
protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
    virtual Bool_t    Write();
			
public:
    PTaggBack();
    virtual ~PTaggBack();
    virtual Bool_t  Init();
    Bool_t InitFreeScalers();
    Bool_t InitBackgroundScalers();

};
#endif
