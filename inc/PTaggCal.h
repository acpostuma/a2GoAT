#ifndef __PTaggCal_h__
#define __PTaggCal_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "PPhysics.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"

class	PTaggCal  : public PPhysics
{
private:
    Int_t nTaggerChannels = 0;

    TH1*	TaggerCurScal;
    TH1*	TaggerSumScal;

    std::vector<TGraph> TaggEnerCalib;
    Int_t numReads = 0;
    Double_t minNMR = 2;
    Double_t maxNMR = 0;
    Double_t preNMR = 2;
    
protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
    virtual Bool_t    Write();
			
public:
    PTaggCal();
    virtual ~PTaggCal();
    virtual Bool_t  Init();

};
#endif
