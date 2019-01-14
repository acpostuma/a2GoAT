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
#include "TLine.h"
#include "TLegend.h"

class	PTaggCal  : public PPhysics
{
private:
    Int_t nTaggerChannels = 0;

    TH1*	TaggerCurScal;
    TH1*	TaggerSumScal;

    std::vector<TGraph> TaggEnerCalib;
    std::vector<TLine> Intersections;

    TString scanName = "Tagger Field Scan";

    Int_t nScalerReads = 0;

    Int_t numReads = 0;
    Int_t readNum = 0;
    Double_t NMRdif = 0.0001;

    Double_t minNMR = 2;
    Double_t maxNMR = 0;
    Double_t NMRmin = 0;
    Double_t NMRmax = 2;

    Int_t preMaxBin;
    Double_t preNMR = 2;
    Double_t preHiPer;
    Double_t preLoPer;
    Int_t newMaxBin;
    Double_t newNMR = 2;
    Double_t newHiPer;
    Double_t newLoPer;

protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
    virtual Bool_t    Write();
			
public:
    PTaggCal();
    virtual ~PTaggCal();
    virtual Bool_t  Init();
    Bool_t InitScanName();
    Bool_t InitNMRScan();
    Bool_t InitNMRLimits();

};
#endif
