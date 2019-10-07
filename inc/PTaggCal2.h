#ifndef __PTaggCal2_h__
#define __PTaggCal2_h__

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
#include "TLine.h"
#include "TLegend.h"

class	PTaggCal2  : public PPhysics
{
private:
    Int_t nTaggerChannels = 0;

    TH1*	TaggerCurScal;
    TH1*	TaggerSumScal;

    std::vector<TH1*> TaggerScalers;
    std::vector<Float_t> NMRReads;
    std::vector<Float_t> NMRReadsAll;

    TGraph nReadsStableNMR;

    std::vector<TGraphErrors> TaggEnerCalib;
    std::vector<TLine> Intersections;

    TString scanName = "Tagger Field Scan";

    Int_t numReads = 0;
    Int_t readNum = 0;
    Int_t readOff = 0;
    Float_t NMRdif = 0.0001;

    Float_t thisMinNMR = 2;
    Float_t thisMaxNMR = 0;

    Float_t scanMinNMR = 2;
    Float_t scanMaxNMR = 0;

    Float_t NMRmin = 0;
    Float_t NMRmax = 2;

    Int_t preMaxBin;
    Float_t preNMR = 2;
    Double_t preHiPer;
    Double_t preLoPer;
    Int_t newMaxBin;
    Float_t newNMR = 2;
    Double_t newHiPer;
    Double_t newLoPer;

    Double_t preCur = 0;

protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
    virtual Bool_t    Write();
			
public:
    PTaggCal2();
    virtual ~PTaggCal2();
    virtual Bool_t  Init();
    Bool_t InitScanName();
    Bool_t InitNMRScan();
    Bool_t InitNMRLimits();

};
#endif
