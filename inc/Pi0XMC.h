#ifndef __Pi0XMC_h__
#define __Pi0XMC_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "PPhysics.h"
#include "TH2.h"
#include "TH3.h"
#include "TRandom3.h"

class	Pi0XMC  : public PPhysics
{
private:
    const Double_t mpi0      = 134.9766;   // MeV

    GH1*	time;
    GH1*        time_all;
    GH1*	time_cut;

    GH1*        theta;
    GH1*        theta_all;
    GH1*        theta_MC;
    GH1*        thetaCM;
    GH1*        thetaCM_all;
    GH1*        thetaCM_MC;

    GH1*        FPD;
    GH1*        FPD_all;

    TH1*	TaggerAccScal;

    string  	config;

    Int_t       evtNum;

protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void    ProcessScalerRead();
    virtual Bool_t  Write();

public:
    Pi0XMC();
    virtual ~Pi0XMC();
    virtual Bool_t  Init();

    void        FillThetaMC(GH1* gHist, GH1* gHistCM, Bool_t TaggerBinning = kFALSE);
    void        FillTheta(const GTreeParticle& tree, GH1* gHist, GH1* gHistCM, Bool_t TaggerBinning = kFALSE);
    void        FillTheta(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, GH1* gHistCM, Bool_t TaggerBinning = kFALSE);

    void        FillFPD(const GTreeParticle& tree, GH1* gHist);
    void        FillFPD(const GTreeParticle& tree, Int_t particle_index, GH1* gHist);

    Double_t    Compute_ThetaCM(Double_t, Double_t, Double_t, Double_t, Double_t);
};
#endif
