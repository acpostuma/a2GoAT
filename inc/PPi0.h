#ifndef __PPi0_h__
#define __PPi0_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "PPhysics.h"
#include "TH2.h"
#include "TH3.h"

class	PPi0  : public PPhysics
{
private:
    GH1*	time;
    GH1*	time_cut;

    GH1*	IM;

    GH1*	IM_gg;
    GH1*	IM_ggg;

    GH1*	IM_rr;
    GH1*	IM_rrr;

    GH1*	IM_gr;
    GH1*	IM_ggr;
    GH1*	IM_grr;

    TH3*	IM_all;

    GH1*	MM;

    GH1*	MM_gg;
    GH1*	MM_ggg;

    GH1*	MM_rr;
    GH1*	MM_rrr;

    GH1*	MM_gr;
    GH1*	MM_ggr;
    GH1*	MM_grr;

    GH2*    MM_IM;

    TH1*	TaggerAccScal;

    string  	config;
protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void    ProcessScalerRead();
    virtual Bool_t  Write();

public:
    PPi0();
    virtual ~PPi0();
    virtual Bool_t  Init();

    void	FillMassMissingMass(const GTreeParticle& tree, GH2* gHist, Bool_t TaggerBinning = kFALSE);
    void	FillMassMissingMass(const GTreeParticle& tree, Int_t particle_index, GH2* gHist, Bool_t TaggerBinning = kFALSE);
    void 	FillMassMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH2* gHist, Bool_t TaggerBinning = kFALSE);

};
#endif
