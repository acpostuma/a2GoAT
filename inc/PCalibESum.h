#ifndef __PCalibESum_h__
#define __PCalibESum_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "PPhysics.h"

class	PCalibESum  : public PPhysics
{
private:
    GH1*	time;
    GH1*	time_cut;
    GH1*	time_2g;      
    GH1*	time_2g_cut;   
     
    GH1*	IM;
    GH1*	IM_2g;

    GH1*	MM;
    GH1*	MM_2g;

    GH2*    ESumCalib;

    Int_t   CBESum;

protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
    virtual Bool_t    Write();
			
public:
    PCalibESum();
    virtual ~PCalibESum();
    virtual Bool_t  Init();
    void	FillESum(const GTreeParticle& tree, GH2* gHist, Bool_t TaggerBinning = kFALSE);
    void	FillESum(const GTreeParticle& tree, Int_t particle_index, GH2* gHist, Bool_t TaggerBinning = kFALSE);
    void 	FillESum(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH2* gHist, Bool_t TaggerBinning = kFALSE);

};
#endif
