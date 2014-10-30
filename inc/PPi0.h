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
    TH1*	time;
    TH1*	time_cut;
     
    TH1*	IM;

    TH1*	IM_gg;
    TH1*	IM_ggg;

    TH1*	IM_rr;
    TH1*	IM_rrr;

    TH1*	IM_gr;
    TH1*	IM_ggr;
    TH1*	IM_grr;

    TH3*	IM_all;

    GH1*	MM;

    GH1*	MM_gg;
    GH1*	MM_ggg;

    GH1*	MM_rr;
    GH1*	MM_rrr;

    GH1*	MM_gr;
    GH1*	MM_ggr;
    GH1*	MM_grr;

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
    virtual Bool_t  Init(const char* configfile);

};
#endif
