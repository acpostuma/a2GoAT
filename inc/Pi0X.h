#ifndef __Pi0X_h__
#define __Pi0X_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "PPhysics.h"
#include "TH2.h"
#include "TH3.h"
#include "TRandom3.h"

class	Pi0X  : public PPhysics
{
private:

    GH1*	time;
    GH1*        time_all;
    GH1*	time_cut;

    GH1*        egamma;
    GH1*        egamma_all;

    GH1*        phi;
    GH2*        thetaphi;
    GH1*        phi_all;
    GH2*        thetaphi_all;
    GH1*        Cphi;
    GH2*        Cthetaphi;
    GH1*        Cphi_all;
    GH2*        Cthetaphi_all;

    GH1*        theta;
    GH1*        theta_hel0;
    GH1*        theta_hel1;
    GH1*        theta_hel0CM;
    GH1*        theta_hel1CM;
    GH1*        thetaCB_hel0;
    GH1*        thetaCB_hel1;
    GH1*        thetaCB_hel0CM;
    GH1*        thetaCB_hel1CM;
    GH1*        thetaTAPS_hel0;
    GH1*        thetaTAPS_hel1;
    GH1*        thetaTAPS_hel0CM;
    GH1*        thetaTAPS_hel1CM;
    GH1*        theta_all;
    GH1*        theta_all_hel0;
    GH1*        theta_all_hel1;
    GH1*        Ctheta;
    GH1*        Ctheta_hel0;
    GH1*        Ctheta_hel1;
    GH1*        Ctheta_hel0CM;
    GH1*        Ctheta_hel1CM;
    GH1*        CthetaCB_hel0;
    GH1*        CthetaCB_hel1;
    GH1*        CthetaCB_hel0CM;
    GH1*        CthetaCB_hel1CM;
    GH1*        CthetaTAPS_hel0;
    GH1*        CthetaTAPS_hel1;
    GH1*        CthetaTAPS_hel0CM;
    GH1*        CthetaTAPS_hel1CM;
    GH1*        Ctheta_all;
    GH1*        Ctheta_all_hel0;
    GH1*        Ctheta_all_hel1;

    GH1*        FPD;
    GH1*        FPD_hel0;
    GH1*        FPD_hel1;
    GH1*        FPD_CB_hel0;
    GH1*        FPD_CB_hel1;
    GH1*        FPD_TAPS_hel0;
    GH1*        FPD_TAPS_hel1;
    GH1*        FPD_all;
    GH1*        FPD_all_hel0;
    GH1*        FPD_all_hel1;

    GH1*	IM;
    GH1*	IM_gg;
    GH1*	IM_ggg;
    GH1*	CIM_gg;

    GH1*	MM;
    GH1*	MM_gg;
    GH1*	MM_ggg;
    GH1*	CMM_gg;

    GH1*        helicity;
    GH1*        helicityZE;
    GH1*        helerrors;
    GH1*        errcode;

    TH1*	TaggerAccScal;

    string  	config;

    Int_t       evtNum;

protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void    ProcessScalerRead();
    virtual Bool_t  Write();

public:
    Pi0X();
    virtual ~Pi0X();
    virtual Bool_t  Init();

    void        FillPhotonEnergy(const GTreeParticle& tree, GH1* gHist);
    void        FillPhotonEnergy(const GTreeParticle& tree, Int_t particle_index, GH1* gHist);

    void        FillHelicity(GH1* gHist, GH1* gHist2, GH1* gHist3, GH1* gHist4);

    void        FillAngularDist(const GTreeParticle& tree, GH1* hHist, GH2* ghHist);
    void        FillAngularDist(const GTreeParticle& tree, Int_t particle_index, GH1* hHist, GH2* ghHist);

    void        FillTheta(const GTreeParticle& tree, GH1* gHist, Bool_t TaggerBinning = kFALSE);
    void        FillTheta(const GTreeParticle& tree, GH1* gHist0, GH1* gHist1, Bool_t TaggerBinning = kFALSE);
    void        FillTheta(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning = kFALSE);
    void        FillTheta(const GTreeParticle& tree, Int_t particle_index, GH1* gHist0, GH1* gHist1, Bool_t TaggerBinning = kFALSE);
    //    void        FillTheta(const GTreeParticle& tree, Int_t particle_index, GH1* gHist0, GH1* gCBHist0, GH1* gTAPSHist0, 
    //			  GH1* gHist1, GH1* gCBHist1, GH1* gTAPSHist1, Bool_t TaggerBinning = kFALSE);
    void        FillTheta(const GTreeParticle& tree, Int_t particle_index, 
			  GH1* gHist0, GH1* gCBHist0, GH1* gTAPSHist0, GH1* gHist0CM, GH1* gCBHist0CM, GH1* gTAPSHist0CM,
			  GH1* gHist1, GH1* gCBHist1, GH1* gTAPSHist1, GH1* gHist1CM, GH1* gCBHist1CM, GH1* gTAPSHist1CM,
			  Bool_t TaggerBinning = kFALSE);
    void        FillFPD(const GTreeParticle& tree, GH1* gHist);
    void        FillFPD(const GTreeParticle& tree, GH1* gHist0, GH1* gHist1);
    void        FillFPD(const GTreeParticle& tree, Int_t particle_index, GH1* gHist);
    void        FillFPD(const GTreeParticle& tree, Int_t particle_index, GH1* gHist0, GH1* gCBHist0, GH1* gTAPSHist0, GH1* gHist1, GH1* gCBHist1, GH1* gTAPSHist1);

    Double_t    Compute_ThetaCM(Double_t, Double_t, Double_t, Double_t, Double_t);
};
#endif
