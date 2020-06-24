#ifndef _PAnalyze_h__
#define _PAnalyze_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "PPhysics.h"
#include "TF1.h"
#include "TRandom3.h"

class	PAnalyze  : public PPhysics
{
private:
    TH1*    TaggerAccScal;
    TH1*    LiveTimeScal;

    TH1*    CorrTaggScal;
    TH1*    PolarizeScal;
    TH1*    Helicity;

    TH1*    Tagg_Tm;
    TH1*    Tagg_0;
    TH1*    Tagg_0_R;
    TH1*    Tagg_1;
    TH1*    Tagg_1_R;

    TH1*    Inc_Tm;
    TH2*    Inc_0;
    TH2*    Inc_0_R;
    TH2*    Inc_1;
    TH2*    Inc_1_R;

    TH1*    Pi0_IM_A;
    TH1*    Pi0_IM_E;
    TH1*    Pi0_IM_I;
    TH3*    Pi0_IM_CC;
    TH3*    Pi0_IM_CT;

    TH1*    Pi0_CA;

    TH2*    Pi0_Tm_NE;
    TH2*    Pi0_Tm_NI;
    TH2*    Pi0_Tm_CE;
    TH2*    Pi0_Tm_CI;
    TH2*    Pi0_Tm_WE;
    TH2*    Pi0_Tm_WI;
    TH2*    Pi0_Tm_TE;
    TH2*    Pi0_Tm_TI;

    TH3*    Pi0_OA;
    TH3*    Pi0_OA_R;
    TH3*    Pi0_OA_Cut;
    TH3*    Pi0_OA_Cut_R;

    TH3*    Pi0_MM_NE_0;
    TH3*    Pi0_MM_NE_0_R;
    TH3*    Pi0_MM_NE_1;
    TH3*    Pi0_MM_NE_1_R;
    TH3*    Pi0_Ph_NE_0;
    TH3*    Pi0_Ph_NE_0_R;
    TH3*    Pi0_Ph_NE_1;
    TH3*    Pi0_Ph_NE_1_R;

    TH3*    Pi0_MM_NI_0;
    TH3*    Pi0_MM_NI_0_R;
    TH3*    Pi0_MM_NI_1;
    TH3*    Pi0_MM_NI_1_R;
    TH3*    Pi0_Ph_NI_0;
    TH3*    Pi0_Ph_NI_0_R;
    TH3*    Pi0_Ph_NI_1;
    TH3*    Pi0_Ph_NI_1_R;

    TH3*    Pi0_MM_CE_0;
    TH3*    Pi0_MM_CE_0_R;
    TH3*    Pi0_MM_CE_1;
    TH3*    Pi0_MM_CE_1_R;
    TH3*    Pi0_Ph_CE_0;
    TH3*    Pi0_Ph_CE_0_R;
    TH3*    Pi0_Ph_CE_1;
    TH3*    Pi0_Ph_CE_1_R;

    TH3*    Pi0_MM_CI_0;
    TH3*    Pi0_MM_CI_0_R;
    TH3*    Pi0_MM_CI_1;
    TH3*    Pi0_MM_CI_1_R;
    TH3*    Pi0_Ph_CI_0;
    TH3*    Pi0_Ph_CI_0_R;
    TH3*    Pi0_Ph_CI_1;
    TH3*    Pi0_Ph_CI_1_R;

    TH3*    Pi0_MM_WE_0;
    TH3*    Pi0_MM_WE_0_R;
    TH3*    Pi0_MM_WE_1;
    TH3*    Pi0_MM_WE_1_R;
    TH3*    Pi0_Ph_WE_0;
    TH3*    Pi0_Ph_WE_0_R;
    TH3*    Pi0_Ph_WE_1;
    TH3*    Pi0_Ph_WE_1_R;

    TH3*    Pi0_MM_WI_0;
    TH3*    Pi0_MM_WI_0_R;
    TH3*    Pi0_MM_WI_1;
    TH3*    Pi0_MM_WI_1_R;
    TH3*    Pi0_Ph_WI_0;
    TH3*    Pi0_Ph_WI_0_R;
    TH3*    Pi0_Ph_WI_1;
    TH3*    Pi0_Ph_WI_1_R;

    TH3*    Pi0_MM_TE_0;
    TH3*    Pi0_MM_TE_0_R;
    TH3*    Pi0_MM_TE_1;
    TH3*    Pi0_MM_TE_1_R;
    TH3*    Pi0_Ph_TE_0;
    TH3*    Pi0_Ph_TE_0_R;
    TH3*    Pi0_Ph_TE_1;
    TH3*    Pi0_Ph_TE_1_R;

    TH3*    Pi0_MM_TI_0;
    TH3*    Pi0_MM_TI_0_R;
    TH3*    Pi0_MM_TI_1;
    TH3*    Pi0_MM_TI_1_R;
    TH3*    Pi0_Ph_TI_0;
    TH3*    Pi0_Ph_TI_0_R;
    TH3*    Pi0_Ph_TI_1;
    TH3*    Pi0_Ph_TI_1_R;

    TH3*    Pi0_Re_All;
    TH3*    Pi0_Re_All_R;
    TH3*    Pi0_Re_Det;
    TH3*    Pi0_Re_Det_R;
    TH3*    Pi0_Re_Dif;
    TH3*    Pi0_Re_Dif_R;
    TH3*    Pi0_Re_NoE;
    TH3*    Pi0_Re_NoE_R;

    TH1*    Comp_CA;

    TH2*    Comp_Tm_NE;
    TH2*    Comp_Tm_NI;
    TH2*    Comp_Tm_CE;
    TH2*    Comp_Tm_CI;
    TH2*    Comp_Tm_WE;
    TH2*    Comp_Tm_WI;
    TH2*    Comp_Tm_TE;
    TH2*    Comp_Tm_TI;

    TH3*    Comp_OA;
    TH3*    Comp_OA_R;
    TH3*    Comp_OA_Cut;
    TH3*    Comp_OA_Cut_R;

    TH3*    Comp_MM_NE_0;
    TH3*    Comp_MM_NE_0_R;
    TH3*    Comp_MM_NE_1;
    TH3*    Comp_MM_NE_1_R;
    TH3*    Comp_Ph_NE_0;
    TH3*    Comp_Ph_NE_0_R;
    TH3*    Comp_Ph_NE_1;
    TH3*    Comp_Ph_NE_1_R;

    TH3*    Comp_MM_NI_0;
    TH3*    Comp_MM_NI_0_R;
    TH3*    Comp_MM_NI_1;
    TH3*    Comp_MM_NI_1_R;
    TH3*    Comp_Ph_NI_0;
    TH3*    Comp_Ph_NI_0_R;
    TH3*    Comp_Ph_NI_1;
    TH3*    Comp_Ph_NI_1_R;

    TH3*    Comp_MM_CE_0;
    TH3*    Comp_MM_CE_0_R;
    TH3*    Comp_MM_CE_1;
    TH3*    Comp_MM_CE_1_R;
    TH3*    Comp_Ph_CE_0;
    TH3*    Comp_Ph_CE_0_R;
    TH3*    Comp_Ph_CE_1;
    TH3*    Comp_Ph_CE_1_R;

    TH3*    Comp_MM_CI_0;
    TH3*    Comp_MM_CI_0_R;
    TH3*    Comp_MM_CI_1;
    TH3*    Comp_MM_CI_1_R;
    TH3*    Comp_Ph_CI_0;
    TH3*    Comp_Ph_CI_0_R;
    TH3*    Comp_Ph_CI_1;
    TH3*    Comp_Ph_CI_1_R;

    TH3*    Comp_MM_WE_0;
    TH3*    Comp_MM_WE_0_R;
    TH3*    Comp_MM_WE_1;
    TH3*    Comp_MM_WE_1_R;
    TH3*    Comp_Ph_WE_0;
    TH3*    Comp_Ph_WE_0_R;
    TH3*    Comp_Ph_WE_1;
    TH3*    Comp_Ph_WE_1_R;

    TH3*    Comp_MM_WI_0;
    TH3*    Comp_MM_WI_0_R;
    TH3*    Comp_MM_WI_1;
    TH3*    Comp_MM_WI_1_R;
    TH3*    Comp_Ph_WI_0;
    TH3*    Comp_Ph_WI_0_R;
    TH3*    Comp_Ph_WI_1;
    TH3*    Comp_Ph_WI_1_R;

    TH3*    Comp_MM_TE_0;
    TH3*    Comp_MM_TE_0_R;
    TH3*    Comp_MM_TE_1;
    TH3*    Comp_MM_TE_1_R;
    TH3*    Comp_Ph_TE_0;
    TH3*    Comp_Ph_TE_0_R;
    TH3*    Comp_Ph_TE_1;
    TH3*    Comp_Ph_TE_1_R;

    TH3*    Comp_MM_TI_0;
    TH3*    Comp_MM_TI_0_R;
    TH3*    Comp_MM_TI_1;
    TH3*    Comp_MM_TI_1_R;
    TH3*    Comp_Ph_TI_0;
    TH3*    Comp_Ph_TI_0_R;
    TH3*    Comp_Ph_TI_1;
    TH3*    Comp_Ph_TI_1_R;

    // Do I want these?

    TH2*    Split_OA_E;
    TH3*    MM_CA_OA;
    TH3*    MM_CA_OA_R;

    TH2*    Comp_Tm_N_MM;
    TH2*    Comp_Tm_N_CS;
    TH2*    Comp_Tm_N_MMCS;

    TH3*    Comp_MM_N_C;
    TH3*    Comp_MM_N_C_R;

    TH3*    Comp_CS;
    TH3*    Comp_CS_MM;
    TH3*    Comp_CS_MM_R;
    TH3*    Reco_CS;
    TH3*    Reco_CS_MM;
    TH3*    Reco_CS_MM_R;

    //

    Int_t   verbosity;
    Bool_t  excl_pi0;
    Bool_t  excl_pro;

    Double_t IMCut;
    Double_t MMLoC;
    Double_t MMHiC;
    Double_t OACut;
    Double_t ESCut;

    Bool_t  save_randoms;
    Bool_t   split_search;

    Double_t taps_eff;
    std::vector<Bool_t> ignoreTrack;

    Bool_t  cir_beam;
    Bool_t  lin_beam;
    Double_t beamPol;
    std::vector<Int_t> beamPolTime;
    std::vector<Double_t> beamPolMeas;

    Double_t targPol;
    std::vector<Int_t> targPolTime;
    std::vector<Double_t> targPolMeas;

    Bool_t  firstEvent;

protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
    virtual Bool_t    Write();

public:
    PAnalyze();
    virtual ~PAnalyze();
    virtual Bool_t  Init();
    Bool_t InitVerbosity();
    Bool_t InitExclusivity();
    Bool_t InitInvariantMass();
    Bool_t InitMissingMass();
    Bool_t InitOpeningAngle();
    Bool_t InitEnergySum();
    Bool_t InitSaveRandoms();
    Bool_t InitTAPSEff();
    Bool_t InitBeamPol();
    Bool_t InitTargPol();
    Bool_t InitSplitSearch();
    Double_t TwoBodyAngleToEnergyMin(Double_t eBeam, Double_t mTarg, Double_t mPar1, Double_t mPar2, Double_t tPar1);
    Double_t TwoBodyAngleToEnergyMax(Double_t eBeam, Double_t mTarg, Double_t mPar1, Double_t mPar2, Double_t tPar1);
    Double_t TwoBodyEnergyToAngle(Double_t eBeam, Double_t mTarg, Double_t mPar1, Double_t mPar2, Double_t ePar1);
    Double_t CalcCircBeamPol(Double_t E_e, Double_t P_e, Double_t E_g);

};
#endif
