#ifndef _PAnalyze_h__
#define _PAnalyze_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "PPhysics.h"
#include "TF1.h"

class	PAnalyze  : public PPhysics
{
private:
    TH1*    TaggerAccScal;
    TH1*    LiveTimeScal;

    TH1*    CorrTaggScal;
    TH1*    PolarizeScal;
    TH1*    Helicity;

    TH1*    TaggTime;
    TH1*    TaggHel0;
    TH1*    TaggHel0_R;
    TH1*    TaggHel1;
    TH1*    TaggHel1_R;

    TH1*    IncTime;
    TH1*    IncHel0;
    TH1*    IncHel0_R;
    TH1*    IncHel1;
    TH1*    IncHel1_R;
    TH2*    IncHits;
    TH2*    IncHits_R;

    TH1*    Pi0IM;
    TH1*    Pi0IM_XX;
    TH1*    Pi0IM_NN;
    TH3*    Pi0IM_CBSum;
    TH3*    Pi0IM_CB;
    TH3*    Pi0IM_TAPS;
    TH3*    Pi0IM_CBTAPS;
    TH1*    Coplanarity;
    TH2*    Split_OA_E;
    TH3*    MM_CA_OA;
    TH3*    MM_CA_OA_R;

    TH2*    Pi0_Dt;
    TH2*    Pi0_Dt_NN;

    TH3*    Pi0_IM_MM;
    TH3*    Pi0_IM_MM_R;
    TH3*    Pi0_IM_MM_NN;
    TH3*    Pi0_IM_MM_NN_R;

    TH3*    Pi0_MM;
    TH3*    Pi0_MM_R;
    TH3*    Pi0_MM_IM;
    TH3*    Pi0_MM_IM_R;

    TH3*    Pi0_MM_NN;
    TH3*    Pi0_MM_NN_R;
    TH3*    Pi0_MM_NC;
    TH3*    Pi0_MM_NC_R;
    TH3*    Pi0_MM_CC;
    TH3*    Pi0_MM_CC_R;

    TH3*    Pi0_MM_NNX;
    TH3*    Pi0_MM_NNX_R;
    TH3*    Pi0_MM_NCX;
    TH3*    Pi0_MM_NCX_R;
    TH3*    Pi0_MM_CCX;
    TH3*    Pi0_MM_CCX_R;

    TH3*    Pi0_Re_All;
    TH3*    Pi0_Re_All_R;
    TH3*    Pi0_Re_Det;
    TH3*    Pi0_Re_Det_R;
    TH3*    Pi0_Re_NoE;
    TH3*    Pi0_Re_NoE_R;

    TH2*    Comp_Dt_N_MM;
    TH2*    Comp_Dt_N_CS;
    TH2*    Comp_Dt_N_MMCS;

    TH2*    Comp_Dt_N;
    TH2*    Comp_Dt_C;
    TH2*    Comp_Dt_NN;
    TH2*    Comp_Dt_NC;
    TH2*    Comp_Dt_NT;
    TH2*    Comp_Dt_NW;
    TH2*    Comp_Dt_NNX;
    TH2*    Comp_Dt_NCX;
    TH2*    Comp_Dt_NTX;
    TH2*    Comp_Dt_NWX;

    TH3*    Comp_MM_N_C;
    TH3*    Comp_MM_N_C_R;

    TH3*    Comp_CS;
    TH3*    Comp_CS_MM;
    TH3*    Comp_CS_MM_R;
    TH3*    Reco_CS;
    TH3*    Reco_CS_MM;
    TH3*    Reco_CS_MM_R;

    TH3*    Comp_MM_N_0;
    TH3*    Comp_MM_N_0_R;
    TH3*    Comp_MM_N_1;
    TH3*    Comp_MM_N_1_R;
    TH3*    Comp_Ph_N_0;
    TH3*    Comp_Ph_N_0_R;
    TH3*    Comp_Ph_N_1;
    TH3*    Comp_Ph_N_1_R;

    TH3*    Comp_MM_C_0;
    TH3*    Comp_MM_C_0_R;
    TH3*    Comp_MM_C_1;
    TH3*    Comp_MM_C_1_R;
    TH3*    Comp_Ph_C_0;
    TH3*    Comp_Ph_C_0_R;
    TH3*    Comp_Ph_C_1;
    TH3*    Comp_Ph_C_1_R;

    TH3*    Comp_MM_NN_0;
    TH3*    Comp_MM_NN_0_R;
    TH3*    Comp_MM_NN_1;
    TH3*    Comp_MM_NN_1_R;
    TH3*    Comp_Ph_NN_0;
    TH3*    Comp_Ph_NN_0_R;
    TH3*    Comp_Ph_NN_1;
    TH3*    Comp_Ph_NN_1_R;

    TH3*    Comp_MM_NC_0;
    TH3*    Comp_MM_NC_0_R;
    TH3*    Comp_MM_NC_1;
    TH3*    Comp_MM_NC_1_R;
    TH3*    Comp_Ph_NC_0;
    TH3*    Comp_Ph_NC_0_R;
    TH3*    Comp_Ph_NC_1;
    TH3*    Comp_Ph_NC_1_R;

    TH3*    Comp_MM_NT_0;
    TH3*    Comp_MM_NT_0_R;
    TH3*    Comp_MM_NT_1;
    TH3*    Comp_MM_NT_1_R;
    TH3*    Comp_Ph_NT_0;
    TH3*    Comp_Ph_NT_0_R;
    TH3*    Comp_Ph_NT_1;
    TH3*    Comp_Ph_NT_1_R;

    TH3*    Comp_MM_NW_0;
    TH3*    Comp_MM_NW_0_R;
    TH3*    Comp_MM_NW_1;
    TH3*    Comp_MM_NW_1_R;
    TH3*    Comp_Ph_NW_0;
    TH3*    Comp_Ph_NW_0_R;
    TH3*    Comp_Ph_NW_1;
    TH3*    Comp_Ph_NW_1_R;

    TH3*    Comp_MM_NNX_0;
    TH3*    Comp_MM_NNX_0_R;
    TH3*    Comp_MM_NNX_1;
    TH3*    Comp_MM_NNX_1_R;
    TH3*    Comp_Ph_NNX_0;
    TH3*    Comp_Ph_NNX_0_R;
    TH3*    Comp_Ph_NNX_1;
    TH3*    Comp_Ph_NNX_1_R;

    TH3*    Comp_MM_NCX_0;
    TH3*    Comp_MM_NCX_0_R;
    TH3*    Comp_MM_NCX_1;
    TH3*    Comp_MM_NCX_1_R;
    TH3*    Comp_Ph_NCX_0;
    TH3*    Comp_Ph_NCX_0_R;
    TH3*    Comp_Ph_NCX_1;
    TH3*    Comp_Ph_NCX_1_R;

    TH3*    Comp_MM_NTX_0;
    TH3*    Comp_MM_NTX_0_R;
    TH3*    Comp_MM_NTX_1;
    TH3*    Comp_MM_NTX_1_R;
    TH3*    Comp_Ph_NTX_0;
    TH3*    Comp_Ph_NTX_0_R;
    TH3*    Comp_Ph_NTX_1;
    TH3*    Comp_Ph_NTX_1_R;

    TH3*    Comp_MM_NWX_0;
    TH3*    Comp_MM_NWX_0_R;
    TH3*    Comp_MM_NWX_1;
    TH3*    Comp_MM_NWX_1_R;
    TH3*    Comp_Ph_NWX_0;
    TH3*    Comp_Ph_NWX_0_R;
    TH3*    Comp_Ph_NWX_1;
    TH3*    Comp_Ph_NWX_1_R;

    Int_t   verbosity;
    Bool_t  excl_pi0;
    Bool_t  excl_pro;

    Double_t IMCut;
    Double_t MMLoC;
    Double_t MMHiC;
    Double_t OACut;
    Double_t ESCut;

    Bool_t  save_randoms;

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
    Bool_t InitBeamPol();
    Bool_t InitTargPol();
    Double_t TwoBodyAngleToEnergyMin(Double_t eBeam, Double_t mTarg, Double_t mPar1, Double_t mPar2, Double_t tPar1);
    Double_t TwoBodyAngleToEnergyMax(Double_t eBeam, Double_t mTarg, Double_t mPar1, Double_t mPar2, Double_t tPar1);
    Double_t TwoBodyEnergyToAngle(Double_t eBeam, Double_t mTarg, Double_t mPar1, Double_t mPar2, Double_t ePar1);
    Double_t CalcCircBeamPol(Double_t E_e, Double_t P_e, Double_t E_g);

};
#endif
