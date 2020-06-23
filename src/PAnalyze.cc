#include "PAnalyze.h"

PAnalyze::PAnalyze()
{
    CorrTaggScal = new TH1D("CorrTaggScal", "Livetime corrected tagger scalers;Tagger Channel;LT*Scaler", 352, 0, 352);
    PolarizeScal = new TH1D("PolarizeScal", "Livetime and polarization corrected tagger scalers;Tagger Channel;LT*Pt*Pg*Scaler", 352, 0, 352);

    Helicity = new TH1D("Helicity", "Helicity", 3, 0, 3);

    TaggTime = new TH1D("TaggTime", "Tagger Time;t_{#gamma} (ns)", 1400, -700, 700);
    TaggHel0 = new TH1D("TaggHel0", "Tagger Hits;Tagger Channel", 352, 0, 352);
    TaggHel0_R = new TH1D("TaggHel0_R", "Tagger Hits;Tagger Channel", 352, 0, 352);
    TaggHel1 = new TH1D("TaggHel1", "Tagger Hits;Tagger Channel", 352, 0, 352);
    TaggHel1_R = new TH1D("TaggHel1_R", "Tagger Hits;Tagger Channel", 352, 0, 352);

    IncTime = new TH1D("IncTime", "Tagger - CB Time;t_{Tagger}-t_{CB} (ns)", 1400, -700, 700);
    IncHel0 = new TH1D("IncHel0", "Total Inclusive Hits;Tagger Channel", 352, 0, 352);
    IncHel0_R = new TH1D("IncHel0_R", "Total Inclusive Hits;Tagger Channel", 352, 0, 352);
    IncHel1 = new TH1D("IncHel1", "Total Inclusive Hits;Tagger Channel", 352, 0, 352);
    IncHel1_R = new TH1D("IncHel1_R", "Total Inclusive Hits;Tagger Channel", 352, 0, 352);
    IncHits = new TH2D("IncHits", "Total Inclusive Hits;Tagger Channel;Number of Tracks", 352, 0, 352, 10, 0, 10);
    IncHits_R = new TH2D("IncHits_R", "Total Inclusive Hits;Tagger Channel;Number of Tracks", 352, 0, 352, 10, 0, 10);

    Pi0IM   = new TH1D("Pi0IM", "Neutral Pion Invariant Mass;m_{#gamma#gamma} (MeV)", 400, 0, 400);
    Pi0IM_XX   = new TH1D("Pi0IM_XX", "Neutral Pion Invariant Mass;m_{#gamma#gamma} (MeV)", 400, 0, 400);
    Pi0IM_NN   = new TH1D("Pi0IM_NN", "Neutral Pion Invariant Mass;m_{#gamma#gamma} (MeV)", 400, 0, 400);
    Pi0IM_CBSum = new TH3D("Pi0IM_CBSum", "CB Invariant Mass;E_{#gamma1}+E_{#gamma2} (MeV);#sqrt{2*E_{#gamma1}*E_{#gamma2}} (MeV);m_{#gamma#gamma} (MeV)", 60, 150, 450, 50, 100, 300, 60, 115, 155);
    Pi0IM_CB = new TH3D("Pi0IM_CB", "CB Invariant Mass;E_{#gamma_1};E_{#gamma_2};m_{#gamma#gamma} (MeV)", 80, 0, 800, 80, 0, 800, 200, 0, 200);
    Pi0IM_TAPS = new TH3D("Pi0IM_TAPS", "TAPS Invariant Mass;E_{#gamma_1};E_{#gamma_2};m_{#gamma#gamma} (MeV)", 80, 0, 800, 80, 0, 800, 200, 0, 200);
    Pi0IM_CBTAPS = new TH3D("Pi0IM_CBTAPS", "CB/TAPS Invariant Mass;E_{#gamma_{CB}};E_{#gamma_{TAPS}};m_{#gamma#gamma} (MeV)", 80, 0, 800, 80, 0, 800, 200, 0, 200);
    Coplanarity = new TH1D("Coplanarity","Coplanarity Angle;#phi_{0}-#phi_{1} (deg)",360,0,360);
    Split_OA_E = new TH2D("Split_OA_E", "Split OA vs Energy Ratio;E_{split}/E_{#gamma};#theta_{OA} (deg)", 200, 0, 1, 180, 0, 180);
    MM_CA_OA = new TH3D("MM_CA_OA", "Kinematic Cuts;m_{miss}-m_{targ} (MeV);#phi_{0}-#phi_{1} (deg);#theta_{OA} (deg)", 80, -80, 120, 72, 0, 360, 36, 0, 180);
    MM_CA_OA_R = new TH3D("MM_CA_OA_R", "Kinematic Cuts;m_{miss}-m_{targ} (MeV);#phi_{0}-#phi_{1} (deg);#theta_{OA} (deg)", 80, -80, 120, 72, 0, 360, 36, 0, 180);

    Pi0_Dt = new TH2D("Pi0_Dt", "Pi0 Time;Tagger Channel;t_{#gamma}-t_{#pi^{0}} (ns)", 352, 0, 352, 1400, -700, 700);
    Pi0_Dt_NN = new TH2D("Pi0_Dt_NN", "Pi0 Time;Tagger Channel;t_{#gamma}-t_{#pi^{0}} (ns)", 352, 0, 352, 1400, -700, 700);

    Pi0_IM_MM = new TH3D("Pi0_IM_MM", "Pi0 Missing Mass;Tagger Channel;m_{#gamma#gamma} (MeV);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 100, 0, 200, 80, -80, 120);
    Pi0_IM_MM_R = new TH3D("Pi0_IM_MM_R", "Pi0 Missing Mass;Tagger Channel;m_{#gamma#gamma} (MeV);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 100, 0, 200, 80, -80, 120);
    Pi0_IM_MM_NN = new TH3D("Pi0_IM_MM_NN", "Pi0 Missing Mass;Tagger Channel;m_{#gamma#gamma} (MeV);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 100, 0, 200, 80, -80, 120);
    Pi0_IM_MM_NN_R = new TH3D("Pi0_IM_MM_NN_R", "Pi0 Missing Mass;Tagger Channel;m_{#gamma#gamma} (MeV);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 100, 0, 200, 80, -80, 120);

    Pi0_MM = new TH3D("Pi0_MM", "Pi0 Missing Mass;Tagger Channel;#theta^{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_R = new TH3D("Pi0_MM_R", "Pi0 Missing Mass;Tagger Channel;#theta^{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_IM = new TH3D("Pi0_MM_IM", "Pi0 Missing Mass;Tagger Channel;#theta^{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_IM_R = new TH3D("Pi0_MM_IM_R", "Pi0 Missing Mass;Tagger Channel;#theta^{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);

    Pi0_MM_NN = new TH3D("Pi0_MM_NN", "Pi0 Missing Mass;Tagger Channel;#theta^{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_NN_R = new TH3D("Pi0_MM_NN_R", "Pi0 Missing Mass;Tagger Channel;#theta^{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_NC = new TH3D("Pi0_MM_NC", "Pi0 Missing Mass;Tagger Channel;#theta^{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_NC_R = new TH3D("Pi0_MM_NC_R", "Pi0 Missing Mass;Tagger Channel;#theta^{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_CC = new TH3D("Pi0_MM_CC", "Pi0 Missing Mass;Tagger Channel;#theta^{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_CC_R = new TH3D("Pi0_MM_CC_R", "Pi0 Missing Mass;Tagger Channel;#theta^{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);

    Pi0_MM_NNX = new TH3D("Pi0_MM_NNX", "Pi0 Missing Mass;Tagger Channel;#theta^{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_NNX_R = new TH3D("Pi0_MM_NNX_R", "Pi0 Missing Mass;Tagger Channel;#theta^{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_NCX = new TH3D("Pi0_MM_NCX", "Pi0 Missing Mass;Tagger Channel;#theta^{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_NCX_R = new TH3D("Pi0_MM_NCX_R", "Pi0 Missing Mass;Tagger Channel;#theta^{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_CCX = new TH3D("Pi0_MM_CCX", "Pi0 Missing Mass;Tagger Channel;#theta^{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_CCX_R = new TH3D("Pi0_MM_CCX_R", "Pi0 Missing Mass;Tagger Channel;#theta^{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);

    Pi0_Re_All = new TH3D("Pi0_Re_All", "Pi0 Missing Mass;Recoil Energy (MeV);#theta^{p} (deg);m_{miss}-m_{targ} (MeV)", 350, 0, 350, 36, 0, 180, 80, -80, 120);
    Pi0_Re_All_R = new TH3D("Pi0_Re_All_R", "Pi0 Missing Mass;Recoil Energy (MeV);#theta^{p} (deg);m_{miss}-m_{targ} (MeV)", 350, 0, 350, 36, 0, 180, 80, -80, 120);
    Pi0_Re_Det = new TH3D("Pi0_Re_Det", "Pi0 Missing Mass;Recoil Energy (MeV);#theta^{p} (deg);m_{miss}-m_{targ} (MeV)", 350, 0, 350, 36, 0, 180, 80, -80, 120);
    Pi0_Re_Det_R = new TH3D("Pi0_Re_Det_R", "Pi0 Missing Mass;Recoil Energy (MeV);#theta^{p} (deg);m_{miss}-m_{targ} (MeV)", 350, 0, 350, 36, 0, 180, 80, -80, 120);
    Pi0_Re_NoE = new TH3D("Pi0_Re_NoE", "Pi0 Missing Mass;Recoil Energy (MeV);#theta^{p} (deg);m_{miss}-m_{targ} (MeV)", 350, 0, 350, 36, 0, 180, 80, -80, 120);
    Pi0_Re_NoE_R = new TH3D("Pi0_Re_NoE_R", "Pi0 Missing Mass;Recoil Energy (MeV);#theta^{p} (deg);m_{miss}-m_{targ} (MeV)", 350, 0, 350, 36, 0, 180, 80, -80, 120);

    Comp_Dt_N_MM = new TH2D("Comp_Dt_N_MM", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);
    Comp_Dt_N_CS = new TH2D("Comp_Dt_N_CS", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);
    Comp_Dt_N_MMCS = new TH2D("Comp_Dt_N_MMCS", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);

    Comp_Dt_N = new TH2D("Comp_Dt_N", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);
    Comp_Dt_C = new TH2D("Comp_Dt_C", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);
    Comp_Dt_NN = new TH2D("Comp_Dt_NN", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);
    Comp_Dt_NC = new TH2D("Comp_Dt_NC", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);
    Comp_Dt_NT = new TH2D("Comp_Dt_NT", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);
    Comp_Dt_NW = new TH2D("Comp_Dt_NW", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);
    Comp_Dt_NNX = new TH2D("Comp_Dt_NNX", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);
    Comp_Dt_NCX = new TH2D("Comp_Dt_NCX", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);
    Comp_Dt_NTX = new TH2D("Comp_Dt_NTX", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);
    Comp_Dt_NWX = new TH2D("Comp_Dt_NWX", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);

    Comp_MM_N_C = new TH3D("Comp_MM_N_C", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_N_C_R = new TH3D("Comp_MM_N_C_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);

    Comp_CS = new TH3D("Comp_CS", "Compton Cluster Size;E^{#gamma} (MeV);#theta_{#gamma} (deg);Cluster Size", 300, 0, 300, 36, 0, 180, 12, 0, 12);
    Comp_CS_MM = new TH3D("Comp_CS_MM", "Compton Cluster Size Cut by MM;E^{#gamma} (MeV);#theta_{#gamma} (deg);Cluster Size", 300, 0, 300, 36, 0, 180, 12, 0, 12);
    Comp_CS_MM_R = new TH3D("Comp_CS_MM_R", "Compton Cluster Size Cut by MM;E^{#gamma} (MeV);#theta_{#gamma} (deg);Cluster Size", 300, 0, 300, 36, 0, 180, 12, 0, 12);
    Reco_CS = new TH3D("Reco_CS", "Recoil Cluster Size;E^{p} (MeV);#theta_{p} (deg);Cluster Size", 300, 0, 300, 36, 0, 180, 12, 0, 12);
    Reco_CS_MM = new TH3D("Reco_CS_MM", "Recoil Cluster Size Cut by MM;E^{p} (MeV);#theta_{p} (deg);Cluster Size", 300, 0, 300, 36, 0, 180, 12, 0, 12);
    Reco_CS_MM_R = new TH3D("Reco_CS_MM_R", "Recoil Cluster Size Cut by MM;E^{p} (MeV);#theta_{p} (deg);Cluster Size", 300, 0, 300, 36, 0, 180, 12, 0, 12);

    Comp_MM_N_0 = new TH3D("Comp_MM_N_0", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_N_0_R = new TH3D("Comp_MM_N_0_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_N_1 = new TH3D("Comp_MM_N_1", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_N_1_R = new TH3D("Comp_MM_N_1_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_Ph_N_0 = new TH3D("Comp_Ph_N_0", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_N_0_R = new TH3D("Comp_Ph_N_0_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_N_1 = new TH3D("Comp_Ph_N_1", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_N_1_R = new TH3D("Comp_Ph_N_1_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Comp_MM_C_0 = new TH3D("Comp_MM_C_0", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_C_0_R = new TH3D("Comp_MM_C_0_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_C_1 = new TH3D("Comp_MM_C_1", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_C_1_R = new TH3D("Comp_MM_C_1_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_Ph_C_0 = new TH3D("Comp_Ph_C_0", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_C_0_R = new TH3D("Comp_Ph_C_0_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_C_1 = new TH3D("Comp_Ph_C_1", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_C_1_R = new TH3D("Comp_Ph_C_1_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Comp_MM_NN_0 = new TH3D("Comp_MM_NN_0", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NN_0_R = new TH3D("Comp_MM_NN_0_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NN_1 = new TH3D("Comp_MM_NN_1", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NN_1_R = new TH3D("Comp_MM_NN_1_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_Ph_NN_0 = new TH3D("Comp_Ph_NN_0", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NN_0_R = new TH3D("Comp_Ph_NN_0_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NN_1 = new TH3D("Comp_Ph_NN_1", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NN_1_R = new TH3D("Comp_Ph_NN_1_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Comp_MM_NC_0 = new TH3D("Comp_MM_NC_0", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NC_0_R = new TH3D("Comp_MM_NC_0_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NC_1 = new TH3D("Comp_MM_NC_1", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NC_1_R = new TH3D("Comp_MM_NC_1_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_Ph_NC_0 = new TH3D("Comp_Ph_NC_0", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NC_0_R = new TH3D("Comp_Ph_NC_0_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NC_1 = new TH3D("Comp_Ph_NC_1", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NC_1_R = new TH3D("Comp_Ph_NC_1_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Comp_MM_NT_0 = new TH3D("Comp_MM_NT_0", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NT_0_R = new TH3D("Comp_MM_NT_0_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NT_1 = new TH3D("Comp_MM_NT_1", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NT_1_R = new TH3D("Comp_MM_NT_1_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_Ph_NT_0 = new TH3D("Comp_Ph_NT_0", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NT_0_R = new TH3D("Comp_Ph_NT_0_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NT_1 = new TH3D("Comp_Ph_NT_1", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NT_1_R = new TH3D("Comp_Ph_NT_1_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Comp_MM_NW_0 = new TH3D("Comp_MM_NW_0", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NW_0_R = new TH3D("Comp_MM_NW_0_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NW_1 = new TH3D("Comp_MM_NW_1", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NW_1_R = new TH3D("Comp_MM_NW_1_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_Ph_NW_0 = new TH3D("Comp_Ph_NW_0", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NW_0_R = new TH3D("Comp_Ph_NW_0_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NW_1 = new TH3D("Comp_Ph_NW_1", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NW_1_R = new TH3D("Comp_Ph_NW_1_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Comp_MM_NNX_0 = new TH3D("Comp_MM_NNX_0", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NNX_0_R = new TH3D("Comp_MM_NNX_0_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NNX_1 = new TH3D("Comp_MM_NNX_1", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NNX_1_R = new TH3D("Comp_MM_NNX_1_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_Ph_NNX_0 = new TH3D("Comp_Ph_NNX_0", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NNX_0_R = new TH3D("Comp_Ph_NNX_0_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NNX_1 = new TH3D("Comp_Ph_NNX_1", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NNX_1_R = new TH3D("Comp_Ph_NNX_1_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Comp_MM_NCX_0 = new TH3D("Comp_MM_NCX_0", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NCX_0_R = new TH3D("Comp_MM_NCX_0_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NCX_1 = new TH3D("Comp_MM_NCX_1", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NCX_1_R = new TH3D("Comp_MM_NCX_1_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_Ph_NCX_0 = new TH3D("Comp_Ph_NCX_0", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NCX_0_R = new TH3D("Comp_Ph_NCX_0_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NCX_1 = new TH3D("Comp_Ph_NCX_1", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NCX_1_R = new TH3D("Comp_Ph_NCX_1_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Comp_MM_NTX_0 = new TH3D("Comp_MM_NTX_0", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NTX_0_R = new TH3D("Comp_MM_NTX_0_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NTX_1 = new TH3D("Comp_MM_NTX_1", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NTX_1_R = new TH3D("Comp_MM_NTX_1_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_Ph_NTX_0 = new TH3D("Comp_Ph_NTX_0", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NTX_0_R = new TH3D("Comp_Ph_NTX_0_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NTX_1 = new TH3D("Comp_Ph_NTX_1", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NTX_1_R = new TH3D("Comp_Ph_NTX_1_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Comp_MM_NWX_0 = new TH3D("Comp_MM_NWX_0", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NWX_0_R = new TH3D("Comp_MM_NWX_0_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NWX_1 = new TH3D("Comp_MM_NWX_1", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NWX_1_R = new TH3D("Comp_MM_NWX_1_R", "Compton Missing Mass;Tagger Channel;#theta^{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_Ph_NWX_0 = new TH3D("Comp_Ph_NWX_0", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NWX_0_R = new TH3D("Comp_Ph_NWX_0_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NWX_1 = new TH3D("Comp_Ph_NWX_1", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NWX_1_R = new TH3D("Comp_Ph_NWX_1_R", "Compton Phi Distribution;Tagger Channel;#theta^{#gamma} (deg);#phi^{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    verbosity = 0;
    excl_pi0 = false;
    excl_pro = false;

    IMCut = 134.98;
    MMLoC = 800.00;
    MMHiC = 938.27;
    OACut = 180;
    ESCut = 0;

    save_randoms = false;

    cir_beam = false;
    lin_beam = false;
    beamPol = 1;
    targPol = 1;

    firstEvent = true;
}

PAnalyze::~PAnalyze()
{
}

Bool_t	PAnalyze::Init()
{
    cout << "Initialising physics analysis..." << endl;
    cout << "--------------------------------------------------" << endl << endl;

    if(!InitBackgroundCuts()) return kFALSE;
    if(!InitTargetMass()) return kFALSE;
    if(!InitVerbosity()) return kFALSE;
    if(!InitExclusivity()) return kFALSE;
    if(!InitInvariantMass()) return kFALSE;
    if(!InitMissingMass()) return kFALSE;
    if(!InitOpeningAngle()) return kFALSE;
    if(!InitEnergySum()) return kFALSE;
    if(!InitSaveRandoms()) return kFALSE;
    if(!InitBeamPol()) return kFALSE;
    if(!InitTargPol()) return kFALSE;

    if(!PPhysics::Init()) return kFALSE;

    TaggerAccScal = GetScalerHist("TaggerAccScal");
    if(!TaggerAccScal)
    {
        cout << "No tagger scaler histogram available" << endl;
        return kFALSE;
    }
    LiveTimeScal = GetScalerHist("LiveTimeScal");
    if(!LiveTimeScal)
    {
        cout << "No live time histogram available" << endl;
        return kFALSE;
    }

    cout << "--------------------------------------------------" << endl;
    return kTRUE;
}

Bool_t	PAnalyze::Start()
{
    if(!IsAcquFile())
    {
        cout << "ERROR: Input File is not an Acqu file." << endl;
        return kFALSE;
    }
    SetAsGoATFile();

    TraverseValidEvents();

    //GoosyTagger(TaggerAccScal);
    //GoosyVuprom(TaggerAccScal);
    GoosyNewFPDRecabled(TaggerAccScal);

    return kTRUE;
}

void	PAnalyze::ProcessEvent()
{
    if (firstEvent)
    {
        Int_t fileTime = GetSetupParameters()->GetTimeStamp();
        cout << endl << "Time Stamp for this run: " << fileTime << endl;
        Int_t timeDiff = fileTime;

        for (UInt_t i = 0; i<beamPolTime.size(); i++)
        {
            if (TMath::Abs(beamPolTime.at(i) - fileTime) < timeDiff)
            {
                beamPol = beamPolMeas.at(i);
                timeDiff = TMath::Abs(beamPolTime.at(i) - fileTime);
            }
        }
        for (UInt_t i = 0; i<targPolTime.size(); i+=2)
        {
            if ((targPolTime.at(i) < fileTime) && (fileTime < targPolTime.at(i+1)))
            {
                Double_t tau = ((targPolTime.at(i)-targPolTime.at(i+1))/TMath::Log(targPolMeas.at(i+1)/targPolMeas.at(i)));
                targPol = targPolMeas.at(i)*TMath::Exp((targPolTime.at(i)-fileTime)/tau);
                break;
            }
        }

        cout << "E Beam polarization for this run: " << 100*beamPol << "%" << endl;
        cout << "Target polarization for this run: " << 100*targPol << "%" << endl << endl;
        firstEvent = false;
    }

    if (GetEventNumber()%100000 == 0) cout << "Event " << GetEventNumber() << endl;

    if (GetTrigger()->GetEnergySum() < ESCut) return;

    Float_t event_weight = 1;
    if (IsMCFile()) event_weight = GetTrigger()->GetMCWeight();

    Bool_t hel = 0;
    if (cir_beam && !(IsMCFile()))
    {
        if(!(GetTrigger()->HasHelicity())) return;
        hel = GetTrigger()->GetHelicity();
    }

    for (Int_t i=0; i<(GetTrigger()->GetNErrors()); i++)
    {
        if (!(IsMCFile()) && (cir_beam || (GetTrigger()->GetErrorModuleID())[i] >= 0 || (GetTrigger()->GetErrorModuleIndex())[i] >= 0))
        {
            Helicity->Fill(2);
            return;
        }
    }

    Helicity->Fill(hel);

    Int_t tagg_chan, clus_size;
    Double_t comp_time, pi0_time, tagg_time, tagg_energy, time, energy, theta, theta_cm, phi, phi_cm, copl_ang, open_ang, split_energy, m_miss;
    Double_t reco_energy, reco_theta;
    Int_t reco_size;
    TLorentzVector beam, target, compton, compton_cm, pi0, missing, ptot;
    TVector3 lab_to_cm, recoil, split;
    TLorentzVector track_0, track_1, miss_0, miss_1;
    target = GetTarget();

    Bool_t comp_N, comp_C, comp_NN, comp_NC, comp_NT, comp_NW, comp_NNX, comp_NCX, comp_NTX, comp_NWX, comp_X;
    comp_N = comp_C = comp_NN = comp_NC = comp_NT = comp_NW = comp_NNX = comp_NCX = comp_NTX = comp_NWX = comp_X = false;
    Bool_t fill_Ph_0, fill_Ph_1, fill_MM_0, fill_MM_1, cut_CA, cut_OA;
    fill_Ph_0 = fill_Ph_1 = fill_MM_0 = fill_MM_1 = cut_CA = cut_OA = false;
    split_energy = 0;
    Bool_t pi0_NN, pi0_NC, pi0_CC, pi0_NNX, pi0_NCX, pi0_CCX, pi0_X;
    pi0_NN = pi0_NC = pi0_CC = pi0_NNX = pi0_NCX = pi0_CCX = pi0_X = false;

    if (!lin_beam)
    {
        if(hel) fill_MM_1 = true;
        else fill_MM_0 = true;
    }

    //if(GetDecodeDoubles()) GetTagger()->DecodeDoubles();

    Int_t nTagg = GetTagger()->GetNTagged();
    //Int_t nDoub = GetTagger()->GetNDouble();

    Double_t minIM = 0;
    Bool_t IM_CB = false;
    Bool_t IM_TAPS = false;
    Bool_t IM_CBTAPS = false;
    Double_t e0, e1;

    //////////////////////////////////////////////////
    // Initial pi0 stuff
    //////////////////////////////////////////////////
    if (GetTracks()->GetNTracks() == 2 && GetTracks()->GetClusterEnergy(0) != 0 && GetTracks()->GetClusterEnergy(1) != 0)
    {
        pi0 = GetTracks()->GetVector(0) + GetTracks()->GetVector(1);
        minIM = pi0.M();
        pi0_time = 0.5*(GetTracks()->GetTime(0) + GetTracks()->GetTime(1));
        IM_CB = (GetTracks()->HasCB(0) && GetTracks()->HasCB(1));
        IM_TAPS = (GetTracks()->HasTAPS(0) && GetTracks()->HasTAPS(1));
        IM_CBTAPS = ((GetTracks()->HasCB(0) && GetTracks()->HasTAPS(1)) || (GetTracks()->HasTAPS(0) && GetTracks()->HasCB(1)));
        if (GetTracks()->HasTAPS(0) && GetTracks()->HasCB(1))
        {
            e0 = GetTracks()->GetClusterEnergy(1);
            e1 = GetTracks()->GetClusterEnergy(0);
        }
        else
        {
            e0 = GetTracks()->GetClusterEnergy(0);
            e1 = GetTracks()->GetClusterEnergy(1);
        }
        if (GetTracks()->IsNeutral(0) && GetTracks()->IsNeutral(1)) pi0_NN = true;
        else if (GetTracks()->IsNeutral(0) || GetTracks()->IsNeutral(1)) pi0_NC = true;
        else pi0_CC = true;
    }
    else if (GetTracks()->GetNTracks() == 3)
    {
        if (GetTracks()->GetClusterEnergy(0) != 0 && GetTracks()->GetClusterEnergy(1) != 0)
        {
            ptot = GetTracks()->GetVector(0) + GetTracks()->GetVector(1);
            if(TMath::Abs(134.98 - ptot.M()) < TMath::Abs(134.98 - minIM))
            {
                pi0 = GetTracks()->GetVector(0) + GetTracks()->GetVector(1);
                minIM = pi0.M();
                pi0_time = 0.5*(GetTracks()->GetTime(0) + GetTracks()->GetTime(1));
                IM_CB = (GetTracks()->HasCB(0) && GetTracks()->HasCB(1));
                IM_TAPS = (GetTracks()->HasTAPS(0) && GetTracks()->HasTAPS(1));
                IM_CBTAPS = ((GetTracks()->HasCB(0) && GetTracks()->HasTAPS(1)) || (GetTracks()->HasTAPS(0) && GetTracks()->HasCB(1)));
                if (GetTracks()->HasTAPS(0) && GetTracks()->HasCB(1))
                {
                    e0 = GetTracks()->GetClusterEnergy(1);
                    e1 = GetTracks()->GetClusterEnergy(0);
                }
                else
                {
                    e0 = GetTracks()->GetClusterEnergy(0);
                    e1 = GetTracks()->GetClusterEnergy(1);
                }
                recoil = GetTracks()->GetUnitVector(2);
                reco_energy = GetTracks()->GetClusterEnergy(2);
                copl_ang = TMath::Abs(pi0.Phi()*TMath::RadToDeg()-GetTracks()->GetPhi(2));
                cut_CA = (TMath::Abs(180.0-copl_ang) < OACut);
                pi0_NNX = pi0_NCX = pi0_CCX = false;
                if (GetTracks()->IsNeutral(0) && GetTracks()->IsNeutral(1)) pi0_NNX = true;
                else if (GetTracks()->IsNeutral(0) || GetTracks()->IsNeutral(1)) pi0_NCX = true;
                else pi0_CCX = true;
            }
        }
        if (GetTracks()->GetClusterEnergy(0) != 0 && GetTracks()->GetClusterEnergy(2) != 0)
        {
            ptot = GetTracks()->GetVector(0) + GetTracks()->GetVector(2);
            if(TMath::Abs(134.98 - ptot.M()) < TMath::Abs(134.98 - minIM))
            {
                pi0 = GetTracks()->GetVector(0) + GetTracks()->GetVector(2);
                minIM = pi0.M();
                pi0_time = 0.5*(GetTracks()->GetTime(0) + GetTracks()->GetTime(2));
                IM_CB = (GetTracks()->HasCB(0) && GetTracks()->HasCB(2));
                IM_TAPS = (GetTracks()->HasTAPS(0) && GetTracks()->HasTAPS(2));
                IM_CBTAPS = ((GetTracks()->HasCB(0) && GetTracks()->HasTAPS(2)) || (GetTracks()->HasTAPS(0) && GetTracks()->HasCB(2)));
                if (GetTracks()->HasTAPS(0) && GetTracks()->HasCB(2))
                {
                    e0 = GetTracks()->GetClusterEnergy(2);
                    e1 = GetTracks()->GetClusterEnergy(0);
                }
                else
                {
                    e0 = GetTracks()->GetClusterEnergy(0);
                    e1 = GetTracks()->GetClusterEnergy(2);
                }
                recoil = GetTracks()->GetUnitVector(1);
                reco_energy = GetTracks()->GetClusterEnergy(1);
                copl_ang = TMath::Abs(pi0.Phi()*TMath::RadToDeg()-GetTracks()->GetPhi(1));
                cut_CA = (TMath::Abs(180.0-copl_ang) < OACut);
                pi0_NNX = pi0_NCX = pi0_CCX = false;
                if (GetTracks()->IsNeutral(0) && GetTracks()->IsNeutral(2)) pi0_NNX = true;
                else if (GetTracks()->IsNeutral(0) || GetTracks()->IsNeutral(2)) pi0_NCX = true;
                else pi0_CCX = true;
            }
        }
        if (GetTracks()->GetClusterEnergy(1) != 0 && GetTracks()->GetClusterEnergy(2) != 0)
        {
            ptot = GetTracks()->GetVector(1) + GetTracks()->GetVector(2);
            if(TMath::Abs(134.98 - ptot.M()) < TMath::Abs(134.98 - minIM))
            {
                pi0 = GetTracks()->GetVector(1) + GetTracks()->GetVector(2);
                minIM = pi0.M();
                pi0_time = 0.5*(GetTracks()->GetTime(1) + GetTracks()->GetTime(2));
                IM_CB = (GetTracks()->HasCB(1) && GetTracks()->HasCB(2));
                IM_TAPS = (GetTracks()->HasTAPS(1) && GetTracks()->HasTAPS(2));
                IM_CBTAPS = ((GetTracks()->HasCB(1) && GetTracks()->HasTAPS(2)) || (GetTracks()->HasTAPS(1) && GetTracks()->HasCB(2)));
                if (GetTracks()->HasTAPS(1) && GetTracks()->HasCB(2))
                {
                    e0 = GetTracks()->GetClusterEnergy(2);
                    e1 = GetTracks()->GetClusterEnergy(1);
                }
                else
                {
                    e0 = GetTracks()->GetClusterEnergy(1);
                    e1 = GetTracks()->GetClusterEnergy(2);
                }
                recoil = GetTracks()->GetUnitVector(0);
                reco_energy = GetTracks()->GetClusterEnergy(0);
                copl_ang = TMath::Abs(pi0.Phi()*TMath::RadToDeg()-GetTracks()->GetPhi(0));
                cut_CA = (TMath::Abs(180.0-copl_ang) < OACut);
                pi0_NNX = pi0_NCX = pi0_CCX = false;
                if (GetTracks()->IsNeutral(1) && GetTracks()->IsNeutral(2)) pi0_NNX = true;
                else if (GetTracks()->IsNeutral(1) || GetTracks()->IsNeutral(2)) pi0_NCX = true;
                else pi0_CCX = true;
            }
        }
    }

    if (minIM > 0)
    {
        if (TMath::Abs(134.98 - minIM) <= IMCut) pi0_X = pi0_NN || pi0_NC || pi0_CC || pi0_NNX || pi0_NCX || pi0_CCX;

        Pi0IM->Fill(minIM,event_weight);
        if (GetTracks()->GetNTracks() == 2)
        {
            Pi0IM_XX->Fill(minIM,event_weight);
            if (pi0_NN) Pi0IM_NN->Fill(minIM,event_weight);
        }
        if (IM_CB)
        {
            Pi0IM_CBSum->Fill(e0+e1,TMath::Sqrt(2*e0*e1),minIM,event_weight);
            Pi0IM_CB->Fill(e0,e1,minIM,event_weight);
        }
        if (IM_TAPS) Pi0IM_TAPS->Fill(e0,e1,minIM,event_weight);
        if (IM_CBTAPS) Pi0IM_CBTAPS->Fill(e0,e1,minIM,event_weight);
    }

    //////////////////////////////////////////////////
    // Initial compton stuff
    //////////////////////////////////////////////////
    if (!pi0_X && GetTracks()->GetNTracks() == 1 && GetTracks()->GetClusterEnergy(0) != 0)
    {
        compton = GetTracks()->GetVector(0);
        comp_time = GetTracks()->GetTime(0);
        energy = compton.E();
        theta = compton.Theta()*TMath::RadToDeg();
        phi = compton.Phi()*TMath::RadToDeg();
        clus_size = GetTracks()->GetClusterSize(0);
        if (lin_beam)
        {
            if (phi < -90 || (phi > 0 && phi < 90)) fill_MM_0 = true;
            else fill_MM_1 = true;
        }
        if (GetTracks()->IsNeutral(0)) comp_N = true;
        else comp_C = true;
    }
    else if (!pi0_X && GetTracks()->GetNTracks() == 2 && (GetTracks()->GetClusterEnergy(0) != 0 || GetTracks()->GetClusterEnergy(1) != 0))
    {
        copl_ang = TMath::Abs(GetTracks()->GetPhi(0)-GetTracks()->GetPhi(1));
        Coplanarity->Fill(copl_ang,event_weight);

        if (TMath::Abs(180.0-copl_ang) < OACut)
        {
            if (GetTracks()->IsNeutral(0) && (GetTracks()->GetClusterEnergy(0) > GetTracks()->GetClusterEnergy(1)))
            {
                compton = GetTracks()->GetVector(0);
                comp_time = GetTracks()->GetTime(0);
                energy = compton.E();
                theta = compton.Theta()*TMath::RadToDeg();
                phi = compton.Phi()*TMath::RadToDeg();
                clus_size = GetTracks()->GetClusterSize(0);
                if (lin_beam)
                {
                    if (phi < -90 || (phi > 0 && phi < 90)) fill_MM_0 = true;
                    else fill_MM_1 = true;
                }
                recoil = GetTracks()->GetUnitVector(1);
                reco_energy = GetTracks()->GetClusterEnergy(1);
                reco_theta = recoil.Theta()*TMath::RadToDeg();
                reco_size = GetTracks()->GetClusterSize(1);
                if (GetTracks()->IsNeutral(1)) comp_NN = true;
                else if (GetTracks()->GetClusterEnergy(1) > 0 && GetTracks()->HasCB(1)) comp_NC = true;
                else if (GetTracks()->GetClusterEnergy(1) > 0 && GetTracks()->HasTAPS(1)) comp_NT = true;
                else comp_NW = true;
            }
            else if (GetTracks()->IsNeutral(1) && (GetTracks()->GetClusterEnergy(1) > GetTracks()->GetClusterEnergy(0)))
            {
                compton = GetTracks()->GetVector(1);
                comp_time = GetTracks()->GetTime(1);
                energy = compton.E();
                theta = compton.Theta()*TMath::RadToDeg();
                phi = compton.Phi()*TMath::RadToDeg();
                clus_size = GetTracks()->GetClusterSize(1);
                if (lin_beam)
                {
                    if (phi < -90 || (phi > 0 && phi < 90)) fill_MM_0 = true;
                    else fill_MM_1 = true;
                }
                recoil = GetTracks()->GetUnitVector(0);
                reco_energy = GetTracks()->GetClusterEnergy(0);
                reco_theta = recoil.Theta()*TMath::RadToDeg();
                reco_size = GetTracks()->GetClusterSize(0);
                if (GetTracks()->IsNeutral(0)) comp_NN = true;
                else if (GetTracks()->GetClusterEnergy(0) > 0 && GetTracks()->HasCB(0)) comp_NC = true;
                else if (GetTracks()->GetClusterEnergy(0) > 0 && GetTracks()->HasTAPS(0)) comp_NT = true;
                else comp_NW = true;
            }
        }
    }
    else if (!pi0_X && GetTracks()->GetNTracks() == 3 && (GetTracks()->GetClusterEnergy(0) != 0 || GetTracks()->GetClusterEnergy(1) != 0 || GetTracks()->GetClusterEnergy(2) != 0))
    {
        if (GetTracks()->IsNeutral(0) && (GetTracks()->GetClusterEnergy(0) > GetTracks()->GetClusterEnergy(1)) && (GetTracks()->GetClusterEnergy(0) > GetTracks()->GetClusterEnergy(2)))
        {
            if (TMath::Abs(180.0-TMath::Abs(GetTracks()->GetPhi(0)-GetTracks()->GetPhi(1))) < TMath::Abs(180.0-TMath::Abs(GetTracks()->GetPhi(0)-GetTracks()->GetPhi(2))))
            {
                copl_ang = TMath::Abs(GetTracks()->GetPhi(0)-GetTracks()->GetPhi(1));
                Coplanarity->Fill(copl_ang,event_weight);

                if (TMath::Abs(180.0-copl_ang) < OACut)
                {
                    compton = GetTracks()->GetVector(0);
                    comp_time = GetTracks()->GetTime(0);
                    energy = compton.E();
                    theta = compton.Theta()*TMath::RadToDeg();
                    phi = compton.Phi()*TMath::RadToDeg();
                    clus_size = GetTracks()->GetClusterSize(0);
                    if (lin_beam)
                    {
                        if (phi < -90 || (phi > 0 && phi < 90)) fill_MM_0 = true;
                        else fill_MM_1 = true;
                    }
                    recoil = GetTracks()->GetUnitVector(1);
                    reco_energy = GetTracks()->GetClusterEnergy(1);
                    reco_theta = recoil.Theta()*TMath::RadToDeg();
                    reco_size = GetTracks()->GetClusterSize(1);
                    if (GetTracks()->IsNeutral(1)) comp_NNX = true;
                    else if (GetTracks()->GetClusterEnergy(1) > 0 && GetTracks()->HasCB(1)) comp_NCX = true;
                    else if (GetTracks()->GetClusterEnergy(1) > 0 && GetTracks()->HasTAPS(1)) comp_NTX = true;
                    else comp_NWX = true;
                    split = GetTracks()->GetUnitVector(2);
                    split_energy = GetTracks()->GetClusterEnergy(2);
                }
            }
            else if (TMath::Abs(180.0-TMath::Abs(GetTracks()->GetPhi(0)-GetTracks()->GetPhi(2))) < TMath::Abs(180.0-TMath::Abs(GetTracks()->GetPhi(0)-GetTracks()->GetPhi(1))))
            {
                copl_ang = TMath::Abs(GetTracks()->GetPhi(0)-GetTracks()->GetPhi(2));
                Coplanarity->Fill(copl_ang,event_weight);

                if (TMath::Abs(180.0-copl_ang) < OACut)
                {
                    compton = GetTracks()->GetVector(0);
                    comp_time = GetTracks()->GetTime(0);
                    energy = compton.E();
                    theta = compton.Theta()*TMath::RadToDeg();
                    phi = compton.Phi()*TMath::RadToDeg();
                    clus_size = GetTracks()->GetClusterSize(0);
                    if (lin_beam)
                    {
                        if (phi < -90 || (phi > 0 && phi < 90)) fill_MM_0 = true;
                        else fill_MM_1 = true;
                    }
                    recoil = GetTracks()->GetUnitVector(2);
                    reco_energy = GetTracks()->GetClusterEnergy(2);
                    reco_theta = recoil.Theta()*TMath::RadToDeg();
                    reco_size = GetTracks()->GetClusterSize(2);
                    if (GetTracks()->IsNeutral(2)) comp_NNX = true;
                    else if (GetTracks()->GetClusterEnergy(2) > 0 && GetTracks()->HasCB(2)) comp_NCX = true;
                    else if (GetTracks()->GetClusterEnergy(2) > 0 && GetTracks()->HasTAPS(2)) comp_NTX = true;
                    else comp_NWX = true;
                    split = GetTracks()->GetUnitVector(1);
                    split_energy = GetTracks()->GetClusterEnergy(1);
                }
            }
        }
        else if (GetTracks()->IsNeutral(1) && (GetTracks()->GetClusterEnergy(1) > GetTracks()->GetClusterEnergy(0)) && (GetTracks()->GetClusterEnergy(1) > GetTracks()->GetClusterEnergy(2)))
        {
            if (TMath::Abs(180.0-TMath::Abs(GetTracks()->GetPhi(1)-GetTracks()->GetPhi(0))) < TMath::Abs(180.0-TMath::Abs(GetTracks()->GetPhi(1)-GetTracks()->GetPhi(2))))
            {
                copl_ang = TMath::Abs(GetTracks()->GetPhi(1)-GetTracks()->GetPhi(0));
                Coplanarity->Fill(copl_ang,event_weight);

                if (TMath::Abs(180.0-copl_ang) < OACut)
                {
                    compton = GetTracks()->GetVector(1);
                    comp_time = GetTracks()->GetTime(1);
                    energy = compton.E();
                    theta = compton.Theta()*TMath::RadToDeg();
                    phi = compton.Phi()*TMath::RadToDeg();
                    clus_size = GetTracks()->GetClusterSize(1);
                    if (lin_beam)
                    {
                        if (phi < -90 || (phi > 0 && phi < 90)) fill_MM_0 = true;
                        else fill_MM_1 = true;
                    }
                    recoil = GetTracks()->GetUnitVector(0);
                    reco_energy = GetTracks()->GetClusterEnergy(0);
                    reco_theta = recoil.Theta()*TMath::RadToDeg();
                    reco_size = GetTracks()->GetClusterSize(0);
                    if (GetTracks()->IsNeutral(0)) comp_NNX = true;
                    else if (GetTracks()->GetClusterEnergy(0) > 0 && GetTracks()->HasCB(0)) comp_NCX = true;
                    else if (GetTracks()->GetClusterEnergy(0) > 0 && GetTracks()->HasTAPS(0)) comp_NTX = true;
                    else comp_NWX = true;
                    split = GetTracks()->GetUnitVector(2);
                    split_energy = GetTracks()->GetClusterEnergy(2);
                }
            }
            else if (TMath::Abs(180.0-TMath::Abs(GetTracks()->GetPhi(1)-GetTracks()->GetPhi(2))) < TMath::Abs(180.0-TMath::Abs(GetTracks()->GetPhi(1)-GetTracks()->GetPhi(0))))
            {
                copl_ang = TMath::Abs(GetTracks()->GetPhi(1)-GetTracks()->GetPhi(2));
                Coplanarity->Fill(copl_ang,event_weight);

                if (TMath::Abs(180.0-copl_ang) < OACut)
                {
                    compton = GetTracks()->GetVector(1);
                    comp_time = GetTracks()->GetTime(1);
                    energy = compton.E();
                    theta = compton.Theta()*TMath::RadToDeg();
                    phi = compton.Phi()*TMath::RadToDeg();
                    clus_size = GetTracks()->GetClusterSize(1);
                    if (lin_beam)
                    {
                        if (phi < -90 || (phi > 0 && phi < 90)) fill_MM_0 = true;
                        else fill_MM_1 = true;
                    }
                    recoil = GetTracks()->GetUnitVector(2);
                    reco_energy = GetTracks()->GetClusterEnergy(2);
                    reco_theta = recoil.Theta()*TMath::RadToDeg();
                    reco_size = GetTracks()->GetClusterSize(2);
                    if (GetTracks()->IsNeutral(2)) comp_NNX = true;
                    else if (GetTracks()->GetClusterEnergy(2) > 0 && GetTracks()->HasCB(2)) comp_NCX = true;
                    else if (GetTracks()->GetClusterEnergy(2) > 0 && GetTracks()->HasTAPS(2)) comp_NTX = true;
                    else comp_NWX = true;
                    split = GetTracks()->GetUnitVector(0);
                    split_energy = GetTracks()->GetClusterEnergy(0);
                }
            }
        }
        else if (GetTracks()->IsNeutral(2) && (GetTracks()->GetClusterEnergy(2) > GetTracks()->GetClusterEnergy(0)) && (GetTracks()->GetClusterEnergy(2) > GetTracks()->GetClusterEnergy(1)))
        {
            if (TMath::Abs(180.0-TMath::Abs(GetTracks()->GetPhi(2)-GetTracks()->GetPhi(0))) < TMath::Abs(180.0-TMath::Abs(GetTracks()->GetPhi(2)-GetTracks()->GetPhi(1))))
            {
                copl_ang = TMath::Abs(GetTracks()->GetPhi(2)-GetTracks()->GetPhi(0));
                Coplanarity->Fill(copl_ang,event_weight);

                if (TMath::Abs(180.0-copl_ang) < OACut)
                {
                    compton = GetTracks()->GetVector(2);
                    comp_time = GetTracks()->GetTime(2);
                    energy = compton.E();
                    theta = compton.Theta()*TMath::RadToDeg();
                    phi = compton.Phi()*TMath::RadToDeg();
                    clus_size = GetTracks()->GetClusterSize(2);
                    if (lin_beam)
                    {
                        if (phi < -90 || (phi > 0 && phi < 90)) fill_MM_0 = true;
                        else fill_MM_1 = true;
                    }
                    recoil = GetTracks()->GetUnitVector(0);
                    reco_energy = GetTracks()->GetClusterEnergy(0);
                    reco_theta = recoil.Theta()*TMath::RadToDeg();
                    reco_size = GetTracks()->GetClusterSize(0);
                    if (GetTracks()->IsNeutral(0)) comp_NNX = true;
                    else if (GetTracks()->GetClusterEnergy(0) > 0 && GetTracks()->HasCB(0)) comp_NCX = true;
                    else if (GetTracks()->GetClusterEnergy(0) > 0 && GetTracks()->HasTAPS(0)) comp_NTX = true;
                    else comp_NWX = true;
                    split = GetTracks()->GetUnitVector(1);
                    split_energy = GetTracks()->GetClusterEnergy(1);
                }
            }
            else if (TMath::Abs(180.0-TMath::Abs(GetTracks()->GetPhi(2)-GetTracks()->GetPhi(1))) < TMath::Abs(180.0-TMath::Abs(GetTracks()->GetPhi(2)-GetTracks()->GetPhi(0))))
            {
                copl_ang = TMath::Abs(GetTracks()->GetPhi(2)-GetTracks()->GetPhi(1));
                Coplanarity->Fill(copl_ang,event_weight);

                if (TMath::Abs(180.0-copl_ang) < OACut)
                {
                    compton = GetTracks()->GetVector(2);
                    comp_time = GetTracks()->GetTime(2);
                    energy = compton.E();
                    theta = compton.Theta()*TMath::RadToDeg();
                    phi = compton.Phi()*TMath::RadToDeg();
                    clus_size = GetTracks()->GetClusterSize(2);
                    if (lin_beam)
                    {
                        if (phi < -90 || (phi > 0 && phi < 90)) fill_MM_0 = true;
                        else fill_MM_1 = true;
                    }
                    recoil = GetTracks()->GetUnitVector(1);
                    reco_energy = GetTracks()->GetClusterEnergy(1);
                    reco_theta = recoil.Theta()*TMath::RadToDeg();
                    reco_size = GetTracks()->GetClusterSize(1);
                    if (GetTracks()->IsNeutral(1)) comp_NNX = true;
                    else if (GetTracks()->GetClusterEnergy(1) > 0 && GetTracks()->HasCB(1)) comp_NCX = true;
                    else if (GetTracks()->GetClusterEnergy(1) > 0 && GetTracks()->HasTAPS(1)) comp_NTX = true;
                    else comp_NWX = true;
                    split = GetTracks()->GetUnitVector(0);
                    split_energy = GetTracks()->GetClusterEnergy(0);
                }
            }
        }
    }
    comp_X = comp_N || comp_C || comp_NN || comp_NC || comp_NT || comp_NW || comp_NNX || comp_NCX || comp_NTX || comp_NWX;

    if (split_energy > 0) Split_OA_E->Fill(split_energy/compton.E(),(TMath::RadToDeg()*compton.Angle(split)));
    if (comp_X && !comp_C) Comp_CS->Fill(energy, theta, clus_size, event_weight);
    if (comp_NN || comp_NC || comp_NT || comp_NNX || comp_NCX || comp_NTX) Reco_CS->Fill(reco_energy, reco_theta, reco_size, event_weight);

    Int_t num_clust = 0;
    Double_t avg_time = 0;
    for (Int_t i = 0; i < GetTracks()->GetNTracks(); i++)
    {
        if (GetTracks()->GetClusterEnergy(i) != 0)
        {
            avg_time += GetTracks()->GetTime(i);
            num_clust++;
        }
    }
    avg_time = avg_time/num_clust;

    //////////////////////////////////////////////////
    // Tagger loop
    //////////////////////////////////////////////////
    for (Int_t i = 0; i < nTagg; i++)
    {
        if (RejectTagged(i)) continue;

        tagg_chan = GetTagger()->GetTaggedChannel(i);
        tagg_time = GetTagger()->GetTaggedTime(i);
        tagg_energy = GetTagger()->GetTaggedEnergy(i);
        beam = TLorentzVector(0., 0., tagg_energy, tagg_energy);

        TaggTime->Fill(tagg_time);
        if (hel)
        {
            if (GHistBGSub::IsPrompt(tagg_time)) TaggHel1->Fill(tagg_chan);
            else if (GHistBGSub::IsRandom(tagg_time)) TaggHel1_R->Fill(tagg_chan);
        }
        else
        {
            if (GHistBGSub::IsPrompt(tagg_time)) TaggHel0->Fill(tagg_chan);
            else if (GHistBGSub::IsRandom(tagg_time)) TaggHel0_R->Fill(tagg_chan);
        }

        time = tagg_time - avg_time;
        IncTime->Fill(time);
        if (hel)
        {
            if (GHistBGSub::IsPrompt(time)) IncHel1->Fill(tagg_chan);
            else if (GHistBGSub::IsRandom(time)) IncHel1_R->Fill(tagg_chan);
        }
        else
        {
            if (GHistBGSub::IsPrompt(time)) IncHel0->Fill(tagg_chan);
            else if (GHistBGSub::IsRandom(time)) IncHel0_R->Fill(tagg_chan);
        }
        if (GHistBGSub::IsPrompt(time)) IncHits->Fill(tagg_chan,GetTracks()->GetNTracks());
        else if (GHistBGSub::IsRandom(time)) IncHits_R->Fill(tagg_chan,GetTracks()->GetNTracks());

        if (minIM > 0)
        {
            time = tagg_time - pi0_time;
            missing = beam + target - pi0;
            m_miss = missing.M()-target.M();
            Pi0_Dt->Fill(tagg_chan,time,event_weight);

            if (GHistBGSub::IsPrompt(time))
            {
                Pi0_MM->Fill(tagg_chan,pi0.Theta()*TMath::RadToDeg(),m_miss,event_weight);
                Pi0_IM_MM->Fill(tagg_chan,minIM,m_miss,event_weight);
            }
            else if (GHistBGSub::IsRandom(time))
            {
                Pi0_MM_R->Fill(tagg_chan,pi0.Theta()*TMath::RadToDeg(),m_miss,event_weight);
                Pi0_IM_MM_R->Fill(tagg_chan,minIM,m_miss,event_weight);
            }
        }

        if (pi0_X)
        {
            time = tagg_time - pi0_time;
            missing = beam + target - pi0;
            m_miss = missing.M()-target.M();

            if (GHistBGSub::IsPrompt(time))
            {
                Pi0_MM_IM->Fill(tagg_chan,pi0.Theta()*TMath::RadToDeg(),m_miss,event_weight);
                Pi0_Re_All->Fill(missing.E()-missing.M(),missing.Theta()*TMath::RadToDeg(),m_miss,event_weight);
            }
            else if (GHistBGSub::IsRandom(time))
            {
                Pi0_MM_IM_R->Fill(tagg_chan,pi0.Theta()*TMath::RadToDeg(),m_miss,event_weight);
                Pi0_Re_All_R->Fill(missing.E()-missing.M(),missing.Theta()*TMath::RadToDeg(),m_miss,event_weight);
            }

            if (pi0_NN)
            {
                Pi0_Dt_NN->Fill(tagg_chan,time,event_weight);
                if (GHistBGSub::IsPrompt(time))
                {
                    Pi0_MM_NN->Fill(tagg_chan,pi0.Theta()*TMath::RadToDeg(),m_miss,event_weight);
                    Pi0_IM_MM_NN->Fill(tagg_chan,minIM,m_miss,event_weight);
                }
                else if (GHistBGSub::IsRandom(time))
                {
                    Pi0_MM_NN_R->Fill(tagg_chan,pi0.Theta()*TMath::RadToDeg(),m_miss,event_weight);
                    Pi0_IM_MM_NN_R->Fill(tagg_chan,minIM,m_miss,event_weight);
                }
            }
            else if (pi0_NC)
            {
                if (GHistBGSub::IsPrompt(time)) Pi0_MM_NC->Fill(tagg_chan,pi0.Theta()*TMath::RadToDeg(),m_miss,event_weight);
                else if (GHistBGSub::IsRandom(time)) Pi0_MM_NC_R->Fill(tagg_chan,pi0.Theta()*TMath::RadToDeg(),m_miss,event_weight);
            }
            else if (pi0_CC)
            {
                if (GHistBGSub::IsPrompt(time)) Pi0_MM_CC->Fill(tagg_chan,pi0.Theta()*TMath::RadToDeg(),m_miss,event_weight);
                else if (GHistBGSub::IsRandom(time)) Pi0_MM_CC_R->Fill(tagg_chan,pi0.Theta()*TMath::RadToDeg(),m_miss,event_weight);
            }
            else if (pi0_NNX)
            {
                if (GHistBGSub::IsPrompt(time)) Pi0_MM_NNX->Fill(tagg_chan,pi0.Theta()*TMath::RadToDeg(),m_miss,event_weight);
                else if (GHistBGSub::IsRandom(time)) Pi0_MM_NNX_R->Fill(tagg_chan,pi0.Theta()*TMath::RadToDeg(),m_miss,event_weight);

            }
            else if (pi0_NCX)
            {
                if (GHistBGSub::IsPrompt(time)) Pi0_MM_NCX->Fill(tagg_chan,pi0.Theta()*TMath::RadToDeg(),m_miss,event_weight);
                else if (GHistBGSub::IsRandom(time)) Pi0_MM_NCX_R->Fill(tagg_chan,pi0.Theta()*TMath::RadToDeg(),m_miss,event_weight);
            }
            else if (pi0_CCX)
            {
                if (GHistBGSub::IsPrompt(time)) Pi0_MM_CCX->Fill(tagg_chan,pi0.Theta()*TMath::RadToDeg(),m_miss,event_weight);
                else if (GHistBGSub::IsRandom(time)) Pi0_MM_CCX_R->Fill(tagg_chan,pi0.Theta()*TMath::RadToDeg(),m_miss,event_weight);
            }

            if (cut_CA)
            {
                open_ang = (TMath::RadToDeg()*missing.Angle(recoil));
                if (open_ang < OACut)
                {
                    if (reco_energy > 0)
                    {
                        if (GHistBGSub::IsPrompt(time)) Pi0_Re_Det->Fill(missing.E()-missing.M(),missing.Theta()*TMath::RadToDeg(),m_miss,event_weight);
                        else if (GHistBGSub::IsRandom(time)) Pi0_Re_Det_R->Fill(missing.E()-missing.M(),missing.Theta()*TMath::RadToDeg(),m_miss,event_weight);
                    }
                    else
                    {
                        if (GHistBGSub::IsPrompt(time)) Pi0_Re_NoE->Fill(missing.E()-missing.M(),missing.Theta()*TMath::RadToDeg(),m_miss,event_weight);
                        else if (GHistBGSub::IsRandom(time)) Pi0_Re_NoE_R->Fill(missing.E()-missing.M(),missing.Theta()*TMath::RadToDeg(),m_miss,event_weight);
                    }
                }
            }
        }

        else if(comp_X)
        {
            time = tagg_time - comp_time;
            missing = beam + target - compton;
            m_miss = missing.M()-target.M();
            if (!hel) fill_Ph_0 = (missing.M() >= MMLoC && missing.M() < MMHiC);
            else fill_Ph_1 = (missing.M() >= MMLoC && missing.M() < MMHiC);

            ptot = beam + target;
            lab_to_cm = -ptot.BoostVector();
            compton_cm = compton;
            compton_cm.Boost(lab_to_cm);
            theta_cm = compton_cm.Theta()*TMath::RadToDeg();
            phi_cm = compton_cm.Phi()*TMath::RadToDeg();

            if (!comp_N && !comp_C)
            {
                open_ang = (TMath::RadToDeg()*missing.Angle(recoil));
                cut_OA = (open_ang < OACut);
                if (GHistBGSub::IsPrompt(time)) MM_CA_OA->Fill(m_miss,copl_ang,open_ang,event_weight);
                else if (GHistBGSub::IsRandom(time)) MM_CA_OA_R->Fill(m_miss,copl_ang,open_ang,event_weight);
            }

            if (comp_N)
            {
                if(TMath::Abs(m_miss) < 15) Comp_Dt_N_MM->Fill(tagg_chan,time,event_weight);
                if(GetTracks()->GetClusterSize(0) > 1) Comp_Dt_N_CS->Fill(tagg_chan,time,event_weight);
                if((TMath::Abs(m_miss) < 15) && (GetTracks()->GetClusterSize(0) > 1)) Comp_Dt_N_MMCS->Fill(tagg_chan,time,event_weight);
                Comp_Dt_N->Fill(tagg_chan,time,event_weight);
                if (GHistBGSub::IsPrompt(time))
                {
                    if (GetTracks()->GetClusterSize(0) > 1) Comp_MM_N_C->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (TMath::Abs(m_miss) < 15) Comp_CS_MM->Fill(energy, theta, clus_size, event_weight);
                    if (fill_MM_0) Comp_MM_N_0->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_N_1->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_N_0->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_N_1->Fill(tagg_chan,theta,phi,event_weight);
                }
                else if (GHistBGSub::IsRandom(time))
                {
                    if(GetTracks()->GetClusterSize(0) > 1) Comp_MM_N_C_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (TMath::Abs(m_miss) < 15) Comp_CS_MM_R->Fill(energy, theta, clus_size, event_weight);
                    if (fill_MM_0) Comp_MM_N_0_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_N_1_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_N_0_R->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_N_1_R->Fill(tagg_chan,theta,phi,event_weight);
                }
            }
            else if (comp_C)
            {
                Comp_Dt_C->Fill(tagg_chan,time,event_weight);
                if (GHistBGSub::IsPrompt(time))
                {
                    if (fill_MM_0) Comp_MM_C_0->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_C_1->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_C_0->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_C_1->Fill(tagg_chan,theta,phi,event_weight);
                }
                else if (GHistBGSub::IsRandom(time))
                {
                    if (fill_MM_0) Comp_MM_C_0_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_C_1_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_C_0_R->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_C_1_R->Fill(tagg_chan,theta,phi,event_weight);
                }
            }
            else if (comp_NN && cut_OA)
            {
                Comp_Dt_NN->Fill(tagg_chan,time,event_weight);
                if (GHistBGSub::IsPrompt(time))
                {
                    if (TMath::Abs(m_miss) < 15)
                    {
                        Comp_CS_MM->Fill(energy, theta, clus_size, event_weight);
                        Reco_CS_MM->Fill(reco_energy, reco_theta, reco_size, event_weight);
                    }
                    if (fill_MM_0) Comp_MM_NN_0->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_NN_1->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_NN_0->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_NN_1->Fill(tagg_chan,theta,phi,event_weight);
                }
                else if (GHistBGSub::IsRandom(time))
                {
                    if (TMath::Abs(m_miss) < 15)
                    {
                        Comp_CS_MM_R->Fill(energy, theta, clus_size, event_weight);
                        Reco_CS_MM_R->Fill(reco_energy, reco_theta, reco_size, event_weight);
                    }
                    if (fill_MM_0) Comp_MM_NN_0_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_NN_1_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_NN_0_R->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_NN_1_R->Fill(tagg_chan,theta,phi,event_weight);
                }
            }
            else if (comp_NC && cut_OA)
            {
                Comp_Dt_NC->Fill(tagg_chan,time,event_weight);
                if (GHistBGSub::IsPrompt(time))
                {
                    if (TMath::Abs(m_miss) < 15)
                    {
                        Comp_CS_MM->Fill(energy, theta, clus_size, event_weight);
                        Reco_CS_MM->Fill(reco_energy, reco_theta, reco_size, event_weight);
                    }
                    if (fill_MM_0) Comp_MM_NC_0->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_NC_1->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_NC_0->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_NC_1->Fill(tagg_chan,theta,phi,event_weight);
                }
                else if (GHistBGSub::IsRandom(time))
                {
                    if (TMath::Abs(m_miss) < 15)
                    {
                        Comp_CS_MM_R->Fill(energy, theta, clus_size, event_weight);
                        Reco_CS_MM_R->Fill(reco_energy, reco_theta, reco_size, event_weight);
                    }
                    if (fill_MM_0) Comp_MM_NC_0_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_NC_1_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_NC_0_R->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_NC_1_R->Fill(tagg_chan,theta,phi,event_weight);
                }
            }
            else if (comp_NT && cut_OA)
            {
                Comp_Dt_NT->Fill(tagg_chan,time,event_weight);
                if (GHistBGSub::IsPrompt(time))
                {
                    if (TMath::Abs(m_miss) < 15)
                    {
                        Comp_CS_MM->Fill(energy, theta, clus_size, event_weight);
                        Reco_CS_MM->Fill(reco_energy, reco_theta, reco_size, event_weight);
                    }
                    if (fill_MM_0) Comp_MM_NT_0->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_NT_1->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_NT_0->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_NT_1->Fill(tagg_chan,theta,phi,event_weight);
                }
                else if (GHistBGSub::IsRandom(time))
                {
                    if (TMath::Abs(m_miss) < 15)
                    {
                        Comp_CS_MM_R->Fill(energy, theta, clus_size, event_weight);
                        Reco_CS_MM_R->Fill(reco_energy, reco_theta, reco_size, event_weight);
                    }
                    if (fill_MM_0) Comp_MM_NT_0_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_NT_1_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_NT_0_R->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_NT_1_R->Fill(tagg_chan,theta,phi,event_weight);
                }
            }
            else if (comp_NW && cut_OA)
            {
                Comp_Dt_NW->Fill(tagg_chan,time,event_weight);
                if (GHistBGSub::IsPrompt(time))
                {
                    if (TMath::Abs(m_miss) < 15) Comp_CS_MM->Fill(energy, theta, clus_size, event_weight);
                    if (fill_MM_0) Comp_MM_NW_0->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_NW_1->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_NW_0->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_NW_1->Fill(tagg_chan,theta,phi,event_weight);
                }
                else if (GHistBGSub::IsRandom(time))
                {
                    if (TMath::Abs(m_miss) < 15) Comp_CS_MM_R->Fill(energy, theta, clus_size, event_weight);
                    if (fill_MM_0) Comp_MM_NW_0_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_NW_1_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_NW_0_R->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_NW_1_R->Fill(tagg_chan,theta,phi,event_weight);
                }
            }
            else if (comp_NNX && cut_OA)
            {
                Comp_Dt_NNX->Fill(tagg_chan,time,event_weight);
                if (GHistBGSub::IsPrompt(time))
                {
                    if (TMath::Abs(m_miss) < 15)
                    {
                        Comp_CS_MM->Fill(energy, theta, clus_size, event_weight);
                        Reco_CS_MM->Fill(reco_energy, reco_theta, reco_size, event_weight);
                    }
                    if (fill_MM_0) Comp_MM_NNX_0->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_NNX_1->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_NNX_0->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_NNX_1->Fill(tagg_chan,theta,phi,event_weight);
                }
                else if (GHistBGSub::IsRandom(time))
                {
                    if (TMath::Abs(m_miss) < 15)
                    {
                        Comp_CS_MM_R->Fill(energy, theta, clus_size, event_weight);
                        Reco_CS_MM_R->Fill(reco_energy, reco_theta, reco_size, event_weight);
                    }
                    if (fill_MM_0) Comp_MM_NNX_0_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_NNX_1_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_NNX_0_R->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_NNX_1_R->Fill(tagg_chan,theta,phi,event_weight);
                }
            }
            else if (comp_NCX && cut_OA)
            {
                Comp_Dt_NCX->Fill(tagg_chan,time,event_weight);
                if (GHistBGSub::IsPrompt(time))
                {
                    if (TMath::Abs(m_miss) < 15)
                    {
                        Comp_CS_MM->Fill(energy, theta, clus_size, event_weight);
                        Reco_CS_MM->Fill(reco_energy, reco_theta, reco_size, event_weight);
                    }
                    if (fill_MM_0) Comp_MM_NCX_0->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_NCX_1->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_NCX_0->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_NCX_1->Fill(tagg_chan,theta,phi,event_weight);
                }
                else if (GHistBGSub::IsRandom(time))
                {
                    if (TMath::Abs(m_miss) < 15)
                    {
                        Comp_CS_MM_R->Fill(energy, theta, clus_size, event_weight);
                        Reco_CS_MM_R->Fill(reco_energy, reco_theta, reco_size, event_weight);
                    }
                    if (fill_MM_0) Comp_MM_NCX_0_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_NCX_1_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_NCX_0_R->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_NCX_1_R->Fill(tagg_chan,theta,phi,event_weight);
                }
            }
            else if (comp_NTX && cut_OA)
            {
                Comp_Dt_NTX->Fill(tagg_chan,time,event_weight);
                if (GHistBGSub::IsPrompt(time))
                {
                    if (TMath::Abs(m_miss) < 15)
                    {
                        Comp_CS_MM->Fill(energy, theta, clus_size, event_weight);
                        Reco_CS_MM->Fill(reco_energy, reco_theta, reco_size, event_weight);
                    }
                    if (fill_MM_0) Comp_MM_NTX_0->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_NTX_1->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_NTX_0->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_NTX_1->Fill(tagg_chan,theta,phi,event_weight);
                }
                else if (GHistBGSub::IsRandom(time))
                {
                    if (TMath::Abs(m_miss) < 15)
                    {
                        Comp_CS_MM_R->Fill(energy, theta, clus_size, event_weight);
                        Reco_CS_MM_R->Fill(reco_energy, reco_theta, reco_size, event_weight);
                    }
                    if (fill_MM_0) Comp_MM_NTX_0_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_NTX_1_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_NTX_0_R->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_NTX_1_R->Fill(tagg_chan,theta,phi,event_weight);
                }
            }
            else if (comp_NWX && cut_OA)
            {
                Comp_Dt_NWX->Fill(tagg_chan,time,event_weight);
                if (GHistBGSub::IsPrompt(time))
                {
                    if (TMath::Abs(m_miss) < 15) Comp_CS_MM->Fill(energy, theta, clus_size, event_weight);
                    if (fill_MM_0) Comp_MM_NWX_0->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_NWX_1->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_NWX_0->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_NWX_1->Fill(tagg_chan,theta,phi,event_weight);
                }
                else if (GHistBGSub::IsRandom(time))
                {
                    if (TMath::Abs(m_miss) < 15) Comp_CS_MM_R->Fill(energy, theta, clus_size, event_weight);
                    if (fill_MM_0) Comp_MM_NWX_0_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    else if (fill_MM_1) Comp_MM_NWX_1_R->Fill(tagg_chan,theta,m_miss,event_weight);
                    if (fill_Ph_0) Comp_Ph_NWX_0_R->Fill(tagg_chan,theta,phi,event_weight);
                    else if (fill_Ph_1) Comp_Ph_NWX_1_R->Fill(tagg_chan,theta,phi,event_weight);
                }
            }
        }
    }
}

Bool_t 	PAnalyze::InitVerbosity()
{
    Int_t sc1;
    string config = ReadConfig("Verbosity");
    if(sscanf( config.c_str(), "%d\n", &sc1) == 1)
    {
        cout << "Setting verbosity: " << sc1 << endl << endl;
        verbosity = sc1;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Verbosity not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

Bool_t 	PAnalyze::InitExclusivity()
{
    Int_t sc1, sc2;
    string config = ReadConfig("Exclusivity");
    if(sscanf( config.c_str(), "%d%d\n", &sc1, &sc2) == 2)
    {
        cout << "Setting exclusivity: pi0 = " << sc1 << ", p = " << sc2 << endl << endl;
        excl_pi0 = sc1;
        excl_pro = sc2;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Exclusivity not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

Bool_t 	PAnalyze::InitInvariantMass()
{
    Double_t sc1;
    string config = ReadConfig("Invariant-Mass-Cut");
    if(sscanf( config.c_str(), "%lf\n", &sc1) == 1)
    {
        cout << "Setting invariant mass width: +/-" << sc1 << " MeV " << endl << endl;
        IMCut = sc1;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Invariant mass cut not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

Bool_t 	PAnalyze::InitMissingMass()
{
    Double_t sc1, sc2;
    string config = ReadConfig("Missing-Mass-Cut");
    if(sscanf( config.c_str(), "%lf%lf\n", &sc1, &sc2) == 2)
    {
        cout << "Setting missing mass cut: " << sc1 << " to " << sc2 << " MeV " << endl << endl;
        MMLoC = sc1;
        MMHiC = sc2;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Missing mass cut not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

Bool_t 	PAnalyze::InitOpeningAngle()
{
    Double_t sc1;
    string config = ReadConfig("Opening-Angle-Cut");
    if(sscanf( config.c_str(), "%lf\n", &sc1) == 1)
    {
        cout << "Setting opening angle: " << sc1 << " deg " << endl << endl;
        OACut = sc1;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Opening angle cut not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

Bool_t 	PAnalyze::InitEnergySum()
{
    Double_t sc1;
    string config = ReadConfig("Energy-Sum-Cut");
    if(sscanf( config.c_str(), "%lf\n", &sc1) == 1)
    {
        cout << "Setting energy sum: " << sc1 << " MeV " << endl << endl;
        ESCut = sc1;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Energy sum cut not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

Bool_t 	PAnalyze::InitSaveRandoms()
{
    Int_t sc1;
    string config = ReadConfig("Save-Randoms");
    if(sscanf( config.c_str(), "%d\n", &sc1) == 1)
    {
        cout << "Save random hists: " << sc1 << endl << endl;
        save_randoms = sc1;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Save randoms not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

Bool_t  PAnalyze::InitBeamPol()
{
    string config;
    Int_t instance = 0;
    Int_t time;
    Double_t meas;

    do
    {
        config = ReadConfig("Beam-Polarization",instance);

        if(sscanf( config.c_str(), "%d %lf\n", &time, &meas) == 2)
        {
            beamPolTime.push_back(time);
            beamPolMeas.push_back(meas);
            cout << "Beam polarization measurement: " << 100*meas << "% at " << time << endl;
            instance++;
        }
        else if(strcmp(config.c_str(), "nokey") != 0)
        {
            cout << "Beam polarization not set correctly" << endl;
            return kFALSE;
        }
    } while (strcmp(config.c_str(), "nokey") != 0);

    if(instance)
    {
        cir_beam = true;
        cout << endl;
    }

    Int_t sc1;
    config = ReadConfig("Linear-Pol");
    if(sscanf( config.c_str(), "%d\n", &sc1) == 1)
    {
        cout << "Use linear polarization: " << sc1 << endl << endl;
        lin_beam = sc1;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Linear polarization not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;
}

Bool_t  PAnalyze::InitTargPol()
{
    string config;
    Int_t instance = 0;
    Int_t time1, time2;
    Double_t meas1, meas2;

    do
    {
        config = ReadConfig("Target-Polarization",instance);

        if(sscanf( config.c_str(), "%d %lf %d %lf\n", &time1, &meas1, &time2, &meas2) == 4)
        {
            targPolTime.push_back(time1);
            targPolMeas.push_back(meas1);
            targPolTime.push_back(time2);
            targPolMeas.push_back(meas2);
            cout << "Target polarization measurement: " << endl;
            cout << "\tInitial: " << 100*meas1 << "% at " << time1 << endl;
            cout << "\tFinal:   " << 100*meas2 << "% at " << time2 << endl;
            instance++;
        }
        else if(strcmp(config.c_str(), "nokey") != 0)
        {
            cout << "Target polarization not set correctly" << endl;
            return kFALSE;
        }
    } while (strcmp(config.c_str(), "nokey") != 0);

    if(instance) cout << endl;

    return kTRUE;
}

Double_t  PAnalyze::TwoBodyAngleToEnergyMin(Double_t eBeam, Double_t mTarg, Double_t mPar1, Double_t mPar2, Double_t tPar1)
{
    Double_t dX = ((eBeam*(mPar1-mTarg))-(0.5*(TMath::Power((mTarg-mPar1),2.0)))+(0.5*mPar2*mPar2));

    Double_t dA = ((TMath::Power(eBeam*TMath::Cos(tPar1*(TMath::DegToRad())),2.0))-(TMath::Power((eBeam+mTarg),2.0)));

    Double_t dB = ((2.0*mPar1*(TMath::Power(eBeam*TMath::Cos(tPar1*(TMath::DegToRad())),2.0)))-((2.0*dX*(eBeam+mTarg))));

    Double_t dC = (-dX*dX);

    Double_t dY = ((dB*dB)-(4.0*dA*dC));

    if (dA==0 || dY<0) return 0;

    Double_t dP = ((-dB+(TMath::Sqrt(dY)))/(2.0*dA));
    Double_t dM = ((-dB-(TMath::Sqrt(dY)))/(2.0*dA));
    return TMath::Min(dP,dM);
}

Double_t  PAnalyze::TwoBodyAngleToEnergyMax(Double_t eBeam, Double_t mTarg, Double_t mPar1, Double_t mPar2, Double_t tPar1)
{
    Double_t dX = ((eBeam*(mPar1-mTarg))-(0.5*(TMath::Power((mTarg-mPar1),2.0)))+(0.5*mPar2*mPar2));

    Double_t dA = ((TMath::Power(eBeam*TMath::Cos(tPar1*(TMath::DegToRad())),2.0))-(TMath::Power((eBeam+mTarg),2.0)));

    Double_t dB = ((2.0*mPar1*(TMath::Power(eBeam*TMath::Cos(tPar1*(TMath::DegToRad())),2.0)))-((2.0*dX*(eBeam+mTarg))));

    Double_t dC = (-dX*dX);

    Double_t dY = ((dB*dB)-(4.0*dA*dC));

    if (dA==0 || dY<0) return 0;

    Double_t dP = ((-dB+(TMath::Sqrt(dY)))/(2.0*dA));
    Double_t dM = ((-dB-(TMath::Sqrt(dY)))/(2.0*dA));
    return TMath::Max(dP,dM);
}

Double_t  PAnalyze::TwoBodyEnergyToAngle(Double_t eBeam, Double_t mTarg, Double_t mPar1, Double_t mPar2, Double_t ePar1)
{
    return ((TMath::RadToDeg())*(TMath::ACos(((ePar1*(eBeam+mTarg))+(eBeam*(mPar1-mTarg))-(0.5*(TMath::Power((mTarg-mPar1),2.0)))+(0.5*mPar2*mPar2))/(eBeam*TMath::Sqrt((ePar1*ePar1)+(2.0*ePar1*mPar1))))));
}

Double_t PAnalyze::CalcCircBeamPol(Double_t E_e, Double_t P_e, Double_t E_g)
{
    Double_t P_g = P_e*(((4*E_g*E_e)-(E_g*E_g))/((4*E_e*E_e)-(4*E_g*E_e)+(3*E_g*E_g)));

    return P_g;
}

void	PAnalyze::ProcessScalerRead()
{
    PPhysics::ProcessScalerRead();
}

Bool_t	PAnalyze::Write()
{
    Double_t livetime = ((LiveTimeScal->GetBinContent(2))/(LiveTimeScal->GetBinContent(1)));
    Double_t counts;
    Double_t photonPol = 1;
    for (Int_t i=0; i<352; i++)
    {
        counts = TaggerAccScal->GetBinContent(i+1);
        if(beamPol < 1) photonPol = CalcCircBeamPol(450,beamPol,GetSetupParameters()->GetTaggerPhotonEnergy(i));
        CorrTaggScal->SetBinContent(i+1,(livetime*counts));
        PolarizeScal->SetBinContent(i+1,(livetime*targPol*photonPol*counts));
    }

    Double_t ratio = GHistBGSub::GetBackgroundSumbtractionFactor();

    if(save_randoms && !(IsMCFile()))
    {
        TH1D *TaggHel0_P = (TH1D*)TaggHel0->Clone("TaggHel0_P");
        TH1D *TaggHel1_P = (TH1D*)TaggHel1->Clone("TaggHel1_P");

        TH1D *IncHel0_P = (TH1D*)IncHel0->Clone("IncHel0_P");
        TH1D *IncHel1_P = (TH1D*)IncHel1->Clone("IncHel1_P");
        TH2D *IncHits_P = (TH2D*)IncHits->Clone("IncHits_P");

        TH3D *MM_CA_OA_P = (TH3D*)MM_CA_OA->Clone("MM_CA_OA_P");

        TH3D *Pi0_IM_MM_P = (TH3D*)Pi0_IM_MM->Clone("Pi0_IM_MM_P");
        TH3D *Pi0_IM_MM_NN_P = (TH3D*)Pi0_IM_MM_NN->Clone("Pi0_IM_MM_NN_P");

        TH3D *Pi0_MM_P = (TH3D*)Pi0_MM->Clone("Pi0_MM_P");
        TH3D *Pi0_MM_IM_P = (TH3D*)Pi0_MM_IM->Clone("Pi0_MM_IM_P");

        TH3D *Pi0_MM_NN_P = (TH3D*)Pi0_MM_NN->Clone("Pi0_MM_NN_P");
        TH3D *Pi0_MM_NC_P = (TH3D*)Pi0_MM_NC->Clone("Pi0_MM_NC_P");
        TH3D *Pi0_MM_CC_P = (TH3D*)Pi0_MM_CC->Clone("Pi0_MM_CC_P");

        TH3D *Pi0_MM_NNX_P = (TH3D*)Pi0_MM_NNX->Clone("Pi0_MM_NNX_P");
        TH3D *Pi0_MM_NCX_P = (TH3D*)Pi0_MM_NCX->Clone("Pi0_MM_NCX_P");
        TH3D *Pi0_MM_CCX_P = (TH3D*)Pi0_MM_CCX->Clone("Pi0_MM_CCX_P");

        TH3D *Pi0_Re_All_P = (TH3D*)Pi0_Re_All->Clone("Pi0_Re_All_P");
        TH3D *Pi0_Re_Det_P = (TH3D*)Pi0_Re_Det->Clone("Pi0_Re_Det_P");
        TH3D *Pi0_Re_NoE_P = (TH3D*)Pi0_Re_NoE->Clone("Pi0_Re_NoE_P");

        TH3D *Comp_MM_N_C_P = (TH3D*)Comp_MM_N_C->Clone("Comp_MM_N_C_P");

        TH3D *Comp_CS_MM_P = (TH3D*)Comp_CS_MM->Clone("Comp_CS_MM_P");
        TH3D *Reco_CS_MM_P = (TH3D*)Reco_CS_MM->Clone("Reco_CS_MM_P");

        TH3D *Comp_MM_N_0_P = (TH3D*)Comp_MM_N_0->Clone("Comp_MM_N_0_P");
        TH3D *Comp_MM_N_1_P = (TH3D*)Comp_MM_N_1->Clone("Comp_MM_N_1_P");
        TH3D *Comp_Ph_N_0_P = (TH3D*)Comp_Ph_N_0->Clone("Comp_Ph_N_0_P");
        TH3D *Comp_Ph_N_1_P = (TH3D*)Comp_Ph_N_1->Clone("Comp_Ph_N_1_P");

        TH3D *Comp_MM_C_0_P = (TH3D*)Comp_MM_C_0->Clone("Comp_MM_C_0_P");
        TH3D *Comp_MM_C_1_P = (TH3D*)Comp_MM_C_1->Clone("Comp_MM_C_1_P");
        TH3D *Comp_Ph_C_0_P = (TH3D*)Comp_Ph_C_0->Clone("Comp_Ph_C_0_P");
        TH3D *Comp_Ph_C_1_P = (TH3D*)Comp_Ph_C_1->Clone("Comp_Ph_C_1_P");

        TH3D *Comp_MM_NN_0_P = (TH3D*)Comp_MM_NN_0->Clone("Comp_MM_NN_0_P");
        TH3D *Comp_MM_NN_1_P = (TH3D*)Comp_MM_NN_1->Clone("Comp_MM_NN_1_P");
        TH3D *Comp_Ph_NN_0_P = (TH3D*)Comp_Ph_NN_0->Clone("Comp_Ph_NN_0_P");
        TH3D *Comp_Ph_NN_1_P = (TH3D*)Comp_Ph_NN_1->Clone("Comp_Ph_NN_1_P");

        TH3D *Comp_MM_NC_0_P = (TH3D*)Comp_MM_NC_0->Clone("Comp_MM_NC_0_P");
        TH3D *Comp_MM_NC_1_P = (TH3D*)Comp_MM_NC_1->Clone("Comp_MM_NC_1_P");
        TH3D *Comp_Ph_NC_0_P = (TH3D*)Comp_Ph_NC_0->Clone("Comp_Ph_NC_0_P");
        TH3D *Comp_Ph_NC_1_P = (TH3D*)Comp_Ph_NC_1->Clone("Comp_Ph_NC_1_P");

        TH3D *Comp_MM_NT_0_P = (TH3D*)Comp_MM_NT_0->Clone("Comp_MM_NT_0_P");
        TH3D *Comp_MM_NT_1_P = (TH3D*)Comp_MM_NT_1->Clone("Comp_MM_NT_1_P");
        TH3D *Comp_Ph_NT_0_P = (TH3D*)Comp_Ph_NT_0->Clone("Comp_Ph_NT_0_P");
        TH3D *Comp_Ph_NT_1_P = (TH3D*)Comp_Ph_NT_1->Clone("Comp_Ph_NT_1_P");

        TH3D *Comp_MM_NW_0_P = (TH3D*)Comp_MM_NW_0->Clone("Comp_MM_NW_0_P");
        TH3D *Comp_MM_NW_1_P = (TH3D*)Comp_MM_NW_1->Clone("Comp_MM_NW_1_P");
        TH3D *Comp_Ph_NW_0_P = (TH3D*)Comp_Ph_NW_0->Clone("Comp_Ph_NW_0_P");
        TH3D *Comp_Ph_NW_1_P = (TH3D*)Comp_Ph_NW_1->Clone("Comp_Ph_NW_1_P");

        TH3D *Comp_MM_NNX_0_P = (TH3D*)Comp_MM_NNX_0->Clone("Comp_MM_NNX_0_P");
        TH3D *Comp_MM_NNX_1_P = (TH3D*)Comp_MM_NNX_1->Clone("Comp_MM_NNX_1_P");
        TH3D *Comp_Ph_NNX_0_P = (TH3D*)Comp_Ph_NNX_0->Clone("Comp_Ph_NNX_0_P");
        TH3D *Comp_Ph_NNX_1_P = (TH3D*)Comp_Ph_NNX_1->Clone("Comp_Ph_NNX_1_P");

        TH3D *Comp_MM_NCX_0_P = (TH3D*)Comp_MM_NCX_0->Clone("Comp_MM_NCX_0_P");
        TH3D *Comp_MM_NCX_1_P = (TH3D*)Comp_MM_NCX_1->Clone("Comp_MM_NCX_1_P");
        TH3D *Comp_Ph_NCX_0_P = (TH3D*)Comp_Ph_NCX_0->Clone("Comp_Ph_NCX_0_P");
        TH3D *Comp_Ph_NCX_1_P = (TH3D*)Comp_Ph_NCX_1->Clone("Comp_Ph_NCX_1_P");

        TH3D *Comp_MM_NTX_0_P = (TH3D*)Comp_MM_NTX_0->Clone("Comp_MM_NTX_0_P");
        TH3D *Comp_MM_NTX_1_P = (TH3D*)Comp_MM_NTX_1->Clone("Comp_MM_NTX_1_P");
        TH3D *Comp_Ph_NTX_0_P = (TH3D*)Comp_Ph_NTX_0->Clone("Comp_Ph_NTX_0_P");
        TH3D *Comp_Ph_NTX_1_P = (TH3D*)Comp_Ph_NTX_1->Clone("Comp_Ph_NTX_1_P");

        TH3D *Comp_MM_NWX_0_P = (TH3D*)Comp_MM_NWX_0->Clone("Comp_MM_NWX_0_P");
        TH3D *Comp_MM_NWX_1_P = (TH3D*)Comp_MM_NWX_1->Clone("Comp_MM_NWX_1_P");
        TH3D *Comp_Ph_NWX_0_P = (TH3D*)Comp_Ph_NWX_0->Clone("Comp_Ph_NWX_0_P");
        TH3D *Comp_Ph_NWX_1_P = (TH3D*)Comp_Ph_NWX_1->Clone("Comp_Ph_NWX_1_P");

        TaggHel0_P->Write();
        TaggHel1_P->Write();

        IncHel0_P->Write();
        IncHel1_P->Write();
        IncHits_P->Write();

        MM_CA_OA_P->Write();

        Pi0_IM_MM_P->Write();
        Pi0_IM_MM_NN_P->Write();

        Pi0_MM_P->Write();
        Pi0_MM_IM_P->Write();

        Pi0_MM_NN_P->Write();
        Pi0_MM_NC_P->Write();
        Pi0_MM_CC_P->Write();

        Pi0_MM_NNX_P->Write();
        Pi0_MM_NCX_P->Write();
        Pi0_MM_CCX_P->Write();

        Pi0_Re_All_P->Write();
        Pi0_Re_Det_P->Write();
        Pi0_Re_NoE_P->Write();

        Comp_MM_N_C_P->Write();

        Comp_CS_MM_P->Write();
        Reco_CS_MM_P->Write();

        Comp_MM_N_0_P->Write();

        Comp_MM_C_0_P->Write();

        Comp_MM_NN_0_P->Write();

        Comp_MM_NC_0_P->Write();

        Comp_MM_NT_0_P->Write();

        Comp_MM_NW_0_P->Write();

        Comp_MM_NNX_0_P->Write();

        Comp_MM_NCX_0_P->Write();

        Comp_MM_NTX_0_P->Write();

        Comp_MM_NWX_0_P->Write();

        if(cir_beam || lin_beam)
        {
            Comp_MM_N_1_P->Write();
            Comp_Ph_N_0_P->Write();
            Comp_Ph_N_1_P->Write();

            Comp_MM_C_1_P->Write();
            Comp_Ph_C_0_P->Write();
            Comp_Ph_C_1_P->Write();

            Comp_MM_NN_1_P->Write();
            Comp_Ph_NN_0_P->Write();
            Comp_Ph_NN_1_P->Write();

            Comp_MM_NC_1_P->Write();
            Comp_Ph_NC_0_P->Write();
            Comp_Ph_NC_1_P->Write();

            Comp_MM_NT_1_P->Write();
            Comp_Ph_NT_0_P->Write();
            Comp_Ph_NT_1_P->Write();

            Comp_MM_NW_1_P->Write();
            Comp_Ph_NW_0_P->Write();
            Comp_Ph_NW_1_P->Write();

            Comp_MM_NNX_1_P->Write();
            Comp_Ph_NNX_0_P->Write();
            Comp_Ph_NNX_1_P->Write();

            Comp_MM_NCX_1_P->Write();
            Comp_Ph_NCX_0_P->Write();
            Comp_Ph_NCX_1_P->Write();

            Comp_MM_NTX_1_P->Write();
            Comp_Ph_NTX_0_P->Write();
            Comp_Ph_NTX_1_P->Write();

            Comp_MM_NWX_1_P->Write();
            Comp_Ph_NWX_0_P->Write();
            Comp_Ph_NWX_1_P->Write();
        }

        delete TaggHel0_P;
        delete TaggHel1_P;

        delete IncHel0_P;
        delete IncHel1_P;
        delete IncHits_P;

        delete MM_CA_OA_P;

        delete Pi0_IM_MM_P;
        delete Pi0_IM_MM_NN_P;

        delete Pi0_MM_P;
        delete Pi0_MM_IM_P;

        delete Pi0_MM_NN_P;
        delete Pi0_MM_NC_P;
        delete Pi0_MM_CC_P;

        delete Pi0_MM_NNX_P;
        delete Pi0_MM_NCX_P;
        delete Pi0_MM_CCX_P;

        delete Pi0_Re_All_P;
        delete Pi0_Re_Det_P;
        delete Pi0_Re_NoE_P;

        delete Comp_MM_N_C_P;

        delete Comp_CS_MM_P;
        delete Reco_CS_MM_P;

        delete Comp_MM_N_0_P;
        delete Comp_MM_N_1_P;
        delete Comp_Ph_N_0_P;
        delete Comp_Ph_N_1_P;

        delete Comp_MM_C_0_P;
        delete Comp_MM_C_1_P;
        delete Comp_Ph_C_0_P;
        delete Comp_Ph_C_1_P;

        delete Comp_MM_NN_0_P;
        delete Comp_MM_NN_1_P;
        delete Comp_Ph_NN_0_P;
        delete Comp_Ph_NN_1_P;

        delete Comp_MM_NC_0_P;
        delete Comp_MM_NC_1_P;
        delete Comp_Ph_NC_0_P;
        delete Comp_Ph_NC_1_P;

        delete Comp_MM_NT_0_P;
        delete Comp_MM_NT_1_P;
        delete Comp_Ph_NT_0_P;
        delete Comp_Ph_NT_1_P;

        delete Comp_MM_NW_0_P;
        delete Comp_MM_NW_1_P;
        delete Comp_Ph_NW_0_P;
        delete Comp_Ph_NW_1_P;

        delete Comp_MM_NNX_0_P;
        delete Comp_MM_NNX_1_P;
        delete Comp_Ph_NNX_0_P;
        delete Comp_Ph_NNX_1_P;

        delete Comp_MM_NCX_0_P;
        delete Comp_MM_NCX_1_P;
        delete Comp_Ph_NCX_0_P;
        delete Comp_Ph_NCX_1_P;

        delete Comp_MM_NTX_0_P;
        delete Comp_MM_NTX_1_P;
        delete Comp_Ph_NTX_0_P;
        delete Comp_Ph_NTX_1_P;

        delete Comp_MM_NWX_0_P;
        delete Comp_MM_NWX_1_P;
        delete Comp_Ph_NWX_0_P;
        delete Comp_Ph_NWX_1_P;
    }
    if(!(IsMCFile()))
    {
        TaggHel0->Sumw2();
        TaggHel0->Add(TaggHel0_R,-ratio);
        TaggHel1->Sumw2();
        TaggHel1->Add(TaggHel1_R,-ratio);

        IncHel0->Sumw2();
        IncHel0->Add(IncHel0_R,-ratio);
        IncHel1->Sumw2();
        IncHel1->Add(IncHel1_R,-ratio);
        IncHits->Sumw2();
        IncHits->Add(IncHits_R,-ratio);

        MM_CA_OA->Sumw2();
        MM_CA_OA->Add(MM_CA_OA_R,-ratio);

        Pi0_IM_MM->Sumw2();
        Pi0_IM_MM->Add(Pi0_IM_MM_R,-ratio);
        Pi0_IM_MM_NN->Sumw2();
        Pi0_IM_MM_NN->Add(Pi0_IM_MM_NN_R,-ratio);

        Pi0_MM->Sumw2();
        Pi0_MM->Add(Pi0_MM_R,-ratio);
        Pi0_MM_IM->Sumw2();
        Pi0_MM_IM->Add(Pi0_MM_IM_R,-ratio);

        Pi0_MM_NN->Sumw2();
        Pi0_MM_NN->Add(Pi0_MM_NN_R,-ratio);
        Pi0_MM_NC->Sumw2();
        Pi0_MM_NC->Add(Pi0_MM_NC_R,-ratio);
        Pi0_MM_CC->Sumw2();
        Pi0_MM_CC->Add(Pi0_MM_CC_R,-ratio);

        Pi0_MM_NNX->Sumw2();
        Pi0_MM_NNX->Add(Pi0_MM_NNX_R,-ratio);
        Pi0_MM_NCX->Sumw2();
        Pi0_MM_NCX->Add(Pi0_MM_NCX_R,-ratio);
        Pi0_MM_CCX->Sumw2();
        Pi0_MM_CCX->Add(Pi0_MM_CCX_R,-ratio);

        Pi0_Re_All->Sumw2();
        Pi0_Re_All->Add(Pi0_Re_All_R,-ratio);
        Pi0_Re_Det->Sumw2();
        Pi0_Re_Det->Add(Pi0_Re_Det_R,-ratio);
        Pi0_Re_NoE->Sumw2();
        Pi0_Re_NoE->Add(Pi0_Re_NoE_R,-ratio);

        Comp_MM_N_C->Sumw2();
        Comp_MM_N_C->Add(Comp_MM_N_C_R,-ratio);

        Comp_CS_MM->Sumw2();
        Comp_CS_MM->Add(Comp_CS_MM_R,-ratio);
        Reco_CS_MM->Sumw2();
        Reco_CS_MM->Add(Reco_CS_MM_R,-ratio);

        Comp_MM_N_0->Sumw2();
        Comp_MM_N_0->Add(Comp_MM_N_0_R,-ratio);
        Comp_MM_N_1->Sumw2();
        Comp_MM_N_1->Add(Comp_MM_N_1_R,-ratio);
        Comp_Ph_N_0->Sumw2();
        Comp_Ph_N_0->Add(Comp_Ph_N_0_R,-ratio);
        Comp_Ph_N_1->Sumw2();
        Comp_Ph_N_1->Add(Comp_Ph_N_1_R,-ratio);

        Comp_MM_C_0->Sumw2();
        Comp_MM_C_0->Add(Comp_MM_C_0_R,-ratio);
        Comp_MM_C_1->Sumw2();
        Comp_MM_C_1->Add(Comp_MM_C_1_R,-ratio);
        Comp_Ph_C_0->Sumw2();
        Comp_Ph_C_0->Add(Comp_Ph_C_0_R,-ratio);
        Comp_Ph_C_1->Sumw2();
        Comp_Ph_C_1->Add(Comp_Ph_C_1_R,-ratio);

        Comp_MM_NN_0->Sumw2();
        Comp_MM_NN_0->Add(Comp_MM_NN_0_R,-ratio);
        Comp_MM_NN_1->Sumw2();
        Comp_MM_NN_1->Add(Comp_MM_NN_1_R,-ratio);
        Comp_Ph_NN_0->Sumw2();
        Comp_Ph_NN_0->Add(Comp_Ph_NN_0_R,-ratio);
        Comp_Ph_NN_1->Sumw2();
        Comp_Ph_NN_1->Add(Comp_Ph_NN_1_R,-ratio);

        Comp_MM_NC_0->Sumw2();
        Comp_MM_NC_0->Add(Comp_MM_NC_0_R,-ratio);
        Comp_MM_NC_1->Sumw2();
        Comp_MM_NC_1->Add(Comp_MM_NC_1_R,-ratio);
        Comp_Ph_NC_0->Sumw2();
        Comp_Ph_NC_0->Add(Comp_Ph_NC_0_R,-ratio);
        Comp_Ph_NC_1->Sumw2();
        Comp_Ph_NC_1->Add(Comp_Ph_NC_1_R,-ratio);

        Comp_MM_NT_0->Sumw2();
        Comp_MM_NT_0->Add(Comp_MM_NT_0_R,-ratio);
        Comp_MM_NT_1->Sumw2();
        Comp_MM_NT_1->Add(Comp_MM_NT_1_R,-ratio);
        Comp_Ph_NT_0->Sumw2();
        Comp_Ph_NT_0->Add(Comp_Ph_NT_0_R,-ratio);
        Comp_Ph_NT_1->Sumw2();
        Comp_Ph_NT_1->Add(Comp_Ph_NT_1_R,-ratio);

        Comp_MM_NW_0->Sumw2();
        Comp_MM_NW_0->Add(Comp_MM_NW_0_R,-ratio);
        Comp_MM_NW_1->Sumw2();
        Comp_MM_NW_1->Add(Comp_MM_NW_1_R,-ratio);
        Comp_Ph_NW_0->Sumw2();
        Comp_Ph_NW_0->Add(Comp_Ph_NW_0_R,-ratio);
        Comp_Ph_NW_1->Sumw2();
        Comp_Ph_NW_1->Add(Comp_Ph_NW_1_R,-ratio);

        Comp_MM_NNX_0->Sumw2();
        Comp_MM_NNX_0->Add(Comp_MM_NNX_0_R,-ratio);
        Comp_MM_NNX_1->Sumw2();
        Comp_MM_NNX_1->Add(Comp_MM_NNX_1_R,-ratio);
        Comp_Ph_NNX_0->Sumw2();
        Comp_Ph_NNX_0->Add(Comp_Ph_NNX_0_R,-ratio);
        Comp_Ph_NNX_1->Sumw2();
        Comp_Ph_NNX_1->Add(Comp_Ph_NNX_1_R,-ratio);

        Comp_MM_NCX_0->Sumw2();
        Comp_MM_NCX_0->Add(Comp_MM_NCX_0_R,-ratio);
        Comp_MM_NCX_1->Sumw2();
        Comp_MM_NCX_1->Add(Comp_MM_NCX_1_R,-ratio);
        Comp_Ph_NCX_0->Sumw2();
        Comp_Ph_NCX_0->Add(Comp_Ph_NCX_0_R,-ratio);
        Comp_Ph_NCX_1->Sumw2();
        Comp_Ph_NCX_1->Add(Comp_Ph_NCX_1_R,-ratio);

        Comp_MM_NTX_0->Sumw2();
        Comp_MM_NTX_0->Add(Comp_MM_NTX_0_R,-ratio);
        Comp_MM_NTX_1->Sumw2();
        Comp_MM_NTX_1->Add(Comp_MM_NTX_1_R,-ratio);
        Comp_Ph_NTX_0->Sumw2();
        Comp_Ph_NTX_0->Add(Comp_Ph_NTX_0_R,-ratio);
        Comp_Ph_NTX_1->Sumw2();
        Comp_Ph_NTX_1->Add(Comp_Ph_NTX_1_R,-ratio);

        Comp_MM_NWX_0->Sumw2();
        Comp_MM_NWX_0->Add(Comp_MM_NWX_0_R,-ratio);
        Comp_MM_NWX_1->Sumw2();
        Comp_MM_NWX_1->Add(Comp_MM_NWX_1_R,-ratio);
        Comp_Ph_NWX_0->Sumw2();
        Comp_Ph_NWX_0->Add(Comp_Ph_NWX_0_R,-ratio);
        Comp_Ph_NWX_1->Sumw2();
        Comp_Ph_NWX_1->Add(Comp_Ph_NWX_1_R,-ratio);
    }
    if(save_randoms && !(IsMCFile()))
    {
        TaggHel0_R->Write();
        TaggHel1_R->Write();

        IncHel0_R->Write();
        IncHel1_R->Write();
        IncHits_R->Write();

        MM_CA_OA_R->Write();

        Pi0_IM_MM_R->Write();
        Pi0_IM_MM_NN_R->Write();

        Pi0_MM_R->Write();
        Pi0_MM_IM_R->Write();

        Pi0_MM_NN_R->Write();
        Pi0_MM_NC_R->Write();
        Pi0_MM_CC_R->Write();

        Pi0_MM_NNX_R->Write();
        Pi0_MM_NCX_R->Write();
        Pi0_MM_CCX_R->Write();

        Pi0_Re_All_R->Write();
        Pi0_Re_Det_R->Write();
        Pi0_Re_NoE_R->Write();

        Comp_MM_N_C_R->Write();

        Comp_CS_MM_R->Write();
        Reco_CS_MM_R->Write();

        Comp_MM_N_0_R->Write();

        Comp_MM_C_0_R->Write();

        Comp_MM_NN_0_R->Write();

        Comp_MM_NC_0_R->Write();

        Comp_MM_NT_0_R->Write();

        Comp_MM_NW_0_R->Write();

        Comp_MM_NNX_0_R->Write();

        Comp_MM_NCX_0_R->Write();

        Comp_MM_NTX_0_R->Write();

        Comp_MM_NWX_0_R->Write();

        if(cir_beam || lin_beam)
        {
            Comp_MM_N_1_R->Write();
            Comp_Ph_N_0_R->Write();
            Comp_Ph_N_1_R->Write();

            Comp_MM_C_1_R->Write();
            Comp_Ph_C_0_R->Write();
            Comp_Ph_C_1_R->Write();

            Comp_MM_NN_1_R->Write();
            Comp_Ph_NN_0_R->Write();
            Comp_Ph_NN_1_R->Write();

            Comp_MM_NC_1_R->Write();
            Comp_Ph_NC_0_R->Write();
            Comp_Ph_NC_1_R->Write();

            Comp_MM_NT_1_R->Write();
            Comp_Ph_NT_0_R->Write();
            Comp_Ph_NT_1_R->Write();

            Comp_MM_NW_1_R->Write();
            Comp_Ph_NW_0_R->Write();
            Comp_Ph_NW_1_R->Write();

            Comp_MM_NNX_1_R->Write();
            Comp_Ph_NNX_0_R->Write();
            Comp_Ph_NNX_1_R->Write();

            Comp_MM_NCX_1_R->Write();
            Comp_Ph_NCX_0_R->Write();
            Comp_Ph_NCX_1_R->Write();

            Comp_MM_NTX_1_R->Write();
            Comp_Ph_NTX_0_R->Write();
            Comp_Ph_NTX_1_R->Write();

            Comp_MM_NWX_1_R->Write();
            Comp_Ph_NWX_0_R->Write();
            Comp_Ph_NWX_1_R->Write();
        }
    }

    delete TaggHel0_R;
    delete TaggHel1_R;

    delete IncHel0_R;
    delete IncHel1_R;
    delete IncHits_R;

    delete MM_CA_OA_R;

    delete Pi0_IM_MM_R;
    delete Pi0_IM_MM_NN_R;

    delete Pi0_MM_R;
    delete Pi0_MM_IM_R;

    delete Pi0_MM_NN_R;
    delete Pi0_MM_NC_R;
    delete Pi0_MM_CC_R;

    delete Pi0_MM_NNX_R;
    delete Pi0_MM_NCX_R;
    delete Pi0_MM_CCX_R;

    delete Pi0_Re_All_R;
    delete Pi0_Re_Det_R;
    delete Pi0_Re_NoE_R;

    delete Comp_MM_N_C_R;

    delete Comp_CS_MM_R;
    delete Reco_CS_MM_R;

    delete Comp_MM_N_0_R;
    delete Comp_MM_N_1_R;
    delete Comp_Ph_N_0_R;
    delete Comp_Ph_N_1_R;

    delete Comp_MM_C_0_R;
    delete Comp_MM_C_1_R;
    delete Comp_Ph_C_0_R;
    delete Comp_Ph_C_1_R;

    delete Comp_MM_NN_0_R;
    delete Comp_MM_NN_1_R;
    delete Comp_Ph_NN_0_R;
    delete Comp_Ph_NN_1_R;

    delete Comp_MM_NC_0_R;
    delete Comp_MM_NC_1_R;
    delete Comp_Ph_NC_0_R;
    delete Comp_Ph_NC_1_R;

    delete Comp_MM_NT_0_R;
    delete Comp_MM_NT_1_R;
    delete Comp_Ph_NT_0_R;
    delete Comp_Ph_NT_1_R;

    delete Comp_MM_NW_0_R;
    delete Comp_MM_NW_1_R;
    delete Comp_Ph_NW_0_R;
    delete Comp_Ph_NW_1_R;

    delete Comp_MM_NNX_0_R;
    delete Comp_MM_NNX_1_R;
    delete Comp_Ph_NNX_0_R;
    delete Comp_Ph_NNX_1_R;

    delete Comp_MM_NCX_0_R;
    delete Comp_MM_NCX_1_R;
    delete Comp_Ph_NCX_0_R;
    delete Comp_Ph_NCX_1_R;

    delete Comp_MM_NTX_0_R;
    delete Comp_MM_NTX_1_R;
    delete Comp_Ph_NTX_0_R;
    delete Comp_Ph_NTX_1_R;

    delete Comp_MM_NWX_0_R;
    delete Comp_MM_NWX_1_R;
    delete Comp_Ph_NWX_0_R;
    delete Comp_Ph_NWX_1_R;

    if(!cir_beam && !lin_beam)
    {
        delete Comp_MM_N_1;
        delete Comp_Ph_N_0;
        delete Comp_Ph_N_1;

        delete Comp_MM_C_1;
        delete Comp_Ph_C_0;
        delete Comp_Ph_C_1;

        delete Comp_MM_NN_1;
        delete Comp_Ph_NN_0;
        delete Comp_Ph_NN_1;

        delete Comp_MM_NC_1;
        delete Comp_Ph_NC_0;
        delete Comp_Ph_NC_1;

        delete Comp_MM_NT_1;
        delete Comp_Ph_NT_0;
        delete Comp_Ph_NT_1;

        delete Comp_MM_NW_1;
        delete Comp_Ph_NW_0;
        delete Comp_Ph_NW_1;

        delete Comp_MM_NNX_1;
        delete Comp_Ph_NNX_0;
        delete Comp_Ph_NNX_1;

        delete Comp_MM_NCX_1;
        delete Comp_Ph_NCX_0;
        delete Comp_Ph_NCX_1;

        delete Comp_MM_NTX_1;
        delete Comp_Ph_NTX_0;
        delete Comp_Ph_NTX_1;

        delete Comp_MM_NWX_1;
        delete Comp_Ph_NWX_0;
        delete Comp_Ph_NWX_1;
    }

    // Write all GH1's and TObjects defined in this class
    if(!(GTreeManager::Write())) return false;

    return true;
}
