#include "PPi0.h"

PPi0::PPi0()
{ 
    time 	        = new GH1("time",      "time",     1400, -700, 700);
    time_cut 	= new GH1("time_cut",  "time_cut", 1400, -700, 700);

    IM 		= new GH1("IM", 	"GoAT - #pi^{0}#rightarrowX;m_{#pi^{0}} (MeV)",                    100, 110, 160);

    IM_gg 	= new GH1("IM_gg", 	"GoAT - #pi^{0}#rightarrow#gamma+#gamma;m_{#pi^{0}} (MeV)",        100, 110, 160);
    IM_ggg 	= new GH1("IM_ggg", 	"GoAT - #pi^{0}#rightarrow#gamma+#gamma+#gamma;m_{#pi^{0}} (MeV)", 100, 110, 160);

    IM_rr 	= new GH1("IM_rr", 	"GoAT - #pi^{0}#rightarrowr+r;m_{#pi^{0}} (MeV)",                  100, 110, 160);
    IM_rrr 	= new GH1("IM_rrr", 	"GoAT - #pi^{0}#rightarrowr+r+r;m_{#pi^{0}} (MeV)",                100, 110, 160);

    IM_gr 	= new GH1("IM_gr", 	"GoAT - #pi^{0}#rightarrow#gamma+r;m_{#pi^{0}} (MeV)",             100, 110, 160);
    IM_ggr 	= new GH1("IM_ggr", 	"GoAT - #pi^{0}#rightarrow#gamma+#gamma+r;m_{#pi^{0}} (MeV)",      100, 110, 160);
    IM_grr 	= new GH1("IM_grr", 	"GoAT - #pi^{0}#rightarrow#gamma+r+r;m_{#pi^{0}} (MeV)",           100, 110, 160);

    IM_all 	= new TH3D("IM_all", 	"GoAT - #pi^{0}#rightarrowX;N_{#gamma};N_{r};m_{#pi^{0}} (MeV)",   10, 0, 10, 10, 0, 10, 100, 110, 160);

    MM 		= new GH1("MM", 	"GoAT - #pi^{0}#rightarrowX;m_{miss} (MeV)",                    400, 800, 1200);

    MM_gg 	= new GH1("MM_gg", 	"GoAT - #pi^{0}#rightarrow#gamma+#gamma;m_{miss} (MeV)",        400, 800, 1200);
    MM_ggg 	= new GH1("MM_ggg", 	"GoAT - #pi^{0}#rightarrow#gamma+#gamma+#gamma;m_{miss} (MeV)", 400, 800, 1200);

    MM_rr 	= new GH1("MM_rr", 	"GoAT - #pi^{0}#rightarrowr+r;m_{miss} (MeV)",                  400, 800, 1200);
    MM_rrr 	= new GH1("MM_rrr", 	"GoAT - #pi^{0}#rightarrowr+r+r;m_{miss} (MeV)",                400, 800, 1200);

    MM_gr 	= new GH1("MM_gr", 	"GoAT - #pi^{0}#rightarrow#gamma+r;m_{miss} (MeV)",             400, 800, 1200);
    MM_ggr 	= new GH1("MM_ggr", 	"GoAT - #pi^{0}#rightarrow#gamma+#gamma+r;m_{miss} (MeV)",      400, 800, 1200);
    MM_grr 	= new GH1("MM_grr", 	"GoAT - #pi^{0}#rightarrow#gamma+r+r;m_{miss} (MeV)",           400, 800, 1200);

    MM_IM 	= new GH2("MM_IM", 	"GoAT - #pi^{0}#rightarrow#gamma+r+r;m_{#pi^{0}} (MeV);m_{miss} (MeV)",          100, 110, 160, 400, 800, 1200);

    TaggerAccScal = new TH1D("TaggerAccScal","TaggerAccScal",352,0,352);
}

PPi0::~PPi0()
{
}

Bool_t	PPi0::Init()
{
    cout << "Initialising physics analysis..." << endl;
    cout << "--------------------------------------------------" << endl << endl;

    if(!InitBackgroundCuts()) return kFALSE;
    if(!InitTargetMass()) return kFALSE;
    if(!InitTaggerChannelCuts()) return kFALSE;
    if(!InitTaggerScalers()) return kFALSE;
    cout << "--------------------------------------------------" << endl;
    return kTRUE;
}

Bool_t	PPi0::Start()
{
    if(!IsGoATFile()){
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();

    TraverseValidEvents();

    return kTRUE;
}

void	PPi0::ProcessEvent()
{

    // Time diff (tagger - pi0)
    FillTime(*GetNeutralPions(),0,time);
    FillTimeCut(*GetNeutralPions(),0,time_cut);

    // Any decays
    FillMass(*GetNeutralPions(),0,IM);
    FillMissingMass(*GetNeutralPions(),0,MM);
    FillMassMissingMass(*GetNeutralPions(),0,MM_IM);
    
    // 2 photon decay
    if      ((GetNeutralPions()->GetNSubPhotons(0) == 2) && (GetNeutralPions()->GetNSubRootinos(0) == 0)){
        FillMass(*GetNeutralPions(),0,IM_gg);
        FillMissingMass(*GetNeutralPions(),0,MM_gg);
    }

    // 3 photon decay
    else if ((GetNeutralPions()->GetNSubPhotons(0) == 3) && (GetNeutralPions()->GetNSubRootinos(0) == 0)){
        FillMass(*GetNeutralPions(),0,IM_ggg);
        FillMissingMass(*GetNeutralPions(),0,MM_ggg);
    }

    // 2 rootino decay
    else if ((GetNeutralPions()->GetNSubPhotons(0) == 0) && (GetNeutralPions()->GetNSubRootinos(0) == 2)){
        FillMass(*GetNeutralPions(),0,IM_rr);
        FillMissingMass(*GetNeutralPions(),0,MM_rr);
    }

    // 3 rootino decay
    else if ((GetNeutralPions()->GetNSubPhotons(0) == 0) && (GetNeutralPions()->GetNSubRootinos(0) == 3)){
        FillMass(*GetNeutralPions(),0,IM_rrr);
        FillMissingMass(*GetNeutralPions(),0,MM_rrr);
    }

    // 1 photon and 1 rootino decay
    else if ((GetNeutralPions()->GetNSubPhotons(0) == 1) && (GetNeutralPions()->GetNSubRootinos(0) == 1)){
        FillMass(*GetNeutralPions(),0,IM_gr);
        FillMissingMass(*GetNeutralPions(),0,MM_gr);
    }

    // 2 photon and 1 rootino decay
    else if ((GetNeutralPions()->GetNSubPhotons(0) == 2) && (GetNeutralPions()->GetNSubRootinos(0) == 1)){
        FillMass(*GetNeutralPions(),0,IM_ggr);
        FillMissingMass(*GetNeutralPions(),0,MM_ggr);
    }

    // 1 photon and 2 rootino decay
    else if ((GetNeutralPions()->GetNSubPhotons(0) == 1) && (GetNeutralPions()->GetNSubRootinos(0) == 2)){
        FillMass(*GetNeutralPions(),0,IM_grr);
        FillMissingMass(*GetNeutralPions(),0,MM_grr);
    }

    // 3D histogram showing invariant mass for different photon rootino decays
    IM_all->Fill(GetNeutralPions()->GetNSubPhotons(0),GetNeutralPions()->GetNSubRootinos(0),GetNeutralPions()->Particle(0).M());

}  

void	PPi0::ProcessScalerRead()
{
    // Fill Tagger Scalers
    FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal);
}

Bool_t	PPi0::Write()
{
    // Write all GH1's and TObjects defined in this class
    GTreeManager::Write();
}

void PPi0::FillMassMissingMass(const GTreeParticle& tree, GH2* gHist, Bool_t TaggerBinning)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
    {
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
        {
            FillMassMissingMass(tree, i, j, gHist, TaggerBinning);
        }
    }
}

void PPi0::FillMassMissingMass(const GTreeParticle& tree, Int_t particle_index, GH2* gHist, Bool_t TaggerBinning)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
    {
        FillMassMissingMass(tree, particle_index, i, gHist, TaggerBinning);
    }
}

void PPi0::FillMassMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH2* gHist, Bool_t TaggerBinning)
{
    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < GetTC_cut_min()) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > GetTC_cut_max()) return;

    // calc particle time diff
    Double_t time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);

    // calc missing p4
    TLorentzVector missingp4 = CalcMissingP4(tree, particle_index,tagger_index);

    // Fill GH1
    if(TaggerBinning)   gHist->Fill(tree.GetMass(particle_index), missingp4.M(),time, GetTagger()->GetTaggedChannel(tagger_index));
    else gHist->Fill(tree.GetMass(particle_index), missingp4.M(),time);
}
