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
  
  TaggerAccScal = new TH1D("TaggerAccScal","TaggerAccScal",352,0,352);
}

PPi0::~PPi0()
{
}

Bool_t	PPi0::Init(const char* configfile)
{
  cout << "Initialising physics analysis..." << endl;
  cout << "--------------------------------------------------" << endl << endl;
  if(configfile) SetConfigFile(configfile);

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
  FillTime(*neutralPions,0,time);
  FillTimeCut(*neutralPions,0,time_cut);
  
  // Any decays
  FillMass(*neutralPions,0,IM);
  FillMissingMass(*neutralPions,0,MM);
    
  // 2 photon decay
  if      ((neutralPions->GetNSubPhotons(0) == 2) && (neutralPions->GetNSubRootinos(0) == 0)){
    FillMass(*neutralPions,0,IM_gg);
    FillMissingMass(*neutralPions,0,MM_gg);
  }

  // 3 photon decay
  else if ((neutralPions->GetNSubPhotons(0) == 3) && (neutralPions->GetNSubRootinos(0) == 0)){
    FillMass(*neutralPions,0,IM_ggg);
    FillMissingMass(*neutralPions,0,MM_ggg);
  }

  // 2 rootino decay
  else if ((neutralPions->GetNSubPhotons(0) == 0) && (neutralPions->GetNSubRootinos(0) == 2)){
    FillMass(*neutralPions,0,IM_rr);
    FillMissingMass(*neutralPions,0,MM_rr);
  }

  // 3 rootino decay
  else if ((neutralPions->GetNSubPhotons(0) == 0) && (neutralPions->GetNSubRootinos(0) == 3)){
    FillMass(*neutralPions,0,IM_rrr);
    FillMissingMass(*neutralPions,0,MM_rrr);
  }

  // 1 photon and 1 rootino decay
  else if ((neutralPions->GetNSubPhotons(0) == 1) && (neutralPions->GetNSubRootinos(0) == 1)){
    FillMass(*neutralPions,0,IM_gr);
    FillMissingMass(*neutralPions,0,MM_gr);
  }

  // 2 photon and 1 rootino decay
  else if ((neutralPions->GetNSubPhotons(0) == 2) && (neutralPions->GetNSubRootinos(0) == 1)){
    FillMass(*neutralPions,0,IM_ggr);
    FillMissingMass(*neutralPions,0,MM_ggr);
  }

  // 1 photon and 2 rootino decay
  else if ((neutralPions->GetNSubPhotons(0) == 1) && (neutralPions->GetNSubRootinos(0) == 2)){
    FillMass(*neutralPions,0,IM_grr);
    FillMissingMass(*neutralPions,0,MM_grr);
  }

  // 3D histogram showing invariant mass for different photon rootino decays
  IM_all->Fill(neutralPions->GetNSubPhotons(0),neutralPions->GetNSubRootinos(0),neutralPions->Particle(0).M());

}  

void	PPi0::ProcessScalerRead()
{
  // Fill Tagger Scalers
  FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal);
}

Bool_t	PPi0::Write()
{
  // Write some TH1s
  GTreeManager::Write(TaggerAccScal);
  GTreeManager::Write(IM_all);

  // Write all GH1's easily
  GTreeManager::Write();
}
