#include "Pi0X.h"

Pi0X::Pi0X()
{ 
  helicity          = new GH1("helicity",          "helicity",    2, 0, 2);
  helicityZE        = new GH1("helicityZE",        "helicityZE",  2, 0, 2);
  helerrors         = new GH1("helicityerr",       "helicityerr", 15, 0, 15);
  errcode           = new GH1("errcode",           "errcode",     11, 0, 11);

  time 	            = new GH1("time",              "time",       1400, -700, 700);
  time_all          = new GH1("time_all",          "time_all",   1400, -700, 700);
  time_cut          = new GH1("time_cut",          "time_cut",   1400, -700, 700);

  egamma            = new GH1("egamma",            "GOAT - E_{#gamma}",               1000, 0, 1600);
  egamma_all        = new GH1("egamma_all",        "GoAT - E_{#gamma} all tracks",    1000, 0, 1600);

  FPD               = new GH1("FPD",               "GoAT - FPD hits (p-r)",                            352, 0, 352);
  FPD_hel0          = new GH1("FPD_hel0",          "GoAT - FPD hits (p-r) helicity 0",                 352, 0, 352);
  FPD_CB_hel0       = new GH1("FPD_CB_hel0",       "GoAT - FPD hits (p-r) helicity 0, CB trigger",     352, 0, 352);
  FPD_TAPS_hel0     = new GH1("FPD_TAPS_hel0",     "GoAT - FPD hits (p-r) helicity 0, TAPS trigger",   352, 0, 352);
  FPD_hel1          = new GH1("FPD_hel1",          "GoAT - FPD hits (p-r) helicity 1",                 352, 0, 352);
  FPD_CB_hel1       = new GH1("FPD_CB_hel1",       "GoAT - FPD hits (p-r) helicity 1, CB trigger",     352, 0, 352);
  FPD_TAPS_hel1     = new GH1("FPD_TAPS_hel1",     "GoAT - FPD hits (p-r) helicity 1, TAPS trigger",   352, 0, 352);
  FPD_all           = new GH1("FPD_all",           "GoAT - FPD hits all tracks (p-r)",                 352, 0, 352);
  FPD_all_hel0      = new GH1("FPD_all_hel0",      "GoAT - FPD hits all tracks (p-r) hel 0",           352, 0, 352);
  FPD_all_hel1      = new GH1("FPD_all_hel1",      "GoAT - FPD hits all tracks (p-r) hel 1",           352, 0, 352);

  phi               = new GH1("phi",               "GOAT - #phi distribution",                         360, 0, 360);
  thetaphi          = new GH2("thetaphi",          "GOAT - #theta vs #phi prompt",                     180, 0, 180, 360, 0, 360);
  phi_all           = new GH1("phi_all",           "GOAT - #phi distribution all tracks",              360, 0, 360);
  thetaphi_all      = new GH2("thetaphi_all",      "GOAT - #theta vs #phi prompt all tracks",          180, 0, 180, 360, 0, 360);
  theta             = new GH1("theta",             "GOAT - #theta distribution",                       180, 0, 180);
  theta_all         = new GH1("theta_all",         "GOAT - #theta distribution all tracks",            180, 0, 180);
  theta_hel0        = new GH1("theta_hel0",        "GOAT - #theta distribution - hel0",                180, 0, 180);
  theta_hel1        = new GH1("theta_hel1",        "GOAT - #theta distribution - hel1",                180, 0, 180);
  theta_hel0CM      = new GH1("theta_hel0CM",      "GOAT - #theta CM distribution - hel0",             180, 0, 180);
  theta_hel1CM      = new GH1("theta_hel1CM",      "GOAT - #theta CM distribution - hel1",             180, 0, 180);
  thetaCB_hel0      = new GH1("thetaCB_hel0",      "GOAT - #theta distribution - hel0, CB trigg",      180, 0, 180);
  thetaCB_hel1      = new GH1("thetaCB_hel1",      "GOAT - #theta distribution - hel1, CB trigg",      180, 0, 180);
  thetaCB_hel0CM    = new GH1("thetaCB_hel0CM",    "GOAT - #theta CM distribution - hel0, CB trigg",   180, 0, 180);
  thetaCB_hel1CM    = new GH1("thetaCB_hel1CM",    "GOAT - #theta CM distribution - hel1, CB trigg",   180, 0, 180);
  thetaTAPS_hel0    = new GH1("thetaTAPS_hel0",    "GOAT - #theta distribution - hel0, TAPS trigg",    180, 0, 180);
  thetaTAPS_hel1    = new GH1("thetaTAPS_hel1",    "GOAT - #theta distribution - hel1, TAPS trigg",    180, 0, 180);
  thetaTAPS_hel0CM  = new GH1("thetaTAPS_hel0CM",  "GOAT - #theta CM distribution - hel0, TAPS trigg", 180, 0, 180);
  thetaTAPS_hel1CM  = new GH1("thetaTAPS_hel1CM",  "GOAT - #theta CM distribution - hel1, TAPS trigg", 180, 0, 180);
  theta_all_hel0    = new GH1("theta_all_hel0",    "GOAT - #theta distribution all tracks - hel0",     180, 0, 180);
  theta_all_hel1    = new GH1("theta_all_hel1",    "GOAT - #theta distribution all tracks - hel1",     180, 0, 180);

  Cphi               = new GH1("Cphi",               "GOAT + cond - #phi distribution",                         360, 0, 360);
  Cthetaphi          = new GH2("Cthetaphi",          "GOAT + cond - #theta vs #phi prompt",                     180, 0, 180, 360, 0, 360);
  Cphi_all           = new GH1("Cphi_all",           "GOAT + cond - #phi distribution all tracks",              360, 0, 360);
  Cthetaphi_all      = new GH2("Cthetaphi_all",      "GOAT + cond - #theta vs #phi prompt all tracks",          180, 0, 180, 360, 0, 360);
  Ctheta             = new GH1("Ctheta",             "GOAT + cond - #theta distribution",                       180, 0, 180);
  Ctheta_all         = new GH1("Ctheta_all",         "GOAT + cond - #theta distribution all tracks",            180, 0, 180);
  Ctheta_hel0        = new GH1("Ctheta_hel0",        "GOAT + cond - #theta distribution - hel0",                180, 0, 180);
  Ctheta_hel1        = new GH1("Ctheta_hel1",        "GOAT + cond - #theta distribution - hel1",                180, 0, 180);
  Ctheta_hel0CM      = new GH1("Ctheta_hel0CM",      "GOAT + cond - #theta CM distribution - hel0",             180, 0, 180);
  Ctheta_hel1CM      = new GH1("Ctheta_hel1CM",      "GOAT + cond - #theta CM distribution - hel1",             180, 0, 180);
  CthetaCB_hel0      = new GH1("CthetaCB_hel0",      "GOAT + cond - #theta distribution - hel0, CB trigg",      180, 0, 180);
  CthetaCB_hel1      = new GH1("CthetaCB_hel1",      "GOAT + cond - #theta distribution - hel1, CB trigg",      180, 0, 180);
  CthetaCB_hel0CM    = new GH1("CthetaCB_hel0CM",    "GOAT + cond - #theta CM distribution - hel0, CB trigg",   180, 0, 180);
  CthetaCB_hel1CM    = new GH1("CthetaCB_hel1CM",    "GOAT + cond - #theta CM distribution - hel1, CB trigg",   180, 0, 180);
  CthetaTAPS_hel0    = new GH1("CthetaTAPS_hel0",    "GOAT + cond - #theta distribution - hel0, TAPS trigg",    180, 0, 180);
  CthetaTAPS_hel1    = new GH1("CthetaTAPS_hel1",    "GOAT + cond - #theta distribution - hel1, TAPS trigg",    180, 0, 180);
  CthetaTAPS_hel0CM  = new GH1("CthetaTAPS_hel0CM",  "GOAT + cond - #theta CM distribution - hel0, TAPS trigg", 180, 0, 180);
  CthetaTAPS_hel1CM  = new GH1("CthetaTAPS_hel1CM",  "GOAT + cond - #theta CM distribution - hel1, TAPS trigg", 180, 0, 180);
  Ctheta_all_hel0    = new GH1("Ctheta_all_hel0",    "GOAT + cond - #theta distribution all tracks - hel0",     180, 0, 180);
  Ctheta_all_hel1    = new GH1("Ctheta_all_hel1",    "GOAT + cond - #theta distribution all tracks - hel1",     180, 0, 180);

  IM 	  	  = new GH1("IM", 	       "GoAT - #pi^{0}#rightarrowX;m_{#pi^{0}} (MeV)",                    100, 110, 160);
  IM_gg   	  = new GH1("IM_gg", 	       "GoAT - #pi^{0}#rightarrow#gamma+#gamma;m_{#pi^{0}} (MeV)",        100, 110, 160);
  IM_ggg 	  = new GH1("IM_ggg", 	       "GoAT - #pi^{0}#rightarrow#gamma+#gamma+#gamma;m_{#pi^{0}} (MeV)", 100, 110, 160);
  CIM_gg   	  = new GH1("CIM_gg", 	       "GoAT + cond - #pi^{0}#rightarrow#gamma+#gamma;m_{#pi^{0}} (MeV)", 100, 110, 160);

  MM 	  	  = new GH1("MM", 	       "GoAT - #pi^{0}#rightarrowX;m_{miss} (MeV)",                    400, 600, 2000);
  MM_gg 	  = new GH1("MM_gg", 	       "GoAT - #pi^{0}#rightarrow#gamma+#gamma;m_{miss} (MeV)",        400, 600, 2000);
  MM_ggg 	  = new GH1("MM_ggg", 	       "GoAT - #pi^{0}#rightarrow#gamma+#gamma+#gamma;m_{miss} (MeV)", 400, 600, 2000);
  CMM_gg 	  = new GH1("CMM_gg", 	       "GoAT + cond - #pi^{0}#rightarrow#gamma+#gamma;m_{miss} (MeV)", 400, 600, 2000);

  TaggerAccScal   = new TH1D("TaggerAccScal",   "TaggerAccScal", 352, 0, 352);

}

Pi0X::~Pi0X()
{
}

Bool_t Pi0X::Init()
{
  cout << "Initialising physics analysis..." << endl;
  cout << "--------------------------------------------------" << endl << endl;

  if(!InitBackgroundCuts()) return kFALSE;
  if(!InitTargetMass()) return kFALSE;
  if(!InitTaggerChannelCuts()) return kFALSE;
  if(!InitTaggerScalers()) return kFALSE;
  cout << "--------------------------------------------------" << endl;

  evtNum = 0;
  return kTRUE;
}

Bool_t Pi0X::Start()
{
  cout << "=== Pi0X::Start() " << endl;
  if(!IsGoATFile()){
    cout << "ERROR: Input File is not a GoAT file." << endl;
    return kFALSE;
  }
  SetAsPhysicsFile();

  TraverseValidEvents();

  return kTRUE;
}

void Pi0X::ProcessEvent()
{
  //  cout << "===> NEW event " << endl;
  // Time diff (tagger - pi0)
  FillTime(*GetNeutralPions(),0,time);
  FillTime(*GetNeutralPions(), time_all);
  FillTimeCut(*GetNeutralPions(),0,time_cut);

  FillHelicity(helicityZE, helerrors, errcode, helicity);

  // Photon energy
  FillPhotonEnergy(*GetNeutralPions(), 0, egamma);
  FillPhotonEnergy(*GetNeutralPions(), egamma_all);

  // FPD hits
  FillFPD(*GetNeutralPions(), 0, FPD);
  FillFPD(*GetNeutralPions(), FPD_all);
  FillFPD(*GetNeutralPions(), 0, FPD_hel0, FPD_CB_hel0, FPD_TAPS_hel0, FPD_hel1, FPD_CB_hel1, FPD_TAPS_hel1);
  FillFPD(*GetNeutralPions(), FPD_all_hel0, FPD_all_hel1);

  // Any decays
  FillMass(*GetNeutralPions(),0,IM);
  FillMissingMass(*GetNeutralPions(),0,MM, kTRUE);

  // 2 photon decay
  //  if ((GetNeutralPions()->GetNSubPhotons(0)==2) && (GetNeutralPions()->GetNSubRootinos(0) == 0)) 
  //  if ((GetNeutralPions()->GetNSubPhotons(0)==2))
    {
      // Angular distributions
      FillAngularDist(*GetNeutralPions(), 0, phi, thetaphi);
      FillAngularDist(*GetNeutralPions(), phi_all, thetaphi_all);
      FillTheta(*GetNeutralPions(), 0, theta, kTRUE);
      //    FillTheta(*GetNeutralPions(), 0, theta_hel0, thetaCB_hel0, thetaTAPS_hel0, theta_hel1, thetaCB_hel1, thetaTAPS_hel1, kTRUE);
      FillTheta(*GetNeutralPions(), 0, 
		theta_hel0, thetaCB_hel0, thetaTAPS_hel0, theta_hel0CM, thetaCB_hel0CM, thetaTAPS_hel0CM,
		theta_hel1, thetaCB_hel1, thetaTAPS_hel1, theta_hel1CM, thetaCB_hel1CM, thetaTAPS_hel1CM,
		kTRUE);
      FillTheta(*GetNeutralPions(), theta_all, kTRUE);
      FillTheta(*GetNeutralPions(), theta_all_hel0, theta_all_hel1, kTRUE);
      
      // mass and missing mass
      FillMass(*GetNeutralPions(),0,IM_gg);
      FillMissingMass(*GetNeutralPions(),0,MM_gg, kTRUE);
    }
  if (GetNeutralPions()->GetNSubPhotons(0) == GetNeutralPions()->GetNSubParticles(0) )
    {
      // Angular distributions
      FillAngularDist(*GetNeutralPions(), 0, Cphi, Cthetaphi);
      FillAngularDist(*GetNeutralPions(), Cphi_all, Cthetaphi_all);
      FillTheta(*GetNeutralPions(), 0, Ctheta, kTRUE);
      FillTheta(*GetNeutralPions(), 0, 
		Ctheta_hel0, CthetaCB_hel0, CthetaTAPS_hel0, Ctheta_hel0CM, CthetaCB_hel0CM, CthetaTAPS_hel0CM,
		Ctheta_hel1, CthetaCB_hel1, CthetaTAPS_hel1, Ctheta_hel1CM, CthetaCB_hel1CM, CthetaTAPS_hel1CM,
		kTRUE);
      FillTheta(*GetNeutralPions(), Ctheta_all, kTRUE);
      FillTheta(*GetNeutralPions(), Ctheta_all_hel0, Ctheta_all_hel1, kTRUE);
      
      // mass and missing mass
      FillMass(*GetNeutralPions(),0,CIM_gg);
      FillMissingMass(*GetNeutralPions(),0,CMM_gg, kTRUE);
    }
  
//   // 3 photon decay
//   else if ((GetNeutralPions()->GetNSubPhotons(0) == 3) && (GetNeutralPions()->GetNSubRootinos(0) == 0)){
//     FillMass(*GetNeutralPions(),0,IM_ggg);
//     FillMissingMass(*GetNeutralPions(),0,MM_ggg, kTRUE);
//   }

  evtNum++;
}  

void	Pi0X::ProcessScalerRead()
{
  // Fill Tagger Scalers
  FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal);
}

Bool_t	Pi0X::Write()
{
  cout << "Pi0X::Write() " << endl;
  // Write all GH1's and TObjects defined in this class
  return GTreeManager::Write();
}


void Pi0X::FillPhotonEnergy(const GTreeParticle& tree, GH1* gHist)
{
  for (Int_t i = 0; i < tree.GetNParticles(); i++)
    {
      for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
	  // Is tagger channel rejected by user?
	  if(GetTagger()->GetTaggedChannel(j) < GetTC_cut_min()) continue;
	  if(GetTagger()->GetTaggedChannel(j) > GetTC_cut_max()) continue;

	  // calc particle time diff
	  Double_t time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);

	  Double_t gammae = GetTagger()->GetTaggedEnergy(j);
	  gHist->Fill(gammae, time);
	}
    }
}


void Pi0X::FillPhotonEnergy(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
  for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
    {
      // Is tagger channel rejected by user?
      if(GetTagger()->GetTaggedChannel(j) < GetTC_cut_min()) continue;
      if(GetTagger()->GetTaggedChannel(j) > GetTC_cut_max()) continue;

      Double_t time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
      Double_t gammae = GetTagger()->GetTaggedEnergy(j);
      gHist->Fill(gammae, time);
    }
}


void Pi0X::FillHelicity(GH1* gHist, GH1* gHist2, GH1* gHist3, GH1* gHist4)
{
  Int_t nerror = 0;
  Bool_t hel;
  Bool_t ErrFlag = kFALSE;

  nerror = GetTrigger()->GetNErrors();
  gHist2->Fill(nerror);

  hel = GetTrigger()->GetHelicity();
  if (nerror==0) {
    if (hel==kFALSE)
      gHist->Fill(0);
    else if (hel==kTRUE)
      gHist->Fill(1);
  }

  const Int_t* errcode;
  for (int i=0; i<nerror; i++) {
    errcode = GetTrigger()->GetErrorCode();
    gHist3->Fill(*errcode);
    if (*errcode==9 || *errcode==10) {
      ErrFlag = kTRUE;
      break;
    }
  }
  if (!ErrFlag && nerror>0) {
    if (hel==kFALSE)
      gHist4->Fill(0);
    else if (hel==kTRUE)
      gHist4->Fill(1);
  }
}


void Pi0X::FillAngularDist(const GTreeParticle& tree, GH1* hHist, GH2* ghHist)
{
  // hHist   --> phi hist
  // ghHist  --> 2D theta-phi prompt hist
  Int_t nerror = GetTrigger()->GetNErrors();
  for (Int_t i = 0; i < tree.GetNParticles(); i++)
    {
      for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
	  // Is tagger channel rejected by user?
	  if(GetTagger()->GetTaggedChannel(j) < GetTC_cut_min()) continue;
	  if(GetTagger()->GetTaggedChannel(j) > GetTC_cut_max()) continue;
	  
	  // calc particle time diff
	  Double_t time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);

	  // Fill theta/phi histograms 
	  if (nerror==0) {
	    hHist->Fill((TVector2::Phi_0_2pi(tree.GetPhiRad(i)))*TMath::RadToDeg(), time);
	    ghHist->Fill(tree.GetTheta(i), (TVector2::Phi_0_2pi(tree.GetPhiRad(i)))*TMath::RadToDeg(), time);
	  }
	}
    }
}


void Pi0X::FillAngularDist(const GTreeParticle& tree, Int_t particle_index, GH1* hHist, GH2* ghHist)
{
  // hHist   --> phi hist
  // ghHist  --> 2D theta-phi prompt hist
  Int_t nerror = GetTrigger()->GetNErrors();

  for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
    {
      // Is tagger channel rejected by user?
      if(GetTagger()->GetTaggedChannel(j) < GetTC_cut_min()) continue;
      if(GetTagger()->GetTaggedChannel(j) > GetTC_cut_max()) continue;
      
      // calc particle time diff
      Double_t time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
      
      // Fill theta/phi histograms 
      if (nerror==0) {
	hHist->Fill((TVector2::Phi_0_2pi(tree.GetPhiRad(particle_index)))*TMath::RadToDeg(), time);
	ghHist->Fill(tree.GetTheta(particle_index), (TVector2::Phi_0_2pi(tree.GetPhiRad(particle_index)))*TMath::RadToDeg(), time);
      }
    }
}


void Pi0X::FillTheta(const GTreeParticle& tree, GH1* gHist, Bool_t TaggerBinning)
{
  Int_t nerror = GetTrigger()->GetNErrors();

  for (Int_t i = 0; i < tree.GetNParticles(); i++)
    {
      for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
	  // Is tagger channel rejected by user?
	  if(GetTagger()->GetTaggedChannel(j) < GetTC_cut_min()) continue;
	  if(GetTagger()->GetTaggedChannel(j) > GetTC_cut_max()) continue;
	  
	  // calc particle time diff
	  Double_t time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);

	  // Fill theta/phi histograms 
	  if (nerror==0) {
	    if (TaggerBinning) gHist->Fill(tree.GetTheta(i), time,  GetTagger()->GetTaggedChannel(j));
	    else gHist->Fill(tree.GetTheta(i), time);
	  }
	}
    }
}


void Pi0X::FillTheta(const GTreeParticle& tree, GH1* gHist0, GH1* gHist1, Bool_t TaggerBinning)
{
  Int_t nerror = GetTrigger()->GetNErrors();
  Bool_t helicity = GetTrigger()->GetHelicity();

  for (Int_t i = 0; i < tree.GetNParticles(); i++)
    {
      for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
	  // Is tagger channel rejected by user?
	  if(GetTagger()->GetTaggedChannel(j) < GetTC_cut_min()) continue;
	  if(GetTagger()->GetTaggedChannel(j) > GetTC_cut_max()) continue;
	  
	  // calc particle time diff
	  Double_t time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);

	  // Fill theta/phi histograms 
	  if (nerror==0 && helicity==kFALSE) {
	    if (TaggerBinning)
	      gHist0->Fill(tree.GetTheta(i), time, GetTagger()->GetTaggedChannel(j));
	    else 
	      gHist0->Fill(tree.GetTheta(i), time);
	  }
	  else if (nerror==0 && helicity==kTRUE) {
	    if (TaggerBinning)
	      gHist1->Fill(tree.GetTheta(i), time,  GetTagger()->GetTaggedChannel(j));
	    else 
	      gHist1->Fill(tree.GetTheta(i), time);
	  }
	}
    }
}


void Pi0X::FillTheta(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning)
{
  Int_t nerror = GetTrigger()->GetNErrors();
  for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
    {
      // Is tagger channel rejected by user?
      if(GetTagger()->GetTaggedChannel(j) < GetTC_cut_min()) continue;
      if(GetTagger()->GetTaggedChannel(j) > GetTC_cut_max()) continue;
      
      // calc particle time diff
      Double_t time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
      
      // Fill theta/phi histograms
      if (nerror==0) {
	if (TaggerBinning) 
	  gHist->Fill(tree.GetTheta(particle_index), time, GetTagger()->GetTaggedChannel(j));
	else
	  gHist->Fill(tree.GetTheta(particle_index), time);
      }
    }
}


void Pi0X::FillTheta(const GTreeParticle& tree, Int_t particle_index, GH1* gHist0, GH1* gHist1, Bool_t TaggerBinning)
{
  Int_t nerror = GetTrigger()->GetNErrors();
  Bool_t helicity = GetTrigger()->GetHelicity();

  for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
    {
      // Is tagger channel rejected by user?
      if(GetTagger()->GetTaggedChannel(j) < GetTC_cut_min()) continue;
      if(GetTagger()->GetTaggedChannel(j) > GetTC_cut_max()) continue;
      
      // calc particle time diff
      Double_t time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
      
      // Fill theta/phi histograms 
      if (nerror==0 && helicity==kFALSE) {
	if (TaggerBinning)
	  gHist0->Fill(tree.GetTheta(particle_index), time, GetTagger()->GetTaggedChannel(j));
	else
	  gHist0->Fill(tree.GetTheta(particle_index), time);
      }
      else if (nerror==0 && helicity==kTRUE) {
	if (TaggerBinning)
	  gHist1->Fill(tree.GetTheta(particle_index), time, GetTagger()->GetTaggedChannel(j));
	else
	  gHist1->Fill(tree.GetTheta(particle_index), time);
      }
    }
}


void Pi0X::FillTheta(const GTreeParticle& tree, Int_t particle_index, GH1* gHist0, GH1* gCBHist0, GH1* gTAPSHist0, GH1* gHist0CM, GH1* gCBHist0CM, GH1* gTAPSHist0CM, GH1* gHist1, GH1* gCBHist1, GH1* gTAPSHist1, GH1* gHist1CM, GH1* gCBHist1CM, GH1* gTAPSHist1CM, Bool_t TaggerBinning)
{
  Int_t nerror = GetTrigger()->GetNErrors();
  Bool_t helicity = GetTrigger()->GetHelicity();

  const Int_t *tp = GetTrigger()->GetTriggerPattern();

  Double_t mom      = TMath::Sqrt(tree.GetTotalEnergy(particle_index)*tree.GetTotalEnergy(particle_index) - tree.GetMass(particle_index)*tree.GetMass(particle_index));

  for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
    {
      // Is tagger channel rejected by user?
      if(GetTagger()->GetTaggedChannel(j) < GetTC_cut_min()) continue;
      if(GetTagger()->GetTaggedChannel(j) > GetTC_cut_max()) continue;
      
      // calc particle time diff
      Double_t time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
      Double_t theta_cm = Compute_ThetaCM(GetTagger()->GetTaggedEnergy(j), mom, tree.GetTotalEnergy(particle_index), tree.GetThetaRad(particle_index), TVector2::Phi_0_2pi(tree.GetPhiRad(particle_index)));
      //     cout << "egamma " << GetTagger()->GetTaggedEnergy(j) << ", mom " << mom << ", ene " << tree.GetTotalEnergy(particle_index) << ", theta " << tree.GetTheta(particle_index) << ", phi " << TVector2::Phi_0_2pi(tree.GetPhiRad(particle_index))*TMath::RadToDeg() << ", thetacm " << theta_cm << endl;
      // Fill theta/phi histograms 

      // all triggers
      if (nerror==0 && helicity==kFALSE) {
	if (TaggerBinning) {
	  gHist0->Fill(tree.GetTheta(particle_index), time, GetTagger()->GetTaggedChannel(j));
	  gHist0CM->Fill(theta_cm, time, GetTagger()->GetTaggedChannel(j));
	}
	else {
	  gHist0->Fill(tree.GetTheta(particle_index), time);
	  gHist0CM->Fill(theta_cm, time);
	}
      }
      else if (nerror==0 && helicity==kTRUE) {
	if (TaggerBinning) {
	  gHist1->Fill(tree.GetTheta(particle_index), time, GetTagger()->GetTaggedChannel(j));
	  gHist1CM->Fill(theta_cm, time, GetTagger()->GetTaggedChannel(j));
	}
	else {
	  gHist1->Fill(tree.GetTheta(particle_index), time);
	  gHist1CM->Fill(theta_cm, time);
	}
      }

      // CB trigger 
      if (tp[0] == 0) {
	if (nerror==0 && helicity==kFALSE) {
	  if (TaggerBinning) {
	    gCBHist0->Fill(tree.GetTheta(particle_index), time, GetTagger()->GetTaggedChannel(j));
	    gCBHist0CM->Fill(theta_cm, time, GetTagger()->GetTaggedChannel(j));
	  }
	  else {
	    gCBHist0->Fill(tree.GetTheta(particle_index), time);
	    gCBHist0CM->Fill(theta_cm, time);
	  }
	}
	else if (nerror==0 && helicity==kTRUE) {
	  if (TaggerBinning) {
	    gCBHist1->Fill(tree.GetTheta(particle_index), time, GetTagger()->GetTaggedChannel(j));
	    gCBHist1CM->Fill(theta_cm, time, GetTagger()->GetTaggedChannel(j));
	  }
	  else {
	    gCBHist1->Fill(tree.GetTheta(particle_index), time);
	    gCBHist1CM->Fill(theta_cm, time);
	  }
	}
      }

      // TAPS trigger
      else if (tp[0] == 19) {
	if (nerror==0 && helicity==kFALSE) {
	  if (TaggerBinning) {
	    gTAPSHist0->Fill(tree.GetTheta(particle_index), time, GetTagger()->GetTaggedChannel(j));
	    gTAPSHist0CM->Fill(theta_cm, time, GetTagger()->GetTaggedChannel(j));
	  }
	  else {
	    gTAPSHist0->Fill(tree.GetTheta(particle_index), time);
	    gTAPSHist0CM->Fill(theta_cm, time);
	  }
	}
	else if (nerror==0 && helicity==kTRUE) {
	  if (TaggerBinning) {
	    gTAPSHist1->Fill(tree.GetTheta(particle_index), time, GetTagger()->GetTaggedChannel(j));
	    gTAPSHist1CM->Fill(theta_cm, time, GetTagger()->GetTaggedChannel(j));
	  }
	  else {
	    gTAPSHist1->Fill(tree.GetTheta(particle_index), time);
	    gTAPSHist1CM->Fill(theta_cm, time);
	  }
	}
      }
    }
}


void Pi0X::FillFPD(const GTreeParticle& tree, GH1* gHist)
{

  for (Int_t i = 0; i < tree.GetNParticles(); i++)
    {
      for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
	  // Is tagger channel rejected by user?
	  if(GetTagger()->GetTaggedChannel(j) < GetTC_cut_min()) continue;
	  if(GetTagger()->GetTaggedChannel(j) > GetTC_cut_max()) continue;
	  
	  // calc particle time diff
	  Double_t time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);

	  // Fill FPD histograms 
	  gHist->Fill(GetTagger()->GetTaggedChannel(j), time);
	}
    }
}

void Pi0X::FillFPD(const GTreeParticle& tree, GH1* gHist0, GH1* gHist1)
{
  Int_t nerror = GetTrigger()->GetNErrors();
  Bool_t helicity = GetTrigger()->GetHelicity();

  for (Int_t i = 0; i < tree.GetNParticles(); i++)
    {
      for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
	  // Is tagger channel rejected by user?
	  if(GetTagger()->GetTaggedChannel(j) < GetTC_cut_min()) continue;
	  if(GetTagger()->GetTaggedChannel(j) > GetTC_cut_max()) continue;
	  
	  // calc particle time diff
	  Double_t time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);

	  // Fill FPD histograms 
	  if (nerror==0 && helicity==kFALSE)
	    gHist0->Fill(GetTagger()->GetTaggedChannel(j), time);
	  else if (nerror==0 && helicity==kTRUE)
	    gHist1->Fill(GetTagger()->GetTaggedChannel(j), time);
	}
    }
}


void Pi0X::FillFPD(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
  for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
    {
      // Is tagger channel rejected by user?
      if(GetTagger()->GetTaggedChannel(j) < GetTC_cut_min()) continue;
      if(GetTagger()->GetTaggedChannel(j) > GetTC_cut_max()) continue;
      
      // calc particle time diff
      Double_t time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
      
      // Fill FPD histograms 
      gHist->Fill(GetTagger()->GetTaggedChannel(j), time);
    }
}


void Pi0X::FillFPD(const GTreeParticle& tree, Int_t particle_index, GH1* gHist0, GH1* gCBHist0, GH1* gTAPSHist0, GH1* gHist1, GH1* gCBHist1, GH1* gTAPSHist1)
{
  Int_t nerror = GetTrigger()->GetNErrors();
  Bool_t helicity = GetTrigger()->GetHelicity();
  const Int_t *tp = GetTrigger()->GetTriggerPattern();

  for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
    {
      // Is tagger channel rejected by user?
      if(GetTagger()->GetTaggedChannel(j) < GetTC_cut_min()) continue;
      if(GetTagger()->GetTaggedChannel(j) > GetTC_cut_max()) continue;

      // calc particle time diff
      Double_t time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
      
      // Fill FPD histograms 
      if (nerror==0 && helicity==kFALSE) {
	gHist0->Fill(GetTagger()->GetTaggedChannel(j), time);
	if (tp[0] == 0)
	  gCBHist0->Fill(GetTagger()->GetTaggedChannel(j), time);
	else if (tp[0] == 19)
	  gTAPSHist0->Fill(GetTagger()->GetTaggedChannel(j), time);
      }
      else if (nerror==0 && helicity==kTRUE) {
      	gHist1->Fill(GetTagger()->GetTaggedChannel(j), time);   
	if (tp[0] == 0)
	  gCBHist1->Fill(GetTagger()->GetTaggedChannel(j), time);
	else if (tp[0] == 19)
	  gTAPSHist1->Fill(GetTagger()->GetTaggedChannel(j), time);
      }
  
    }

}


Double_t Pi0X::Compute_ThetaCM(Double_t Egamma, Double_t mom, Double_t ene, Double_t theta, Double_t phi)
{
  //  static const Double_t mproton  = 938.272013; // GeV
  static const Double_t mneutron = 939.565346; // GeV
  Double_t theta_cm;

  Double_t px = mom * TMath::Sin(theta) * TMath::Cos(phi);
  Double_t py = mom * TMath::Sin(theta) * TMath::Sin(phi);
  Double_t pz = mom * TMath::Cos(theta);

  Double_t v;
  v = Egamma / (Egamma + mneutron);

  Double_t gamma = 1./TMath::Sqrt(1 - v*v);
  Double_t px_cm = px;
  Double_t py_cm = py;
  Double_t pz_cm = -gamma*v*ene + gamma*pz;
  Double_t mom_cm = TMath::Sqrt(px_cm*px_cm + py_cm*py_cm + pz_cm*pz_cm);

  theta_cm = TMath::ACos((gamma*mom*TMath::Cos(theta) - gamma*v*ene) / mom_cm);
  
  return theta_cm * TMath::RadToDeg();

}
