#include "Pi0XMC.h"

Pi0XMC::Pi0XMC()
{ 
  time 	               = new GH1("time",                 "time",       1400, -700, 700);
  time_all             = new GH1("time_all",             "time_all",   1400, -700, 700);
  time_cut             = new GH1("time_cut",             "time_cut",   1400, -700, 700);

  FPD                  = new GH1("FPD",                  "GoAT - FPD hits (p-r)",                         352, 0, 352);
  FPD_all              = new GH1("FPD_all",              "GoAT - FPD hits all tracks (p-r)",              352, 0, 352);

  theta                = new GH1("theta",                "GOAT - #theta distribution",                    180, 0, 180);
  theta_all            = new GH1("theta_all",            "GOAT - #theta distribution all tracks",         180, 0, 180);
  thetaCM              = new GH1("thetaCM",              "GOAT - #theta CM distribution",                 180, 0, 180);
  thetaCM_all          = new GH1("thetaCM_all",          "GOAT - #theta CM distribution all tracks",      180, 0, 180);
  theta_MC             = new GH1("theta_MC",             "GOAT - #theta MC distribution",                 180, 0, 180);
  thetaCM_MC           = new GH1("thetaCM_MC",           "GOAT - #theta MC CM distribution",              180, 0, 180);

  TaggerAccScal   = new TH1D("TaggerAccScal",   "TaggerAccScal", 352, 0, 352);

}

Pi0XMC::~Pi0XMC()
{
}

Bool_t Pi0XMC::Init()
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

Bool_t Pi0XMC::Start()
{
  cout << "=== Pi0XMC::Start() " << endl;
  if(!IsGoATFile()){
    cout << "ERROR: Input File is not a GoAT file." << endl;
    return kFALSE;
  }
  SetAsPhysicsFile();

  TraverseValidEvents();

  return kTRUE;
}

void Pi0XMC::ProcessEvent()
{

  FillThetaMC(theta_MC, thetaCM_MC, kTRUE);

  // Time diff (tagger - pi0)
  FillTime(*GetNeutralPions(),0,time);
  FillTime(*GetNeutralPions(), time_all);
  FillTimeCut(*GetNeutralPions(),0,time_cut);

  // FPD hits
  FillFPD(*GetNeutralPions(), 0, FPD);
  FillFPD(*GetNeutralPions(), FPD_all);

  // Angular distributions
  FillTheta(*GetNeutralPions(), 0, theta, thetaCM, kTRUE);
  FillTheta(*GetNeutralPions(), theta_all, thetaCM_all, kTRUE);

  evtNum++;
}  

void	Pi0XMC::ProcessScalerRead()
{
  // Fill Tagger Scalers
  FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal);
}

Bool_t	Pi0XMC::Write()
{
  cout << "Pi0XMC::Write() " << endl;
  // Write all GH1's and TObjects defined in this class
  return GTreeManager::Write();
}



void Pi0XMC::FillTheta(const GTreeParticle& tree, GH1* gHist, GH1* gHistCM, Bool_t TaggerBinning)
{
  //  for (Int_t i = 0; i < tree.GetNParticles(); i++)
  for (Int_t i = 0; i < 1; i++)
    {
      if (GetTagger()->GetNTagged()==0) continue;
      //      for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
      for (Int_t j = 0; j < 1; j++)
	{
	  // Is tagger channel rejected by user?
	  if(GetTagger()->GetTaggedChannel(j) < GetTC_cut_min()) continue;
	  if(GetTagger()->GetTaggedChannel(j) > GetTC_cut_max()) continue;
	  
	  // calc particle time diff
	  Double_t time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);
	  Double_t mom = TMath::Sqrt(tree.GetTotalEnergy(i)*tree.GetTotalEnergy(i) - mpi0*mpi0);
	  Double_t theta_cm = Compute_ThetaCM(GetTagger()->GetTaggedEnergy(j), mom, tree.GetTotalEnergy(i), tree.GetThetaRad(i), TVector2::Phi_0_2pi(tree.GetPhiRad(i))); // deg
	  // Fill theta histograms 
	  if (TaggerBinning) {
	    gHist->Fill(tree.GetTheta(i), time,  GetTagger()->GetTaggedChannel(j));
	    gHistCM->Fill(theta_cm, time, GetTagger()->GetTaggedChannel(j));
	  }
	  else {
	    gHist->Fill(tree.GetTheta(i), time);
	    gHistCM->Fill(theta_cm, time);
	  }
	}
    }
}



void Pi0XMC::FillTheta(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, GH1* gHistCM, Bool_t TaggerBinning)
{
  //  for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
  for (Int_t j = 0; j < 1; j++)
    {
      // Is tagger channel rejected by user?
      if(GetTagger()->GetTaggedChannel(j) < GetTC_cut_min()) continue;
      if(GetTagger()->GetTaggedChannel(j) > GetTC_cut_max()) continue;
      
      // calc particle time diff
      Double_t time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
      Double_t mom = TMath::Sqrt(tree.GetTotalEnergy(particle_index)*tree.GetTotalEnergy(particle_index) - mpi0*mpi0);
      Double_t theta_cm = Compute_ThetaCM(GetTagger()->GetTaggedEnergy(j), mom, tree.GetTotalEnergy(particle_index), tree.GetThetaRad(particle_index), TVector2::Phi_0_2pi(tree.GetPhiRad(particle_index))); // deg

      // Fill theta histograms 
      if (TaggerBinning) {
	gHist->Fill(tree.GetTheta(particle_index), time, GetTagger()->GetTaggedChannel(j));
	gHistCM->Fill(theta_cm, time, GetTagger()->GetTaggedChannel(j));
      }
      else {
	gHist->Fill(tree.GetTheta(particle_index), time);
	gHistCM->Fill(theta_cm, time);
      }
    }
}



void Pi0XMC::FillThetaMC(GH1* gHist, GH1* gHistCM, Bool_t TaggerBinning)
{
  //  for (Int_t i = 0; i < geant->GetNTrueParticles(); i++) --> i = 0 for pi0
    {
      //      for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
      for (Int_t j = 0; j < 1; j++)
	{
	  // Is tagger channel rejected by user?
 	  if(GetTagger()->GetTaggedChannel(j) < GetTC_cut_min()) continue;
 	  if(GetTagger()->GetTaggedChannel(j) > GetTC_cut_max()) continue;
	  
	  Double_t time = 0.;
	  Double_t Ebeam = geant->GetBeam().E()*1e3;
	  Double_t Energy = (geant->GetTrueVector(0)).Energy()*1e3; // MeV
	  Double_t Momentum = (geant->GetTrueVector(0)).P()*1e3; // MeV
	  TVector3 p(geant->GetTrueVector(0).Px()*1e3,geant->GetTrueVector(0).Py()*1e3, geant->GetTrueVector(0).Pz()*1e3);
	  //	  cout << "MOMENTUM " << p.X() << "\t " << p.Y() << "\t" << p.Z() << "\t" << p.Mag() << "\t" << p.Theta()*TMath::RadToDeg() << "\t" << TVector2::Phi_0_2pi(p.Phi())*TMath::RadToDeg() << endl;
	  Double_t ThetaRad = (geant->GetTrueVector(0)).Theta(); // Rad
	  Double_t PhiRad = TVector2::Phi_0_2pi((geant->GetTrueVector(0)).Phi()); // Rad
	  Double_t thetaMC_CM = Compute_ThetaCM(Ebeam, Momentum, Energy, ThetaRad, PhiRad); // deg

	  // Fill theta histograms 
	  if (TaggerBinning) { 
	    gHist->Fill(ThetaRad*TMath::RadToDeg(), time, GetTagger()->GetTaggedChannel(j));
	    gHistCM->Fill(thetaMC_CM, time, GetTagger()->GetTaggedChannel(j));
	  }
	  else {
	    gHist->Fill(ThetaRad*TMath::RadToDeg(), time);
	    gHistCM->Fill(thetaMC_CM*TMath::RadToDeg(), time);
	  }
	}
    }
}


void Pi0XMC::FillFPD(const GTreeParticle& tree, GH1* gHist)
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



void Pi0XMC::FillFPD(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
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


Double_t Pi0XMC::Compute_ThetaCM(Double_t Egamma, Double_t mom, Double_t ene, Double_t theta, Double_t phi)
{
  static const Double_t mproton  = 938.272013; // GeV
  Double_t theta_cm;

  Double_t px = mom * TMath::Sin(theta) * TMath::Cos(phi);
  Double_t py = mom * TMath::Sin(theta) * TMath::Sin(phi);
  Double_t pz = mom * TMath::Cos(theta);

  Double_t v;
  v = Egamma / (Egamma + mproton);

  Double_t gamma = 1./TMath::Sqrt(1 - v*v);
  Double_t px_cm = px;
  Double_t py_cm = py;
  Double_t pz_cm = -gamma*v*ene + gamma*pz;
  Double_t mom_cm = TMath::Sqrt(px_cm*px_cm + py_cm*py_cm + pz_cm*pz_cm);

  theta_cm = TMath::ACos((gamma*mom*TMath::Cos(theta) - gamma*v*ene) / mom_cm);
  
  return theta_cm * TMath::RadToDeg();

}
