#include "TotalInclusive.h"

TotalInclusive::TotalInclusive()
{ 
  time 	           = new GH1("time",              "time",       1400, -700, 700);
  time_all         = new GH1("time_all",          "time_all",   1400, -700, 700);
  time_cut         = new GH1("time_cut",          "time_cut",   1400, -700, 700);

  FPD_TI_hel0      = new GH1("FPD_TI_hel0",      "GoAT - FPD hits 4TI all tracks (p-r) hel 0",   352, 0, 352);
  FPD_TI_hel1      = new GH1("FPD_TI_hel1",      "GoAT - FPD hits 4TI all tracks (p-r) hel 1",   352, 0, 352);
  FPD_TI_CB_hel0   = new GH1("FPD_TI_CB_hel0",   "GoAT - FPD hits 4TI CB trigger (p-r) hel 0",   352, 0, 352);
  FPD_TI_CB_hel1   = new GH1("FPD_TI_CB_hel1",   "GoAT - FPD hits 4TI CB trigger (p-r) hel 1",   352, 0, 352);
  FPD_TI_TAPS_hel0 = new GH1("FPD_TI_TAPS_hel0", "GoAT - FPD hits 4TI TAPS trigger (p-r) hel 0", 352, 0, 352);
  FPD_TI_TAPS_hel1 = new GH1("FPD_TI_TAPS_hel1", "GoAT - FPD hits 4TI TAPS trigger (p-r) hel 1", 352, 0, 352);

  CBtrigg          = new GH1("Esum_CBtrigg",     "GoAT - Esum in case of CB trigger",            200, 0, 1500);
  TAPStrigg        = new GH1("Esum_TAPStrigg",   "GoAT - Esum in case of TAPS trigger",          200, 0, 1500);
  MIXtrigg         = new GH1("Esum_MIXtrigg",    "GoAT - Esum in case of MIX (0, 1, 16) trigger",200, 0, 1500);

  TaggerAccScal    = new TH1D("TaggerAccScal",   "TaggerAccScal", 352, 0, 352);
}

TotalInclusive::~TotalInclusive()
{
}

Bool_t TotalInclusive::Init()
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

Bool_t TotalInclusive::Start()
{
  cout << "=== TotalInclusive::Start() " << endl;
  if(!IsGoATFile()){
    cout << "ERROR: Input File is not a GoAT file." << endl;
    return kFALSE;
  }
  SetAsPhysicsFile();

  TraverseValidEvents();

  return kTRUE;
}

void TotalInclusive::ProcessEvent()
{
  //  cout << "New event " << endl;
  // Time diff (tagger - pi0)
  FillTime(*GetNeutralPions(),0,time);
  FillTime(*GetNeutralPions(), time_all);
  FillTimeCut(*GetNeutralPions(),0,time_cut);

  // FPD hits
  FillFPD_TI(*GetTracks(), FPD_TI_hel0, FPD_TI_hel1, 
	     FPD_TI_CB_hel0, FPD_TI_CB_hel1, 
	     FPD_TI_TAPS_hel0, FPD_TI_TAPS_hel1);

  // check Esum
  FillEsum(*GetTracks(), CBtrigg, TAPStrigg, MIXtrigg);


  evtNum++;
}  

void	TotalInclusive::ProcessScalerRead()
{
  // Fill Tagger Scalers
  FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal);
}

Bool_t	TotalInclusive::Write()
{
  cout << "TotalInclusive::Write() " << endl;
  // Write all GH1's and TObjects defined in this class
  return GTreeManager::Write();
}


void TotalInclusive::FillFPD_TI(const GTreeTrack& tree, GH1* gHist0, GH1* gHist1, GH1* gCBHist0, GH1* gCBHist1, GH1* gTAPSHist0, GH1* gTAPSHist1)
{
  Int_t nerror = GetTrigger()->GetNErrors();
  Bool_t helicity = GetTrigger()->GetHelicity();

  Double_t time = -999999;

  const Int_t *tp = GetTrigger()->GetTriggerPattern();

  if (tree.GetNTracks() > 0)
    {
      for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
	  // Is tagger channel rejected by user?
	  if(GetTagger()->GetTaggedChannel(j) < GetTC_cut_min()) continue;
	  if(GetTagger()->GetTaggedChannel(j) > GetTC_cut_max()) continue;
	  
	  // calc particle time diff
	  for (Int_t i = 0; i < tree.GetNTracks(); i++) {
	    time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);
	    if (GHistBGSub::IsPrompt(time)) break;
	  }

	  // Fill FPD histograms 
	  // all triggers
	  if (nerror==0 && helicity==kFALSE)
	    gHist0->Fill(GetTagger()->GetTaggedChannel(j), time);
	  else if (nerror==0 && helicity==kTRUE)
	    gHist1->Fill(GetTagger()->GetTaggedChannel(j), time);
	  // CB trigger 
	  if (tp[0] == 0) {
	    if (nerror==0 && helicity==kFALSE)
	      gCBHist0->Fill(GetTagger()->GetTaggedChannel(j), time);
	    else if (nerror==0 && helicity==kTRUE)
	      gCBHist1->Fill(GetTagger()->GetTaggedChannel(j), time);
	  }
	  // TAPS trigger
	  else if (tp[0] == 19) {
	    if (nerror==0 && helicity==kFALSE)
	      gTAPSHist0->Fill(GetTagger()->GetTaggedChannel(j), time);
	    else if (nerror==0 && helicity==kTRUE)
	      gTAPSHist1->Fill(GetTagger()->GetTaggedChannel(j), time);
	  }
	}
    }
}

void TotalInclusive::FillEsum(const GTreeTrack& tree, GH1* gHist0, GH1* gHist1, GH1* gHist2)
{
  //  cout << "===> New event " << endl;
  Int_t nerror = GetTrigger()->GetNErrors();

  Double_t Esum = 0;

  const Int_t *tp = GetTrigger()->GetTriggerPattern();

  for (Int_t i=0; i<tree.GetNTracks(); i++) {
    Esum += tree.GetClusterEnergy(i);
  }

  if (nerror==0 && tree.GetNTracks()>0) {
    if (tp[0]==0) // CB trigger
      gHist0->Fill(Esum);
    else if (tp[0]==19) //TAPS trigger
      gHist1->Fill(Esum);
  }

  Int_t counter = 0;
  for (Int_t k=0; k<GetTrigger()->GetNTriggerPattern(); k++) {
    if (tp[k]==0 || tp[k]==1 || tp[k]==16)
      counter++;
  }
  if (counter==3) {
    //    cout << "Filling mix histo " << endl;
    gHist2->Fill(Esum);
  }

}

