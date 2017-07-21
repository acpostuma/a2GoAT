#include "PCalibESum.h"

PCalibESum::PCalibESum()
{ 
    time 	= new GH1("time", 	"time", 	1400, -700, 700);
    time_cut 	= new GH1("time_cut", 	"time_cut", 	1400, -700, 700);

    time_2g 	= new GH1("time_2g",	"time_2g", 	1400, -700, 700);
    time_2g_cut = new GH1("time_2g_cut","time_2g_cut", 	1400, -700, 700);

    IM 		= new GH1("IM", 	"IM", 		400,   0, 400);
    IM_2g 	= new GH1("IM_2g", 	"IM_2g", 	400,   0, 400);
  
    MM		= new GH1("MM", 	"MM", 	 	400,   800, 1200);     
    MM_2g	= new GH1("MM_2g", 	"MM_2g", 	400,   800, 1200);

    ESumCalib = new GH2("ESumCalib", "ESumCalib", 410, 0, 4100, 100, 0, 1000);
}

PCalibESum::~PCalibESum()
{
}

Bool_t	PCalibESum::Init()
{
	cout << "Initialising physics analysis..." << endl;
	cout << "--------------------------------------------------" << endl << endl;

	if(!InitBackgroundCuts()) return kFALSE;
    if(!InitTargetMass()) return kFALSE;

    if(!PPhysics::Init()) return kFALSE;

    cout << "--------------------------------------------------" << endl;
	return kTRUE;
}

Bool_t	PCalibESum::Start()
{
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();

    GetDetectorHits()->SetNewBranchAddress("CBESum",&CBESum);

    TraverseValidEvents();

    return kTRUE;
}

void	PCalibESum::ProcessEvent()
{
    // Uncomment the following line to decode double hits in the tagger
    //GetTagger()->DecodeDoubles();

    // fill time diff (tagger - pi0), all pi0
    FillTime(*GetNeutralPions(),time);
    FillTimeCut(*GetNeutralPions(),time_cut);
	
	// fill missing mass, all pi0
    FillMissingMass(*GetNeutralPions(),MM);
	
	// fill invariant mass, all pi0
    FillMass(*GetNeutralPions(),IM);
		
    // Some neutral decays
    for (Int_t i = 0; i < GetNeutralPions()->GetNParticles(); i++)
    {
        // Fill MM for 2 photon decay
        if (GetNeutralPions()->GetNSubParticles(i) == 2)
        {
            // fill time diff (tagger - pi0), this pi0
            FillTime(*GetNeutralPions(),i,time_2g);
            FillTimeCut(*GetNeutralPions(),i,time_2g_cut);
			
            // fill missing mass, this pi0
            FillMissingMass(*GetNeutralPions(),i,MM_2g);
            
            // fill invariant mass, this pi0
            FillMass(*GetNeutralPions(),i,IM_2g);

            // fill esum vs pi0 energy
            FillESum(*GetNeutralPions(),i,ESumCalib);
        }

    }

}

void PCalibESum::FillESum(const GTreeParticle& tree, GH2* gHist, Bool_t TaggerBinning)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
    {
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
        {
            FillESum(tree, i, j, gHist, TaggerBinning);
        }
    }
}

void PCalibESum::FillESum(const GTreeParticle& tree, Int_t particle_index, GH2* gHist, Bool_t TaggerBinning)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
    {
        FillESum(tree, particle_index, i, gHist, TaggerBinning);
    }
}

void PCalibESum::FillESum(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH2* gHist, Bool_t TaggerBinning)
{
    if(RejectTagged(tagger_index)) return;

    // calc particle time diff
    Double_t time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);

    // calc missing p4
    TLorentzVector missingp4 = CalcMissingP4(tree, particle_index,tagger_index);

    // Cut on missing mass
    if(TMath::Abs(missingp4.M()-938.27)<30)
    {
        // Fill GH2
        if(TaggerBinning)   gHist->Fill(CBESum,tree.GetTotalEnergy(particle_index),time, GetTagger()->GetTaggedChannel(tagger_index));
        else gHist->Fill(CBESum,tree.GetTotalEnergy(particle_index),time);
    }
}

void	PCalibESum::ProcessScalerRead()
{
    PPhysics::ProcessScalerRead();
}

Bool_t	PCalibESum::Write()
{
    // Write all GH1's and TObjects defined in this class
    return GTreeManager::Write();
}
