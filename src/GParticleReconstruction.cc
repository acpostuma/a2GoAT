#include "GParticleReconstruction.h"


using namespace std;


GParticleReconstruction::GParticleReconstruction()  :
    CBTime(0),
    TAPSTime(0),
    CBTimeAfterCut(0),
    TAPSTimeAfterCut(0),
    DoScalerCorrection(kFALSE),
    DoTrigger(kFALSE),
    E_Sum(50),
    multiplicity(1)
{
    CBTimeCut[0]    = -1000000.0;
    CBTimeCut[1]    = 1000000.0;
    TAPSTimeCut[0]    = -1000000.0;
    TAPSTimeCut[1]    = 1000000.0;
}

GParticleReconstruction::~GParticleReconstruction()
{
    if(CBTime)              delete  CBTime;
    if(TAPSTime)            delete  TAPSTime;
    if(CBTimeAfterCut)      delete  CBTimeAfterCut;
    if(TAPSTimeAfterCut)    delete  TAPSTimeAfterCut;
}

Bool_t  GParticleReconstruction::Trigger()
{
    if(trigger->GetESum() < E_Sum)
        return kFALSE;
    if(trigger->GetMult() < multiplicity)
        return kFALSE;
    return kTRUE;
}

void  GParticleReconstruction::ProcessEvent()
{
    if(DoTrigger)
    {
        if(!Trigger())
            return;
    }
    if(rawEvent->GetNCB() != 2 && rawEvent->GetNCB() != 6 && rawEvent->GetNCB() != 10)
        return;
    if(rawEvent->GetNTAPS()>1)
        return;

    for(Int_t i=0; i<rawEvent->GetNParticles(); i++)
    {
        if(rawEvent->GetApparatus(i) == GTreeRawEvent::APPARATUS_CB)
            CBTime->Fill(rawEvent->GetTime(i));
        if(rawEvent->GetApparatus(i) == GTreeRawEvent::APPARATUS_TAPS)
            TAPSTime->Fill(rawEvent->GetTime(i));
    }

    for(Int_t i=0; i<rawEvent->GetNParticles(); i++)
    {
        if(rawEvent->GetApparatus(i) == GTreeRawEvent::APPARATUS_CB)
        {
            if(rawEvent->GetTime(i)<CBTimeCut[0] || rawEvent->GetTime(i)>CBTimeCut[1])
                return;
        }
        if(rawEvent->GetApparatus(i) == GTreeRawEvent::APPARATUS_TAPS)
        {
            if(rawEvent->GetTime(i)<TAPSTimeCut[0] || rawEvent->GetTime(i)>TAPSTimeCut[1])
                return;
        }
    }

    GCorrectScalers::ProcessEvent();

    photons->Clear();
    protons->Clear();
    for(Int_t i=0; i<rawEvent->GetNParticles(); i++)
    {
        if(rawEvent->GetApparatus(i) == GTreeRawEvent::APPARATUS_CB)
        {
            CBTimeAfterCut->Fill(rawEvent->GetTime(i));
            photons->AddParticle(rawEvent->GetVector(i), i);
        }
        if(rawEvent->GetApparatus(i) == GTreeRawEvent::APPARATUS_TAPS)
        {
            TAPSTimeAfterCut->Fill(rawEvent->GetTime(i));
            protons->AddParticle(rawEvent->GetVector(i, 938.272046), i);
        }
    }
    photons->Fill();
    protons->Fill();
}

Bool_t  GParticleReconstruction::Process()
{
    file_out->cd();
    CBTime              = new TH1D("CBTimeOR", "CBTimeOR", 10000, -1000, 1000);
    TAPSTime            = new TH1D("TAPSTimeOR", "TAPSTimeOR", 10000, -1000, 1000);
    CBTimeAfterCut      = new TH1D("CBTimeOR_Cut", "CBTimeOR_Cut", 10000, -1000, 1000);
    TAPSTimeAfterCut    = new TH1D("TAPSTimeOR_Cut", "TAPSTimeOR_Cut", 10000, -1000, 1000);

    if(DoScalerCorrection)
    {
        if(!GCorrectScalers::Process())  return kFALSE;

        if(!Write(CBTime))  return kFALSE;
        if(!Write(TAPSTime))  return kFALSE;
        if(!Write(CBTimeAfterCut))  return kFALSE;
        if(!Write(TAPSTimeAfterCut))  return kFALSE;
        return kTRUE;
    }

    scalers->Clone();

    file_out->cd();
    taggerTime  = new TH1D("TaggerTimeOR", "TaggerTimeOR", 10000, -1000, 1000);
    accepted    = new TH1I("Accepted", "Events with correct scalers (all=0,accepted=1,rejected=2)", 3, 0, 3);
    accepted->SetBinContent(1, rawEvent->GetNEntries());
    accepted->SetBinContent(2, rawEvent->GetNEntries());
    accepted->SetBinContent(3, 0);

    TraverseEntries(0, rawEvent->GetNEntries()+1);

    if(!Write())    return kFALSE;
    if(!Write(taggerTime))  return kFALSE;
    if(!Write(accepted))  return kFALSE;
    if(!Write(CBTime))  return kFALSE;
    if(!Write(TAPSTime))  return kFALSE;
    if(!Write(CBTimeAfterCut))  return kFALSE;
    if(!Write(TAPSTimeAfterCut))  return kFALSE;
    return kTRUE;
}
