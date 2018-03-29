#include "PTaggCal.h"

PTaggCal::PTaggCal()
{
}

PTaggCal::~PTaggCal()
{
}

Bool_t	PTaggCal::Init()
{
    cout << "Initialising tagger calibration analysis..." << endl;
	cout << "--------------------------------------------------" << endl << endl;

    if(!PPhysics::Init()) return kFALSE;

    TaggerCurScal = GetScalerHist("TaggerCurScal");
    if(!TaggerCurScal)
    {
        cout << "No tagger scaler histogram available" << endl;
        return kFALSE;
    }
    TaggerSumScal = (TH1*)TaggerCurScal->Clone("TaggerSumScal");

    for(Int_t i=0; i<=329; i++)
    {
        TGraph temp(0);
        temp.SetName(Form("TaggEnerCalib_%d",i));
        if(i%10) temp.SetMarkerColor(i%10);
	else temp.SetMarkerColor(11);
        TaggEnerCalib.push_back(temp);
    }
    
    cout << "--------------------------------------------------" << endl;
	return kTRUE;
}

Bool_t	PTaggCal::Start()
{
    if(!IsAcquFile())
    {
        cout << "ERROR: Input File is not an Acqu file." << endl;
        return kFALSE;
    }
    SetAsGoATFile();

    TraverseValidEvents();

    return kTRUE;
}

void	PTaggCal::ProcessEvent()
{
}

void	PTaggCal::ProcessScalerRead()
{
    PPhysics::ProcessScalerRead();

    GoosyNewFPD(TaggerCurScal);

    Float_t nmr = GetScalers()->GetNMR();
    if(nmr < minNMR) minNMR = nmr;
    if(nmr > maxNMR) maxNMR = nmr;

    if(TMath::Abs(nmr-preNMR) < 0.0001)
    {
        TaggerSumScal->Add(TaggerCurScal,1);
        numReads++;
    }
    else
    {
        if(numReads > 12)
        {
            TaggerSumScal->Scale(1.0/numReads);
            Int_t maxbin = TaggerSumScal->GetMaximumBin();
            if(TaggerSumScal->GetBinContent(maxbin) > 100)
            {
                Double_t meanchan = (TaggerSumScal->GetMean()-0.5);
                if(TMath::Abs(meanchan-(maxbin-1)) < 1)
                {
                    TaggEnerCalib.at(0).SetPoint(TaggEnerCalib.at(0).GetN(),preNMR,meanchan);
                    TaggEnerCalib.at(329).SetPoint(TaggEnerCalib.at(329).GetN(),preNMR,TaggerSumScal->Integral(maxbin-1,maxbin+1));
                }
                if(maxbin>1 && maxbin<328)
                {
                    cout << numReads << " reads, NMR = " << preNMR << ", Max chan = " << maxbin-1 << "\t";
                    cout << TaggerSumScal->GetBinContent(maxbin-1) << "\t" << TaggerSumScal->GetBinContent(maxbin) << "\t" << TaggerSumScal->GetBinContent(maxbin+1) << endl;
                }

                if(maxbin>1) TaggEnerCalib.at(maxbin-1).SetPoint(TaggEnerCalib.at(maxbin-1).GetN(),preNMR,TaggerSumScal->GetBinContent(maxbin-1)/TaggerSumScal->Integral());
                TaggEnerCalib.at(maxbin).SetPoint(TaggEnerCalib.at(maxbin).GetN(),preNMR,TaggerSumScal->GetBinContent(maxbin)/TaggerSumScal->Integral());
                if(maxbin<328) TaggEnerCalib.at(maxbin+1).SetPoint(TaggEnerCalib.at(maxbin+1).GetN(),preNMR,TaggerSumScal->GetBinContent(maxbin+1)/TaggerSumScal->Integral());
            }
        }

        preNMR = nmr;
        TaggerSumScal->Reset();
        numReads = 0;
    }
    
    TaggerCurScal->Reset();
}

Bool_t	PTaggCal::Write()
{    

    // Write all GH1's and TObjects defined in this class
    if(!(GTreeManager::Write())) return false;

    TCanvas *c1 = new TCanvas("c1","c1");
    TaggEnerCalib.at(0).Draw("A*L");
    c1->SaveAs("TaggEnerCalibMax.pdf");
    
    TCanvas *c2 = new TCanvas("c2","c2");
    TGraph temp(2);
    temp.SetPoint(0,minNMR,0);
    temp.SetPoint(1,maxNMR,1);
    temp.Draw("AP");
    for(Int_t i=1; i<=328; i++) if(TaggEnerCalib.at(i).GetN()) TaggEnerCalib.at(i).Draw("*L");
    c2->SaveAs("TaggEnerCalibPer.pdf");

    TCanvas *c3 = new TCanvas("c3","c3");
    TaggEnerCalib.at(329).Draw("A*L");
    c3->SaveAs("TaggEnerCalibRate.pdf");
    
    TFile f1("TaggEnerCalib.root","RECREATE");
    for(Int_t i=0; i<=329; i++) if(TaggEnerCalib.at(i).GetN()) TaggEnerCalib.at(i).Write();
    c1->Write();
    c2->Write();
    f1.Close();
    
    delete c1;
    delete c2;
    delete c3;
    
    return true;
}
