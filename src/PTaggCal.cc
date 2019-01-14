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

    if(!InitScanName()) return kFALSE;
    if(!InitNMRScan()) return kFALSE;
    if(!InitNMRLimits()) return kFALSE;

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

Bool_t 	PTaggCal::InitScanName()
{
    char name[256];
    string config = ReadConfig("Scan-Name");
    if(strcmp(config.c_str(), "nokey") == 0)
    {
        cout << "Using default scan name: " << scanName << endl << endl;
    }
    else if(sscanf( config.c_str(), "%[^\n]s", name) == 1)
    {
        cout << "Naming NMR scan: " << name << endl << endl;
        scanName = name;
    }
    else
    {
        cout << "Scan name not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

Bool_t 	PTaggCal::InitNMRScan()
{
    Int_t n;
    Double_t d;
    string config = ReadConfig("NMR-Scan");
    if(strcmp(config.c_str(), "nokey") == 0)
    {
        cout << "Using default NMR scan settings: " << numReads << " reads with difference of " << NMRdif << endl << endl;
    }
    else if(sscanf( config.c_str(), "%d%lf\n", &n, &d) == 2)
    {
        cout << "Setting NMR scan: " << n << " reads with difference of " << d << endl << endl;
        numReads = n;
        NMRdif = d;
    }
    else
    {
        cout << "NMR scan not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

Bool_t 	PTaggCal::InitNMRLimits()
{
    Double_t min, max;
    string config = ReadConfig("NMR-Limits");
    if(strcmp(config.c_str(), "nokey") == 0)
    {
        cout << "Using default NMR limits: " << NMRmin << " to " << NMRmax << endl << endl;
    }
    else if(sscanf( config.c_str(), "%lf%lf\n", &min, &max) == 2)
    {
        cout << "Setting NMR limits: " << min << " to " << max << endl << endl;
        NMRmin = min;
        NMRmax = max;
    }
    else
    {
        cout << "NMR limits not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

void	PTaggCal::ProcessScalerRead()
{
    PPhysics::ProcessScalerRead();

    //nScalerReads++;
    //if(nScalerReads>)

    TaggerSumScal->SetBins(TaggerCurScal->GetNbinsX(),0,TaggerCurScal->GetNbinsX());

    Float_t nmr = GetScalers()->GetNMR();
    if(nmr < NMRmin || nmr > NMRmax) return;

    if(nmr < minNMR) minNMR = nmr;
    if(nmr > maxNMR) maxNMR = nmr;

    if(TMath::Abs(nmr-newNMR) < NMRdif)
    {
        TaggerSumScal->Add(TaggerCurScal,1);
        readNum++;
    }
    else
    {
        if(readNum > numReads)
        {
            //GoosyNewFPD(TaggerSumScal);
            GoosyNewFPDRecabled(TaggerSumScal);
            TaggerSumScal->Scale(1.0/readNum);
            newMaxBin = TaggerSumScal->GetMaximumBin();
            //cout << "Max channel = " << newMaxBin-1 << endl;
            if(TaggerSumScal->GetBinContent(newMaxBin) > 100)
            {
                Double_t meanchan = (TaggerSumScal->GetMean()-0.5);
                if(TMath::Abs(meanchan-(newMaxBin-1)) < 1)
                {
                    TaggEnerCalib.at(0).SetPoint(TaggEnerCalib.at(0).GetN(),newNMR,meanchan);
                    TaggEnerCalib.at(329).SetPoint(TaggEnerCalib.at(329).GetN(),newNMR,TaggerSumScal->Integral(newMaxBin-1,newMaxBin+1));
                }
                /*
                if(newMaxBin>1 && newMaxBin<328)
                {
                    cout << readNum << " reads, NMR = " << newNMR << ", Max chan = " << newMaxBin-1 << "\t";
                    cout << TaggerSumScal->GetBinContent(newMaxBin-1) << "\t" << TaggerSumScal->GetBinContent(newMaxBin) << "\t" << TaggerSumScal->GetBinContent(newMaxBin+1) << endl;
                }
                */
                if(newMaxBin>1) TaggEnerCalib.at(newMaxBin-1).SetPoint(TaggEnerCalib.at(newMaxBin-1).GetN(),newNMR,TaggerSumScal->GetBinContent(newMaxBin-1)/TaggerSumScal->Integral());
                TaggEnerCalib.at(newMaxBin).SetPoint(TaggEnerCalib.at(newMaxBin).GetN(),newNMR,TaggerSumScal->GetBinContent(newMaxBin)/TaggerSumScal->Integral());
                if(newMaxBin<328) TaggEnerCalib.at(newMaxBin+1).SetPoint(TaggEnerCalib.at(newMaxBin+1).GetN(),newNMR,TaggerSumScal->GetBinContent(newMaxBin+1)/TaggerSumScal->Integral());

                newHiPer = TaggerSumScal->GetBinContent(newMaxBin)/TaggerSumScal->Integral();
                newLoPer = TaggerSumScal->GetBinContent(newMaxBin+1)/TaggerSumScal->Integral();
                if(newMaxBin==(preMaxBin-1))
                {
                    Double_t intersect = ((((newHiPer-newLoPer)*preNMR)+((preHiPer-preLoPer)*newNMR))/((newHiPer-newLoPer)+(preHiPer-preLoPer)));
                    cout << "Intersection found at " << intersect << " between channels " << newMaxBin-1 << " and " << preMaxBin-1 << endl;
                    TLine temp(intersect,0,intersect,1);
                    temp.SetLineColor(2);
                    Intersections.push_back(temp);

                }
                preMaxBin = newMaxBin;
                preNMR = newNMR;
                preHiPer = TaggerSumScal->GetBinContent(newMaxBin)/TaggerSumScal->Integral();
                preLoPer = TaggerSumScal->GetBinContent(newMaxBin-1)/TaggerSumScal->Integral();
            }
        }

        newNMR = nmr;
        TaggerSumScal->Reset();
        readNum = 0;
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
    
    TLegend l1(0.85,0.92,0.95,0.95);
    l1.SetHeader("Chan (JA)");
    Int_t nChan = 0;
    scanName.Append(";NMR (T);Hits in Channel / Hits in Ladder");

    TCanvas *c2 = new TCanvas("c2","c2");
    TGraph temp(2);
    temp.SetPoint(0,minNMR,0);
    temp.SetPoint(1,maxNMR,1);
    temp.SetTitle(scanName);
    temp.Draw("AP");
    for(Int_t i=1; i<=328; i++) if(TaggEnerCalib.at(i).GetN())
    {
        nChan++;
        l1.AddEntry(&TaggEnerCalib.at(i),Form("%d (%d)",i-1,329-i),"p");
        TaggEnerCalib.at(i).Draw("*L");
    }
    for(UInt_t i=0; i<Intersections.size(); i++) Intersections.at(i).Draw();
    if(nChan<30)
    {
        l1.SetY1(0.92-0.03*nChan);
        l1.Draw();
    }
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
