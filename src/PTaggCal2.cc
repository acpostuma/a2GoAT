#include "PTaggCal2.h"

PTaggCal2::PTaggCal2()
{
}

PTaggCal2::~PTaggCal2()
{
}

Bool_t	PTaggCal2::Init()
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

    for(Int_t i=0; i<numReads; i++)
    {
        //TH1* temp = (TH1*)TaggerCurScal->Clone(Form("TaggerScalers_%i",i));
        TaggerScalers.push_back((TH1*)TaggerCurScal->Clone(Form("TaggerScalers_%i",i)));
        NMRReads.push_back(0);
    }

    for(Int_t i=0; i<=369; i++)
    {
        TGraphErrors temp(0);
        temp.SetName(Form("TaggEnerCalib_%d",i));
        if(i%10) temp.SetMarkerColor(i%10);
        else temp.SetMarkerColor(11);
        TaggEnerCalib.push_back(temp);
    }
    
    nReadsStableNMR.SetName("nReadsStableNMR");

    cout << "--------------------------------------------------" << endl;
	return kTRUE;
}

Bool_t	PTaggCal2::Start()
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

void	PTaggCal2::ProcessEvent()
{
}

Bool_t 	PTaggCal2::InitScanName()
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

Bool_t 	PTaggCal2::InitNMRScan()
{
    Int_t n;
    Float_t d;
    string config = ReadConfig("NMR-Scan");
    if(strcmp(config.c_str(), "nokey") == 0)
    {
        cout << "Using default NMR scan settings: " << numReads << " reads with difference of " << NMRdif << endl << endl;
    }
    else if(sscanf( config.c_str(), "%d%f\n", &n, &d) == 2)
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

Bool_t 	PTaggCal2::InitNMRLimits()
{
    Float_t min, max;
    string config = ReadConfig("NMR-Limits");
    if(strcmp(config.c_str(), "nokey") == 0)
    {
        cout << "Using default NMR limits: " << NMRmin << " to " << NMRmax << endl << endl;
    }
    else if(sscanf( config.c_str(), "%f%f\n", &min, &max) == 2)
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

void	PTaggCal2::ProcessScalerRead()
{
    PPhysics::ProcessScalerRead();

    Float_t cur = GetScalers()->GetCur();

    Float_t nmr = GetScalers()->GetNMR();

    if(nmr < NMRmin || nmr > NMRmax) return;

    if(nmr < scanMinNMR) scanMinNMR = nmr;
    if(nmr > scanMaxNMR) scanMaxNMR = nmr;

    Float_t avgNMR = 0;
    for(Int_t i=0; i<numReads; i++) avgNMR += NMRReads.at(i);
    avgNMR = avgNMR/numReads;

    cout << cur << "\t" << nmr << "\t" << TaggerCurScal->GetBinContent(TaggerCurScal->GetMaximumBin()) << "\t" << TaggerCurScal->Integral() << endl;

    Float_t sigNMR = 0;
    if(numReads > 1)
    {
        for(Int_t i=0; i<numReads; i++) sigNMR += TMath::Power((avgNMR-NMRReads.at(i)),2);
        sigNMR = TMath::Sqrt(sigNMR/(numReads-1));
    }
    //readOff = readNum%numReads;

    //if((TMath::Abs(nmr-avgNMR) < NMRdif) && (cur == preCur))
    if((cur != 0) && (cur == preCur))
    {
        for(Int_t i=0; i<(TaggerCurScal->GetNbinsX()); i++) TaggerScalers.at(readNum%numReads)->SetBinContent(i+1,TaggerCurScal->GetBinContent(i+1));
        NMRReads.at(readNum%numReads) = nmr;
        readNum++;
    }
    else if(readNum >= numReads)
    {
        if(sigNMR < NMRdif)
        {
            if((cur == preCur) && (preCur != 0)) cout << "Value found before change for dipole current of " << preCur << endl;

            TaggerSumScal->Reset();
            TaggerSumScal->SetBins(TaggerCurScal->GetNbinsX(),0,TaggerCurScal->GetNbinsX());

            for(Int_t i=0; i<numReads; i++)
            {
                TaggerSumScal->Add(TaggerScalers.at(i),1);
            }

            //GoosyNewFPD(TaggerSumScal);
            //GoosyNewFPDRecabled(TaggerSumScal);
            GoosyNewFPDRecabledYoke(TaggerSumScal);
            //TaggerSumScal->Scale(1.0/numReads);
            newMaxBin = TaggerSumScal->GetMaximumBin();
            cout << "Max channel = " << newMaxBin-1 << endl;
            if(TaggerSumScal->GetBinContent(newMaxBin) > 100)
            {
                Double_t meanchan = (TaggerSumScal->GetMean()-0.5);
                if(TMath::Abs(meanchan-(newMaxBin-1)) < 1)
                {
                    TaggEnerCalib.at(0).SetPoint(TaggEnerCalib.at(0).GetN(),avgNMR,meanchan);
                    TaggEnerCalib.at(0).SetPointError(TaggEnerCalib.at(0).GetN()-1,sigNMR,TaggerSumScal->GetRMS());
                    TaggEnerCalib.at(369).SetPoint(TaggEnerCalib.at(369).GetN(),avgNMR,TaggerSumScal->Integral(newMaxBin-1,newMaxBin+1));
                    TaggEnerCalib.at(369).SetPointError(TaggEnerCalib.at(369).GetN()-1,sigNMR,TMath::Sqrt(TaggerSumScal->Integral(newMaxBin-1,newMaxBin+1)));
                }

                if(newMaxBin>1 && newMaxBin<368)
                {
                    cout << readNum << " reads, NMR = " << avgNMR << ", Max chan = " << newMaxBin-1 << "\t";
                    cout << TaggerSumScal->GetBinContent(newMaxBin-1) << "\t" << TaggerSumScal->GetBinContent(newMaxBin) << "\t" << TaggerSumScal->GetBinContent(newMaxBin+1) << "\t";
                    cout << TaggerSumScal->Integral() << endl;
                }

                Double_t value, error;
                if(newMaxBin>1)
                {
                    value = TaggerSumScal->GetBinContent(newMaxBin-1)/TaggerSumScal->Integral(newMaxBin-2,newMaxBin);
                    TaggEnerCalib.at(newMaxBin-1).SetPoint(TaggEnerCalib.at(newMaxBin-1).GetN(),avgNMR,value);
                    error = value*TMath::Sqrt((1/TaggerSumScal->GetBinContent(newMaxBin-1))-(1/TaggerSumScal->Integral(newMaxBin-2,newMaxBin)));
                    TaggEnerCalib.at(newMaxBin-1).SetPointError(TaggEnerCalib.at(newMaxBin-1).GetN()-1,sigNMR,error);
                }
                value = TaggerSumScal->GetBinContent(newMaxBin)/TaggerSumScal->Integral(newMaxBin-1,newMaxBin+1);
                TaggEnerCalib.at(newMaxBin).SetPoint(TaggEnerCalib.at(newMaxBin).GetN(),avgNMR,value);
                error = value*TMath::Sqrt((1/TaggerSumScal->GetBinContent(newMaxBin))-(1/TaggerSumScal->Integral(newMaxBin-1,newMaxBin+1)));
                TaggEnerCalib.at(newMaxBin).SetPointError(TaggEnerCalib.at(newMaxBin).GetN()-1,sigNMR,error);
                if(newMaxBin<368)
                {
                    value = TaggerSumScal->GetBinContent(newMaxBin+1)/TaggerSumScal->Integral(newMaxBin,newMaxBin+2);
                    TaggEnerCalib.at(newMaxBin+1).SetPoint(TaggEnerCalib.at(newMaxBin+1).GetN(),avgNMR,value);
                    error = value*TMath::Sqrt((1/TaggerSumScal->GetBinContent(newMaxBin+1))-(1/TaggerSumScal->Integral(newMaxBin,newMaxBin+2)));
                    TaggEnerCalib.at(newMaxBin+1).SetPointError(TaggEnerCalib.at(newMaxBin+1).GetN()-1,sigNMR,error);
                }

                newHiPer = TaggerSumScal->GetBinContent(newMaxBin)/TaggerSumScal->Integral();
                newLoPer = TaggerSumScal->GetBinContent(newMaxBin+1)/TaggerSumScal->Integral();
                if(newMaxBin==(preMaxBin-1))
                {
                    Double_t intersect = ((((newHiPer-newLoPer)*preNMR)+((preHiPer-preLoPer)*avgNMR))/((newHiPer-newLoPer)+(preHiPer-preLoPer)));
                    cout << "Intersection found at " << intersect << " between channels " << newMaxBin-1 << " and " << preMaxBin-1 << endl;
                    TLine temp(intersect,0,intersect,1);
                    temp.SetLineColor(2);
                    Intersections.push_back(temp);

                }
                preMaxBin = newMaxBin;
                preNMR = avgNMR;
                preHiPer = TaggerSumScal->GetBinContent(newMaxBin)/TaggerSumScal->Integral();
                preLoPer = TaggerSumScal->GetBinContent(newMaxBin-1)/TaggerSumScal->Integral();
            }
        }
        else if((cur != preCur) && (preCur != 0)) cout << "No value found for dipole current of " << preCur << "(NMR = " << avgNMR << " +/- " << sigNMR << " T)" << endl;
        else cout << "Not sure how I got here!" << endl;

        cout << preCur << endl << "\tAverage NMR for " << numReads << " reads = " << avgNMR << " +/- " << sigNMR << endl;

        avgNMR = NMRReadsAll.at(readNum-1);
        for(Int_t i=(readNum-2); i>=0; i--)
        {
            avgNMR = avgNMR*(readNum-i-1);
            avgNMR += NMRReadsAll.at(i);
            avgNMR = avgNMR/(readNum-i);

            sigNMR = 0;
            for(Int_t j=(readNum-1); j>=i; j--) sigNMR += TMath::Power((avgNMR-NMRReadsAll.at(j)),2);
            sigNMR = TMath::Sqrt(sigNMR/(readNum-i-1));

            if(sigNMR > NMRdif)
            {
                avgNMR = avgNMR*(readNum-i);
                avgNMR -= NMRReadsAll.at(i);
                avgNMR = avgNMR/(readNum-i-1);

                sigNMR = 0;
                for(Int_t j=(readNum-1); j>i; j--) sigNMR += TMath::Power((avgNMR-NMRReadsAll.at(j)),2);
                sigNMR = TMath::Sqrt(sigNMR/(readNum-i-2));

                if((readNum-i-1) < 100) nReadsStableNMR.SetPoint(nReadsStableNMR.GetN(),preCur,readNum-i-1);
                cout << "\tAverage NMR for " << readNum-i-1 << " reads = " << avgNMR << " +/- " << sigNMR << endl;
                break;
            }
        }

        NMRReadsAll.resize(0);
        readNum = 0;
    }
    else
    {
        NMRReadsAll.resize(0);
        readNum = 0;
    }

    NMRReadsAll.push_back(nmr);
    preCur = cur;
    TaggerCurScal->Reset();
}

Bool_t	PTaggCal2::Write()
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
    temp.SetPoint(0,scanMinNMR,0);
    temp.SetPoint(1,scanMaxNMR,1);
    temp.SetTitle(scanName);
    temp.Draw("AP");
    for(Int_t i=1; i<=368; i++) if(TaggEnerCalib.at(i).GetN())
    {
        nChan++;
        l1.AddEntry(&TaggEnerCalib.at(i),Form("%d (%d)",i-1,369-i),"p");
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
    TaggEnerCalib.at(369).Draw("A*L");
    c3->SaveAs("TaggEnerCalibRate.pdf");
    
    TFile f1("TaggEnerCalib.root","RECREATE");
    for(Int_t i=0; i<=369; i++) if(TaggEnerCalib.at(i).GetN()) TaggEnerCalib.at(i).Write();
    nReadsStableNMR.Write();
    c1->Write();
    c2->Write();
    f1.Close();
    
    delete c1;
    delete c2;
    delete c3;
    
    return true;
}
