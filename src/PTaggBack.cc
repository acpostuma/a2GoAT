#include "PTaggBack.h"

PTaggBack::PTaggBack()
{
    FreeScalers = true;
}

PTaggBack::~PTaggBack()
{
}

Bool_t	PTaggBack::Init()
{
    cout << "Initialising tagging efficiency analysis..." << endl;
	cout << "--------------------------------------------------" << endl << endl;

    if(!InitBackgroundCuts()) return kFALSE;
    if(!InitFreeScalers()) return kFALSE;
    if(!InitBackgroundScalers()) return kFALSE;

    if(!PPhysics::Init()) return kFALSE;

    TaggerAccScal = GetScalerHist("TaggerAccScal");
    if(!TaggerAccScal)
    {
        cout << "No tagger scaler histogram available" << endl;
        return kFALSE;
    }
    TaggerPreScal = (TH1*)TaggerAccScal->Clone("TaggerPreScal");
    TaggerCurScal = (TH1*)TaggerAccScal->Clone("TaggerCurScal");
    TaggerFirScal = (TH1*)TaggerAccScal->Clone("TaggerFirScal");
    LiveTimeScal = GetScalerHist("LiveTimeScal");
    if(!LiveTimeScal)
    {
        cout << "No live time histogram available" << endl;
        return kFALSE;
    }
    
    cout << "--------------------------------------------------" << endl;
	return kTRUE;
}

Bool_t	PTaggBack::Start()
{
    if(!IsAcquFile())
    {
        cout << "ERROR: Input File is not an Acqu file." << endl;
        return kFALSE;
    }
    SetAsGoATFile();

    TraverseValidEvents();

    //GoosyTagger(TaggerAccScal);
    //GoosyVuprom(TaggerAccScal);
    GoosyNewFPD(TaggerAccScal);

    return kTRUE;
}

void	PTaggBack::ProcessEvent()
{
}

Bool_t 	PTaggBack::InitFreeScalers()
{
    Int_t fs;
    string config = ReadConfig("Free-Running-Scal");
    if(strcmp(config.c_str(), "nokey") == 0)
    {
        cout << "Assuming scalers are free running! " << endl << endl;
        FreeScalers = true;
    }
    else if(sscanf( config.c_str(), "%d\n", &fs) == 1)
    {
        cout << "Setting free running scalers: " << fs << endl << endl;
        FreeScalers = (Bool_t)fs;
    }
    else
    {
        cout << "Free running scalers not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

Bool_t 	PTaggBack::InitBackgroundScalers()
{
    Int_t sn, s1, s2;
    string config = ReadConfig("Background-Scalers");
    if(strcmp(config.c_str(), "nokey") == 0)
    {
        cout << "No background scalers selection! " << endl << endl;
    }
    else if(sscanf( config.c_str(), "%d%d%d\n", &sn, &s1, &s2) == 3)
    {
        cout << "Setting background scalers: " << s1 << " to " << s2 << " every " << sn << " reads" << endl << endl;
        nScalRead = sn;
        back_bin_lo = s1;
        back_bin_hi = s2;
    }
    else
    {
        cout << "Background scalers not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

void	PTaggBack::ProcessScalerRead()
{
    PPhysics::ProcessScalerRead();

    if(clock_pre==0)
    {
        TaggerFirScal->Add(TaggerAccScal);
        GoosyNewFPD(TaggerFirScal);
    }

    nScalIter++;
    if(nScalIter==nScalRead)
    {
        TaggerCurScal->Add(TaggerAccScal,TaggerPreScal,1,-1);
        TaggerPreScal->Add(TaggerCurScal,1);
        GoosyNewFPD(TaggerCurScal);

        Double_t clock_now = LiveTimeScal->GetBinContent(1)/1e6;
        Double_t clock_new = clock_now-clock_pre;
        Double_t counts = TaggerCurScal->Integral(back_bin_lo,back_bin_hi);
        Double_t countrate = counts/clock_new;
        Double_t error = countrate*TMath::Sqrt((1/counts)+(1/(1e6*clock_new)));
        clock_pre = clock_now;
        cout << clock_now << "\t" << clock_new << "\t" << countrate << endl;
        Int_t nPoints = TaggBackground.GetN();
        TaggBackground.SetPoint(nPoints,clock_now,countrate);
        TaggBackground.SetPointError(nPoints,0,error);

        nScalIter=0;
    }
}

Bool_t	PTaggBack::Write()
{    
    // Write all GH1's and TObjects defined in this class
    if(!(GTreeManager::Write())) return false;

    TCanvas *c1 = new TCanvas("c1","c1");
    TaggBackground.Draw("AP");
    if(back_bin_lo==back_bin_hi) TaggBackground.SetTitle(Form("Tagger Channel %d",back_bin_lo-1));
    else TaggBackground.SetTitle(Form("Tagger Channels %d - %d",back_bin_lo-1,back_bin_hi-1));
    TaggBackground.GetXaxis()->SetTitle("Time (s)");
    TaggBackground.GetYaxis()->SetTitle("Background Rate (Hz)");

    Double_t xmin, xmax, ymin, ymax;
    TaggBackground.GetPoint(0,xmin,ymax);
    TaggBackground.GetPoint(TaggBackground.GetN()-1,xmax,ymin);
    cout << xmin << "\t" << ymin << "\t" << xmax << "\t" << ymax << endl;

    TF1 *f1 = new TF1("f1","[0]+[1]*TMath::Exp(-x/[2])",0,xmax);
    f1->SetParameter(0,ymin);
    f1->SetParLimits(0,0.5*ymin,1.5*ymin);
    f1->SetParameter(1,ymax-ymin);
    f1->SetParLimits(1,0.5*(ymax-ymin),1.5*(ymax-ymin));
    f1->SetParameter(2,1700);
    TaggBackground.Fit("f1");
    if(back_bin_lo==back_bin_hi) c1->SaveAs(Form("TaggBack_Chan%d.pdf",back_bin_lo-1));
    else c1->SaveAs(Form("TaggBack_Chan%d-%d.pdf",back_bin_lo-1,back_bin_hi-1));
    delete c1;
    delete f1;
    
    return true;
}
