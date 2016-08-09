#include <TH1.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
using namespace std;

void plotRMS1(string fileName)
{
	TCanvas can("can","can",1200,800);
	TFile f(fileName.c_str());
	gStyle->SetOptStat(0);
	TTree *t = (TTree*)f.Get("RMS");
	t->Draw("rms:chid>>hist(8257,0,8256,1000,0,6)");
	TH1F *h = (TH1F*)gDirectory->Get("hist");
	h->SetTitle(Form("Event %s;Channel ID;ADC RMS",(fileName.substr(19,fileName.length()-24)).c_str()));
	h->Draw("p0");
	can.Print(Form("rms_%s.png",(fileName.substr(19,fileName.length()-24)).c_str() ));
}

void plotFFT(string fileName)
{
	TCanvas can("can","can",1200,800);
	can.SetLogy();
	TFile f(fileName.c_str());
	gStyle->SetOptStat(0);
	gStyle->SetLegendBorderSize(0); 
	
	TH1F *UFFT = (TH1F*)f.Get("UFFT_mag");
	TH1F *UFFT_long = (TH1F*)f.Get("UlongwireFFT_mag");
	TH1F *mis = (TH1F*)f.Get("UFFT_mag_misconf");
	TH1F *VFFT = (TH1F*)f.Get("VFFT_mag");
	TH1F *VFFT_long = (TH1F*)f.Get("VlongwireFFT_mag");
	TH1F *WFFT = (TH1F*)f.Get("WFFT_mag");
	
	TTree* rms = (TTree*)f.Get("RMS");
	Int_t UN = 0, UlN = 0, misN = 0, VN = 0, VlN = 0, WN = 0;
	Bool_t chstat = 1;
	rms->SetBranchAddress("chstat",&chstat);
	for(Int_t c = 0; c<8256;c++)
	{
		rms->GetEntry(c);
		UN += (Int_t)chstat*(c<2400);
		UlN += (Int_t)chstat*(c>=673 && c<1729);
		misN += (Int_t)chstat*((c>=2016 && c<=2095) || (c>=2192 && c<=2302) || (c>=2352 && c<=2383));
		VN += (Int_t)chstat*(c>=2400 && c<4800);
		VlN += (Int_t)chstat*(c>=3073 && c<4129);
		WN += (Int_t)chstat*(c>=4800);
	}
	
	UFFT->Add(mis,1);
	UFFT->Scale(UN/UFFT->Integral());
	UFFT_long->Scale(UlN/UFFT_long->Integral());
	mis->Scale(misN/mis->Integral());
	VFFT->Scale(VN/VFFT->Integral());
	VFFT_long->Scale(VlN/VFFT_long->Integral());
	WFFT->Scale(WN/WFFT->Integral());
	
	/*UFFT->Add(mis,1);
	UFFT->Scale(2400/UFFT->Integral());
	UFFT_long->Scale(1056/UFFT_long->Integral());
	mis->Scale(223/mis->Integral());
	VFFT->Scale(2400/VFFT->Integral());
	VFFT_long->Scale(1056/VFFT_long->Integral());
	WFFT->Scale(3456/WFFT->Integral());*/
	
	UFFT->SetMarkerColor(kBlue);UFFT->SetFillColor(kBlue);UFFT->SetLineColor(kBlue);
	VFFT->SetMarkerColor(kRed);VFFT->SetFillColor(kRed);VFFT->SetLineColor(kRed);
	WFFT->SetMarkerColor(kCyan);WFFT->SetFillColor(kCyan);WFFT->SetLineColor(kCyan);
	VFFT_long->SetMarkerColor(kOrange);VFFT_long->SetFillColor(kOrange);VFFT_long->SetLineColor(kOrange);
	UFFT_long->SetMarkerColor(kViolet);UFFT_long->SetFillColor(kViolet);UFFT_long->SetLineColor(kViolet);
	mis->SetMarkerColor(kGreen);mis->SetFillColor(kGreen);mis->SetLineColor(kGreen);
	
	WFFT->SetMaximum(10);
	WFFT->SetMinimum(0.01);
	WFFT->GetXaxis()->SetRange(0,4800);
	WFFT->SetTitle(Form("Event %s;Frequency (kHz);FFT Magnitude",(fileName.substr(19,fileName.length()-24)).c_str()));
	for(Int_t bin = 500; bin < 4800; bin+=500) WFFT->GetXaxis()->SetBinLabel(bin,Form("%i",bin*208/1000));
	WFFT->GetXaxis()->LabelsOption("h");
	
	WFFT->Draw("p0");
	VFFT->Draw("p0 same");
	UFFT->Draw("p0 same");
	UFFT_long->Draw("p0 same");
	VFFT_long->Draw("p0 same");
	mis->Draw("p0 same");
	
	TLegend *leg = new TLegend(0.7,0.7,0.89,0.89);
   	leg->AddEntry(UFFT,"U","l");
   	leg->AddEntry(UFFT_long,"U (Long Wires)","l");
   	leg->AddEntry(mis,"U (Misconfigured)","l");
	leg->AddEntry(VFFT,"V","l");
   	leg->AddEntry(VFFT_long,"V (Long Wires)","l");
	leg->AddEntry(WFFT,"Y","l");
   	leg->Draw();
	
	//h->SetTitle(Form("Event %s;Channel ID;ADC RMS",(fileName.substr(19,fileName.length()-24)).c_str()));
	can.Print(Form("FFT_%s.png",(fileName.substr(19,fileName.length()-24)).c_str()));
}

void plotRMS(const char * eventsFileName)
{
	string rootFileName;
	ifstream file(eventsFileName);
	
	if(file.is_open()){
		while(!file.eof()){
			getline(file,rootFileName);
			if(rootFileName != "") {plotRMS1(rootFileName); plotFFT(rootFileName);}
		}
		file.close();
	}
}
