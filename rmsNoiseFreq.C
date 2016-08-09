#include <iostream>
#include <fstream>
#include "TH2.h"
#include "TH1.h"
#include "TTree.h"
#include "TF1.h"

/* Calculate RMS for particular channel (signal removal) */
Bool_t GetR(UShort_t c, Double_t& R, TH2F *h, UShort_t nticks)
{
	TH1F chan("chan","Noise RMS", 800, -100,100);
	Double_t RMS = 0;
	Double_t xq;
	Double_t par[3]; 
	
	for(UShort_t t = 1; t <= nticks; t++)
		chan.Fill(h->GetBinContent(c,t));
		
	xq = 0.5-.34;
	chan.GetQuantiles(1,&par[0],&xq);
	xq = 0.5;
	chan.GetQuantiles(1,&par[1],&xq);
	xq = 0.5+.34;
	chan.GetQuantiles(1,&par[2],&xq);
	
	RMS = sqrt((pow(par[2]-par[1],2)+pow(par[1]-par[0],2))/2);
	
	if(RMS < .1){R = 0; return 0;}
	
	R = RMS;
	return 1;
}

void rmsNoiseFreq(const char * filename)
{
	TCanvas can("can","can",400,400);
	TFile *data = new TFile(filename);
	TH2F *hu = (TH2F*)data->Get("hu_raw");
	TH2F *hv = (TH2F*)data->Get("hv_raw");
	TH2F *hw = (TH2F*)data->Get("hw_raw");
	TH2F *hu_orig = (TH2F*)data->Get("hu_orig");
	TH2F *hv_orig = (TH2F*)data->Get("hv_orig");
	TH2F *hw_orig = (TH2F*)data->Get("hw_orig");
	
	TTree *run = (TTree*)data->Get("Trun");
	run->Draw("eventNo");run->Draw("runNo");run->Draw("subRunNo");
	TLeaf *eventNo = (TLeaf*)run->GetLeaf("eventNo");
	TLeaf *runNo = (TLeaf*)run->GetLeaf("runNo");
	TLeaf *subRunNo = (TLeaf*)run->GetLeaf("subRunNo");
	
	Int_t eNo = 0, rNo = 0, sRNo = 0;
	if(eventNo) eNo = (Int_t)eventNo->GetValue();
	else cerr<<"Can't get event number. Set to zero"<<std::endl;
	if(runNo) rNo = (Int_t)runNo->GetValue();
	else cerr<<"Can't get run number. Set to zero"<<std::endl;
	if(subRunNo) sRNo = (Int_t)subRunNo->GetValue();
	else cerr<<"Can't get subrun number. Set to zero"<<std::endl;
	
	TH2F *h;
	UShort_t nticks = (UShort_t)hu->GetNbinsY();
	
	stringstream rootfName, textfName;
	rootfName<<"noiseData_"<<setfill('0')<<setw(7)<<rNo<<"-"<<setw(5)<<sRNo<<"_"<<setw(5)<<eNo<<".root";
	textfName<<"noiseData_"<<setfill('0')<<setw(7)<<rNo<<"-"<<setw(5)<<sRNo<<"_"<<setw(5)<<eNo<<".txt";
	TFile f(rootfName.str().c_str(),"RECREATE","RMS Tree File");

	Double_t rms;
	Bool_t chstat, outlier = 0, chirping = 0;
	Int_t chid, c, dead = 0;
	Int_t outlierChannels[10000] = {};
	
	TH1F rmsU1("rmsU1","rmsU1",673,1,673);
	TH1F rmsCUflat("rmsCUflat","rmsCUflat",1057,673,1729);
	TH1F rmsU2("rmsU2","rmsU2",673,1729,2401);
	TH1F rmsV1("rmsV1","rmsV1",673,2401,3073);
	TH1F rmsCVflat("rmsCVflat","rmsCVflat",1057,3073,4129);
	TH1F rmsV2("rmsV2","rmsV2",673,4129,4801);
	TH1F rmsCW("rmsCU","rmsCU",3456,4801,8257);
	
	TTree *tree = new TTree("RMS","RMS for all channels");
	tree->Branch("rms",&rms,"rms/D",60000);
	tree->Branch("chid",&chid,"chid/s");
	tree->Branch("chstat",&chstat,"chstat/O");
	tree->Branch("outlier_flag",&outlier,"outlier_flag/O");
	
	TTree *Trun = run->CloneTree();
	
	
	for(chid = 1; chid <= 8256; chid++)
	{
		h = hu; c = chid;
		if(c > 4800) {h = hw; c = c - 4800;}
		else if(c > 2400) {h = hv; c = c - 2400;}
		
		chstat = GetR(c,rms,h,nticks);
		tree->Fill();
		if(chstat == 0) dead++;
		else
		{
			if(chid<673) rmsU1.Fill(chid,rms);
			if(chid>672 && c<1729) rmsCUflat.Fill(chid,rms);
			if(chid>1728 && chid<2401) rmsU2.Fill(chid,rms);
			if(chid>2400 && chid<3073) rmsV1.Fill(chid,rms);
			if(chid>3072 && chid<4129) rmsCVflat.Fill(chid,rms);
			if(chid>4128 && chid<4801) rmsV2.Fill(chid,rms);
			if(chid>4800) rmsCW.Fill(chid,rms);
		}
	}
	
	tree->Draw("rms>>rmsUflat","(chid>672) && (chid<1729) && (chstat == 1)");
	tree->Draw("rms>>rmsVflat","(chid>3072) && (chid<4129) && (chstat == 1)");

	tree->Draw("rms>>rmsW","(chid>4800) && (chstat == 1)");
	tree->Draw("rms>>miscon1","(chid>2015) && (chid<2096) && (chstat == 1)");
	tree->Draw("rms>>miscon2","(chid>2191) && (chid<2303) && (chstat == 1)");
	tree->Draw("rms>>miscon3","(chid>2351) && (chid<2384) && (chstat == 1)");
	
	TH1F *rmsUflat = (TH1F*)gDirectory->Get("rmsUflat");
	TH1F *rmsVflat = (TH1F*)gDirectory->Get("rmsVflat");
	
	TH1F *rmsW = (TH1F*)gDirectory->Get("rmsW");
	TH1F *miscon1 = (TH1F*)gDirectory->Get("miscon1");
	TH1F *miscon2 = (TH1F*)gDirectory->Get("miscon2");
	TH1F *miscon3 = (TH1F*)gDirectory->Get("miscon3");
	
	/*W Plane*/
	Double_t meanW = rmsW->GetMean();
	Double_t sigW  = rmsW->GetStdDev();
	tree->Draw("rms>>oW",Form("(chid>4800) && (chstat == 1) && (fabs(rms-%f)>3*%f)",meanW,sigW));
	TH1F *oW = (TH1F*)gDirectory->Get("oW");
	Int_t outW = oW->GetEntries();
	
	/*Long Wire U Plane*/
	Double_t meanUflat = rmsUflat->GetMean();
	Double_t sigUflat  = rmsUflat->GetStdDev();
	tree->Draw("rms>>oUflat",Form("(chid>672) && (chid<1729) && (chstat == 1) && (fabs(rms-%f)>3*%f)",meanUflat,sigUflat));
	TH1F *oUflat = (TH1F*)gDirectory->Get("oUflat");
	Int_t outUflat = oUflat->GetEntries();
	
	/*Long Wire V Plane*/
	Double_t meanVflat = rmsVflat->GetMean();
	Double_t sigVflat  = rmsVflat->GetStdDev();
	tree->Draw("rms>>oVflat",Form("(chid>3072) && (chid<4129) && (chstat == 1) && (fabs(rms-%f)>3*%f)",meanVflat,sigVflat));
	TH1F *oVflat = (TH1F*)gDirectory->Get("oVflat");
	Int_t outVflat = oVflat->GetEntries();
	
	
	/*Varying Length Wire Induction Planes*/
	TF1* fitU1, *fitU2, *fitV1, *fitV2;
	TH1F res("res","res",100,-20,20);
	
	/*U Plane Upstream*/
	rmsU1.Fit("pol1");
	fitU1 = (TF1*)rmsU1.GetFunction("pol1");
	
	for(Int_t c = 1; c<673; c++)
		if(rmsU1.GetBinContent(c) > 0) res.Fill(fitU1->Eval(c)-rmsU1.GetBinContent(c));
		
	Double_t U1p0 = fitU1->GetParameter(0);
	Double_t U1p1 = fitU1->GetParameter(1);
	Double_t U1resRMS = res.GetRMS();
	Int_t U1out = 0;
	

	for(Int_t c = 1; c<673; c++)
		if(rmsU1.GetBinContent(c)>0 && fabs(fitU1->Eval(c)-rmsU1.GetBinContent(c))>3*U1resRMS) U1out++;
	res.Reset();
	
	/*U Plane Downstream*/
	rmsU2.Fit("pol1");
	fitU2 = (TF1*)rmsU2.GetFunction("pol1");
	
	for(Int_t c = 1729; c<2401; c++)
		if(rmsU2.GetBinContent(c-1728) > 0) res.Fill(fitU2->Eval(c)-rmsU2.GetBinContent(c-1728));
		
	Double_t U2p0 = fitU2->GetParameter(0);
	Double_t U2p1 = fitU2->GetParameter(1);
	Double_t U2resRMS = res.GetRMS();
	Int_t U2out = 0;
	
	for(Int_t c = 1729; c<2401; c++)
		if(rmsU2.GetBinContent(c-1728)>0 && fabs(fitU2->Eval(c)-rmsU2.GetBinContent(c-1728))>3*U2resRMS) U2out++;
	res.Reset();
	
	/*V Plane Upstream*/
	rmsV1.Fit("pol1");
	fitV1 = (TF1*)rmsV1.GetFunction("pol1");
	
	for(Int_t c = 2401; c<3073; c++)
		if(rmsV1.GetBinContent(c-2400) > 0) res.Fill(fitV1->Eval(c)-rmsV1.GetBinContent(c-2400));
		
	Double_t V1p0 = fitV1->GetParameter(0);
	Double_t V1p1 = fitV1->GetParameter(1);
	Double_t V1resRMS = res.GetRMS();
	Int_t V1out = 0;
	
	for(Int_t c = 2401; c<3073; c++)
		if(rmsV1.GetBinContent(c-2400)>0 && fabs(fitV1->Eval(c)-rmsV1.GetBinContent(c-2400))>3*V1resRMS) V1out++;
	res.Reset();
	
	/*V Plane Downstream*/
	rmsV2.Fit("pol1");
	fitV2 = (TF1*)rmsV2.GetFunction("pol1");
	
	for(Int_t c = 4129; c<4801; c++)
		if(rmsV2.GetBinContent(c-4128) > 0) res.Fill(fitV2->Eval(c)-rmsV2.GetBinContent(c-4128));
		
	Double_t V2p0 = fitV2->GetParameter(0);
	Double_t V2p1 = fitV2->GetParameter(1);
	Double_t V2resRMS = res.GetRMS();
	Int_t V2out = 0;
	
	for(Int_t c = 4129; c<4801; c++)
		if(rmsV2.GetBinContent(c-4128)>0 && fabs(fitV2->Eval(c)-rmsV2.GetBinContent(c-4128))>3*V2resRMS) V2out++;
	res.Reset();
	
	Double_t miscon1RMS = miscon1->GetRMS();
	Double_t miscon2RMS = miscon2->GetRMS();
	Double_t miscon3RMS = miscon3->GetRMS();
	Double_t miscon1Mean = miscon1->GetMean();
	Double_t miscon2Mean = miscon2->GetMean();
	Double_t miscon3Mean = miscon3->GetMean();
	
	for(Int_t c = 1; c<=8256; c++)
	{
		if(c<673 && rmsU1.GetBinContent(c)>0 && fabs(fitU1->Eval(c)-rmsU1.GetBinContent(c))>3*U1resRMS) outlierChannels[c-1] = 1;
		if(c>=673 && c<1729 && rmsCUflat.GetBinContent(c-672)>0 && fabs(meanUflat-rmsCUflat.GetBinContent(c-672))>3*sigUflat) outlierChannels[c-1] = 1;
		if(c>=1729 && c<2401 && rmsU2.GetBinContent(c-1728)>0 && fabs(fitU2->Eval(c)-rmsU2.GetBinContent(c-1728))>3*U2resRMS) outlierChannels[c-1] = 1;
		if(c>=2401 && c<3073 && rmsV1.GetBinContent(c-2400)>0 && fabs(fitV1->Eval(c)-rmsV1.GetBinContent(c-2400))>3*V1resRMS) outlierChannels[c-1] = 1;
		if(c>=3073 && c<4129 && rmsCVflat.GetBinContent(c-3072)>0 && fabs(meanVflat-rmsCVflat.GetBinContent(c-3072))>3*sigVflat) outlierChannels[c-1] = 1;
		if(c>=4129 && c<4801 && rmsV2.GetBinContent(c-4128)>0 && fabs(fitV2->Eval(c)-rmsV2.GetBinContent(c-4128))>3*V2resRMS) outlierChannels[c-1] = 1;
		if(c>=4801 && rmsCW.GetBinContent(c-4800)>0 && fabs(meanW-rmsCW.GetBinContent(c-4800))>3*sigW) outlierChannels[c-1] = 1;
	}
	
	Int_t nbins = hu_orig->GetNbinsY();
	TH1 *hum = new TH1F("UFFT_mag","UFFT_mag",nbins,1,nbins+1), *huml = new TH1F("UlongwireFFT_mag","UlongwireFFT_mag",nbins,1,nbins+1),
	    *humMis = new TH1F("UFFT_mag_misconf","UFFT_mag_misconf",nbins,1,nbins+1), 
	    *hvm = new TH1F("VFFT_mag","VFFT_mag",nbins,1,nbins+1), *hvml = new TH1F("VlongwireFFT_mag","VlongwireFFT_mag",nbins,1,nbins+1),
	    *hwm = new TH1F("WFFT_mag","WFFT_mag",nbins,1,nbins+1);
	    
	hum->GetXaxis()->SetTitle("Frequency"); hum->GetYaxis()->SetTitle("Magnitude");
	huml->GetXaxis()->SetTitle("Frequency"); huml->GetYaxis()->SetTitle("Magnitude");
	humMis->GetXaxis()->SetTitle("Frequency"); humMis->GetYaxis()->SetTitle("Magnitude");
	hvm->GetXaxis()->SetTitle("Frequency"); hvm->GetYaxis()->SetTitle("Magnitude");
	hvml->GetXaxis()->SetTitle("Frequency"); hvml->GetYaxis()->SetTitle("Magnitude");
	hwm->GetXaxis()->SetTitle("Frequency"); hwm->GetYaxis()->SetTitle("Magnitude");
	
	TH1 *ChanTemp = new TH1F("ChanTemp","ChanTemp",nbins, 1, nbins+1);
	TH1 *magTemp = 0;
	
	for(Int_t c = 1; c<2401; c++)
	{		
		for(Int_t t = 1; t <=nbins; t++)
			ChanTemp->Fill(t,hu_orig->GetBinContent(c,t));
		magTemp = ChanTemp->FFT(0,"MAG");
		
		for(Int_t i = 0; i< nbins;i++){
			Double_t rho = magTemp->GetBinContent(i+1);
			if(i==0) rho=0;
			if(!((c>=2016 && c<=2095) || (c>=2192 && c<=2302) || (c>=2352 && c<=2383))) hum->SetBinContent(i+1,rho+hum->GetBinContent(i+1));
			else humMis->SetBinContent(i+1,rho+humMis->GetBinContent(i+1));
			if(c>=673 && c<1729) huml->SetBinContent(i+1,rho+huml->GetBinContent(i+1));
		}
		ChanTemp->Reset(); delete magTemp;
	}
	for(Int_t c = 2401; c<4801; c++)
	{		
		for(Int_t t = 1; t <=nbins; t++)
			ChanTemp->Fill(t,hv_orig->GetBinContent(c-2400,t));
		magTemp = ChanTemp->FFT(0,"MAG");
		
		for(Int_t i = 0; i< nbins;i++){
			Double_t rho = magTemp->GetBinContent(i+1);
			if(i==0) rho=0;
			hvm->SetBinContent(i+1,rho+hvm->GetBinContent(i+1));
			if(c>=3073 && c<4129) hvml->SetBinContent(i+1,rho+hvml->GetBinContent(i+1));
		}
		ChanTemp->Reset(); delete magTemp;
	}
	for(Int_t c = 4801; c<8257; c++)
	{		
		for(Int_t t = 1; t <=nbins; t++)
			ChanTemp->Fill(t,hw_orig->GetBinContent(c-4800,t));
		magTemp = ChanTemp->FFT(0,"MAG");
		
		for(Int_t i = 0; i< nbins;i++){
			Double_t rho = magTemp->GetBinContent(i+1);
			if(i==0) rho=0;
			hwm->SetBinContent(i+1,rho+hwm->GetBinContent(i+1));
		}
		ChanTemp->Reset(); delete magTemp;
	}
	
	tree->Reset();
	for(chid = 1; chid <= 8256; chid++)
	{
		h = hu; c = chid;
		if(c > 4800) {h = hw; c = c - 4800;}
		else if(c > 2400) {h = hv; c = c - 2400;}
		
		chstat = GetR(c,rms,h,nticks);
		outlier = outlierChannels[chid-1];
		tree->Fill();
	}
	
	f.Write();
	
	f.Delete("rmsU1;1");
	f.Delete("rmsCUflat;1");
	f.Delete("rmsU2;1");
	f.Delete("rmsV1;1");
	f.Delete("rmsCVflat;1");
	f.Delete("rmsV2;1");
	f.Delete("rmsCU;1");
	f.Delete("rmsUflat;1");
	f.Delete("rmsVflat;1");
	f.Delete("rmsW;1");
	f.Delete("miscon1;1");
	f.Delete("miscon2;1");
	f.Delete("miscon3;1");
	f.Delete("oW;1");
	f.Delete("oVflat;1");
	f.Delete("oUflat;1");
	f.Delete("res;1");
	f.Delete("ChanTemp;1");
	f.Close();
	
	ofstream outFile;
	outFile.open(textfName.str().c_str());
	outFile<<rNo<<"\t"<<sRNo<<"\t"<<eNo<<std::endl;
	outFile<<U1p0<<"\t"<<U1p1<<std::endl;
	outFile<<U1resRMS<<std::endl;
	outFile<<U1out<<std::endl;
	
	outFile<<meanUflat<<std::endl;
	outFile<<sigUflat<<std::endl;
	outFile<<outUflat<<std::endl;

	outFile<<U2p0<<"\t"<<U2p1<<std::endl;
	outFile<<U2resRMS<<std::endl;
	outFile<<U2out<<std::endl;

	outFile<<V1p0<<"\t"<<V1p1<<std::endl;
	outFile<<V1resRMS<<std::endl;
	outFile<<V1out<<std::endl;
	
	outFile<<meanVflat<<std::endl;
	outFile<<sigVflat<<std::endl;
	outFile<<outVflat<<std::endl;

	outFile<<V2p0<<"\t"<<V2p1<<std::endl;
	outFile<<V2resRMS<<std::endl;
	outFile<<V2out<<std::endl;
	
	outFile<<meanW<<std::endl;
	outFile<<sigW<<std::endl;
	outFile<<outW<<std::endl;
	
	outFile<<miscon1Mean<<"\t"<<miscon1RMS<<std::endl;
	outFile<<miscon2Mean<<"\t"<<miscon2RMS<<std::endl;
	outFile<<miscon3Mean<<"\t"<<miscon3RMS<<std::endl;
	
	outFile<<dead<<std::endl;
	
	outFile.close();
	can.Close();
}
