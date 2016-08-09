#include <iostream>
#include <fstream>
#include <TTree.h>
#include <TH1.h>
using namespace std;

void plotDeadChannelsRun(const char * fileName, const Int_t nEvents)
{
	TCanvas can("can","can",1200,800);
	gStyle->SetOptStat(0);
	gStyle->SetErrorX(0.0001);
	gStyle->SetLegendBorderSize(0); 
	string rootFileName;
	ifstream file(fileName);
	TTree *tree;
	TFile *f = 0;
	Int_t fileCount = 0;
	vector<TTree*> eventTrees(nEvents);
	vector<TFile*> eventFiles(nEvents);
	vector<Int_t> run(nEvents), subrun(nEvents), event(nEvents);
	
	//Save trees to eventTrees
	if(file.is_open()){
		while(!file.eof() && fileCount<nEvents){
			getline(file,rootFileName);
			eventFiles.at(fileCount) = TFile::Open(rootFileName.c_str(),"update");
			tree = (TTree*)eventFiles.at(fileCount)->Get("RMS");
			eventTrees.at(fileCount) = tree;
			if(eventTrees.at(fileCount)->GetEntry(0) == -1) cerr<<"Could not read TTree RMS in file "<<rootFileName.c_str()<<endl;
			tree = (TTree*)eventFiles.at(fileCount)->Get("Trun");
			tree->Draw("runNo");tree->Draw("subRunNo");tree->Draw("eventNo");
			run.at(fileCount) = tree->GetLeaf("runNo")->GetValue();
			subrun.at(fileCount) = tree->GetLeaf("subRunNo")->GetValue();
			event.at(fileCount) = tree->GetLeaf("eventNo")->GetValue();
			eventTrees.at(fileCount)->SetDirectory(0);
			delete eventFiles.at(fileCount++);
		}
		file.close();
	}	
	if(fileCount != nEvents) cerr<<"Mismatch: number of events and number of input root files"<<endl;
	
	//Create and set branch addresses
	Double_t rms;
	UShort_t chid;
	Bool_t chstat, outlier_flag, chirping_flag;
	
	for(Int_t i = 0; i < nEvents; i++){
		eventTrees.at(i)->SetBranchAddress("rms",&rms);
		eventTrees.at(i)->SetBranchAddress("chid",&chid);
		eventTrees.at(i)->SetBranchAddress("chstat",&chstat);
		eventTrees.at(i)->SetBranchAddress("outlier_flag",&outlier_flag);
		//eventTrees.at(i)->SetBranchAddress("chirping_flag",&chirping_flag);
		}
	
	//Declare relevant plots and parameters
	Int_t firstRun = *min_element(run.begin(),run.end()), lastRun = *max_element(run.begin(),run.end());
	vector<TH1F*> deadChannels(3);
	deadChannels.at(0) = new TH1F("deadChannelsU","Dead Channels in U Plane",(lastRun-firstRun+1),firstRun,lastRun+1);
	deadChannels.at(1) = new TH1F("deadChannelsV","Dead Channels in V Plane",(lastRun-firstRun+1),firstRun,lastRun+1);
	deadChannels.at(2) = new TH1F("deadChannelsW","Dead Channels in W Plane",(lastRun-firstRun+1),firstRun,lastRun+1);
	
	Double_t norm = 0;
	//Fill Plots
	for(Int_t n = 0; n < nEvents; n++){
		norm = count(run.begin(),run.end(),run.at(n));
		if (norm == 0.0) cout<<n<<" : "<<run.at(n)<<endl;
		deadChannels.at(0)->Fill(run.at(n), eventTrees.at(n)->Draw("","chstat==0&&chid<2401")/norm);
		deadChannels.at(1)->Fill(run.at(n), eventTrees.at(n)->Draw("","chstat==0&&chid<4801&&chid>2400")/norm);
		deadChannels.at(2)->Fill(run.at(n), eventTrees.at(n)->Draw("","chstat==0&&chid>4800")/norm);
		/*for(Int_t c = 0; c < 8256; c++){
			chstat = 1;
			eventTrees.at(n)->GetEntry(c);
			if(chstat == 1 || norm == 0) continue;
			if(c < 2401) deadChannels.at(0)->Fill(run.at(n), 1/norm);
			if(c>2400 && c<4801) deadChannels.at(1)->Fill(run.at(n), 1/norm);
			if(c>4800) deadChannels.at(2)->Fill(run.at(n), 1/norm);
		}*/
	}

	deadChannels.at(0)->SetMarkerStyle(kFullDotMedium);
	deadChannels.at(1)->SetMarkerStyle(kFullDotMedium);
	deadChannels.at(2)->SetMarkerStyle(kFullDotMedium);
	deadChannels.at(0)->SetMarkerColor(kBlue);
	deadChannels.at(1)->SetMarkerColor(kRed);
	deadChannels.at(2)->SetMarkerColor(kCyan);
	deadChannels.at(0)->SetTitle(";Run Number;Number of Dead Channels");
	deadChannels.at(0)->SetMaximum(850);
	deadChannels.at(0)->GetYaxis()->SetTitleOffset(1.5);
	deadChannels.at(0)->Draw("P HIST");
	deadChannels.at(1)->Draw("P HIST same");
	deadChannels.at(2)->Draw("P HIST same");
	
	TLegend *leg = new TLegend(0.7,0.75,0.89,0.89);
   	leg->AddEntry(deadChannels.at(0),"U Plane","p");
   	leg->AddEntry(deadChannels.at(1),"V Plane","p");
   	leg->AddEntry(deadChannels.at(2),"Y Plane","p");
   	leg->Draw();
   
	can.Print(Form("deadChannels_%i-%i_%i.png",firstRun,lastRun,nEvents));	
	TFile outF(Form("deadChannels_%i-%i_%i.root",firstRun,lastRun,nEvents),"RECREATE");
	deadChannels.at(0)->Write();
	deadChannels.at(1)->Write();
	deadChannels.at(2)->Write();
	outF.Close();
}

void plotDeadChannelsTime(const char * fileName, const Int_t nEvents)
{
	TCanvas can("can","can",1200,800);
	gStyle->SetOptStat(0);
	gStyle->SetErrorX(0.0001);
	gStyle->SetLegendBorderSize(0); 
	string rootFileName;
	ifstream file(fileName);
	TTree *tree;
	TFile *f = 0;
	Int_t fileCount = 0;
	vector<TTree*> eventTrees(nEvents);
	vector<TFile*> eventFiles(nEvents);
	vector<Int_t> run(nEvents), subrun(nEvents), event(nEvents), eventTime(nEvents);
	
	//Save trees to eventTrees
	if(file.is_open()){
		while(!file.eof() && fileCount<nEvents){
			getline(file,rootFileName);
			eventFiles.at(fileCount) = TFile::Open(rootFileName.c_str(),"update");
			tree = (TTree*)eventFiles.at(fileCount)->Get("RMS");
			eventTrees.at(fileCount) = tree;
			if(eventTrees.at(fileCount)->GetEntry(0) == -1) cerr<<"Could not read TTree RMS in file "<<rootFileName.c_str()<<endl;
			tree = (TTree*)eventFiles.at(fileCount)->Get("Trun");
			tree->Draw("runNo");tree->Draw("subRunNo");tree->Draw("eventNo");tree->Draw("eventTime");
			run.at(fileCount) = tree->GetLeaf("runNo")->GetValue();
			subrun.at(fileCount) = tree->GetLeaf("subRunNo")->GetValue();
			event.at(fileCount) = tree->GetLeaf("eventNo")->GetValue();
			eventTime.at(fileCount) = tree->GetLeaf("eventTime")->GetValue();
			eventTrees.at(fileCount)->SetDirectory(0);
			delete eventFiles.at(fileCount++);
		}
		file.close();
	}	
	if(fileCount != nEvents) cerr<<"Mismatch: number of events and number of input root files"<<endl;
	
	//Create and set branch addresses
	Double_t rms;
	UShort_t chid;
	Bool_t chstat, outlier_flag, chirping_flag;
	
	for(Int_t i = 0; i < nEvents; i++){
		eventTrees.at(i)->SetBranchAddress("rms",&rms);
		eventTrees.at(i)->SetBranchAddress("chid",&chid);
		eventTrees.at(i)->SetBranchAddress("chstat",&chstat);
		eventTrees.at(i)->SetBranchAddress("outlier_flag",&outlier_flag);
		//eventTrees.at(i)->SetBranchAddress("chirping_flag",&chirping_flag);
		}
	
	//Declare relevant plots and parameters // 2/10/16 @ 5:30 (UTC)
	Int_t firstRun = *min_element(run.begin(),run.end()), lastRun = *max_element(run.begin(),run.end());
	Double_t firstTime = *min_element(eventTime.begin(),eventTime.end()), lastTime = *max_element(eventTime.begin(),eventTime.end());
	
	vector<TH1F*> deadChannels(3);
	deadChannels.at(0) = new TH1F("deadChannelsU","Dead Channels in U Plane",(lastTime-firstTime+1),0,lastTime-firstTime+1);
	deadChannels.at(1) = new TH1F("deadChannelsV","Dead Channels in V Plane",(lastTime-firstTime+1),0,lastTime-firstTime+1);
	deadChannels.at(2) = new TH1F("deadChannelsW","Dead Channels in W Plane",(lastTime-firstTime+1),0,lastTime-firstTime+1);
	vector<TH1F*> deadChannelsFeb(3), deadChannelsAp(3), deadChannelsJun(3);
	
	Double_t norm = 0;
	//Fill Plots
	for(Int_t n = 0; n < nEvents; n++){
		for(Int_t c = 0; c < 8256; c++){
			eventTrees.at(n)->GetEntry(c);
			if(c < 2401) deadChannels.at(0)->Fill(eventTime.at(n)-firstTime, deadChannels.at(0)->GetBinContent(eventTime.at(n)-firstTime)+(1-chstat)/2400);
			if(c>2400 && c<4801) deadChannels.at(1)->Fill(eventTime.at(n)-firstTime, deadChannels.at(1)->GetBinContent(eventTime.at(n)-firstTime)+(1-chstat)/2400);
			if(c>4800) deadChannels.at(2)->Fill(eventTime.at(n)-firstTime, deadChannels.at(2)->GetBinContent(eventTime.at(n)-firstTime)+(1-chstat)/3456);
			}
		}

	deadChannels.at(0)->SetMarkerStyle(kFullDotMedium);
	deadChannels.at(1)->SetMarkerStyle(kFullDotMedium);
	deadChannels.at(2)->SetMarkerStyle(kFullDotMedium);
	deadChannels.at(0)->SetMarkerColor(kBlue);
	deadChannels.at(1)->SetMarkerColor(kRed);
	deadChannels.at(2)->SetMarkerColor(kCyan);
	deadChannels.at(0)->SetTitle(";Time Since 2-10-16 @ 5:30 (UTC);Average Proportion of Dead Channels");
	deadChannels.at(0)->SetMaximum(0.5);
	deadChannels.at(0)->GetYaxis()->SetTitleOffset(1.5);
	deadChannels.at(0)->Draw("P HIST");
	deadChannels.at(1)->Draw("P HIST same");
	deadChannels.at(2)->Draw("P HIST same");
	
	deadChannelsFeb.at(0) = (TH1F*)deadChannels.at(0)->Clone();
	deadChannelsAp.at(0) = (TH1F*)deadChannels.at(0)->Clone();
	deadChannelsJun.at(0) = (TH1F*)deadChannels.at(0)->Clone();
	deadChannelsFeb.at(1) = (TH1F*)deadChannels.at(0)->Clone();
	deadChannelsAp.at(1) = (TH1F*)deadChannels.at(0)->Clone();
	deadChannelsJun.at(1) = (TH1F*)deadChannels.at(0)->Clone();
	deadChannelsFeb.at(2) = (TH1F*)deadChannels.at(0)->Clone();
	deadChannelsAp.at(2) = (TH1F*)deadChannels.at(0)->Clone();
	deadChannelsJun.at(2) = (TH1F*)deadChannels.at(0)->Clone();
	deadChannelsFeb.at(0)->SetFillColor(33);
	deadChannelsAp.at(0)->SetFillColor(33);
	deadChannelsJun.at(0)->SetFillColor(33);
	deadChannelsFeb.at(1)->SetFillColor(33);
	deadChannelsAp.at(1)->SetFillColor(33);
	deadChannelsJun.at(1)->SetFillColor(33);
	deadChannelsFeb.at(2)->SetFillColor(33);
	deadChannelsAp.at(2)->SetFillColor(33);
	deadChannelsJun.at(2)->SetFillColor(33);
	
	deadChannelsFeb.at(0)->GetXaxis()->SetRange(0,1456790400);
	deadChannelsAp.at(0)->GetXaxis()->SetRange(1459468800,1462060800);
	deadChannelsJun.at(0)->GetXaxis()->SetRange(1464739200,1467331200);
	deadChannelsFeb.at(1)->GetXaxis()->SetRange(0,1456790400);
	deadChannelsAp.at(1)->GetXaxis()->SetRange(1459468800,1462060800);
	deadChannelsJun.at(1)->GetXaxis()->SetRange(1464739200,1467331200);
	deadChannelsFeb.at(2)->GetXaxis()->SetRange(0,1456790400);
	deadChannelsAp.at(2)->GetXaxis()->SetRange(1459468800,1462060800);
	deadChannelsJun.at(2)->GetXaxis()->SetRange(1464739200,1467331200);
	
	deadChannelsFeb.at(0)->Draw("P HIST same");
	deadChannelsFeb.at(1)->Draw("P HIST same");
	deadChannelsFeb.at(2)->Draw("P HIST same");
	deadChannelsAp.at(0)->Draw("P HIST same");
	deadChannelsAp.at(1)->Draw("P HIST same");
	deadChannelsAp.at(2)->Draw("P HIST same");
	deadChannelsAp.at(0)->Draw("P HIST same");
	deadChannelsAp.at(1)->Draw("P HIST same");
	deadChannelsAp.at(2)->Draw("P HIST same");
	
	TLegend *leg = new TLegend(0.7,0.7,0.89,0.89);
   	leg->AddEntry(deadChannels.at(0),"U Plane","p");
   	leg->AddEntry(deadChannels.at(1),"V Plane","p");
   	leg->AddEntry(deadChannels.at(2),"Y Plane","p");
   	leg->Draw();
   
	can.Print(Form("deadChannelsT_%i-%i_%i.png",firstRun,lastRun,nEvents));	
	TFile outF(Form("deadChannelsT_%i-%i_%i.root",firstRun,lastRun,nEvents),"RECREATE");
	deadChannels.at(0)->Write();
	deadChannels.at(1)->Write();
	deadChannels.at(2)->Write();
	outF.Close();
}

void plotTotDeadChannels(const char * fileName, const Int_t nEvents)
{
	TCanvas can("can","can",1200,800);
	gStyle->SetOptStat(0);
	gStyle->SetErrorX(0.0001);
	string rootFileName;
	ifstream file(fileName);
	TTree *tree;
	TFile *f = 0;
	Int_t fileCount = 0;
	vector<TTree*> eventTrees(nEvents);
	vector<TFile*> eventFiles(nEvents);
	vector<Int_t> run(nEvents), subrun(nEvents), event(nEvents);
	
	//Save trees to eventTrees
	if(file.is_open()){
		while(!file.eof() && fileCount<nEvents){
			getline(file,rootFileName);
			eventFiles.at(fileCount) = TFile::Open(rootFileName.c_str(),"update");
			tree = (TTree*)eventFiles.at(fileCount)->Get("RMS");
			eventTrees.at(fileCount) = tree;
			if(eventTrees.at(fileCount)->GetEntry(0) == -1) cerr<<"Could not read TTree RMS in file "<<rootFileName.c_str()<<endl;
			tree = (TTree*)eventFiles.at(fileCount)->Get("Trun");
			tree->Draw("runNo");tree->Draw("subRunNo");tree->Draw("eventNo");
			run.at(fileCount) = tree->GetLeaf("runNo")->GetValue();
			subrun.at(fileCount) = tree->GetLeaf("subRunNo")->GetValue();
			event.at(fileCount) = tree->GetLeaf("eventNo")->GetValue();
			eventTrees.at(fileCount)->SetDirectory(0);
			delete eventFiles.at(fileCount++);
		}
		file.close();
	}	
	if(fileCount != nEvents) cerr<<"Mismatch: number of events and number of input root files"<<endl;
	
	//Create and set branch addresses
	Double_t rms;
	UShort_t chid;
	Bool_t chstat, outlier_flag, chirping_flag;
	
	for(Int_t i = 0; i < nEvents; i++){
		eventTrees.at(i)->SetBranchAddress("rms",&rms);
		eventTrees.at(i)->SetBranchAddress("chid",&chid);
		eventTrees.at(i)->SetBranchAddress("chstat",&chstat);
		eventTrees.at(i)->SetBranchAddress("outlier_flag",&outlier_flag);
		//eventTrees.at(i)->SetBranchAddress("chirping_flag",&chirping_flag);
		}
	
	//Declare relevant plots and parameters
	Int_t firstRun = *min_element(run.begin(),run.end()), lastRun = *max_element(run.begin(),run.end());
	TH1F *deadChannels = new TH1F("deadChannels","Dead Channels",(lastRun-firstRun+1),firstRun,lastRun+1);
	
	Double_t norm = 0;
	//Fill Plots
	for(Int_t n = 0; n < nEvents; n++){
		norm = count(run.begin(),run.end(),run.at(n));
		for(Int_t c = 0; c < 8256; c++){
			eventTrees.at(n)->GetEntry(c);
			deadChannels->Fill(run.at(n), deadChannels->GetBinContent(run.at(n))+(1-chstat)/norm);
			}
		}

	deadChannels->SetMarkerStyle(kFullDotMedium);
	deadChannels->SetTitle(";Run Number;Average Total Number of Dead Channels");
	deadChannels->GetYaxis()->SetTitleOffset(1.5);
	deadChannels->SetMinimum(1400);
	deadChannels->Draw("P HIST");
   
	can.Print(Form("totalDeadChannels_%i-%i_%i.png",firstRun,lastRun,nEvents));	
}

void plotTextMetrics(const char * fileName)
{
	TCanvas can("can","can",1200,800);
	gStyle->SetOptStat(0);
	gStyle->SetErrorX(0.0001);
	TFile f(fileName);
	TTree *metricTree = (TTree*)f.Get("metricTree");
	Int_t run, firstRun, lastRun;
	metricTree->SetBranchAddress("run",&run);
	Int_t nEntries = metricTree->GetEntries();
	metricTree->GetEntry(0); firstRun = run;
	metricTree->GetEntry(nEntries-1); lastRun = run;
	
	//U Plane Upstream
	metricTree->Draw("U1p0:run>>U1p0");
	TH1F *U1p0= (TH1F*)gDirectory->Get("U1p0");
	U1p0->SetTitle("Upstream U Plane;Run Number;Zero Order Linear Fit Parameter");
	U1p0->SetMarkerStyle(kFullDotSmall);
	U1p0->Draw("P HIST");
	can.Print(Form("U1p0_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("U1p1:run>>U1p1");
	TH1F *U1p1= (TH1F*)gDirectory->Get("U1p1");
	U1p1->SetTitle("Upstream U Plane;Run Number;First Order Linear Fit Parameter");
	U1p1->SetMarkerStyle(kFullDotSmall);
	U1p1->Draw("P HIST");
	can.Print(Form("U1p1_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("U1rms:run>>U1rms");
	TH1F *U1rms= (TH1F*)gDirectory->Get("U1rms");
	U1rms->SetTitle("Upstream U Plane;Run Number;Fit Residuals RMS");
	U1rms->SetMarkerStyle(kFullDotSmall);
	U1rms->Draw("P HIST");
	can.Print(Form("U1rms_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("U1outlier:run>>U1out");
	TH1F *U1out= (TH1F*)gDirectory->Get("U1out");
	U1out->SetTitle("Upstream U Plane;Run Number;Outliers of Fit");
	U1out->SetMarkerStyle(kFullDotSmall);
	U1out->Draw("P HIST");
	can.Print(Form("U1outlier_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	//U Plane Long Wire
	metricTree->Draw("Umean:run>>Umean");
	TH1F *Umean= (TH1F*)gDirectory->Get("Umean");
	Umean->SetTitle("Long Wire U Plane;Run Number;Average");
	Umean->SetMarkerStyle(kFullDotSmall);
	Umean->Draw("P HIST");
	can.Print(Form("Umean_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("Urms:run>>Urms");
	TH1F *Urms= (TH1F*)gDirectory->Get("Urms");
	Urms->SetTitle("Long Wire U Plane;Run Number;RMS");
	Urms->SetMarkerStyle(kFullDotSmall);
	Urms->Draw("P HIST");
	can.Print(Form("Urms_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("Uoutlier:run>>Uout");
	TH1F *Uout= (TH1F*)gDirectory->Get("Uout");
	Uout->SetTitle("Long Wire U Plane;Run Number;Outliers");
	Uout->SetMarkerStyle(kFullDotSmall);
	Uout->Draw("P HIST");
	can.Print(Form("Uoutlier_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	//U Plane Downstream
	metricTree->Draw("U2p0:run>>U2p0");
	TH1F *U2p0= (TH1F*)gDirectory->Get("U2p0");
	U2p0->SetTitle("Downstream U Plane;Run Number;Zero Order Linear Fit Parameter");
	U2p0->SetMarkerStyle(kFullDotSmall);
	U2p0->Draw("P HIST");
	can.Print(Form("U2p0_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("U2p1:run>>U2p1");
	TH1F *U2p1= (TH1F*)gDirectory->Get("U2p1");
	U2p1->SetTitle("Downstream U Plane;Run Number;First Order Linear Fit Parameter");
	U2p1->SetMarkerStyle(kFullDotSmall);
	U2p1->Draw("P HIST");
	can.Print(Form("U2p1_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("U2rms:run>>U2rms");
	TH1F *U2rms= (TH1F*)gDirectory->Get("U2rms");
	U2rms->SetTitle("Downstream U Plane;Run Number;Fit Residuals RMS");
	U2rms->SetMarkerStyle(kFullDotSmall);
	U2rms->Draw("P HIST");
	can.Print(Form("U2rms_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("U2outlier:run>>U2out");
	TH1F *U2out= (TH1F*)gDirectory->Get("U2out");
	U2out->SetTitle("Downstream U Plane;Run Number;Outliers of Fit");
	U2out->SetMarkerStyle(kFullDotSmall);
	U2out->Draw("P HIST");
	can.Print(Form("U2outlier_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	//V Plane Upstream
	metricTree->Draw("V1p0:run>>V1p0");
	TH1F *V1p0= (TH1F*)gDirectory->Get("V1p0");
	V1p0->SetTitle("Upstream V Plane;Run Number;Zero Order Linear Fit Parameter");
	V1p0->SetMarkerStyle(kFullDotSmall);
	V1p0->Draw("P HIST");
	can.Print(Form("V1p0_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("V1p1:run>>V1p1");
	TH1F *V1p1= (TH1F*)gDirectory->Get("V1p1");
	V1p1->SetTitle("Upstream V Plane;Run Number;First Order Linear Fit Parameter");
	V1p1->SetMarkerStyle(kFullDotSmall);
	V1p1->Draw("P HIST");
	can.Print(Form("V1p1_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("V1rms:run>>V1rms");
	TH1F *V1rms= (TH1F*)gDirectory->Get("V1rms");
	V1rms->SetTitle("Upstream V Plane;Run Number;Fit Residuals RMS");
	V1rms->SetMarkerStyle(kFullDotSmall);
	V1rms->Draw("P HIST");
	can.Print(Form("V1rms_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("V1outlier:run>>V1out");
	TH1F *V1out= (TH1F*)gDirectory->Get("V1out");
	V1out->SetTitle("Upstream V Plane;Run Number;Outliers of Fit");
	V1out->SetMarkerStyle(kFullDotSmall);
	V1out->Draw("P HIST");
	can.Print(Form("V1outlier_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	//V Plane Long Wire
	metricTree->Draw("Vmean:run>>Vmean");
	TH1F *Vmean= (TH1F*)gDirectory->Get("Vmean");
	Vmean->SetTitle("Long Wire V Plane;Run Number;Average");
	Vmean->SetMarkerStyle(kFullDotSmall);
	Vmean->Draw("P HIST");
	can.Print(Form("Vmean_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("Vrms:run>>Vrms");
	TH1F *Vrms= (TH1F*)gDirectory->Get("Vrms");
	Vrms->SetTitle("Long Wire V Plane;Run Number;RMS");
	Vrms->SetMarkerStyle(kFullDotSmall);
	Vrms->Draw("P HIST");
	can.Print(Form("Vrms_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("Voutlier:run>>Vout");
	TH1F *Vout= (TH1F*)gDirectory->Get("Vout");
	Vout->SetTitle("Long Wire V Plane;Run Number;Outliers");
	Vout->SetMarkerStyle(kFullDotSmall);
	Vout->Draw("P HIST");
	can.Print(Form("Voutlier_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	//V Plane Downstream
	metricTree->Draw("V2p0:run>>V2p0");
	TH1F *V2p0= (TH1F*)gDirectory->Get("V2p0");
	V2p0->SetTitle("Downstream V Plane;Run Number;Zero Order Linear Fit Parameter");
	V2p0->SetMarkerStyle(kFullDotSmall);
	V2p0->Draw("P HIST");
	can.Print(Form("V2p0_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("V2p1:run>>V2p1");
	TH1F *V2p1= (TH1F*)gDirectory->Get("V2p1");
	V2p1->SetTitle("Downstream V Plane;Run Number;First Order Linear Fit Parameter");
	V2p1->SetMarkerStyle(kFullDotSmall);
	V2p1->Draw("P HIST");
	can.Print(Form("V12p1_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("V2rms:run>>V2rms");
	TH1F *V2rms= (TH1F*)gDirectory->Get("V2rms");
	V2rms->SetTitle("Downstream V Plane;Run Number;Fit Residuals RMS");
	V2rms->SetMarkerStyle(kFullDotSmall);
	V2rms->Draw("P HIST");
	can.Print(Form("V2rms_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("V2outlier:run>>V2out");
	TH1F *V2out= (TH1F*)gDirectory->Get("V2out");
	V2out->SetTitle("Downstream V Plane;Run Number;Outliers of Fit");
	V2out->SetMarkerStyle(kFullDotSmall);
	V2out->Draw("P HIST");
	can.Print(Form("V2outlier_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	//W Plane
	metricTree->Draw("Wmean:run>>Wmean");
	TH1F *Wmean= (TH1F*)gDirectory->Get("Wmean");
	Wmean->SetTitle("W Plane;Run Number;Average");
	Wmean->SetMarkerStyle(kFullDotSmall);
	Wmean->Draw("P HIST");
	can.Print(Form("Wmean_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("Wrms:run>>Wrms");
	TH1F *Wrms= (TH1F*)gDirectory->Get("Wrms");
	Wrms->SetTitle("W Plane;Run Number;RMS");
	Wrms->SetMarkerStyle(kFullDotSmall);
	Wrms->Draw("P HIST");
	can.Print(Form("Wrms_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("Woutlier:run>>Wout");
	TH1F *Wout= (TH1F*)gDirectory->Get("Wout");
	Wout->SetTitle("Long Wire W Plane;Run Number;Outliers");
	Wout->SetMarkerStyle(kFullDotSmall);
	Wout->Draw("P HIST");
	can.Print(Form("Woutlier_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	
	//Misconfigured Channels
	metricTree->Draw("mis1mean:run>>mis1mean");
	TH1F *mis1mean= (TH1F*)gDirectory->Get("mis1mean");
	mis1mean->SetTitle("Channels 2016-2095;Run Number;Average");
	mis1mean->SetMarkerStyle(kFullDotSmall);
	mis1mean->Draw("P HIST");
	can.Print(Form("mis1mean_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("mis1rms:run>>mis1rms");
	TH1F *mis1rms= (TH1F*)gDirectory->Get("mis1rms");
	mis1rms->SetTitle("Channels 2016-2095;Run Number;RMS");
	mis1rms->SetMarkerStyle(kFullDotSmall);
	mis1rms->Draw("P HIST");
	can.Print(Form("mis1rms_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("mis2mean:run>>mis2mean");
	TH1F *mis2mean= (TH1F*)gDirectory->Get("mis2mean");
	mis2mean->SetTitle("Channels 2192-2302;Run Number;Average");
	mis2mean->SetMarkerStyle(kFullDotSmall);
	mis2mean->Draw("P HIST");
	can.Print(Form("mis2mean_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("mis2rms:run>>mis2rms");
	TH1F *mis2rms= (TH1F*)gDirectory->Get("mis2rms");
	mis2rms->SetTitle("Channels 2192-2302;Run Number;RMS");
	mis2rms->SetMarkerStyle(kFullDotSmall);
	mis2rms->Draw("P HIST");
	can.Print(Form("mis2rms_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("mis3mean:run>>mis3mean");
	TH1F *mis3mean= (TH1F*)gDirectory->Get("mis3mean");
	mis3mean->SetTitle("Channels 2352-2383;Run Number;Average");
	mis3mean->SetMarkerStyle(kFullDotSmall);
	mis3mean->Draw("P HIST");
	can.Print(Form("mis3mean_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("mis3rms:run>>mis3rms");
	TH1F *mis3rms= (TH1F*)gDirectory->Get("mis3rms");
	mis3rms->SetTitle("Channels 2352-2383;Run Number;RMS");
	mis3rms->SetMarkerStyle(kFullDotSmall);
	mis3rms->Draw("P HIST");
	can.Print(Form("mis3rms_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	metricTree->Draw("numDead:run>>numDead");
	TH1F *numDead= (TH1F*)gDirectory->Get("numDead");
	numDead->SetTitle(";Run Number;Number of Dead Channels");
	numDead->SetMarkerStyle(kFullDotSmall);
	numDead->Draw("P HIST");
	can.Print(Form("numDead_%i-%i_%i.png",firstRun,lastRun,nEntries));
	
	TFile ff(Form("metricHists_%i-%i_%i.root",firstRun,lastRun,nEntries),"RECREATE");
	U1p0->Write();
	U1p1->Write();
	U1rms->Write();
	U1out->Write();
	Umean->Write();
	Urms->Write();
	Uout->Write();
	U2p0->Write();
	U2p1->Write();
	U2rms->Write();
	U2out->Write();
	V1p0->Write();
	V1p1->Write();
	V1rms->Write();
	V1out->Write();
	Vmean->Write();
	Vrms->Write();
	Vout->Write();
	V2p0->Write();
	V2p1->Write();
	V2rms->Write();
	V2out->Write();
	Wmean->Write();
	Wrms->Write();
	Wout->Write();
	mis1mean->Write();
	mis1rms->Write();
	mis2mean->Write();
	mis2rms->Write();
	mis3mean->Write();
	mis3rms->Write();
	numDead->Write();
	ff.Close();
}

void plotFFTMetrics(const char * FFTfileName)
{
	TCanvas can("can","can",1200,800);
	gStyle->SetOptStat(0);
	gStyle->SetLegendBorderSize(0); 
	
	TFile f(FFTfileName);
	vector<TTree*> trees(6);
	trees.at(0) = (TTree*)f.Get("UFFT_mag");
	trees.at(1) = (TTree*)f.Get("UlongwireFFT_mag");
	trees.at(2) = (TTree*)f.Get("UFFT_mag_misconf");
	trees.at(3) = (TTree*)f.Get("VFFT_mag");
	trees.at(4) = (TTree*)f.Get("VlongwireFFT_mag");
	trees.at(5) = (TTree*)f.Get("WFFT_mag");
	TTree *stats = (TTree*)f.Get("runStats");
	nEvents = (Int_t)stats->GetEntries();
	
	Double_t mlow, mmid, mhigh, mco, m36, m108, m900;
	for(Int_t i = 0; i<6; i++)
	{
		trees.at(i)->SetBranchAddress("mlow",&mlow);
		trees.at(i)->SetBranchAddress("mmid",&mmid);
		trees.at(i)->SetBranchAddress("mhigh",&mhigh);
		trees.at(i)->SetBranchAddress("mco",&mco);
		trees.at(i)->SetBranchAddress("m36",&m36);
		trees.at(i)->SetBranchAddress("m108",&m108);
		trees.at(i)->SetBranchAddress("m900",&m900);
	}
	
	stats->GetEntry(0);
	Int_t firstRun = stats->GetLeaf("run")->GetValue(0);
	stats->GetEntry(nEvents-1);
	Int_t lastRun = stats->GetLeaf("run")->GetValue(0);
	
	Int_t run = 0;
	stats->SetBranchAddress("run",&run);
	
	vector<TH1F*> hlow(6);
	vector<TH1F*> hmid(6);
	vector<TH1F*> hhigh(6);
	vector<TH1F*> hco(6);
	vector<TH1F*> h36(6);
	vector<TH1F*> h108(6);
	vector<TH1F*> h900(6);
	
	for(Int_t i = 0; i<6; i++)
	{
		hlow.at(i) = new TH1F(Form("hlow%i",i),"",lastRun-firstRun+1,firstRun,lastRun);
		hmid.at(i) = new TH1F(Form("hmid%i",i),"",lastRun-firstRun+1,firstRun,lastRun);
		hhigh.at(i) = new TH1F(Form("hhigh%i",i),"",lastRun-firstRun+1,firstRun,lastRun);
		hco.at(i) = new TH1F(Form("hco%i",i),"",lastRun-firstRun+1,firstRun,lastRun);
		h36.at(i) = new TH1F(Form("h36%i",i),"",lastRun-firstRun+1,firstRun,lastRun);
		h108.at(i) = new TH1F(Form("h108%i",i),"",lastRun-firstRun+1,firstRun,lastRun);
		h900.at(i) = new TH1F(Form("h900%i",i),"",lastRun-firstRun+1,firstRun,lastRun);
	}
	
	Int_t Nrun, r;
	for(Int_t n = 0; n<nEvents; n++)
	{
		stats->GetEntry(n);
		r = run;
		Nrun = stats->Draw("",Form("run==%i",r),"goff"); 
		for(Int_t i = 0; i< 6; i++)
		{
			trees.at(i)->GetEntry(n);
			hlow.at(i)->Fill(r,mlow/Nrun);
			hmid.at(i)->Fill(r,mmid/Nrun);
			hhigh.at(i)->Fill(r,mhigh/Nrun);
			hco.at(i)->Fill(r,mco/Nrun);
			h36.at(i)->Fill(r,m36/Nrun);
			h108.at(i)->Fill(r,m108/Nrun);
			h900.at(i)->Fill(r,m900/Nrun);
		}
	}
	
	hlow.at(0)->SetMarkerColor(kBlue);hlow.at(0)->SetFillColor(kBlue);hlow.at(0)->SetLineColor(kBlue);
	hmid.at(0)->SetMarkerColor(kBlue);hmid.at(0)->SetFillColor(kBlue);hmid.at(0)->SetLineColor(kBlue);
	hhigh.at(0)->SetMarkerColor(kBlue);hhigh.at(0)->SetFillColor(kBlue);hhigh.at(0)->SetLineColor(kBlue);
	hco.at(0)->SetMarkerColor(kBlue);hco.at(0)->SetFillColor(kBlue);hco.at(0)->SetLineColor(kBlue);
	h36.at(0)->SetMarkerColor(kBlue);h36.at(0)->SetFillColor(kBlue);h36.at(0)->SetLineColor(kBlue);
	h108.at(0)->SetMarkerColor(kBlue);h108.at(0)->SetFillColor(kBlue);h108.at(0)->SetLineColor(kBlue);
	h900.at(0)->SetMarkerColor(kBlue);h900.at(0)->SetFillColor(kBlue);h900.at(0)->SetLineColor(kBlue);
	
	hlow.at(1)->SetMarkerColor(kViolet);hlow.at(1)->SetFillColor(kViolet);hlow.at(1)->SetLineColor(kViolet);
	hmid.at(1)->SetMarkerColor(kViolet);hmid.at(1)->SetFillColor(kViolet);hmid.at(1)->SetLineColor(kViolet);
	hhigh.at(1)->SetMarkerColor(kViolet);hhigh.at(1)->SetFillColor(kViolet);hhigh.at(1)->SetLineColor(kViolet);
	hco.at(1)->SetMarkerColor(kViolet);hco.at(1)->SetFillColor(kViolet);hco.at(1)->SetLineColor(kViolet);
	h36.at(1)->SetMarkerColor(kViolet);h36.at(1)->SetFillColor(kViolet);h36.at(1)->SetLineColor(kViolet);
	h108.at(1)->SetMarkerColor(kViolet);h108.at(1)->SetFillColor(kViolet);h108.at(1)->SetLineColor(kViolet);
	h900.at(1)->SetMarkerColor(kViolet);h900.at(1)->SetFillColor(kViolet);h900.at(1)->SetLineColor(kViolet);
	
	hlow.at(2)->SetMarkerColor(kGreen);hlow.at(2)->SetFillColor(kGreen);hlow.at(2)->SetLineColor(kGreen);
	hmid.at(2)->SetMarkerColor(kGreen);hmid.at(2)->SetFillColor(kGreen);hmid.at(2)->SetLineColor(kGreen);
	hhigh.at(2)->SetMarkerColor(kGreen);hhigh.at(2)->SetFillColor(kGreen);hhigh.at(2)->SetLineColor(kGreen);
	hco.at(2)->SetMarkerColor(kGreen);hco.at(2)->SetFillColor(kGreen);hco.at(2)->SetLineColor(kGreen);
	h36.at(2)->SetMarkerColor(kGreen);h36.at(2)->SetFillColor(kGreen);h36.at(2)->SetLineColor(kGreen);
	h108.at(2)->SetMarkerColor(kGreen);h108.at(2)->SetFillColor(kGreen);h108.at(2)->SetLineColor(kGreen);
	h900.at(2)->SetMarkerColor(kGreen);h900.at(2)->SetFillColor(kGreen);h900.at(2)->SetLineColor(kGreen);
	
	hlow.at(3)->SetMarkerColor(kRed);hlow.at(3)->SetFillColor(kRed);hlow.at(3)->SetLineColor(kRed);
	hmid.at(3)->SetMarkerColor(kRed);hmid.at(3)->SetFillColor(kRed);hmid.at(3)->SetLineColor(kRed);
	hhigh.at(3)->SetMarkerColor(kRed);hhigh.at(3)->SetFillColor(kRed);hhigh.at(3)->SetLineColor(kRed);
	hco.at(3)->SetMarkerColor(kRed);hco.at(3)->SetFillColor(kRed);hco.at(3)->SetLineColor(kRed);
	h36.at(3)->SetMarkerColor(kRed);h36.at(3)->SetFillColor(kRed);h36.at(3)->SetLineColor(kRed);
	h108.at(3)->SetMarkerColor(kRed);h108.at(3)->SetFillColor(kRed);h108.at(3)->SetLineColor(kRed);
	h900.at(3)->SetMarkerColor(kRed);h900.at(3)->SetFillColor(kRed);h900.at(3)->SetLineColor(kRed);
	
	hlow.at(4)->SetMarkerColor(kOrange);hlow.at(4)->SetFillColor(kOrange);hlow.at(4)->SetLineColor(kOrange);
	hmid.at(4)->SetMarkerColor(kOrange);hmid.at(4)->SetFillColor(kOrange);hmid.at(4)->SetLineColor(kOrange);
	hhigh.at(4)->SetMarkerColor(kOrange);hhigh.at(4)->SetFillColor(kOrange);hhigh.at(4)->SetLineColor(kOrange);
	hco.at(4)->SetMarkerColor(kOrange);hco.at(4)->SetFillColor(kOrange);hco.at(4)->SetLineColor(kOrange);
	h36.at(4)->SetMarkerColor(kOrange);h36.at(4)->SetFillColor(kOrange);h36.at(4)->SetLineColor(kOrange);
	h108.at(4)->SetMarkerColor(kOrange);h108.at(4)->SetFillColor(kOrange);h108.at(4)->SetLineColor(kOrange);
	h900.at(4)->SetMarkerColor(kOrange);h900.at(4)->SetFillColor(kOrange);h900.at(4)->SetLineColor(kOrange);
	
	hlow.at(5)->SetMarkerColor(kCyan);hlow.at(5)->SetFillColor(kCyan);hlow.at(5)->SetLineColor(kCyan);
	hmid.at(5)->SetMarkerColor(kCyan);hmid.at(5)->SetFillColor(kCyan);hmid.at(5)->SetLineColor(kCyan);
	hhigh.at(5)->SetMarkerColor(kCyan);hhigh.at(5)->SetFillColor(kCyan);hhigh.at(5)->SetLineColor(kCyan);
	hco.at(5)->SetMarkerColor(kCyan);hco.at(5)->SetFillColor(kCyan);hco.at(5)->SetLineColor(kCyan);
	h36.at(5)->SetMarkerColor(kCyan);h36.at(5)->SetFillColor(kCyan);h36.at(5)->SetLineColor(kCyan);
	h108.at(5)->SetMarkerColor(kCyan);h108.at(5)->SetFillColor(kCyan);h108.at(5)->SetLineColor(kCyan);
	h900.at(5)->SetMarkerColor(kCyan);h900.at(5)->SetFillColor(kCyan);h900.at(5)->SetLineColor(kCyan);
	
	for(Int_t i = 0;i<6;i++){
		hlow.at(i)->SetMarkerStyle(kFullDotMedium);
		hmid.at(i)->SetMarkerStyle(kFullDotMedium);
		hhigh.at(i)->SetMarkerStyle(kFullDotMedium);
		hco.at(i)->SetMarkerStyle(kFullDotMedium);
		h36.at(i)->SetMarkerStyle(kFullDotMedium);
		h108.at(i)->SetMarkerStyle(kFullDotMedium);
		h900.at(i)->SetMarkerStyle(kFullDotMedium);
	}
	TLegend *leg = new TLegend(0.7,0.7,0.89,0.89);
   	leg->AddEntry(hlow.at(0),"U","l");
   	leg->AddEntry(hlow.at(1),"U (Long Wires)","l");
   	leg->AddEntry(hlow.at(2),"U (Misconfigured)","l");
	leg->AddEntry(hlow.at(3),"V","l");
   	leg->AddEntry(hlow.at(4),"V (Long Wires)","l");
	leg->AddEntry(hlow.at(5),"Y","l");
	leg->SetFillStyle(0);
	
	vector<Int_t> mlowMax(6), mmidMax(6), mhighMax(6), mcoMax(6), m36Max(6), m108Max(6), m900Max(6);
	
	for(Int_t i = 0; i<6;i++)
	{
		mlowMax.at(i) = hlow.at(i)->GetMaximum();
		mmidMax.at(i) = hmid.at(i)->GetMaximum();
		mhighMax.at(i) = hhigh.at(i)->GetMaximum();
		mcoMax.at(i) = hco.at(i)->GetMaximum();
		m36Max.at(i) = h36.at(i)->GetMaximum();
		m108Max.at(i) = h108.at(i)->GetMaximum();
		m900Max.at(i) = h900.at(i)->GetMaximum();
	}
	
	//hlow.at(0)->SetMaximum(*max_element(mlowMax.begin(),mlowMax.end()));
	hlow.at(0)->SetMaximum(40);
	hlow.at(0)->SetTitle(";Run Number;Average of 42-104 kHz Magnitude");
	hlow.at(0)->Draw("P HIST");
	for(Int_t i = 1; i<6; i++) hlow.at(i)->Draw("P HIST SAME");
	leg->Draw();
	can.Print(Form("hlow_%i-%i_%i.png",firstRun,lastRun,nEvents));
	
	//hmid.at(0)->SetMaximum(*max_element(mmidMax.begin(),mmidMax.end()));
	hmid.at(0)->SetMaximum(20);
	hmid.at(0)->SetTitle(";Run Number;Average of 108-208 kHz Magnitude");
	hmid.at(0)->Draw("P HIST");
	for(Int_t i = 1; i<6; i++) hmid.at(i)->Draw("P HIST SAME");
	leg->Draw();
	can.Print(Form("hmid_%i-%i_%i.png",firstRun,lastRun,nEvents));
	
	//hhigh.at(0)->SetMaximum(*max_element(mhighMax.begin(),mhighMax.end()));
	hhigh.at(0)->SetMaximum(3);
	hhigh.at(0)->SetTitle(";Run Number;Average of 667-833 kHz Magnitude");
	hhigh.at(0)->Draw("P HIST");
	for(Int_t i = 1; i<6; i++) hhigh.at(i)->Draw("P HIST SAME");
	leg->Draw();
	can.Print(Form("hhigh_%i-%i_%i.png",firstRun,lastRun,nEvents));
	
	//hco.at(0)->SetMaximum(*max_element(mcoMax.begin(),mcoMax.end()));
	hco.at(0)->SetMaximum(80);
	hco.at(0)->SetTitle(";Run Number;Average of 10-30 kHz Magnitude");
	hco.at(0)->Draw("P HIST");
	for(Int_t i = 1; i<6; i++) hco.at(i)->Draw("P HIST SAME");
	leg->Draw();
	can.Print(Form("hco_%i-%i_%i.png",firstRun,lastRun,nEvents));
	
	//h36.at(0)->SetMaximum(*max_element(m36Max.begin(),m36Max.end()));
	h36.at(0)->SetMaximum(100);
	h36.at(0)->SetTitle(";Run Number;Average of 36+/-0.83 kHz Magnitude");
	h36.at(0)->Draw("P HIST");
	for(Int_t i = 1; i<6; i++) h36.at(i)->Draw("P HIST SAME");
	leg->Draw();
	can.Print(Form("h36_%i-%i_%i.png",firstRun,lastRun,nEvents));
	
	//h108.at(0)->SetMaximum(*max_element(m108Max.begin(),m108Max.end()));
	h108.at(0)->SetMaximum(35);
	h108.at(0)->SetTitle(";Run Number;Average of 108+/-0.83 kHz Magnitude");
	h108.at(0)->Draw("P HIST");
	for(Int_t i = 1; i<6; i++) h108.at(i)->Draw("P HIST SAME");
	leg->Draw();
	can.Print(Form("h108_%i-%i_%i.png",firstRun,lastRun,nEvents));
	
	//h900.at(0)->SetMaximum(*max_element(m900Max.begin(),m900Max.end()));
	h900.at(0)->SetMaximum(5);
	h900.at(0)->SetTitle(";Run Number;Average of 900+/-42 kHz Magnitude");
	h900.at(0)->Draw("P HIST");
	for(Int_t i = 1; i<6; i++) h900.at(i)->Draw("P HIST SAME");
	leg->Draw();
	can.Print(Form("h900_%i-%i_%i.png",firstRun,lastRun,nEvents));
	
	f.Close();
}
