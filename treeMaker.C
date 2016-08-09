#include <TH1.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
using namespace std;

void getStats(TH1F * h, Double_t &mlow, Double_t &mmid, Double_t &mhigh, Double_t &mco, Double_t &m36, Double_t &m108, Double_t &m900, Int_t N)
{
	mlow = 0; mmid = 0; mhigh = 0; mco = 0; m36 = 0; m108 = 0; m900 = 0;
	for(Int_t bin = 200; bin<500;bin++)
		mlow += h->GetBinContent(bin)/300./N/208.3;
	for(Int_t bin = 500; bin<1000;bin++)
		mmid += h->GetBinContent(bin)/500./N/208.3;
	for(Int_t bin = 3200; bin<4000;bin++)
		mhigh += h->GetBinContent(bin)/800./N/208.3;
	for(Int_t bin = 48; bin<145;bin++)
		mco += h->GetBinContent(bin)/97./N/208.3;
	for(Int_t bin = 168; bin<177;bin++)
		m36 += h->GetBinContent(bin)/9./N/208.3;
	for(Int_t bin = 514;bin<523;bin++)
		m108 += h->GetBinContent(518)/9./N/208.3;
	for(Int_t bin = 4120; bin<4521;bin++)
		m900 += h->GetBinContent(bin)/400./N/208.3;
}

void makeFFTMetrics(const char * eventsFileName, const Int_t nEvents)
{
	TCanvas can("can","can",700,700);
	TFile FFTFile(Form("FFTMetrics_%i.root",nEvents),"RECREATE","Trees Containing FFT Metrics");
	TFile *f; TTree *RUN;
	
	vector<TTree*> trees(6);
	trees.at(0) = new TTree("UFFT_mag","U Plane FFT Magnitude Metrics");
	trees.at(1) = new TTree("UlongwireFFT_mag","U Plane (long wire) FFT Magnitude Metrics");
	trees.at(2) = new TTree("UFFT_mag_misconf","U Plane Misconfigured Channels FFT Magnitude Metrics");
	trees.at(3) = new TTree("VFFT_mag","V Plane FFT Magnitude Metrics");
	trees.at(4) = new TTree("VlongwireFFT_mag","V Plane (long wire) FFT Magnitude Metrics");
	trees.at(5) = new TTree("WFFT_mag","W Plane FFT Magnitude Metrics");
	TTree* runStats = new TTree("runStats","Run subrub event numbers for corresponding tree elements");
	
	Double_t mlow[6] = {}, mmid[6] = {}, mhigh[6] = {}, mco[6] = {}, m36[6] = {}, m108[6] = {}, m900[6] = {};
	for(Int_t i = 0; i < 6; i++)
	{
		trees.at(i)->Branch("mlow",&mlow[i],"mlow/D",60000);
		trees.at(i)->Branch("mmid",&mmid[i],"mmid/D",60000);
		trees.at(i)->Branch("mhigh",&mhigh[i],"mhigh/D",60000);
		trees.at(i)->Branch("mco",&mco[i],"mco/D",60000);
		trees.at(i)->Branch("m36",&m36[i],"m36/D",60000);
		trees.at(i)->Branch("m108",&m108[i],"m108/D",60000);
		trees.at(i)->Branch("m900",&m900[i],"m900/D",60000);
	}
	Int_t run, subrun, event, eventTime;
	runStats->Branch("run",&run,"run/I");
	runStats->Branch("subrun",&subrun,"subrun/I");
	runStats->Branch("event",&event,"event/I");
	//runStats->Branch("eventTime",&eventTime,"eventTime/D");
	
	string rootFileName;
	ifstream file(eventsFileName);
	Int_t fileCount = 0;
	Int_t NU, NUl, Nmis, NV, NVl, NW;
	TTree* RMS;
	
	if(file.is_open()){
		while(!file.eof() && fileCount<nEvents){
			getline(file,rootFileName);
			f = TFile::Open(rootFileName.c_str());
			RMS = (TTree*)f->Get("RMS");
			NU = RMS->Draw("","chid<2400 && chstat==1","goff");
			NUl = RMS->Draw("","chid>=673 && chid<1729 && chstat==1","goff");
			Nmis = RMS->Draw("","((chid>=2016 && chid<=2095) || (chid>=2192 && chid<=2302) || (chid>=2352 && chid<=2383)) && chstat==1","goff");
			NV = RMS->Draw("","chid>=2400 && chid<4800 && chstat==1","goff");
			NVl = RMS->Draw("","chid>=3073 && chid<4129 && chstat==1","goff");
			NW = RMS->Draw("","chid>4800 && chstat==1","goff");
			
			getStats((TH1F*)f->Get("UFFT_mag"),mlow[0],mmid[0],mhigh[0],mco[0],m36[0],m108[0],m900[0],NU);
			trees.at(0)->Fill();
			getStats((TH1F*)f->Get("UlongwireFFT_mag"),mlow[1],mmid[1],mhigh[1],mco[1],m36[1],m108[1],m900[1],NUl);
			trees.at(1)->Fill();
			getStats((TH1F*)f->Get("UFFT_mag_misconf"),mlow[2],mmid[2],mhigh[2],mco[2],m36[2],m108[2],m900[2],Nmis);
			trees.at(2)->Fill();
			getStats((TH1F*)f->Get("VFFT_mag"),mlow[3],mmid[3],mhigh[3],mco[3],m36[3],m108[3],m900[3],NV);
			trees.at(3)->Fill();
			getStats((TH1F*)f->Get("VlongwireFFT_mag"),mlow[4],mmid[4],mhigh[4],mco[4],m36[4],m108[4],m900[4],NVl);
			trees.at(4)->Fill();
			getStats((TH1F*)f->Get("WFFT_mag"),mlow[5],mmid[5],mhigh[5],mco[5],m36[5],m108[5],m900[5],NW);
			trees.at(5)->Fill();
			
			RUN = (TTree*)f->Get("Trun");
			RUN->Draw("runNo");RUN->Draw("subRunNo");RUN->Draw("eventNo");//RUN->Draw("eventTime");
			run = RUN->GetLeaf("runNo")->GetValue();
			subrun = RUN->GetLeaf("subRunNo")->GetValue();
			event = RUN->GetLeaf("eventNo")->GetValue();
			//eventTime = RUN->GetLeaf("eventTime")->GetValue();
			runStats->Fill();
			can.Close();
			f->Close();
			fileCount++;
		}
		file.close();
	}
	if(fileCount != nEvents) cerr<<"Mismatch in number of files read"<<endl;
	
	FFTFile.Write();
	can.Close();
}

void chirpID(const char * fileName, const Int_t start, const Int_t nEvents = 500)
{
	TCanvas can("can","can",700,700);
	string rootFileName;
	ifstream file(fileName);
	TTree *tree;
	TFile *f = 0;
	Int_t fileCount = 0, fileAt = 0;
	vector<TTree*> eventTrees(nEvents);
	vector<TFile*> eventFiles(nEvents);
	vector<string> eventFileNames(nEvents);
	vector<TBranch*> eventBranches(nEvents);
	vector<Bool_t> chstats(nEvents);
	Bool_t chstat = 0, chirping_flag = 0;
	
	if(file.is_open()){
		while(!file.eof() && fileCount<nEvents){
			if(fileAt<start) {fileAt++; continue;}
			getline(file,rootFileName);
			eventFileNames.at(fileCount) = rootFileName;
			eventFiles.at(fileCount) = TFile::Open(rootFileName.c_str(),"update");
			tree = (TTree*)eventFiles.at(fileCount)->Get("RMS");
			eventTrees.at(fileCount) = tree;
			eventTrees.at(fileCount)->SetDirectory(0);
			delete eventFiles.at(fileCount);
			eventTrees.at(fileCount++)->GetEntry(0);
		}
		file.close();
	}	
	if(fileCount != nEvents) cerr<<"Mismatch: number of events and number of input root files"<<endl;
	
	cout<<"Files have been read"<<endl;
	
	for(Int_t i = 0; i < nEvents; i++){
		if(i%(nEvents/5) == 0) cout<<"Branch addresses are being set for event "<<i<<endl;
		if(eventTrees.at(i)->GetBranch("chirping_flag")) {
  			eventTrees.at(i)->SetBranchAddress("chirping_flag",&chirping_flag);
			eventBranches.at(i) = eventTrees.at(i)->GetBranch("chirping_flag");
			eventBranches.at(i)->Reset();
			}
		else eventBranches.at(i) = eventTrees.at(i)->Branch("chirping_flag",&chirping_flag,"chirping_flag/O");
		eventTrees.at(i)->SetBranchAddress("chstat",&chstat);
		//eventTrees.at(i)->SetEntries(8256);
	}

	for(Int_t i = 0; i < 8256; i++){
		if(i%500 == 0) cout<<"Channel "<<i<<" being evaluated for chirping"<<endl;
		for(Int_t j = 0; j < nEvents; j++){
			eventTrees.at(j)->GetEntry(i);
			chstats.at(j) = chstat;
		}	
		for(Int_t j = 0; j < nEvents; j++){
			if(chstats.at(j) == 1) continue;
			eventTrees.at(j)->GetEntry(i);
			if(find(chstats.begin(), chstats.end(), 1) != chstats.end()) chirping_flag = 1;
			else chirping_flag = 0;
			eventBranches.at(j)->Fill();
		}
	}
	
	for(Int_t n = 0; n < nEvents; n++){
		eventFiles.at(n) = TFile::Open(eventFileNames.at(n).c_str(),"update");
		eventFiles.at(n)->Delete("RMS;1");
		eventTrees.at(n)->Write();
		delete eventFiles.at(n);
	}	
}

void addTime(const char * fileName, const char * timeInfo, const Int_t nEvents)
{
	TCanvas can("can","can",700,700);
	string rootFileName;
	ifstream file(fileName);
	string line;
	ifstream timeFile(timeInfo);
	TTree *tree;
	TFile *f = 0;
	Int_t fileCount = 0;
	vector<TTree*> eventTrees(nEvents);
	vector<TFile*> eventFiles(nEvents);
	vector<string> eventFileNames(nEvents);
	vector<TBranch*> eventBranches(nEvents);
	Int_t runNo, subRunNo, eventNo;
	Double_t eventTime;
	vector<Int_t> RUN(nEvents), SUBRUN(nEvents), EVENT(nEvents);
	vector<Double_t> TIME(nEvents);
	
	if(file.is_open()){
		while(!file.eof() && fileCount<nEvents){
			getline(file,rootFileName);
			eventFileNames.at(fileCount) = rootFileName;
			eventFiles.at(fileCount) = TFile::Open(rootFileName.c_str(),"update");
			tree = (TTree*)eventFiles.at(fileCount)->Get("Trun");
			eventTrees.at(fileCount) = tree;
			eventTrees.at(fileCount)->Draw("runNo");
			eventTrees.at(fileCount)->Draw("subRunNo");
			eventTrees.at(fileCount)->Draw("eventNo");
			eventTrees.at(fileCount)->SetDirectory(0);
			delete eventFiles.at(fileCount);
			fileCount++;
		}
		file.close();
	}	
	if(fileCount != nEvents) cerr<<"Mismatch: number of events and number of input root files"<<endl;
	
	cout<<"Files have been read"<<endl;
	
	for(Int_t i = 0; i < nEvents; i++){
		if(i%(nEvents/5) == 0) cout<<"Branch addresses are being set for event "<<i<<endl;
		if(eventTrees.at(i)->GetBranch("eventTime")) {
  			eventTrees.at(i)->SetBranchAddress("eventTime",&eventTime);
			eventBranches.at(i) = eventTrees.at(i)->GetBranch("eventTime");
			eventBranches.at(i)->Reset();
			}
		else eventBranches.at(i) = eventTrees.at(i)->Branch("eventTime",&eventTime,"eventTime/D");
		//eventTrees.at(i)->SetBranchAddress("runNo",&runNo);
		//eventTrees.at(i)->SetBranchAddress("subRunNo",&subRunNo);
		//eventTrees.at(i)->SetBranchAddress("eventNo",&eventNo);
	}
	
	Int_t eventCount = 0;
	if(timeFile.is_open()){
		while (getline(timeFile, line))
    		{
        		istringstream ss(line);
			ss >> RUN.at(eventCount) >> SUBRUN.at(eventCount) >> EVENT.at(eventCount) >> TIME.at(eventCount);
			eventCount++;
		}
	}
	
	
	std::vector<Int_t>::iterator runIt;
	Int_t e = -1;
	
	for(Int_t i = 0; i< nEvents; i++)
	{
		//cout<<eventTrees.at(i)->GetEntry(0)<<endl;
		eventTrees.at(i)->Draw("runNo"); eventTrees.at(i)->Draw("subRunNo"); eventTrees.at(i)->Draw("eventNo");
		runNo = eventTrees.at(i)->GetLeaf("runNo")->GetValue();
		subRunNo = eventTrees.at(i)->GetLeaf("subRunNo")->GetValue();
		eventNo = eventTrees.at(i)->GetLeaf("eventNo")->GetValue();
		
		runIt = find (RUN.begin(), RUN.end(), runNo);
		for(Int_t j = runIt-RUN.begin(); j< nEvents; j++) 
		{	
			if(RUN.at(j) == runNo && SUBRUN.at(j) == subRunNo && EVENT.at(j) == eventNo){ e = j; break;}
		}
		if(e == -1) {cerr<<"Couldn't find time information for event "<<runNo<<"\t"<<subRunNo<<"\t"<<eventNo<<endl; continue;}
		eventTime = TIME.at(e);
		eventBranches.at(i)->Fill();
		e = -1;
	}
	
	for(Int_t n = 0; n < nEvents; n++){
		eventFiles.at(n) = TFile::Open(eventFileNames.at(n).c_str(),"update");
		eventFiles.at(n)->Delete("Trun;1");
		eventTrees.at(n)->Write();
		delete eventFiles.at(n);
	}	
}

void makeMetricTree(const char * eventFileName, const Int_t nEvents)
{
	string textFileName, line;
	ifstream inputFile(eventFileName);
	Int_t fileCount = 0, lineCount = 0;
	Int_t i1,i2,i3;
	Double_t d1, d2, eventTime;
	TFile *rf;
	TTree *ft;
	
	TFile metricRootFile(Form("RMSMetrics_%i.root",nEvents),"RECREATE","Trees Containing RMS Metrics");
	TTree *metricTree = new TTree("metricTree","Noise Metrics");
	//Run stats
	TBranch *brun = metricTree->Branch("run",&i1,"run/I");
	TBranch *bsubrun = metricTree->Branch("subrun",&i2,"subrun/I");
	TBranch *bevent = metricTree->Branch("event",&i3,"event/I");
	TBranch *beventTime = metricTree->Branch("eventTime",&eventTime,"eventTime/D");
	//U plane
	TBranch *bU1p0 = metricTree->Branch("U1p0",&d1,"U1p0/D");			//0th order parameter in linear fit of rms in upstream U plane channels
	TBranch *bU1p1 = metricTree->Branch("U1p1",&d2,"U1p1/D");			//1st order parameter in linear fit of rms in upstream U plane channels
	TBranch *bU1rms = metricTree->Branch("U1rms",&d1,"U1rms/D");			//RMS of residuals of fit in upstream U plane channels
	TBranch *bU1outlier = metricTree->Branch("U1outlier",&i1,"U1outlier/I");	//Number of outliers in upstream U plane channels based on 3RMS cut
	TBranch *bUmean = metricTree->Branch("Umean",&d1,"Umean/D");			//Average of RMS in U plane long wires
	TBranch *bUrms = metricTree->Branch("Urms",&d1,"Urms/D");			//Standard deviation of RMS in U plane long wires
	TBranch *bUoutlier = metricTree->Branch("Uoutlier",&i1,"Uoutlier/I");		//Number of outliers in U plane long wires with 3RMS cut
	TBranch *bU2p0 = metricTree->Branch("U2p0",&d1,"U2p0/D");			//0th order parameter in linear fit of rms in downstream U plane channels
	TBranch *bU2p1 = metricTree->Branch("U2p1",&d2,"U2p1/D");			//1st order parameter in linear fit of rms in downstream U plane channels
	TBranch *bU2rms = metricTree->Branch("U2rms",&d1,"U2rms/D");			//RMS of residuals of fit in downstream U plane channels
	TBranch *bU2outlier = metricTree->Branch("U2outlier",&i1,"U2outlier/I");	//Number of outliers in downstream U plane channels based on 2RMS cut
	//V plane
	TBranch *bV1p0 = metricTree->Branch("V1p0",&d1,"V1p0/D");			//0th order parameter in linear fit of rms in upstream U plane channels
	TBranch *bV1p1 = metricTree->Branch("V1p1",&d2,"V1p1/D");			//1st order parameter in linear fit of rms in upstream U plane channels
	TBranch *bV1rms = metricTree->Branch("V1rms",&d1,"V1rms/D");			//RMS of residuals of fit in upstream U plane channels
	TBranch *bV1outlier = metricTree->Branch("V1outlier",&i1,"V1outlier/I");	//Number of outliers in upstream U plane channels based on 3RMS cut
	TBranch *bVmean = metricTree->Branch("Vmean",&d1,"Vmean/D");			//Average of RMS in U plane long wires
	TBranch *bVrms = metricTree->Branch("Vrms",&d1,"Vrms/D");			//Standard deviation of RMS in U plane long wires
	TBranch *bVoutlier = metricTree->Branch("Voutlier",&i1,"Voutlier/I");		//Number of outliers in U plane long wires with 3RMS cut
	TBranch *bV2p0 = metricTree->Branch("V2p0",&d1,"V2p0/D");			//0th order parameter in linear fit of rms in downstream U plane channels
	TBranch *bV2p1 = metricTree->Branch("V2p1",&d2,"V2p1/D");			//1st order parameter in linear fit of rms in downstream U plane channels
	TBranch *bV2rms = metricTree->Branch("V2rms",&d1,"V2rms/D");			//RMS of residuals of fit in downstream U plane channels
	TBranch *bV2outlier = metricTree->Branch("V2outlier",&i1,"V2outlier/I");	//Number of outliers in downstream U plane channels based on 2RMS cut
	//W plane
	TBranch *bWmean = metricTree->Branch("Wmean",&d1,"Wmean/D");			//Average of RMS in W plane long wires
	TBranch *bWrms = metricTree->Branch("Wrms",&d1,"Wrms/D");			//Standard deviation of RMS in W plane long wires
	TBranch *bWoutlier = metricTree->Branch("Woutlier",&i1,"Woutlier/I");		//Number of outliers in W plane long wires with 3RMS cut
	//Misconfigured Channels
	TBranch *bmis1mean = metricTree->Branch("mis1mean",&d1,"mis1mean/D");		//Average of RMS in channels 2016-2095
	TBranch *bmis1rms = metricTree->Branch("mis1rms",&d2,"mis1rms/D");		//Standard deviation of RMS in channels 2016-2095
	TBranch *bmis2mean = metricTree->Branch("mis2mean",&d1,"mis2mean/D");		//Average of RMS in channels 2192-2303
	TBranch *bmis2rms = metricTree->Branch("mis2rms",&d2,"mis2rms/D");		//Standard deviation of RMS in channels 2192-2303
	TBranch *bmis3mean = metricTree->Branch("mis3mean",&d1,"mis3mean/D");		//Average of RMS in channels 2352-2383
	TBranch *bmis3rms = metricTree->Branch("mis3rms",&d2,"mis3rms/D");		//Standard deviation of RMS in channels 2352-2383
	//Number of dead channels
	TBranch *bdead = metricTree->Branch("numDead",&i1,"numDead/I");
		
	metricTree->SetEntries(nEvents);
		 	
	if(inputFile.is_open()){
		while(getline(inputFile, textFileName)){
			string rootFileName = textFileName.substr(0,textFileName.length()-3) + "root";
			rf = TFile::Open(rootFileName.c_str());
			rt = (TTree*)rf->Get("Trun");
			rt->Draw("eventTime");
			eventTime = rt->GetLeaf("eventTime")->GetValue();
			delete rf;
			beventTime->Fill();
			ifstream metricFile(textFileName);
			if(metricFile.is_open()){
				while (getline(metricFile, line))
    				{
        				istringstream ss(line);
					
					//Fill run stats
					if(lineCount == 0){
					ss >> i1 >> i2 >> i3;
					brun->Fill(); bsubrun->Fill(); bevent->Fill();}
					//Fill U plane
					//Upstream
					if(lineCount == 1){
					ss >> d1 >> d2;
					bU1p0->Fill(); bU1p1->Fill();}
					if(lineCount == 2){
					ss >> d1;
					bU1rms->Fill();}
					if(lineCount == 3){
					ss >> i1;
					bU1outlier->Fill();}
					//Longwire
					if(lineCount == 4){
					ss >> d1;
					bUmean->Fill();}
					if(lineCount == 5){
					ss >> d1;
					bUrms->Fill();}
					if(lineCount == 6){
					ss >> i1;
					bUoutlier->Fill();}
					//Downstream
					if(lineCount == 7){
					ss >> d1 >> d2;
					bU2p0->Fill(); bU2p1->Fill();}
					if(lineCount == 8){
					ss >> d1;
					bU2rms->Fill();}
					if(lineCount == 9){
					ss >> i1;
					bU2outlier->Fill();}
					
					//Fill V plane
					//Upstream
					if(lineCount == 10){
					ss >> d1 >> d2;
					bV1p0->Fill(); bV1p1->Fill();}
					if(lineCount == 11){
					ss >> d1;
					bV1rms->Fill();}
					if(lineCount == 12){
					ss >> i1;
					bV1outlier->Fill();}
					//Longwire
					if(lineCount == 13){
					ss >> d1;
					bVmean->Fill();}
					if(lineCount == 14){
					ss >> d1;
					bVrms->Fill();}
					if(lineCount == 15){
					ss >> i1;
					bVoutlier->Fill();}
					//Downstream
					if(lineCount == 16){
					ss >> d1 >> d2;
					bV2p0->Fill(); bV2p1->Fill();}
					if(lineCount == 17){
					ss >> d1;
					bV2rms->Fill();}
					if(lineCount == 18){
					ss >> i1;
					bV2outlier->Fill();}
					//Fill W Plan
					if(lineCount == 19){
					ss >> d1;
					bWmean->Fill();}
					if(lineCount == 20){
					ss >> d1;
					bWrms->Fill();}
					if(lineCount == 21){
					ss >> i1;
					bWoutlier->Fill();}
					//Fill misconfigured channels
					if(lineCount == 22){
					ss >> d1 >> d2;
					bmis1mean->Fill(); bmis1rms->Fill();}
					if(lineCount == 23){
					ss >> d1 >> d2;
					bmis2mean->Fill(); bmis2rms->Fill();}
					if(lineCount == 24){
					ss >> d1 >> d2;
					bmis3mean->Fill(); bmis3rms->Fill();}
					if(lineCount == 25){
					ss >> i1;
					bdead->Fill();}
					
					lineCount++;
    				}
			}
			else cerr<<"Couldn't access "<<metricFile<<endl;
			fileCount++;
			lineCount++;
			if(lineCount < 26) cerr<<"Only  "<<lineCount<<" out of 26 lines read in "<<metricFile<<endl;
			lineCount = 0;
		}
	}
	else cerr<<"Couldn't access "<<inputFile<<endl;
	if(fileCount!=nEvents) cerr<<"Mismatch: number of events and number of input root files"<<endl;
	metricRootFile.Write();
}
