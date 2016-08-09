#include <TH1.h>
#include <TTree.h>
using namespace std;

void getStats(TH1F * h, Double_t &mlow, Double_t &mmid, Double_t &mhigh, Double_t &mco, Double_t &m36, Double_t &m108, Double_t &m900)
{
	mlow = 0; mmid = 0; mhigh = 0; mco = 0; m36 = 0; m108 = 0; m900 = 0;
	for(Int_t bin = 200; bin<500;bin++)
		mlow += h->GetBinContent(bin)/300.;
	for(Int_t bin = 500; bin<1000;bin++)
		mmid += h->GetBinContent(bin)/500.;
	for(Int_t bin = 3200; bin<4000;bin++)
		mhigh += h->GetBinContent(bin)/800.;
	for(Int_t bin = 48; bin<145;bin++)
		mco += h->GetBinContent(bin)/97.;
	m36 = h->GetBinContent(172);
	m108 = h->GetBinContent(518);
	for(Int_t bin = 4120; bin<4521;bin++)
		m900 += h->GetBinContent(bin)/400.;
}

void makeFFTMetrics(const char * eventsFileName, const Int_t nEvents)
{
	TCanvas can("can","can",700,700);
	TFile FFTFile("FFTMetrics.root","RECREATE","Trees Containing FFT Metrics");
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
	Double_t run, subrun, event;
	runStats->Branch("run",&run,"run/I");
	runStats->Branch("subrun",&subrun,"subrun/I");
	runStats->Branch("event",&event,"event/I");
	
	string rootFileName;
	ifstream file(eventsFileName);
	Int_t fileCount = 0;
	
	if(file.is_open()){
		while(!file.eof() && fileCount<nEvents){
			getline(file,rootFileName);
			f = TFile::Open(rootFileName.c_str());
			getStats((TH1F*)f->Get("UFFT_mag"),mlow[0],mmid[0],mhigh[0],mco[0],m36[0],m108[0],m900[0]);
			trees.at(0)->Fill();
			getStats((TH1F*)f->Get("UlongwireFFT_mag"),mlow[1],mmid[1],mhigh[1],mco[1],m36[1],m108[1],m900[1]);
			trees.at(1)->Fill();
			getStats((TH1F*)f->Get("UFFT_mag_misconf"),mlow[2],mmid[2],mhigh[2],mco[2],m36[2],m108[2],m900[2]);
			trees.at(2)->Fill();
			getStats((TH1F*)f->Get("VFFT_mag"),mlow[3],mmid[3],mhigh[3],mco[3],m36[3],m108[3],m900[3]);
			trees.at(3)->Fill();
			getStats((TH1F*)f->Get("VlongwireFFT_mag"),mlow[4],mmid[4],mhigh[4],mco[4],m36[4],m108[4],m900[4]);
			trees.at(4)->Fill();
			getStats((TH1F*)f->Get("WFFT_mag"),mlow[5],mmid[5],mhigh[5],mco[5],m36[5],m108[5],m900[5]);
			trees.at(5)->Fill();
			
			RUN = (TTree*)f->Get("Trun");
			//RUN->Draw("runNo");RUN->Draw("subRunNo");RUN->Draw("eventNo");
			run = RUN->GetLeaf("runNo")->GetValue();
			subrun = RUN->GetLeaf("subRunNo")->GetValue();
			event = RUN->GetLeaf("eventNo")->GetValue();
			runStats->Fill();
			can.Close();
	
			fileCount++;
		}
		file.close();
	}
	if(fileCount != nEvents) cerr<<"Mismatch in number of files read"<<endl;
	
	FFTFile.Write();
}
