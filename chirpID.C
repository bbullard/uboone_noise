#include <iostream>
#include <fstream>
#include <TTree.h>
using namespace std;

void chirpID(const char * fileName, const Int_t nEvents)
{
	TCanvas can("can","can",700,700);
	string rootFileName;
	ifstream file(fileName);
	TTree *tree;
	TFile *f = 0;
	Int_t fileCount = 0;
	vector<TTree*> eventTrees(nEvents);
	vector<TFile*> eventFiles(nEvents);
	vector<TBranch*> eventBranches(nEvents);
	vector<Bool_t> chstats(nEvents);
	Bool_t chstat = 0, chirping_flag = 0, doTheThingWinThePoints = 1;
	
	if(file.is_open()){
		while(!file.eof() && fileCount<nEvents){
			getline(file,rootFileName);
			eventFiles.at(fileCount) = TFile::Open(rootFileName.c_str(),"update");
			tree = (TTree*)eventFiles.at(fileCount)->Get("RMS");
			eventTrees.at(fileCount++) = tree;
		}
		file.close();
	}	
	if(fileCount != nEvents) cerr<<"Mismatch: number of events and number of input root files"<<endl;
	
	for(Int_t i = 0; i < nEvents; i++){
		eventBranches.at(i) = eventTrees.at(i)->Branch("chirping_flag",&chirping_flag,"chirping_flag/O");
		eventTrees.at(i)->SetBranchAddress("chstat",&chstat);
	}

	for(Int_t i = 0; i < 8256; i++){
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
		eventFiles.at(n)->Delete("RMS;1");
		eventFiles.at(n)->Write();
		delete eventFiles.at(n);
	}	
}
