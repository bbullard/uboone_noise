#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
using namespace std;

void getDeadChannels(const char * eventsFileName)
{
	string rootFileName;
	ifstream file(eventsFileName);
	TFile *f;
	TTree *t, *tr;
	
	if(file.is_open()){
		while(!file.eof()){
			getline(file,rootFileName);
			f = TFile::Open(rootFileName.c_str());
			t = (TTree*)f->Get("RMS");
			tr = (TTree*)f->Get("Trun");
			tr->GetEntry(0);
			cout<<"Run "<<tr->GetLeaf("runNo")->GetValue(0)<<", subrun "<<tr->GetLeaf("subRunNo")->GetValue(0)<<", event "<<tr->GetLeaf("eventNo")->GetValue(0)<<endl;
			cout<<"Dead channels (U): "<<t->Draw("","chstat==0&&chid<2401")<<endl;
			cout<<"Dead channels (V): "<<t->Draw("","chstat==0&&chid<4801&&chid>2400")<<endl;
			cout<<"Dead channels (W): "<<t->Draw("","chstat==0&&chid>4800")<<endl<<endl;
			delete f;
		}
		file.close();
	}
}
