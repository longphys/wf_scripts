// This is a script to pre-process waveforms from a ROOT file, 
// applying an ARC (Amplitude Ratio Correction) 
// and saving the results in a new ROOT file.

#include "TH1.h"
#include "TH2.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCutG.h"
#include "TLegend.h"
#include "TSystem.h"

#include <iostream>

void draw_wf()
{
	TFile* f = new TFile("/home/long/scripts/wf_scripts/BC404_graphs.root","read");
    TTree* t = (TTree*)f->Get("results");
    
    std::vector<double> *height = nullptr;
    std::vector<double> bin;

    t->SetBranchAddress("height", &height);

    for(int i = 0; i<2000; i++)
    {
        bin.push_back(i);
    }
    t->GetEntry(0);

    TGraph* gr = new TGraph(height->size(),&bin[0],&(height->at(0)));

    gr->SetMarkerStyle(7);
    gr->SetLineColor(kRed);
    gr->Draw("ALP");
}
