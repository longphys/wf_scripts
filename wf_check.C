// This script just checks the pre-processed waveforms

#include "TH1.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"

#include <iostream>

const int n_wf = 512;
double n_wf_height = 1000.;

void wf_check()
{
	auto timer = new TStopwatch();
	timer->Start();

    TFile* file_array = new TFile("/home/long/data/wf_files/output/wf_array.root", "read");
	TTree* tree_array = (TTree*)file_array->Get("wf_array");
    tree_array->Print();

	double wf_array_n_gamma[n_wf];
	double wf_array_gamma[n_wf];
	tree_array->SetBranchAddress("wf_n_gamma", wf_array_n_gamma);
	tree_array->SetBranchAddress("wf_gamma", wf_array_gamma);
	
	TCanvas* canvas_1 = new TCanvas("canvas_1", "canvas_1", 1400, 700);
	canvas_1->Divide(2,1);

	const int n = 20;
	std::cout << "Number of entries: " << n << "\n";
	
	// Histogram min and max
	int wf_min = 0;
	int wf_max = n_wf;

    
    for (int i = 0; i < n; i++)
    {
        TH1D* hist_wf_n_gamma = new TH1D("hist_wf_n_gamma", "hist_wf_n_gamma", n_wf, wf_min, wf_max);
        TH1D* hist_wf_gamma = new TH1D("hist_wf_gamma", "hist_wf_gamma", n_wf, wf_min, wf_max);
        tree_array->GetEntry(i);
        for (int j = 0; j < n_wf; j++)
        {
            hist_wf_n_gamma->SetBinContent(j, wf_array_n_gamma[j]);
            hist_wf_gamma->SetBinContent(j, wf_array_gamma[j]);
        }
        canvas_1->cd(1);
        hist_wf_n_gamma->SetName(Form("hist_n_gamma_%d", i));
        hist_wf_n_gamma->Draw("same");
        canvas_1->cd(2);
        hist_wf_gamma->SetName(Form("hist_gamma_%d", i));
        hist_wf_gamma->Draw("same");
        canvas_1->Modified();
        canvas_1->Update();
        gSystem->ProcessEvents();
        
        sleep(1);
    }

	std::cout << "time: " << timer->RealTime() << " seconds \n";
}
