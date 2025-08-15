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

int ApplyARC(TH1D* hist, double arc_k, int arc_tau_d) 
{
    // Returns trigger bin index based on ARC maximum
    int nbins = hist->GetNbinsX();
    double max_arc = -1e9;
    int trigger_index = 1;

    for (int i = arc_tau_d+1; i <= nbins; ++i)
    {
        double v1 = hist->GetBinContent(i - arc_tau_d);
        double v2 = hist->GetBinContent(i);
        double arc = v1 - arc_k * v2;
        if (arc > max_arc)
        {
            max_arc = arc;
            trigger_index = i;
        }
	}
    return trigger_index;  // bin index of aligned pulse
}

const int n_wf = 512;
double n_wf_height = 1000.;

void wf_pre_analyse_1()
{

	auto timer = new TStopwatch();
	timer->Start();

	//Input local files
	//TFile* file_n_gamma = new TFile("/mnt/c/Users/Long/Desktop/data/wf_files/input/stilbene_neutrons.root", "read");
	//TFile* file_gamma = new TFile("/mnt/c/Users/Long/Desktop/data/wf_files/input/stilbene_cs137.root", "read");

	TFile* file_n_gamma = new TFile("~/data/wf_files/input/stilbene_neutrons.root", "read");
	TFile* file_gamma = new TFile("~/data/wf_files/input/stilbene_cs137.root", "read");
	
	//Input files from link
	//TFile* file_n_gamma = TFile::Open("https://zenodo.org/records/16795081/files/stilbene_neutrons.root?download=1");
	//TFile* file_gamma = TFile::Open("https://zenodo.org/records/16795081/files/stilbene_cs137.root?download=1");

	TTree* tree_n_gamma = (TTree*)file_n_gamma->Get("tt");
	TTree* tree_gamma = (TTree*)file_gamma->Get("tt");
	//tree->Print();
	
	UShort_t time_n_gamma[n_wf];
	UShort_t time_gamma[n_wf];
	tree_n_gamma->SetBranchAddress("wf", time_n_gamma);
	tree_gamma->SetBranchAddress("wf", time_gamma);
	
	TCanvas* canvas_1 = new TCanvas("canvas_1", "canvas_1", 1400, 700);
	canvas_1->Divide(2,1);
	
	// TCanvas* canvas_2 = new TCanvas("canvas_2", "canvas_2", 1400, 700);
	// canvas_2->Divide(2,1);
	
	TCanvas* canvas_3 = new TCanvas("canvas_3", "canvas_3", 1500, 450);
	canvas_3->Divide(3,1);
	//const int n = tree_gamma->GetEntries();
	//const int n = tree_n_gamma->GetEntries();
	const int n = 10000;
	// const int n = 10;
	
	std::cout << "Number of entries: " << n << "\n";
	
	// Histogram min and max
	int wf_min = 0;
	int wf_max = n_wf;
	
	// Viewing window min and max
	//int wf_min_view = 50;
	//int wf_max_view = 200;
	int wf_min_view = wf_min;
	int wf_max_view = wf_max;
	
	int wf_charge_total_min = 90;
	int wf_charge_total_max = 130;
	
	int wf_charge_tail_min = 115;
	
	/*
	TCutG *cutg = new TCutG("mycut", 0);
	
	//! neutron anomaly
	cutg->SetPoint(0, 2500.0, 0.56);
	cutg->SetPoint(1, 26000.0, 0.36);
	cutg->SetPoint(2, 20500.0, 0.3);
	cutg->SetPoint(3, 2210.0, 0.455);
	cutg->SetPoint(4, 1950.0, 0.51);
	cutg->SetPoint(5, 2500.0, 0.56);
  	*/
  	
  	/*
  	std::ifstream infile_cut_1_cut_1("cut_1.txt");
    if (!infile_cut_1.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
    }

    std::vector<double> vx_cut_1, vy_cut_2;
    double x_cut_1, y_cut_1;

    while (infile_cut_1 >> x_cut_1 >> y_cut_1) {
        vx_cut_1.push_back(x_cut_1);
        vy_cut_2.push_back(y_cut_1);
    }

    infile_cut_1.close();
    */
    
    std::vector<Double_t> vx_cut_1{
      578.6806737851293, 626.220044572261, 3383.503550225918, 11275.03910088983, 28721.98817976728, 46121.3978878576, 56294.82323630386, 55724.35078685827, 45075.53173054069, 30433.40552810403,
      19927.20458414786, 11132.42098852844, 5998.16894351818, 4239.212224394294, 2955.649213141731, 2004.861797399089, 1529.468089527769, 1386.849977166374, 1054.074381656449, 673.7594153593936,
      673.7594153593936, 578.6806737851293
   };
   std::vector<Double_t> vy_cut_1{
      0.5250000018626453, 0.2398058228556391, 0.107524267486432, 0.06140776194487352, 0.05655339294049888, 0.1063106752353383, 0.1718446567943951, 0.2313106770979836, 0.1572815497812713, 0.1148058209929937,
      0.1099514519886192, 0.1354368892615857, 0.1742718412965824, 0.209466016578298, 0.2847087361461039, 0.3490291254540669, 0.4084951457576555, 0.4861650498276487, 0.5274271863648325, 0.5298543708670198,
      0.5298543708670198, 0.5250000018626453
   };

    int n_cut_1 = vx_cut_1.size();
    if (n_cut_1 < 3) {
        std::cerr << "Error: Not enough points to form a cut!" << std::endl;
    }

    TCutG* cut_1 = new TCutG("cut_1", n_cut_1);
    for (int i = 0; i < n_cut_1; ++i) {
        cut_1->SetPoint(i, vx_cut_1[i], vy_cut_1[i]);
    }
    cut_1->SetLineColor(kRed);
	
	TH1D* hist_spectrum_n_gamma = new TH1D("spectrum_n_gamma", "spectrum", 500, 0., 60000.);
	hist_spectrum_n_gamma->GetXaxis()->SetTitle("Amplitude (Channels)");
	hist_spectrum_n_gamma->GetYaxis()->SetTitle("Count/Channel");
	hist_spectrum_n_gamma->GetXaxis()->CenterTitle();
	hist_spectrum_n_gamma->GetYaxis()->CenterTitle();
	
	TH1D* hist_spectrum_gamma = new TH1D("spectrum_gamma", "spectrum", 500, 0., 60000.);
	hist_spectrum_gamma->GetXaxis()->SetTitle("Amplitude (Channels)");
	hist_spectrum_gamma->GetYaxis()->SetTitle("Count/Channel");
	hist_spectrum_gamma->GetXaxis()->CenterTitle();
	hist_spectrum_gamma->GetYaxis()->CenterTitle();
	
	TH2D* hist_Q_ratio_n_gamma = new TH2D("Q-ratio Map n-gamma", "Q-ratio Map n-gamma", 500, 0., 60000., 800, 0., 0.8);
	hist_Q_ratio_n_gamma->GetXaxis()->SetTitle("Charge (a. unit)");
	hist_Q_ratio_n_gamma->GetYaxis()->SetTitle("Q-ratio");
	hist_Q_ratio_n_gamma->GetXaxis()->CenterTitle();
	hist_Q_ratio_n_gamma->GetYaxis()->CenterTitle();
	
	TH2D* hist_Q_ratio_gamma = new TH2D("Q-ratio Map gamma", "Q-ratio Map gamma", 500, 0., 60000., 800, 0., 0.8);
	hist_Q_ratio_gamma->GetXaxis()->SetTitle("Charge (a. unit)");
	hist_Q_ratio_gamma->GetYaxis()->SetTitle("Q-ratio");
	hist_Q_ratio_gamma->GetXaxis()->CenterTitle();
	hist_Q_ratio_gamma->GetYaxis()->CenterTitle();
	
	TGraph* graph_Q_ratio_n_gamma = new TGraph();
	
	TGraph* graph_Q_ratio_gamma = new TGraph();
	
	graph_Q_ratio_n_gamma->SetMarkerColor(kRed);
	graph_Q_ratio_gamma->SetMarkerColor(kBlue);
	graph_Q_ratio_n_gamma->SetMarkerStyle(1);
	graph_Q_ratio_gamma->SetMarkerStyle(1);
	
	TMultiGraph* graph_both = new TMultiGraph();
	graph_both->GetXaxis()->SetTitle("Charge (a. unit)");
	graph_both->GetYaxis()->SetTitle("Q-ratio");
	graph_both->GetXaxis()->CenterTitle();
	graph_both->GetYaxis()->CenterTitle();
	graph_both->Add(graph_Q_ratio_n_gamma, "AP");
	graph_both->Add(graph_Q_ratio_gamma, "AP");
	
	//TFile* file_save_array = new TFile("/mnt/c/Users/Long/Desktop/data/wf_files/output/wf_array.root", "recreate");
	TFile* file_save_array = new TFile("~/data/wf_files/output/wf_array.root", "recreate");
	
	TTree* tree_wf_array = new TTree("wf_array", "Tree of n_wf-element array");

	double wf_array_n_gamma[n_wf];
	double wf_array_gamma[n_wf];
	tree_wf_array->Branch("wf_n_gamma", wf_array_n_gamma, Form("wf_n_gamma[%d]/D", n_wf));
	tree_wf_array->Branch("wf_gamma", wf_array_gamma, Form("wf_gamma[%d]/D", n_wf));
	
	TString name;
	for (int i = 0; i < n; i++)
	{	
		if (i%10000==0)
		// if (i%50000==0)
		{
			std::cout << "Entry No. " << i << "\n";
			
			canvas_1->cd(1)->Modified();
			canvas_1->cd(1)->Update();
			canvas_1->cd(2)->Modified();
			canvas_1->cd(2)->Update();
			
			canvas_3->cd(1)->Modified();
			canvas_3->cd(1)->Update();
			canvas_3->cd(2)->Modified();
			canvas_3->cd(2)->Update();
			canvas_3->cd(3)->Modified();
			canvas_3->cd(3)->Update();
			
			gSystem->ProcessEvents();
		}
    	name = Form("hist_temp_n_gamma_%d",i);
		TH1D* hist_temp_n_gamma = new TH1D(name, name, wf_max-wf_min, wf_min, wf_max);
		tree_n_gamma->GetEntry(i);
		
    	name = Form("hist_temp_gamma_%d",i);
		TH1D* hist_temp_gamma = new TH1D(name, name, wf_max-wf_min, wf_min, wf_max);
		tree_gamma->GetEntry(i);
		
		// Get and fill waveforms
		for (int j = 0; j < wf_max-wf_min; j++)
		{
			double inverse_height = 4096-time_n_gamma[j+wf_min];
			hist_temp_n_gamma->SetBinContent(j+1, inverse_height);
			//hist_temp_n_gamma->SetBinContent(j+1, time_n_gamma[j+wf_min]);
		}
		
		for (int j = 0; j < wf_max-wf_min; j++)
		{
			double inverse_height = 4096-time_gamma[j+wf_min];
			hist_temp_gamma->SetBinContent(j+1, inverse_height);
			//hist_temp_gamma->SetBinContent(j+1, time_gamma[j+wf_min]);
		}
		
		int maximum_n_gamma = hist_temp_n_gamma->GetMaximum();
		double scale_factor_n_gamma = n_wf_height/(maximum_n_gamma);
		hist_temp_n_gamma->Scale(scale_factor_n_gamma, "noSW2");
		
		int maximum_gamma = hist_temp_gamma->GetMaximum();
		double scale_factor_gamma = n_wf_height/(maximum_gamma);
		hist_temp_gamma->Scale(scale_factor_gamma, "noSW2");
		
		// Apply ARC
		// std::cout << "ARC_bin_index = " << ApplyARC(hist_temp, 0.95, 2) << "\n";
		
		int arc_trigger_bin_n_gamma = ApplyARC(hist_temp_n_gamma, 0.95, 2);
		int arc_trigger_bin_gamma = ApplyARC(hist_temp_gamma, 0.95, 2);
		int arc_bin_align_to = 100;
		int shift_n_gamma = arc_bin_align_to - arc_trigger_bin_n_gamma;
		int shift_gamma = arc_bin_align_to - arc_trigger_bin_gamma;
		
    	name = Form("hist_temp_aligned_n_gamma_%d",i); 
		TH1D* hist_temp_aligned_n_gamma = 
		new TH1D(
				name, 
				name, 
				hist_temp_n_gamma->GetNbinsX(), 
				hist_temp_n_gamma->GetXaxis()->GetXmin(), 
				hist_temp_n_gamma->GetXaxis()->GetXmax());
				
    	name = Form("hist_temp_aligned_gamma_%d",i); 
		TH1D* hist_temp_aligned_gamma = 
		new TH1D(
				name, 
				name, 
				hist_temp_gamma->GetNbinsX(), 
				hist_temp_gamma->GetXaxis()->GetXmin(), 
				hist_temp_gamma->GetXaxis()->GetXmax());
				
		for(int j = 0; j < hist_temp_aligned_n_gamma->GetNbinsX(); j++)
		{
			int newbin = j + shift_n_gamma;
			if (newbin >= 1 && newbin <= hist_temp_aligned_n_gamma->GetNbinsX())
			{
        		hist_temp_aligned_n_gamma->SetBinContent(newbin, hist_temp_n_gamma->GetBinContent(j));
			}
		}
		
		for(int j = 0; j < hist_temp_aligned_gamma->GetNbinsX(); j++)
		{
			int newbin = j + shift_gamma;
			if (newbin >= 1 && newbin <= hist_temp_aligned_gamma->GetNbinsX())
			{
        		hist_temp_aligned_gamma->SetBinContent(newbin, hist_temp_gamma->GetBinContent(j));
			}
		}
		
		// Save waveforms as arrays of 512 elements into a file
		for(int j = 0; j < wf_max; j++)
		{	
			wf_array_n_gamma[j] = hist_temp_aligned_n_gamma->GetBinContent(j);
			wf_array_gamma[j] = hist_temp_aligned_gamma->GetBinContent(j);
			//std::cout << "wf_array_n_gamma[" << j << "] = " << wf_array_n_gamma[j] << "\n";
			//std::cout << "wf_array_gamma[" << j << "] = " << wf_array_gamma[j] << "\n";
		}
		tree_wf_array->Fill();
		
		// Calculation of the charge
		hist_temp_aligned_n_gamma->Scale(1/scale_factor_n_gamma, "noSW2");
		hist_temp_aligned_gamma->Scale(1/scale_factor_gamma, "noSW2");
		double charge_total_n_gamma = 0.;
		double charge_tail_n_gamma = 0.;
		double charge_total_gamma = 0.;
		double charge_tail_gamma = 0.;
		
		for(int j = wf_charge_total_min; j < wf_charge_total_max; j++)
		{
			charge_total_n_gamma += hist_temp_aligned_n_gamma->GetBinContent(j+1);
			charge_total_gamma += hist_temp_aligned_gamma->GetBinContent(j+1);
			if(j < wf_charge_tail_min)
			{}
			else
			{
				charge_tail_n_gamma += hist_temp_aligned_n_gamma->GetBinContent(j+1);
				charge_tail_gamma += hist_temp_aligned_gamma->GetBinContent(j+1);
			}
		}
		double q_ratio_n_gamma = charge_tail_n_gamma/charge_total_n_gamma;
		double q_ratio_gamma = charge_tail_gamma/charge_total_gamma;
		
		//std::cout << "charge total = " << charge_total << "\n";
		//std::cout << "charge tail/charge total = " << q_ratio << "\n";
			
		// Analyze waveforms and fill spectrum
		hist_Q_ratio_n_gamma->Fill(charge_total_n_gamma, q_ratio_n_gamma);
		graph_Q_ratio_n_gamma->AddPoint(charge_total_n_gamma, q_ratio_n_gamma);
        if (cut_1->IsInside(charge_total_n_gamma, q_ratio_n_gamma))
        { //! if (x = energy, y = timediff is inside cutg) return 1
			hist_spectrum_n_gamma->Fill(charge_total_n_gamma);
        }
        else
        {}
        
        if (cut_1->IsInside(charge_total_gamma, q_ratio_gamma))
        { //! if (x = energy, y = timediff is inside cutg) return 1
		    if (i < 280000)
		    {
				hist_spectrum_gamma->Fill(charge_total_gamma);
			}
        }
        //if (i < tree_gamma->GetEntriesFast())
        //if (i < tree_gamma->GetEntriesFast() & i > 300000)
        //if (i < 300000 & i > 250000)
        if (i < 280000)
        {
			//hist_spectrum_gamma->Fill(charge_total_gamma);
			hist_Q_ratio_gamma->Fill(charge_total_gamma, q_ratio_gamma);
			graph_Q_ratio_gamma->AddPoint(charge_total_gamma, q_ratio_gamma);
		}
		
		// Draw Histograms
		if (i == 0)
		{
			canvas_1->cd(1);
			hist_spectrum_n_gamma->Draw();
		}
		if (i == 0)
		{
			canvas_1->cd(2);
			hist_spectrum_gamma->Draw();
		}
		
		// // Clone for hist_temp or hist_temp_aligned;
		// //TH1D* hist_temp_clone = (TH1D*)hist_temp->Clone();
		// TH1D* hist_temp_clone_n_gamma = (TH1D*)hist_temp_aligned_n_gamma->Clone();
    	// name = Form("hist_temp_clone_n_gamma_%d",i);
    	// hist_temp_clone_n_gamma->SetName(name);
    	// hist_temp_clone_n_gamma->SetTitle(name);
		// hist_temp_clone_n_gamma->SetDirectory(0);
		// hist_temp_clone_n_gamma->GetXaxis()->SetRangeUser(wf_min_view, wf_max_view);
		// hist_temp_clone_n_gamma->GetYaxis()->SetRangeUser(0, 1100);
		
		// TH1D* hist_temp_clone_gamma = (TH1D*)hist_temp_aligned_gamma->Clone();
    	// name = Form ("hist_temp_clone_gamma_%d",i);
    	// hist_temp_clone_gamma->SetName(name);
    	// hist_temp_clone_gamma->SetTitle(name);
		// hist_temp_clone_gamma->SetDirectory(0);
		// hist_temp_clone_gamma->GetXaxis()->SetRangeUser(wf_min_view, wf_max_view);
		// hist_temp_clone_gamma->GetYaxis()->SetRangeUser(0, 1100);
		
		// canvas_2->cd(1);
		
		// //hist_temp->Draw();
		// sleep(1);
		// int colorIndex = i % 50 + 1;  // ROOT has colors 1–50 (looping)
		
        // hist_temp_clone_n_gamma->SetLineColorAlpha(colorIndex, 0.05);  // Faint line
        // hist_temp_clone_n_gamma->SetLineWidth(1);
		// hist_temp_clone_n_gamma->Draw("same");
		
		
		// canvas_2->cd(2);
        // hist_temp_clone_gamma->SetLineColorAlpha(colorIndex, 0.05);  // Faint line
        // hist_temp_clone_gamma->SetLineWidth(1);
		// hist_temp_clone_gamma->Draw("same");
		
		// canvas_2->Modified();
		// canvas_2->Update();
		
		hist_temp_n_gamma->SetDirectory(0);
		hist_temp_gamma->SetDirectory(0);
		hist_temp_aligned_n_gamma->SetDirectory(0);
		hist_temp_aligned_gamma->SetDirectory(0);
		
		delete hist_temp_n_gamma;
		delete hist_temp_gamma;
		delete hist_temp_aligned_n_gamma;
		delete hist_temp_aligned_gamma;
		
		if (i == 0)
		{	
			canvas_3->cd(1);
			hist_Q_ratio_n_gamma->Draw();
			cut_1->Draw("same");
			canvas_3->cd(2);
			hist_Q_ratio_gamma->Draw();
			cut_1->Draw("same");
			canvas_3->cd(3);
			graph_both->Draw("a");
			
			TLegend *legend_1 = new TLegend(0.75, 0.6, 0.98, 0.75);
			//legend_1->SetHeader("test", "C");
			legend_1->SetBorderSize(2);
			legend_1->AddEntry(hist_Q_ratio_n_gamma, "some neutron source", "p");
			legend_1->AddEntry(hist_Q_ratio_gamma, "Cs137", "p");
			legend_1->Draw();
		}
		graph_both->GetXaxis()->SetRangeUser(
				hist_Q_ratio_n_gamma->GetXaxis()->GetXmin(), 
				hist_Q_ratio_n_gamma->GetXaxis()->GetXmax());
		graph_both->GetYaxis()->SetRangeUser(
				hist_Q_ratio_n_gamma->GetYaxis()->GetXmin(), 
				hist_Q_ratio_n_gamma->GetYaxis()->GetXmax());
	}
	
	tree_wf_array->Write();
	file_save_array->Close();
	
	//TFile* file_save_spectrum = new TFile("/mnt/c/Users/Long/Desktop/data/wf_files/output/wf_out.root" , "recreate");
	TFile* file_save_spectrum = new TFile("~/data/wf_files/output/wf_out.root" , "recreate");
	hist_spectrum_n_gamma->Write("spectrum_n_gamma");
	hist_spectrum_gamma->Write("spectrum_gamma");
	
	hist_Q_ratio_n_gamma->Write("q_ratio_n_gamma");
	hist_Q_ratio_gamma->Write("q_ratio_gamma");
	
	graph_both->Write("q_ratio_both");
	
	file_save_spectrum->Write();
	file_save_spectrum->Close();
	
	std::cout << "time: " << timer->RealTime() << " seconds \n";
}
