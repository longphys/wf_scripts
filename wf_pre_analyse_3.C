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

// #define ORIGINAL_WF
// #define DRAW_EACH_WF
#define APPLY_FIXED_AVERAGE
// #define APPLY_MOVING_AVERAGE
// #define NORMALIZE_WF
#define HEIGHT_THRESHOLD
double height_threshold_n_gamma = 10.;
double height_threshold_gamma = 10.;
// #define USE_SIGNAL_HEIGHT // Use the pulse height to draw the 2D histograms instead of the integral.

//! Short window
const int n_wf = 2000;
//! Long window
// const int n_wf = 2500;
double n_wf_height = 1000.;
double view_hist_height = 200.;

// Histogram min and max
int wf_min = 0;
int wf_max = n_wf;

// Viewing window min and max
//int wf_min_view = 50;
//int wf_max_view = 200;
int wf_min_view = wf_min;
#ifdef APPLY_FIXED_AVERAGE
int wf_max_view = (int)wf_max/4.;
#endif
#ifdef APPLY_MOVING_AVERAGE
int wf_max_view = wf_max;
#endif


void wf_pre_analyse_3()
{
	auto timer = new TStopwatch();
	timer->Start();

	//! Input local files

	//! Short window
	TFile* file_n_gamma = new TFile("~/data/wf_files/input/root_files/bc404_pu_c13_20250723-0002.root", "read");
	
	//! Long window
	// TFile* file_n_gamma = new TFile("~/data/wf_files/input/root_files/bc404_pu_c13_200ps.root", "read");

	//! Short window
	TFile* file_gamma = new TFile("~/data/wf_files/input/root_files/bc404_na22_20250723.root", "read");

	//! Long window
	// TFile* file_gamma = new TFile("~/data/wf_files/input/root_files/bc404_na22_20250725-0005.root", "read");
	// TFile* file_gamma = new TFile("~/data/wf_files/input/root_files/bc404_cs137_20250725-0006.root", "read");
	
	
	//! Input files from link
	//TFile* file_n_gamma = TFile::Open("https://zenodo.org/records/16795081/files/stilbene_neutrons.root?download=1");
	//TFile* file_gamma = TFile::Open("https://zenodo.org/records/16795081/files/stilbene_cs137.root?download=1");
	
	// TFile* file_n_gamma = TFile::Open("https://zenodo.org/records/16906274/files/bc404_pu_c13_200ps.root?download=1");
	// TFile* file_n_gamma = TFile::Open("https://zenodo.org/records/16906274/files/bc404_pu_c13_1.root?download=1");
	// TFile* file_gamma = TFile::Open("https://zenodo.org/records/16906274/files/bc404_na22.root?download=1");

	TTree* tree_n_gamma = (TTree*)file_n_gamma->Get("wf");
	TTree* tree_gamma = (TTree*)file_gamma->Get("wf");
	//tree->Print();
	
	std::vector<float>* time_n_gamma = nullptr;
	std::vector<float>* time_gamma   = nullptr;
	std::vector<float>* voltage_n_gamma = nullptr;
	std::vector<float>* voltage_gamma   = nullptr;
	// UShort_t time_n_gamma[n_wf];
	// UShort_t time_gamma[n_wf];
	tree_n_gamma->SetBranchAddress("Time", &time_n_gamma);
	tree_gamma->SetBranchAddress("Time", &time_gamma);
	tree_n_gamma->SetBranchAddress("Voltage", &voltage_n_gamma);
	tree_gamma->SetBranchAddress("Voltage", &voltage_gamma);
	
	TCanvas* canvas_2 = new TCanvas("canvas_2", "canvas_2", 1400, 700);
	canvas_2->Divide(2,1);
	
	//const int n = tree_gamma->GetEntries();
	//const int n = tree_n_gamma->GetEntries();

	//! Short window
	const int n = 28000;
	//! Long window
	// const int n = 50000;
	// const int n = 100;
	
	std::cout << "Number of entries: " << n << "\n";

	TFile* outFile = new TFile("BC404_wf.root","recreate");
	TTree* outTree = new TTree("tree","tree");

	std::vector<double> h_ave_n_gamma;
	std::vector<double> h_ave_gamma;

	double integral_total_n_gamma;
	double integral_total_gamma;
	std::vector<double> integral_tail_n_gamma;
	std::vector<double> integral_tail_gamma;
	
	double min_total;
	std::vector<double> min_tail;
	double max;
	
	outTree->Branch("t_n_gamma", &time_n_gamma); // Time values of a waveform
	outTree->Branch("t_gamma", &time_gamma);
	outTree->Branch("v_n_gamma", &voltage_n_gamma); // Voltage values of a waveform
	outTree->Branch("v_gamma", &voltage_gamma);
	outTree->Branch("h_ave_n_gamma", &h_ave_n_gamma); // Averaged heights of a waveform (voltage)
	outTree->Branch("h_ave_gamma", &h_ave_gamma);
	outTree->Branch("integral_total_n_gamma", &integral_total_n_gamma); // Total/tail integral of a waveform, for different ranges.
	outTree->Branch("integral_total_gamma", &integral_total_gamma);
	outTree->Branch("integral_tail_n_gamma", &integral_tail_n_gamma);
	outTree->Branch("integral_tail_gamma", &integral_tail_gamma);
	outTree->Branch("min_total", &min_total); // Total/tail ranges to take integral
	outTree->Branch("min_tail", &min_tail); 
	outTree->Branch("max", &max); 
	
	//! Short window
	int wf_charge_total_min = 145;
	int wf_charge_total_max = 500;
	int wf_charge_tail_min = 200;
	int wf_charge_tail_min_final = 450;
	double baseline = 188.;
	
	//! Long window
	// int wf_charge_total_min = 1200;
	// int wf_charge_total_max = 2450;
	// int wf_charge_tail_min = 1250;
	// int wf_charge_tail_min_final = 2400;
	// double baseline = 185.;

	int wf_charge_tail_min_incre = 1;
	// int wf_charge_tail_min_incre = 50;
	int times = (wf_charge_tail_min_final-wf_charge_tail_min)/wf_charge_tail_min_incre;
	// std::cout << "times = " << times << "\n";

	TString name;
	for (int i = 0; i < n; i++)
	// for (int i = 9705; i < 9706; i++)
	// for (int i = 1000; i < 1001; i++)
	// for (int i = 8781; i < 8782; i++)
	{	
		if (i%1000==0)
		// if (i%50000==0)
		{
			std::cout << "Entry No. " << i << "\n";
		}
    	name = Form("hist_temp_n_gamma_%d",i);
		TH1D* hist_temp_n_gamma = new TH1D(name, name, wf_max-wf_min, wf_min, wf_max);
		tree_n_gamma->GetEntry(i);
		
    	name = Form("hist_temp_gamma_%d",i);
		TH1D* hist_temp_gamma = new TH1D(name, name, wf_max-wf_min, wf_min, wf_max);
		// TH1D* hist_temp_gamma = new TH1D(name, name, wf_max-wf_min, wf_min+1500, wf_max+1500);
		tree_gamma->GetEntry(i);

		hist_temp_n_gamma->SetDirectory(0);
		hist_temp_gamma->SetDirectory(0);
		
		// Get and fill waveforms
		double inverse_height_a;
		double inverse_height_b;

		#ifdef ORIGINAL_WF
		for (int j = 0; j < wf_max-wf_min; j++)
		{
			inverse_height_a = (*voltage_n_gamma)[j+wf_min]; //! Original waveform
			inverse_height_b = (*voltage_gamma)[j+wf_min];
			hist_temp_n_gamma->SetBinContent(j+1, inverse_height_a);
			hist_temp_gamma->SetBinContent(j+1, inverse_height_b);
		}
		#else
		for (int j = 0; j < wf_max-wf_min; j++)
		{
			double inverse_height_a = baseline-(*voltage_n_gamma)[j+wf_min];
			// double inverse_height_a = 183.-(*voltage_n_gamma)[j+wf_min];
			
			if(inverse_height_a < 0)
			{
				hist_temp_n_gamma->SetBinContent(j+1, 1.);
			}
			else
			{
				hist_temp_n_gamma->SetBinContent(j+1, inverse_height_a);
			}

			double inverse_height_b = baseline-(*voltage_gamma)[j+wf_min];
			// double inverse_height_b = 183.-(*voltage_gamma)[j+wf_min];
			// std::cout << "inverse_height_b = " << inverse_height_b << "\n";
			if(inverse_height_b < 0)
			{
				hist_temp_gamma->SetBinContent(j+1, 1.);
			}
			else
			{
				hist_temp_gamma->SetBinContent(j+1, inverse_height_b);
			}
		}
		#endif

		int maximum_n_gamma = hist_temp_n_gamma->GetMaximum();
		int maximum_gamma = hist_temp_gamma->GetMaximum();

		// Apply average
		#ifdef APPLY_FIXED_AVERAGE
		for (int j = 0; j < wf_max-wf_min; j+=4) // fixed width (4 bins) average
		{
			double average_a = 
			(
				hist_temp_n_gamma->GetBinContent(j)+
				hist_temp_n_gamma->GetBinContent(j+1)+
				hist_temp_n_gamma->GetBinContent(j+2)+
				hist_temp_n_gamma->GetBinContent(j+3)
			)/4.0;
			hist_temp_n_gamma->SetBinContent(j, average_a);
			hist_temp_n_gamma->SetBinContent(j+1, average_a);
			hist_temp_n_gamma->SetBinContent(j+2, average_a);
			hist_temp_n_gamma->SetBinContent(j+3, average_a);
			
			// if (j<1996)
			// {
			// 	for(int t = 0; t <4; t++)
			// 	{
					h_ave_n_gamma.push_back(average_a);
			// 	}
			// }

			double average_b = 
			(
				hist_temp_gamma->GetBinContent(j)+
				hist_temp_gamma->GetBinContent(j+1)+
				hist_temp_gamma->GetBinContent(j+2)+
				hist_temp_gamma->GetBinContent(j+3)
			)/4.0;
			hist_temp_gamma->SetBinContent(j, average_b);
			hist_temp_gamma->SetBinContent(j+1, average_b);
			hist_temp_gamma->SetBinContent(j+2, average_b);
			hist_temp_gamma->SetBinContent(j+3, average_b);
			
			// if (j<1996)
			// {
			// 	for(int t = 0; t <4; t++)
			// 	{
					h_ave_gamma.push_back(average_b);
			// 	}
			// }
		}

		hist_temp_n_gamma->Clear();
		hist_temp_gamma->Clear();
		
		hist_temp_n_gamma->SetBins(h_ave_gamma.size(),0.,(double)h_ave_gamma.size());
		hist_temp_gamma->SetBins(h_ave_gamma.size(),0.,(double)h_ave_gamma.size());
		
		for(int new_i = 0; new_i < h_ave_gamma.size(); new_i++)
		{
			hist_temp_n_gamma->SetBinContent(new_i, h_ave_n_gamma.at(new_i));
			hist_temp_gamma->SetBinContent(new_i, h_ave_gamma.at(new_i));
		}

		#endif
		
		#ifdef APPLY_MOVING_AVERAGE
		for (int j = 0; j < wf_max-wf_min; j++) // moving average
		{
			double average_a = 
			(
				hist_temp_n_gamma->GetBinContent(j)+
				hist_temp_n_gamma->GetBinContent(j+1)+
				hist_temp_n_gamma->GetBinContent(j+2)+
				hist_temp_n_gamma->GetBinContent(j+3)
			)/4.0;
			hist_temp_n_gamma->SetBinContent(j, average_a);
			
			h_ave_n_gamma.push_back(average_a);

			double average_b = 
			(
				hist_temp_gamma->GetBinContent(j)+
				hist_temp_gamma->GetBinContent(j+1)+
				hist_temp_gamma->GetBinContent(j+2)+
				hist_temp_gamma->GetBinContent(j+3)
			)/4.0;
			hist_temp_gamma->SetBinContent(j, average_b);
			
			h_ave_gamma.push_back(average_b);
		}
		#endif
		
		#ifdef NORMALIZE_WF
		double scale_factor_n_gamma = n_wf_height/(maximum_n_gamma);
		hist_temp_n_gamma->Scale(scale_factor_n_gamma, "noSW2");
		
		double scale_factor_gamma = n_wf_height/(maximum_gamma);
		hist_temp_gamma->Scale(scale_factor_gamma, "noSW2");
		#endif

		double charge_total_n_gamma = 0.;
		double charge_tail_n_gamma = 0.;
		double charge_total_gamma = 0.;
		double charge_tail_gamma = 0.;
		
		double q_ratio_n_gamma;
		double q_ratio_gamma;
		
		int safe_bins = hist_temp_n_gamma->GetNbinsX();

		// fixed bins intergration
		for(int j = wf_charge_total_min; j < wf_charge_total_max; j++)
		{
			if(j > safe_bins)
			{break;}
			else
			{
				charge_total_n_gamma = charge_total_n_gamma + hist_temp_n_gamma->GetBinContent(j);
				charge_total_gamma = charge_total_gamma + hist_temp_gamma->GetBinContent(j);
			}
		}
		integral_total_n_gamma = charge_total_n_gamma;
		integral_total_gamma = charge_total_gamma;
		
		min_total = wf_charge_total_min;
		max = wf_charge_total_max;

		int temp = wf_charge_tail_min;
		for (int iii = 0; iii < times; iii++)
		{
			wf_charge_tail_min += wf_charge_tail_min_incre;
			for(int j = wf_charge_tail_min; j < wf_charge_total_max; j++)
			{
				if(j > safe_bins)
				{break;}
				else
				{							
					charge_tail_n_gamma = charge_tail_n_gamma + hist_temp_n_gamma->GetBinContent(j);
					charge_tail_gamma = charge_tail_gamma + hist_temp_gamma->GetBinContent(j);
				}
			}
			// std::cout << "charge_tail_n_gamma = " << charge_tail_n_gamma << "\n";
			// std::cout << "charge_tail_gamma = " << charge_tail_gamma << "\n";
			// sleep(1);
			integral_tail_n_gamma.push_back(charge_tail_n_gamma);
			integral_tail_gamma.push_back(charge_tail_gamma);
			min_tail.push_back(wf_charge_tail_min);

			charge_tail_n_gamma = 0.;
			charge_tail_gamma = 0.;
		}
		wf_charge_tail_min = temp;
		
		// std::cout << "integral_tail_n_gamma[0] = " << integral_tail_n_gamma[0] << "\n";
		// std::cout << "integral_tail_gamma[0] = " << integral_tail_gamma[0] << "\n";
		// sleep(1);
		
		q_ratio_n_gamma = charge_tail_n_gamma/charge_total_n_gamma;
		q_ratio_gamma = charge_tail_gamma/charge_total_gamma;

		#ifdef USE_SIGNAL_HEIGHT
		charge_total_n_gamma = maximum_n_gamma;
		charge_total_gamma = maximum_gamma;
		#endif

		#ifdef HEIGHT_THRESHOLD
		if(maximum_n_gamma > height_threshold_n_gamma)
		{
			// hist_Q_ratio_n_gamma->Fill(charge_total_n_gamma, q_ratio_n_gamma);
			// graph_Q_ratio_n_gamma->AddPoint(charge_total_n_gamma, q_ratio_n_gamma);
		}
		if(maximum_gamma > height_threshold_gamma)
		{
			// hist_Q_ratio_gamma->Fill(charge_total_gamma, q_ratio_gamma);
			// graph_Q_ratio_gamma->AddPoint(charge_total_gamma, q_ratio_gamma);
		}
		#else
		hist_Q_ratio_n_gamma->Fill(charge_total_n_gamma, q_ratio_n_gamma);
		graph_Q_ratio_n_gamma->AddPoint(charge_total_n_gamma, q_ratio_n_gamma);
		hist_Q_ratio_gamma->Fill(charge_total_gamma, q_ratio_gamma);
		graph_Q_ratio_gamma->AddPoint(charge_total_gamma, q_ratio_gamma);
		#endif

		#ifdef DRAW_EACH_WF
			//! Turn this on to check each waveforms
			//? Clone for hist_temp or hist_temp_aligned;
			TH1D* hist_temp_clone_n_gamma = (TH1D*)hist_temp_n_gamma->Clone();
			name = Form("hist_temp_clone_n_gamma_%d",i);
			hist_temp_clone_n_gamma->SetName(name);
			hist_temp_clone_n_gamma->SetTitle(name);
			hist_temp_clone_n_gamma->SetDirectory(0);
			hist_temp_clone_n_gamma->GetXaxis()->SetRangeUser(wf_min_view, wf_max_view);
			hist_temp_clone_n_gamma->GetYaxis()->SetRangeUser(0, view_hist_height);
			
			TH1D* hist_temp_clone_gamma = (TH1D*)hist_temp_gamma->Clone();
			name = Form ("hist_temp_clone_gamma_%d",i);
			hist_temp_clone_gamma->SetName(name);
			hist_temp_clone_gamma->SetTitle(name);
			hist_temp_clone_gamma->SetDirectory(0);
			hist_temp_clone_gamma->GetXaxis()->SetRangeUser(wf_min_view, wf_max_view);
			hist_temp_clone_gamma->GetYaxis()->SetRangeUser(0, view_hist_height);
			
			canvas_2->cd(1);
			
			//hist_temp->Draw();
			// sleep(1);
			int colorIndex = i % 50 + 1;  // ROOT has colors 1–50 (looping)
			
			hist_temp_clone_n_gamma->SetLineColorAlpha(colorIndex, 0.05);  // Faint line
			hist_temp_clone_n_gamma->SetLineWidth(1);
			hist_temp_clone_n_gamma->Draw("same");
			
			
			canvas_2->cd(2);
			hist_temp_clone_gamma->SetLineColorAlpha(colorIndex, 0.05);  // Faint line
			hist_temp_clone_gamma->SetLineWidth(1);
			hist_temp_clone_gamma->Draw("same");
			
			canvas_2->Modified();
			canvas_2->Update();
			//! Up to here
		#endif

		delete hist_temp_n_gamma;
		delete hist_temp_gamma;
		
		outTree->Fill();
		h_ave_n_gamma.clear();
		h_ave_gamma.clear();
		
		integral_tail_n_gamma.clear();
		integral_tail_gamma.clear();
		min_tail.clear();
	}

	outTree->Write();
	outFile->Close();

	std::cout << "wf_charge_total_min = " << wf_charge_total_min << "\n";
	std::cout << "wf_charge_total_max = " << wf_charge_total_max << "\n";
	std::cout << "wf_charge_tail_min increases from = " << wf_charge_tail_min << " to "<< wf_charge_tail_min_final <<"\n";
	std::cout << "with increment = " << wf_charge_tail_min_incre << "\n";

	std::cout << "time: " << timer->RealTime() << " seconds \n";
	
}
