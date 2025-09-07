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

// #define ORIGINAL_WF //! only for checking the waveform, all analysis will fail.
// #define DRAW_EACH_WF
#define APPLY_FIXED_AVERAGE
// #define APPLY_MOVING_AVERAGE
#define NORMALIZE_WF

//! Short window
const int n_wf = 2000;
//! Long window
// const int n_wf = 2500;
double n_wf_height = 200.;
double view_hist_height = 500.;

// Histogram min and max
int wf_min = 0;
int wf_max = n_wf;

// Viewing window min and max
int wf_min_view = wf_min;
#ifdef APPLY_FIXED_AVERAGE
int wf_max_view = (int)wf_max/4.;
#endif
#ifdef APPLY_MOVING_AVERAGE
int wf_max_view = wf_max;
#endif


void wf_pre_analyse_5()
{
	auto timer = new TStopwatch();
	timer->Start();

	//! Input local files

	//! Short window
	// TFile* file_n_gamma = new TFile("~/data/wf_files/input/root_files/bc404_pu_c13_20250723-0002.root", "read");
	TFile* file_n_gamma = new TFile("~/data/wf_files/input/root_files/stilbene_pu_c13_20250723-0004.root", "read");
	
	//! Long window
	// TFile* file_n_gamma = new TFile("~/data/wf_files/input/root_files/bc404_pu_c13_200ps.root", "read");

	//! Short window
	// TFile* file_gamma = new TFile("~/data/wf_files/input/root_files/bc404_na22_20250723.root", "read");
	TFile* file_gamma = new TFile("~/data/wf_files/input/root_files/stilbene_na22_20250723-0003.root", "read");

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
	
	#ifdef DRAW_EACH_WF
	TCanvas* canvas_2 = new TCanvas("canvas_2", "canvas_2", 1400, 700);
	canvas_2->Divide(2,1);
	#endif
	
	//const int n = tree_gamma->GetEntries();
	//const int n = tree_n_gamma->GetEntries();

	//! Short window
	// const int n = 28000; // plastic
	const int n = 24000; // stilbene
	//! Long window
	// const int n = 50000;
	// const int n = 100;
	
	std::cout << "Number of entries: " << n << "\n";

	TFile* outFile = new TFile("stilbene_wf.root","recreate");
	// TFile* outFile = new TFile("plastic_wf.root","recreate");
	// TFile* outFile = new TFile("wf_test.root","recreate");
	TTree* outTree = new TTree("tree","tree");

	std::vector<double> h_ave_n_gamma;
	std::vector<double> h_ave_gamma;

	double front_n_gamma; // from the maxima to the start of the signal
	double front_gamma;

	std::vector<double> moving_front_n_gamma; // from a moving channel to the start of the signal
	std::vector<double> moving_front_gamma;
	
	double tail_n_gamma; // from the maxima to the end of the signal
	double tail_gamma;

	std::vector<double> moving_tail_n_gamma; // from a moving channel to the end of the signal
	std::vector<double> moving_tail_gamma;
	
	double maximum_n_gamma;
	double maximum_gamma;
	double maximum_n_gamma_pos;
	double maximum_gamma_pos;
	
	outTree->Branch("t_n_gamma", &time_n_gamma); // Time values of a waveform
	outTree->Branch("t_gamma", &time_gamma);
	outTree->Branch("v_n_gamma", &voltage_n_gamma); // Voltage values of a waveform
	outTree->Branch("v_gamma", &voltage_gamma);
	outTree->Branch("h_ave_n_gamma", &h_ave_n_gamma); // Averaged heights of a waveform
	outTree->Branch("h_ave_gamma", &h_ave_gamma);
	
	outTree->Branch("front_n_gamma", &front_n_gamma);
	outTree->Branch("front_gamma", &front_gamma);
	
	outTree->Branch("moving_front_n_gamma", &moving_front_n_gamma);
	outTree->Branch("moving_front_gamma", &moving_front_gamma);
	
	outTree->Branch("tail_n_gamma", &tail_n_gamma);
	outTree->Branch("tail_gamma", &tail_gamma);

	outTree->Branch("moving_tail_n_gamma", &moving_tail_n_gamma);
	outTree->Branch("moving_tail_gamma", &moving_tail_gamma);

	outTree->Branch("maximum_n_gamma", &maximum_n_gamma); // Maximum height of a waveform
	outTree->Branch("maximum_gamma", &maximum_gamma); 
	outTree->Branch("maximum_n_gamma_pos", &maximum_n_gamma_pos); // Position of the Maximum height of a waveform
	outTree->Branch("maximum_gamma_pos", &maximum_gamma_pos); 
	
	//! Short window
	// number of bins away from the maxima.
	int front_start = 100;  // The start of the signal (number of channels away from the maxima)
	int front_end = 0; // Moving end of the front (number of channels away from the maxima)
	int tail_end = 350;  // The end of the signal (number of channels away from the maxima)
	int tail_start = 0; // Moving start of the tail (number of channels away from the maxima)
	int incre_front = 1; // Increment for the front
	int incre_tail = 1; // Increment for the tail

	double baseline = 188.;
	
	//! Long window

	// number of times for the integral scan to run
	int times_front = (front_start-front_end)/incre_front;
	int times_tail = (tail_end-tail_start)/incre_tail;
	std::cout << "times_front = " << times_front << "\n";
	std::cout << "times_tail = " << times_tail << "\n";

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
			
			h_ave_n_gamma.push_back(average_a);

			double average_b = 
			(
				hist_temp_gamma->GetBinContent(j)+
				hist_temp_gamma->GetBinContent(j+1)+
				hist_temp_gamma->GetBinContent(j+2)+
				hist_temp_gamma->GetBinContent(j+3)
			)/4.0;
			
			h_ave_gamma.push_back(average_b);
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
		
		maximum_n_gamma = hist_temp_n_gamma->GetMaximum();
		maximum_gamma = hist_temp_gamma->GetMaximum();
		maximum_n_gamma_pos = hist_temp_n_gamma->GetMaximumBin();
		maximum_gamma_pos = hist_temp_gamma->GetMaximumBin();

		#ifdef NORMALIZE_WF
		double scale_factor_n_gamma = n_wf_height/(maximum_n_gamma);
		hist_temp_n_gamma->Scale(scale_factor_n_gamma, "noSW2");
		
		double scale_factor_gamma = n_wf_height/(maximum_gamma);
		hist_temp_gamma->Scale(scale_factor_gamma, "noSW2");
		#endif

		front_n_gamma = 0.;
		front_gamma = 0.;
		tail_n_gamma = 0.;
		tail_gamma = 0.;
		
		// Calculation of the front
		for(int j = maximum_n_gamma_pos-front_start; j < maximum_n_gamma_pos-front_end; j++)
		{
			if(j < 0)
			{continue;}
			else
			{
				front_n_gamma += hist_temp_n_gamma->GetBinContent(j);
			}
		}
		for(int j = maximum_gamma_pos-front_start; j < maximum_gamma_pos-front_end; j++)
		{
			if(j < 0)
			{continue;}
			else
			{
				front_gamma += hist_temp_gamma->GetBinContent(j);
			}
		}

		// Calculation of the tail		
		for(int j = maximum_n_gamma_pos+tail_start; j < maximum_n_gamma_pos+tail_end; j++)
		{
			if (j > hist_temp_n_gamma->GetNbinsX())
			{break;}
			else
			{
				tail_n_gamma += hist_temp_n_gamma->GetBinContent(j);
			}
		}
		for(int j = maximum_gamma_pos+tail_start; j < maximum_gamma_pos+tail_end; j++)
		{
			if (j > hist_temp_gamma->GetNbinsX())
			{break;}
			else
			{
				tail_gamma += hist_temp_gamma->GetBinContent(j);
			}
		}

		// Calculation of the moving front
		for (int iii = 0; iii < times_front; iii++)
		{
			double moving_front_n_gamma_temp = 0.;
			double moving_front_gamma_temp = 0.;
			for(int j = maximum_n_gamma_pos-front_start; j < maximum_n_gamma_pos-(front_end+iii); j++)
			{
				if(j < 0)
				{continue;}
				else
				{							
					moving_front_n_gamma_temp += hist_temp_n_gamma->GetBinContent(j);
				}
			}			
			for(int j = maximum_gamma_pos-front_start; j < maximum_gamma_pos-(front_end+iii); j++)
			{
				if(j < 0)
				{continue;}
				else
				{							
					moving_front_gamma_temp += hist_temp_gamma->GetBinContent(j);
				}
			}
			// std::cout << "moving_front_n_gamma_temp["<< iii <<"] = " << moving_front_n_gamma_temp << "\n";
			// std::cout << "moving_front_gamma_temp["<< iii <<"] = " << moving_front_gamma_temp << "\n";
			moving_front_n_gamma.push_back(moving_front_n_gamma_temp);
			moving_front_gamma.push_back(moving_front_gamma_temp);
		}		

		// Calculation of the moving tail
		for (int iii = 0; iii < times_tail; iii++)
		{
			double moving_tail_n_gamma_temp = 0.;
			double moving_tail_gamma_temp = 0.;
			for(int j = maximum_n_gamma_pos+(tail_start+iii); j < maximum_n_gamma_pos+tail_end; j++)
			{
				if (j > hist_temp_n_gamma->GetNbinsX())
				{break;}
				else
				{							
					moving_tail_n_gamma_temp += hist_temp_n_gamma->GetBinContent(j);
				}
			}			
			for(int j = maximum_gamma_pos+(tail_start+iii); j < maximum_gamma_pos+tail_end; j++)
			{
				if (j > hist_temp_gamma->GetNbinsX())
				{break;}
				else
				{							
					moving_tail_gamma_temp += hist_temp_gamma->GetBinContent(j);
				}
			}
			// std::cout << "moving_tail_n_gamma_temp["<< iii <<"] = " << moving_tail_n_gamma_temp << "\n";
			// std::cout << "moving_tail_gamma_temp["<< iii <<"] = " << moving_tail_gamma_temp << "\n";
			moving_tail_n_gamma.push_back(moving_tail_n_gamma_temp);
			moving_tail_gamma.push_back(moving_tail_gamma_temp);
		}
		
		// std::cout << "moving_tail_n_gamma[0] = " << moving_tail_n_gamma[0] << "\n";
		// std::cout << "moving_tail_gamma[0] = " << moving_tail_gamma[0] << "\n";
		// sleep(1);

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
		
		moving_front_n_gamma.clear();
		moving_front_gamma.clear();
		moving_tail_n_gamma.clear();
		moving_tail_gamma.clear();
	}

	outTree->Write();
	outFile->Close();

	std::cout << "Integral of the front is from -" << front_start << " to -"<< front_end <<" from the signal maxima.\n";
	std::cout << "Integral of the tail is from +" << tail_start << " to +"<< tail_end <<" from the signal maxima.\n";
	std::cout << "Moving integral of the front is from -" << front_start << " to -x (" << front_start << "<x<" << front_end <<"), from the signal maxima.\n";
	std::cout << "with increment = " << incre_front << "\n";
	std::cout << "Moving integral of the tail is from +y (" << tail_start << "<y<" << tail_end << ") to +" << tail_end << ", from the signal maxima.\n";
	std::cout << "with increment = " << incre_tail << "\n";

	std::cout << "time: " << timer->RealTime() << " seconds \n";
	
}
