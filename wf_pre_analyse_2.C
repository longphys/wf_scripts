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

#define ORIGINAL_WF
#define DRAW_EACH_WF
#define APPLY_AVERAGE
// #define APPLY_ARC
#define UN_NORMALIZE_WF
#define HEIGHT_THRESHOLD
double height_threshold_n_gamma = 10.;
double height_threshold_gamma = 10.;
#define USE_SIGNAL_HEIGHT // Use the pulse height to draw the 2D histograms instead of the integral.

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
int wf_max_view = wf_max;

void wf_pre_analyse_2()
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
	// UShort_t time_n_gamma[n_wf];
	// UShort_t time_gamma[n_wf];
	tree_n_gamma->SetBranchAddress("Voltage", &time_n_gamma);
	tree_gamma->SetBranchAddress("Voltage", &time_gamma);
	
	TCanvas* canvas_1 = new TCanvas("canvas_1", "canvas_1", 1400, 700);
	canvas_1->Divide(2,1);
	
	TCanvas* canvas_2 = new TCanvas("canvas_2", "canvas_2", 1400, 700);
	canvas_2->Divide(2,1);
	
	TCanvas* canvas_3 = new TCanvas("canvas_3", "canvas_3", 1500, 450);
	canvas_3->Divide(3,1);
	//const int n = tree_gamma->GetEntries();
	//const int n = tree_n_gamma->GetEntries();

	//! Short window
	const int n = 28000;
	//! Long window
	// const int n = 50000;
	// const int n = 100;
	
	std::cout << "Number of entries: " << n << "\n";

	TH1D* hist_spectrum_n_gamma = new TH1D("spectrum_n_gamma", "spectrum", 500, 0., 500000.);
	hist_spectrum_n_gamma->GetXaxis()->SetTitle("Amplitude (Channels)");
	hist_spectrum_n_gamma->GetYaxis()->SetTitle("Count/Channel");
	hist_spectrum_n_gamma->GetXaxis()->CenterTitle();
	hist_spectrum_n_gamma->GetYaxis()->CenterTitle();
	
	TH1D* hist_spectrum_gamma = new TH1D("spectrum_gamma", "spectrum", 500, 0., 500000.);
	hist_spectrum_gamma->GetXaxis()->SetTitle("Amplitude (Channels)");
	hist_spectrum_gamma->GetYaxis()->SetTitle("Count/Channel");
	hist_spectrum_gamma->GetXaxis()->CenterTitle();
	hist_spectrum_gamma->GetYaxis()->CenterTitle();
	
	#ifdef USE_SIGNAL_HEIGHT
	TH2D* hist_Q_ratio_n_gamma = new TH2D("Q-ratio Map n-gamma", "Q-ratio Map n-gamma", 100, 0., 300., 800, 0., 1);
	TH2D* hist_Q_ratio_gamma = new TH2D("Q-ratio Map gamma", "Q-ratio Map gamma", 100, 0., 300., 800, 0., 1);
	hist_Q_ratio_n_gamma->SetTitle("Q-ratio Map n_gamma;Pulse height (a. unit);Q-ratio");
	hist_Q_ratio_gamma->SetTitle("Q-ratio Map gamma;Pulse height (a. unit);Q-ratio");
	#else
	TH2D* hist_Q_ratio_n_gamma = new TH2D("Q-ratio Map n-gamma", "Q-ratio Map n-gamma", 500, 0., 50000., 800, 0., 1);
	TH2D* hist_Q_ratio_gamma = new TH2D("Q-ratio Map gamma", "Q-ratio Map gamma", 500, 0., 50000., 800, 0., 1);
	hist_Q_ratio_n_gamma->SetTitle("Q-ratio Map n_gamma;Charge (a. unit);Q-ratio");
	hist_Q_ratio_gamma->SetTitle("Q-ratio Map gamma;Charge (a. unit);Q-ratio");
	#endif
	hist_Q_ratio_n_gamma->GetXaxis()->CenterTitle();
	hist_Q_ratio_n_gamma->GetYaxis()->CenterTitle();
	
	hist_Q_ratio_gamma->GetXaxis()->CenterTitle();
	hist_Q_ratio_gamma->GetYaxis()->CenterTitle();
	
	TGraph* graph_Q_ratio_n_gamma = new TGraph();
	
	TGraph* graph_Q_ratio_gamma = new TGraph();
	
	graph_Q_ratio_n_gamma->SetMarkerColor(kRed);
	graph_Q_ratio_gamma->SetMarkerColor(kBlue);
	#ifdef USE_SIGNAL_HEIGHT
	graph_Q_ratio_n_gamma->SetMarkerStyle(3);
	graph_Q_ratio_gamma->SetMarkerStyle(3);
	#else
	graph_Q_ratio_n_gamma->SetMarkerStyle(1);
	graph_Q_ratio_gamma->SetMarkerStyle(1);
	#endif
	
	TMultiGraph* graph_both = new TMultiGraph();
	#ifdef USE_SIGNAL_HEIGHT
	graph_both->GetXaxis()->SetTitle("Pulse height (a. unit)");
	#else
	graph_both->GetXaxis()->SetTitle("Charge (a. unit)");
	#endif
	graph_both->GetYaxis()->SetTitle("Q-ratio");
	graph_both->GetXaxis()->CenterTitle();
	graph_both->GetYaxis()->CenterTitle();
	graph_both->Add(graph_Q_ratio_n_gamma, "AP");
	graph_both->Add(graph_Q_ratio_gamma, "AP");

	TFile* outFile = new TFile("BC404_graphs.root","recreate");
	TTree* outTree = new TTree("results","results");
	std::vector<double> height;
	outTree->Branch("height", &height);
	
	//! Short window
	int wf_charge_total_min = 600;
	int wf_charge_total_max = 2000;
	int wf_charge_tail_min = 1000;
	int wf_charge_tail_min_final = 1300;
	double baseline = 188.;
	
	//! Long window
	// int wf_charge_total_min = 1200;
	// int wf_charge_total_max = 2450;
	// int wf_charge_tail_min = 1250;
	// int wf_charge_tail_min_final = 2400;
	// double baseline = 185.;

	int wf_charge_tail_min_incre = 10;
	// int wf_charge_tail_min_incre = 50;
	int times = (wf_charge_tail_min_final-wf_charge_tail_min)/wf_charge_tail_min_incre;
	// std::cout << "times = " << times << "\n";
	std::cout << "wf_charge_total_min = " << wf_charge_total_min << "\n";
	std::cout << "wf_charge_total_max = " << wf_charge_total_max << "\n";
	std::cout << "wf_charge_tail_min will increase from = " << wf_charge_tail_min << " to "<< wf_charge_tail_min_final <<"\n";
	std::cout << "with increment = " << wf_charge_tail_min_incre << "\n";

	TString name;
for (int ii = 0; ii < times; ii++)
{
	std::cout << "iteration: " << ii+1 << "/" << times << "\n";
	wf_charge_tail_min += wf_charge_tail_min_incre;
std::cout << "wf_charge_tail_min = " << wf_charge_tail_min << "\n";
	// for (int i = 0; i < n; i++)
	// for (int i = 9705; i < 9706; i++)
	// for (int i = 1000; i < 1001; i++)
	for (int i = 8781; i < 8782; i++)
	{	
		if (i%5000==0)
		// if (i%50000==0)
		{
			// std::cout << "Entry No. " << i << "\n";
			
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
		// TH1D* hist_temp_gamma = new TH1D(name, name, wf_max-wf_min, wf_min+1500, wf_max+1500);
		tree_gamma->GetEntry(i);
		
		// Get and fill waveforms
		double inverse_height_a;
		double inverse_height_b;

		#ifdef ORIGINAL_WF
		for (int j = 0; j < wf_max-wf_min; j++)
		{
			inverse_height_a = (*time_n_gamma)[j+wf_min]; //! Original waveform
			inverse_height_b = (*time_gamma)[j+wf_min];
			hist_temp_n_gamma->SetBinContent(j+1, inverse_height_a);
			hist_temp_gamma->SetBinContent(j+1, inverse_height_b);
		}
		#else
		for (int j = 0; j < wf_max-wf_min; j++)
		{
			double inverse_height_a = baseline-(*time_n_gamma)[j+wf_min];
			// double inverse_height_a = 183.-(*time_n_gamma)[j+wf_min];
			
			if(inverse_height_a < 0)
			{
				hist_temp_n_gamma->SetBinContent(j+1, 1.);
			}
			else
			{
				hist_temp_n_gamma->SetBinContent(j+1, inverse_height_a);
			}

			double inverse_height_b = baseline-(*time_gamma)[j+wf_min];
			// double inverse_height_b = 183.-(*time_gamma)[j+wf_min];
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
		#ifdef APPLY_AVERAGE
		// for (int j = 0; j < wf_max-wf_min; j+=4)
		for (int j = 0; j < wf_max-wf_min; j++)
		{
			double average_a = 
			(
				hist_temp_n_gamma->GetBinContent(j)+
				hist_temp_n_gamma->GetBinContent(j+1)+
				hist_temp_n_gamma->GetBinContent(j+2)+
				hist_temp_n_gamma->GetBinContent(j+3)
			)/4.0;
			hist_temp_n_gamma->SetBinContent(j, average_a);
			// hist_temp_n_gamma->SetBinContent(j+1, average_a);
			// hist_temp_n_gamma->SetBinContent(j+2, average_a);
			// hist_temp_n_gamma->SetBinContent(j+3, average_a);
			
			if (j<1996)
			{
				// for(int t = 0; t <4; t++)
				for(int t = 0; t <1; t++)
				{
					height.push_back(average_a);
				}
			}

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
		}
		#endif
		
		double scale_factor_n_gamma = n_wf_height/(maximum_n_gamma);
		hist_temp_n_gamma->Scale(scale_factor_n_gamma, "noSW2");
		
		double scale_factor_gamma = n_wf_height/(maximum_gamma);
		hist_temp_gamma->Scale(scale_factor_gamma, "noSW2");

		
		int shift_n_gamma = 0;
		int shift_gamma = 0;

		// Apply ARC
		#ifdef APPLY_ARC
		// std::cout << "ARC_bin_index = " << ApplyARC(hist_temp, 0.95, 2) << "\n";
		
		int arc_trigger_bin_n_gamma = ApplyARC(hist_temp_n_gamma, 0.95, 2);
		int arc_trigger_bin_gamma = ApplyARC(hist_temp_gamma, 0.95, 2);

		// std::cout << "arc_trigger_bin_n_gamma = " << arc_trigger_bin_n_gamma << "\n";
		// std::cout << "arc_trigger_bin_gamma = " << arc_trigger_bin_gamma << "\n";
		// sleep(1);

		int arc_bin_align_to = 600;
		shift_n_gamma = arc_bin_align_to - arc_trigger_bin_n_gamma;
		shift_gamma = arc_bin_align_to - arc_trigger_bin_gamma;
		#endif
		
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

		// Calculation of the charge
		#ifdef UN_NORMALIZE_WF
			//! Turn this off to enable amplitude normalization
			hist_temp_aligned_n_gamma->Scale(1/scale_factor_n_gamma, "noSW2");
			hist_temp_aligned_gamma->Scale(1/scale_factor_gamma, "noSW2");
		#endif

		double charge_total_n_gamma = 0.;
		double charge_tail_n_gamma = 0.;
		double charge_total_gamma = 0.;
		double charge_tail_gamma = 0.;
		
		// fixed bins intergration
		for(int j = wf_charge_total_min; j < wf_charge_total_max; j++)
		{
			charge_total_n_gamma = charge_total_n_gamma + hist_temp_aligned_n_gamma->GetBinContent(j+1);
			charge_total_gamma = charge_total_gamma + hist_temp_aligned_gamma->GetBinContent(j+1);
			if(j > wf_charge_tail_min)
			{
				charge_tail_n_gamma = charge_tail_n_gamma + hist_temp_aligned_n_gamma->GetBinContent(j+1);
				charge_tail_gamma = charge_tail_gamma + hist_temp_aligned_gamma->GetBinContent(j+1);
			}
		}

		double q_ratio_n_gamma = charge_tail_n_gamma/charge_total_n_gamma;
		double q_ratio_gamma = charge_tail_gamma/charge_total_gamma;
		
		//std::cout << "charge total = " << charge_total << "\n";
		//std::cout << "charge tail/charge total = " << q_ratio << "\n";
			
		// Analyze waveforms and fill spectrum

		if(ii == 0)
		{
			hist_spectrum_n_gamma->Fill(charge_total_n_gamma);
			hist_spectrum_gamma->Fill(charge_total_gamma);
		}
		#ifdef USE_SIGNAL_HEIGHT
		charge_total_n_gamma = maximum_n_gamma;
		charge_total_gamma = maximum_gamma;
		#endif
		#ifdef HEIGHT_THRESHOLD
		if(maximum_n_gamma > height_threshold_n_gamma)
		{
			hist_Q_ratio_n_gamma->Fill(charge_total_n_gamma, q_ratio_n_gamma);
			graph_Q_ratio_n_gamma->AddPoint(charge_total_n_gamma, q_ratio_n_gamma);
		}
		if(maximum_gamma > height_threshold_gamma)
		{
			hist_Q_ratio_gamma->Fill(charge_total_gamma, q_ratio_gamma);
			graph_Q_ratio_gamma->AddPoint(charge_total_gamma, q_ratio_gamma);
		}
		#else
		hist_Q_ratio_n_gamma->Fill(charge_total_n_gamma, q_ratio_n_gamma);
		graph_Q_ratio_n_gamma->AddPoint(charge_total_n_gamma, q_ratio_n_gamma);
		hist_Q_ratio_gamma->Fill(charge_total_gamma, q_ratio_gamma);
		graph_Q_ratio_gamma->AddPoint(charge_total_gamma, q_ratio_gamma);
		#endif
		// Draw Histograms
		if (i == 0 && ii == 0)
		{
			canvas_1->cd(1);
			hist_spectrum_n_gamma->Draw();
			canvas_1->cd(2);
			hist_spectrum_gamma->Draw();
		}

		#ifdef DRAW_EACH_WF
			//! Turn this on to check each waveforms
			//? Clone for hist_temp or hist_temp_aligned;
			TH1D* hist_temp_clone_n_gamma = (TH1D*)hist_temp_aligned_n_gamma->Clone();
			name = Form("hist_temp_clone_n_gamma_%d",i);
			hist_temp_clone_n_gamma->SetName(name);
			hist_temp_clone_n_gamma->SetTitle(name);
			hist_temp_clone_n_gamma->SetDirectory(0);
			hist_temp_clone_n_gamma->GetXaxis()->SetRangeUser(wf_min_view, wf_max_view);
			hist_temp_clone_n_gamma->GetYaxis()->SetRangeUser(0, view_hist_height);
			
			TH1D* hist_temp_clone_gamma = (TH1D*)hist_temp_aligned_gamma->Clone();
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
			// cut_1->Draw("same");
			canvas_3->cd(2);
			hist_Q_ratio_gamma->Draw();
			// cut_1->Draw("same");
			canvas_3->cd(3);
			graph_both->Draw("a");
			
			TLegend *legend_1 = new TLegend(0.6, 0.2, 0.9, 0.4);
			//legend_1->SetHeader("test", "C");
			legend_1->SetBorderSize(2);
			legend_1->SetTextSize(0.025);
			legend_1->AddEntry(hist_Q_ratio_n_gamma, "some neutron source", "p");
			legend_1->AddEntry(hist_Q_ratio_gamma, "some gamma source", "p");
			legend_1->Draw();
		}
		graph_both->GetXaxis()->SetRangeUser(
				hist_Q_ratio_n_gamma->GetXaxis()->GetXmin(), 
				hist_Q_ratio_n_gamma->GetXaxis()->GetXmax());
		graph_both->GetYaxis()->SetRangeUser(
				hist_Q_ratio_n_gamma->GetYaxis()->GetXmin(), 
				hist_Q_ratio_n_gamma->GetYaxis()->GetXmax());
				
				// break;
				outTree->Fill();
				height.clear();
	}

	if(ii==0)
	{
		hist_spectrum_n_gamma->Write();
		hist_spectrum_gamma->Write();
	}
	
	hist_Q_ratio_n_gamma->SetName(Form("n_gamma_%d_%d_tail_%d_%d", wf_charge_total_min, wf_charge_total_max, wf_charge_tail_min, wf_charge_total_max));
	hist_Q_ratio_gamma->SetName(Form("gamma_%d_%d_tail_%d_%d", wf_charge_total_min, wf_charge_total_max, wf_charge_tail_min, wf_charge_total_max));
	graph_both->SetName(Form("both_%d_%d_tail_%d_%d", wf_charge_total_min, wf_charge_total_max, wf_charge_tail_min, wf_charge_total_max));

	hist_Q_ratio_n_gamma->Write();
	hist_Q_ratio_gamma->Write();
	graph_both->Write();

	hist_Q_ratio_n_gamma->Reset();
	hist_Q_ratio_gamma->Reset();
	// Clear the TMultiGraph
	graph_Q_ratio_n_gamma->Set(0);
	graph_Q_ratio_gamma->Set(0);
	//
}
outTree->Write();
outFile->Close();
	std::cout << "time: " << timer->RealTime() << " seconds \n";
}
