#include <TVector2.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <iostream>
#include <map>
#include <algorithm>

using std::fstream;
using std::string;
using std::getline;
using std::istringstream;
using std::vector;
using std::pair;
using std::make_pair;
using std::cout;
using std::endl;
using std::stof;
using std::map;
using std::cerr;
using std::max;

void analyse(){

    TFile* file = new TFile("~/data/wf_files/input/root_files/56V_1250MS.root", "read");
    TTree* tree = (TTree*)file->Get("wf");

    Int_t eventNr = 0;
    
	std::vector<float>* vTime = nullptr;
	std::vector<float>* vVoltage_1 = nullptr;
	std::vector<float>* vVoltage_2 = nullptr;

    float baseline = -150.;

    double integral_1;
    double integral_2;

    tree->SetBranchAddress("event", &eventNr);
    
    tree->SetBranchAddress("Time", &vTime);
    tree->SetBranchAddress("Voltage_1", &vVoltage_1);
    tree->SetBranchAddress("Voltage_2", &vVoltage_2);
    tree->SetBranchAddress("Integral_1", &integral_1);
    tree->SetBranchAddress("Integral_2", &integral_2);
    
    TH1D* h_baseline_1 = new TH1D("h_baseline_1", "h_baseline_1", 100, -2.,30.);
    h_baseline_1->SetTitle("Baseline Channel 1; Voltage Difference (mV); Counts");
    TH1D* h_baseline_2 = new TH1D("h_baseline_2", "h_baseline_2", 100, -2.,30.);
    h_baseline_2->SetTitle("Baseline Channel 2; Voltage Difference (mV); Counts");

    TH1D* h_rise_time_1 = new TH1D("h_rise_time_1", "h_rise_time_1", 50, 0.,5.);
    h_rise_time_1->SetTitle("Rise Time Channel 1; Time (ns); Counts");
    TH1D* h_rise_time_2 = new TH1D("h_rise_time_2", "h_rise_time_2", 50, 0.,5.);
    h_rise_time_2->SetTitle("Rise Time Channel 2; Time (ns); Counts");

    TH1D* h_diff_rise_time = new TH1D("h_diff_rise_time", "h_diff_rise_time", 100, -5.,5.);
    h_diff_rise_time->SetTitle("Difference in Rise Time; Time Difference (ns); Counts");

    TH1D* h_start_1 = new TH1D("h_start_1", "h_start_1", 100, -2., 2.);
    h_start_1->SetTitle("Signal Start Time Channel 1; Time (ns); Counts");
    TH1D* h_start_2 = new TH1D("h_start_2", "h_start_2", 100, -2., 2.);
    h_start_2->SetTitle("Signal Start Time Channel 2; Time (ns); Counts");

    // float min_set = 5000.;
    // float max_set = 10000.;
    // float min_set = 15000.;
    // float max_set = 500000.;
    
    float min_set = 10000.;
    float max_set = 500000.;

    float max_1_total = 0.;
    float max_2_total = 0.;

    for (int i = 0; i<tree->GetEntries(); i++)
    // for (int i = 0; i<100; i++)
    {
        std::cout << "entry no. " << i << "\n";
        tree->GetEntry(i);

        float start_1 = 0.;
        float start_2 = 0.;
        float max_1 = -1e9;
        float max_2 = -1e9;
        float max_1_time = -1e9;
        float max_2_time = -1e9;

        float trigger_thresh = 3.5;

        if (integral_1 > min_set && integral_1 < max_set && integral_2 > min_set && integral_2 < max_set)
        {

            //! Baseline calculation
            for (int j = 0; j<(*vTime).size()-1600; j++)
            // for (int j = 0; j<(*vTime).size(); j++)
            {
                // std::cout << "(*vVoltage_1)[" << j << "] - baseline = " << (*vVoltage_1)[j] - baseline << "\n";
                // std::cout << "(*vVoltage_2)[" << j << "] - baseline = " << (*vVoltage_2)[j] - baseline << "\n";
                h_baseline_1->Fill((*vVoltage_1)[j] - baseline);
                h_baseline_2->Fill((*vVoltage_2)[j] - baseline);
            }

            //! Timestamp of the start of the signal
            for (int j = 0; j<(*vTime).size(); j++)
            {
                if( ((*vVoltage_1)[j] - baseline) > trigger_thresh )
                {
                    start_1 = (*vTime)[j];
                    // std::cout << "trigger is at " << start << "\n";
                    h_start_1->Fill(start_1);
                    break;
                }
            }

            for (int j = 0; j<(*vTime).size(); j++)
            {
                if( ((*vVoltage_2)[j] - baseline) > trigger_thresh )
                {
                    start_2 = (*vTime)[j];
                    h_start_2->Fill(start_2);
                    // std::cout << "trigger is at " << start << "\n";
                    break;
                }
            }
            
            //! Time stamp of the maximum of the signal
            for (int j = 0; j<(*vTime).size(); j++)
            {            
                if( (*vVoltage_1)[j] > max_1 )
                {
                    // std::cout << "(*vVoltage_1)[" << j << "] = " << (*vVoltage_1)[j] << "\n";
                    max_1 = (*vVoltage_1)[j];
                    max_1_time = (*vTime)[j];
                    // std::cout << "max_1 = " << max_1 << "\n";
                }
                
                if( (*vVoltage_2)[j] > max_2 )
                {
                    // std::cout << "(*vVoltage_2)[" << j << "] = " << (*vVoltage_2)[j] << "\n";
                    max_2 = (*vVoltage_2)[j];
                    max_2_time = (*vTime)[j];
                    // std::cout << "max_2 = " << max_2 << "\n";
                }
            }
            max_1_total += max_1;
            max_2_total += max_2;
            
            float rise_time_1 = (max_1_time-start_1)/2.;
            float rise_time_2 = (max_2_time-start_2)/2.;
            h_rise_time_1->Fill(rise_time_1);
            h_rise_time_2->Fill(rise_time_2);
            h_diff_rise_time->Fill( rise_time_1 - rise_time_2 );
            // sleep(5);
        }
    }

    std::cout << "Average max_1 = " << max_1_total/tree->GetEntries() - baseline << "\n";
    float SNR_1 = (max_1_total/tree->GetEntries() - baseline)/3.0;
    std::cout << "Signal-to-noise ratio = " << SNR_1 << "\n";
    std::cout << "Average max_2 = " << max_2_total/tree->GetEntries() - baseline << "\n";
    float SNR_2 = (max_2_total/tree->GetEntries() - baseline)/3.0;
    std::cout << "Signal-to-noise ratio = " << SNR_2 << "\n";
    
    std::cout << "Jitter_1 = " << 1.69/SNR_1 << " ns\n";
    std::cout << "Jitter_2 = " << 1.57/SNR_2 << " ns\n";
    TCanvas* c1 = new TCanvas("c1", "c1");
    c1->Divide(2,1);
    c1->cd(1);
    h_baseline_1->Draw();
    c1->cd(2);
    h_baseline_2->Draw();

    TCanvas* c2 = new TCanvas("c2", "c2");
    c2->Divide(2,1);
    c2->cd(1);
    h_rise_time_1->Draw();
    c2->cd(2);
    h_rise_time_2->Draw();

    TCanvas* c3 = new TCanvas("c3", "c3");
    c3->cd();
    h_diff_rise_time->Draw();
    h_diff_rise_time->SetStats(0);
    c3->SetLogy();
    h_diff_rise_time->Fit("gaus");
    double sigma = h_diff_rise_time->GetFunction("gaus")->GetParameter(2);
    cout << "Timing resolution (sigma) = " << sigma/sqrt(2) << " ns" << endl;

    TCanvas* c4 = new TCanvas("c4", "c4");
    c4->Divide(2,1);
    c4->cd(1);
    h_start_1->Draw();
    c4->cd(2);
    h_start_2->Draw();

    TCanvas* c5 = new TCanvas("c5", "c5");
    TGraph* time_res_by_integral = new TGraph();
    time_res_by_integral->SetTitle("Time Resolution vs Integral; Integral (a. unit); Time Resolution (ns)");
    time_res_by_integral->SetPoint(0, 7500., 0.478);
    time_res_by_integral->SetPoint(1, 12500., 0.331);
    time_res_by_integral->SetPoint(2, 17500., 0.298);
    time_res_by_integral->SetPoint(3, 22500., 0.240);
    time_res_by_integral->SetPoint(4, 27500., 0.221);
    time_res_by_integral->SetPoint(5, 32500., 0.250);
    // time_res_by_integral->SetPoint(6, 37500., 0.079);
    time_res_by_integral->SetMarkerStyle(21);
    c5->cd();
    time_res_by_integral->Draw("AP");
}