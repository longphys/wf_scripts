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

void refill(){
    TFile* f_mea = new TFile("/home/long/scripts/wf_scripts/picoscope/out_bc404_na22.root", "read");
    TTree* t_mea = (TTree*)f_mea->Get("wf");

    double integral;
    t_mea->SetBranchAddress("integral", &integral);

    double lowest_value = -9700.;

    int entries = t_mea->GetEntries();

    TFile* f_out = new TFile("/home/long/scripts/wf_scripts/picoscope/final_out_na22.root", "recreate");
    TTree* t_out = new TTree("Events", "Events");
    double amp;
    t_out->Branch("Amplitude", &amp);

    for(int i = 0; i<entries; i++)
    {
        t_mea->GetEntry(i);
        amp = integral-lowest_value;
        t_out->Fill();
    }

    t_out->Write();

    f_out->Close();
    f_mea->Close();

}