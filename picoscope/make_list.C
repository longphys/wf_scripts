#include <fstream>
#include <iostream>
#include <iomanip>
#include "TString.h"

void make_list() {
    // --- settings ---
    TString out_txt = "/home/long/scripts/wf_scripts/picoscope/lists/list_csi_na22_08102025.txt";
    TString base_dir = "/home/long/data/wf_files/picoscope/new/csi_na22_08102025/csi_na22_08102025/";
    TString prefix   = "csi_na22_08102025_";
    int first_idx    = 1;
    // int last_idx     = 32760;  // change to however many files you have
    int last_idx     = 50000;  // change to however many files you have
    TString ext      = ".csv";

    // --- open output file ---
    std::ofstream outfile(out_txt.Data());
    if (!outfile.is_open()) {
        std::cerr << "Error: Cannot open output file " << out_txt << std::endl;
        return;
    }

    // --- write file list ---
    for (int i = first_idx; i <= last_idx; i++) {
        outfile << base_dir
                << prefix
                << std::setw(5) << std::setfill('0') << i
                << ext << std::endl;
    }

    outfile.close();
    std::cout << "Wrote list to " << out_txt << std::endl;
}
