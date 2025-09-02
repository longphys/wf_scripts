#include <fstream>
#include <iostream>
#include <iomanip>
#include "TString.h"

void make_list() {
    // --- settings ---
    TString out_txt = "/home/long/scripts/wf_scripts/picoscope/lists/list_bc404_cs137_20250725-0006.txt";
    TString base_dir = "/home/long/data/wf_files/input/picoscope_bc404_cs137/20250725-0006/";
    TString prefix   = "20250725-0006_";
    int first_idx    = 1;
    int last_idx     = 100000;  // change to however many files you have
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
