#include <fstream>
#include <iostream>
#include <iomanip>
#include "TString.h"

void make_list() {
    // --- settings ---
    TString out_txt = "/home/user/data/list_new.txt";
    // TString base_dir = "/home/user/data/20251030_56V_1250MS/";
    // TString base_dir = "/home/user/data/20251030_58V_1250MS_128/";
    TString base_dir = "/home/user/data/20251030_58V_1250MS_132/";
    // TString prefix   = "20251030_56V_1250MS_";
    // TString prefix   = "20251030_58V_1250MS_128_";
    TString prefix   = "20251030_58V_1250MS_132_";
    int first_idx    = 1;
    // int last_idx     = 32760;  // change to however many files you have
    int last_idx     = 5063;  // change to however many files you have
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
                // << std::setw(5) << std::setfill('0') << i
                << std::setw(4) << std::setfill('0') << i
                << ext << std::endl;
    }

    outfile.close();
    std::cout << "Wrote list to " << out_txt << std::endl;
}
