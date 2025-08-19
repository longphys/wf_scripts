#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

int read_csv_2() {

    // Open output ROOT file
    TFile* outputFile = new TFile("~/data/wf_files/output/wf_picoscope_pu_c13.root", "recreate");
    TTree* tree = new TTree("tt", "Waveforms");

    // Variables for branches
    double time = 0;
    double voltage = 0;
    int fileIndex = 0;

    tree->Branch("fileID", &fileIndex, "fileID/I");
    tree->Branch("time", &time, "time/D");
    tree->Branch("voltage", &voltage, "voltage/D");

    // Read the list of files
    ifstream list("list_bc404_pu_c13.txt");

    vector<string> files;
    string filename;
    while (getline(list, filename)) {
        if (!filename.empty())
            files.push_back(filename);
    }
    list.close();

    // Process each file
    for (size_t f = 0; f < files.size(); f++) {
        ifstream infile(files[f]);
        if (!infile.is_open()) {
            cerr << "Could not open file: " << files[f] << endl;
            continue;
        }

        cout << "Processing " << files[f] << endl;

        // Skip first 3 lines
        string line;
        for (int i = 0; i < 3 && getline(infile, line); i++);

        // Read lines
        while (getline(infile, line)) {
            if (line.empty()) continue;

            istringstream ss(line);
            string col1, col2;

            if (!getline(ss, col1, ',')) continue;
            if (!getline(ss, col2, ',')) continue;

            try {
                time = stod(col1);
                voltage = stod(col2);
            } catch (...) {
                continue; // skip malformed
            }

            fileIndex = static_cast<int>(f);
            tree->Fill();
        }
    }

    // Save ROOT file
    outputFile->cd();
    tree->Write();
    outputFile->Close();

    return 0;
}
