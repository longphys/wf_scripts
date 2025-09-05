#include <TVector2.h>
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

// Code parameters
#define LINES_TO_SKIP 3
//#define BASELINE 180
#define DELIM ','  // We use ',' to separate values
// #define INLIST "/home/long/scripts/wf_scripts/picoscope/lists/list_bc404_cs137_20250725-0006.txt"
// #define INLIST "/home/long/scripts/wf_scripts/picoscope/lists/list_bc404_pu_c13_20250723-0002.txt"
// #define INLIST "/home/long/scripts/wf_scripts/picoscope/lists/list_stilbene_pu_c13_20250723-0004.txt"
#define INLIST "/home/long/scripts/wf_scripts/picoscope/lists/list_stilbene_na22_20250723-0003.txt"
// #define OUTROOT "/home/long/data/wf_files/input/root_files/bc404_cs137_20250725-0006.root"
// #define OUTROOT "/home/long/data/wf_files/input/root_files/stilbene_pu_c13_20250723-0004.root"
#define OUTROOT "/home/long/data/wf_files/input/root_files/stilbene_na22_20250723-0003.root"
// #define OUTROOT "./test_wf.root"

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

int read_csv_1() {
    // 1. Read file list
    fstream inList(INLIST);
    string line;
    vector<string> vWfFileList;
    int counter = 0;

    while (getline(inList, line))
    {
        #ifdef FILES_TO_READ 
            if (counter == FILES_TO_READ)
            {
                break;
            }
        #endif /*FILES_TO_READ*/

        vWfFileList.push_back(line);
        counter++;
    }
    inList.close();


    size_t eventNr = 0;
    vector<float> vTime;
    vector<float> vVoltage;
    double integral;


    TFile* outputFile = new TFile(OUTROOT, "recreate");
    TTree* tree = new TTree("wf", "wf");
    tree->Branch("event", &eventNr, "event/I");
    tree->Branch("Time", &vTime);
    tree->Branch("Voltage", &vVoltage);
    tree->Branch("Integral", &integral, "Integral/D");

    #ifdef DEBUG
        TCanvas* cnCommon = new TCanvas("cnCommon", "cnCommon", 700, 700);
    #endif

     std:: cout  <<"processed values: " ;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 3. Process each file
    for (size_t i = 0; i < vWfFileList.size(); i++) {
    
        fstream archive(vWfFileList[i]);
        std:: cout <<"processed file: "<< vWfFileList[i] << std:: endl ;
        if (!archive.is_open()) {
            cerr << "Error opening: " << vWfFileList[i] << endl;
            //continue;
            return 1;
        }

        // Reset variables for each event
        vTime.clear();
        vVoltage.clear();
        int skip = 0;
        eventNr = i;
        
        // 4. Read headings (skip initial lines)
        for (int j = 0; j < LINES_TO_SKIP; j++) {
            if (!getline(archive, line)) break;
        }

        // 5. Process data
        map<string, size_t> indiceColumns;
        size_t idx_time = 0, idx_voltage = 1; // Default values
        
        if (getline(archive, line)) {
            // Parse headers
            istringstream ss(line);
            string column;
            vector<string> columns;
            
            while (getline(ss, column, DELIM)) {
                column.erase(remove(column.begin(), column.end(), ' '), column.end());
                columns.push_back(column);
                indiceColumns[column] = columns.size() - 1;
            }
            
            // Identify column indexes
            if (indiceColumns.find("Time") != indiceColumns.end() && 
                indiceColumns.find("ChannelA") != indiceColumns.end()) {
                idx_time = indiceColumns["Time"];
                idx_voltage = indiceColumns["ChannelA"];
            }
        }

        // 6. Read data
        while (getline(archive, line)) {
            istringstream ss(line);
            string value;
            vector<string> fila_actual;
            
            while (getline(ss, value, DELIM)) {
                fila_actual.push_back(value);
            }
            
            if (fila_actual.size() > max(idx_time, idx_voltage)) {
                try {
                    float time = stof(fila_actual[idx_time]);
                    float voltage = stof(fila_actual[idx_voltage]);
                
                    //printf("%f %f\n", time, voltage);
                    vTime.push_back(time);
                    vVoltage.push_back(voltage);   
                                  
                } 
                catch (...) {
                    
                }
            }
            for (int iii = 0; iii <(int)vVoltage.size(); iii++)
            {
                integral+=190.-vVoltage.at(iii);
            }            
        }
        archive.close();

        // 8. Fill tree
        tree->Fill();
        integral = 0.;
    }

    // 9. Save and close
    outputFile->Write();
    outputFile->Close();

    #ifdef DEBUG
        delete cnCommon;
    #endif

    return 0;
}