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
#define DELIM ','  // We use ',' to separate values
#define INLIST "/home/long/data/list_new.txt"
// #define INLIST "/home/long/scripts/wf_scripts/picoscope/lists/list_csi_co60_08102025.txt"
#define OUTROOT "~/data/wf_files/input/root_files/56V_1250MS.root"
// #define OUTROOT "~/data/58V_1250MS_128.root"
// #define OUTROOT "~/data/58V_1250MS_132.root"
// #define OUTROOT "/home/long/data/wf_files/input/root_files/csi_co60_08102025.root"

//! REMEMBER TO CHANGE THE BASELINE
// double baseline = 300.; //gagg
// double baseline = -135.; //plastic
double baseline = -150.; //plastic

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
    vector<float> vVoltage_1;
    vector<float> vVoltage_2;
    double integral_1;
    double integral_2;

    TFile* outputFile = new TFile(OUTROOT, "recreate");
    
    TTree* tree = new TTree("wf", "wf");
    tree->Branch("event", &eventNr, "event/I");
    tree->Branch("Time", &vTime);
    tree->Branch("Voltage_1", &vVoltage_1);
    tree->Branch("Voltage_2", &vVoltage_2);
    tree->Branch("Integral_1", &integral_1, "Integral_1/D");
    tree->Branch("Integral_2", &integral_2, "Integral_2/D");

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
        vVoltage_1.clear();
        vVoltage_2.clear();
        int skip = 0;
        eventNr = i;
        
        // 4. Read headings (skip initial lines)
        for (int j = 0; j < LINES_TO_SKIP; j++) {
            if (!getline(archive, line)) break;
        }

        // 5. Process data
        map<string, size_t> indiceColumns;
        size_t idx_time = 0, idx_voltage_1 = 1, idx_voltage_2 = 2; // Default values
        
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
                indiceColumns.find("ChannelA") != indiceColumns.end() &&
                indiceColumns.find("ChannelB") != indiceColumns.end()) {
                idx_time = indiceColumns["Time"];
                idx_voltage_1 = indiceColumns["ChannelA"];
                idx_voltage_2 = indiceColumns["ChannelB"];
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
            
            // if (fila_actual.size() > max(idx_time, idx_voltage_1, idx_voltage_2)) {
            if (fila_actual.size() > max(idx_time, idx_voltage_2)) {
                try {
                    float time = stof(fila_actual[idx_time]);
                    float voltage_1 = stof(fila_actual[idx_voltage_1]);
                    float voltage_2 = stof(fila_actual[idx_voltage_2]);
                
                    // printf("%f %f %f\n", time, voltage_1, voltage_2);
                    // sleep(1);
                    vTime.push_back(time);
                    vVoltage_1.push_back(voltage_1);
                    vVoltage_2.push_back(voltage_2);   
                                  
                } 
                catch (...) {
                    
                }
            }     
        }
        archive.close();

        for (int iii = 0; iii <(int)vVoltage_1.size(); iii++)
        // for (int iii = 1500; iii <(int)vVoltage.size()-1400; iii++)
        // for (int iii = 1500; iii <2500; iii++)
        // for (int iii = 3700; iii <10000; iii++)
        {
            // integral_1+=baseline-vVoltage_1.at(iii);
            // integral_2+=baseline-vVoltage_2.at(iii);
            integral_1+=vVoltage_1.at(iii)-baseline;
            integral_2+=vVoltage_2.at(iii)-baseline;
        }       

        // 8. Fill tree
        tree->Fill();
        integral_1 = 0.;
        integral_2 = 0.;
    }

    // 9. Save and close
    // outputFile->Write();
    tree->Write();
    outputFile->Close();

    #ifdef DEBUG
        delete cnCommon;
    #endif

    return 0;
}