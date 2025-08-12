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
#define INLIST "list_bc404_na22.txt"
#define OUTROOT "out_bc404_na22.root"

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

// Function to calculate the integral (trapezoidal method)
double calculateIntegral(const vector<float>& time, const vector<float>& voltage) {
    if (time.size() < 2 || voltage.size() < 2 || time.size() != voltage.size()) 
        return 0.0;
    
    double integral = 0.0;
    for (size_t i = 1; i < time.size(); ++i) {
        double dt = time[i] - time[i-1];
        double v_prom = -(voltage[i] + voltage[i-1]) / 2.0;
        integral += v_prom * dt;
    }
    return integral;
}

int prueba4() {
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
    double integral = 0.;
    vector<float> vTraceX;
    vector<float> vTraceY;


    TFile* outputFile = new TFile(OUTROOT, "recreate");
    TTree* tree = new TTree("wf", "wf");
    tree->Branch("event", &eventNr, "event/I");
    tree->Branch("traceX", &vTraceX);
    tree->Branch("traceY", &vTraceY);
    tree->Branch("integral", &integral, "integral/D");

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
        vTraceX.clear();
        vTraceY.clear();
        int skip = 0;
        integral = 0.0;
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
                    vTraceX.push_back(time);
                    vTraceY.push_back(voltage);   
                                  
                } 
                catch (...) {
                    
                }
            }            
        }
        archive.close();
        integral = calculateIntegral(vTraceX, vTraceY);
        //printf("%f\n", integral);

        // 8. Fill tree
        tree->Fill();
        
    }

    // 9. Save and close
    outputFile->Write();
    outputFile->Close();

    #ifdef DEBUG
        delete cnCommon;
    #endif

    return 0;
}