void analyse(){

    TFile* file = new TFile(~/data/56V_1250MS.root, "recreate");
    TTree* tree = (TTree*)file->Get("wf", "wf");

    size_t eventNr = 0;
    vector<float> vTime;
    vector<float> vVoltage_1;
    vector<float> vVoltage_2;
    double integral_1;
    double integral_2;

    tree->SetBranchAddress("event", &eventNr);
    tree->SetBranchAddress("Time", &vTime);
    tree->SetBranchAddress("Voltage_1", &vVoltage_1);
    tree->SetBranchAddress("Voltage_2", &vVoltage_2);
    tree->SetBranchAddress("Integral_1", &integral_1);
    tree->SetBranchAddress("Integral_2", &integral_2);

    for (int i = 0; i<tree->GetEntries(); i++)
    {
        tree->GetEntry(i);
        double 
        for (int j = 0; j<vVoltage_1.size(); j++)
        {

        }
    }
}