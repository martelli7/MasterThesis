#include <RtypesCore.h>
#include <TFile.h>
#include <TTree.h>
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>

#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <map>

void filterTree()
{
    TFile *inputFile = new TFile("filteredFinal.root", "READ");
    TTree *inputtree = (TTree*)inputFile->Get("filteredPhoton");

    // Create the output file (it will overwrite the file if it exists)
    TFile* outputFile = new TFile("output.root", "RECREATE");

    // Clone the tree structure (this copies the branches and the structure)
    TTree* outputtree = inputtree->CloneTree(0);

    // std::unique_ptr<TFile> outputFile( TFile::Open("output.root", "UPDATE") );
    // auto outputtree = outputFile->Get<TTree>("filteredPhoton");

    // Branches
    float in_posX;
    float in_posY;
    float in_posZ;
    float in_localPosX;
    float in_localPosY;
    float in_localPosZ;
    float in_edep;

    int in_photonID;
    int in_eventID;
    int in_crystalID;
    int in_rsectorID;
    int in_nCrystalCompton;
    int in_nCrystalRayleigh;

    char in_processName;

    int myin_eventID;

    // Set input branches
    inputtree->SetBranchAddress("posX", &in_posX);
    inputtree->SetBranchAddress("posY", &in_posY);
    inputtree->SetBranchAddress("posZ", &in_posZ);
    inputtree->SetBranchAddress("localPosX", &in_localPosX);
    inputtree->SetBranchAddress("localPosY", &in_localPosY);
    inputtree->SetBranchAddress("localPosZ", &in_localPosZ);
    inputtree->SetBranchAddress("edep", &in_edep);
    inputtree->SetBranchAddress("photonID", &in_photonID);
    inputtree->SetBranchAddress("eventID", &in_eventID);
    inputtree->SetBranchAddress("crystalID", &in_crystalID);
    inputtree->SetBranchAddress("rsectorID", &in_rsectorID);
    inputtree->SetBranchAddress("nCrystalCompton", &in_nCrystalCompton);
    inputtree->SetBranchAddress("nCrystalRayleigh", &in_nCrystalRayleigh);
    inputtree->SetBranchAddress("processName", &in_processName);

    // Counters
    Long64_t iEntry = 0;
    Long64_t iEntry2 = 0;

    bool alreadyFilled = 0;

    Long64_t nEntries = inputtree->GetEntries();

    // std::cout <<nEntries << std::endl;
    inputtree->SetBranchStatus("*",0);
    inputtree->SetBranchStatus("eventID",1);
    inputtree->SetBranchStatus("processName",1);

    while (iEntry < 1000)
    {
        inputtree->GetEntry(iEntry);

        // std::cout << iEntry << std::endl;

        myin_eventID = in_eventID;

        // if(iEntry > 0)
        // {
            if(in_processName=='c')
            {
                inputtree->GetEntry(iEntry+1);

                if(in_processName=='p' && myin_eventID==in_eventID)
                {

                    inputtree->SetBranchStatus("*",1);
                    // inputtree->SetBranchStatus("eventID",1);
                    // inputtree->SetBranchStatus("processName",1);
                    inputtree->GetEntry(iEntry);
                    outputtree->Fill();
                    // outputtree->Write("",TObject::kOverwrite);

                    inputtree->GetEntry(iEntry+1);
                    outputtree->Fill();

                    inputtree->SetBranchStatus("*",0);
                    inputtree->SetBranchStatus("eventID",1);
                    inputtree->SetBranchStatus("processName",1);

                    // std::cout << iEntry << std::endl;

                    iEntry++;
                }
                else
                {
                    bool a=0;
                    while(a==0)
                    {

                        iEntry++;

                        inputtree->GetEntry(iEntry);
                        if (in_eventID!=myin_eventID)
                        {
                            a = 1;
                        }
                    }
                }
            }
            else
            {
                int b=0;

                while(b==0)
                {
                    iEntry++;

                    inputtree->GetEntry(iEntry);
                    if (in_eventID!=myin_eventID)
                    {
                        b = 1;
                    }
                }
            }
    }
        // else
        // {
        //     if(in_processName=='c')
        //     {
        //         inputtree->GetEntry(iEntry+1);
        //     }
        // }
        outputtree->Write("",TObject::kOverwrite);
        outputFile->Close();

        delete outputFile;
}




void createEmptyRootFile() {
    // Apri il file esistente
    TFile *inputFile = new TFile("filteredFinal.root", "READ");
    TTree *tree = (TTree*)inputFile->Get("filteredPhoton"); // Sostituisci con il nome effettivo del TTree

    // Crea il nuovo file
    TFile *outputFile = new TFile("output.root", "RECREATE");
    TTree *newTree = tree->CloneTree(0); // Copia solo la struttura senza dati

    // Salva il nuovo file
    newTree->Write();
    outputFile->Close();
    inputFile->Close();

    delete inputFile;
    delete outputFile;
}



int main()
{
    ROOT::EnableImplicitMT(16); // multi-threading

    // createEmptyRootFile();
    filterTree();

    return 0;
}
