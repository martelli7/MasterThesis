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




void checkRootFile(std::string inputFile)
{
    // Apri il file esistente
    TFile *checkFile = new TFile(inputFile.c_str(), "READ");
    TTree *checktree = (TTree*)checkFile->Get("filteredPhoton"); // Sostituisci con il nome effettivo del TTree

    // Branches
    float in_posX;
    float in_posY;
    float in_posZ;
    float in_localPosX;
    float in_localPosY;
    float in_localPosZ;
    float in_edep; //MeV

    int in_photonID;
    int in_eventID;
    int in_crystalID;
    int in_rsectorID;
    int in_nCrystalCompton;
    int in_nCrystalRayleigh;

    char in_processName[10];

    // Set input branches
    checktree->SetBranchAddress("posX", &in_posX);
    checktree->SetBranchAddress("posY", &in_posY);
    checktree->SetBranchAddress("posZ", &in_posZ);
    checktree->SetBranchAddress("localPosX", &in_localPosX);
    checktree->SetBranchAddress("localPosY", &in_localPosY);
    checktree->SetBranchAddress("localPosZ", &in_localPosZ);
    checktree->SetBranchAddress("edep", &in_edep);
    checktree->SetBranchAddress("photonID", &in_photonID);
    checktree->SetBranchAddress("eventID", &in_eventID);
    checktree->SetBranchAddress("crystalID", &in_crystalID);
    checktree->SetBranchAddress("rsectorID", &in_rsectorID);
    checktree->SetBranchAddress("nCrystalCompton", &in_nCrystalCompton);
    checktree->SetBranchAddress("nCrystalRayleigh", &in_nCrystalRayleigh);
    checktree->SetBranchAddress("processName", &in_processName);

    int myin_eventID;
    char myin_processName[10];
    float myin_edep;

    Long64_t nEntries = checktree->GetEntries(); //10956826
    Long64_t iEntry = 0;

    // Enable branches
    checktree->SetBranchStatus("*",0);
    checktree->SetBranchStatus("eventID",1);
    checktree->SetBranchStatus("processName",1);
    checktree->SetBranchStatus("edep",1);

    Long64_t goodEv = 0;
    Long64_t bad1Ev = 0;
    Long64_t bad2Ev = 0;
    Long64_t badEdep = 0;

    // while(iEntry < nEntries)
    // {
    //     checktree->GetEntry(iEntry);
    //
    //     myin_eventID = in_eventID;
    //     std::copy(std::begin(in_processName), std::end(in_processName), std::begin(myin_processName));
    //
    //     if(myin_processName[0]=='c') // First hit->Compton
    //     {
    //         if(in_edep > 0.34066666)
    //         {
    //             // std::cout << "Bad Energy occurred: " << std::endl;
    //             // std::cout << in_processName << ", " << myin_eventID << ", " << in_edep << ", " << iEntry << std::endl;
    //             badEdep++;
    //         }
    //
    //         checktree->GetEntry(iEntry+1);
    //
    //         if(in_processName[0]=='p' && myin_eventID==in_eventID) // Second hit->Photoelectric
    //         {
    //             goodEv++;
    //         }
    //         else
    //         {
    //             std::cout << "2 Bad Events occurred: " << std::endl;
    //             std::cout << myin_processName << ", " << myin_eventID << ", " << iEntry << std::endl;
    //             std::cout << in_processName << ", " << in_eventID << ", " << iEntry+1 << std::endl;
    //             bad2Ev++;
    //         }
    //
    //         iEntry=iEntry+2;
    //     }
    //     else
    //     {
    //         std::cout << "1 Bad Event occurred: " << std::endl;
    //         std::cout << myin_processName << ", " << myin_eventID << ", " << iEntry << std::endl;
    //         bad1Ev++;
    //         iEntry++;
    //     }
    // }

    while(iEntry < nEntries)
    {
        checktree->GetEntry(iEntry);

        myin_eventID = in_eventID;
        myin_edep = in_edep;
        if(myin_edep>0.3406)
        {
            std::copy(std::begin(in_processName), std::end(in_processName), std::begin(myin_processName));

            checktree->GetEntry(iEntry+1);

            float tot_energy = in_edep+myin_edep;

            std::cout << myin_edep << " + " << in_edep << " = " << tot_energy << std::endl;
        }


        iEntry=iEntry+2;
    }

     std::cout << "Good Events = " << goodEv << ", including Bad Edep = " << badEdep <<  std::endl;
     std::cout << "Bad Single Events = " << bad1Ev << std::endl;
     std::cout << "Bad Double Events = " << bad2Ev << std::endl;


    // Salva il nuovo file
    checkFile->Close();

    delete checkFile;
}


int main(int argc, char* argv[])
{
    ROOT::EnableImplicitMT(16); // multi-threading

    std::string inputFile = argv[1];

    checkRootFile(inputFile);

    return 0;
}
