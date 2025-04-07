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

void filterTree(const char* inputF, const char* outputF)
{
    TFile *inputFile = new TFile(inputF, "READ");
    TTree *inputtree = (TTree*)inputFile->Get("filteredPhoton");

    // Create the output file (it will overwrite the file if it exists)
    TFile* outputFile = new TFile(outputF, "RECREATE");

    // Clone the tree structure (this copies the branches and the structure)
    TTree* outputtree = inputtree->CloneTree(0);

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

    char in_processName[15];

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

    Long64_t iEntry = 0;
    Long64_t nEntries = inputtree->GetEntries(); // 27561311

    //Enable branches
    inputtree->SetBranchStatus("*",0);
    inputtree->SetBranchStatus("eventID",1);
    inputtree->SetBranchStatus("processName",1);

    while (iEntry < nEntries)
    {
        inputtree->GetEntry(iEntry);

        myin_eventID = in_eventID;

        if(iEntry==122 || iEntry==125)
        {
            std::cout << "Prima dei cicli: " << in_processName << ", " <<in_eventID << std::endl;
        }

        if(iEntry==nEntries-1)
        {
            std::cout << "Ultimo valore: " << in_processName << ", " <<in_eventID << std::endl;
        }

        if(in_processName[0]=='c')
        {
            if(iEntry==122 || iEntry==125 || iEntry==nEntries-1)
            {
                std::cout << "dentro if \'c\' the first hit: " << in_processName << ", " <<in_eventID << std::endl;
            }

            inputtree->GetEntry(iEntry+1);

            if(in_processName[0]=='p' && myin_eventID==in_eventID)
            {

                if(iEntry==122 || iEntry==125|| iEntry==nEntries-1)
                {
                    std::cout << "dentro if \'p\' the second hit: " << in_processName << ", " <<in_eventID << std::endl;
                }

                inputtree->SetBranchStatus("*",1);

                inputtree->GetEntry(iEntry);
                outputtree->Fill();

                inputtree->GetEntry(iEntry+1);
                outputtree->Fill();

                inputtree->SetBranchStatus("*",0);
                inputtree->SetBranchStatus("eventID",1);
                inputtree->SetBranchStatus("processName",1);

                iEntry++;
            }
            else
            {
                if(iEntry==122 || iEntry==125 || iEntry==nEntries-1)
                {
                    std::cout << "dentro else (not \'p\' the second hit): " << in_processName << ", " <<in_eventID << std::endl;
                }

                bool a=0;
                while(a==0)
                {
                    iEntry++;
                    if(iEntry >= nEntries)
                    {
                        break;
                    }

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
            int c = 0;
            if(iEntry==122 || iEntry==125 || iEntry==nEntries-1)
            {
                std::cout << "dentro else (not \'c\' the first hit): " << in_processName << ", " <<in_eventID << std::endl;
                c = 1;
            }

            int b=0;

            while(b==0)
            {
                iEntry++;
                if(iEntry >= nEntries)
                {
                    break;
                }

                inputtree->GetEntry(iEntry);

                if (in_eventID!=myin_eventID)
                {
                    b = 1;
                }
            }

            if(c==1)
            {
                std::cout << "dopo while->prossimo evento: " << in_processName << ", " <<in_eventID << ", " << iEntry << std::endl;
                c=0;
            }
        }
    }

    std::cout << "Fine" << std::endl;
    outputtree->Write("",TObject::kOverwrite);
    outputFile->Close();
    inputFile->Close();

    delete outputFile;
    delete inputFile;


}



int main(int argc, char* argv[])
{
    ROOT::EnableImplicitMT(16); // multi-threading

    const char* inputF = argv[1];
    const char* outputF = argv[2];

    filterTree(inputF, outputF);

    return 0;
}
