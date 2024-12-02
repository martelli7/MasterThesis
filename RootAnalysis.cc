// $Id: RootAnalysis.cc 1 2024-12-02 martelli $
/**
 * @file   RootAnalysis.cc
 *
 * @date   02 Dec 2024
 * @author martelli
 *
 * @brief  Simple program for opening a previous .root file
 */


//ROOT
#include "TTree.h"
#include "TFile.h"
#include <RtypesCore.h>
#include <TBranch.h>
#include <TObjArray.h>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
//C++
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <variant>


using namespace std;
using DataType = std::variant<int, float, double>;

//Needs the filename as argument
int main (int argc, char** argv)
{

    const char* filename = argv[1];

    TFile *file = TFile::Open(filename, "READ");

    //Check the correct opening
    if(!file || file->IsZombie())
    {
        cerr << "Error. Cannot open file named: " << argv[1] << endl;

        //If cannot open the file, end the program
        return 1;
    }
    //Which file I am opening
    else cout << "Opening file named: " << file->GetName() << endl;

    //Display the file content
    cout << "File content: " << endl;
    file->ls();


    //**** Writing vectors with file's branches ****//

    //Get the pointer to the TTree where hits are stored
    TTree* tree = nullptr;
    tree = (TTree*)file->Get("Hits;-32768");


    //Check the correct opening
    if(!tree)
    {
        cerr << "Error. Cannot open the TTree." << endl;

        //If cannot open the TTree, end the program
        return 1;
    }
    //Which TTree I am opening
    else cout << "Opening TTree named: " << tree->GetName() << endl;

    //Print TTree content
    tree->Print();

    Long64_t nEntries = tree->GetEntries(); //number of entries inside each Branch (=35651444)

    //Retrieve the TTree's branches: look at this example https://root.cern/doc/v632/tree1_8C_source.html
    TObjArray* branches = tree->GetListOfBranches(); //get the branches list
    TBranch* branch = nullptr;

    int nBranches = branches->GetEntries(); //number of branches inside the TTree (=42)

    Int_t PDGEncoding, trackID, parentID,   gantryID, rsectorID, moduleID, submoduleID, crystalID, layerID, photonID, nPhantomCompton, nCrystalCompton, nCrystalRayleigh, primaryID,   sourceID, eventID, runID,   volumeID[10];
    Double_t trackLocalTime, time;
    Float_t edep, stepLength, trackLength, posX, posY, posZ, localPosX, localPosY, localPosZ, momDirX, momDirY, momDirZ,   sourcePosX, sourcePosY, sourcePosZ,   axialPos, rotationAngle;
    Char_t processName, comptVolName, RayleighVolName;

    tree->SetBranchAddress("PDGEncoding", &PDGEncoding);
    tree->SetBranchAddress("trackID", &trackID);
    tree->SetBranchAddress("parentID", &parentID);

    tree->SetBranchAddress("trackLocalTime", &trackLocalTime);
    tree->SetBranchAddress("time", &time);

    tree->SetBranchAddress("edep", &edep);
    tree->SetBranchAddress("stepLength", &stepLength);
    tree->SetBranchAddress("trackLength", &trackLength);
    tree->SetBranchAddress("posX", &posX);
    tree->SetBranchAddress("posY", &posY);
    tree->SetBranchAddress("posZ", &posZ);
    tree->SetBranchAddress("localPosX", &localPosX);
    tree->SetBranchAddress("localPosY", &localPosY);
    tree->SetBranchAddress("localPosZ", &localPosZ);
    tree->SetBranchAddress("momDirX", &momDirX);
    tree->SetBranchAddress("momDirY", &momDirY);
    tree->SetBranchAddress("momDirZ", &momDirZ);

    tree->SetBranchAddress("gantryID", &gantryID);
    tree->SetBranchAddress("rsectorID", &rsectorID);
    tree->SetBranchAddress("moduleID", &moduleID);
    tree->SetBranchAddress("submoduleID", &submoduleID);
    tree->SetBranchAddress("crystalID", &crystalID);
    tree->SetBranchAddress("layerID", &layerID);
    tree->SetBranchAddress("photonID", &photonID);
    tree->SetBranchAddress("nPhantomCompton", &nPhantomCompton);
    tree->SetBranchAddress("nCrystalCompton", &nCrystalCompton);
    tree->SetBranchAddress("nCrystalRayleigh", &nCrystalRayleigh);
    tree->SetBranchAddress("primaryID", &primaryID);

    tree->SetBranchAddress("sourcePosX", &sourcePosX);
    tree->SetBranchAddress("sourcePosY", &sourcePosY);
    tree->SetBranchAddress("sourcePosZ", &sourcePosZ);

    tree->SetBranchAddress("sourceID", &sourceID);
    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("runID", &runID);

    tree->SetBranchAddress("axialPos", &axialPos);
    tree->SetBranchAddress("rotationAngle", &rotationAngle);

    tree->SetBranchAddress("volumeID", volumeID);

    tree->SetBranchAddress("processName", &processName);
    tree->SetBranchAddress("comptVolName", &comptVolName);
    tree->SetBranchAddress("RayleighVolName", &RayleighVolName);


    //Retrieve the events!
    for (Long64_t i = 0; i < nBranches; i++)
    {
        tree->GetEntry(i); //it wants a Long64_t, returns the entries of all the branches of the i-th event
    }
    tree->ResetBranchAddresses();

    //Debug stuff
    cout << nBranches << endl;

    //Closing and deleting
    file->Close();
    delete file;
    file = 0;

    if(file == 0) cout << "File closed correctly." << endl;
    else cerr << "Error. Cannot close the file." << endl;

    return 0;
}









