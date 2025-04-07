#include <ROOT/RDataFrame.hxx>
#include "TFile.h"
#include <RtypesCore.h>
#include <TMath.h>
#include <cmath>
#include <cstring>


void filterDataFrame(const char* fileName, const char* treeName, const char* fileName1, const char* fileName2)
{
    // Open the RootDataFrame
    ROOT::RDataFrame df(treeName, fileName);

    // Filtering
    auto df_filtered1 = df.Filter("photonID==1", "Select the first Photon only")
        .Filter("nCrystalCompton < 2","Select rows with only one Compton")
        //.Count();   //Count() does not work with Snapshot()
        .Snapshot("filteredPhoton", "filtered_Phys4_Photon1.root", {"posX", "posY", "posZ", "edep", "localPosX", "localPosY", "localPosZ", "photonID", "eventID", "parentID", "crystalID", "rsectorID", "nCrystalCompton", "nCrystalRayleigh", "processName"});

    auto df_filtered2 = df.Filter("photonID==2", "Select the first Photon only")
        .Filter("nCrystalCompton < 2","Select rows with only one Compton")
        //.Count();   //Count() does not work with Snapshot()
        .Snapshot("filteredPhoton", "filtered_Phys4_Photon2.root", {"posX", "posY", "posZ", "edep", "localPosX", "localPosY", "localPosZ", "photonID", "eventID", "parentID", "crystalID", "rsectorID", "nCrystalCompton", "nCrystalRayleigh", "processName"});
}


int main(int argc, char* argv[])
{
    ROOT::EnableImplicitMT(); // multi-threading

    const char* fileName = argv[1];
    const char* fileName1 = argv[2];
    const char* fileName2 = argv[3];
    const char* treeName = "Hits";

    filterDataFrame(fileName,treeName,fileName1,fileName2);

    return 0;
}
