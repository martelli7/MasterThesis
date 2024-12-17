// $Id: RootAnalysis.cc 2 2024-12-12 martelli $
/**
 * @file   RootAnalysis.cc
 *
 * @date   12 Dec 2024
 * @author martelli
 *
 * @brief  Simple program for opening a previous .root file
 */


//ROOT
#include "TTree.h"
#include "TFile.h"
#include <ROOT/RVec.hxx>
#include <RtypesCore.h>
#include <TBranch.h>
#include <TMathBase.h>
#include <TObjArray.h>
#include "TH1.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
//C++
#include <iostream>
#include <vector>
#include <string>
#include <map>
//#include <variant>


// commands for ROOT terminal:
// TTree *tree = (TTree*)_file0->Get("Hits;-32768");
// tree->Draw("parentID>>(200,0,200)");
// tree->Draw("posX:posY","nCrystalCompton==1 & eventID==137","colz");

/*const std::vector<Float_t>& posX,
 * const std::vector<Float_t>& posY,
 * const std::vector<Float_t>& posZ) {
 *
 * root [5] tree->Scan("eventID:posX:posY:posZ","eventID==99")
 ************************************************************
 *    Row   *   eventID *      posX *      posY *      posZ *
 ************************************************************
 *   121232 *        99 * -45.14472 * 105.52783 * -4.894348 *
 *   121233 *        99 * -51.04073 * 98.252044 * -3.312198 *
 *   121234 *        99 * 40.279876 * -94.15602 * 4.3669276 *
 ************************************************************
 *
 *
 */


int main()
{
    ROOT::EnableImplicitMT(); // multi-threading

    auto fileName = "Tesi_NN_0.root";
    auto treeName = "Hits;-32768";

    // apriamo il RootDataFrame
    ROOT::RDataFrame df(treeName, fileName);



    // Filtering
    auto df_filtered = df.Filter("nCrystalCompton == 1", "Events with just 1 Compton Scattering")
    .Filter("parentID==0","Photons which are not daughters") //domanda: gli eventi con parentID = 2,3,4 cosa sono...?
    .Filter("posX!=0 && posY!=0 && posZ!=0", "Events with photons")
    .Filter("photonID==1","Select the first photon only")
    //.Filter()   //possibile concatenare più filtri
    //.Count();   //Count() non funziona usando anche Snapshot()
    .Snapshot("filteredHits", "filtered_file_photon1.root", {"posX", "posY", "posZ", "photonID", "eventID", "parentID", "nCrystalCompton"});
    //il metodo Snapshot() crea un nuovo file con i filtri e le colonne selezionate, ora il file è molto leggero

    //controlla quante entrate hanno superato il filto
    //std::cout << *df_filtered << " entries passed all filters" << std::endl;

    // Histo2D posX:posY -> funziona!
    /*auto myHist1 = df.Histo2D({"histName", "posX:posY{nCrystalCompton==1}", 64u, -150., 150., 32u, -150., 150.}, "posX", "posY");
     *
     *    TCanvas c;
     *    myHist1->Draw("COLZ"); //scala a colori
     *    c.SaveAs("posXposY.root");
     */

    //----------------------------------------------------------//
    // NOTA: questa parte non funziona, ma secondo me l'idea è utilizzabile
    // mi dà un errore di memoria a riga 113 e 124 quando cerco di accedere agli array... forse il problema è risolvibile semplicemente dichiarandoli di grandezza nota e fine
    // ora apriamo il nuovo file root e analizziamo le posizioni


    std::unique_ptr<TFile> myFile( TFile::Open("filtered_file_photon1.root") );
    auto tree = myFile->Get<TTree>("filteredHits");

    //std::vector<Float_t> posX, posY, posZ;
    //Float_t posX[10], posY[10], posZ[10];
    Float_t *posX, *posY, *posZ;
    Int_t *photonID, *eventID;
    Int_t *row;
    tree->SetBranchAddress("posX", &posX);
    tree->SetBranchAddress("posY", &posY);
    tree->SetBranchAddress("posZ", &posZ);
    tree->SetBranchAddress("photonID", &photonID);
    tree->SetBranchAddress("eventID", &eventID);

    for (Long64_t iEntry = 0; tree->LoadTree(iEntry) >= 0; ++iEntry) {
        // Load the data for the given tree entry
        tree->GetEntry(iEntry);
        row[iEntry] = iEntry;
    }

    Long64_t i = 0;
    Float_t dx1, dy1, dz1, dx2, dy2, dz2;
    Float_t distance;

    TH1D *hDistanceScattering = new TH1D("hDistanceScattering", "Scattering Distance;Distance [units]];Counts", 100, 0, 50);

    while(i < sizeof(row)) {
        // Retrieve eventID entries to check the same-event-photons
        if(eventID[i] == eventID[i+1] && eventID[i+1] == eventID[i+2])
        {
            // This event has a scattering event!
            // have to check which couple (i, i+1 or i+1, i+2) is actually a scattering

            // First couple
            dx1 = TMath::Abs(posX[i] - posX[i+1]);
            dy1 = TMath::Abs(posY[i] - posY[i+1]);
            dz1 = TMath::Abs(posZ[i] - posZ[i+1]);

            // Second couple
            dx2 = TMath::Abs(posX[i+1] - posX[i+2]);
            dy2 = TMath::Abs(posY[i+1] - posY[i+2]);
            dz2 = TMath::Abs(posZ[i+1] - posZ[i+2]);

            // First couple condition
            if(dx1 < 30)
            {
                distance = TMath::Sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
                hDistanceScattering->Fill(distance);
            }
            // Second couple condition
            else
            {
                distance = TMath::Sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
                hDistanceScattering->Fill(distance);
            }
            // Now we skip to the following two events (already processed)
            i = i+3;
        }
        else
        {
            i++; // Scattering not found, we keep scrolling each eventID
        }
    }

    TCanvas c1;
    hDistanceScattering->Draw();
    c1.SaveAs("distanceScatteringPhoton1.root");

    delete [] eventID;
    delete [] posX;
    delete [] posY;
    delete [] posZ;
    delete [] photonID;
    delete hDistanceScattering;


    return 0;
}









