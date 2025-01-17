// $Id: RootAnalysis.cc 4 2025-1-16 martelli $
/**
 * @file   RootAnalysis.cc
 *
 * @date   16 Jan 2025
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

    // auto fileName = "Tesi_NN_0.root";
    // auto treeName = "Hits;-32768";
    //
    // // apriamo il RootDataFrame
    // ROOT::RDataFrame df(treeName, fileName);
    //
    //
    //
    // // Filtering
    // auto df_filtered = df.Filter("nCrystalCompton == 1", "Events with just 1 Compton Scattering")
    // .Filter("parentID==0","Photons which are not daughters")
    // .Filter("posX!=0 && posY!=0 && posZ!=0", "Events with photons")
    // //.Filter("photonID==1","Select the first photon only")
    // //.Filter()   //possibile concatenare più filtri
    // //.Count();   //Count() non funziona usando anche Snapshot()
    // .Snapshot("filteredHits", "filtered_file2.root", {"posX", "posY", "posZ", "edep", "photonID", "eventID", "parentID", "crystalID", "nCrystalCompton"});
    // //il metodo Snapshot() crea un nuovo file con i filtri e le colonne selezionate, ora il file è molto leggero

    //controlla quante entrate hanno superato il filto
    //std::cout << *df_filtered << " entries passed all filters" << std::endl;

    // Histo2D posX:posY -> funziona!
    //auto myHist1 = df.Histo2D({"histName", "posX:posY{nCrystalCompton==1}", 64u, -150., 150., 32u, -150., 150.}, "posX", "posY");
    /*
     *    TCanvas c;
     *    myHist1->Draw("COLZ"); //scala a colori
     *    c.SaveAs("posXposY.root");
     */

    //----------------------------------------------------------//

    std::unique_ptr<TFile> myFile( TFile::Open("filtered_file2.root") );
    auto tree = myFile->Get<TTree>("filteredHits");

    Long64_t nEntries = tree->GetEntries(); //number of entries inside each Branch (=9677397)

    // Branches
    Float_t posX;
    Float_t posY;
    Float_t posZ;
    Float_t edep;
    Int_t photonID;
    Int_t eventID;
    Int_t crystalID;

    // My arrays
    Float_t *myposX = new Float_t[3];
    Float_t *myposY = new Float_t[3];
    Float_t *myposZ = new Float_t[3];
    Float_t *myedep = new Float_t[3];
    Int_t *myphotonID = new Int_t[3];
    Int_t *myeventID = new Int_t[3];;
    Int_t *mycrystalID = new Int_t[3];

    //tree->SetBranchStatus("*",0); //disable all branches

    tree->SetBranchAddress("posX", &posX);
    tree->SetBranchAddress("posY", &posY);
    tree->SetBranchAddress("posZ", &posZ);
    tree->SetBranchAddress("edep", &edep);
    tree->SetBranchAddress("photonID", &photonID);
    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("crystalID", &crystalID);

    Long64_t iEntry = 0;
    Long64_t iEntry2 = 0;
    Long64_t iEntry3 = 0;

    Float_t dx1, dy1, dz1, dx2, dy2, dz2;
    Int_t dxc1, dxc2, dyc1, dyc2;
    Float_t distance, crydistance;

    // Histograms
    TH1D *hDistanceScattering = new TH1D("hDistanceScattering", "Scattering Distance;Distance [mm];Counts", 1000, 0, 30);
    TH1D *hDistanceCrystal = new TH1D("hDistanceCrystal", "Crystal Scattering Distance;Distance [crystals];Counts", 1000, 0, 30);
    TH1D *hDistanceZ = new TH1D("hDistanceZ", "Scattering Distance along z-axis;Distance [mm];Counts", 1000, 0, 30);

    TH2D *hScatterEnergy = new TH2D("hScatterEnergy", "Distance[mm] vs 1st Photon Energy[keV]", 80, 0, 40, 17, 0, 511);

    //for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry)
    while (iEntry < nEntries)
    {
        // Load the data for the given tree entry
        tree->GetEntry(iEntry);

        if ( crystalID == 119 ) // Select events on the central crystal only
        {
            // Transfer the current leaf to my arrays
            myposX[0] = posX;
            myposY[0] = posY;
            myposZ[0] = posZ;
            myedep[0] = edep*1000; // Mev-->keV
            myphotonID[0] = photonID;
            myeventID[0] = eventID;
            mycrystalID[0] = crystalID;

            iEntry2 = iEntry+1;
            iEntry3 = iEntry+2;

            // Call the second leaf
            tree->GetEntry(iEntry2);

            myposX[1] = posX;
            myposY[1] = posY;
            myposZ[1] = posZ;
            myedep[1] = edep*1000;
            myphotonID[1] = photonID;
            myeventID[1] = eventID;
            mycrystalID[1] = crystalID;

            // Call the third leaf
            tree->GetEntry(iEntry3);

            myposX[2] = posX;
            myposY[2] = posY;
            myposZ[2] = posZ;
            myedep[2] = edep*1000;
            myphotonID[2] = photonID;
            myeventID[2] = eventID;
            mycrystalID[2] = crystalID;

            // Retrieve eventID entries to check the same-event-photons
            if( myeventID[0] == myeventID[1] && myeventID[1] == myeventID[2] )
            {
                // This event has a scattering event!
                // have to check which couple (0,1 or 1,2) is actually a scattering

                // First couple
                dx1 = TMath::Abs( myposX[0] - myposX[1]);
                dy1 = TMath::Abs( myposY[0] - myposY[1]);
                dz1 = TMath::Abs( myposZ[0] - myposZ[1]);

                // Second couple
                dx2 = TMath::Abs( myposX[1] - myposX[2]);
                dy2 = TMath::Abs( myposY[1] - myposY[2]);
                dz2 = TMath::Abs( myposZ[1] - myposZ[2]);

                // Compute crystal distance...
                // First couple
                dxc1 = TMath::Abs(mycrystalID[0]/16 - mycrystalID[1]/16);
                dyc1 = TMath::Abs(mycrystalID[0]%16 - mycrystalID[1]%16);

                // Second couple
                dxc2 = TMath::Abs(mycrystalID[1]/16 - mycrystalID[2]/16);
                dyc2 = TMath::Abs(mycrystalID[1]%16 - mycrystalID[2]%16);

                // First couple condition
                if(dx1 < 30 && dz1 > 3)
                {
                    distance = TMath::Sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);

                    crydistance = TMath::Sqrt(dxc1*dxc1 + dyc1*dyc1);

                    // Display some outputs...
                    //std::cout << myphotonID[i] << ", " << i << ", "<< myeventID[i] << ", (" << dx1 << "," << dy1 << "," << dz1 << ") " << "(" << dxc1 << "," << dyc1 << ")" <<  std::endl;

                    // Filling histograms
                    hDistanceScattering->Fill(distance);
                    hDistanceCrystal->Fill(crydistance);
                    hDistanceZ->Fill(dz1);
                    hScatterEnergy->Fill(distance, myedep[0]);
                }
                // Second couple condition
                else if(dz2 > 3)
                {
                    distance = TMath::Sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);

                    crydistance = TMath::Sqrt(dxc2*dxc2 + dyc2*dyc2);

                    // Display some outputs...
                    //std::cout << myphotonID[i] << ", " << i << ", "<< myeventID[i] << ", (" << dx2 << "," << dy2 << "," << dz2 << ") " << "(" << dxc1 << "," << dyc1 << ")" <<  std::endl;

                    // Filling histograms
                    hDistanceScattering->Fill(distance);
                    hDistanceCrystal->Fill(crydistance);
                    hDistanceZ->Fill(dz2);
                    hScatterEnergy->Fill(distance, myedep[1]);
                }
                // Now we skip the following 4th event (3 events processed in one time)
                iEntry = iEntry+3;
            }
            else if(myeventID[0] == myeventID[1] && myeventID[1] != myeventID[2])
            {
                // This is a scattering + second photon is lost

                // First couple
                dx1 = TMath::Abs( myposX[0] - myposX[1]);
                dy1 = TMath::Abs( myposY[0] - myposY[1]);
                dz1 = TMath::Abs( myposZ[0] - myposZ[1]);

                // Compute crystal distance...
                // First couple
                dxc1 = TMath::Abs(mycrystalID[0]/16 - mycrystalID[1]/16);
                dyc1 = TMath::Abs(mycrystalID[0]%16 - mycrystalID[1]%16);

                crydistance = TMath::Sqrt(dxc1*dxc1 + dyc1*dyc1);

                // First couple condition
                if(dx1 < 30 && dz1 > 3)
                {
                    distance = TMath::Sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);

                    // Display some outputs...
                    //std::cout << myphotonID[i] << ", " << i << ", "<< myeventID[i] << ", (" << dx1 << "," << dy1 << "," << dz1 << ")"<< std::endl;

                    // Filling histograms
                    hDistanceScattering->Fill(distance);
                    hDistanceCrystal->Fill(crydistance);
                    hDistanceZ->Fill(dz1);
                    hScatterEnergy->Fill(distance, myedep[0]);
                }
                // Now we skip to the 3rd event (2 events processed in one time)
                iEntry = iEntry+2;
            }
            else if(myeventID[0] != myeventID[1] && myeventID[1] == myeventID[2])
            {
                // This is a scattering + second photon is lost

                // Second couple
                dx2 = TMath::Abs( myposX[1] - myposX[2]);
                dy2 = TMath::Abs( myposY[1] - myposY[2]);
                dz2 = TMath::Abs( myposZ[1] - myposZ[2]);

                // Second couple
                dxc2 = TMath::Abs(mycrystalID[1]/16 - mycrystalID[2]/16);
                dyc2 = TMath::Abs(mycrystalID[1]%16 - mycrystalID[2]%16);

                if (dx2 < 30 && dz2 > 3)
                {
                    distance = TMath::Sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);

                    crydistance = TMath::Sqrt(dxc2*dxc2 + dyc2*dyc2);

                    // Display some outputs...
                    //std::cout << myphotonID[i] << ", " << i << ", "<< myeventID[i] << ", (" << dx2 << "," << dy2 << "," << dz2 << ")"<< std::endl;

                    // Filling histograms
                    hDistanceScattering->Fill(distance);
                    hDistanceCrystal->Fill(crydistance);
                    hDistanceZ->Fill(dz2);
                    hScatterEnergy->Fill(distance, myedep[1]);
                }
                // Now we skip to the 3rd event (2 events processed in one time)
                iEntry = iEntry+2;
            }
            else
            {
                iEntry++; // Scattering not found, we keep scrolling each eventID (just 1 event processed)
            }
        } //close if crystalID
        else
        {
            iEntry++; // The current event does not take place in the central crystal, we keep scrolling each eventID
        }
    } //close while

    // TCanvas c1;
    // hDistanceScattering->Draw();
    // c1.SaveAs("distanceScatteringPhotons2.root");
    // c1.SaveAs("distanceScatteringPhotons2.png");
    //
    // TCanvas c2;
    // hDistanceCrystal->Draw();
    // c2.SaveAs("distanceCrystalScattering2.root");
    // c2.SaveAs("distanceCrystalScattering2.png");
    //
    // TCanvas c3;
    // hDistanceZ->Draw();
    // c3.SaveAs("distanceScatteringZ2.root");
    // c3.SaveAs("distanceScatteringZ2.png");

    TCanvas c4;
    hScatterEnergy->Draw("colz");
    c4.SaveAs("scatterEnergy.root");
    c4.SaveAs("scatterEnergy.png");

    return 0;
}









