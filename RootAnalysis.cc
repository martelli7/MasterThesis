// $Id: RootAnalysis.cc 4 2025-1-24 martelli $
/**
 * @file   RootAnalysis.cc
 *
 * @date   24 Jan 2025
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
#include <cmath>
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

// Physics constants
const Double_t h = TMath::H();     // Planck [J*s]
const Double_t c = TMath::C();     // light speed [m/s]
const Double_t me = 9.10938356e-31;    // m electron [kg]
const Double_t Qe = TMath::Qe(); // q electron [eV->J]
const Double_t en511 = 511*1000*Qe; // [J]
const Double_t lambdaIn = h*c/en511; // lambda of 511keV [m]

// this function takes in the energy deposit by the scattered 511 and returns its scatter angle by reversing Compton formula
Float_t scatterAngle (Double_t edep)
{
    Double_t enGamma2 = (511*1000 - edep*1000)*Qe;
    Double_t lambdaOut = h*c/enGamma2;

    Double_t angleR = TMath::ACos( 1 + (lambdaIn - lambdaOut)*me*c/h  ); // Rad

    Double_t angleD = angleR*TMath::RadToDeg(); // Deg

    return angleD;
}



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
    // .Snapshot("filteredHits", "root/filtered_file2.root", {"posX", "posY", "posZ", "edep", "photonID", "eventID", "parentID", "crystalID", "nCrystalCompton"});
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
    Float_t distance, distanceXY, crydistance, sAngle;

    // Histograms
    TH1D *hDistanceScattering = new TH1D("hDistanceScattering", "Scattering Distance;Distance [mm];Counts", 1000, 0, 30);
    TH1D *hDistanceCrystal = new TH1D("hDistanceCrystal", "Crystal Scattering Distance;Distance [crystals];Counts", 1000, 0, 30);
    TH1D *hDistanceZ = new TH1D("hDistanceZ", "Scattering Distance along z-axis;Distance [mm];Counts", 1000, 0, 30);
    TH1D *hScatterAngle = new TH1D("hScatterAngle", "Scattering Angle;Angle[#degree];Counts", 380, -5,180);

    TH2D *hScatterEnergy0 = new TH2D("hScatterEnergy0", "#Deltaz = [0#;3)mm;xy-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
    TH2D *hScatterEnergy1 = new TH2D("hScatterEnergy1", "#Deltaz = [3#;6)mm;xy-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
    TH2D *hScatterEnergy2 = new TH2D("hScatterEnergy2", "#Deltaz = [6#;9)mm;xy-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
    TH2D *hScatterEnergy3 = new TH2D("hScatterEnergy3", "#Deltaz = [9#;12)mm;xy-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
    TH2D *hScatterEnergy4 = new TH2D("hScatterEnergy4", "#Deltaz = [12#;15)mm;xy-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);

    //TH2D *hScatterEnergyZ = new TH2D("hScatterEnergyZ", "#DeltaZ[mm] vs 1st Photon Energy[keV]", 5, 0, 15, 17, 0, 511);

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
                if(dx1 < 30)
                {
                    distance = TMath::Sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
                    distanceXY = TMath::Sqrt(dx1 * dx1 + dy1 * dy1);

                    crydistance = TMath::Sqrt(dxc1*dxc1 + dyc1*dyc1);

                    sAngle = scatterAngle(myedep[0]);

                    // Display some outputs...
                    //std::cout << myphotonID[i] << ", " << i << ", "<< myeventID[i] << ", (" << dx1 << "," << dy1 << "," << dz1 << ") " << "(" << dxc1 << "," << dyc1 << ")" <<  std::endl;

                    // Filling histograms
                    hDistanceScattering->Fill(distance);
                    hDistanceCrystal->Fill(crydistance);
                    hDistanceZ->Fill(dz1);
                    hScatterAngle->Fill(sAngle);

                    if (dz1 >= 0 && dz1 < 3)
                    {
                        hScatterEnergy0->Fill(distanceXY, myedep[0]);
                    }
                    else if(dz1 >= 3 && dz1 < 6)
                    {
                        hScatterEnergy1->Fill(distanceXY, myedep[0]);
                    }
                    else if (dz1 >= 6 && dz1 < 9)
                    {
                        hScatterEnergy2->Fill(distanceXY, myedep[0]);
                    }
                    else if (dz1 >= 9 && dz1 < 12)
                    {
                        hScatterEnergy3->Fill(distanceXY, myedep[0]);
                    }
                    else if (dz1 >= 12 && dz1 < 15)
                    {
                        hScatterEnergy4->Fill(distanceXY, myedep[0]);
                    }
                }
                // Second couple condition
                else
                {
                    distance = TMath::Sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
                    distanceXY = TMath::Sqrt(dx2 * dx2 + dy2 * dy2);

                    crydistance = TMath::Sqrt(dxc2*dxc2 + dyc2*dyc2);

                    sAngle = scatterAngle(myedep[1]);

                    // Display some outputs...
                    //std::cout << myphotonID[i] << ", " << i << ", "<< myeventID[i] << ", (" << dx2 << "," << dy2 << "," << dz2 << ") " << "(" << dxc1 << "," << dyc1 << ")" <<  std::endl;

                    // Filling histograms
                    hDistanceScattering->Fill(distance);
                    hDistanceCrystal->Fill(crydistance);
                    hDistanceZ->Fill(dz2);
                    hScatterAngle->Fill(sAngle);

                    if (dz1 >= 0 && dz1 < 3)
                    {
                        hScatterEnergy0->Fill(distanceXY, myedep[1]);
                    }
                    else if(dz1 >= 3 && dz1 < 6)
                    {
                        hScatterEnergy1->Fill(distanceXY, myedep[1]);
                    }
                    else if (dz2 >= 6 && dz2 < 9)
                    {
                        hScatterEnergy2->Fill(distanceXY, myedep[1]);
                    }
                    else if (dz2 >= 9 && dz2 < 12)
                    {
                        hScatterEnergy3->Fill(distanceXY, myedep[1]);
                    }
                    else if (dz2 >= 12 && dz2 < 15)
                    {
                        hScatterEnergy4->Fill(distanceXY, myedep[1]);
                    }
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
                if(dx1 < 30)
                {
                    distance = TMath::Sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
                    distanceXY = TMath::Sqrt(dx1 * dx1 + dy1 * dy1);

                    sAngle = scatterAngle(myedep[0]);

                    // Display some outputs...
                    //std::cout << myphotonID[i] << ", " << i << ", "<< myeventID[i] << ", (" << dx1 << "," << dy1 << "," << dz1 << ")"<< std::endl;

                    // Filling histograms
                    hDistanceScattering->Fill(distance);
                    hDistanceCrystal->Fill(crydistance);
                    hDistanceZ->Fill(dz1);
                    hScatterAngle->Fill(sAngle);

                    if (dz1 >= 0 && dz1 < 3)
                    {
                        hScatterEnergy0->Fill(distanceXY, myedep[0]);
                    }
                    else if(dz1 >= 3 && dz1 < 6)
                    {
                        hScatterEnergy1->Fill(distanceXY, myedep[0]);
                    }
                    else if (dz1 >= 6 && dz1 < 9)
                    {
                        hScatterEnergy2->Fill(distanceXY, myedep[0]);
                    }
                    else if (dz1 >= 9 && dz1 < 12)
                    {
                        hScatterEnergy3->Fill(distanceXY, myedep[0]);
                    }
                    else if (dz1 >= 12 && dz1 < 15)
                    {
                        hScatterEnergy4->Fill(distanceXY, myedep[0]);
                    }
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

                if (dx2 < 30)
                {
                    distance = TMath::Sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
                    distanceXY = TMath::Sqrt(dx2 * dx2 + dy2 * dy2);

                    crydistance = TMath::Sqrt(dxc2*dxc2 + dyc2*dyc2);

                    sAngle = scatterAngle(myedep[1]);

                    // Display some outputs...
                    //std::cout << myphotonID[i] << ", " << i << ", "<< myeventID[i] << ", (" << dx2 << "," << dy2 << "," << dz2 << ")"<< std::endl;

                    // Filling histograms
                    hDistanceScattering->Fill(distance);
                    hDistanceCrystal->Fill(crydistance);
                    hDistanceZ->Fill(dz2);
                    hScatterAngle->Fill(sAngle);

                    if (dz1 >= 0 && dz1 < 3)
                    {
                        hScatterEnergy0->Fill(distanceXY, myedep[1]);
                    }
                    else if(dz1 >= 3 && dz1 < 6)
                    {
                        hScatterEnergy1->Fill(distanceXY, myedep[1]);
                    }
                    else if (dz2 >= 6 && dz2 < 9)
                    {
                        hScatterEnergy2->Fill(distanceXY, myedep[1]);
                    }
                    else if (dz2 >= 9 && dz2 < 12)
                    {
                        hScatterEnergy3->Fill(distanceXY, myedep[1]);
                    }
                    else if (dz2 >= 12 && dz2 < 15)
                    {
                        hScatterEnergy4->Fill(distanceXY, myedep[1]);
                    }
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

    // Save root and png files

    TFile *f = new TFile("root/histograms.root", "RECREATE");
    f->cd();

    hScatterEnergy0->SetOption("colz");
    hScatterEnergy1->SetOption("colz");
    hScatterEnergy2->SetOption("colz");
    hScatterEnergy3->SetOption("colz");
    hScatterEnergy4->SetOption("colz");

    hDistanceScattering->Write("hDistanceScattering");
    hDistanceCrystal->Write("hDistanceCrystal");
    hDistanceZ->Write("hDistanceZ");
    hScatterAngle->Write("hScatterAngle");

    hScatterEnergy0->Write("hScatterEnergy0");
    hScatterEnergy1->Write("hScatterEnergy1");
    hScatterEnergy2->Write("hScatterEnergy2");
    hScatterEnergy3->Write("hScatterEnergy3");
    hScatterEnergy4->Write("hScatterEnergy4");

    f->Close();

    TCanvas c1;
    c1.Clear();

    hDistanceScattering->Draw();
    c1.SaveAs("png/hDistanceScattering.png");
    c1.Clear();

    hDistanceCrystal->Draw();
    c1.SaveAs("png/hDistanceCrystal.png");
    c1.Clear();

    hDistanceZ->Draw();
    c1.SaveAs("png/hDistanceZ.png");
    c1.Clear();

    hScatterAngle->Draw();
    c1.SaveAs("png/hScatterAngle.png");
    c1.Clear();

    hScatterEnergy0->Draw("colz");
    c1.SaveAs("png/hScatterEnergy0.png");
    c1.Clear();

    hScatterEnergy1->Draw("colz");
    c1.SaveAs("png/hScatterEnergy1.png");
    c1.Clear();

    hScatterEnergy2->Draw("colz");
    c1.SaveAs("png/hScatterEnergy2.png");
    c1.Clear();

    hScatterEnergy3->Draw("colz");
    c1.SaveAs("png/hScatterEnergy3.png");
    c1.Clear();

    hScatterEnergy4->Draw("colz");
    c1.SaveAs("png/hScatterEnergy4.png");
    c1.Clear();

    return 0;
}









