// $Id: Analysis.cc 6 2025-3-15 martelli $
/**
 * @file
 * @brief Analysis
 */

#include "../include/Analysis.hh"
#include <ROOT/RDataFrame.hxx>
#include "TFile.h"
#include <RtypesCore.h>
#include <TMath.h>
#include <cmath>

class TFile;


// **** Constructor **** //
Analysis::Analysis(bool &filter, bool &analyze, bool &scatterHisto, bool &scatterEnergy, bool &axisCone) :
doFilter(filter),
doAnalyze(analyze),
doScatterHisto(scatterHisto),
doScatterEnergy(scatterEnergy),
doAxisCone(axisCone)
{
    if(doAnalyze)
    {
        if(doScatterHisto)
        {
            hDistanceScattering = new TH1D("hDistanceScattering", "Scattering Distance;Distance [mm];Counts", 1000, 0, 30);
            hDistanceCrystal = new TH1D("hDistanceCrystal", "Scattering Crystal Distance;Distance [crystals];Counts", 1000, 0, 30);
            hDistanceZ = new TH1D("hDistanceZ", "Scattering Distance along y- axis (i.e. in depht);Distance [mm];Counts", 1000, -30, 30);
            hBetaAngle = new TH1D("hBetaAngle", "Distribution of #beta ;#beta[degree];Counts", 380, 0,180);
        }

        if(doScatterEnergy)
        {
            hScatterEnergy0 = new TH2D("hScatterEnergy0", "#Deltay = [0#;3)mm;xz-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
            hScatterEnergy1 = new TH2D("hScatterEnergy1", "#Deltay = [3#;6)mm;xz-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
            hScatterEnergy2 = new TH2D("hScatterEnergy2", "#Deltay = [6#;9)mm;xz-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
            hScatterEnergy3 = new TH2D("hScatterEnergy3", "#Deltay = [9#;12)mm;xz-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
            hScatterEnergy4 = new TH2D("hScatterEnergy4", "#Deltay = [12#;15)mm;xz-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
        }

        if(doAxisCone)
        {
            hAxisThetaPhi = new TH2D("hAxisThetaPhi", "Angles of the Cone Axis;#theta[deg];#phi[deg]", 1800, 0, 180, 1800, 0, 180);
            hAxisBetaPhi = new TH2D("hAxisBetaPhi", "Distribution #beta vs #phi;#beta[deg];#phi[deg]", 1800, 0, 180, 1800, 0, 180);

            hAxisBetaTheta = new TH2D("hAxisBetaTheta", "Distribution #beta vs #theta;#beta[deg];#theta[deg]", 1800, 0, 180, 1800, 0, 180);

            hAxisDifference = new TH1D("hAxisDifference", "Distribution of the difference #theta - #beta;#theta - #beta[deg];Counts", 3600, -360,360);

            hThetaAngle = new TH1D("hThetaAngle", "Distribution #theta ;#theta[deg];Counts", 380, 0,180);
            hPhiAngle = new TH1D("hPhiAngle", "Distribution of #phi ;#phi[deg];Counts", 380, 0,180);

            hEnergy0 =  new TH1D("hEnergy0", "Energy = [0,100);#theta - #beta[deg];Counts", 500, -100, 100);
            hEnergy1 =  new TH1D("hEnergy1", "Energy = [100,200);#theta - #beta[deg];Counts", 500, -100, 100);
            hEnergy2 =  new TH1D("hEnergy2", "Energy = [200,300];#theta - #beta[deg];Counts", 500, -100, 100);
            hEnergy3 =  new TH1D("hEnergy3", "Energy = [300,340];#theta - #beta[deg];Counts", 500, -100, 100);
        }
    }
}

Analysis::~Analysis()
{}

// **** Function called in the main file **** //
void Analysis::Run(const char* file, const char* tree)
{
    if(doFilter)
    {
        FilterRootFile(file, tree);
    }

    if(doAnalyze)
    {
        Analyze(file, tree);
    }
}

// **** beta is the scattering angle calculated from the energy deposition (should be equal to theta angle) **** //
float Analysis::betaAngle (float &energydeposit, int &count1, int &count2, int &count3, int &countNaN)
{
    double enGamma2 = en511 - energydeposit*1000*Qe;
    double lambdaOut = h*c/enGamma2;

    float angleR = TMath::ACos( 2 - (511/(511-energydeposit)) );

    if(energydeposit > 340.6666666 && energydeposit <= 360)
    {
        count1++;
        // std::cout << energydeposit << std::endl;
    }
    else if (energydeposit > 360 && energydeposit <= 400)
    {
        count2++;
    }
    else if(energydeposit > 400)
    {
        count3++;
    }

    if(angleR != angleR) // angleR is Nan
    {
        countNaN++;
        // std::cout << energydeposit << std::endl;
    }

    return angleR*TMath::RadToDeg(); // [Deg]
}


// **** theta is the scattering angle calculated from the dot product of the Y' (~R1 vector) local axis and the Cone Axis **** //
float Analysis::thetaAngle(float &x1, float &y1, float &z1, float &x2, float &y2, float &z2, float &norm1)
{
    float dx = x2 - x1;
    float dy = y2 - y1;
    float dz = z2 - z1;

    float dotR1Delta = x1*dx + y1*dy + z1*dz;
    float normDelta = norm(dx, dy, dz);
    float normR1 = norm(x1, y1, z1);

    float thetaR = TMath::ACos(dotR1Delta/(normDelta*normR1));

    return thetaR*TMath::RadToDeg(); // [deg]
}



// **** SIGMA, has to be checked **** //
float Analysis::sigmaThetaCone (float &dx, float &dy)
{
    float rxy2 = TMath::Power(dx,2) + TMath::Power(dy,2);

    float sigmaThetaR = TMath::Sqrt(TMath::Power(dx,2)*TMath::Power(sigmaY,2) + TMath::Power(dy,2)*TMath::Power(sigmaX,2)) / rxy2; // [Rad]

    return sigmaThetaR*TMath::RadToDeg(); // [Deg]
}

float Analysis::sigmaPhiCone (float &dx, float &dy, float &dz)
{
    float rxy2 = TMath::Power(dx,2) + TMath::Power(dy,2);

    float sigmaRxy = TMath::Sqrt(TMath::Power(dx,2)*TMath::Power(sigmaX,2) + TMath::Power(dy,2)*TMath::Power(sigmaY,2)/rxy2);

    float sigmaPhiR = TMath::Sqrt(TMath::Power(dz,2)*TMath::Power(sigmaRxy,2) + rxy2*TMath::Power(sigmaZ,2)) / (TMath::Power(dz,2)+rxy2); // [Rad]

    return sigmaPhiR*TMath::RadToDeg(); // [Deg]
}

void Analysis::thisEvent(float &x1, float &y1, float &z1, float &x2, float &y2, float &z2, float &locx1, float &locy1, float &locz1, float &locx2, float &locy2, float &locz2, float &edep1, float &edep2, int &crysID1, int &crysID2, int &count1, int &count2, int &count3, int &countNaN, int &countDiff)
{
    // Compute distance
    float dx = x2 - x1;
    float dy = y2 - y1;
    float dz = z2 - z1;
    float dxlocal = locx2 - locx1;

    float distance = norm(dx, dy, dz);
    float distanceXY = TMath::Sqrt(distance*distance - dxlocal * dxlocal);

    // Compute crystal distance

    float dxc = TMath::Abs(crysID1/16 - crysID2/16);
    float dyc = TMath::Abs(crysID1%16 - crysID2%16);

    float crydistance = TMath::Sqrt(dxc*dxc + dyc*dyc);

    // if((511-edep1)/511 > 0.33)
    // {
    float betaCone = betaAngle(edep1, count1, count2, count3, countNaN);
    float thetaCone = thetaAngle(x1, y1, z1, x2, y2, z2, distance);
    float phiCone = TMath::ACos(dz/distance)*TMath::RadToDeg();
    // }

    // Filling histograms
    if(doScatterHisto)
    {
        hDistanceScattering->Fill(distance);
        hDistanceCrystal->Fill(crydistance);
        hDistanceZ->Fill(dz);
        hBetaAngle->Fill(betaCone);
    }

    if(doAxisCone)
    {
        hAxisThetaPhi->Fill(thetaCone, phiCone);
        hAxisBetaPhi->Fill(betaCone, phiCone);
        hAxisBetaTheta->Fill(betaCone, thetaCone);
        hAxisDifference->Fill(thetaCone - betaCone);

        hThetaAngle->Fill(thetaCone);
        hPhiAngle->Fill(phiCone);

        if(edep1 >= 0 && edep1 < 100)
        {
            hEnergy0->Fill(thetaCone - betaCone);
        }
        else if(edep1 >= 100 && edep1 < 200)
        {
            hEnergy1->Fill(thetaCone - betaCone);
        }
        else if(edep1 >= 200 && edep1 < 300)
        {
            hEnergy2->Fill(thetaCone - betaCone);
        }
        else if(edep1 >= 300 && edep1 <= 340)
        {
            hEnergy3->Fill(thetaCone - betaCone);
        }

        if (thetaCone - betaCone > 5 || thetaCone - betaCone < - 5)
        {
            countDiff++;
        }
    }

    if(doScatterEnergy)
    {
        if (dxlocal >= 0 && dxlocal < 3)
        {
            hScatterEnergy0->Fill(distanceXY, edep1);
        }
        else if(dxlocal >= 3 && dxlocal < 6)
        {
            hScatterEnergy1->Fill(distanceXY, edep1);
        }
        else if (dxlocal >= 6 && dxlocal < 9)
        {
            hScatterEnergy2->Fill(distanceXY, edep1);
        }
        else if (dxlocal >= 9 && dxlocal < 12)
        {
            hScatterEnergy3->Fill(distanceXY, edep1);
        }
        else if (dxlocal >= 12 && dxlocal < 15)
        {
            hScatterEnergy4->Fill(distanceXY, edep1);
        }
    }

}


// It returns the norm of a vector
float Analysis::norm (float &x1, float &y1, float &z1)
{
    return TMath::Sqrt(x1*x1 + y1*y1 + z1*z1);
}


// **** **** functions called in the run() function **** **** //

void Analysis::FilterRootFile(const char* fileName, const char* treeName)
{
    // Open the RootDataFrame
    ROOT::RDataFrame df(treeName, fileName);

    // Filtering
    auto df_filtered = df.Filter("nCrystalCompton <= 1", "Rows (not events!!!) with nCrystalCompton = 0 or 1")
    .Filter("parentID==0","Photons which are not sons")
    .Filter("posX!=0 && posY!=0 && posZ!=0", "Events with photons")
    //.Filter()
    //.Count();   //Count() does not work with Snapshot()
    .Snapshot("filteredHits", "filtered_file4.root", {"posX", "posY", "posZ", "edep", "localPosX", "localPosY", "localPosZ", "photonID", "eventID", "parentID", "crystalID", "rsectorID", "nCrystalCompton", "nCrystalRayleigh"});
}





// **** **** this is the most important function **** **** //
void Analysis::Analyze(const char* filteredFileName, const char* filteredTreeName)
{
    std::unique_ptr<TFile> myFile( TFile::Open(filteredFileName) );
    auto tree = myFile->Get<TTree>(filteredTreeName);

    Long64_t nEntries = tree->GetEntries(); //number of entries inside each Branch (=9677397)

    // Counters
    int count1 = 0;
    int count2 = 0;
    int count3 = 0;
    int countNaN = 0; //NaN values
    int countDiff = 0; //

    // My arrays
    myposX = new float[3];
    myposY = new float[3];
    myposZ = new float[3];
    mylocalPosX = new float[3];
    mylocalPosY = new float[3];
    mylocalPosZ = new float[3];
    myedep = new float[3];

    myphotonID = new int[3];
    myeventID = new int[3];
    mycrystalID = new int[3];
    myrsectorID = new int[3];
    mynCrystalCompton = new int[3];
    mynCrystalRayleigh = new int[3];

    tree->SetBranchAddress("posX", &posX);
    tree->SetBranchAddress("posY", &posY);
    tree->SetBranchAddress("posZ", &posZ);
    tree->SetBranchAddress("localPosX", &localPosX);
    tree->SetBranchAddress("localPosY", &localPosY);
    tree->SetBranchAddress("localPosZ", &localPosZ);
    tree->SetBranchAddress("edep", &edep);
    tree->SetBranchAddress("photonID", &photonID);
    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("crystalID", &crystalID);
    tree->SetBranchAddress("rsectorID", &rsectorID);
    tree->SetBranchAddress("nCrystalCompton", &nCrystalCompton);
    tree->SetBranchAddress("nCrystalRayleigh", &nCrystalRayleigh);

    iEntry = 0;
    iEntry2 = 0;
    iEntry3 = 0;

    while (iEntry < nEntries)
    {
        // Load the data for the given tree entry
        tree->GetEntry(iEntry);

        int row = crystalID / 16;
        int col = crystalID % 16;

        if ( row >= 3 && row <= 11 && col >= 3 && row <= 11 ) // Select events on the central 9x9 matrix
        {
            // Transfer the current leaf to my arrays
            myposX[0] = posX;
            myposY[0] = posY;
            myposZ[0] = posZ;
            mylocalPosX[0] = localPosX;
            myedep[0] = edep*1000; // Mev-->keV

            myphotonID[0] = photonID;
            myeventID[0] = eventID;
            mycrystalID[0] = crystalID;
            myrsectorID[0] = rsectorID;
            mynCrystalCompton[0] = nCrystalCompton;
            mynCrystalRayleigh[0] = nCrystalRayleigh;

            iEntry2 = iEntry+1;
            iEntry3 = iEntry+2;

            // Call the second leaf
            tree->GetEntry(iEntry2);

            myposX[1] = posX;
            myposY[1] = posY;
            myposZ[1] = posZ;
            mylocalPosX[1] = localPosX;
            myedep[1] = edep*1000; // Mev-->keV

            myphotonID[1] = photonID;
            myeventID[1] = eventID;
            mycrystalID[1] = crystalID;
            myrsectorID[1] = rsectorID;
            mynCrystalCompton[1] = nCrystalCompton;
            mynCrystalRayleigh[1] = nCrystalRayleigh;

            // Call the third leaf
            tree->GetEntry(iEntry3);

            myposX[2] = posX;
            myposY[2] = posY;
            myposZ[2] = posZ;
            mylocalPosX[2] = localPosX;
            myedep[2] = edep*1000; // Mev-->keV

            myphotonID[2] = photonID;
            myeventID[2] = eventID;
            mycrystalID[2] = crystalID;
            myrsectorID[2] = rsectorID;
            mynCrystalCompton[2] = nCrystalCompton;
            mynCrystalRayleigh[2] = nCrystalRayleigh;

            // Retrieve eventID entries to check the same-event-photons + check there are no Rayleigh
            if( myeventID[0] == myeventID[1] && myeventID[1] == myeventID[2] && mynCrystalRayleigh[0]==0 && mynCrystalRayleigh[1]==0 && mynCrystalRayleigh[2]==0 )
            {
                if( mynCrystalCompton[0]==1 && mynCrystalCompton[1]==1 && mynCrystalCompton[2]==0 && myrsectorID[0]==myrsectorID[1] )
                {
                    //**** This is a Compton + Photoelectric ****//
                    thisEvent(myposX[0], myposY[0], myposZ[0], myposX[1], myposY[1], myposZ[1], mylocalPosX[0], mylocalPosY[0], mylocalPosZ[0], mylocalPosX[1], mylocalPosY[1], mylocalPosZ[1], myedep[0], myedep[1], mycrystalID[0], mycrystalID[1], count1, count2, count3, countNaN, countDiff);
                }
                else if ( mynCrystalCompton[0]==0 && mynCrystalCompton[1]==1 && mynCrystalCompton[2]==1 && myrsectorID[1]==myrsectorID[2] )
                {
                    //**** This is a Compton + Photoelectric ****//
                    thisEvent(myposX[1], myposY[1], myposZ[1], myposX[2], myposY[2], myposZ[2], mylocalPosX[1], mylocalPosY[1], mylocalPosZ[1], mylocalPosX[2], mylocalPosY[2], mylocalPosZ[2], myedep[1], myedep[2], mycrystalID[1], mycrystalID[2], count1, count2, count3, countNaN, countDiff);
                }
                // Now we skip the following 4th row (3 rows processed in one time)
                iEntry = iEntry+3;
            }
            else if( myeventID[0] == myeventID[1] && myeventID[1] != myeventID[2] && mynCrystalRayleigh[0]==0 && mynCrystalRayleigh[1]==0 && mynCrystalRayleigh[2]==0 )
            {
                if( mynCrystalCompton[0]==1 && mynCrystalCompton[1]==1 && myrsectorID[0]==myrsectorID[1] )
                {
                    //**** This is a Compton + Second photon is lost ****//
                    thisEvent(myposX[0], myposY[0], myposZ[0], myposX[1], myposY[1], myposZ[1], mylocalPosX[0], mylocalPosY[0], mylocalPosZ[0], mylocalPosX[1], mylocalPosY[1], mylocalPosZ[1], myedep[0], myedep[1], mycrystalID[0], mycrystalID[1], count1, count2, count3, countNaN, countDiff);
                }
                // Now we skip to the 3rd row (2 rows processed in one time)
                iEntry = iEntry+2;
            }
            else if( myeventID[0] != myeventID[1] && myeventID[1] == myeventID[2] && mynCrystalRayleigh[1]==0 && mynCrystalRayleigh[2]==0 )
            {
                if ( mynCrystalCompton[1]==1 && mynCrystalCompton[2]==1 && myrsectorID[1]==myrsectorID[2] )
                {
                    //**** This is a Compton + Second photon is lost ****//
                    thisEvent(myposX[1], myposY[1], myposZ[1], myposX[2], myposY[2], myposZ[2], mylocalPosX[1], mylocalPosY[1], mylocalPosZ[1], mylocalPosX[2], mylocalPosY[2], mylocalPosZ[2], myedep[1], myedep[2], mycrystalID[1], mycrystalID[2], count1, count2, count3, countNaN, countDiff);
                }
                // Now we skip to the 3rd row (2 rows processed in one time)
                iEntry = iEntry+2;
            }
            else
            {
                iEntry++; // Compton not found, we keep scrolling each row (just 1 row processed)
            }
        } //close if crystalID
        else
        {
            iEntry++; // The current event does not take place in the central matrix, we keep scrolling each row
        }
    } //close while

    // Save root and png files

    std::cout << count1 << ", " << count2 << ", " << count3 << "; " << countDiff << std::endl;
    std::cout << "NaN values are = " << countNaN << std::endl;

    f = new TFile("root/histograms.root", "RECREATE");
    f->cd();

    if(doScatterEnergy)
    {
        hScatterEnergy0->SetOption("colz");
        hScatterEnergy1->SetOption("colz");
        hScatterEnergy2->SetOption("colz");
        hScatterEnergy3->SetOption("colz");
        hScatterEnergy4->SetOption("colz");
    }
    if(doAxisCone)
    {
        hAxisThetaPhi->SetOption("colz");
        hAxisBetaPhi->SetOption("colz");
        hAxisBetaTheta->SetOption("colz");
    }

    if(doScatterHisto)
    {
        hDistanceScattering->Write("hDistanceScattering");
        hDistanceCrystal->Write("hDistanceCrystal");
        hDistanceZ->Write("hDistanceZ");
        hBetaAngle->Write("hBetaAngle");
    }

    if(doScatterEnergy)
    {
        hScatterEnergy0->Write("hScatterEnergy0");
        hScatterEnergy1->Write("hScatterEnergy1");
        hScatterEnergy2->Write("hScatterEnergy2");
        hScatterEnergy3->Write("hScatterEnergy3");
        hScatterEnergy4->Write("hScatterEnergy4");
    }

    if(doAxisCone)
    {
        hAxisThetaPhi->Write("hAxisThetaPhi");
        hAxisBetaPhi->Write("hAxisBetaPhi");
        hAxisBetaTheta->Write("hAxisBetaTheta");
        hAxisDifference->Write("hAxisDifference");

        hThetaAngle->Write("hThetaAngle");
        hPhiAngle->Write("hPhiAngle");

        hEnergy0->Write("hEnergy0");
        hEnergy1->Write("hEnergy1");
        hEnergy2->Write("hEnergy2");
        hEnergy3->Write("hEnergy3");
    }

    f->Close();

    TCanvas c1;
    c1.Clear();

    if(doScatterHisto)
    {
        hDistanceScattering->Draw();
        c1.SaveAs("png/hDistanceScattering.png");
        c1.Clear();

        hDistanceCrystal->Draw();
        c1.SaveAs("png/hDistanceCrystal.png");
        c1.Clear();

        hDistanceZ->Draw();
        c1.SaveAs("png/hDistanceZ.png");
        c1.Clear();

        hBetaAngle->Draw();
        c1.SaveAs("png/hBetaAngle.png");
        c1.Clear();
    }

    if(doScatterEnergy)
    {
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
    }

    if(doAxisCone)
    {
        hAxisThetaPhi->Draw("colz");
        c1.SaveAs("png/hAxisThetaPhi.png");
        c1.Clear();

        hAxisBetaPhi->Draw("colz");
        c1.SaveAs("png/hAxisBetaPhi.png");
        c1.Clear();

        hAxisBetaTheta->Draw("colz");
        c1.SaveAs("png/hAxisBetaTheta.png");
        c1.Clear();

        hAxisDifference->Draw();
        c1.SaveAs("png/hAxisDifference.png");
        c1.Clear();
    }


}


















