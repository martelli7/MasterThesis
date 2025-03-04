// $Id: Analysis.cc 1 2025-3-3 martelli $
/**
 * @file
 * @brief Analysis
 */

#include "Analysis.hh"
#include <ROOT/RDataFrame.hxx>
#include "TFile.h"
#include <cmath>

class TFile;

//**** 'singleton' is used to store the address of memory of 'analysis' ****//
// Analysis* Analysis::singleton = 0;

//**** GetIstance() method creates the pointer to the 'analysis' variable of Analysis type and returns its address of memory ****//
// Analysis* Analysis::GetInstance()
// {
//     if ( singleton == 0 )
//     {
//         static Analysis analysis; //**** this object is created only once, i.e. when 'singleton' is zero ****//
//         singleton = &analysis;
//     }
//     return singleton;
// }

Analysis::~Analysis()
{}

//**** Constructor ****//
Analysis::Analysis(bool &filter, bool &analyze, bool &scatterHisto, bool &scatterEnergy, bool &axisCone) :
doFilter(filter),
doAnalyze(analyze),
doScatterHisto(scatterHisto),
doScatterEnergy(scatterEnergy),
doAxisCone(axisCone)
{}

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

double Analysis::scatterAngle (float &energydeposit)
{
    double enGamma2 = (511*1000 - energydeposit*1000)*Qe;
    double lambdaOut = h*c/enGamma2;

    double angleR = TMath::ACos( 1 + (lambdaIn - lambdaOut)*me*c/h  ); // [Rad]

    // return angleR*TMath::RadToDeg(); // [Deg]
    return angleR; // [Rad]
}

float Analysis::sigmaThetaCone (float &dx, float &dy)
{
    float rxy2 = TMath::Power(dx,2) + TMath::Power(dy,2);

    float sigmaThetaR = TMath::Sqrt(TMath::Power(dx,2)*TMath::Power(sigmaY,2) + TMath::Power(dy,2)*TMath::Power(sigmaX,2)) / rxy2; // [Rad]

    float sigmaThetaD = sigmaThetaR*TMath::RadToDeg(); // [Deg]

    return sigmaThetaD;
}

float Analysis::sigmaPhiCone (float &dx, float &dy, float &dz)
{
    float rxy2 = TMath::Power(dx,2) + TMath::Power(dy,2);

    float sigmaRxy = TMath::Sqrt(TMath::Power(dx,2)*TMath::Power(sigmaX,2) + TMath::Power(dy,2)*TMath::Power(sigmaY,2)/rxy2);

    float sigmaPhiR = TMath::Sqrt(TMath::Power(dz,2)*TMath::Power(sigmaRxy,2) + rxy2*TMath::Power(sigmaZ,2)) / (TMath::Power(dz,2)+rxy2); // [Rad]

    float sigmaPhiD = sigmaPhiR*TMath::RadToDeg(); // [Deg]

    return sigmaPhiD;
}

// It returns the norm of a vector
float Analysis::norm (float &x1, float &y1, float &z1)
{
    return TMath::Sqrt(x1*x1 + y1*y1 + z1*z1);
}

// Theta is the angle between the Y' local axis and the Cone Axis
// double thetaAngle(float &x1, float &y1, float &z1, float &x2, float &y2, float &z2, float &normDelta)
// {
//     float dotR1R2 = x1*x2 + y1*y2 + z1*z2;
//     float normR1 = norm(x1, y1, z1);
//     // float dx = x2 - x1;
//     // float dy = y2 - y1;
//     // float dz = z2 - z1;
//
//
//     // float normDelta1 = norm(dx, dy, dz);
//
//     return dotR1R2/(normDelta*normR1) - normR1/normDelta; // cos
//     // return TMath::ACos(dotR1R2/(normDelta*normR1) - normR1/normDelta)*TMath::RadToDeg(); // [deg]
// }
float Analysis::thetaAngle(float &x1, float &y1, float &z1, float &x2, float &y2, float &z2, float &norm1)
{
    // float normR1 = norm(x1, y1, z1);
    // float normDelta = norm(dx,dy,dz);
    //
    // float dotR1Delta = (dx/normDelta)*(x1/normR1) + (dy/normDelta)*(y1/normR1) + (dz/normDelta)*(z1/normR1);

    float dx = x2 - x1;
    float dy = y2 - y1;
    float dz = z2 - z1;


    float dotR1Delta = x1*dx + y1*dy + z1*dz;
    float normDelta = norm(dx, dy, dz);
    float normR1 = norm(x1, y1, z1);

    return dotR1Delta/(normDelta*normR1); // cos
    // return TMath::ACos(dotR1R2/(normDelta*normR1) - normR1/normDelta)*TMath::RadToDeg(); // [deg]
}



void Analysis::FilterRootFile(const char* fileName, const char* treeName)
{
    // Open the RootDataFrame
    ROOT::RDataFrame df(treeName, fileName);

    // Filtering
    auto df_filtered = df.Filter("nCrystalCompton <= 1", "Rows (not events!!!) with nCrystalCompton = 0 or 1")
        .Filter("parentID==0","Photons which are not daughters")
        .Filter("posX!=0 && posY!=0 && posZ!=0", "Events with photons")
        //.Filter("photonID==1","Select the first photon only") //non usare
        //.Filter()   //possibile concatenare più filtri
        //.Count();   //Count() does not work with Snapshot()
        .Snapshot("filteredHits", "filtered_file4.root", {"posX", "posY", "posZ", "edep", "localPosX", "localPosY", "localPosZ", "photonID", "eventID", "parentID", "crystalID", "rsectorID", "nCrystalCompton", "nCrystalRayleigh"});
}

void Analysis::Analyze(const char* filteredFileName, const char* filteredTreeName)
{
    std::unique_ptr<TFile> myFile( TFile::Open(filteredFileName) );
    auto tree = myFile->Get<TTree>(filteredTreeName);

    Long64_t nEntries = tree->GetEntries(); //number of entries inside each Branch (=9677397)

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

    if(doScatterHisto)
    {
        hDistanceScattering = new TH1D("hDistanceScattering", "Scattering Distance;Distance [mm];Counts", 1000, 0, 30);
        hDistanceCrystal = new TH1D("hDistanceCrystal", "Crystal Scattering Distance;Distance [crystals];Counts", 1000, 0, 30);
        hDistanceZ = new TH1D("hDistanceZ", "Scattering Distance along x-local axis (depht);Distance [mm];Counts", 1000, -10, 30);
        hScatterAngle = new TH1D("hScatterAngle", "Scattering Angle;Angle[#degree];Counts", 380, -5,180);
    }

    if(doScatterEnergy)
    {
        hScatterEnergy0 = new TH2D("hScatterEnergy0", "#Deltaz = [0#;3)mm;xy-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
        hScatterEnergy1 = new TH2D("hScatterEnergy1", "#Deltaz = [3#;6)mm;xy-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
        hScatterEnergy2 = new TH2D("hScatterEnergy2", "#Deltaz = [6#;9)mm;xy-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
        hScatterEnergy3 = new TH2D("hScatterEnergy3", "#Deltaz = [9#;12)mm;xy-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
        hScatterEnergy4 = new TH2D("hScatterEnergy4", "#Deltaz = [12#;15)mm;xy-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
    }

    if(doAxisCone)
    {
        hAxisCone = new TH2D("hAxisCone", "Angles of the Cone Axis;#theta[deg];#phi[deg]", 2000, -1, 1, 2000, -1, 1);
        hAxisCone1 = new TH2D("hAxisCone1", "Angles of the Cone Axis;#beta[deg];#phi[deg]", 2000, -1, 1, 2000, -1, 1);
    }

    while (iEntry < nEntries)
    {
        // Load the data for the given tree entry
        tree->GetEntry(iEntry);

        if ( crystalID == 119 ) // Select events on the central crystal only ################ errore logico
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

            // Retrieve eventID entries to check the same-event-photons
            if( myeventID[0] == myeventID[1] && myeventID[1] == myeventID[2] )
            {
                if( mynCrystalRayleigh[0]==0 && mynCrystalRayleigh[1]==0 && mynCrystalRayleigh[2]==0 ) // no rayleigh scattering
                {
                    // This event has a scattering event!
                    // have to check which couple (0,1 or 1,2) is actually a scattering

                    // First couple
                    dx1 = myposX[1] - myposX[0]; // these are the coordinates of the second point in the SR of the detector, i.e. the #Delta vector
                    dy1 = myposY[1] - myposY[0];
                    dz1 = myposZ[1] - myposZ[0];

                    dx1local = mylocalPosX[1] - mylocalPosX[0]; // in the GATE SR of the detector, x-axis is the depth!!!

                    dx1abs = TMath::Abs( myposX[1] - myposX[0] );
                    dy1abs = TMath::Abs( myposY[1] - myposY[0] );
                    dz1abs = TMath::Abs( myposZ[1] - myposZ[0] );

                    // Second couple
                    dx2 = myposX[2] - myposX[1];
                    dy2 = myposY[2] - myposY[1];
                    dz2 = myposZ[2] - myposZ[1];

                    dx2local = mylocalPosX[2] - mylocalPosX[1];

                    dx2abs = TMath::Abs( myposX[2] - myposX[1]);
                    dy2abs = TMath::Abs( myposY[2] - myposY[1]);
                    dz2abs = TMath::Abs( myposZ[2] - myposZ[1]);

                    // Compute crystal distance...
                    // First couple
                    dxc1 = TMath::Abs(mycrystalID[0]/16 - mycrystalID[1]/16);
                    dyc1 = TMath::Abs(mycrystalID[0]%16 - mycrystalID[1]%16);

                    // Second couple
                    dxc2 = TMath::Abs(mycrystalID[1]/16 - mycrystalID[2]/16);
                    dyc2 = TMath::Abs(mycrystalID[1]%16 - mycrystalID[2]%16);

                    // First couple condition
                    if( mynCrystalCompton[0]==1 && mynCrystalCompton[1]==1 && mynCrystalCompton[2]==0 )
                    {
                        distance = TMath::Sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
                        // distance1 = norm(dx1, dy1, dz1);
                        // std::cout << "differenza = " << distance1 - distance << std::endl;
                        // distanceXY = TMath::Sqrt(dx1 * dx1 + dy1 * dy1);
                        distanceXY = TMath::Sqrt(distance*distance - dx1local * dx1local);

                        crydistance = TMath::Sqrt(dxc1*dxc1 + dyc1*dyc1);

                        sAngle = TMath::Cos(scatterAngle(myedep[0]));
                        // sAngle = scatterAngle(myedep[0]);

                        // Display some outputs...
                        //std::cout << myphotonID[i] << ", " << i << ", "<< myeventID[i] << ", (" << dx1 << "," << dy1 << "," << dz1 << ") " << "(" << dxc1 << "," << dyc1 << ")" <<  std::endl;

                        // Filling histograms
                        if(doScatterHisto)
                        {
                            hDistanceScattering->Fill(distance);
                            hDistanceCrystal->Fill(crydistance);
                            hDistanceZ->Fill(dz1);
                            hScatterAngle->Fill(sAngle);
                        }

                        if(doAxisCone)
                        {
                            // Compute Angles of the Axis of the Cone!!!
                            thetaCone = thetaAngle(myposX[0], myposY[0], myposZ[0], myposX[1], myposY[1], myposZ[1], distance);
                            // thetaCone = thetaAngle(myposX[0], myposY[0], myposZ[0], dx1, dy1, dz1);
                            // alfaCone = TMath::ACos(dz1/distance)*TMath::RadToDeg();
                            alfaCone = dz1/distance;

                            // std::cout << thetaCone << ", " <<alfaCone << std::endl;

                            if (thetaCone > -2)
                            {
                                hAxisCone->Fill(thetaCone, alfaCone);
                            }

                            if (sAngle > -2)
                            {
                                hAxisCone1->Fill(sAngle, alfaCone);
                            }
                        }

                        if(doScatterEnergy)
                        {
                            if (dx1local >= 0 && dx1local < 3)
                            {
                                hScatterEnergy0->Fill(distanceXY, myedep[0]);
                            }
                            else if(dx1local >= 3 && dx1local < 6)
                            {
                                hScatterEnergy1->Fill(distanceXY, myedep[0]);
                            }
                            else if (dx1local >= 6 && dx1local < 9)
                            {
                                hScatterEnergy2->Fill(distanceXY, myedep[0]);
                            }
                            else if (dx1local >= 9 && dx1local < 12)
                            {
                                hScatterEnergy3->Fill(distanceXY, myedep[0]);
                            }
                            else if (dx1local >= 12 && dx1local < 15)
                            {
                                hScatterEnergy4->Fill(distanceXY, myedep[0]);
                            }
                        }
                    }
                    //Second couple condition
                    else if ( mynCrystalCompton[0]==0 && mynCrystalCompton[1]==1 && mynCrystalCompton[2]==1 )
                    {
                        distance = TMath::Sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
                        // distanceXY = TMath::Sqrt(dx2 * dx2 + dy2 * dy2);
                        distanceXY = TMath::Sqrt(distance*distance - dx2local * dx2local);

                        crydistance = TMath::Sqrt(dxc2*dxc2 + dyc2*dyc2);

                        sAngle = TMath::Cos(scatterAngle(myedep[1]));
                        // sAngle = scatterAngle(myedep[1]);

                        // Display some outputs...
                        //std::cout << myphotonID[i] << ", " << i << ", "<< myeventID[i] << ", (" << dx2 << "," << dy2 << "," << dz2 << ") " << "(" << dxc1 << "," << dyc1 << ")" <<  std::endl;


                        // Filling histograms
                        if(doScatterHisto)
                        {
                            hDistanceScattering->Fill(distance);
                            hDistanceCrystal->Fill(crydistance);
                            hDistanceZ->Fill(dz2);
                            hScatterAngle->Fill(sAngle);
                        }

                        if(doAxisCone)
                        {
                            // Compute Angles of the Axis of the Cone!!!
                            thetaCone = thetaAngle(myposX[1], myposY[1], myposZ[1], myposX[2], myposY[2], myposZ[2], distance);
                            // thetaCone = thetaAngle(myposX[1], myposY[1], myposZ[1], dx2, dy2, dz2);
                            // alfaCone = TMath::ACos(dz2/distance)*TMath::RadToDeg();
                            alfaCone = dz2/distance;
                            if (thetaCone > -2)
                            {
                                hAxisCone->Fill(thetaCone, alfaCone);
                            }

                            if (sAngle > -2)
                            {
                                hAxisCone1->Fill(sAngle, alfaCone);
                            }
                        }

                        if(doScatterEnergy)
                        {
                            if (dx2local >= 0 && dx2local < 3)
                            {
                                hScatterEnergy0->Fill(distanceXY, myedep[1]);
                            }
                            else if(dx2local >= 3 && dx2local < 6)
                            {
                                hScatterEnergy1->Fill(distanceXY, myedep[1]);
                            }
                            else if (dx2local >= 6 && dx2local < 9)
                            {
                                hScatterEnergy2->Fill(distanceXY, myedep[1]);
                            }
                            else if (dx2local >= 9 && dx2local < 12)
                            {
                                hScatterEnergy3->Fill(distanceXY, myedep[1]);
                            }
                            else if (dx2local >= 12 && dx2local < 15)
                            {
                                hScatterEnergy4->Fill(distanceXY, myedep[1]);
                            }
                        }
                    }
                }
                // If there is a rayleigh scattering or more than 1 compton we skipped to the 4th event without processing
                // Now we skip the following 4th event (3 events processed in one time)
                iEntry = iEntry+3;
            }
            else if( myeventID[0] == myeventID[1] && myeventID[1] != myeventID[2] )
            {
                if( mynCrystalRayleigh[0]==0 && mynCrystalRayleigh[1]==0 )
                {
                    //**** This is a scattering + second photon is lost ****//

                    // First couple
                    dx1 = myposX[1] - myposX[0];
                    dy1 = myposY[1] - myposY[0];
                    dz1 = myposZ[1] - myposZ[0];

                    dx1local = mylocalPosX[1] - mylocalPosX[0];

                    dx1abs = TMath::Abs( myposX[1] - myposX[0] );
                    dy1abs = TMath::Abs( myposY[1] - myposY[0] );
                    dz1abs = TMath::Abs( myposZ[1] - myposZ[0] );

                    // Compute crystal distance...
                    // First couple
                    dxc1 = TMath::Abs(mycrystalID[0]/16 - mycrystalID[1]/16);
                    dyc1 = TMath::Abs(mycrystalID[0]%16 - mycrystalID[1]%16);

                    crydistance = TMath::Sqrt(dxc1*dxc1 + dyc1*dyc1);

                    // First couple condition
                    if( mynCrystalCompton[0]==0 && mynCrystalCompton[1]==1 )
                    {
                        distance = TMath::Sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
                        // distanceXY = TMath::Sqrt(dx1 * dx1 + dy1 * dy1);
                        distanceXY = TMath::Sqrt(distance*distance - dx1local * dx1local);

                        sAngle = TMath::Cos(scatterAngle(myedep[0]));
                        // sAngle = scatterAngle(myedep[0]);

                        // Display some outputs...
                        //std::cout << myphotonID[i] << ", " << i << ", "<< myeventID[i] << ", (" << dx1 << "," << dy1 << "," << dz1 << ")"<< std::endl;

                        // Filling histograms
                        if(doScatterHisto)
                        {
                            hDistanceScattering->Fill(distance);
                            hDistanceCrystal->Fill(crydistance);
                            hDistanceZ->Fill(dz1);
                            hScatterAngle->Fill(sAngle);
                        }

                        if(doAxisCone)
                        {
                            // Compute Angles of the Axis of the Cone!!!
                            thetaCone = thetaAngle(myposX[0], myposY[0], myposZ[0], myposX[1], myposY[1], myposZ[1], distance);
                            // thetaCone = thetaAngle(myposX[0], myposY[0], myposZ[0], dx1, dy1, dz1);
                            // alfaCone = TMath::ACos(dz1/distance)*TMath::RadToDeg();
                            alfaCone = dz1/distance;

                            if (thetaCone > -2)
                            {
                                hAxisCone->Fill(thetaCone, alfaCone);

                            }

                            if (sAngle > -2)
                            {
                                hAxisCone1->Fill(sAngle, alfaCone);
                            }
                        }

                        if(doScatterEnergy)
                        {
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
                    }
                }
                // If there is a rayleigh scattering we skip to the 3rd event without processing
                // Now we skip to the 3rd event (2 events processed in one time)
                iEntry = iEntry+2;
            }
            else if( myeventID[0] != myeventID[1] && myeventID[1] == myeventID[2] )
            {
                if( mynCrystalRayleigh[1]==0 && mynCrystalRayleigh[2]==0 )
                {
                    // This is a scattering + second photon is lost

                    // Second couple
                    dx2 = myposX[2] - myposX[1];
                    dy2 = myposY[2] - myposY[1];
                    dz2 = myposZ[2] - myposZ[1];

                    dx2abs = TMath::Abs( myposX[2] - myposX[1] );
                    dy2abs = TMath::Abs( myposY[2] - myposY[1] );
                    dz2abs = TMath::Abs( myposZ[2] - myposZ[1] );

                    // Second couple
                    dxc2 = TMath::Abs(mycrystalID[1]/16 - mycrystalID[2]/16);
                    dyc2 = TMath::Abs(mycrystalID[1]%16 - mycrystalID[2]%16);

                    // Second couple condition
                    if ( mynCrystalCompton[1]==1 && mynCrystalCompton[2]==1 )
                    {
                        distance = TMath::Sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
                        // distanceXY = TMath::Sqrt(dx2 * dx2 + dy2 * dy2);
                        distanceXY = TMath::Sqrt(distance*distance - dx2local * dx2local);

                        crydistance = TMath::Sqrt(dxc2*dxc2 + dyc2*dyc2);

                        sAngle = TMath::Cos(scatterAngle(myedep[1]));
                        // sAngle = scatterAngle(myedep[1]);

                        // Display some outputs...
                        //std::cout << myphotonID[i] << ", " << i << ", "<< myeventID[i] << ", (" << dx2 << "," << dy2 << "," << dz2 << ")"<< std::endl;

                        // Compute Angles of the Axis of the Cone!!!
                        thetaCone = thetaAngle(myposX[1], myposY[1], myposZ[1], myposX[2], myposY[2], myposZ[2], distance);
                        // thetaCone = thetaAngle(myposX[1], myposY[1], myposZ[1], dx2, dy2, dz2);
                        // alfaCone = TMath::ACos(dz2/distance)*TMath::RadToDeg();
                        alfaCone = dz2/distance;

                        // Filling histograms
                        if(doScatterHisto)
                        {
                            hDistanceScattering->Fill(distance);
                            hDistanceCrystal->Fill(crydistance);
                            hDistanceZ->Fill(dz2);
                            hScatterAngle->Fill(sAngle);
                        }

                        if(doAxisCone)
                        {
                            if (thetaCone > -2)
                            {
                                hAxisCone->Fill(thetaCone, alfaCone);
                            }

                            if (sAngle > -2)
                            {
                                hAxisCone1->Fill(sAngle, alfaCone);
                            }
                        }

                        if(doScatterEnergy)
                        {
                            if (dz2 >= 0 && dz2 < 3)
                            {
                                hScatterEnergy0->Fill(distanceXY, myedep[1]);
                            }
                            else if(dz2 >= 3 && dz2 < 6)
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
                    }
                }
                // If there is a rayleigh scattering we skip to the 3rd event without processing
                // Now we skip to the 3rd event (2 events processed in one time)
                iEntry = iEntry+2;
            }
            else
            {
                iEntry++; // Compton not found, we keep scrolling each eventID (just 1 event processed)
            }
        } //close if crystalID
        else
        {
            iEntry++; // The current event does not take place in the central crystal, we keep scrolling each eventID
        }
    } //close while

    // Save root and png files

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
        hAxisCone->SetOption("colz");
        hAxisCone1->SetOption("colz");
    }

    if(doScatterHisto)
    {
        hDistanceScattering->Write("hDistanceScattering");
        hDistanceCrystal->Write("hDistanceCrystal");
        hDistanceZ->Write("hDistanceZ");
        hScatterAngle->Write("hScatterAngle");
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
        hAxisCone->Write("hAxisCone");
        hAxisCone1->Write("hAxisCone1");
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

        hScatterAngle->Draw();
        c1.SaveAs("png/hScatterAngle.png");
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
        hAxisCone->Draw("colz");
        c1.SaveAs("png/hAxisCone.png");
        c1.Clear();

        hAxisCone1->Draw("colz");
        c1.SaveAs("png/hAxisCone1.png");
        c1.Clear();
    }


}


















