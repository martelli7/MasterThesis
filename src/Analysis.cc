// $Id: Analysis.cc 7 2025-4-7 martelli $
/**
 * @file
 * @brief Analysis
 */

#include "../include/Analysis.hh"
#include <ROOT/RDataFrame.hxx>
#include <Math/Vector3D.h>
#include "Math/Vector3Dfwd.h"
#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/CoordinateSystemTags.h"
#include "Math/GenVector/PositionVector3D.h"
#include <Math/Rotation3D.h>
#include <Math/Transform3D.h>
#include "TFile.h"
#include <RtypesCore.h>
#include <TMath.h>
#include <cmath>
#include <cstdlib>
#include <cstring>

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
            hBetaAngle = new TH1D("hBetaAngle", "Distribution of #beta;#beta[degree];Counts", 380, 0,180);
        }

        if(doScatterEnergy)
        {
            hScatterEnergy0 = new TH2D("hScatterEnergy0", "#Deltay = [0#;3)mm;xz-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
            hScatterEnergy1 = new TH2D("hScatterEnergy1", "#Deltay = [3#;6)mm;xz-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
            hScatterEnergy2 = new TH2D("hScatterEnergy2", "#Deltay = [6#;9)mm;xz-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
            hScatterEnergy3 = new TH2D("hScatterEnergy3", "#Deltay = [9#;12)mm;xz-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
            hScatterEnergy4 = new TH2D("hScatterEnergy4", "#Deltay = [12#;15)mm;xz-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
            hScatterEnergy5 = new TH2D("hScatterEnergy5", "#Deltay < 0 mm;xz-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
            hScatterEnergy6 = new TH2D("hScatterEnergy6", "#Deltay >= 15 mm;xz-distance[mm];Energy[keV]", 80, 0, 40, 17, 0, 511);
        }

        if(doAxisCone)
        {
            hAxisThetaPhi = new TH2D("hAxisThetaPhi", "Angles of the Cone Axis;#theta[deg];#phi[deg]", 1800, 0, 180, 1800, 0, 180);
            hAxisBetaPhi = new TH2D("hAxisBetaPhi", "Distribution #beta vs #phi;#beta[deg];#phi[deg]", 1800, 0, 180, 1800, 0, 180);

            hAxisBetaTheta = new TH2D("hAxisBetaTheta", "Distribution #beta vs #theta;#beta[deg];#theta[deg]", 1800, 0, 180, 1800, 0, 180);
            hAxisBetaThetaOg = new TH2D("hAxisBetaThetaOg", "Distribution #beta vs #theta (Og);#beta[deg];#theta[deg]", 1800, 0, 180, 1800, 0, 180);

            hAxisPhiGamma = new TH2D("hAxisPhiGamma", "Distribution #phi vs #gamma;#phi[deg];#gamma[deg]", 1800, 0, 180, 1800, 0, 180);

            hAxisPhiATan2 = new TH2D("hAxisPhiATan2", "Distribution #phi vs ATan2;#phi[deg];ATan2", 380, -90, +90, 720, -180, 180);

            hAxisDifference = new TH1D("hAxisDifference", "Distribution of the difference #theta - #beta;#theta - #beta[deg];Counts", 3600, -360,360);
            hAxisDifferenceOg = new TH1D("hAxisDifferenceOg", "Distribution of the difference #theta - #beta (Og);#theta - #beta[deg];Counts", 3600, -360,360);

            hThetaAngle = new TH1D("hThetaAngle", "Distribution #theta ;#theta[deg];Counts", 380, 0,180);
            hThetaAngleOg = new TH1D("hThetaAngleOg", "Distribution #theta (Og) ;#theta[deg];Counts", 380, 0,180);
            hPhiAngle = new TH1D("hPhiAngle", "Distribution of #phi ;#phi[deg];Counts", 380, -90,90);
            // hThetaAngle = new TH1D("hThetaAngle", "Distribution #theta ;#theta[deg];Counts", 720, 0,360);
            // hPhiAngle = new TH1D("hPhiAngle", "Distribution of #phi ;#phi[deg];Counts", 720, 0,360);
            hGammaAngle = new TH1D("hGammaAngle", "Distribution of #gamma ;#gamma[deg];Counts", 380, 0,180);
            hAlfaAngle = new TH1D("hAlfaAngle", "Sum ;#alfa + #theta[deg];Counts", 380, 0,180);

            hPhiATan2 = new TH1D("hPhiATan2", "Distribution of ATan2 ;ATan2[deg];Counts", 720, -180,+180);

            hEnergy0 =  new TH1D("hEnergy0", "Energy = [0,100);#theta - #beta[deg];Counts", 500, -100, 100);
            hEnergy1 =  new TH1D("hEnergy1", "Energy = [100,200);#theta - #beta[deg];Counts", 500, -100, 100);
            hEnergy2 =  new TH1D("hEnergy2", "Energy = [200,300];#theta - #beta[deg];Counts", 500, -100, 100);
            hEnergy3 =  new TH1D("hEnergy3", "Energy = [300,340];#theta - #beta[deg];Counts", 500, -100, 100);


            hEnergySpectrum =  new TH1D("hEnergySpectrum", "Energy deposition;Energy [keV];Counts", 1040, 0, 520);
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

// It returns the norm of a vector
float Analysis::norm (float &compx, float &compy, float &compz)
{
    return TMath::Sqrt(compx*compx + compy*compy + compz*compz);
}


// **** beta is the scattering angle calculated from the energy deposition (should be equal to theta angle) **** //
float Analysis::betaAngle (float &energydeposit, int &countNaN)
{
    double enGamma2 = en511 - energydeposit*1000*Qe;
    double lambdaOut = h*c/enGamma2;

    float angleR = TMath::ACos( 2 - (511/(511-energydeposit)) );

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

void Analysis::thisEvent(float &x1, float &y1, float &z1, float &x2, float &y2, float &z2, float &locx1, float &locy1, float &locz1, float &locx2, float &locy2, float &locz2, float &edep1, float &edep2, int &crysID1, int &crysID2, int &countNaN, int &countDiff, int &countNaNPhi1, int &countNaNPhi2)
{
    // Compute distance
    float dx = x2 - x1;
    float dy = y2 - y1;
    float dz = z2 - z1;
    float dxlocal = locx2 - locx1;

    // using CoordSystem = ROOT::Math::Cartesian3D<double>;
    // using PositionTag = ROOT::Math::DefaultCoordinateSystemTag;

    ROOT::Math::XYZVector vecR1(x1,y1,z1);
    auto versR1 = vecR1.Unit();

    ROOT::Math::XYZVector vecR2(x2,y2,z2);

    ROOT::Math::XYZVector oldZ(0,0,1);

    auto vecDelta = vecR2 - vecR1;

    // std::cout <<vecDelta.X() << ", " <<vecDelta.Y() << ", " <<vecDelta.Z() << std::endl;

    // New SR

    // auto newZ = oldZ;

    auto newY = versR1;
    auto newX = versR1.Cross(oldZ).Unit();
    // auto newZ = oldZ;
    auto newZ = newX.Cross(versR1).Unit();

    // Correct axis for the computation of the angles Theta, Phi
    auto newwX = newZ;
    auto newwY = newX;
    auto newwZ = newY;

    ROOT::Math::Rotation3D rotMatrix(newwX, newwY, newwZ);

    ROOT::Math::XYZVector rotatedDelta = rotMatrix*vecDelta;

    // std::cout <<rotatedDelta.X() << ", " <<rotatedDelta.Y() << ", " <<rotatedDelta.Z() << std::endl;
    // std::cout << "\n" << std::endl;

    // New Components

    double newXcomp = vecDelta.Dot(newwX);
    double newYcomp = vecDelta.Dot(newwY);
    double newZcomp = vecDelta.Dot(newwZ);
    //
    // ROOT::Math::XYZVector newDelta(newXcomp,newYcomp,newZcomp);

    double thetaCone = TMath::RadToDeg()*rotatedDelta.Theta();
    double phiCone = TMath::RadToDeg()*rotatedDelta.Phi();

    // std::cout <<TMath::RadToDeg()*rotatedDelta.Theta()<< ", " <<TMath::RadToDeg()*rotatedDelta.Phi() << std::endl;
    // std::cout << "\n" << std::endl;

    float ATang2Phi = TMath::ATan2(newYcomp, newXcomp)*TMath::RadToDeg();


    // ROOT::Math::PositionVector3D<CoordSystem, PositionTag> vecR1;
    // vecR1.SetCoordinates(x1,y1,z1);
    //
    // ROOT::Math::PositionVector3D<CoordSystem, PositionTag> vecR2;
    // vecR2.SetCoordinates(x2,y2,z2);
    //
    // ROOT::Math::PositionVector3D<CoordSystem, PositionTag> axisZ;
    // axisZ.SetCoordinates(0,0,1);

    float distance = norm(dx, dy, dz);
    // std::cout <<distance << ", " << sqrt(rotatedDelta.X()*rotatedDelta.X() + rotatedDelta.Y()*rotatedDelta.Y() + rotatedDelta.Z()*rotatedDelta.Z()) << ", " <<sqrt(vecDelta.X()*vecDelta.X() + vecDelta.Y()*vecDelta.Y() + vecDelta.Z()*vecDelta.Z()) << std::endl; // solo per capire se almeno
    // std::cout << "\n" << std::endl;

    float distanceXY = TMath::Sqrt(distance*distance - dxlocal * dxlocal);

    // Compute crystal distance

    float dxc = TMath::Abs(crysID1/16 - crysID2/16);
    float dyc = TMath::Abs(crysID1%16 - crysID2%16);

    float crydistance = TMath::Sqrt(dxc*dxc + dyc*dyc);

    // float vectx = y1*dz - z1*dy; // R1 x Delta
    // float vecty = z1*dx - x1*dz;
    // float vectz = x1*dy - y1*dx;

    float vectx = y1; // R1 x asse z' (y' x z') -> ASSE X'
    float vecty = -x1;
    float vectz = 0;

    float normAsseX = norm(vectx, vecty, vectz);

    float compz1 = (-x1*z1)/(pow(norm(x1,y1,z1),2));
    float compz2 = (-y1*z1)/(pow(norm(x1,y1,z1),2));
    float compz3 = (pow(x1,2)+pow(y1,2))/(pow(norm(x1,y1,z1),2));

    float normz = norm(compz1, compz2, compz3);

    // if((511-edep1)/511 > 0.33)
    // {
    float betaCone = betaAngle(edep1, countNaN);
    float thetaConeOg = thetaAngle(x1, y1, z1, x2, y2, z2, distance);
    //float phiCone = TMath::ACos((norm(x1,y1,z1)*dz)/(norm(vectx,vecty,vectz)))*TMath::RadToDeg();
    //float phiCone = TMath::ACos(dz/(distance*TMath::Sqrt(1-pow((x1*dx + y1*dy + z1*dz)/(distance*norm(x1,y1,z1)),2))))*TMath::RadToDeg();
    //float phiCone = TMath::ACos(dz/(distance*TMath::Sin(thetaCone*TMath::DegToRad())))*TMath::RadToDeg();
    //float phiCone = acos((compz1*dx + compz3*dz)/(normz*distance*sin(thetaCone*TMath::RadToDeg())))*TMath::RadToDeg();
    //float phiCone = TMath::ASin((-vectx*dx -vecty*dy)/(normAsseX*distance*TMath::Sin(thetaCone*TMath::DegToRad())))*TMath::RadToDeg();
    float alfaCone = TMath::ACos(dz/(distance))*TMath::RadToDeg();

    float gammaCone = TMath::ACos((vectx*dx + vecty*dy)/(normAsseX*distance*TMath::Sin(thetaCone*TMath::DegToRad())))*TMath::RadToDeg();

    // float ATang2Phi = TMath::ATan2(-vecty*dy, -vectx*dx)*TMath::RadToDeg();
    // }


    //if(gammaCone > 89 && gammaCone < 91)
    if(phiCone != phiCone)
    {
        // std::cout << "Gamma = " << gammaCone << "; Phi = " << phiCone << "; Norm(R1) = " << norm(x1,y1,z1) << ", DeltaZ = " << dz << ", Distance = " << distance << ", beta = " << betaCone << ", thetaCone = " << thetaCone << ", edep = " << edep1 << std::endl;
        // std::cout << "Gamma + Phi = " << abs(gammaCone) + abs(phiCone) << std::endl;
        // std::cout << "Sin(Phi1) = " << (-vectx*dx -vecty*dy)/(normAsseX*distance*TMath::Sin(thetaCone*TMath::DegToRad())) << std::endl;
        countNaNPhi1++;
    }

    float phiCone2 = TMath::ACos(dz/(distance*TMath::Sin(thetaCone*TMath::DegToRad())))*TMath::RadToDeg();

    if(phiCone2!=phiCone2)
    {
        countNaNPhi2++;
    }


    hEnergySpectrum->Fill(edep1);

    if(locx2 - locx1 < 0)
    {
        // std::cout << locx1 << ", " << locx2 << std::endl;
    }
    if(edep1 < 340.666)
    {
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
            hAxisBetaThetaOg->Fill(betaCone, thetaConeOg);
            hAxisDifference->Fill(thetaCone - betaCone);
            hAxisDifferenceOg->Fill(thetaConeOg - betaCone);

            hThetaAngle->Fill(thetaCone);
            hThetaAngleOg->Fill(thetaConeOg);
            hPhiAngle->Fill(phiCone);
            hGammaAngle->Fill(gammaCone);

            hAxisPhiATan2->Fill(phiCone, ATang2Phi);
            hPhiATan2->Fill(ATang2Phi);
            //hAlfaAngle->Fill(thetaCone+alfaCone);

            hAxisPhiGamma->Fill(phiCone, gammaCone);


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
            else if (dxlocal >= 15)
            {
                hScatterEnergy6->Fill(distanceXY, edep1);
            }
            else if (dxlocal < 0)
            {
                hScatterEnergy5->Fill(distanceXY, edep1);
            }
        }
    }

}





// **** **** functions called in the run() function **** **** //

void Analysis::FilterRootFile(const char* fileName, const char* treeName)
{
    // Open the RootDataFrame
    ROOT::RDataFrame df(treeName, fileName);

    // Filtering
    auto df_filtered0 = df.Filter("photonID==1", "Select the first Photon only")
    .Filter("nCrystalCompton < 2","Select rows with only one Compton")
    //.Count();   //Count() does not work with Snapshot()
    .Snapshot("filteredPhoton", "filtered_Photon1.root", {"posX", "posY", "posZ", "edep", "localPosX", "localPosY", "localPosZ", "photonID", "eventID", "parentID", "crystalID", "rsectorID", "nCrystalCompton", "nCrystalRayleigh", "processName"});

    auto df_filtered1 = df.Filter("photonID==2", "Select the first Photon only")
    .Filter("nCrystalCompton < 2","Select rows with only one Compton")
    //.Count();   //Count() does not work with Snapshot()
    .Snapshot("filteredPhoton", "filtered_Photon2.root", {"posX", "posY", "posZ", "edep", "localPosX", "localPosY", "localPosZ", "photonID", "eventID", "parentID", "crystalID", "rsectorID", "nCrystalCompton", "nCrystalRayleigh", "processName"});
}





// **** **** this is the most important function **** **** //
void Analysis::Analyze(const char* filteredFileName, const char* filteredTreeName)
{
    TFile *myFile = new TFile(filteredFileName, "READ");
    TTree *tree = (TTree*)myFile->Get(filteredTreeName);

    Long64_t nEntries = tree->GetEntries(); // 10956826

    // Counters
    int count1 = 0;
    int count2 = 0;
    int count3 = 0;
    int countNaN = 0; //NaN values
    int countDiff = 0;
    int countNaNPhi1 = 0;
    int countNaNPhi2 = 0;

    // My arrays
    myposX = new float[2];
    myposY = new float[2];
    myposZ = new float[2];
    mylocalPosX = new float[2];
    mylocalPosY = new float[2];
    mylocalPosZ = new float[2];
    myedep = new float[2];

    myphotonID = new int[2];
    myeventID = new int[2];
    mycrystalID = new int[2];
    myrsectorID = new int[2];
    mynCrystalCompton = new int[2];
    mynCrystalRayleigh = new int[2];

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
    // tree->SetBranchAddress("processName", &processName);

    iEntry = 0;

    while (iEntry < nEntries)
    {
        tree->GetEntry(iEntry);

        myposX[0] = posX;
        myposY[0] = posY;
        myposZ[0] = posZ;
        mylocalPosX[0] = localPosX;
        myedep[0] = edep*1000; // Mev-->keV

        myeventID[0] = eventID;
        mycrystalID[0] = crystalID;

        tree->GetEntry(iEntry+1);

        myposX[1] = posX;
        myposY[1] = posY;
        myposZ[1] = posZ;
        mylocalPosX[1] = localPosX;
        myedep[1] = edep*1000; // Mev-->keV

        myeventID[1] = eventID;
        mycrystalID[1] = crystalID;
        myrsectorID[1] = rsectorID;

        thisEvent(myposX[0], myposY[0], myposZ[0], myposX[1], myposY[1], myposZ[1], mylocalPosX[0], mylocalPosY[0], mylocalPosZ[0], mylocalPosX[1], mylocalPosY[1], mylocalPosZ[1], myedep[0], myedep[1], mycrystalID[0], mycrystalID[1], countNaN, countDiff, countNaNPhi1, countNaNPhi2);

        if(myedep[0] > 340.666666 && myedep[0] <= 360)
        {
            count1++;
        }
        else if (myedep[0] > 360 && myedep[0] <= 400)
        {
            count2++;
        }
        else if(myedep[0] > 400)
        {
            count3++;
        }

        iEntry = iEntry+2;
    }

    // Save root and png files

    std::cout << "340.666-360, 360-400, >400; theta-beta > |5|" << std::endl;
    std::cout << count1 << ", " << count2 << ", " << count3 << "; " << countDiff << std::endl;
    std::cout << "NaN values are = " << countNaN << std::endl;
    std::cout << "Phi1 NaN values are = " << countNaNPhi1 << std::endl;
    std::cout << "Phi2 NaN values are = " << countNaNPhi2 << std::endl;

    f = new TFile("root/histograms.root", "RECREATE");
    f->cd();

    if(doScatterEnergy)
    {
        hScatterEnergy0->SetOption("colz");
        hScatterEnergy1->SetOption("colz");
        hScatterEnergy2->SetOption("colz");
        hScatterEnergy3->SetOption("colz");
        hScatterEnergy4->SetOption("colz");
        hScatterEnergy5->SetOption("colz");
        hScatterEnergy6->SetOption("colz");
    }
    if(doAxisCone)
    {
        hAxisThetaPhi->SetOption("colz");
        hAxisBetaPhi->SetOption("colz");
        hAxisBetaTheta->SetOption("colz");
        hAxisBetaThetaOg->SetOption("colz");
        hAxisPhiGamma->SetOption("colz");

        hAxisPhiATan2->SetOption("colz");
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
        hScatterEnergy5->Write("hScatterEnergy5");
        hScatterEnergy6->Write("hScatterEnergy6");
    }

    if(doAxisCone)
    {
        hAxisThetaPhi->Write("hAxisThetaPhi");
        hAxisBetaPhi->Write("hAxisBetaPhi");
        hAxisBetaTheta->Write("hAxisBetaTheta");
        hAxisBetaThetaOg->Write("hAxisBetaThetaOg");
        hAxisDifference->Write("hAxisDifference");
        hAxisDifferenceOg->Write("hAxisDifferenceOg");

        hAxisPhiGamma->Write("hAxisPhiGamma");

        hAxisPhiATan2->Write("hAxisPhiATan2");

        hThetaAngle->Write("hThetaAngle");
        hThetaAngleOg->Write("hThetaAngleOg");
        hPhiAngle->Write("hPhiAngle");
        hGammaAngle->Write("hGammaAngle");

        hPhiATan2->Write("hPhiATan2");

        hEnergy0->Write("hEnergy0");
        hEnergy1->Write("hEnergy1");
        hEnergy2->Write("hEnergy2");
        hEnergy3->Write("hEnergy3");

        hEnergySpectrum->Write("hEnergySpectrum");
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

        hScatterEnergy5->Draw("colz");
        c1.SaveAs("png/hScatterEnergy5.png");
        c1.Clear();

        hScatterEnergy6->Draw("colz");
        c1.SaveAs("png/hScatterEnergy6.png");
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

        hThetaAngle->Draw();
        c1.SaveAs("png/hThetaAngle.png");
        c1.Clear();

        hPhiAngle->Draw();
        c1.SaveAs("png/hPhiAngle.png");
        c1.Clear();

        hEnergy0->Draw();
        c1.SaveAs("png/hEnergy0.png");
        c1.Clear();

        hEnergy1->Draw();
        c1.SaveAs("png/hEnergy1.png");
        c1.Clear();

        hEnergy2->Draw();
        c1.SaveAs("png/hEnergy2.png");
        c1.Clear();

        hEnergy3->Draw();
        c1.SaveAs("png/hEnergy3.png");
        c1.Clear();
    }

    myFile->Close();
    delete myFile;

}


















