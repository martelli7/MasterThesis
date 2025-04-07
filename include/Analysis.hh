// $Id: Analysis.hh 7 2025-4-7 martelli $
/**
 * @file
 * @brief Analysis
 */

/*
 * Analysis.hh
 *
 *  Created on: 9 Feb 2010
 *      Author: adotti
 */

#include <cstddef>
#ifndef ANALYSIS_HH
#define ANALYSIS_HH 1

//ROOT
#include "TTree.h"
#include "TFile.h"
#include <ROOT/RVec.hxx>
#include <RtypesCore.h>
#include <TBranch.h>
#include <TMathBase.h>
#include <TObjArray.h>
#include "TH1.h"
#include <TGraphErrors.h>
#include "TMath.h"

#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>

//C++
#include <cmath>
#include <string>

// Physics constants
const double h = TMath::H();     // Planck [J*s]
const double c = TMath::C();     // light speed [m/s]
const double me = 9.10938356e-31;    // m electron [kg]
const double Qe = TMath::Qe(); // q electron [eV->J]
const double en511 = 511*1000*Qe; // [J]
const double lambdaIn = h*c/en511; // lambda of 511keV [m]

// X, Y, Z Sigma
const float sigmaX = 10; // [mm]
const float sigmaY = 10; // [mm]
const float sigmaZ = 10; // [mm]


class TH2D;
class TFile;

class Analysis {

public:

    Analysis(bool &filter, bool &analyze, bool &scatterHisto, bool &scatterEnergy, bool &axisCone); //constructor

    ~Analysis(); //destructor

    void Run(const char* file, const char* tree);

    void FilterRootFile(const char* fileName, const char* treeName); //to filter a big root file in a smaller one

    void Analyze(const char* filteredFileName, const char* filteredTreeName);

private:

    // this function takes in the energy deposit by the scattered 511 and returns its scatter angle by reversing Compton formula
    float betaAngle (float &energydeposit, int &countNaN);

    // It returns the norm of a vector
    float norm (float &x1, float &y1, float &z1);

    float thetaAngle (float &x1, float &y1, float &z1, float &x2, float &y2, float &z2, float &norm1);

    // SIGMA
    // this function takes in dx and dy and returns the sigma of the Theta Angle defined between the Y' local axis and the Cone Axis (delta vector)
    float sigmaThetaCone (float &dx, float &dy);

    float sigmaPhiCone (float &dx, float &dy, float &dz);


    void thisEvent(float &x1, float &y1, float &z1, float &x2, float &y2, float &z2, float &locx1, float &locy1, float &locz1, float &locx2, float &locy2, float &locz2, float &edep1, float &edep2, int &crysID1, int &crysID2, int &countNaN, int &countDiff, int &countNaNPhi1, int &countNaNPhi2);

    // Booleans
    bool doFilter, doAnalyze, doScatterHisto, doScatterEnergy, doAxisCone;

    // Branches
    float posX;
    float posY;
    float posZ;
    float localPosX;
    float localPosY;
    float localPosZ;
    float edep;

    int photonID;
    int eventID;
    int crystalID;
    int rsectorID;
    int nCrystalCompton;
    int nCrystalRayleigh;

    // My arrays
    float *myposX;
    float *myposY;
    float *myposZ;
    float *mylocalPosX;
    float *mylocalPosY;
    float *mylocalPosZ;
    float *myedep;

    int *myphotonID;
    int *myeventID;
    int *mycrystalID;
    int *myrsectorID;
    int *mynCrystalCompton;
    int *mynCrystalRayleigh;

    // Counters
    Long64_t iEntry;

    // ROOT objects
    TFile *f;

    // ROOT objects
    TH1D *hDistanceScattering = nullptr;
    TH1D *hDistanceCrystal = nullptr;
    TH1D *hDistanceZ = nullptr;
    TH1D *hBetaAngle = nullptr;

    TH2D *hScatterEnergy0 = nullptr;
    TH2D *hScatterEnergy1 = nullptr;
    TH2D *hScatterEnergy2 = nullptr;
    TH2D *hScatterEnergy3 = nullptr;
    TH2D *hScatterEnergy4 = nullptr;
    TH2D *hScatterEnergy5 = nullptr;
    TH2D *hScatterEnergy6 = nullptr;

    TH2D *hAxisThetaPhi = nullptr;
    TH2D *hAxisBetaPhi = nullptr;
    TH2D *hAxisBetaTheta = nullptr;
    TH2D *hAxisBetaThetaOg = nullptr;
    TH2D *hAxisPhiGamma = nullptr;
    TH2D *hAxisPhiATan2 = nullptr;


    TH1D *hAxisDifference = nullptr;
    TH1D *hAxisDifferenceOg = nullptr;
    TH1D *hThetaAngle = nullptr;
    TH1D *hThetaAngleOg = nullptr;
    TH1D *hPhiAngle = nullptr;
    TH1D *hGammaAngle = nullptr;
    TH1D *hAlfaAngle = nullptr;

    TH1D *hPhiATan2 = nullptr;

    TH1D *hEnergy0 =  nullptr;
    TH1D *hEnergy1 =  nullptr;
    TH1D *hEnergy2 =  nullptr;
    TH1D *hEnergy3 =  nullptr;

    TH1D *hEnergySpectrum = nullptr;
};

#endif
