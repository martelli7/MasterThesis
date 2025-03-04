// $Id: Analysis.hh 1 2025-3-3 martelli $
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
    // static Analysis* GetInstance(); //use it to access to the analysis class

    Analysis(bool &filter, bool &analyze, bool &scatterHisto, bool &scatterEnergy, bool &axisCone); //constructor

    ~Analysis(); //destructor

    void Run(const char* file, const char* tree);

    void FilterRootFile(const char* fileName, const char* treeName); //to filter a big root file in a smaller one

    void Analyze(const char* filteredFileName, const char* filteredTreeName);

private:
    // static Analysis* singleton;

    // this function takes in the energy deposit by the scattered 511 and returns its scatter angle by reversing Compton formula
    double scatterAngle (float &energydeposit);

    // this function takes in dx and dy and returns the sigma of the Theta Angle defined between the Y' local axis and the Cone Axis (delta vector)
    float sigmaThetaCone (float &dx, float &dy);

    float sigmaPhiCone (float &dx, float &dy, float &dz);

    // It returns the norm of a vector
    float norm (float &x1, float &y1, float &z1);

    float thetaAngle(float &x1, float &y1, float &z1, float &x2, float &y2, float &z2, float &norm1);


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
    float *myposX = new float[3];
    float *myposY = new float[3];
    float *myposZ = new float[3];
    float *mylocalPosX = new float[3];
    float *mylocalPosY = new float[3];
    float *mylocalPosZ = new float[3];
    float *myedep = new float[3];

    int *myphotonID = new int[3];
    int *myeventID = new int[3];
    int *mycrystalID = new int[3];
    int *myrsectorID = new int[3];
    int *mynCrystalCompton = new int[3];
    int *mynCrystalRayleigh = new int[3];

    // Variables
    Long64_t iEntry;
    Long64_t iEntry2;
    Long64_t iEntry3;

    float dx1, dy1, dz1, dx2, dy2, dz2;
    float dx1abs, dy1abs, dz1abs, dx2abs, dy2abs, dz2abs, dx1local, dx2local;
    int dxc1, dxc2, dyc1, dyc2;
    float distance, distanceXY, crydistance, distance1;        // distance is the norm of the DELTA VECTOR
    double sAngle;
    double thetaCone, alfaCone, thetaSigmaCone, alfaSigmaCone;

    // ROOT objects
    TFile *f;

    TH1D *hDistanceScattering;
    TH1D *hDistanceCrystal;
    TH1D *hDistanceZ;
    TH1D *hScatterAngle;

    TH2D *hScatterEnergy0;
    TH2D *hScatterEnergy1;
    TH2D *hScatterEnergy2;
    TH2D *hScatterEnergy3;
    TH2D *hScatterEnergy4;

    TH2D *hAxisCone;
    TH2D *hAxisCone1;

};

#endif
