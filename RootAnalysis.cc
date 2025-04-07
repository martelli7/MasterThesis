// $Id: RootAnalysis.cc 6 2025-3-3 martelli $
/**
 * @file   RootAnalysis.cc
 *
 * @date   3 Mar 2025
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
#include <TGraphErrors.h>
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

#include "Analysis.hh"
#include "include/Analysis.hh"


// commands for ROOT terminal:
// TTree *tree = (TTree*)_file0->Get("Hits;-32768");
// tree->Draw("parentID>>(200,0,200)");
// tree->Draw("posX:posY","nCrystalCompton==1 & eventID==137","colz");

/*const std::vector<Float_t>& posX,
 * const std::vector<Float_t>& posY,
 * const std::vector<Float_t>& posZ) {
 *
 * root [6]  tree->Scan("eventID:rsectorID:photonID:processName:nCrystalCompton:nCrystalRayleigh:posX")
 ************************************************************************************
 *    Row   *   eventID * rsectorID *  photonID * nCrystalC * nCrystalR *      posX *
 ************************************************************************************
 *        0 *         4 *         4 *         2 *         1 *         0 * -33.72324 *
 *        1 *         4 *         4 *         2 *         1 *         0 * -33.86138 *
 * unica coppia compton
 *
 *        2 *        16 *        10 *         1 *         0 *         0 * 58.599784 *
 *        3 *        16 *         4 *         2 *         0 *         1 * -56.17806 *
 *        4 *        16 *         4 *         2 *         1 *         1 * -60.42935 *
 *        5 *        16 *         4 *         2 *         2 *         1 * -63.02994 *
 *        6 *        16 *         4 *         2 *         2 *         1 * -64.34118 *
 * fotoelettrico + righe 3,4 -> rayleigh, righe 4,5 -> primo compton, righe 5,6 -> secondo compton
 *
 *        7 *        17 *         8 *         1 *         1 *         0 * -56.12685 *
 *        8 *        17 *         8 *         1 *         1 *         0 * -56.98064 *
 *        9 *        17 *         2 *         2 *         0 *         0 * 58.198852 *
 * unica coppia compton + fotoelettrico
 *
 *       10 *        22 *         8 *         1 *         1 *         0 * -31.88264 *
 *       11 *        22 *         8 *         1 *         2 *         0 * -32.00119 *
 *       12 *        22 *         8 *         1 *         3 *         0 * -39.27561 *
 *       13 *        22 *         8 *         1 *         3 *         0 * -38.80887 *
 *       14 *        22 *         2 *         2 *         1 *         0 * 30.226743 *
 *
 *       15 *        24 *         0 *         2 *         1 *         0 * 105.85574 *
 *       16 *        24 *         0 *         2 *         1 *         0 * 107.73904 *
 * unica coppia compton
 *
 *       17 *        25 *         1 *         1 *         1 *         0 * 93.452529 *
 *       18 *        25 *         1 *         1 *         1 *         1 * 93.679504 *
 *       19 *        25 *         1 *         1 *         1 *         2 * 93.836738 *
 *       20 *        25 *         7 *         2 *         0 *         1 * -93.80313 *
 *       21 *        25 *         7 *         2 *         0 *         1 * -96.04290 *
 *
 *       22 *        27 *         8 *         1 *         0 *         1 * -34.67410 *
 *
 *       23 *        30 *         3 *         1 *         1 *         0 * -21.99345 *
 *       24 *        30 *         3 *         1 *         2 *         0 * -21.85355 *
 *       25 *        30 *         3 *         1 *         2 *         0 * -21.34117 *
 *       26 *        30 *         9 *         2 *         0 *         0 *  22.45331 *
 *
 */


int main(int argc, char* argv[])
{

    ROOT::EnableImplicitMT(); // multi-threading

    bool filter = 0;
    bool analyze = 1;
    bool scatterHisto = 1;
    bool scatterEnergy = 1;
    bool axisCone = 1;

    const char* fileName = argv[1];
    const char* treeName = "filteredPhoton";


    Analysis ansys(filter, analyze, scatterHisto, scatterEnergy, axisCone);
    ansys.Run(fileName, treeName);

    return 0;
}









