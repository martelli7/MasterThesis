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
#include <RtypesCore.h>
#include <TBranch.h>
#include <TObjArray.h>
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
    //.Filter()   //possibile concatenare più filtri
    //.Count();   //Count() non funziona usando anche Snapshot()
    .Snapshot("filteredHits", "filtered_file3.root", {"posX", "posY", "posZ", "photonID", "eventID", "parentID", "nCrystalCompton"});
    //il metodo Snapshot() crea un nuovo file con i filtri e le colonne selezionate, ora il file è molto leggero

    //controlla quante entrate hanno superato il filto
    //std::cout << *df_filtered << " entries passed all filters" << std::endl;

    // Histo2D posX:posY
    /*auto myHist1 = df.Histo2D({"histName", "posX:posY{nCrystalCompton==1}", 64u, -150., 150., 32u, -150., 150.}, "posX", "posY");
     *
     *    TCanvas c;
     *    myHist1->Draw("COLZ"); //scala a colori
     *    c.SaveAs("posXposY.root");
     */



    return 0;
}









