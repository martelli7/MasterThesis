#include <TFile.h>
#include <TTree.h>

// command to compile:  g++ -o copyTree copyTree.cc `root-config --cflags --libs`;

void createEmptyRootFile() {
    // Apri il file esistente
    TFile *inputFile = new TFile("filteredFinal.root", "READ");
    TTree *tree = (TTree*)inputFile->Get("filteredPhoton"); // Sostituisci con il nome effettivo del TTree

    // Crea il nuovo file
    TFile *outputFile = new TFile("output.root", "RECREATE");
    TTree *newTree = tree->CloneTree(0); // Copia solo la struttura senza dati

    // Salva il nuovo file
    newTree->Write();
    outputFile->Close();
    inputFile->Close();

    delete inputFile;
    delete outputFile;
}

int main()
{
    createEmptyRootFile();

    return 0;
}
