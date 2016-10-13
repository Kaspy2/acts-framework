#include "TFile.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TTree.h"

/// This root script reads in two histograms with the names 'hist1Name' and
/// 'hist2Name' from the root file 'inFile' and draws them normalized with
/// different colors in one Canvas

void
compareDistributions(std::string inFile1,
                     std::string hist1Name,
                     std::string inFile2,
                     std::string hist2Name)
{
    std::cout << "Opening file: " << inFile1 << std::endl;
    TFile inputFile1(inFile1.c_str());
    std::cout << "Opening file: " << inFile2 << std::endl;
    TFile inputFile2(inFile2.c_str());
    std::cout << "Comparing Histograms: " << hist1Name << " & " << hist2Name
    << std::endl;
    
    TH1F* h1 = (TH1F*)inputFile1.Get(hist1Name.c_str());
    TH1F* h2 = (TH1F*)inputFile2.Get(hist2Name.c_str());
    
    h1->DrawNormalized();
    h2->SetLineColor(2);
    h2->DrawNormalized("same");
    
    TLegend* leg = new TLegend(0.72, 0.71, 0.99, 0.95);
    leg->AddEntry(h1, hist1Name.c_str());
    leg->AddEntry(h2, hist2Name.c_str());
    leg->Draw();
    
    h1->SetDirectory(0);
    h2->SetDirectory(0);
    
    inputFile1.Close();
    inputFile2.Close();
}
