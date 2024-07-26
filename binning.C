#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>

void log_binned_histogram(){
  
    //TFile *file = TFile::Open("/lustrehome/mbossa/Nuses/Builds/JobsOutput/NUSES_wt_CaloHERDgeantinoPow_1000-10000000/rootOutput/NUSES_wt_CaloHERDgeantinoPow_1000-10000000_10000-evt-0.root");
       TFile *file = TFile::Open("NUSES_wt_CaloHERDgeantinoPow_1000-10000000_10000-evt-0.root");

    TTree *tree = (TTree*)file->Get("Primary");
    

    const char* variable = "PrimaryParticleEnergy";
    int n_bins = 50; 

    double x_min = 1.0; 
    double x_max = 1000.0; 

    double log_min = TMath::Log10(x_min);
    double log_max = TMath::Log10(x_max);
    double bin_edges[n_bins + 1];
    for (int i = 0; i <= n_bins; ++i) {
        bin_edges[i] = TMath::Power(10, log_min + i * (log_max - log_min) / n_bins);
    }

    TH1D *hist = new TH1D("hist", "Istogramma con binning logaritmico;Variabile;Conteggio", n_bins, bin_edges);
    tree->Draw(TString::Format("%s>>hist", variable), "", "goff");

    TCanvas *c = new TCanvas("c", "Istogramma con binning logaritmico", 800, 600);
    c->SetLogx(); // Impostazione dell'asse x in scala logaritmica
    hist->Draw();

    c->SaveAs("log_binned_histogram.png");
    file->Close();
    delete file;
    delete hist;
    delete c;
}
