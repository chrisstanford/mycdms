// std includes
#include <algorithm>
#include <iterator>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <ctime>
#include <map>

// ROOT includes
#include "TSystem.h"
#include "TLegend.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPaveText.h"
#include "TFitResult.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TDirectory.h"
#include "TPaletteAxis.h"
#include "TCutG.h"

#include "mycolors.h"

using namespace std;

// std::random_device rd;  //Will be used to obtain a seed for the random number engine
// std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
// std::uniform_real_distribution<> dis_uni(0,1);
// Pileup
double event_rate=100; // 100 Hz
double separation_thresh=0.e-6; // 20 us
double prob_pileup=event_rate*separation_thresh;

//poisson_distribution<int> dis_pileup(1./event_rate); // 100 Hz

// int main(int argc, char* argv[]) {  
//   TString infilename = argv[1]; // filename 
//   TFile infile = TFile::Open(infile);
//   TTree* events = (TTree*)infile->Get("events");
//   const int N = events->GetEntries();
int spectrumdraw(TString suffix) {
  const int nbins=500;

  double xbinedges[nbins+1] = {0.,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.,1.02,1.04,1.06,1.08,1.1,1.12,1.14,1.16,1.18,1.2,1.22,1.24,1.26,1.28,1.3,1.32,1.34,1.36,1.38,1.4,1.42,1.44,1.46,1.48,1.5,1.52,1.54,1.56,1.58,1.6,1.62,1.64,1.66,1.68,1.7,1.72,1.74,1.76,1.78,1.8,1.82,1.84,1.86,1.88,1.9,1.92,1.94,1.96,1.98,2.,2.02,2.04,2.06,2.08,2.1,2.12,2.14,2.16,2.18,2.2,2.22,2.24,2.26,2.28,2.3,2.32,2.34,2.36,2.38,2.4,2.42,2.44,2.46,2.48,2.5,2.52,2.54,2.56,2.58,2.6,2.62,2.64,2.66,2.68,2.7,2.72,2.74,2.76,2.78,2.8,2.82,2.84,2.86,2.88,2.9,2.92,2.94,2.96,2.98,3.,3.02,3.04,3.06,3.08,3.1,3.12,3.14,3.16,3.18,3.2,3.22,3.24,3.26,3.28,3.3,3.32,3.34,3.36,3.38,3.4,3.42,3.44,3.46,3.48,3.5,3.52,3.54,3.56,3.58,3.6,3.62,3.64,3.66,3.68,3.7,3.72,3.74,3.76,3.78,3.8,3.82,3.84,3.86,3.88,3.9,3.92,3.94,3.96,3.98,4.,4.02,4.04,4.06,4.08,4.1,4.12,4.14,4.16,4.18,4.2,4.22,4.24,4.26,4.28,4.3,4.32,4.34,4.36,4.38,4.4,4.42,4.44,4.46,4.48,4.5,4.52,4.54,4.56,4.58,4.6,4.62,4.64,4.66,4.68,4.7,4.72,4.74,4.76,4.78,4.8,4.82,4.84,4.86,4.88,4.9,4.92,4.94,4.96,4.98,5.,5.02,5.04,5.06,5.08,5.1,5.12,5.14,5.16,5.18,5.2,5.22,5.24,5.26,5.28,5.3,5.32,5.34,5.36,5.38,5.4,5.42,5.44,5.46,5.48,5.5,5.52,5.54,5.56,5.58,5.6,5.62,5.64,5.66,5.68,5.7,5.72,5.74,5.76,5.78,5.8,5.82,5.84,5.86,5.88,5.9,5.92,5.94,5.96,5.98,6.,6.02,6.04,6.06,6.08,6.1,6.12,6.14,6.16,6.18,6.2,6.22,6.24,6.26,6.28,6.3,6.32,6.34,6.36,6.38,6.4,6.42,6.44,6.46,6.48,6.5,6.52,6.54,6.56,6.58,6.6,6.62,6.64,6.66,6.68,6.7,6.72,6.74,6.76,6.78,6.8,6.82,6.84,6.86,6.88,6.9,6.92,6.94,6.96,6.98,7.,7.02,7.04,7.06,7.08,7.1,7.12,7.14,7.16,7.18,7.2,7.22,7.24,7.26,7.28,7.3,7.32,7.34,7.36,7.38,7.4,7.42,7.44,7.46,7.48,7.5,7.52,7.54,7.56,7.58,7.6,7.62,7.64,7.66,7.68,7.7,7.72,7.74,7.76,7.78,7.8,7.82,7.84,7.86,7.88,7.9,7.92,7.94,7.96,7.98,8.,8.02,8.04,8.06,8.08,8.1,8.12,8.14,8.16,8.18,8.2,8.22,8.24,8.26,8.28,8.3,8.32,8.34,8.36,8.38,8.4,8.42,8.44,8.46,8.48,8.5,8.52,8.54,8.56,8.58,8.6,8.62,8.64,8.66,8.68,8.7,8.72,8.74,8.76,8.78,8.8,8.82,8.84,8.86,8.88,8.9,8.92,8.94,8.96,8.98,9.,9.02,9.04,9.06,9.08,9.1,9.12,9.14,9.16,9.18,9.2,9.22,9.24,9.26,9.28,9.3,9.32,9.34,9.36,9.38,9.4,9.42,9.44,9.46,9.48,9.5,9.52,9.54,9.56,9.58,9.6,9.62,9.64,9.66,9.68,9.7,9.72,9.74,9.76,9.78,9.8,9.82,9.84,9.86,9.88,9.9,9.92,9.94,9.96,9.98,10.};
  // data spectrum
  double data[nbins] = {0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,6.0e+00,1.5e+01,2.5e+01,5.7e+01,1.2e+02,1.430e+02,2.150e+02,3.160e+02,3.580e+02,4.420e+02,5.360e+02,6.170e+02,7.520e+02,9.450e+02,1.126e+03,1.355e+03,1.615e+03,1.906e+03,2.323e+03,2.747e+03,3.065e+03,3.524e+03,4.068e+03,4.594e+03,5.063e+03,5.667e+03,5.945e+03,6.318e+03,6.661e+03,6.697e+03,6.841e+03,6.985e+03,7.076e+03,7.428e+03,7.633e+03,8.498e+03,8.762e+03,9.365e+03,9.746e+03,9.650e+03,9.372e+03,9.042e+03,8.046e+03,7.240e+03,5.861e+03,4.718e+03,3.710e+03,2.741e+03,2.084e+03,1.420e+03,9.750e+02,6.820e+02,4.870e+02,3.510e+02,2.620e+02,1.940e+02,1.780e+02,1.530e+02,1.410e+02,1.130e+02,1.050e+02,1.260e+02,1.120e+02,1.060e+02,1.120e+02,9.5e+01,1.080e+02,8.8e+01,8.7e+01,8.3e+01,8.4e+01,8.2e+01,7.4e+01,9.4e+01,6.5e+01,7.5e+01,8.7e+01,6.5e+01,6.6e+01,7.5e+01,6.3e+01,5.8e+01,7.2e+01,7.4e+01,6.5e+01,6.5e+01,5.6e+01,6.9e+01,4.1e+01,8.4e+01,8.6e+01,5.9e+01,5.0e+01,4.2e+01,5.5e+01,4.5e+01,3.3e+01,3.8e+01,2.3e+01,3.3e+01,2.6e+01,1.6e+01,2.5e+01,1.5e+01,1.7e+01,1.1e+01,1.1e+01,1.3e+01,1.6e+01,1.1e+01,1.3e+01,1.4e+01,1.0e+01,1.2e+01,1.1e+01,1.2e+01,1.7e+01,8.0e+00,8.0e+00,1.4e+01,8.0e+00,1.4e+01,1.2e+01,6.0e+00,6.0e+00,9.0e+00,9.0e+00,7.0e+00,6.0e+00,2.0e+00,1.5e+01,9.0e+00,8.0e+00,1.5e+01,8.0e+00,7.0e+00,3.0e+00,1.2e+01,9.0e+00,1.7e+01,7.0e+00,6.0e+00,7.0e+00,8.0e+00,9.0e+00,6.0e+00,1.0e+01,8.0e+00,5.0e+00,7.0e+00,1.0e+00,6.0e+00,5.0e+00,6.0e+00,7.0e+00,4.0e+00,3.0e+00,7.0e+00,2.0e+00,6.0e+00,4.0e+00,3.0e+00,4.0e+00,7.0e+00,2.0e+00,2.0e+00,0.0e+00,6.0e+00,0.0e+00,5.0e+00,1.0e+00,4.0e+00,4.0e+00,5.0e+00,3.0e+00,2.0e+00,5.0e+00,2.0e+00,1.0e+00,1.0e+00,3.0e+00,2.0e+00,4.0e+00,3.0e+00,2.0e+00,1.0e+00,1.0e+00,3.0e+00,1.0e+00,2.0e+00,4.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,3.0e+00,3.0e+00,2.0e+00,0.0e+00,2.0e+00,0.0e+00,1.0e+00,2.0e+00,5.0e+00,2.0e+00,1.0e+00,0.0e+00,3.0e+00,0.0e+00,0.0e+00,2.0e+00,1.0e+00,1.0e+00,0.0e+00,2.0e+00,0.0e+00,2.0e+00,1.0e+00,3.0e+00,2.0e+00,3.0e+00,0.0e+00,1.0e+00,1.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,2.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,1.0e+00,1.0e+00,3.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,1.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,1.0e+00,1.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,2.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,1.0e+00,0.0e+00,1.0e+00,2.0e+00,1.0e+00,2.0e+00,0.0e+00,1.0e+00,1.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,1.0e+00,2.0e+00,0.0e+00,2.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,1.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,1.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,1.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,1.0e+00,1.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+0};

  const int nbins_new=154;
  double xbinedges_new[nbins_new];
  double x=0;
  for (int i=0; i<nbins_new; i++) {
    xbinedges_new[i] = x;
    if (x>=0) x+=0.02;
    if (x>=1) x+=0.02;
    if (x>=2) x+=0.02;
    if (x>=3) x+=0.02;
    if (x>=4) x+=0.02;
    if (x>=5) x+=0.02;
    //    cout<<i<<" "<<x<<endl;
  }

  // TH1D* h_data = new TH1D("h_data","",nbins,xbinedges);
  // for (int i=0; i<nbins; i++) h_data->SetBinContent(i,data[i]);
  
  // rebin
  TH1D* h_data = new TH1D("h_data",";e-h pairs",nbins_new-1,xbinedges_new);
  for (int i=0; i<nbins; i++) {
    int newbin = h_data->FindBin(xbinedges[i]+0.01);
    h_data->SetBinContent(newbin,h_data->GetBinContent(newbin)+data[i]);
  }
  h_data->ResetStats();
  h_data->Sumw2();
  h_data->SetLineWidth(2);
  h_data->Scale(1./h_data->Integral(),"width");

  //  TFile* infile = TFile::Open("simspectrum-run1.root");
  TFile* infile = TFile::Open("simspectrum-"+suffix+".root");
  TTree* events = (TTree*)infile->Get("events");
  double energy;
  events->SetBranchAddress("energy",&energy);
  const int N = events->GetEntries();
  //  TH1D* h_sim = new TH1D("hsim","",nbins,0,10);
  TH1D* h_sim = new TH1D("hsim","",nbins_new-1,xbinedges_new);

  TRandom3* rand = new TRandom3();
  rand->SetSeed();
  
  int n_pileup=0;
  // Event loop
  for (int i=0; i<N; i++) {
    events->GetEntry(i);
    // Check for pileup
    if (rand->Uniform(1)<prob_pileup) { // pileup
      double energy1 = energy;
      n_pileup++;
      i++;
      //      cout<<"Pileup "<<n_pileup<<"/"<<i<<" ("<<float(n_pileup)/i<<")"<<endl;
      if (i<N) {
	events->GetEntry(i);
	h_sim->Fill(energy1+energy);
      }
    } else { // no pileup
      h_sim->Fill(energy);
    }
  }
  h_sim->Sumw2();
  h_sim->SetLineWidth(2);
  h_sim->Scale(1./h_sim->Integral(),"width");

  TCanvas* c = new TCanvas("c","c",1600,1100);
  c->cd();
  h_data->Draw();
  h_sim->Draw("same");
  h_sim->SetLineColor(kRed);
  c->SetLogy();

  TLegend* leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->AddEntry(h_data,"Data","l");
  leg->AddEntry(h_sim,"MC","l");
  leg->Draw();

  gStyle->SetOptStat(0);
}

// gStyle->SetOptStat(0);
//   TCanvas *c = new TCanvas("c","c",1000,800);
//   for (int i=0; i<N_prob; i++) {
//     h[i]->Sumw2();
//     h[i]->SetLineColor(mycolors12[i]);
//     cout<<" "<<h[i]->GetEntries()<<"/"<<N<<" have NII="<<NII<<" for P="<<probII[i]<<endl;
//     h[i]->Scale(1./h[i]->Integral()/h[i]->GetBinWidth(1));
//     //h[i]->Scale(1./N/h[i]->GetBinWidth(1));
//   }
//   h[0]->Draw();
//   h[0]->SetTitle(Form("Bulk Leakage with at exactly %d I.I.",NII));
//   //h[0]->SetTitle("Leakage charges distributed uniformly across the crystal");
//   for (int i=1; i<N_prob; i++) {
//     h[i]->Draw("same");
//   }
//   h[0]->GetYaxis()->SetRangeUser(1e-4,1.2*h[N_prob-1]->GetMaximum());

//   TF1* f[N_prob];
//   for (int i=0; i<N_prob; i++) {
//     if (NII==0) f[i] = new TF1(Form("f%d",i),i0exp,0,5,2);
//     else if (NII==1) f[i] = new TF1(Form("f%d",i),i1exp,0,5,2);
//     else if (NII==2) f[i] = new TF1(Form("f%d",i),i2exp,0,5,2);
//     else return 1;
//     f[i]->SetParameter(0,1);
//     f[i]->FixParameter(1,probII[i]);
//     //    f[i]->FixParameter(0,f[i]->GetParameter(0)/f[i]->Integral(0,5));
//     h[i]->Fit(f[i],"R0Q");
//     cout<<"f["<<i<<"] integral: "<<f[i]->Integral(0,5)<<endl;
//     f[i]->SetLineColor(mycolors12[i]);
//     f[i]->Draw("same");
//   }
  

//   //  c->SetLogy();

//   TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
//   leg->AddEntry("","I.I. Probability","");
//   leg->SetNColumns(1);
//   for (int i=0; i<N_prob; i++) leg->AddEntry(h[i],Form("%.3f",probII[i]),"le");
//   leg->Draw();

//   TLegend* leg2 = new TLegend(0.6,0.3,0.9,0.6);
//   leg2->AddEntry("","Model","");
//   leg2->SetNColumns(1);
//   for (int i=0; i<N_prob; i++) leg2->AddEntry(f[i],Form("Chi2/NDF = %.3f",f[i]->GetChisquare()/f[i]->GetNDF()),"l");
//   leg2->Draw();

//   gStyle->SetOptFit(0);
//   c->SaveAs(Form("spectrum_%d_II.pdf",NII));
//   c->SaveAs(Form("spectrum_%d_II.png",NII));

//   // Sum
//   TCanvas *csum = new TCanvas("csum","csum",1000,800);
//   csum->cd();
//   for (int i=0; i<N_prob; i++) {
//     hsum[i]->SetLineColor(mycolors12[i]);
//     hsum[i]->Sumw2();
//     hsum[i]->Scale(1./(N*hsum[i]->GetBinWidth(1)));
//     if (i==0) hsum[i]->Draw();
//     else hsum[i]->Draw("same");
//   }
//   TF1* fsum[N_prob];
//   for (int i=0; i<N_prob; i++) {
//     fsum[i] = new TF1(Form("fsum%d",i),isum,0,2,2);
//     fsum[i]->FixParameter(0,1);
//     fsum[i]->FixParameter(1,probII[i]);
//     //    f[i]->FixParameter(0,f[i]->GetParameter(0)/f[i]->Integral(0,5));
//     h[i]->Fit(fsum[i],"R0");
//     fsum[i]->SetLineColor(mycolors12[i]);
//     fsum[i]->Draw("same");
//     cout<<fsum[i]->Eval(0.5)<<endl;
//   }
//   TLegend* leg3 = new TLegend(0.6,0.3,0.9,0.6);
//   leg3->AddEntry("","Model","");
//   leg3->SetNColumns(1);
//   for (int i=0; i<N_prob; i++) leg3->AddEntry(fsum[i],Form("Chi2/NDF = %.3f",fsum[i]->GetChisquare()/fsum[i]->GetNDF()),"l");
//   leg3->Draw();
//   leg->Draw();
//   csum->SetLogy();
//   //  hsum[0]->GetYaxis()->SetRangeUser(0.9,1.1);
//   csum->SaveAs("spectrum_sum.png");
//   return 0;
// }
