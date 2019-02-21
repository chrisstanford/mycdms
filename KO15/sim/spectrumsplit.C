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

// Pileup
double separation_thresh=20.e-6; // 20 us

//poisson_distribution<int> dis_pileup(1./event_rate); // 100 Hz

// int main(int argc, char* argv[]) {  
//   TString infilename = argv[1]; // filename 
//   TFile infile = TFile::Open(infile);
//   TTree* events = (TTree*)infile->Get("events");
//   const int N = events->GetEntries();
// void formatPlot(TH1D* h, TH1D* h_norm, bool integral=false) {
//   h->ResetStats();
//   h->Sumw2();
//   h->SetLineWidth(2);
//   h->SetLineColor(myGreen);
//   if (integral)  h->Scale(1./h_norm->Integral(),"width");
//   else h->Scale(1./h_norm->GetMaximum(),"width");
// }

void formatPlot(TH1D* h, TH1D* h_norm) {
  h->ResetStats();
  h->Sumw2();
  h->SetLineWidth(2);
  h->SetLineColor(myGreen);
  return;
}


int spectrumsplit(TString suffix) {
  const int nbins=500;

  double xbinedges[nbins+1] = {0.,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.,1.02,1.04,1.06,1.08,1.1,1.12,1.14,1.16,1.18,1.2,1.22,1.24,1.26,1.28,1.3,1.32,1.34,1.36,1.38,1.4,1.42,1.44,1.46,1.48,1.5,1.52,1.54,1.56,1.58,1.6,1.62,1.64,1.66,1.68,1.7,1.72,1.74,1.76,1.78,1.8,1.82,1.84,1.86,1.88,1.9,1.92,1.94,1.96,1.98,2.,2.02,2.04,2.06,2.08,2.1,2.12,2.14,2.16,2.18,2.2,2.22,2.24,2.26,2.28,2.3,2.32,2.34,2.36,2.38,2.4,2.42,2.44,2.46,2.48,2.5,2.52,2.54,2.56,2.58,2.6,2.62,2.64,2.66,2.68,2.7,2.72,2.74,2.76,2.78,2.8,2.82,2.84,2.86,2.88,2.9,2.92,2.94,2.96,2.98,3.,3.02,3.04,3.06,3.08,3.1,3.12,3.14,3.16,3.18,3.2,3.22,3.24,3.26,3.28,3.3,3.32,3.34,3.36,3.38,3.4,3.42,3.44,3.46,3.48,3.5,3.52,3.54,3.56,3.58,3.6,3.62,3.64,3.66,3.68,3.7,3.72,3.74,3.76,3.78,3.8,3.82,3.84,3.86,3.88,3.9,3.92,3.94,3.96,3.98,4.,4.02,4.04,4.06,4.08,4.1,4.12,4.14,4.16,4.18,4.2,4.22,4.24,4.26,4.28,4.3,4.32,4.34,4.36,4.38,4.4,4.42,4.44,4.46,4.48,4.5,4.52,4.54,4.56,4.58,4.6,4.62,4.64,4.66,4.68,4.7,4.72,4.74,4.76,4.78,4.8,4.82,4.84,4.86,4.88,4.9,4.92,4.94,4.96,4.98,5.,5.02,5.04,5.06,5.08,5.1,5.12,5.14,5.16,5.18,5.2,5.22,5.24,5.26,5.28,5.3,5.32,5.34,5.36,5.38,5.4,5.42,5.44,5.46,5.48,5.5,5.52,5.54,5.56,5.58,5.6,5.62,5.64,5.66,5.68,5.7,5.72,5.74,5.76,5.78,5.8,5.82,5.84,5.86,5.88,5.9,5.92,5.94,5.96,5.98,6.,6.02,6.04,6.06,6.08,6.1,6.12,6.14,6.16,6.18,6.2,6.22,6.24,6.26,6.28,6.3,6.32,6.34,6.36,6.38,6.4,6.42,6.44,6.46,6.48,6.5,6.52,6.54,6.56,6.58,6.6,6.62,6.64,6.66,6.68,6.7,6.72,6.74,6.76,6.78,6.8,6.82,6.84,6.86,6.88,6.9,6.92,6.94,6.96,6.98,7.,7.02,7.04,7.06,7.08,7.1,7.12,7.14,7.16,7.18,7.2,7.22,7.24,7.26,7.28,7.3,7.32,7.34,7.36,7.38,7.4,7.42,7.44,7.46,7.48,7.5,7.52,7.54,7.56,7.58,7.6,7.62,7.64,7.66,7.68,7.7,7.72,7.74,7.76,7.78,7.8,7.82,7.84,7.86,7.88,7.9,7.92,7.94,7.96,7.98,8.,8.02,8.04,8.06,8.08,8.1,8.12,8.14,8.16,8.18,8.2,8.22,8.24,8.26,8.28,8.3,8.32,8.34,8.36,8.38,8.4,8.42,8.44,8.46,8.48,8.5,8.52,8.54,8.56,8.58,8.6,8.62,8.64,8.66,8.68,8.7,8.72,8.74,8.76,8.78,8.8,8.82,8.84,8.86,8.88,8.9,8.92,8.94,8.96,8.98,9.,9.02,9.04,9.06,9.08,9.1,9.12,9.14,9.16,9.18,9.2,9.22,9.24,9.26,9.28,9.3,9.32,9.34,9.36,9.38,9.4,9.42,9.44,9.46,9.48,9.5,9.52,9.54,9.56,9.58,9.6,9.62,9.64,9.66,9.68,9.7,9.72,9.74,9.76,9.78,9.8,9.82,9.84,9.86,9.88,9.9,9.92,9.94,9.96,9.98,10.};
  // data spectrum
  double y_data_ls_pos[nbins] = {0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 1.0 , 1.0 , 1.0 , 1.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 4.0 , 2.0 , 1.0 , 2.0 , 2.0 , 3.0 , 1.0 , 4.0 , 4.0 , 6.0 , 3.0 , 2.0 , 6.0 , 8.0 , 18.0 , 24.0 , 43.0 , 65.0 , 106.0 , 125.0 , 189.0 , 230.0 , 302.0 , 362.0 , 390.0 , 446.0 , 415.0 , 410.0 , 400.0 , 308.0 , 313.0 , 230.0 , 182.0 , 140.0 , 99.0 , 85.0 , 33.0 , 17.0 , 10.0 , 16.0 , 8.0 , 6.0 , 7.0 , 8.0 , 10.0 , 11.0 , 3.0 , 5.0 , 7.0 , 3.0 , 3.0 , 4.0 , 2.0 , 2.0 , 3.0 , 0.0 , 3.0 , 5.0 , 3.0 , 6.0 , 7.0 , 13.0 , 16.0 , 24.0 , 35.0 , 52.0 , 64.0 , 111.0 , 162.0 , 197.0 , 247.0 , 317.0 , 353.0 , 415.0 , 443.0 , 467.0 , 413.0 , 418.0 , 365.0 , 314.0 , 215.0 , 186.0 , 135.0 , 96.0 , 77.0 , 43.0 , 40.0 , 25.0 , 23.0 , 4.0 , 10.0 , 7.0 , 9.0 , 6.0 , 6.0 , 7.0 , 4.0 , 15.0 , 5.0 , 6.0 , 4.0 , 4.0 , 3.0 , 6.0 , 7.0 , 2.0 , 9.0 , 4.0 , 11.0 , 9.0 , 12.0 , 21.0 , 27.0 , 35.0 , 56.0 , 52.0 , 91.0 , 123.0 , 136.0 , 189.0 , 228.0 , 232.0 , 276.0 , 282.0 , 270.0 , 231.0 , 254.0 , 233.0 , 222.0 , 175.0 , 139.0 , 118.0 , 65.0 , 56.0 , 34.0 , 29.0 , 19.0 , 15.0 , 9.0 , 8.0 , 4.0 , 9.0 , 3.0 , 8.0 , 4.0 , 8.0 , 6.0 , 11.0 , 8.0 , 7.0 , 3.0 , 5.0 , 8.0 , 9.0 , 3.0 , 11.0 , 13.0 , 10.0 , 12.0 , 15.0 , 13.0 , 21.0 , 16.0 , 33.0 , 37.0 , 63.0 , 63.0 , 70.0 , 103.0 , 103.0 , 115.0 , 139.0 , 117.0 , 157.0 , 127.0 , 123.0 , 114.0 , 110.0 , 78.0 , 57.0 , 74.0 , 37.0 , 35.0 , 21.0 , 13.0 , 10.0 , 4.0 , 5.0 , 3.0 , 2.0 , 6.0 , 5.0 , 6.0 , 3.0 , 3.0 , 4.0 , 6.0 , 3.0 , 6.0 , 6.0 , 0.0 , 1.0 , 6.0 , 1.0 , 3.0 , 5.0 , 8.0 , 6.0 , 8.0 , 7.0 , 9.0 , 9.0 , 18.0 , 22.0 , 23.0 , 37.0 , 42.0 , 54.0 , 40.0 , 43.0 , 53.0 , 58.0 , 68.0 , 52.0 , 63.0 , 39.0 , 34.0 , 30.0 , 24.0 , 17.0 , 13.0 , 16.0 , 5.0 , 6.0 , 5.0 , 4.0 , 5.0 , 0.0 , 3.0 , 1.0 , 0.0 , 1.0 , 2.0 , 0.0 , 3.0 , 0.0 , 0.0 , 2.0 , 3.0 , 0.0 , 1.0 , 1.0 , 2.0 , 3.0 , 3.0 , 0.0 , 3.0 , 2.0 , 4.0 , 4.0 , 9.0 , 9.0 , 19.0 , 10.0 , 11.0 , 26.0 , 14.0 , 22.0 , 28.0 , 17.0 , 13.0 , 13.0 , 12.0 , 14.0 , 11.0 , 8.0 , 3.0 , 7.0 , 5.0 , 3.0 , 5.0 , 1.0 , 2.0 , 1.0 , 1.0 , 1.0 , 1.0 , 0.0 , 1.0 , 2.0 , 1.0 , 1.0 , 3.0 , 0.0 , 2.0 , 0.0 , 0.0 , 3.0 , 0.0 , 2.0 , 1.0 , 3.0 , 0.0 , 0.0 , 3.0 , 2.0 , 1.0 , 3.0 , 3.0 , 2.0 , 3.0 , 5.0 , 7.0 , 4.0 , 6.0 , 5.0 , 9.0 , 6.0 , 4.0 , 2.0 , 2.0 , 4.0 , 2.0 , 6.0 , 3.0 , 2.0 , 1.0 , 0.0 , 1.0 , 2.0 , 0.0 , 1.0 , 1.0 , 1.0 , 1.0 , 0.0 , 0.0 , 0.0 , 1.0 , 1.0 , 0.0 , 0.0 , 0.0 , 1.0 , 1.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 1.0 , 0.0 , 2.0 , 0.0 , 2.0 , 2.0 , 0.0 , 1.0 , 0.0 , 1.0 , 0.0 , 2.0 , 1.0 , 1.0 , 1.0 , 1.0 , 2.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0};


  double y_data_ls_neg[nbins] = {0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 2.0 , 1.0 , 1.0 , 2.0 , 4.0 , 5.0 , 5.0 , 5.0 , 12.0 , 29.0 , 52.0 , 90.0 , 126.0 , 258.0 , 373.0 , 521.0 , 655.0 , 812.0 , 898.0 , 1014.0 , 1028.0 , 1025.0 , 915.0 , 806.0 , 574.0 , 440.0 , 324.0 , 231.0 , 141.0 , 95.0 , 63.0 , 33.0 , 17.0 , 18.0 , 12.0 , 11.0 , 8.0 , 3.0 , 13.0 , 10.0 , 3.0 , 5.0 , 4.0 , 6.0 , 10.0 , 7.0 , 7.0 , 9.0 , 7.0 , 5.0 , 7.0 , 7.0 , 1.0 , 10.0 , 12.0 , 10.0 , 22.0 , 31.0 , 48.0 , 84.0 , 132.0 , 203.0 , 321.0 , 419.0 , 565.0 , 685.0 , 843.0 , 958.0 , 1034.0 , 1000.0 , 993.0 , 829.0 , 733.0 , 620.0 , 513.0 , 391.0 , 244.0 , 154.0 , 100.0 , 77.0 , 39.0 , 33.0 , 23.0 , 16.0 , 15.0 , 15.0 , 7.0 , 12.0 , 13.0 , 8.0 , 15.0 , 11.0 , 11.0 , 8.0 , 7.0 , 7.0 , 7.0 , 9.0 , 9.0 , 5.0 , 16.0 , 9.0 , 5.0 , 9.0 , 13.0 , 13.0 , 36.0 , 54.0 , 63.0 , 134.0 , 183.0 , 257.0 , 361.0 , 426.0 , 453.0 , 576.0 , 605.0 , 665.0 , 632.0 , 638.0 , 542.0 , 508.0 , 423.0 , 347.0 , 246.0 , 202.0 , 116.0 , 75.0 , 60.0 , 42.0 , 29.0 , 21.0 , 21.0 , 15.0 , 8.0 , 6.0 , 12.0 , 13.0 , 7.0 , 6.0 , 10.0 , 11.0 , 7.0 , 9.0 , 2.0 , 11.0 , 11.0 , 5.0 , 6.0 , 11.0 , 9.0 , 11.0 , 9.0 , 21.0 , 17.0 , 36.0 , 34.0 , 62.0 , 83.0 , 112.0 , 138.0 , 172.0 , 230.0 , 281.0 , 285.0 , 306.0 , 324.0 , 313.0 , 290.0 , 267.0 , 243.0 , 174.0 , 173.0 , 135.0 , 89.0 , 70.0 , 67.0 , 37.0 , 27.0 , 22.0 , 11.0 , 9.0 , 6.0 , 6.0 , 4.0 , 4.0 , 9.0 , 6.0 , 8.0 , 7.0 , 8.0 , 7.0 , 4.0 , 4.0 , 7.0 , 7.0 , 6.0 , 8.0 , 4.0 , 13.0 , 11.0 , 9.0 , 11.0 , 16.0 , 16.0 , 15.0 , 29.0 , 48.0 , 45.0 , 67.0 , 79.0 , 96.0 , 106.0 , 131.0 , 129.0 , 119.0 , 110.0 , 118.0 , 100.0 , 86.0 , 83.0 , 52.0 , 48.0 , 45.0 , 33.0 , 25.0 , 15.0 , 14.0 , 6.0 , 12.0 , 6.0 , 6.0 , 2.0 , 4.0 , 1.0 , 4.0 , 6.0 , 1.0 , 3.0 , 2.0 , 2.0 , 5.0 , 2.0 , 3.0 , 4.0 , 3.0 , 2.0 , 4.0 , 1.0 , 7.0 , 4.0 , 9.0 , 8.0 , 10.0 , 11.0 , 13.0 , 18.0 , 14.0 , 29.0 , 32.0 , 38.0 , 27.0 , 42.0 , 27.0 , 27.0 , 35.0 , 33.0 , 36.0 , 27.0 , 14.0 , 25.0 , 13.0 , 9.0 , 12.0 , 8.0 , 6.0 , 4.0 , 2.0 , 3.0 , 2.0 , 2.0 , 0.0 , 2.0 , 1.0 , 4.0 , 0.0 , 2.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 2.0 , 0.0 , 4.0 , 2.0 , 2.0 , 0.0 , 0.0 , 2.0 , 2.0 , 3.0 , 5.0 , 5.0 , 11.0 , 8.0 , 12.0 , 12.0 , 11.0 , 10.0 , 7.0 , 10.0 , 11.0 , 5.0 , 11.0 , 13.0 , 8.0 , 6.0 , 5.0 , 3.0 , 8.0 , 6.0 , 2.0 , 4.0 , 1.0 , 1.0 , 1.0 , 1.0 , 2.0 , 1.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 1.0 , 1.0 , 1.0 , 0.0 , 1.0 , 1.0 , 1.0 , 0.0 , 0.0 , 1.0 , 1.0 , 2.0 , 1.0 , 3.0 , 0.0 , 0.0 , 3.0 , 1.0 , 2.0 , 6.0 , 2.0 , 3.0 , 2.0 , 3.0 , 2.0 , 4.0 , 4.0 , 1.0 , 3.0 , 3.0 , 1.0 , 0.0 , 0.0 , 1.0 , 0.0 , 1.0 , 1.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 2.0 , 0.0 , 3.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0};
  
  double y_data_bg_pos[nbins] = {0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 5.0 , 35.0 , 168.0 , 538.0 , 1146.0 , 1643.0 , 1964.0 , 1907.0 , 1611.0 , 1301.0 , 1022.0 , 792.0 , 675.0 , 550.0 , 559.0 , 524.0 , 541.0 , 567.0 , 625.0 , 661.0 , 662.0 , 751.0 , 845.0 , 945.0 , 1067.0 , 1342.0 , 1544.0 , 1905.0 , 2133.0 , 2523.0 , 2824.0 , 2995.0 , 3144.0 , 3325.0 , 3345.0 , 3396.0 , 3366.0 , 3423.0 , 3349.0 , 3446.0 , 3314.0 , 3129.0 , 2914.0 , 2729.0 , 2337.0 , 1976.0 , 1533.0 , 1203.0 , 959.0 , 643.0 , 455.0 , 295.0 , 197.0 , 138.0 , 101.0 , 79.0 , 53.0 , 45.0 , 34.0 , 39.0 , 34.0 , 30.0 , 33.0 , 26.0 , 34.0 , 37.0 , 33.0 , 36.0 , 28.0 , 21.0 , 30.0 , 22.0 , 28.0 , 26.0 , 29.0 , 25.0 , 16.0 , 29.0 , 24.0 , 22.0 , 11.0 , 26.0 , 25.0 , 16.0 , 16.0 , 16.0 , 15.0 , 20.0 , 19.0 , 18.0 , 16.0 , 14.0 , 19.0 , 18.0 , 8.0 , 7.0 , 5.0 , 9.0 , 5.0 , 5.0 , 5.0 , 4.0 , 7.0 , 2.0 , 5.0 , 1.0 , 2.0 , 3.0 , 3.0 , 1.0 , 0.0 , 5.0 , 4.0 , 3.0 , 2.0 , 3.0 , 2.0 , 2.0 , 2.0 , 3.0 , 0.0 , 4.0 , 0.0 , 2.0 , 1.0 , 3.0 , 1.0 , 2.0 , 2.0 , 4.0 , 3.0 , 3.0 , 0.0 , 0.0 , 1.0 , 2.0 , 3.0 , 0.0 , 1.0 , 3.0 , 0.0 , 2.0 , 1.0 , 5.0 , 3.0 , 0.0 , 2.0 , 0.0 , 2.0 , 2.0 , 2.0 , 2.0 , 1.0 , 1.0 , 1.0 , 1.0 , 0.0 , 1.0 , 1.0 , 0.0 , 0.0 , 2.0 , 2.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 1.0 , 1.0 , 0.0 , 0.0 , 0.0 , 1.0 , 2.0 , 1.0 , 2.0 , 1.0 , 0.0 , 1.0 , 1.0 , 0.0 , 0.0 , 0.0 , 1.0 , 2.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 2.0 , 1.0 , 0.0 , 1.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 2.0 , 1.0 , 1.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 2.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0};
  
  double y_data_bg_neg[nbins] = {0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 4.0 , 5.0 , 6.0 , 19.0 , 26.0 , 34.0 , 33.0 , 32.0 , 32.0 , 26.0 , 35.0 , 28.0 , 42.0 , 54.0 , 45.0 , 78.0 , 109.0 , 167.0 , 276.0 , 323.0 , 447.0 , 650.0 , 907.0 , 1140.0 , 1526.0 , 1832.0 , 2171.0 , 2412.0 , 2629.0 , 2630.0 , 2605.0 , 2756.0 , 2752.0 , 2930.0 , 3238.0 , 3597.0 , 4082.0 , 4630.0 , 4854.0 , 5074.0 , 4865.0 , 4451.0 , 4122.0 , 3374.0 , 2729.0 , 2060.0 , 1578.0 , 1086.0 , 714.0 , 524.0 , 319.0 , 193.0 , 138.0 , 108.0 , 84.0 , 80.0 , 59.0 , 52.0 , 50.0 , 36.0 , 37.0 , 45.0 , 37.0 , 36.0 , 38.0 , 33.0 , 38.0 , 30.0 , 36.0 , 36.0 , 35.0 , 32.0 , 19.0 , 26.0 , 27.0 , 31.0 , 28.0 , 24.0 , 23.0 , 19.0 , 23.0 , 25.0 , 22.0 , 20.0 , 32.0 , 25.0 , 28.0 , 38.0 , 31.0 , 19.0 , 31.0 , 29.0 , 23.0 , 23.0 , 7.0 , 17.0 , 15.0 , 14.0 , 9.0 , 7.0 , 5.0 , 6.0 , 5.0 , 5.0 , 3.0 , 4.0 , 6.0 , 4.0 , 4.0 , 8.0 , 3.0 , 2.0 , 2.0 , 4.0 , 5.0 , 3.0 , 4.0 , 3.0 , 3.0 , 3.0 , 3.0 , 2.0 , 2.0 , 6.0 , 2.0 , 2.0 , 2.0 , 0.0 , 5.0 , 0.0 , 3.0 , 2.0 , 5.0 , 4.0 , 4.0 , 2.0 , 3.0 , 5.0 , 4.0 , 7.0 , 2.0 , 6.0 , 4.0 , 3.0 , 3.0 , 1.0 , 2.0 , 2.0 , 3.0 , 0.0 , 2.0 , 0.0 , 3.0 , 3.0 , 3.0 , 3.0 , 2.0 , 3.0 , 3.0 , 1.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 2.0 , 0.0 , 1.0 , 0.0 , 0.0 , 2.0 , 2.0 , 0.0 , 2.0 , 2.0 , 1.0 , 1.0 , 0.0 , 1.0 , 1.0 , 2.0 , 0.0 , 1.0 , 0.0 , 1.0 , 2.0 , 1.0 , 1.0 , 1.0 , 0.0 , 1.0 , 1.0 , 3.0 , 0.0 , 1.0 , 0.0 , 2.0 , 1.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 1.0 , 1.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 2.0 , 0.0 , 2.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 1.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 2.0 , 1.0 , 0.0 , 0.0 , 2.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 1.0 , 1.0 , 0.0 , 1.0 , 0.0 , 1.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 2.0 , 0.0 , 1.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 1.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0};

  // make new bins
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

  // make histograms
  TH1D* h_data_bg_pos = new TH1D("h_data_bg_pos","Background (+140V);e-h pairs",nbins_new-1,xbinedges_new);
  TH1D* h_data_bg_neg = new TH1D("h_data_bg_neg","Background (-140V);e-h pairs",nbins_new-1,xbinedges_new);
  TH1D* h_data_ls_pos = new TH1D("h_data_ls_pos","Laser (+140V);e-h pairs",nbins_new-1,xbinedges_new);
  TH1D* h_data_ls_neg = new TH1D("h_data_ls_neg","Laser (-140V);e-h pairs",nbins_new-1,xbinedges_new);
  
  // rebin
  for (int i=0; i<nbins; i++) {
    int newbin = h_data_bg_pos->FindBin(xbinedges[i]+0.01);
    h_data_bg_pos->SetBinContent(newbin,h_data_bg_pos->GetBinContent(newbin)+y_data_bg_pos[i]);
    h_data_bg_neg->SetBinContent(newbin,h_data_bg_neg->GetBinContent(newbin)+y_data_bg_neg[i]);
    h_data_ls_pos->SetBinContent(newbin,h_data_ls_pos->GetBinContent(newbin)+y_data_ls_pos[i]);
    h_data_ls_neg->SetBinContent(newbin,h_data_ls_neg->GetBinContent(newbin)+y_data_ls_neg[i]);
  }
  formatPlot(h_data_bg_pos,h_data_bg_pos);
  formatPlot(h_data_bg_neg,h_data_bg_neg);
  formatPlot(h_data_ls_pos,h_data_ls_pos);
  formatPlot(h_data_ls_neg,h_data_ls_neg);

  // h_data_bg_pos->Scale(1./h_data_bg_pos->Integral(),"width");
  // h_data_bg_neg->Scale(1./h_data_bg_neg->Integral(),"width");
  // h_data_ls_pos->Scale(1./h_data_ls_pos->Integral(),"width");
  // h_data_ls_neg->Scale(1./h_data_ls_neg->Integral(),"width");

  
  TFile* infile = TFile::Open("simspectrum-"+suffix+".root");
  TTree* events = (TTree*)infile->Get("events");
  double energy;
  int type;
  bool laser;
  int pol;
  double t;
  events->SetBranchAddress("energy",&energy);
  //  events->SetBranchAddress("type",&type);
  events->SetBranchAddress("laser",&laser);
  events->SetBranchAddress("pol",&pol);
  events->SetBranchAddress("t",&t);
  const int N = events->GetEntries();
  //  TH1D* h_sim = new TH1D("hsim","",nbins,0,10);

  TH1D* h_simu_bg_pos = new TH1D("h_simu_bg_pos",";e-h pairs",nbins_new-1,xbinedges_new);
  TH1D* h_simu_bg_neg = new TH1D("h_simu_bg_neg",";e-h pairs",nbins_new-1,xbinedges_new);
  TH1D* h_simu_ls_pos = new TH1D("h_simu_ls_pos",";e-h pairs",nbins_new-1,xbinedges_new);
  TH1D* h_simu_ls_neg = new TH1D("h_simu_ls_neg",";e-h pairs",nbins_new-1,xbinedges_new);

  TH1D* h_simu_bg_pos_pu = new TH1D("h_simu_bg_pos_pu",";e-h pairs",nbins_new-1,xbinedges_new);
  TH1D* h_simu_bg_neg_pu = new TH1D("h_simu_bg_neg_pu",";e-h pairs",nbins_new-1,xbinedges_new);
  TH1D* h_simu_ls_pos_pu = new TH1D("h_simu_ls_pos_pu",";e-h pairs",nbins_new-1,xbinedges_new);
  TH1D* h_simu_ls_neg_pu = new TH1D("h_simu_ls_neg_pu",";e-h pairs",nbins_new-1,xbinedges_new);

  // Event loop
  for (int i=0; i<N; i++) {
    events->GetEntry(i);
    bool L = laser;
    double T = t;
    double E = energy;
    // check for pileup with the next two pulses
    if (i<N-1) {
      events->GetEntry(i+1);
      if (fabs(t-T)<separation_thresh) {
	cout<<setprecision(8)<<"pileup1 "<<T<<" "<<t<<" "<<E<<" "<<energy<<endl;
	E+=energy;
	L = L || laser;
	i++;
	if (i<N-1) {
	  events->GetEntry(i+1);
	  if (fabs(t-T)<separation_thresh) {
	    E+=energy;
	    L = L || laser;
	    i++;
	    cout<<setprecision(8)<<"pileup2 "<<T<<" "<<t<<" "<<energy<<endl;
	    
	  }
	}
	if (L && pol>0) h_simu_ls_pos_pu->Fill(E);
	if (L && pol<0) h_simu_ls_neg_pu->Fill(E);
	if (!L && pol>0) h_simu_bg_pos_pu->Fill(E);
	if (!L && pol<0) h_simu_bg_neg_pu->Fill(E);
      }
    }
    if (L && pol>0) h_simu_ls_pos->Fill(E);
    if (L && pol<0) h_simu_ls_neg->Fill(E);
    if (!L && pol>0) h_simu_bg_pos->Fill(E);
    if (!L && pol<0) h_simu_bg_neg->Fill(E);
    
  }
  
  formatPlot(h_simu_bg_pos_pu,h_simu_bg_pos);
  formatPlot(h_simu_bg_neg_pu,h_simu_bg_neg);
  formatPlot(h_simu_ls_pos_pu,h_simu_ls_pos);
  formatPlot(h_simu_ls_neg_pu,h_simu_ls_neg);

  // h_simu_bg_pos_pu->Scale(h_simu_bg_pos_pu->GetBinWidth(h_simu_bg_pos_pu->FindBin(0.9))*h_data_bg_pos->GetBinContent(h_data_bg_pos->FindBin(0.9))/h_simu_bg_pos->GetBinContent(h_simu_bg_pos->FindBin(0.9)),"width");
  // h_simu_bg_neg_pu->Scale(h_simu_bg_neg_pu->GetBinWidth(h_simu_bg_neg_pu->FindBin(0.9))*h_data_bg_neg->GetBinContent(h_data_bg_neg->FindBin(0.9))/h_simu_bg_neg->GetBinContent(h_simu_bg_neg->FindBin(0.9)),"width");
  // h_simu_ls_pos_pu->Scale(h_simu_ls_pos_pu->GetBinWidth(h_simu_ls_pos_pu->FindBin(0.9))*h_data_ls_pos->GetBinContent(h_data_ls_pos->FindBin(0.9))/h_simu_ls_pos->GetBinContent(h_simu_ls_pos->FindBin(0.9)),"width");
  // h_simu_ls_neg_pu->Scale(h_simu_ls_neg_pu->GetBinWidth(h_simu_ls_neg_pu->FindBin(0.9))*h_data_ls_neg->GetBinContent(h_data_ls_neg->FindBin(0.9))/h_simu_ls_neg->GetBinContent(h_simu_ls_neg->FindBin(0.9)),"width");
  
  formatPlot(h_simu_bg_pos,h_simu_bg_pos);
  formatPlot(h_simu_bg_neg,h_simu_bg_neg);
  formatPlot(h_simu_ls_pos,h_simu_ls_pos);
  formatPlot(h_simu_ls_neg,h_simu_ls_neg);

  cout<<h_data_bg_pos->GetBinContent(h_data_bg_pos->FindBin(0.9))<<endl;
  cout<<h_simu_bg_pos->GetBinContent(h_simu_bg_pos->FindBin(0.9))<<endl;
  // h_simu_bg_pos->Scale(h_simu_bg_pos->GetBinWidth(h_simu_bg_pos->FindBin(0.9))*h_data_bg_pos->GetBinContent(h_data_bg_pos->FindBin(0.9))/h_simu_bg_pos->GetBinContent(h_simu_bg_pos->FindBin(0.9)),"width");
  // h_simu_bg_neg->Scale(h_simu_bg_neg->GetBinWidth(h_simu_bg_neg->FindBin(0.9))*h_data_bg_neg->GetBinContent(h_data_bg_neg->FindBin(0.9))/h_simu_bg_neg->GetBinContent(h_simu_bg_neg->FindBin(0.9)),"width");
  // h_simu_ls_pos->Scale(h_simu_ls_pos->GetBinWidth(h_simu_ls_pos->FindBin(0.9))*h_data_ls_pos->GetBinContent(h_data_ls_pos->FindBin(0.9))/h_simu_ls_pos->GetBinContent(h_simu_ls_pos->FindBin(0.9)),"width");
  // h_simu_ls_neg->Scale(h_simu_ls_neg->GetBinWidth(h_simu_ls_neg->FindBin(0.9))*h_data_ls_neg->GetBinContent(h_data_ls_neg->FindBin(0.9))/h_simu_ls_neg->GetBinContent(h_simu_ls_neg->FindBin(0.9)),"width");
  cout<<h_simu_bg_pos->GetBinContent(h_simu_bg_pos->FindBin(0.9))<<endl;
  
  TCanvas* c = new TCanvas("c","c",1600,1100);
  c->Divide(2,2);
  c->cd(1);
  h_data_ls_pos->Draw();
  h_simu_ls_pos->Draw("same");
  h_simu_ls_pos->SetLineColor(myRed);
  h_simu_ls_pos_pu->Draw("same");
  h_simu_ls_pos_pu->SetLineColor(myOrange);
  gPad->SetLogy();

  TLegend* leg = new TLegend(0.6,0.6,0.85,0.85);
  leg->AddEntry(h_data_bg_pos,"Data","l");
  leg->AddEntry(h_simu_bg_pos,"MC","l");
  leg->AddEntry(h_simu_bg_pos_pu,"Pilup (MC)","l");
  leg->Draw();
  
  c->cd(2);
  h_data_ls_neg->Draw();
  h_simu_ls_neg->Draw("same");
  h_simu_ls_neg->SetLineColor(myRed);
  h_simu_ls_neg_pu->Draw("same");
  h_simu_ls_neg_pu->SetLineColor(myOrange);
  gPad->SetLogy();

  c->cd(3);
  h_data_bg_pos->Draw();
  h_simu_bg_pos->Draw("same");
  h_simu_bg_pos->SetLineColor(myRed);
  h_simu_bg_pos_pu->Draw("same");
  h_simu_bg_pos_pu->SetLineColor(myOrange);
  gPad->SetLogy();

  c->cd(4);
  h_data_bg_neg->Draw();
  h_simu_bg_neg->Draw("same");
  h_simu_bg_neg->SetLineColor(myRed);
  h_simu_bg_neg_pu->Draw("same");
  h_simu_bg_neg_pu->SetLineColor(myOrange);
  gPad->SetLogy();

  // TLegend* leg = new TLegend(0.6,0.6,0.8,0.8);
  // leg->AddEntry(h_data,"Data","l");
  // leg->AddEntry(h_background,"MC","l");
  // leg->Draw();

  gStyle->SetOptStat(0);
  return 0;
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
