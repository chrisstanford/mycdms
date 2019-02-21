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

//  gRandom->GetSeed();

std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_uni(0,1);

double UnitStep(double x) {
  if (x>0) return 1;
  else return 0;
}

double i0(double x) {
  return UnitStep(x) - UnitStep(x-1);
}
double i0exp(double* xp, double* par) {
  double x = xp[0];
  double N = par[0];
  double l = par[1];
  return N * l*exp(-l*x) * i0(x);
}
double i0expsurf(double* xp, double* par) {
  return 1;
}
double i1(double x) {
  return ((-2 + x)*UnitStep(-2 + x) - 2*(-1 + x)*UnitStep(-1 + x)+x*UnitStep(x));
}
double i1exp(double* xp, double* par) {
  double x = xp[0];
  double N = par[0];
  double l = par[1];
  return N * l*exp(-l*x) * i1(x);
}
double i1expsurf(double* xp, double* par) {
  double x = xp[0];
  double N = par[0];
  double l = par[1];
  return N*exp(-l*x)*i0(x-1);
}
double i2(double x) {
  return (-(1./2)*pow(-3+x,2)*UnitStep(-3+x)+1./2*pow(-2+x,2)*UnitStep(-2+x)-1./2*pow(-1+x,2)*UnitStep(-1+x)-2*(-(1./2)*pow(-2+x,2)*UnitStep(-2+x)+1./2*pow(-1+x,2)*UnitStep(-1+x))+1./2*pow(x,2)*UnitStep(x));
}
double i2exp(double* xp, double* par) {
  double x = xp[0];
  double N = par[0];
  double l = par[1];
  return N * l*exp(-l*x) * i2(x);
}
double i2expsurf(double* xp, double* par) {
  double x = xp[0];
  double N = par[0];
  double l = par[1];
  return N*exp(-l*x)*i1(x-1);
}
double isum(double* xp, double* par) {
  double x = xp[0];
  double N = par[0];
  double l = par[1]; // mean free path
  return N*exp(-l*x)*(// mu*(1-exp(-1/mu))*
			      i0(x) + // 0 collisions
			      l*i1(x)); // 1 collisions;
  
}


double energyFromCharge(double x, int N_steps, double probII, double &impacts, double maximpacts) {
  // MFP method
  if (x>1) return 0; // sanity check
  if (impacts>maximpacts) return 0.; // Don't bother anymore

  // Optimization attempt (not working)
  // std::uniform_real_distribution<> dis_uni_cut(exp((x-1)*probII),1);
  // double nextCollision;
  // if (impacts<maximpacts) nextCollision= x-log(dis_uni_cut(gen))/probII;
  // else nextCollision= x-log(dis_uni(gen))/probII;

  double nextCollision = x-log(dis_uni(gen))/probII;
  if (nextCollision>1) return 1-x; // charge makes it to 1 without colliding
  else { // impact occured
    impacts++;
    return
      nextCollision-x // energy from this charge until its next collision
      + energyFromCharge(nextCollision,N_steps,probII,impacts,maximpacts) // energy from this charge after its collision
      + energyFromCharge(nextCollision,N_steps,probII,impacts,maximpacts);// energy from impact charge
  }
		       
  // // Step method
  // double E = 1-x; // initial energy  
  // for (double xi=x; xi<1; xi+=1./N_steps) {
  //   // See if you knock out a charge
  //   if (dis_uni(gen)<(probII/N_steps)) {
  //     impacts +=1;
  //     E += energyFromCharge(xi, N_steps, probII, impacts, maximpacts);
  //   }
  //   // Generate only 1 additional charge
  //   // if (dis_uni(gen)<(probII/N_steps)) {
  //   //   impacts +=1;
  //   //   return E + 1-xi;
  //   // }
  // }
  //  return E;
}

int main(int argc, char* argv[]) {
  const int N = atoi(argv[1]); // number of events
  const int NII = atoi(argv[2]); // number of impact ionizations
  const int maxNII = atoi(argv[3]); // max number of impact ionizations to keep track of

  TRandom3* myRNG = new TRandom3();
  gRandom = myRNG;

  int N_steps = 100;
  double hist_max = 4.5;
  const int N_prob = 5;
  double probII[N_prob] = {0.01,0.1,1.,3.,10.}; // probablility of knocking out a charge over the full length of the crystal
  //  double probII[N_prob] = {5.}; // probablility of knocking out a charge over the full length of the crystal
  TH1D* h[N_prob];
  TH1D* hsum[N_prob];
  for (int i=0; i<N_prob; i++) {
    h[i] = new TH1D(Form("h%d",i),";Total distance travelled by charges",N_steps*hist_max/5,0,hist_max);
    hsum[i] = new TH1D(Form("hsum%d",i),";Total distance travelled by charges",N_steps*hist_max/5,0,hist_max);
  }
  // Event loop
  for (int i=0; i<N_prob; i++) {
    cout<<" Generating "<<N<<" events with P="<<probII[i]<<" and taking events with exactly "<<NII<<" impact ionizations"<<endl;
    int Nloop = N;
    if (probII[i]<0.05 && NII>1) Nloop = 10*N;
    for (int j=0; j<Nloop; j++) {
      if (j%1000000==0) cout<<" "<<j<<"/"<<Nloop<<endl;
      // Generate leakage charge position
      double x = dis_uni(gen);
      // Get energy from charge transport
      double impacts = 0;
      //      double E = energyFromCharge(x, N_steps, probII[i], impacts, NII);
      double E = energyFromCharge(x, N_steps, probII[i], impacts, maxNII);
      hsum[i]->Fill(E);
      // Require x number of impacts
      if (impacts!=NII) continue;
      // Fill histogram
      h[i]->Fill(E);
    }
  }
  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c","c",1000,800);
  for (int i=0; i<N_prob; i++) {
    h[i]->Sumw2();
    h[i]->SetLineColor(mycolors12[i]);
    cout<<" "<<h[i]->GetEntries()<<"/"<<N<<" have NII="<<NII<<" for P="<<probII[i]<<endl;
    h[i]->Scale(1./h[i]->Integral()/h[i]->GetBinWidth(1));
    //h[i]->Scale(1./N/h[i]->GetBinWidth(1));
  }
  h[0]->Draw();
  h[0]->SetTitle(Form("Bulk Leakage with at exactly %d I.I.",NII));
  //h[0]->SetTitle("Leakage charges distributed uniformly across the crystal");
  for (int i=1; i<N_prob; i++) {
    h[i]->Draw("same");
  }
  h[0]->GetYaxis()->SetRangeUser(1e-4,1.2*h[N_prob-1]->GetMaximum());

  TF1* f[N_prob];
  for (int i=0; i<N_prob; i++) {
    if (NII==0) f[i] = new TF1(Form("f%d",i),i0exp,0,5,2);
    else if (NII==1) f[i] = new TF1(Form("f%d",i),i1exp,0,5,2);
    else if (NII==2) f[i] = new TF1(Form("f%d",i),i2exp,0,5,2);
    else return 1;
    f[i]->SetParameter(0,1);
    f[i]->FixParameter(1,probII[i]);
    //    f[i]->FixParameter(0,f[i]->GetParameter(0)/f[i]->Integral(0,5));
    h[i]->Fit(f[i],"R0Q");
    cout<<"f["<<i<<"] integral: "<<f[i]->Integral(0,5)<<endl;
    f[i]->SetLineColor(mycolors12[i]);
    f[i]->Draw("same");
  }
  

  //  c->SetLogy();

  TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry("","I.I. Probability","");
  leg->SetNColumns(1);
  for (int i=0; i<N_prob; i++) leg->AddEntry(h[i],Form("%.3f",probII[i]),"le");
  leg->Draw();

  TLegend* leg2 = new TLegend(0.6,0.3,0.9,0.6);
  leg2->AddEntry("","Model","");
  leg2->SetNColumns(1);
  for (int i=0; i<N_prob; i++) leg2->AddEntry(f[i],Form("Chi2/NDF = %.3f",f[i]->GetChisquare()/f[i]->GetNDF()),"l");
  leg2->Draw();

  gStyle->SetOptFit(0);
  c->SaveAs(Form("spectrum_%d_II.pdf",NII));
  c->SaveAs(Form("spectrum_%d_II.png",NII));

  // Sum
  TCanvas *csum = new TCanvas("csum","csum",1000,800);
  csum->cd();
  for (int i=0; i<N_prob; i++) {
    hsum[i]->SetLineColor(mycolors12[i]);
    hsum[i]->Sumw2();
    hsum[i]->Scale(1./(N*hsum[i]->GetBinWidth(1)));
    if (i==0) hsum[i]->Draw();
    else hsum[i]->Draw("same");
  }
  TF1* fsum[N_prob];
  for (int i=0; i<N_prob; i++) {
    fsum[i] = new TF1(Form("fsum%d",i),isum,0,2,2);
    fsum[i]->FixParameter(0,1);
    fsum[i]->FixParameter(1,probII[i]);
    //    f[i]->FixParameter(0,f[i]->GetParameter(0)/f[i]->Integral(0,5));
    h[i]->Fit(fsum[i],"R0");
    fsum[i]->SetLineColor(mycolors12[i]);
    fsum[i]->Draw("same");
    cout<<fsum[i]->Eval(0.5)<<endl;
  }
  TLegend* leg3 = new TLegend(0.6,0.3,0.9,0.6);
  leg3->AddEntry("","Model","");
  leg3->SetNColumns(1);
  for (int i=0; i<N_prob; i++) leg3->AddEntry(fsum[i],Form("Chi2/NDF = %.3f",fsum[i]->GetChisquare()/fsum[i]->GetNDF()),"l");
  leg3->Draw();
  leg->Draw();
  csum->SetLogy();
  //  hsum[0]->GetYaxis()->SetRangeUser(0.9,1.1);
  csum->SaveAs("spectrum_sum.png");
}
