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
#include <limits>

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
#include "spectrum.h"
#include "mycolors.h"

using namespace std;


int setVariable(string varFile, string varName, string varVal) {
  ifstream filein(varFile); //File to read from
  ofstream fileout("_fileout.sty"); //Temporary file
  if(!filein || !fileout) {
    cout << "Error opening files!" << endl;
    return 1;
  }
  string strTemp;
  bool found = false;
  while(filein >> strTemp) {
    string varNameTemp = strTemp.substr(0, strTemp.find('=')-1); 
    if(varNameTemp == varName) {
      fileout << varName << "=" << varVal << "\n";
      found = true;
    } else {
      fileout << strTemp << "\n";
    }
    //if(found) break; // Update only first occurance
   }
  if (!found) {
    fileout << varName << "=" << varVal << "\n";
  }
  rename("_fileout.sty",varFile.c_str());
  return 0;
}

string getVariable(string varFile, string varName) {
  ifstream filein(varFile); //File to read from
  if(!filein) {
    cout << "Error opening file!" << endl;
    return 0;
  }
  string strTemp;
  bool readNextLine = false;
  while(filein >> strTemp) {
    string varNameTemp = strTemp.substr(0, strTemp.find('='));     
    //    cout<<varNameTemp<<" "<<varName<<endl;
    if(varNameTemp == varName) {
      cout<<varName<<" ";
      int i1 = strTemp.find('=')+1;
      int i2 = min(strTemp.size(),min(strTemp.find(' '),strTemp.find('#')));
      cout<<"|"<<strTemp.substr(i1,i2)<<"|"<<endl;
      return strTemp.substr(i1,i2);
    }
  }
  return "0";
}

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



std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_uni(0,1);

double Spectrum::EnergyFromCharge(double x, int sign) {
  if (num_impacts_ee>max_impacts) return 0.; // Don't bother anymore
  if (num_impacts_eh>max_impacts) return 0.; // Don't bother anymore
  if (num_impacts_he>max_impacts) return 0.; // Don't bother anymore
  if (num_impacts_hh>max_impacts) return 0.; // Don't bother anymore
  if (x>1) return 0.; // sanity check
  if (x<0) return 0.; // sanity check

  // Start with values beyond bounds of the crystal (0->1)
  double delta_impact_ee = numeric_limits<double>::max(); // distance between where this e started and the place it bumps out another e
  double delta_impact_eh = numeric_limits<double>::max();
  double delta_impact_he = numeric_limits<double>::max();
  double delta_impact_hh = numeric_limits<double>::max();
  double delta_trap_e = numeric_limits<double>::max();
  double delta_trap_h = numeric_limits<double>::max();

  // Calculate the next interaction
  if (lambda_impact_ee>0 && sign<0) delta_impact_ee = -log(dis_uni(gen))/lambda_impact_ee;
  if (lambda_impact_eh>0 && sign<0) delta_impact_eh = -log(dis_uni(gen))/lambda_impact_eh;
  if (lambda_impact_he>0 && sign>0) delta_impact_he = -log(dis_uni(gen))/lambda_impact_he;
  if (lambda_impact_hh>0 && sign>0) delta_impact_hh = -log(dis_uni(gen))/lambda_impact_hh;
  if (lambda_trap_e>0 && sign<0) delta_trap_e = -log(dis_uni(gen))/lambda_trap_e;
  if (lambda_trap_h>0 && sign>0) delta_trap_h = -log(dis_uni(gen))/lambda_trap_h;

  const int num_interaction_types = 6;
  double deltas[num_interaction_types] = {delta_impact_ee,delta_impact_eh,delta_impact_he,delta_impact_hh,delta_trap_e,delta_trap_h};
  // Find earliest interaction
  int min_index = 0;
  for(int i=1; i<num_interaction_types; i++) {
    if(fabs(deltas[i]) < fabs(deltas[min_index]))
      min_index = i;     
  }
  // Compute next step. The smallest delta determines the process that will occur.
  if (sign<0 && x+deltas[min_index]>1) return 1-x;// charge makes it to edge without colliding
  else if (sign>0 && x-deltas[min_index]<0) return x;// charge makes it to edge without colliding
  else if (min_index==0 && sign<0) { // e impact ionizes and produces e
    num_impacts_ee++;
    return
      delta_impact_ee // energy from this charge until its impact
      + EnergyFromCharge(x+delta_impact_ee,-1) // energy from this charge after its impact
      + EnergyFromCharge(x+delta_impact_ee,-1);// energy from impact charge
  } else if (min_index==1 && sign<0) {// e impact ionizes and produces h
    num_impacts_eh++;
    return
      delta_impact_eh // energy from this charge until its impact
      + EnergyFromCharge(x+delta_impact_eh,-1) // energy from this charge after its impact
      + EnergyFromCharge(x+delta_impact_eh,+1);// energy from impact charge
  } else if (min_index==2 && sign>0) {// h impact ionizes and produces e
    num_impacts_he++;
    return
      delta_impact_he // energy from this charge until its impact
      + EnergyFromCharge(x-delta_impact_he,+1) // energy from this charge after its impact
      + EnergyFromCharge(x-delta_impact_he,-1);// energy from impact charge
  } else if (min_index==3 && sign>0) {// h impact ionizes and produces h
    num_impacts_hh++;
    return
      delta_impact_hh // energy from this charge until its impact
      + EnergyFromCharge(x-delta_impact_ee,+1) // energy from this charge after its impact
      + EnergyFromCharge(x-delta_impact_ee,+1);// energy from impact charge
  } else if (min_index==4 && sign<0) {// e gets trapped
    if (x+delta_trap_e>1) return 1-x; // charge makes it to edge without trapping
    else {
      num_traps_e++;
      return delta_trap_e; // energy from this charge until its trap
    }    
  } else if (min_index==5 && sign>0) {// h gets trapped
    if (x-delta_trap_e<0) return x; // charge makes it to edge without trapping
    else {
      num_traps_h++;
      return delta_trap_h; // energy from this charge until its trap
    }    
  } else if (sign<0) {// these last two should never actually happen
    return 1-x;
  } else {
    return x;
  } 
  return 0.;
}

void Spectrum::Reset() {
  energy=0.;
  start_position=0.;
  num_impacts_ee=0;
  num_impacts_eh=0;
  num_impacts_he=0;
  num_impacts_hh=0;
  num_traps_e=0;
  num_traps_h=0;
  return;
}

Spectrum::Spectrum(string cfg) {
  LoadConfig(cfg);
}

void Spectrum::LoadConfig(string cfg) {
  lambda_impact_ee = stod(getVariable(cfg,"lambda_impact_ee"));
  lambda_impact_eh = stod(getVariable(cfg,"lambda_impact_eh"));
  lambda_impact_he = stod(getVariable(cfg,"lambda_impact_he"));
  lambda_impact_hh = stod(getVariable(cfg,"lambda_impact_hh"));
  lambda_trap_e = stod(getVariable(cfg,"lambda_trap_e"));
  lambda_trap_h = stod(getVariable(cfg,"lambda_trap_h"));
  lambda_bulk_e_possion = stod(getVariable(cfg,"lambda_bulk_e_poisson"));
  lambda_surf_e_possion = stod(getVariable(cfg,"lambda_surf_e_poisson"));
  frac_surf_e = stod(getVariable(cfg,"frac_surf_e"));
  frac_surf_h = stod(getVariable(cfg,"frac_surf_h"));
  frac_bulk_e = stod(getVariable(cfg,"frac_bulk_e"));
  frac_bulk_e_poisson = stod(getVariable(cfg,"frac_bulk_e_poisson"));
  frac_bulk_h = stod(getVariable(cfg,"frac_bulk_h"));
  resolution = stod(getVariable(cfg,"resolution"));
  thresh = stod(getVariable(cfg,"thresh"));
  thresh_res = stod(getVariable(cfg,"thresh_res"));
  max_impacts = stoi(getVariable(cfg,"max_impacts"));
  suffix = getVariable(cfg,"suffix").c_str();
  frac_laser = stod(getVariable(cfg,"frac_laser"));
  bin_n_laser = stoi(getVariable(cfg,"bin_n_laser"));
  bin_p_laser = stod(getVariable(cfg,"bin_p_laser"));
  
}

// double energyFromCharge(double x, double lambda_impact, double lambda_trap, int &num_impacts, int &num_traps, int max_impacts) {
//   if (num_impacts>max_impacts) return 0.; // Don't bother anymore
//   if (x>1) return 0; // sanity check

//   // MFP method
//   double next_impact = x-log(dis_uni(gen))/lambda_impact;
//   double next_trap = x-log(dis_uni(gen))/lambda_trap;
//   if (next_impact<next_trap) { // charge impact ionizes before being trapped.
//     if (next_impact>1) return 1-x; // charge makes it to 1 without colliding
//     else { // impact occured
//       num_impacts++;
//       return
// 	next_impact-x // energy from this charge until its impact
// 	+ energyFromCharge(next_impact,lambda_impact,lambda_trap,num_impacts,num_traps,max_impacts) // energy from this charge after its impact
// 	+ energyFromCharge(next_impact,lambda_impact,lambda_trap,num_impacts,num_traps,max_impacts);// energy from impact charge
//     }
//   } else { // charge is trapped before impact ionizing
//     if (next_trap>1) return 1-x; // charge makes it to 1 without trapping
//     else {
//       num_traps++;
//       return next_trap-x; // energy from this charge until its trap
//     }
//   }
  
//   // // Step method (impact only)
//   // double E = 1-x; // initial energy  
//   // const int N_steps = 100+100*lambda_impact;
//   // for (double xi=x; xi<1; xi+=1./N_steps) {
//   //   // See if you knock out a charge
//   //   if (dis_uni(gen)<(lambda_impact/N_steps)) {
//   //     num_impacts +=1;
//   //     E += energyFromCharge(xi, lambda_impact, lambda_trap, num_impacts, num_traps, max_impacts);
//   //   }
//   //   // Generate only 1 additional charge
//   //   // if (dis_uni(gen)<(lambda_impact/N_steps)) {
//   //   //   impacts +=1;
//   //   //   return E + 1-xi;
//   //   // }
//   // }
//   // return E;
// }

void Spectrum::Run(int N) {
  TRandom3* myRNG = new TRandom3();
  gRandom = myRNG;
  std::normal_distribution<> dis_res(1,resolution);
  std::normal_distribution<> dis_res_0(0,resolution);
  std::normal_distribution<> dis_thresh(thresh,thresh_res);

  TString outfilename = Form("simspectrum-%s.root",suffix.c_str());
  TFile* outfile = new TFile(outfilename,"RECREATE");
  TTree* outtree = new TTree("events","events");
  int event_type = 0; // 0 = leakage, 1 = possion background
  outtree->Branch("energy",&energy);
  outtree->Branch("num_impacts_ee",&num_impacts_ee);
  outtree->Branch("num_impacts_eh",&num_impacts_eh);
  outtree->Branch("num_impacts_he",&num_impacts_he);
  outtree->Branch("num_impacts_hh",&num_impacts_hh);
  outtree->Branch("num_traps_e",&num_traps_e);
  outtree->Branch("num_traps_h",&num_traps_h);
  outtree->Branch("start_position",&start_position);
  outtree->Branch("event_type",&event_type);

  // Event loops

  // surf e
  poisson_distribution<int> dis_surf_e(lambda_surf_e_possion);
  event_type = 1;
  for (int i=0; i<N*frac_surf_e; i++) {
    if (i%1000000==0) cout<<" "<<i<<"/"<<N*frac_surf_e<<endl;
    // Reset variables
    Reset();
    // Generate leakage charge position
    start_position = 0.;
    // Get number of charges
    int n_charges=1+dis_surf_e(gen);
    // Get energy from charge transport
    for (int c=0; c<n_charges; c++)
      energy += EnergyFromCharge(start_position, -1);
    // Apply resolution
    energy += dis_res_0(gen);
    //    energy *= dis_res(gen);
    // Apply threshhold
    if (energy<dis_thresh(gen)) continue;
    // Fill
    outtree->Fill();
  }
  // surf h
  event_type = 2;
  for (int i=0; i<N*frac_surf_h; i++) {
    if (i%1000000==0) cout<<" "<<i<<"/"<<N*frac_surf_h<<endl;
    // Reset variables
    Reset();
    // Generate leakage charge position
    start_position = 1.;
    // Get energy from charge transport
    energy = EnergyFromCharge(start_position, +1);
    // Apply resolution
    energy += dis_res_0(gen);
    //    energy *= dis_res(gen);
    // Apply threshhold
    if (energy<dis_thresh(gen)) continue;
    // Fill
    outtree->Fill();
  }
  // bulk e
  event_type = 3;
  for (int i=0; i<N*frac_bulk_e; i++) {
    if (i%1000000==0) cout<<" "<<i<<"/"<<N*frac_bulk_e<<endl;
    // Reset variables
    Reset();
    // Generate leakage charge position
    start_position = dis_uni(gen);
    //    cout<<start_position<<endl;
    // Get energy from charge transport
    energy = EnergyFromCharge(start_position, -1);
    // Apply resolution
    energy += dis_res_0(gen);
    //    energy *= dis_res(gen);
    // Apply threshhold
    if (energy<dis_thresh(gen)) continue;
    // Fill
    outtree->Fill();
  }
  // bulk h
  event_type = 4;
  for (int i=0; i<N*frac_bulk_h; i++) {
    if (i%1000000==0) cout<<" "<<i<<"/"<<N*frac_bulk_h<<endl;
    // Reset variables
    Reset();
    // Generate leakage charge position
    start_position = dis_uni(gen);
    // Get energy from charge transport
    energy = EnergyFromCharge(start_position, +1);
    // Apply resolution
    energy += dis_res_0(gen);
    //    energy *= dis_res(gen);
    // Apply threshhold
    if (energy<dis_thresh(gen)) continue;
    // Fill
    outtree->Fill();
  }
  // poisson backgound 
  event_type = 5;
  poisson_distribution<int> dis_background(lambda_bulk_e_possion);
  for (int i=0; i<N*frac_bulk_e_poisson; i++) {
    if (i%1000000==0) cout<<" "<<i<<"/"<<N*frac_bulk_e_poisson<<endl;
    // Reset variables
    Reset();
    // Generate leakage charge position
    start_position = dis_uni(gen);
    int num_charges = dis_background(gen);
    // Get energy from charge transport
    for (int c=0; c<num_charges; c++) {
      energy += EnergyFromCharge(start_position, -1);
    }
    // Apply resolution
    energy += dis_res_0(gen);
    //    energy *= dis_res(gen);
    // Apply threshhold
    if (energy<dis_thresh(gen)) continue;
    // Fill
    outtree->Fill();
  }
  // laser
  event_type = 0;
  binomial_distribution<int> dis_laser(bin_n_laser,bin_p_laser);
  for (int i=0; i<N*frac_laser; i++) {
    if (i%1000000==0) cout<<" "<<i<<"/"<<N*frac_laser<<endl;
    // Reset variables
    Reset();
    // Generate leakage charge position
    start_position = 0;
    int num_charges = dis_laser(gen);
    // Get energy from charge transport
    for (int c=0; c<num_charges; c++) {
      energy += EnergyFromCharge(start_position, -1);
    }
    // Apply resolution
    energy += dis_res_0(gen);
    // Apply threshhold
    if (energy<dis_thresh(gen)) continue;
    // Fill
    outtree->Fill();
  }
  outfile->cd();
  cout<<"Writing "<<outfilename.Data()<<endl;
  outtree->Write();
  outfile->Close();

  return;
}

int main(int argc, char* argv[]) {
  if (argc<2) {
    cout<<"Invalid command line arguments. Aborting."<<endl;
    return 1;
  }
  int N = atoi(argv[1]); // number of events 
  string cfg = argv[2]; // config filename 
  Spectrum* s = new Spectrum(cfg);
  s->Run(N);
  return 0;
}











// int main(int argc, char* argv[]) {
//   if (argc<4) {
//     cout<<"Invalid command line arguments. Aborting."<<endl;
//     return 1;
//   }
//   string cfg = argv[2]; // config filename 
//   string suffix = argv[3]; // output filename suffix

//   double lambda_impact_ee = stod(getVariable(cfg,"lambda_impact_ee"));
//   double lambda_impact_eh = stod(getVariable(cfg,"lambda_impact_eh"));
//   double lambda_impact_he = stod(getVariable(cfg,"lambda_impact_he"));
//   double lambda_impact_hh = stod(getVariable(cfg,"lambda_impact_hh"));
//   double lambda_trap = stod(getVariable(cfg,"lambda_trap"));
//   double lambda_bulk_e_possion = stod(getVariable(cfg,"lambda_bulk_e_poisson"));
//   double frac_surf_e = stod(getVariable(cfg,"frac_surf_e"));
//   double frac_bulk_e = stod(getVariable(cfg,"frac_bulk_e"));
//   double frac_bulk_e_poisson = stod(getVariable(cfg,"frac_bulk_e_poisson"));
//   double frac_bulk_h = stod(getVariable(cfg,"frac_bulk_h"));
//   double resolution = stod(getVariable(cfg,"resolution"));
//   double thresh = stod(getVariable(cfg,"thresh"));
//   double thresh_res = stod(getVariable(cfg,"thresh_res"));
//   int max_impacts = stoi(getVariable(cfg,"max_impacts"));

//   // double lambda_impact = getVariable(cfg,"lambdaimpact");
//   // double lambda_trap = getVariable(cfg,"lambdatrap");
//   // double surf_frac = getVariable(cfg,"surffrac");
//   // double res = getVariable(cfg,"res");
//   // double thresh = getVariable(cfg,"thresh");
//   // double thresh_res = getVariable(cfg,"threshres");
//   // int max_impacts = int(getVariable(cfg,"maxcharges"));
//   // double background_frac = getVariable(cfg,"backgroundfrac");
//   // double lambda_background = getVariable(cfg,"lambdabackground");

//   TRandom3* myRNG = new TRandom3();
//   gRandom = myRNG;
//   std::normal_distribution<> dis_res(1,res);
//   std::normal_distribution<> dis_thresh(thresh,thresh_res);

//   int N_steps = 100;
//   double hist_max = 4.5;
//   TString outfilename = Form("simspectrum-%s.root",suffix.c_str());
//   TFile* outfile = new TFile(outfilename,"RECREATE");
//   TTree* outtree = new TTree("events","events");
//   double energy = 0;
//   double start_position = 0;
//   int num_impacts = 0;
//   int num_traps = 0;
//   int event_type = 0; // 0 = leakage, 1 = possion background
//   outtree->Branch("energy",&energy);
//   outtree->Branch("num_impacts",&num_impacts);
//   outtree->Branch("num_traps",&num_traps);
//   outtree->Branch("lambda_impact",&lambda_impact);
//   outtree->Branch("lambda_trap",&lambda_trap);
//   outtree->Branch("start_position",&start_position);
//   outtree->Branch("event_type",&event_type);
//   cout<<"Prob of Impact Ionization: "<<lambda_impact<<endl;
//   cout<<"Fraction of Surface Events: "<<surf_frac<<endl;
//   // Event loop
//   event_type = 0;
//   for (int i=0; i<N; i++) {
//     if (i%1000000==0) cout<<" "<<i<<"/"<<N<<endl;
//     // Generate leakage charge position
//     bool second_surf = false; // did a surface charge generate a second surface charge
//     if (i<N*surf_frac)
//       start_position = 0.;
//     else
//       start_position = dis_uni(gen);
//     // Reset variables
//     energy = 0;
//     num_impacts = 0;
//     // Get energy from charge transport
//     energy = energyFromElectron(start_position, lambda_impact, lambda_trap, num_impacts, num_traps, max_impacts);
//     // Apply resolution
//     energy *= dis_res(gen);
//     // Apply threshhold
//     if (energy<dis_thresh(gen)) continue;
//     // Fill
//     outtree->Fill();
//   }
//   // Poisson backgound event loop
//   event_type = 1;
//   poisson_distribution<int> dis_background(2);
//   for (int i=0; i<N*background_frac; i++) {
//     if (i%1000000==0) cout<<" "<<i<<"/"<<N*background_frac<<endl;
//     // Generate leakage charge position
//     start_position = dis_uni(gen);
//     // Reset variables
//     energy = 0;
//     num_impacts = 0;
//     int num_charges = dis_background(gen);
//     // Get energy from charge transport
//     for (int c=0; c<num_charges; c++) {
//       energy += energyFromCharge(start_position, lambda_impact, lambda_trap, num_impacts, num_traps, max_impacts);
//     }
//     // Apply resolution
//     energy *= dis_res(gen);
//     // Apply threshhold
//     if (energy<dis_thresh(gen)) continue;
//     // Fill
//     outtree->Fill();
//   }
  
//   outfile->cd();
//   cout<<"Writing "<<outfilename.Data()<<endl;
//   outtree->Write();
//   outfile->Close();
//   return 1;
// }
