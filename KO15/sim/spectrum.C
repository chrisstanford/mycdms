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
  double next_x = x+sign*pol*deltas[min_index];
  // Compute next step. The smallest delta determines the process that will occur.
  if (sign<0 && next_x>1) { // e makes it to 1 without colliding
    int num_surf_eh = dis_surf_eh1(gen);
    num_impacts_eh += num_surf_eh;
    double E = 1-x; // this charge
    for (int i=0; i<num_surf_eh; i++) E += EnergyFromCharge(1,+1); // hs freed at surface by e 
    return E;
  } else if (sign<0 && next_x<0) { // e makes it to 0 without colliding
    int num_surf_eh = dis_surf_eh0(gen);
    num_impacts_eh += num_surf_eh;
    double E = x; // this charge
    for (int i=0; i<num_surf_eh; i++) E += EnergyFromCharge(0,+1); // hs freed at surface by e 
    return E;
  } else if (sign>0 && next_x<0) { // h makes it to 0 without colliding
    int num_surf_he = dis_surf_he0(gen);
    num_impacts_he += num_surf_he;
    double E = x; // this charge
    for (int i=0; i<num_surf_he; i++) E += EnergyFromCharge(0,-1); // es freed at surface by h
    return E;
  } else if (sign>0 && next_x>1) { // h makes it to 1 without colliding
    int num_surf_he = dis_surf_he1(gen);
    num_impacts_he += num_surf_he;
    double E = 1-x; // this charge
    for (int i=0; i<num_surf_he; i++) E += EnergyFromCharge(1,-1); // es freed at surface by h
    return E;
  } else if (min_index==0 && sign<0) { // e impact ionizes and produces e
    num_impacts_ee++;
    return
      delta_impact_ee // energy from this charge until its impact
      + EnergyFromCharge(next_x,-1) // energy from this charge after its impact
      + EnergyFromCharge(next_x,-1);// energy from impact charge
  } else if (min_index==1 && sign<0) {// e impact ionizes and produces h
    num_impacts_eh++;
    return
      delta_impact_eh // energy from this charge until its impact
      + EnergyFromCharge(next_x,-1) // energy from this charge after its impact
      + EnergyFromCharge(next_x,+1);// energy from impact charge
  } else if (min_index==2 && sign>0) {// h impact ionizes and produces e
    num_impacts_he++;
    return
      delta_impact_he // energy from this charge until its impact
      + EnergyFromCharge(next_x,+1) // energy from this charge after its impact
      + EnergyFromCharge(next_x,-1);// energy from impact charge
  } else if (min_index==3 && sign>0) {// h impact ionizes and produces h
    num_impacts_hh++;
    return
      delta_impact_hh // energy from this charge until its impact
      + EnergyFromCharge(next_x,+1) // energy from this charge after its impact
      + EnergyFromCharge(next_x,+1);// energy from impact charge
  } else if (min_index==4 && sign<0) {// e gets trapped
    num_traps_e++;
    return delta_trap_e; // energy from this charge until its trap
  } else if (min_index==5 && sign>0) {// h gets trapped
    num_traps_h++;
    return delta_trap_h; // energy from this charge until its trap
  } else {// this should never actually happen
    cout<<"Reached end of if block. Check code"<<endl;
    return 0;
  }
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
  laser=false;
  return;
}

Spectrum::Spectrum(string cfg) {
  LoadConfig(cfg);

  gen.seed(rd());
  dis_uni.param(uniform_real_distribution<>::param_type(0,1));
  dis_surf_e.param(poisson_distribution<int>::param_type(lambda_surf_e));
  dis_surf_h.param(poisson_distribution<int>::param_type(lambda_surf_h)); 
  dis_surf_eh0.param(poisson_distribution<int>::param_type(lambda_surf_eh0));
  dis_surf_he0.param(poisson_distribution<int>::param_type(lambda_surf_he0)); 
  dis_surf_eh1.param(poisson_distribution<int>::param_type(lambda_surf_eh1));
  dis_surf_he1.param(poisson_distribution<int>::param_type(lambda_surf_he1));
  dis_laser_pos.param(poisson_distribution<int>::param_type(lambda_laser_pos));
  dis_laser_neg.param(poisson_distribution<int>::param_type(lambda_laser_neg));
  
}

void Spectrum::LoadConfig(string cfg) {
  resolution = stod(getVariable(cfg,"resolution"));
  thresh = stod(getVariable(cfg,"thresh"));
  thresh_res = stod(getVariable(cfg,"thresh_res"));
  max_impacts = stoi(getVariable(cfg,"max_impacts"));
  lambda_trap_e = stod(getVariable(cfg,"lambda_trap_e"));
  lambda_trap_h = stod(getVariable(cfg,"lambda_trap_h"));
  lambda_impact_ee = stod(getVariable(cfg,"lambda_impact_ee"));
  lambda_impact_eh = stod(getVariable(cfg,"lambda_impact_eh"));
  lambda_impact_he = stod(getVariable(cfg,"lambda_impact_he"));
  lambda_impact_hh = stod(getVariable(cfg,"lambda_impact_hh"));
  lambda_surf_e = stod(getVariable(cfg,"lambda_surf_e"));
  lambda_surf_h = stod(getVariable(cfg,"lambda_surf_h"));
  lambda_surf_eh0 = stod(getVariable(cfg,"lambda_surf_eh0"));
  lambda_surf_he0 = stod(getVariable(cfg,"lambda_surf_he0"));
  lambda_surf_eh1 = stod(getVariable(cfg,"lambda_surf_eh1"));
  lambda_surf_he1 = stod(getVariable(cfg,"lambda_surf_he1"));
  lambda_laser_pos = stod(getVariable(cfg,"lambda_laser_pos"));
  lambda_laser_neg = stod(getVariable(cfg,"lambda_laser_neg"));
  rate_laser_pos = stod(getVariable(cfg,"rate_laser_pos"));
  rate_laser_neg = stod(getVariable(cfg,"rate_laser_neg"));
  rate_surf_e_pos = stod(getVariable(cfg,"rate_surf_e_pos"));
  rate_surf_e_neg = stod(getVariable(cfg,"rate_surf_e_neg"));
  rate_surf_h_pos = stod(getVariable(cfg,"rate_surf_h_pos"));
  rate_surf_h_neg = stod(getVariable(cfg,"rate_surf_h_neg"));
  rate_bulk_e_pos = stod(getVariable(cfg,"rate_bulk_e_pos"));
  rate_bulk_e_neg = stod(getVariable(cfg,"rate_bulk_e_neg"));
  rate_bulk_h_pos = stod(getVariable(cfg,"rate_bulk_h_pos"));
  rate_bulk_h_neg = stod(getVariable(cfg,"rate_bulk_h_neg"));
  livedays_pos = stod(getVariable(cfg,"livedays_pos"));
  livedays_neg = stod(getVariable(cfg,"livedays_neg"));
  suffix = getVariable(cfg,"suffix");
  

}

void Spectrum::Run() {
  TRandom3* myRNG = new TRandom3();
  gRandom = myRNG;
  std::normal_distribution<> dis_res(1,resolution);
  std::normal_distribution<> dis_res_0(0,resolution);
  std::normal_distribution<> dis_res_1(0,resolution/2);
  std::normal_distribution<> dis_thresh(thresh,thresh_res);

  TString outfilename = Form("simspectrum-%s.root",suffix.c_str());
  TFile* outfile = new TFile(outfilename,"RECREATE");
  TTree* outtree = new TTree("events","events");
  int event_type = 0; // 0 = leakage, 1 = possion background
  double t = 0;
  outtree->Branch("energy",&energy);
  outtree->Branch("num_impacts_ee",&num_impacts_ee);
  outtree->Branch("num_impacts_eh",&num_impacts_eh);
  outtree->Branch("num_impacts_he",&num_impacts_he);
  outtree->Branch("num_impacts_hh",&num_impacts_hh);
  outtree->Branch("num_traps_e",&num_traps_e);
  outtree->Branch("num_traps_h",&num_traps_h);
  outtree->Branch("start_position",&start_position);
  outtree->Branch("event_type",&event_type);
  outtree->Branch("pol",&pol);
  outtree->Branch("laser",&laser);
  outtree->Branch("t",&t);

  // Loop over polarity
  for (int p=-1; p<=1; p+=2) {
    pol = p;

    double rate_laser  = (pol>0) ? rate_laser_pos : rate_laser_neg;
    double rate_surf_e = (pol>0) ? rate_surf_e_pos : rate_surf_e_neg;
    double rate_surf_h = (pol>0) ? rate_surf_h_pos : rate_surf_h_neg;
    double rate_bulk_e = (pol>0) ? rate_bulk_e_pos : rate_bulk_e_neg;
    double rate_bulk_h = (pol>0) ? rate_bulk_h_pos : rate_bulk_h_neg;
    double livedays    = (pol>0) ? livedays_pos : livedays_neg;

    double next_laser  = (rate_laser >0) ? -1 : numeric_limits<double>::max();
    double next_surf_e = (rate_surf_e>0) ? -1 : numeric_limits<double>::max();
    double next_surf_h = (rate_surf_h>0) ? -1 : numeric_limits<double>::max();
    double next_bulk_e = (rate_bulk_e>0) ? -1 : numeric_limits<double>::max();
    double next_bulk_h = (rate_bulk_h>0) ? -1 : numeric_limits<double>::max();

    int num_laser = livedays*24*60*60*rate_laser; // estimated number of laser events;
    t=0; // start time
    while (t<livedays*24*60*60) {
      int i_laser = t*rate_laser;
      //      cout<<i_laser<<" "<<num_laser<<endl;
      if (num_laser>0 && (i_laser*10%num_laser)==0) cout<<i_laser*100/num_laser<<"%"<<endl;
      
      //      if (fmod(t,24*60*60)<1) cout<<int(t/(24*60*60))<<"/"<<livedays<<endl;
      // Determine next event time and type
      if (rate_laser>0  && t>=next_laser)  next_laser  = t+1./rate_laser;
      if (rate_surf_e>0 && t>=next_surf_e) next_surf_e = t-log(dis_uni(gen))/rate_surf_e;
      if (rate_surf_h>0 && t>=next_surf_h) next_surf_h = t-log(dis_uni(gen))/rate_surf_h;
      if (rate_bulk_e>0 && t>=next_bulk_e) next_bulk_e = t-log(dis_uni(gen))/rate_bulk_e;
      if (rate_bulk_h>0 && t>=next_bulk_h) next_bulk_h = t-log(dis_uni(gen))/rate_bulk_h;
      const int num_event_types = 5;
      double next_times[num_event_types] = {next_laser,next_surf_e,next_surf_h,next_bulk_e,next_bulk_h};
      int min_index = 0;
      for(int i=1; i<num_event_types; i++) {
	if(next_times[i] < next_times[min_index])
	  min_index = i;
      }
      t = next_times[min_index];
      // Reset variables
      Reset();
      laser=false;
      // Generate charge position and calculate energy
      int sign = 0; // e or h?
      int num_charges = 1; // how many charges? default = 1
      if (min_index==0) { // laser
	laser = true;
	start_position = 0.;
	sign = pol;
	num_charges = (pol>0) ? dis_laser_pos(gen) : dis_laser_neg(gen);
      } else if (min_index==1) { // surf e
	start_position = double(pol+1)/2;
	sign = -1;
	num_charges += dis_surf_e(gen);
      } else if (min_index==2) { // surf h
	start_position = double(-pol+1)/2;
	sign = +1;
	num_charges += dis_surf_h(gen);
      } else if (min_index==3) { // bulk e
	start_position = dis_uni(gen);
	sign = -1;
      } else if (min_index==4) { // bulk h
	start_position = dis_uni(gen);
	sign = +1;
      } else continue;
      for (int c=0; c<num_charges; c++) energy += EnergyFromCharge(start_position, sign);
      // Apply resolution
      energy += dis_res_0(gen);
      // test
      // if (energy>0.8) energy += dis_res_0(gen);
      // else energy += dis_res_1(gen);
      // Apply threshhold
      if (energy<dis_thresh(gen)) continue;
      // Fill
      outtree->Fill();
    } // End event loop
  } // End loop over polarity
    
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
  string cfg = argv[1]; // config filename 
  Spectrum* s = new Spectrum(cfg);
  s->Run();
  return 0;
}
