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

//  gRandom->GetSeed();

/* The LoadConfig function loads the configuration file given by filename
   It returns a map of key-value pairs stored in the conifuration file */
// std::map<std::string,std::string> LoadConfig(std::string filename)
// {
//   std::map<std::string,std::string> ans; //A map of key-value pairs in the file
//   std::ifstream input(filename); //The input stream
//   while(input) //Keep on going as long as the file stream is good
//     {
//       std::string key; //The key
//       std::string value; //The value
//       std::getline(input, key, ':'); //Read up to the : delimiter into key
//       std::getline(input, value, '\n'); //Read up to the newline into value
//       std::string::size_type pos1 = value.find_first_of("\""); //Find the first quote in the value
//       std::string::size_type pos2 = value.find_last_of("\""); //Find the last quote in the value
//       if(pos1 != std::string::npos && pos2 != std::string::npos && pos2 > pos1) //Check if the found positions are all valid
//         {
// 	  //value = value.substr(pos1+1,pos2-pos1-1); //Take a substring of the part between the quotes
// 	  //	  ans[key] = value; //Store the result in the map
// 	  //ans.insert(pair<string,string>("key",0));// = value; //Store the result in the map
//         }
//     }
//   input.close(); //Close the file stream
//   return ans; //And return the result
// }

// int setVariable(string varFile, string varName, string varVal) {
//   ifstream filein(varFile); //File to read from
//   ofstream fileout("_fileout.sty"); //Temporary file
//   if(!filein || !fileout) {
//     cout << "Error opening files!" << endl;
//     return 1;
//   }
//   string strTemp;
//   bool found = false;
//   while(filein >> strTemp) {
//     string varNameTemp = strTemp.substr(0, strTemp.find('=')-1); 
//     if(varNameTemp == varName) {
//       fileout << varName << " = " << varVal << "\n";
//       found = true;
//     } else {
//       fileout << strTemp << "\n";
//     }
//     //if(found) break; // Update only first occurance
//    }
//   if (!found) {
//     fileout << varName << "=" << varVal << "\n";
//   }
//   rename("_fileout.sty",varFile.c_str());
//   return 0;
// }

// string getVariable(string varFile, string varName) {
//   ifstream filein(varFile); //File to read from
//   if(!filein) {
//     cout << "Error opening file!" << endl;
//     return 0;
//   }
//   string strTemp;
//   bool readNextLine = false;
//   while(filein >> strTemp) {
//     string varNameTemp = strTemp.substr(0, strTemp.find('='));     
//     cout<<varNameTemp<<endl;
//     cout<<strTemp.substr(strTemp.find('=')+1,strTemp.size())<<endl;
//     if(strTemp == varName) return strTemp.substr(strTemp.find('=')+1,strTemp.size());
//   }
//   return "";
// }


int setVariable(TString varFile, TString varName, TString varVal) {
  varName = "\\newcommand{\\"+varName+"}";
  varVal = "{"+varVal+"}";
  ifstream filein(varFile.Data()); //File to read from
  ofstream fileout("fileout.sty"); //Temporary file
  if(!filein || !fileout) {
    cout << "Error opening files!" << endl;
    return 1;
  }
  string strTemp;
  bool found = false;
  bool skipNextLine = false;
  while(filein >> strTemp) {
    if (skipNextLine) {
      skipNextLine = false;
      continue;
    }
    if(strTemp == varName) {
      fileout << varName << "\n";
      fileout << varVal << "\n";
      skipNextLine = true;
      found = true;
    } else {
      fileout << strTemp << "\n";
    }
    //if(found) break; // Update only first occurance
   }
  if (!found) {
    fileout << varName << "\n";
    fileout << varVal << "\n";
  }
  rename("fileout.sty",varFile.Data());
  return 0;
}

double getVariable(TString varFile, TString varName) {
  varName = "\\newcommand{\\"+varName+"}";
  ifstream filein(varFile.Data()); //File to read from
  if(!filein) {
    cout << "Error opening file!" << endl;
    return 0;
  }
  string strTemp;
  bool readNextLine = false;
  while(filein >> strTemp) {
    if (readNextLine) {
      return stod(strTemp.substr(1,strTemp.length()-2));
      continue;
    }
    if(strTemp == varName) {
      readNextLine = true;
    }
  }
  return 0;
}

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


double energyFromCharge(double x, double lambda_impact, double lambda_trap, int &num_impacts, int &num_traps, int max_impacts) {
  if (num_impacts>max_impacts) return 0.; // Don't bother anymore
  if (x>1) return 0; // sanity check

  // MFP method
  double next_impact = x-log(dis_uni(gen))/lambda_impact;
  double next_trap = x-log(dis_uni(gen))/lambda_trap;
  if (next_impact<next_trap) { // charge impact ionizes before being trapped.
    if (next_impact>1) return 1-x; // charge makes it to 1 without colliding
    else { // impact occured
      num_impacts++;
      return
	next_impact-x // energy from this charge until its impact
	+ energyFromCharge(next_impact,lambda_impact,lambda_trap,num_impacts,num_traps,max_impacts) // energy from this charge after its impact
	+ energyFromCharge(next_impact,lambda_impact,lambda_trap,num_impacts,num_traps,max_impacts);// energy from impact charge
    }
  } else { // charge is trapped before impact ionizing
    if (next_trap>1) return 1-x; // charge makes it to 1 without trapping
    else {
      num_traps++;
      return next_trap-x; // energy from this charge until its trap
    }
  }
  
  // // Step method (impact only)
  // double E = 1-x; // initial energy  
  // const int N_steps = 100+100*lambda_impact;
  // for (double xi=x; xi<1; xi+=1./N_steps) {
  //   // See if you knock out a charge
  //   if (dis_uni(gen)<(lambda_impact/N_steps)) {
  //     num_impacts +=1;
  //     E += energyFromCharge(xi, lambda_impact, lambda_trap, num_impacts, num_traps, max_impacts);
  //   }
  //   // Generate only 1 additional charge
  //   // if (dis_uni(gen)<(lambda_impact/N_steps)) {
  //   //   impacts +=1;
  //   //   return E + 1-xi;
  //   // }
  // }
  // return E;
}

int main(int argc, char* argv[]) {
  if (argc<4) {
    cout<<"Invalid command line arguments. Aborting."<<endl;
    return 1;
  }
  const int N = atoi(argv[1]); // number of events
  string cfg = argv[2]; // config filename 
  string suffix = argv[3]; // config filename 

  // setVariable(cfg,"lambdaimpact","0.03");
  // setVariable(cfg,"surffrac","0.6");
  // setVariable(cfg,"res","0.1");
  // setVariable(cfg,"thresh","0.6");
  // setVariable(cfg,"thresh_res","0.05");
  // setVariable(cfg,"max_charges","10");

  // double lambda_impact = stod(getVariable(cfg,"lambdaimpact"));
  // double lambda_trap = stod(getVariable(cfg,"lambdatrap"));
  // double surf_frac = stod(getVariable(cfg,"surffrac"));
  // double res = stod(getVariable(cfg,"res"));
  // double thresh = stod(getVariable(cfg,"thresh"));
  // double thresh_res = stod(getVariable(cfg,"threshres"));
  // int max_charges = stoi(getVariable(cfg,"maxcharges"));
  // cout<<lambda_impact<<" "<<lambda_trap<<" "<<surf_frac<<" "<<res<<" "<<thresh<<" "<<thresh_res<<" "<<max_charges<<endl;
  // return 0;
  // double background_frac = stod(getVariable(cfg,"backgroundfrac"));
  // double lambda_background = stod(getVariable(cfg,"lambdabackground"));

  double lambda_impact = getVariable(cfg,"lambdaimpact");
  double lambda_trap = getVariable(cfg,"lambdatrap");
  double surf_frac = getVariable(cfg,"surffrac");
  double res = getVariable(cfg,"res");
  double thresh = getVariable(cfg,"thresh");
  double thresh_res = getVariable(cfg,"threshres");
  int max_charges = int(getVariable(cfg,"maxcharges"));
  double background_frac = getVariable(cfg,"backgroundfrac");
  double lambda_background = getVariable(cfg,"lambdabackground");

  TRandom3* myRNG = new TRandom3();
  gRandom = myRNG;
  std::normal_distribution<> dis_res(1,res);
  std::normal_distribution<> dis_thresh(thresh,thresh_res);

  int N_steps = 100;
  double hist_max = 4.5;
  TString outfilename = "simspectrum-"+suffix+".root";
  TFile* outfile = new TFile(outfilename,"RECREATE");
  TTree* outtree = new TTree("events","events");
  double energy = 0;
  double start_position = 0;
  int num_impacts = 0;
  int num_traps = 0;
  int event_type = 0; // 0 = leakage, 1 = possion background
  outtree->Branch("energy",&energy);
  outtree->Branch("num_impacts",&num_impacts);
  outtree->Branch("num_traps",&num_traps);
  outtree->Branch("lambda_impact",&lambda_impact);
  outtree->Branch("lambda_trap",&lambda_trap);
  outtree->Branch("start_position",&start_position);
  outtree->Branch("event_type",&event_type);
  cout<<"Prob of Impact Ionization: "<<lambda_impact<<endl;
  cout<<"Fraction of Surface Events: "<<surf_frac<<endl;
  // Event loop
  event_type = 0;
  for (int i=0; i<N; i++) {
    if (i%1000000==0) cout<<" "<<i<<"/"<<N<<endl;
    // Generate leakage charge position
    bool second_surf = false; // did a surface charge generate a second surface charge
    if (i<N*surf_frac)
      start_position = 0.;
    else
      start_position = dis_uni(gen);
    // Reset variables
    energy = 0;
    num_impacts = 0;
    // Get energy from charge transport
    energy = energyFromCharge(start_position, lambda_impact, lambda_trap, num_impacts, num_traps, max_charges);
    // Apply resolution
    energy *= dis_res(gen);
    // Apply threshhold
    if (energy<dis_thresh(gen)) continue;
    // Fill
    outtree->Fill();
  }
  // Poisson backgound event loop
  event_type = 1;
  poisson_distribution<int> dis_background(2);
  for (int i=0; i<N*background_frac; i++) {
    if (i%1000000==0) cout<<" "<<i<<"/"<<N*background_frac<<endl;
    // Generate leakage charge position
    start_position = dis_uni(gen);
    // Reset variables
    energy = 0;
    num_impacts = 0;
    int num_charges = dis_background(gen);
    // Get energy from charge transport
    for (int c=0; c<num_charges; c++) {
      energy += energyFromCharge(start_position, lambda_impact, lambda_trap, num_impacts, num_traps, max_charges);
    }
    // Apply resolution
    energy *= dis_res(gen);
    // Apply threshhold
    if (energy<dis_thresh(gen)) continue;
    // Fill
    outtree->Fill();
  }
  
  outfile->cd();
  cout<<"Writing "<<outfilename.Data()<<endl;
  outtree->Write();
  outfile->Close();
  return 1;
}
