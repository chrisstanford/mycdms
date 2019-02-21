#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <string>
#include <random>

using namespace std;

class Spectrum {
 public:


  // input
  double resolution;
  double thresh;
  double thresh_res;
  int max_impacts;
  double lambda_trap_e;
  double lambda_trap_h;
  double lambda_impact_ee;
  double lambda_impact_eh;
  double lambda_impact_he;
  double lambda_impact_hh;
  double lambda_surf_e;
  double lambda_surf_h;
  double lambda_surf_eh0;
  double lambda_surf_he0;
  double lambda_surf_eh1;
  double lambda_surf_he1;
  double lambda_laser_pos;
  double lambda_laser_neg;
  double rate_laser_pos;
  double rate_laser_neg;
  double rate_surf_e_pos;
  double rate_surf_e_neg;
  double rate_surf_h_pos;
  double rate_surf_h_neg;
  double rate_bulk_e_pos;
  double rate_bulk_e_neg;
  double rate_bulk_h_pos;
  double rate_bulk_h_neg;
  double livedays_pos;
  double livedays_neg;
  string suffix;

  
  // local variables
  double energy;
  double start_position;
  int num_impacts_ee;
  int num_impacts_eh;
  int num_impacts_he;
  int num_impacts_hh;
  int num_traps_e;
  int num_traps_h;
  int pol; // polarity
  bool laser;
  
  // functions
  Spectrum(string cfg);
  void Reset();
  void LoadConfig(string cfg);
  void Run();
  double EnergyFromCharge(double x, int sign);

  // random
  random_device rd;  //Will be used to obtain a seed for the random number engine
  mt19937 gen; //Standard mersenne_twister_engine seeded with rd()
  uniform_real_distribution<> dis_uni;//(0,1);
  poisson_distribution<int> dis_surf_e;
  poisson_distribution<int> dis_surf_h;
  poisson_distribution<int> dis_surf_eh0;
  poisson_distribution<int> dis_surf_he0;
  poisson_distribution<int> dis_surf_eh1;
  poisson_distribution<int> dis_surf_he1;
  poisson_distribution<int> dis_laser_pos;
  poisson_distribution<int> dis_laser_neg;

};


// int setVariable(TString varFile, TString varName, TString varVal) {
//   varName = "\\newcommand{\\"+varName+"}";
//   varVal = "{"+varVal+"}";
//   ifstream filein(varFile.Data()); //File to read from
//   ofstream fileout("fileout.sty"); //Temporary file
//   if(!filein || !fileout) {
//     cout << "Error opening files!" << endl;
//     return 1;
//   }
//   string strTemp;
//   bool found = false;
//   bool skipNextLine = false;
//   while(filein >> strTemp) {
//     if (skipNextLine) {
//       skipNextLine = false;
//       continue;
//     }
//     if(strTemp == varName) {
//       fileout << varName << "\n";
//       fileout << varVal << "\n";
//       skipNextLine = true;
//       found = true;
//     } else {
//       fileout << strTemp << "\n";
//     }
//     //if(found) break; // Update only first occurance
//    }
//   if (!found) {
//     fileout << varName << "\n";
//     fileout << varVal << "\n";
//   }
//   rename("fileout.sty",varFile.Data());
//   return 0;
// }

// double getVariable(TString varFile, TString varName) {
//   varName = "\\newcommand{\\"+varName+"}";
//   ifstream filein(varFile.Data()); //File to read from
//   if(!filein) {
//     cout << "Error opening file!" << endl;
//     return 0;
//   }
//   string strTemp;
//   bool readNextLine = false;
//   while(filein >> strTemp) {
//     if (readNextLine) {
//       return stod(strTemp.substr(1,strTemp.length()-2));
//       continue;
//     }
//     if(strTemp == varName) {
//       readNextLine = true;
//     }
//   }
//   return 0;
// }




#endif
