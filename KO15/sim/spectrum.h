#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <string>

using namespace std;

class Spectrum {
 public:
  // input
  double lambda_impact_ee;
  double lambda_impact_eh;
  double lambda_impact_he;
  double lambda_impact_hh;
  double lambda_trap_e;
  double lambda_trap_h;
  double lambda_bulk_e_possion;
  double lambda_surf_e_possion;
  double frac_surf_e;
  double frac_surf_h;
  double frac_bulk_e;
  double frac_bulk_e_poisson;
  double frac_bulk_h;
  double resolution;
  double thresh;
  double thresh_res;
  string suffix;
  int max_impacts;

  double frac_laser;
  int bin_n_laser;
  double bin_p_laser;

  
  // local variables
  double energy;
  double start_position;
  int num_impacts_ee;
  int num_impacts_eh;
  int num_impacts_he;
  int num_impacts_hh;
  int num_traps_e;
  int num_traps_h;

  // functions
  Spectrum(string cfg);
  void Reset();
  void LoadConfig(string cfg);
  void Run(int N);
  double EnergyFromCharge(double x, int sign);
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
