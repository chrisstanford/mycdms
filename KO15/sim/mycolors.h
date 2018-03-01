#ifndef MYCOLORS_H
#define MYCOLORS_H
#include "TColor.h"
#include "TString.h"

Int_t mycolors8[8] = {TColor::GetColor("#000000"),
			TColor::GetColor("#e69f00"),
			TColor::GetColor("#56b4e9"),
			TColor::GetColor("#009e73"),
			TColor::GetColor("#f0e442"),
			TColor::GetColor("#0072b2"),
			TColor::GetColor("#cc79a7"),
			TColor::GetColor("#d55e00")};

Int_t mycolors12[12] = {
  TColor::GetColor("#aa0a3c"), //Red
  TColor::GetColor("#fa7850"), //Orange
  TColor::GetColor("#0ab45a"), //Green
  TColor::GetColor("#00a0fa"), //Light Blue
  TColor::GetColor("#8214a0"), //Purple
  TColor::GetColor("#005ac8"), //Dark Blue
  TColor::GetColor("#fa78fa"), //Magenta
  TColor::GetColor("#006e82"), //Teal
  TColor::GetColor("#14d2dc"), //Cyan
  TColor::GetColor("#f0f032"), //Yellow
  TColor::GetColor("#a0fa82"), //Light Green
  TColor::GetColor("#fae6be")}; //Peach

/* struct myColors8 { */
/*   const int nColors; */
  
/* myColors8() : nColors(8) {}; */
/*   TString colors = {"#000000","#e69f00","#56b4e9","#009e73","#f0e442","#0072b2","#d55e00","#cc79a7"}; */
/*   TColor GetColor(int i) { */
    
/*   } */
/* }; */

#endif
