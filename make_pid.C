#include "functions.h"
#include "style.h"
using namespace style;

void make_pid(TString cutg_name = "")
//void make_pid(TString cutg_name = "cut_proton")
{
  zcolor(1);

  //auto tag = "develop.1534.2537e2e";
  //TString tag = "yoffsettest.1543.d50e777";
  //TString tag = "yoffset_cal2.develop.RC.1612.62c6e62";
  //TString tag = "develop.RC.1612.62c6e62.yoffset_cal4";
  TString tag = "develop.RC.1612.62c6e62.cov_test";
  version(tag);

  auto tree = new TChain("tracks");
  for (auto run : {2894,2900,2901,2902,2903,2904,2905})
    tree -> Add(Form("/Users/ejungwoo/data/spirit/summary_track_run%d.%s.root",run,tag.Data()));

  if (!cutg_name.IsNull())
    cutg(Form("data/%s.%s.root",cutg_name.Data(),tag.Data()),cutg_name,"p","dedx");
  else {
    cout << endl;
    cout << "======================================================================" << endl;
    cout << "*** Instructions to make graphical cut" << endl;
    cout << endl;
    cout << " 1. Click [View] -> [Toolbar]" << endl;
    cout << " 2. Click [Graphical Cut]. It is the last icon in the [Toolbar]" << endl;
    cout << " 3. Click to make graphical cut points from the histogram" << endl;
    cout << "    (Double click to finish and link last to first point)" << endl;
    cout << " 4. Right Click the graphical cut line, and press [SetName]" << endl;
    cout << " 5. Set the [cutg_name](ex. \"cut_proton\") from the pop up window" << endl;
    cout << " 6. Type; write([cutg_name]) on the terminal. (ex. write(cut_proton))"<< endl;
    cout << " 6. root 'make_pid.C(\"[cutg_name]\") # to check the cut"<< endl;
    cout << "======================================================================" << endl;
    cout << endl;
  }

  TString name = "pid_";
  if (!cutg_name.IsNull()) name = name + cutg_name;
  else                     name = name + "raw";

  auto hist = tp(name,tree,"dedx:p",TCut(cutg_name),tag+";p/Z (MeV/c);dE/dx (ADC/mm)",400,-500,2000,400,0,500);
  auto cvs = cc2(name);
  hist -> GetZaxis() -> SetRangeUser(0,1000);
  hist -> Draw("colz");
  cvs -> SetLogz();

  save(cvs);
  write(hist);
}
