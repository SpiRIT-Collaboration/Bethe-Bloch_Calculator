#include "style.h"
#include "functions.h"
using namespace style;

void calculate_mass(TString fit_name = "fit_proton")
{
  //auto tag = "develop.1534.2537e2e";
  //TString tag = "yoffsettest.1543.d50e777";
  //TString tag = "yoffset_cal2.develop.RC.1612.62c6e62";
  //TString tag = "develop.RC.1612.62c6e62.yoffset_cal4";
  TString tag = "develop.RC.1612.62c6e62.cov_test";
  version(tag);

  auto tree = new TChain("tracks");
  for (auto run : {2894,2900,2901,2902,2903,2904,2905})
    tree -> Add(Form("/Users/ejungwoo/data/spirit/summary_track_run%d.%s.root",run,tag.Data()));

  Double_t p, dedx;
  tree -> SetBranchAddress("p",&p);
  tree -> SetBranchAddress("dedx",&dedx);

  auto minimizing_mass = new TF1("minimizing_mass", fddedx, -500., 5000., 5, 1);

  auto fitFile = new TFile(Form("data/%s.develop.RC.1612.62c6e62.cov_test.root",fit_name.Data()));
  auto fit = (TF1 *) fitFile -> Get(fit_name);
  minimizing_mass -> SetParameter(0, fit -> GetParameter(0));
  minimizing_mass -> SetParameter(1, fit -> GetParameter(1));

  auto massFile = new TFile(Form("data/mass_%s.%s.root",fit_name.Data(),tag.Data()),"recreate");
  auto massTree = new TTree("mass","");
  Double_t mass;
  massTree -> Branch("mass",&mass);

  auto numEntries = tree -> GetEntries();
  for (auto entry = 0; entry < numEntries; ++entry)
  {
    if (entry%50000 == 0)
      cout << entry << " / " << numEntries << endl;

    tree -> GetEntry(entry);

    auto charge = 1;
    if (p < 0) charge = -1;

    minimizing_mass -> SetParameter(2, p);
    minimizing_mass -> SetParameter(3, charge);
    minimizing_mass -> SetParameter(4, dedx);

    ROOT::Math::RootFinder finderp(ROOT::Math::RootFinder::kBRENT);
    ROOT::Math::WrappedTF1 funcp(*minimizing_mass);
    finderp.SetFunction(funcp, 0.1, 10000.);
    finderp.Solve();
    mass = finderp.Root();
    if (mass == 0)
      continue;

    massTree -> Fill();
  }

  massFile -> cd();
  massTree -> Write();
  massFile -> Print();

  TString name = Form("mass_%s.%s",fit_name.Data(),tag.Data());
  cc2(name);
  auto hist = new TH1D(name,"M",400,0,10000);
  massTree -> Draw(Form("mass>>%s",name.Data()));
  hist -> Write();
}
