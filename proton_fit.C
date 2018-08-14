#include "functions.h"
#include "style.h"
using namespace style;

void proton_fit(TString cutg_name = "cut_proton")
{
  zcolor(1);
  
  bool draw_data_point = false;
  bool draw_raw_data = true;

  //auto tag = "develop.1534.2537e2e";
  //TString tag = "yoffsettest.1543.d50e777";
  //TString tag = "yoffset_cal2.develop.RC.1612.62c6e62";
  //TString tag = "develop.RC.1612.62c6e62.yoffset_cal4";
  TString tag = "develop.RC.1612.62c6e62.cov_test";
  version(tag);

  auto xbin1 = 145;
  auto xbin2 = 260;
  auto dxbin = 5;

  auto ybin1 = 100;
  auto ybin2 = 300;
  auto dybin = 10;



  auto file = new TFile(Form("data/pid_%s.%s.root",cutg_name.Data(),tag.Data()));
  auto hist = (TH2D *) file -> Get(Form("pid_%s",cutg_name.Data()));
  auto graph = new TGraphErrors();

  for (auto xbin=xbin1; xbin<xbin2; xbin+=dxbin) {
    auto bin1 = xbin;
    auto bin2 = xbin + dxbin - 1;
    auto x1 = hist -> GetXaxis() -> GetBinLowEdge(bin1);
    auto x2 = hist -> GetXaxis() -> GetBinUpEdge(bin2);
    auto hp = hist -> ProjectionY(TString(hist->GetName())+Form("_%d(%.2f)_%d(%.2f)",bin1,x1,bin2,x2),bin1,bin2);
    auto fp = fitg(hp);
    auto mean = fp -> GetParameter(1);
    auto sigma = fp -> GetParameter(2);
    graph -> SetPoint(graph->GetN(), (x1+x2)/2., mean);
    graph -> SetPointError(graph->GetN()-1, (x2-x1)/2., sigma);

    if (draw_data_point) {
      cc();
      hp -> Draw();
      fp -> Draw("same");
      hp -> GetXaxis() -> SetRangeUser(mean-5*sigma,mean+5*sigma);
    }
  }

  for (auto ybin=ybin1; ybin<ybin2; ybin+=dybin) {
    auto bin1 = ybin;
    auto bin2 = ybin + dybin - 1;
    auto y1 = hist -> GetYaxis() -> GetBinLowEdge(bin1);
    auto y2 = hist -> GetYaxis() -> GetBinUpEdge(bin2);
    auto hp = hist -> ProjectionX(TString(hist->GetName())+Form("_%d(%.2f)_%d(%.2f)",bin1,y1,bin2,y2),bin1,bin2);
    auto fp = fitg(hp);
    auto mean = fp -> GetParameter(1);
    auto sigma = fp -> GetParameter(2);
    graph -> SetPoint(graph->GetN(), mean, (y1+y2)/2.);
    graph -> SetPointError(graph->GetN()-1, sigma, (y2-y1)/2.);

    if (draw_data_point) {
      cc();
      hp -> Draw();
      fp -> Draw("same");
      hp -> GetXaxis() -> SetRangeUser(mean-5*sigma,mean+5*sigma);
    }
  }

  auto fit_proton = new TF1("fit_proton", fdedx_p, 0., 5000., 2, 1);
  make (fit_proton);
  fit_proton -> SetParameters(-5.85199e-02, 6.21632e+02);
  graph -> Fit(fit_proton,"N");
  auto c1 = fit_proton -> GetParameter(0);
  auto c2 = fit_proton -> GetParameter(1);

  {
    TString name = TString(hist->GetName()) + "_fit_data";
    auto cvs = cc2(name);
    hist -> Draw("colz");
    graph -> Draw("samep");
    fit_proton -> Draw("same");
    cvs -> SetLogz();
    write(fit_proton);
  }

  if (draw_raw_data) {
    TString name = TString(hist->GetName()) + "_toraw";
    auto file_raw = new TFile(Form("data/pid_raw.%s.root",tag.Data()));
    auto hist_raw = (TH2D *) file_raw -> Get("pid_raw");
    auto cvs = cc2(name);
    hist_raw -> Draw("colz");
    cvs -> SetLogz();

    TF1 *fdedx = nullptr;
    fdedx = new TF1("fpim",fdedx_pim, -500., 0, 2, 1); fdedx -> SetParameters(c1,c2); fdedx -> Draw("same");
    fdedx = new TF1("fpi", fdedx_pi,  0., 5000., 2, 1); fdedx -> SetParameters(c1,c2); fdedx -> Draw("same");
    fdedx = new TF1("fp",  fdedx_p,   0., 5000., 2, 1); fdedx -> SetParameters(c1,c2); fdedx -> Draw("same");
    fdedx = new TF1("fd",  fdedx_d,   0., 5000., 2, 1); fdedx -> SetParameters(c1,c2); fdedx -> Draw("same");
    fdedx = new TF1("ft",  fdedx_t,   0., 5000., 2, 1); fdedx -> SetParameters(c1,c2); fdedx -> Draw("same");
  }
}
