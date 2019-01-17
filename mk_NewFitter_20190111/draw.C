#include "MassEstimator.h"
using namespace MassEstimator;

TString particle[] = {"proton","deuteron","triton","3He","alpha","6He","6Li","7Li"};
TF1 *bbfit[8], *lvfit[8];
Double_t minp = 100.;
Double_t maxp = 3000.;

void drawRawFunction()
{
	gStyle->SetOptStat(0);

	auto cvsBB = new TCanvas("cvsBB","",1000,700);
	cvsBB->SetLogz();
	cvsBB->SetLogy();

	auto h1 = new TH1F("h1",";Magnetic rigidity (MeV/c/Z); dE/dx (MeV/cm)",10,minp,maxp);
	h1->Draw();
	h1->GetYaxis()->SetRangeUser(0.001,2.);

	bbfit[0] = new TF1("bbfit_p",fBBdedx_p,minp,maxp,2,1);
	bbfit[1] = new TF1("bbfit_d",fBBdedx_d,minp,maxp,2,1);
	bbfit[2] = new TF1("bbfit_t",fBBdedx_t,minp,maxp,2,1);
	bbfit[3] = new TF1("bbfit_he3",fBBdedx_he3,minp,maxp,2,1);
	bbfit[4] = new TF1("bbfit_al",fBBdedx_al,minp,maxp,2,1);
	bbfit[5] = new TF1("bbfit_he6",fBBdedx_he6,minp,maxp,2,1);
	bbfit[6] = new TF1("bbfit_li6",fBBdedx_li6,minp,maxp,2,1);
	bbfit[7] = new TF1("bbfit_li7",fBBdedx_li7,minp,maxp,2,1);
	lvfit[0] = new TF1("lvfit_p",fLVdedx_p,minp,maxp,3,1);
	lvfit[1] = new TF1("lvfit_d",fLVdedx_d,minp,maxp,3,1);
	lvfit[2] = new TF1("lvfit_t",fLVdedx_t,minp,maxp,3,1);
	lvfit[3] = new TF1("lvfit_he3",fLVdedx_he3,minp,maxp,3,1);
	lvfit[4] = new TF1("lvfit_al",fLVdedx_al,minp,maxp,3,1);
	lvfit[5] = new TF1("lvfit_he6",fLVdedx_he6,minp,maxp,3,1);
	lvfit[6] = new TF1("lvfit_li6",fLVdedx_li6,minp,maxp,3,1);
	lvfit[7] = new TF1("lvfit_li7",fLVdedx_li7,minp,maxp,3,1);

	auto latex = new TLatex();
	latex->SetTextFont(132);
	latex->SetTextSize(0.06);
	latex->DrawLatexNDC(0.22,0.82,"Solid line : Bethe-Bloch");
	latex->DrawLatexNDC(0.22,0.74,"Dotted line: Landau-Vavilov");
	latex->SetTextSize(0.05);
	for(Int_t i=0; i<5; i++){
		bbfit[i]->SetLineColor(i+2);
		bbfit[i]->SetParameters(1.,0.);
		bbfit[i]->SetLineWidth(2);
		bbfit[i]->SetNpx(2000.);
		lvfit[i]->SetLineColor(i+2);
		lvfit[i]->SetParameters(1.,0.,1.);
		lvfit[i]->SetLineWidth(2);
		lvfit[i]->SetLineStyle(2);
		lvfit[i]->SetNpx(2000.);
		latex->SetTextColor(i+2);
		latex->DrawLatexNDC(0.72,0.85-i*0.05,particle[i]);
		bbfit[i]->Draw("same");
		lvfit[i]->Draw("same");
	}
	//cvsBB->SaveAs("BBandLV.eps");

}

void draw()
{
	gStyle->SetOptStat(0);
	gStyle->SetPalette(kRainBow);
	auto f = TFile::Open("goodPID.root");
	TH2* h = nullptr; f->GetObject("h2PID",h);

	auto cvsBB = new TCanvas("cvsBB","",1400,1000);
	cvsBB->SetLogz();
	cvsBB->SetLogy();
	h->Draw("col");
	auto cvsLV = new TCanvas("cvsLV","",1400,1000);
	cvsLV->SetLogz();
	cvsLV->SetLogy();
	h->Draw("col");

	bbfit[0] = new TF1("bbfit_p",fBBdedx_p,minp,maxp,2,1);
	bbfit[1] = new TF1("bbfit_d",fBBdedx_d,minp,maxp,2,1);
	bbfit[2] = new TF1("bbfit_t",fBBdedx_t,minp,maxp,2,1);
	bbfit[3] = new TF1("bbfit_he3",fBBdedx_he3,minp,maxp,2,1);
	bbfit[4] = new TF1("bbfit_al",fBBdedx_al,minp,maxp,2,1);
	bbfit[5] = new TF1("bbfit_he6",fBBdedx_he6,minp,maxp,2,1);
	bbfit[6] = new TF1("bbfit_li6",fBBdedx_li6,minp,maxp,2,1);
	bbfit[7] = new TF1("bbfit_li7",fBBdedx_li7,minp,maxp,2,1);
	lvfit[0] = new TF1("lvfit_p",fLVdedx_p,minp,maxp,3,1);
	lvfit[1] = new TF1("lvfit_d",fLVdedx_d,minp,maxp,3,1);
	lvfit[2] = new TF1("lvfit_t",fLVdedx_t,minp,maxp,3,1);
	lvfit[3] = new TF1("lvfit_he3",fLVdedx_he3,minp,maxp,3,1);
	lvfit[4] = new TF1("lvfit_al",fLVdedx_al,minp,maxp,3,1);
	lvfit[5] = new TF1("lvfit_he6",fLVdedx_he6,minp,maxp,3,1);
	lvfit[6] = new TF1("lvfit_li6",fLVdedx_li6,minp,maxp,3,1);
	lvfit[7] = new TF1("lvfit_li7",fLVdedx_li7,minp,maxp,3,1);

	auto latex = new TLatex();
	latex->SetTextFont(132);
	latex->SetTextSize(0.06);
	cvsBB->cd();
	latex->DrawLatexNDC(0.12,0.92,"Bethe-Bloch");
	cvsLV->cd();
	latex->DrawLatexNDC(0.12,0.92,"Landau-Vavilov");
	TString particle[] = {"proton","deuteron","triton","3He","alpha","6He","6Li","7Li"};
	latex->SetTextSize(0.05);
	for(Int_t i=0; i<8; i++){
		bbfit[i]->SetParameters(11500.,0.);
		bbfit[i]->SetLineWidth(2);
		bbfit[i]->SetNpx(2000);
		bbfit[i]->SetLineColor(i+2);
		lvfit[i]->SetParameters(20000.,0.,1.);
		lvfit[i]->SetLineWidth(2);
		lvfit[i]->SetNpx(2000);
		lvfit[i]->SetLineColor(i+2);
		latex->SetTextColor(i+2);

		cvsBB->cd();
		latex->DrawLatexNDC(0.12,0.85-i*0.04,particle[i]);
		bbfit[i]->Draw("same");
		cvsLV->cd();
		latex->DrawLatexNDC(0.12,0.85-i*0.04,particle[i]);
		lvfit[i]->Draw("same");
	}
	//cvsBB->SaveAs("pidWithBB.png");
	//cvsLV->SaveAs("pidWithLV.png");

}

