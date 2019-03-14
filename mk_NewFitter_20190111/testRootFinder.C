#include "MassEstimator.h"
using namespace MassEstimator;

void testRootFinder()
{
	// parameters for the eq. function.
	Double_t *p = new Double_t[8];
	p[0] = 1.;
	p[1] = 0.;
	p[2] = 1.;
	p[3] = 0.;
	p[4] = 12.;
	p[5] = 1.;

	auto out = new TFile("test.root","recreate");
	auto tree = new TTree("tree","");
	Double_t mom,dedx,mass;
	Double_t calMass[7];  // calculated mass for different options.

	tree->Branch("mom",&mom,"mom/D");
	tree->Branch("dedx",&dedx,"dedx/D");
	tree->Branch("mass",&mass,"mass/D");
	tree->Branch("calMass",calMass,"calMass[7]/D");

	auto finder = new ROOT::Math::RootFinder;
	Double_t dx = 0.01;  // minute displacement for the derivation.
	auto func  = [p](double x){ return LVMassFinderEq(&x,p); };
	auto dfunc = [p,dx](double x){ return LVMassFinderDeriv(&x,p,dx);};
	auto func1d  = new ROOT::Math::Functor1D(func);
	auto gradf1d = new ROOT::Math::GradFunctor1D(func,dfunc);

	for(Int_t i=0; i<50000; i++){

		mom = mass = dedx = -1.;
		for(Int_t j=0; j<7; j++)calMass[j]=-1.;

		while(mom<=0)  mom = gRandom->Landau(200,100);
		while(mass<=0) mass = 938.7*gRandom->Gaus(1.,0.2);
		dedx = fLVdedx(1.,mass,&mom,p);

		p[6] = mom;
		p[7] = dedx;


		// kBRENT
		finder->SetMethod(ROOT::Math::RootFinder::kBRENT);
		finder->SetFunction(*func1d,0.1,10000.);
		finder->Solve();
		calMass[0] = finder->Root();

		// These three make a lot of error messages (when calculation failed).
		// If you want to turn off the error or info. messages, set global variable like below.
		// gErrorIgnoreLevel = kBreak, kSysError, kFatal;
		/*
		// BISECTION
		finder.SetMethod(ROOT::Math::RootFinder::kGSL_BISECTION);
		finder.SetFunction(*func1d,0.1,10000.);
		finder.Solve();
		calMass[1] = finder.Root();

		// FALSE_POS
		finder.SetMethod(ROOT::Math::RootFinder::kGSL_FALSE_POS);
		finder.SetFunction(*func1d,0.1,10000.);
		finder.Solve();
		calMass[2] = finder.Root();

		// GSL_BRENT
		finder.SetMethod(ROOT::Math::RootFinder::kGSL_BRENT);
		finder.SetFunction(*func1d,0.1,10000.);
		finder.Solve();
		calMass[3] = finder.Root();
		*/

		// need derivative !
		// Newton
		finder->SetMethod(ROOT::Math::RootFinder::kGSL_NEWTON);
		finder->SetFunction(*gradf1d,1000.);
		finder->Solve();
		calMass[4] = finder->Root();

		// Secant
		finder->SetMethod(ROOT::Math::RootFinder::kGSL_SECANT);
		finder->SetFunction(*gradf1d,1000.);
		finder->Solve();
		calMass[5] = finder->Root();

		// Steffenson
		finder->SetMethod(ROOT::Math::RootFinder::kGSL_STEFFENSON);
		finder->SetFunction(*gradf1d,1000.);
		finder->Solve();
		calMass[6] = finder->Root();


		tree->Fill();
	}

	out->Write();

}

