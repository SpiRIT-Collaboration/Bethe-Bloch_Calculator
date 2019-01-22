
TString pathReco = "/data/Q18393/production/20181219/data/";
TString recoVer  = "HEAD.1769.ef17b59";

void makeGoodPID()
{
	auto chain = new TChain("cbmsim");
	
	for(Int_t runNo = 2900; runNo<2910; runNo++){
		auto splitNo = 0;
		while(gSystem->IsFileInIncludePath(Form(pathReco+"run%d_s%d.reco."+recoVer+".root",runNo,splitNo))){
			chain->Add(Form(pathReco+"run%d_s%d.reco."+recoVer+".root",runNo,splitNo));
			splitNo++;
		}
	}

	chain->SetBranchStatus("*",0);
	chain->SetBranchStatus("STVertex*",1);
	chain->SetBranchStatus("VATracks*",1);
	TClonesArray *vertexArray  = nullptr;
	TClonesArray *vaTrackArray = nullptr;
	chain->SetBranchAddress("STVertex",&vertexArray);
	chain->SetBranchAddress("VATracks",&vaTrackArray);

	auto outFile = new TFile("goodPID.root","recreate");
	auto h2PID   = new TH2D("h2PID","",500,-500,2500,500,10,1500);

	auto nEvent = chain->GetEntries();
	for(auto i: ROOT::TSeqL(nEvent)){
		chain->GetEntry(i);
		
		if(vertexArray->GetEntries()<=0) continue;
		
		auto v = (STVertex*) vertexArray->At(0);
		auto pos = v->GetPos();
		if(TMath::Abs(pos.Z()+13.3)>3.5) continue;
		if(TMath::Abs(pos.X())>15||TMath::Abs(pos.Y()+226.)>20) continue;

		auto nTrack = vaTrackArray->GetEntries();
		for(auto j: ROOT::TSeqL(nTrack)){
			auto track = (STRecoTrack*) vaTrackArray->At(j);
			
			if(track->GetMomentumTargetPlane().Mag()==0||track->GetCharge()==0) continue;
			Double_t dist = (pos-track->GetPOCAVertex()).Mag();
			Int_t ndf     = track->GetClusterIDArray()->size();
			if(dist>10||ndf<20)continue;
			
			TVector3 rig = track->GetMomentumTargetPlane();
			if(rig.Z()<0) rig = -rig;
			Double_t momTarget = rig.Mag();
			Double_t theta = rig.Theta()*180./TMath::Pi();
			Double_t phi   = rig.Phi()*180./TMath::Pi();

			if(theta<30&&TMath::Abs(phi)<40)
				h2PID->Fill(momTarget/track->GetCharge(),track->GetdEdxWithCut(0.,0.7));
		}

			
	}

	outFile->Write();
	outFile->Close();

}
