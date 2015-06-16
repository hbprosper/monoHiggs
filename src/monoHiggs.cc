// ---------------------------------------------------------------------------
// mono Higgs analysis (pp -> 4 muons + missing ET)
// created: Les Houches 2015 Nic & HBP
// ---------------------------------------------------------------------------
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>

#include "TStyle.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TClonesArray.h"

#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"
#include "monoHiggs.h"
#include "nic.h"

using namespace std;

// ---------------------------------------------------------------------------
const double MUMASS=0.105658;
const double ZMASS=91.19;

// Quadratic background function
Double_t Dbackground(Double_t *x, Double_t *par)
{
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}
// Lorentzian Peak function
Double_t lorentzianPeak(Double_t *x, Double_t *par)
{
  return (0.5*par[0]*par[1]/TMath::Pi()) /
    TMath::Max(1.e-10,(x[0]-par[2])*(x[0]-par[2])+ .25*par[1]*par[1]);
}
// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
  return Dbackground(x,par) + lorentzianPeak(x,&par[4]);
}

struct DiLepton
{
  DiLepton() {}
  DiLepton(TLorentzVector& L1_, TLorentzVector& L2_)
    : L1(L1_), L2(L2_), LL(L1_+L2_), used(false) {}
  ~DiLepton() {}
  TLorentzVector L1;
  TLorentzVector L2;
  TLorentzVector LL;
  bool used;
};

// ---------------------------------------------------------------------------
// inputFile    space delimited list of Delphes files
// sample       e.g.: "4mu_sig" or "4mu_bkg"
// pileup       mean number of pileup events
// ---------------------------------------------------------------------------
void monoHiggs::analysis4mu(string inputFile,
			    string prefix_,
			    int pileup,
			    int numberEvents){
  cout << endl << "\t== monoHiggs::analysis4mu ==" << endl;

  // -----------------------------------------
  // Objects
  //  muons  - PT > 5 GeV
  //           |Eta| < 2.4
  //           relisol < 0.4
  //
  //  Z1     - min({mZ - Zmass})
  //           mZ1 >  40 GeV
  //           mZ1 < 120 GeV
  //
  //  Z2     - max( {Z.pT} )
  //           mZ2 >  12 GeV
  //           mZ2 < 120 GeV
  //
  //  L1     - PT > 20 GeV
  //  L2     - PT > 10 GeV
  //
  // -----------------------------------------
  // try to load libDelphes
  try {
    gSystem->Load("libDelphes");
  }
  catch (...) {
    cout << endl << "** can't load libDelphes" << endl;
    exit(0);
  }
  
  bool useEventWeight=false;
  char filename[256];
  char prefix[80];
  sprintf(prefix, "%s", prefix_.c_str());

  // make sure plots and histos directories exist
  nic::shell("mkdir -p plots; mkdir -p histos");
  
  // create empty root file for analysis results
  sprintf(filename, "histos/%s_PU%d_results.root", prefix, pileup);
  TFile* theFile = new TFile(filename, "RECREATE");
  if ( ! theFile->IsOpen() ) nic::ciao(string("can't create ") + filename);

  // create a chain of input files
  cout << endl << "=> input files:" << endl;
  TChain* chain = new TChain("Delphes");
  std::vector<std::string> inputFiles = nic::split(inputFile);
  for (size_t i=0; i < inputFiles.size(); i++){
    chain->Add(inputFiles[i].c_str());
    std::cout<<"\tadded input file "<<inputFiles[i]<<std::endl;
  }

  // create A Delphes tree reader
  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
   
  // initialize pointers to branches correspoding to Delphes objects
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchMuon  = treeReader->UseBranch("Muon");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchTower = treeReader->UseBranch("Tower");
  
  // create empty histograms
  gStyle->SetTitleFont(22,"X_mod");
  gStyle->SetTitleFont(22,"Y");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  cout << endl << "=> initialize histograms" << endl;
  TH1D *h_nEvent = new TH1D("h_nEvent", "Delphes Eff vs cut", 11, -1, 10);
  TH1F *weight = new TH1F("weight","weight",11,-1,10);

  // histogram for isolation cuts
  int numIsoCuts = 10;  
  TH1F* h_isolation = new TH1F("isolation","", numIsoCuts, 1, numIsoCuts);
  for(int k=0; k < numIsoCuts; k++)
    {
      double relisol = 0.1 + k * 0.1;
      char label[80];
      sprintf(label, "> %3.1f", relisol);
      h_isolation->GetXaxis()->SetBinLabel(k+1, label);
    }
  
  TH1F* h_smweight= new TH1F("smweight","Snowmass Event Weight",200,0,200000);
  TH1F* h_nMuons = new TH1F("nMuons","Number of muons",10,0,10);
  
  TH2F* ETA_PTgen = new TH2F("Eta_PTgen","Eta vs PT",100,-5.,5.,200,0.,100.);
  TH2F* ETA_PT = new TH2F("Eta_PT","Eta vs PT",100,-5.,5.,200,0.,100.);
  TH2F* ETA_PTeff = new TH2F("Eta_PTeff","Eta vs PT",100,-5.,5.,200,0.,100.);
  TH1F* ETA = new TH1F("Eta","Eta",100,-5.,5.);
  TH1F* PT = new TH1F("PT","PT",200,0.,100.);
  TH1F* ETAgen = new TH1F("Etagen","Etagen",100,-5.,5.);
  TH1F* PTgen = new TH1F("PTgen","PTgen",200,0.,100.);
  TH1F* ETAres = new TH1F("ETAres","ETAres",100,-1.,1.);
  TH1F* PTres = new TH1F("PTres","PTres",100,-.15,.15);
  TH1F* ETAeff = new TH1F("Etaeff","Etaeff",100,-5.,5.);
  TH1F* PTeff = new TH1F("PTeff","PTeff",200,0.,100.);
  TH1F* ETAeff2 = new TH1F("Etaeff2"," Eff_{#eta} (#mu^{+} #mu^{-})",
			   100,-5.,5.);
  TH1F* PTeff2 = new TH1F("PTeff2"," Eff_{p_{T}} (#mu^{+} #mu^{-})",
			  200,0.,100.);
  TH1F* ETA2 = new TH1F("Eta2","Etaeff2",100,-5.,5.);
  TH1F* PT2 = new TH1F("PT2","PTeff2",200,0.,100.);
  TH1F* muEta = new TH1F("muEta","muEta",100,-5.,5.);
  TH1F* muPT = new TH1F("muPT","muPT",200,0.,100.);	
   
  TH1F* Z1mass = new TH1F("Z1mass","Z1mass",200,20.,120.);
  TH1F* Z2mass = new TH1F("Z2mass","Z2mass",200,0.,100.);
  TH1F* Hmass = new TH1F("Hmass","Hmass",200,50.,150.);
  TH1F* HmassLarge = new TH1F("HmassLarge","HmassLarge",600,0.,300.);
   
  PTeff2->SetLineColor(3);
  PTeff2->SetLineWidth(2);
  PTeff2->SetTitleOffset(1.1,"X");
  PTeff2->GetXaxis()->SetTitle("p_{T}");
   
  ETAeff2->SetTitleOffset(1.1,"X");
  ETAeff2->GetXaxis()->SetTitle("#eta");
  ETAeff2->SetLineWidth(2);
  ETAeff2->SetLineColor(3);	
   
  ETA_PTeff->GetXaxis()->SetTitle("#eta");
  ETA_PTeff->GetYaxis()->SetTitle("p_{T}");

  // -----------------------------------------
  // EVENT LOOP
  // -----------------------------------------
  cout << endl << "=> begin event loop" << endl;
  bool eff2 = false;
  bool original = false;
  bool efficiencies = true;
  bool masses = true;
  
  TLorentzVector L1;
  TLorentzVector L2;
  TLorentzVector L3;
  TLorentzVector L4;
  TLorentzVector Z1;
  TLorentzVector Z2;
  TLorentzVector H;
  TLorentzVector fourL;

  int index[200];
  int index1[200];
  int badL[200];

  long int numberOfEntries = treeReader->GetEntries();
  if ( numberEvents > 0 ) numberOfEntries = numberEvents;
  std::cout <<"\tanalyzing "<<numberOfEntries<<" events" << std::endl;

   for(Int_t entry = 0; entry < numberOfEntries; entry++){
 
    //Load branches for event
    treeReader->ReadEntry(entry);
    if (entry%100==0) std::cout<<"\tprocessing event "<<entry<<std::endl;
    
    // get event weight
    double smweight = 1.0;
    if ( useEventWeight )
      {
	LHEFEvent* event = (LHEFEvent*)branchEvent->At(0);
	smweight = event->Weight / numberOfEntries;
      }
    h_smweight->Fill(smweight);
     
    // histogram number of reco muons in the event
    h_nMuons->Fill(branchMuon->GetEntriesFast());
    
    // histogram some gen level muon quantities and
    // count number of gen level muons.
    int count = 0;		
    for(int i = 0; i < branchParticle->GetEntriesFast(); i++){
      GenParticle* particle = static_cast<GenParticle*>
	(branchParticle->At(i));
      if(fabs(particle->PID) == 13 && particle->Status == 1){
	ETAgen->Fill(particle->Eta);
	PTgen->Fill(particle->PT);
	ETA_PTgen->Fill(particle->Eta, particle->PT);
	count++;
      }
    }
    // 1. first bin of h_nEvent (with value -1) contains the total
    // event weight
    // 2. all bins of weight contains total count if count > 3
    if( count > 3 ) { 
      h_nEvent->Fill(-1.0, smweight);
      for(int k = -1; k < 10; k++) weight->Fill(k);
    }		
     
    // compare reco and gen level muons and, in particular, histogram
    // the muon deltapT / pT, where deltapT = muon-pT - pT and pT is
    // the true pT.
    for(int i=0; i<branchParticle->GetEntriesFast();i++){
      GenParticle* particle = static_cast<GenParticle*>
	(branchParticle->At(i));
      if(fabs(particle->PID) == 13 && particle->Status == 1){
	for(int k = 0; k < branchMuon->GetEntriesFast();k++){
	  Muon* muon = static_cast<Muon*>(branchMuon->At(k));
	  GenParticle* particle2 = static_cast<GenParticle*>
	    (muon->Particle.GetObject());
	  if(particle2 == particle) {
	    PT2->Fill(particle->PT);
	    ETA2->Fill(particle->Eta);
	    muPT->Fill(muon->PT);
	    muEta->Fill(muon->Eta);
	    ETA_PT->Fill(particle->Eta,particle->PT);
	    PTres->Fill((muon->PT - particle->PT)/particle->PT);
	  } 
	}
      }
    }

    // ----------------------------------------------------------
    // flag bad muons: muon pairs that are too close together in
    // deltaR or that have pT > 5 GeV and lie outside the eta
    // acceptance of the muon system.
    // ----------------------------------------------------------             
    int numberBad = 0;
    for(int i = 0; i < branchMuon->GetEntriesFast();i++){
      Muon* muon = static_cast<Muon*>(branchMuon->At(i));

      if(muon->PT < 5 ||  fabs(muon->Eta) > 2.4){
	badL[numberBad] = i;
	numberBad++;
	continue;
      }
      
      for(int j = i+1; j < branchMuon->GetEntriesFast(); j++){
	Muon* muon1 = static_cast<Muon*>(branchMuon->At(j));
	 
	double dR = nic::deltaR(muon->Eta, muon->Phi,
				muon1->Eta, muon1->Phi);
	if(dR < 0.02){
	  std::cout<<"Dropped closeby muon candidates"<<std::endl;
	  badL[numberBad] = i;
	  badL[numberBad+1] = j;
	  numberBad += 2;
	}
      }
    }

    // get list of good muons
    // Ln1 is the number of good muons
    int Ln1 = 0;
    for(int i = 0; i < branchMuon->GetEntriesFast(); i++){	
      bool good = true;
      for(int k = 0; k < numberBad; k++){
	if(badL[k] == i) {
	  good = false;
	  break;
	}
      }
      if( good ){
	index1[Ln1] = i;
	Ln1++;
      }
    }
  
    // get the list of good isolated muons 
    // Ln is number of isolated muons
    int Ln = 0;
    for(int k = 0; k < Ln1; k++){
      Muon* muon = (Muon*) branchMuon->At(index1[k]);	
      double relisol = nic::leptonIsolation(muon, branchTrack, branchTower); 
      if( relisol < 0.4 ){
	index[Ln] = index1[k];
	Ln++;
      }
    }
    
    // ----------------------------------------------------------
    // CUT 2: number of isolated muons > 3
    // ----------------------------------------------------------
    if ( !(Ln > 3) )  continue;
    // Note: bin 3 (value 1) counts the weighted number of
    // good isolated muons.
    h_nEvent->Fill(0.0, smweight);
    
    // find all opposite sign di-leptons
    vector<DiLepton> dilepton;
    for(int k = 0; k < Ln; k++){
      Muon* muon = static_cast<Muon*>(branchMuon->At(index[k]));
      for(int j = k+1; j< Ln; j++){
	Muon* muon1 = static_cast<Muon*>(branchMuon->At(index[j]));
	if( !(muon->Charge*muon1->Charge < 0) ) continue;

	// we have oppositely charged muons
	L1.SetPtEtaPhiM(muon->PT,muon->Eta,muon->Phi,MUMASS);
	L2.SetPtEtaPhiM(muon1->PT,muon1->Eta,muon1->Phi,MUMASS);

	// QCD suppression
	Z1 = L1 + L2;
	if ( !(Z1.M() > 4) ) continue;
	
	dilepton.push_back(DiLepton(L1, L2));
      }
    }

    // ----------------------------------------------------------
    // CUT 3: number of dileptons > 0
    // ----------------------------------------------------------
    if ( !(dilepton.size() > 0) )  continue;
    h_nEvent->Fill(1.0, smweight);

    // find Z candidate closest to Z pole mass
    int indexZ1 = 0;
    double error = 1e10;
    double largestPT = 0;
    for(size_t k = 0; k < dilepton.size(); k++){
      DiLepton& d = dilepton[k];
      if(fabs(d.LL.M() - ZMASS) < error){
	indexZ1 = k;
	error = fabs(d.LL.M() - ZMASS);
      }
    }
    
    // flag this dilepton as "used"
    dilepton[indexZ1].used = true;
    Z1 = dilepton[indexZ1].LL;
    L1 = dilepton[indexZ1].L1;
    L2 = dilepton[indexZ1].L2;    

    // ----------------------------------------------------------
    // CUT 4: 40 < Z1mass < 120
    // ----------------------------------------------------------
    if( !(Z1.M() > 40) )   continue;
    if( !(Z1.M() < 120) )  continue;
    h_nEvent->Fill(2.0, smweight);

    // ----------------------------------------------------------
    // CUT 5: number of potential Z candidates > 1
    // ----------------------------------------------------------
    if( !( dilepton.size() > 1) )  continue;
    h_nEvent->Fill(3.0, smweight);    

    // this will be true if we have di-Z event
    bool diZevent = false;

    // look for another Z candidate
    for(size_t k = 0; k < dilepton.size(); k++){
      if ( dilepton[k].used ) continue;

      DiLepton& d = dilepton[k];
      if ( ! (d.LL.M() > 0) ) continue;
      
      double dR = nic::deltaR(d.LL.Eta(), d.LL.Phi(),
			      Z1.Eta(), Z1.Phi());
      if ( ! (dR > 0.02) )   continue;

      if(d.LL.Pt() > largestPT){
	largestPT  = d.LL.Pt();
	Z2 = d.LL;
	L3 = d.L1;
	L4 = d.L2;
	diZevent = true;
      }
    }
    // ----------------------------------------------------------
    // CUT 6: number of di-Z bosons == 1
    // ----------------------------------------------------------
    if( ! diZevent )  continue;
    h_nEvent->Fill(4.0, smweight);
    
    // ----------------------------------------------------------
    // CUT 7: 12 < Z2mass < 120
    // ----------------------------------------------------------
    if( !(Z2.M() >  12) ) continue;
    if( !(Z2.M() < 120) ) continue;
    h_nEvent->Fill(5.0, smweight);    
    
    // We now have a di-Z event
    std::cout << " event " << entry
	      << "\tZ1mass=" << Z1.M()
	      << "\tZ2mass=" << Z2.M() << std::endl;   
    // ----------------------------------------------------------
    // sort lepton PTs
    // ----------------------------------------------------------
    vector<double> LPT(4,0);
    LPT[0] = L1.Pt();
    LPT[1] = L2.Pt();
    LPT[2] = L3.Pt();
    LPT[3] = L4.Pt();  
    sort(LPT.begin(), LPT.end());
    double LPT1 = LPT[3];
    double LPT2 = LPT[2];
    
    // ----------------------------------------------------------
    // CUT 8: LPT1 > 20 GeV
    // ----------------------------------------------------------
    if ( !(LPT1 > 20) ) continue;
    h_nEvent->Fill(6.0, smweight);
    
    // ----------------------------------------------------------
    // CUT 9: LPT2 >10 GeV
    // ----------------------------------------------------------    
    if ( !(LPT2 > 10) ) continue;
    h_nEvent->Fill(7.0, smweight);
    
    // compute 4-lepton 4-vector
    H = Z1 + Z2;
    
    // ----------------------------------------------------------
    // CUT 10: m4l > 100 GeV
    // ----------------------------------------------------------    
    if ( !( H.M() > 100) ) continue;
    h_nEvent->Fill(8.0, smweight);
    
    Z1mass->Fill(Z1.M(), smweight);
    Z2mass->Fill(Z2.M(), smweight);
    Hmass->Fill(H.M(), smweight);
    HmassLarge->Fill(H.M(), smweight);
   } // END OF EVENT LOOP

   // -----------------------------------------------------------
   // Finish up
   // -----------------------------------------------------------
   if ( original ){	
     TCanvas* can = new TCanvas("can","",0,0,1000,500);
     can->Divide(3,4);
     ETAeff->Divide(ETA,ETAgen);
     PTeff->Divide(PT,PTgen);

     can->cd(1);
     ETA->Draw();
     can->cd(2);
     PT->Draw();
     can->cd(3);

     can->cd(4);
     ETAgen->Draw();
     can->cd(5);
     PTgen->Draw();
     can->cd(6);

     can->cd(7);
     ETAres->Draw();
     can->cd(8);
     PTres->Draw();
     can->cd(9);

     can->cd(10);
     ETAeff->Draw();
     can->cd(11);
     PTeff->Draw();
     can->cd(12);
   }
   
   if ( eff2 ){
     TCanvas* c3 = new TCanvas("c3");
     c3->Divide(3,1);
     PTeff2->Divide(PT2,PTgen);
     ETAeff2->Divide(ETA2,ETAgen);
     
     c3->cd(1);
     ETAeff2->Draw();
     c3->cd(2);
     PTeff2->Draw();
     c3->cd(3);
   }
   
   if ( efficiencies ){
     TCanvas* c = new TCanvas("c");
     c->Divide(1,4);
     ETA_PTeff->Divide(ETA_PT,ETA_PTgen);
     PTeff2->Divide(PT2,PTgen);
     ETAeff2->Divide(ETA2,ETAgen);
     c->cd(1); PTeff2->Draw();
     c->cd(2); ETAeff2->Draw();
     c->cd(3); ETA_PTeff->Draw("colz");

     TCanvas *cuts = new TCanvas();
     sprintf(filename, "plots/%s_PU%d_pT_Eff.png", prefix, pileup);
     cuts->Print(filename, "png");
     
     sprintf(filename, "plots/%s_PU%d_Eta_Eff.png", prefix, pileup);
     cuts->Print(filename, "png");
     c->cd(4); PTres->Draw();
   }
   
   if ( masses ){
     TCanvas *m1 = new TCanvas();
     TCanvas *m2 = new TCanvas();
     TCanvas *m3 = new TCanvas();
     TCanvas *m4 = new TCanvas();
     
     m1->cd(); 
     Z1mass->Draw(); 
     
     m2->cd(); 
     Z2mass->Draw(); 
     
     m3->cd(); 
     Hmass->Draw();

     m4->cd(); 
     HmassLarge->Draw();     

     TCanvas *cuts = new TCanvas();
     cuts->cd();
     TH1F *h_nEvent1 = new TH1F("h_nEvent1","Delphes Eff vs cut",11,-1,10);
     h_nEvent1->Sumw2();
     h_nEvent1->Divide(h_nEvent, weight);
     
     h_nEvent1->SetLineWidth(2);
     h_nEvent1->SetLineColor(3);
     
     h_nEvent1->GetXaxis()->SetBinLabel(1,"N events");
     h_nEvent1->GetXaxis()->SetBinLabel(2,"n(isolepton) #ge 4");
     h_nEvent1->GetXaxis()->SetBinLabel(3,"n(dilepton) #ge 1");
     h_nEvent1->GetXaxis()->SetBinLabel(4,"40 < Z1(1 SFOS pair)  < 120");
     h_nEvent1->GetXaxis()->SetBinLabel(5,"n(dilepton) #ge 2");
     h_nEvent1->GetXaxis()->SetBinLabel(6,"diZevent");
     h_nEvent1->GetXaxis()->SetBinLabel(7,"12 < Z2(+1 SFOS pair)  < 120");
     h_nEvent1->GetXaxis()->SetBinLabel(8,"pT(l1) > 20 GeV");
     h_nEvent1->GetXaxis()->SetBinLabel(9,"pT(l2) > 10 GeV");
     h_nEvent1->GetXaxis()->SetBinLabel(10,"M(ll) > 100 GeV (QCD)");
     
     cuts->cd(); 
     h_nEvent1->SetMinimum(0.0);
     h_nEvent1->SetMaximum(1.1);
     
     h_nEvent1->Draw();
     sprintf(filename, "plots/%s_PU%d_cutflow.png", prefix, pileup);
     cuts->Print(filename, "png");
   }
   
   theFile->cd();
   theFile->Write();
}
