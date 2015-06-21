// ---------------------------------------------------------------------------
// mono Higgs analysis (pp -> 4 leptons + missing ET)
// created: Les Houches 2015 Nic & HBP
// ---------------------------------------------------------------------------
#include <algorithm>
#include <iostream>
#include <fstream>
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
//const double MUMASS=0.105658;
const double ZMASS=91.19;
const int k4MU   = 1;
const int k4E    = 2;
const int k2E2MU = 3;
const int ELECTRON = 11;
const int MUON = 13;
/** Cone-based isolation variable.
    isoType
    1         CMS muon isolation
    2         CMS electron isolation
    3         CMS photon isolation
    
    11        ATLAS muon isolation
    12        ATLAS electron isolation
    13        ATLAS photon isolation
*/
const int CMS_MU_ISOL=1;
//const int CMS_EL_ISOL=2;
// ---------------------------------------------------------------------------
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

struct TParticle : public TLorentzVector
{
  TParticle() : TLorentzVector() {}
  TParticle(int PID_, 
	    double PT, double Eta, double Phi, double Mass=0)
    : TLorentzVector(),
      PID(PID_),
      ID(abs(PID_)),
      Name(nic::particleName(PID_)),
      Value(std::map<std::string, double>()),
      particle(0)
  {
    SetPtEtaPhiM(PT, Eta, Phi, Mass);
  }
  bool operator<(const TParticle& o) const { return o.Pt() < this->Pt(); }
  
  int PID; // PDG ID
  int ID;  // unsigned PDG ID
  std::string Name;
  std::map<std::string, double> Value;
  GenParticle* particle;
};
std::ostream& operator<<(std::ostream& os, const TParticle& o)
{
  char record[80];
  sprintf(record, "%10d %-12s\t%10.2f %10.3f %10.3f",
	  o.PID, o.Name.c_str(), o.Pt(), o.Eta(), o.Phi());
  os << record;
  return os;
}

void copyData(TClonesArray* electrons,
	      TClonesArray* muons,
	      std::vector<TParticle>& lepton)
{
  if ( electrons )
    {
      for(int i=0; i < electrons->GetEntriesFast(); i++)
	{
	  Electron* p = static_cast<Electron*>(electrons->At(i));
	  int PID = -ELECTRON * p->Charge;
	  lepton.push_back(TParticle(PID, p->PT, p->Eta, p->Phi));
	  lepton.back().particle=static_cast<GenParticle*>(p->Particle.GetObject());
	}
    }

  if ( muons )
    {
      for(int i=0; i < muons->GetEntriesFast(); i++)
	{
	  Muon* p = static_cast<Muon*>(muons->At(i));
	  int PID = -MUON * p->Charge;
	  lepton.push_back(TParticle(PID, p->PT, p->Eta, p->Phi));
	  lepton.back().particle=static_cast<GenParticle*>(p->Particle.GetObject());
	}
    }
  sort(lepton.begin(), lepton.end());
}

bool foundLepton(int finalState, int ID)
{
  bool found = false;
  switch (finalState)
    {
    case k4MU:
      found = ID == 13;
      break;
    case k4E:
      found = ID == 11;
      break;
    case k2E2MU:
      found = (ID == 11) || (ID == 13);
      break;
    }
  return found;
}
// ---------------------------------------------------------------------------
// inputFile    input file (with .root extension) of a filelist
// prefix       e.g.: "sig_4mu" or "bkg_4mu"
// pileup       mean number of pileup events
// numberEvents obvious, no?!
// ---------------------------------------------------------------------------
void monoHiggs::analysis(string inputFile,
			 string prefix_,
			 int pileup,
			 int numberEvents){
  cout << endl << "\t== monoHiggs::analysis ==" << endl;

  // get final state 
  int finalState = 0;
  if      ( prefix_.find("4mu")  != std::string::npos )
    finalState = k4MU;
  else if ( prefix_.find("4e")   != std::string::npos )
    finalState = k4E;
  else if ( prefix_.find("2e2mu")!= std::string::npos )
    finalState = k2E2MU;
  else
    nic::ciao("unrecognized final state (need 4mu, 4e, or 2e2mu)");
  
  // -----------------------------------------
  // analysis switches
  // -----------------------------------------
  bool eff2 = false;
  bool original = false;
  bool efficiencies = true;
  bool masses = true;
  bool useEventWeight=false;
  
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

  // inputFile could be either a root file or a filelist
  std::vector<std::string> inputFiles;
  if ( inputFile.find(".root") != std::string::npos )
    inputFiles.push_back(inputFile);
  else
    {
      // assume this is a filelist
      ifstream inp(inputFile);
      if ( !inp.good() ) nic::ciao("can't open " + inputFile);
      string infilename;
      while (getline(inp, infilename))
	{
	  if ( infilename.substr(0,1) == "#" ) continue;
	  if ( infilename.substr(0,1) == " " ) continue;
	  if ( infilename.substr(0,1) == "\n" ) continue;
	  inputFiles.push_back(infilename);
	}
    }
  for (size_t i=0; i < inputFiles.size(); i++){
    chain->Add(inputFiles[i].c_str());
    std::cout<<"\tadded input file "<<inputFiles[i]<<std::endl;
  }

  // create A Delphes tree reader
  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
   
  // initialize pointers to branches correspoding to Delphes objects
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");

  TClonesArray *branchElectron = 0;
  TClonesArray *branchMuon = 0;
  switch (finalState)
    {
    default:
    case k4MU:
      branchMuon = treeReader->UseBranch("Muon");
      break;
    case k4E:
      branchElectron = treeReader->UseBranch("Electron");
      break;
    case k2E2MU:
      {
	branchElectron = treeReader->UseBranch("Electron");
	branchMuon = treeReader->UseBranch("Muon");
      }
      break;
    }
      
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchTower = treeReader->UseBranch("Tower");
  TClonesArray *branchRho   = treeReader->UseBranch("Rho");
  TClonesArray *branchMET   = treeReader->UseBranch("MissingET");

  // create empty histograms
  gStyle->SetTitleFont(22,"X_mod");
  gStyle->SetTitleFont(22,"Y");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  cout << endl << "=> initialize histograms" << endl;
  TH1F *h_nEvent = new TH1F("h_nEvent", "Cut flow", 10, 0, 10);
  h_nEvent->Sumw2();
  h_nEvent->SetLineWidth(2);
  h_nEvent->SetLineColor(kBlue);
  
  h_nEvent->GetXaxis()->SetBinLabel(1,"N events");
  h_nEvent->GetXaxis()->SetBinLabel(2,"N(isolepton) #ge 4");
  h_nEvent->GetXaxis()->SetBinLabel(3,"N(dilepton) #ge 1");
  h_nEvent->GetXaxis()->SetBinLabel(4,"40 < Z1(1 SFOS pair)  < 120");
  h_nEvent->GetXaxis()->SetBinLabel(5,"N(dilepton) #ge 2");
  h_nEvent->GetXaxis()->SetBinLabel(6,"diZevent");
  h_nEvent->GetXaxis()->SetBinLabel(7,"12 < Z2(+1 SFOS pair)  < 120");
  h_nEvent->GetXaxis()->SetBinLabel(8,"pT(l1) > 20 GeV");
  h_nEvent->GetXaxis()->SetBinLabel(9,"pT(l2) > 10 GeV");
  h_nEvent->GetXaxis()->SetBinLabel(10,"M(ll) > 100 GeV (QCD)");
  
  // histogram for isolation cuts
  int numIsoCuts = 10;  
  TH1F* h_isolation = new TH1F("isolation", "", numIsoCuts, 0, numIsoCuts);
  h_isolation->Sumw2();
  for(int k=0; k < numIsoCuts; k++)
    {
      double detisol = 0.1 + k * 0.1;
      char label[80];
      sprintf(label, "> %3.1f", detisol);
      h_isolation->GetXaxis()->SetBinLabel(k+1, label);
    }
  
  TH1F* h_detisol = new TH1F("detisol", "isolation", 80, -2.0, 2.0);
  TH1F* h_rho = new TH1F("rho", "#rho", 50, 0, 100);
  TH1F* h_nLeptons = new TH1F("nLeptons","Number of leptons", 10, 0, 10);
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
  TH1F* ETAeff2 = new TH1F("Etaeff2"," Eff_{#eta} (#l^{+} #l^{-})",
			   100,-5.,5.);
  TH1F* PTeff2 = new TH1F("PTeff2"," Eff_{p_{T}} (#l^{+} #l^{-})",
			  200,0.,100.);
  TH1F* ETA2 = new TH1F("Eta2","Etaeff2",100,-5.,5.);
  TH1F* PT2 = new TH1F("PT2","PTeff2",200,0.,100.);
  
  TH1F* lepEta = new TH1F("lepEta","lepEta",100,-5.,5.);
  TH1F* lepPT  = new TH1F("lepPT","lepPT",200,0.,100.);	
 
  TH1F* Z1mass = new TH1F("Z1mass","Z1mass", 200, 20.,120.);
  TH1F* Z2mass = new TH1F("Z2mass","Z2mass", 200,  0.,100.);
  TH1F* Hmass = new TH1F("Hmass","Hmass",    200, 50.,150.);
  TH1F* HmassLarge = new TH1F("HmassLarge","HmassLarge", 600, 0., 300.);
  TH1F* MET  = new TH1F("MET","Missing E_{T}", 400, 0., 200.);
  
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
  
  TLorentzVector L1;
  TLorentzVector L2;
  TLorentzVector L3;
  TLorentzVector L4;
  TLorentzVector Z1;
  TLorentzVector Z2;
  TLorentzVector H;

  int index[200];
  int index1[200];
  int badL[200];

  long int numberOfEntries = treeReader->GetEntries();
  if ( numberEvents > 0 ) numberOfEntries = numberEvents;
  std::cout <<"\tanalyzing "<<numberOfEntries<<" events" << std::endl;

  // ----------
  // ===> START
  // ----------
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

    // -----------------------------------------------------
    // first bin of h_nEvent contains the total event weight
    // -----------------------------------------------------
    h_nEvent->Fill(0.1, smweight);

    // -----------------------------------------------------
    // copy data to a standard name objects
    // -----------------------------------------------------
    vector<TParticle> lepton;
    copyData(branchElectron, branchMuon, lepton);
    h_nLeptons->Fill(lepton.size());
    
    // histogram some gen level lepton quantities
    for(int i = 0; i < branchParticle->GetEntriesFast(); i++){
      GenParticle* particle = static_cast<GenParticle*>
	(branchParticle->At(i));
      if ( particle->Status == 1)
	{
	  bool found = foundLepton(finalState, abs(particle->PID));
	  if( found ){
	    ETAgen->Fill(particle->Eta);
	    PTgen->Fill(particle->PT);
	    ETA_PTgen->Fill(particle->Eta, particle->PT);
	  }
	}
    }
    
    // compare reco and gen level leptons and, in particular, histogram
    // the lepton deltapT / pT, where deltapT = lepton-pT - pT and pT is
    // the true pT.
    for(int i=0; i<branchParticle->GetEntriesFast();i++){
      GenParticle* particle = static_cast<GenParticle*>
	(branchParticle->At(i));
      if ( !(particle->Status == 1) ) continue;
      
      bool found = foundLepton(finalState, abs(particle->PID));
      if ( !found ) continue;
      
      // we have either a gen-level muon or an electron
      // come true pT with reco-level pT
      for(size_t k = 0; k < lepton.size(); k++){
	TParticle& reco = lepton[k];
	if ( !(reco.particle == particle) ) continue;
	
	PT2->Fill(particle->PT);
	ETA2->Fill(particle->Eta);
	ETA_PT->Fill(particle->Eta, particle->PT);
	
	lepPT->Fill(reco.Pt());
	lepEta->Fill(reco.Eta());
	double dPt = reco.Pt() - particle->PT;
	PTres->Fill(dPt/particle->PT);
      }
    }

    // ----------------------------------------------------------
    // flag bad leptons: lepton pairs that are too close together in
    // deltaR or that have pT > 5 GeV and lie outside the eta
    // acceptance.
    // ----------------------------------------------------------             
    int numberBad = 0;
    for(size_t i = 0; i < lepton.size(); i++){
      TParticle& pi = lepton[i];

      // cuts differ for electrons and muons
      double PtCut;
      double EtaCut;
      if ( pi.ID == MUON )
	{
	  PtCut  = 5.0;
	  EtaCut = 2.4;
	}
      else
	{
	  PtCut  = 7.0;
	  EtaCut = 2.5;
	}

      if ( (pi.Pt() < PtCut) ||  (fabs(pi.Eta()) > EtaCut) ){
	badL[numberBad] = i;
	numberBad++;
	continue;
      }

      for(size_t j = i+1; j < lepton.size(); j++){
	TParticle& pj = lepton[j];
	double dR = nic::deltaR(pi.Eta(), pi.Phi(),
				pj.Eta(), pj.Phi());
	if(dR < 0.02){
	  std::cout<<"Dropped closeby muon candidates"<<std::endl;
	  badL[numberBad] = i;
	  badL[numberBad+1] = j;
	  numberBad += 2;
	}
      }
    }

    // get list of good leptons
    // Ln1 is the number of good leptons
    int Ln1 = 0;
    for(size_t i = 0; i < lepton.size(); i++){	
      bool good = true;
      for(int k = 0; k < numberBad; k++){
	if(badL[k] == (int)i) {
	  good = false;
	  break;
	}
      }
      if( good ){
	index1[Ln1] = i;
	Ln1++;
      }
    }
  
    // get the list of good isolated leptons
    // Ln is number of isolated leptons
    int Ln = 0;
    for(int k = 0; k < Ln1; k++){
      TParticle& p = lepton[index1[k]]; // make a reference not a copy

      // for now use same isolation for muons and electrons
      int isoType = CMS_MU_ISOL;
      double detisol = nic::leptonIsolation(p.Pt(), p.Eta(), p.Phi(),
					    branchTrack,
					    branchTower,
					    branchRho,
					    isoType);
      for(int j=0; j < numIsoCuts; j++)
	{
	  double isol = 0.1 + j * 0.1;
	  if ( detisol > isol )
	    h_isolation->Fill(isol+0.05, smweight);
	}
  
      double rho = static_cast<Rho*>(branchRho->At(0))->Rho;
      h_rho->Fill(rho, smweight);
      h_detisol->Fill(detisol, smweight);

      // need to optimize cut
      if( detisol < 2.0 ){
	index[Ln] = index1[k];
	Ln++;
      }
    }
    
    // ----------------------------------------------------------
    // CUT 2: number of isolated muons >= 4
    // ----------------------------------------------------------
    if ( !(Ln > 3) )  continue;
    h_nEvent->Fill(1.1, smweight);
    
    // find all opposite sign same flavor di-leptons
    vector<DiLepton> dilepton;
    for(int k = 0; k < Ln; k++){
      TParticle& pk = lepton[index[k]];
      
      for(int j = k+1; j< Ln; j++){
	TParticle& pj = lepton[index[j]];
	
	bool OSSF = (pk.PID + pj.PID)==0;
	if ( ! OSSF ) continue;

	// we have oppositely charged same flavor leptons
	L1 = pk; 
	L2 = pj;

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
    h_nEvent->Fill(2.1, smweight);

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
    h_nEvent->Fill(3.1, smweight);

    // ----------------------------------------------------------
    // CUT 5: number of potential Z candidates > 1
    // ----------------------------------------------------------
    if( !( dilepton.size() > 1) )  continue;
    h_nEvent->Fill(4.1, smweight);    

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
    h_nEvent->Fill(5.1, smweight);
    
    // ----------------------------------------------------------
    // CUT 7: 12 < Z2mass < 120
    // ----------------------------------------------------------
    if( !(Z2.M() >  12) ) continue;
    if( !(Z2.M() < 120) ) continue;
    h_nEvent->Fill(6.1, smweight);    
    
    // We now have a di-Z event
    cout << "event " << entry << endl;
    cout << "\tleptons: " << lepton.size() << endl;
    for(size_t i=0; i < lepton.size(); i++)
      cout << "\t" << lepton[i] << endl;    
    std::cout << "\tZ1mass=" << Z1.M()
	      << "\tZ2mass=" << Z2.M()
	      << std::endl;
    
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
    h_nEvent->Fill(7.1, smweight);
    
    // ----------------------------------------------------------
    // CUT 9: LPT2 >10 GeV
    // ----------------------------------------------------------    
    if ( !(LPT2 > 10) ) continue;
    h_nEvent->Fill(8.1, smweight);
    
    // compute 4-lepton 4-vector
    H = Z1 + Z2;
    
    // ----------------------------------------------------------
    // CUT 10: m4l > 100 GeV
    // ----------------------------------------------------------    
    if ( !( H.M() > 100) ) continue;
    h_nEvent->Fill(9.1, smweight);
    
    Z1mass->Fill(Z1.M(), smweight);
    Z2mass->Fill(Z2.M(), smweight);
    Hmass->Fill(H.M(), smweight);
    HmassLarge->Fill(H.M(), smweight);

    MissingET* met = static_cast<MissingET*>(branchMET->At(0));    
    MET->Fill(met->MET, smweight);
    
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

     // get total weighted event count
     double summedWeight = h_nEvent->GetBinContent(1);

     cout << endl << "Summed weights: " << summedWeight << endl;
     h_nEvent->Scale(1.0/summedWeight);
     
     cuts->cd(); 
     h_nEvent->SetMinimum(0.0);
     h_nEvent->SetMaximum(1.1);
     
     h_nEvent->Draw();
     sprintf(filename, "plots/%s_PU%d_cutflow.png", prefix, pileup);
     cuts->Print(filename, "png");
   }
   
   theFile->cd();
   theFile->Write();
}
