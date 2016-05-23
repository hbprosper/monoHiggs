// ---------------------------------------------------------------------------
// mono Higgs analysis (pp -> 4 leptons + missing ET)
// Created: Les Houches 2015 Nic & HBP
// Updated: 31-Oct-2015 HBP - clean up, add ntuple output
// ---------------------------------------------------------------------------
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <map>

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
#include "HardScatter.h"
#include "LHParticle.h"

#include "monoHiggs.h"
#include "nic.h"
#include "Shrub.h"

using namespace std;

// ---------------------------------------------------------------------------
namespace {
  const double ZMASS=91.19;
  const int k4MU   = 1;
  const int k4E    = 2;
  const int k2E2MU = 3;
  const int kALL   = 4;
  
  const int ELECTRON = 11;
  const int MUON     = 13;
  //  const int GLUON    = 21;
  //  const int PHOTON   = 22;
  const int ZBOSON   = 23;
  const int HBOSON   = 25;
  int  DEBUG = 0;
};

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
struct DiLepton : public LHParticle
{
  DiLepton() : LHParticle() {}
  DiLepton(int PID_, LHParticle& L1_, LHParticle& L2_)
    : LHParticle(PID_, 0, 0, 0, 0),
      L1(L1_), L2(L2_)
  {
    TLorentzVector LL = L1 + L2;
    SetPtEtaPhiM(LL.Pt(), LL.Eta(), LL.Phi(), LL.M());
  }
  ~DiLepton() {}

  bool operator<(const LHParticle& p) const
  {
    return abs(this->M()-ZMASS) < abs(p.M()-ZMASS); 
  }
  
  LHParticle L1;
  LHParticle L2;
};

void copyLeptons(TClonesArray* electrons,
		 TClonesArray* muons,
		 std::vector<LHParticle>& lepton)
{
  if ( electrons )
    for(int i=0; i < electrons->GetEntriesFast(); i++)
      {
	Electron* p = static_cast<Electron*>(electrons->At(i));
	int PID = -ELECTRON * p->Charge;
	lepton.push_back(LHParticle(PID, p->PT, p->Eta, p->Phi));
      }

 if ( muons )
    for(int i=0; i < muons->GetEntriesFast(); i++)
      {
	Muon* p = static_cast<Muon*>(muons->At(i));
	int PID = -MUON * p->Charge;
	lepton.push_back(LHParticle(PID, p->PT, p->Eta, p->Phi));
      }
}

void purgeParticles(std::vector<LHParticle>& particles)
{
  // remove particles flagged as bad from vector
  int c = 0;
  for(size_t i = 0; i < particles.size(); i++)
    {
      if ( particles[i].Bad ) continue;
      particles[c] = particles[i];
      c++;
    }
  particles.resize(c);
}

void filterLeptons(std::vector<LHParticle>& lepton)
{
  // ----------------------------------------------------------
  // flag leptons that are too close together in
  // deltaR or that have pT > 5 GeV and lie outside the eta
  // acceptance.
  // ----------------------------------------------------------
  for(size_t i = 0; i < lepton.size(); i++)
    {
      LHParticle& p = lepton[i]; // get a reference not a copy
      
      // cuts differ for electrons and muons
      double PtCut;
      double EtaCut;
      if ( abs(p.PID) == MUON )
	{
	  PtCut  = 5.0;
	  EtaCut = 2.4;
	}
      else
	{
	  PtCut  = 7.0;
	  EtaCut = 2.5;
	}

      p.Bad = true;
      
      if ( !(p.Pt() > PtCut) ) continue;
      if ( !(abs(p.Eta()) < EtaCut) ) continue;
      
      p.Bad = false;
    }

  for(size_t i = 0; i < lepton.size(); i++)
    {
      LHParticle& pi = lepton[i];
      if ( pi.Bad ) continue;
      
      for(size_t j = i+1; j < lepton.size(); j++)
	{
	  LHParticle& pj = lepton[j];
	  if ( pj.Bad ) continue;

	  double dR = nic::deltaR(pi.Eta(), pi.Phi(),
				  pj.Eta(), pj.Phi());
	  if ( dR < 0.02 )
	    {
	      pi.Bad = true;
	      pj.Bad = true;
	    }
	}
    }
  // purge bad leptons from list of leptons
  purgeParticles(lepton);
}


void
listParticles(TClonesArray* particles)
{
  char record[256];
  sprintf(record,"%4s %8s %-10s %8s %8s %8s %4s %4s %4s %4s %4s",
	  "",   "PID", "name",
	  "PT", "Eta", "Phi",
	  "M1", "M2",
	  "D1", "D2", "Stat");
  cout << record << endl;
  for(int i = 0; i < particles->GetEntriesFast(); i++)
    {
      GenParticle* p = static_cast<GenParticle*>(particles->At(i));
      sprintf(record,"%4d %8d %-10s %8.2f %8.3f %8.3f %4d %4d %4d %4d %4d",
	      i, p->PID, LHParticle::name(p->PID).c_str(),
	      p->PT, p->Eta, p->Phi,
	      p->M1, p->M2,
	      p->D1, p->D2, p->Status);
      cout << record << endl;
    }
}


LHParticle* 
getZbosons(vector<LHParticle>&  particles,
	   vector<LHParticle*>& genZ,
	   vector<LHParticle*>& genL)
{
  LHParticle* H=0;
  for(size_t c = 0; c < particles.size(); c++)
    {
      if      ( particles[c].PID == HBOSON )
	{
	  H = &particles[c];
	}
	else if ( particles[c].PID == ZBOSON )
	  {
	    genZ.push_back( &particles[c] );
	
	    for(size_t i = 0; i < particles[c].Daughters.size(); i++)
	      {
		int d = particles[c].Daughters[i];
		int ID = abs(particles[d].PID);
		if( (ID == ELECTRON) || (ID == MUON) )
		  {
		    genL.push_back( &particles[d] );
		  }
	      }
	  }
      }
  //assert(genZ.size()>1);
  //assert(genL.size()>3);
  return H;
}


// match objects
void 
matchObjects(vector<LHParticle>&  r,
	     vector<LHParticle*>& p,
	     nic::Match& match,
	     float dRcut=0.05)
{
  for(size_t i = 0; i < r.size(); i++) r[i].ID  = -1;
  for(size_t i = 0; i < p.size(); i++) p[i]->ID = -1;

  for(size_t i = 0; i < r.size(); i++)
    for(size_t j = 0; j < p.size(); j++)
      match.add(i, r[i].Eta(),  r[i].Phi(),
		j, p[j]->Eta(), p[j]->Phi());
  match.run();

  for(size_t i = 0; i < match.order.size(); i++)
    {
      if ( !(match.order[i].first < dRcut) ) continue;
      int ii = match.order[i].second.first;
      int jj = match.order[i].second.second;
      r[ii].ID  = jj;
      p[jj]->ID = ii;
    }
}

void
findDileptons(vector<LHParticle>& lepton,
	      vector<DiLepton>& dilepton,
	      int finalState)
{
  vector<DiLepton> candidate;
  for(size_t k = 0; k < lepton.size(); k++)
    {
      LHParticle& pk = lepton[k]; // get a reference (an alias), not a copy
      for(size_t j = k+1; j < lepton.size(); j++)
	{
	  LHParticle& pj = lepton[j]; // get a reference, not a copy
	  // require charges to sum to zero
	  // and flavors to be the same
	  if ( !( (pk.PID + pj.PID)==0) ) continue;
	  
	  // we have oppositely charged same flavor leptons
	  LHParticle pkj = pk + pj;
	  
	  // QCD suppression
	  if ( !(pkj.M() > 10) ) continue;
	  
	  candidate.push_back(DiLepton(ZBOSON, pk, pj));
	}
    }
  if ( candidate.size() == 0 ) return;

  // we have one or more opposite sign same flavor dileptons
  
  // sort dilepton candidates in increasing delta(|mass-ZMASS|)
  // and take Z1 to be the first dilepton
  sort(candidate.begin(), candidate.end());    
  dilepton.push_back(candidate.front());
  
  // look for another Z candidate. if more than one, pick the
  // one with the highest pT among those that are further
  // than DR = 0.05 from Z1.
  //
  // But, if final state is 2e 2mu, make sure dileptons
  // differ in flavor
    
  double largestPT = -1;
  int index = -1;
  
  // get the unique identifiers of leptons that comprise Z1
  int L1UID = dilepton[0].L1.UID;
  int L2UID = dilepton[0].L2.UID;    
  
  // be sure to skip first dilepton, which is Z1    
  for(size_t k = 1; k < candidate.size(); k++)
    {
      // make sure current candidate does not share a lepton
      // with Z1
      if ( candidate[k].L1.UID == L1UID ) continue;
      if ( candidate[k].L1.UID == L2UID ) continue;
      
      if ( candidate[k].L2.UID == L1UID ) continue;
      if ( candidate[k].L2.UID == L2UID ) continue;	

      // current candidate does not share a lepton with Z1
      
      // For 2e2mu, the 2nd dilepton must differ in flavor from the first
      if ( finalState == k2E2MU )
	  {
	    bool sameFlavor =
	      abs(candidate[k].L1.PID) ==
	      abs(candidate[0].L1.PID);
	    if ( sameFlavor ) continue;
	  }

      // require dileptons to be further apart than DR = 0.05
      double dR = nic::deltaR(candidate[k].Eta(), candidate[k].Phi(),
			      candidate[0].Eta(), candidate[0].Phi());
      if ( !(dR > 0.05) ) continue;

      if ( candidate[k].Pt() > largestPT )
	{
	  largestPT = candidate[k].Pt();
	  index = k;
	}
    }
  // If index > -1, then we have found a 2nd dilepton.
  if ( index > -1 ) dilepton.push_back(candidate[index]);
}
// ---------------------------------------------------------------------------
// inputFile    input file (with .root extension) or a filelist
// finalstate   e.g.: "s_4mu" or "b_4mu"
// pileup       mean number of pileup events
// luminosity   integrated luminosity
// numberEvents obvious, no?!
// ---------------------------------------------------------------------------
void monoHiggs::hzz4l(string inputFile,
		      string sample,
		      int    numberEvents,
		      double luminosity,
		      double xsection,
		      int    pileup)
{
  cout << endl << "\t== monoHiggs::hzz4l ==" << endl;
  
  // get final state
  string namen = nic::nameonly(inputFile);
  int finalState = 0;
  if      ( sample.find("4mu")  != std::string::npos )
    {
      finalState = k4MU;
      namen += "_4mu";
    }
  else if ( sample.find("4e")   != std::string::npos )
    {
      finalState = k4E;
      namen += "_4me";
    }
  else if ( sample.find("2e2mu")!= std::string::npos )
    {
      finalState = k2E2MU;
      namen += "_2e2mu";
    }
  else if ( sample.find("all")  != std::string::npos )
    {
      finalState = kALL;
      namen += "_4l";
    }
  else
    nic::ciao("unrecognized final state (need 4mu, 4e, 2e2mu, or all)");

  // set background flag
  bool isBackground = sample.find("b") != std::string::npos;
  if ( isBackground )
    cout << endl << "\t=> background sample" << endl;
  else
    cout << endl << "\t=> signal sample" << endl;
  
  char finalstate[80];
  sprintf(finalstate, "%s", namen.c_str());
   
  // -----------------------------------------
  // analysis switches
  // -----------------------------------------
  bool original = true;
  bool masses   = true;
  bool useEventWeight=false;
  
  // -----------------------------------------
  // Objects
  //  muons  - PT > 5 GeV
  //           |Eta| < 2.4
  //           relisol < 0.4
  //
  //  muons  - PT > 7 GeV
  //           |Eta| < 2.7
  //           relisol < 0.4
  //
  //  Z1     - min({mZ - Zmass})
  //           mZ1 > 10 GeV
  //
  // -----------------------------------------
  // try to load libDelphes
  try
    {
      gSystem->Load("libDelphes");
    }
  catch (...)
    {
      nic::ciao("can't load libDelphes");
    }
  
  // make sure plots and histos directories exist
  nic::shell("mkdir -p plots; mkdir -p histos");
  
  // create empty root file for analysis results
  char filename[256];
  sprintf(filename, "histos/%s_PU%d.root", namen.c_str(), pileup);
  TFile* theFile = new TFile(filename, "RECREATE");
  if ( ! theFile->IsOpen() ) nic::ciao(string("can't create ") + filename);

  // create a chain of input files
  cout << endl << "\t=> input files:" << endl;
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
    std::cout<<"\t=> added input file "<<inputFiles[i]<<std::endl;
  }

  // -----------------------------------------
  // create A Delphes tree reader
  // -----------------------------------------
  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);


  // -----------------------------------------
  // initialize branch pointers
  // -----------------------------------------
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchJet   = treeReader->UseBranch("Jet");     
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchTower = treeReader->UseBranch("Tower");
  TClonesArray *branchRho   = treeReader->UseBranch("Rho");
  TClonesArray *branchMET   = treeReader->UseBranch("MissingET");  
  TClonesArray *branchElectron = 0;
  TClonesArray *branchMuon = 0;
  switch (finalState)
    {
    case k4MU:
      branchMuon = treeReader->UseBranch("Muon");
      break;
    case k4E:
      branchElectron = treeReader->UseBranch("Electron");
      break;
    case k2E2MU:
    case kALL:
    default:
      {
	branchElectron = treeReader->UseBranch("Electron");
	branchMuon = treeReader->UseBranch("Muon");
      }
      break;
    }


  // -----------------------------------------
  // create empty histograms
  // -----------------------------------------
  nic::setStyle();

  vector<TH1F*> h;
  
  cout << endl << "=> initialize histograms" << endl;
  TH1F* h_nEvent = new TH1F("h_nEvent", "Cut flow", 10, 0, 10);
  h.push_back(h_nEvent);
  h_nEvent->Sumw2();
  h_nEvent->GetXaxis()->SetBinLabel(1,"N(events)");
  h_nEvent->GetXaxis()->SetBinLabel(2,"N(lepton) > 2");
  h_nEvent->GetXaxis()->SetBinLabel(3,"N(isolepton) > 2");
  h_nEvent->GetXaxis()->SetBinLabel(4,"N(dilepton) > 0");
  
  TH1F* h_nleptons = new TH1F("nleptons", "", 10, 0, 10);
  h.push_back(h_nleptons);
  h_nleptons->GetXaxis()->SetTitle("#font[12]{n}_{leptons}");
  
  TH1F* h_njets = new TH1F("njets", "", 10, 0., 10.);
  h.push_back(h_njets);
  h_njets->GetXaxis()->SetTitle("#font[12]{n}_{jets}");

  TH1F* h_detisol = new TH1F("detisol", "", 100, -5.0, 5.0);
  h.push_back(h_detisol);
  h_detisol->GetXaxis()->SetTitle("#font[12]{l}_{isolation}");
  
  TH1F* h_rho = new TH1F("rho", "", 100, 0, 100);
  h.push_back(h_rho);
  h_rho->GetXaxis()->SetTitle("#font[12]{#rho}");

  int maxlep=4;
  TH1F* h_Eta[maxlep];
  TH1F* h_PT[maxlep];
  TH1F* h_genEta[maxlep];
  TH1F* h_genPT[maxlep];
  for(int c=0; c < maxlep; c++)
    {
      char name[20];
      
      sprintf(name, "Eta%d", c+1);
      h_Eta[c] = new TH1F(name, "", 100, -5., 5.);
      h.push_back(h_Eta[c]);
      h_Eta[c]->GetXaxis()->SetTitle("#font[12]{#eta}_{reco}");

      sprintf(name, "PT%d", c+1);
      h_PT[c] = new TH1F(name, "", 100, 0., 100.);
      h.push_back(h_PT[c]);
      h_PT[c]->GetXaxis()->SetTitle("#font[12]{p}_{T,reco} (GeV)");

      sprintf(name, "genEta%d", c+1);
      h_genEta[c] = new TH1F(name, "", 100, -5., 5.);
      h.push_back(h_genEta[c]);
      h_genEta[c]->GetXaxis()->SetTitle("#font[12]{#eta}_{gen}");

      sprintf(name, "genPT%d", c+1);
      h_genPT[c]  = new TH1F(name,  "", 100, 0., 100.);
      h.push_back(h_genPT[c]);
      h_genPT[c]->GetXaxis()->SetTitle("#font[12]{p}_{T,gen} (GeV)");
    }
  TH1F* h_Z1mass = new TH1F("Z1mass", "", 200,  0,  200.);
  h.push_back(h_Z1mass);
  h_Z1mass->GetXaxis()->SetTitle("#font[12]{m}_{Z1,reco} (GeV)");
  
  TH1F* h_Z2mass = new TH1F("Z2mass", "", 200,  0., 200.);
  h.push_back(h_Z2mass);
  h_Z2mass->GetXaxis()->SetTitle("#font[12]{m}_{Z2,reco} (GeV)");
  
  TH1F* h_genZ1mass = new TH1F("genZ1mass", "", 200, 0., 200.);
  h.push_back(h_genZ1mass);
  h_genZ1mass->GetXaxis()->SetTitle("#font[12]{m}_{Z1,gen} (GeV)");
    
  TH1F* h_genZ2mass = new TH1F("genZ2mass", "", 200, 0., 200.);
  h.push_back(h_Z2mass);
  h_genZ2mass->GetXaxis()->SetTitle("#font[12]{m}_{Z2,gen} (GeV)");
  
  TH1F* h_Hmass = new TH1F("Hmass", "", 200, 0., 200.);
  h.push_back(h_Hmass);
  h_Hmass->GetXaxis()->SetTitle("#font[12]{m}_{H,reco} (GeV)");
  
  TH1F* h_HmassLarge = new TH1F("HmassLarge", "", 200, 0., 400.);
  h.push_back(h_HmassLarge);
  h_HmassLarge->GetXaxis()->SetTitle("#font[12]{m}_{H,reco} (GeV)");
  
  TH1F* h_genHmass = new TH1F("genHmass", "", 200, 0., 200.);
  h.push_back(h_genHmass);
  h_genHmass->GetXaxis()->SetTitle("#font[12]{m}_{H,gen} (GeV)");
  
  TH1F* h_MET  = new TH1F("MET", "", 200, 0., 200.);
  h.push_back(h_MET);
  h_MET->GetXaxis()->SetTitle("missing #font[12]{E}_{T} (GeV)");
  
  TH1F* h_dZ1 = new TH1F("dZ1", "", 80,  0., 20.);
  h.push_back(h_dZ1);
  h_dZ1->GetXaxis()->SetTitle("#font[12]{m}_{Z1(gen)-Z1(reco)} (GeV)");
  
  TH1F* h_dZ2 = new TH1F("dZ2", "",  80,  0., 20.);
  h.push_back(h_dZ2);
  h_dZ2->GetXaxis()->SetTitle("#font[12]{m}_{Z2(gen)-Z2(reco)} (GeV)");

  TH1F* h_dRlepton = new TH1F("dRlepton", "", 100,  0., 0.5);
  h.push_back(h_dRlepton);
  h_dRlepton->GetXaxis()->SetTitle("#Delta#font[12]{R}"
				   "(#font[12]{l}, #font[12]{l_{gen}}");  

  TH1F* h_maxmatch = new TH1F("maxmatch", "", 10,  0., 10);
  h.push_back(h_maxmatch);
  h_maxmatch->GetXaxis()->SetTitle("max[match-index]");

  TH1F* h_nummatch = new TH1F("nummatch", "", 10,  0., 10);
  h.push_back(h_nummatch);
  h_nummatch->GetXaxis()->SetTitle("match multiplicity");

  TH1F* h_costhetastar = new TH1F("costhetastar", "", 100, -1, 1);
  h.push_back(h_costhetastar);
  h_costhetastar->GetXaxis()->SetTitle("cos(#Theta_{*})");

  TH1F* h_costheta1 = new TH1F("costheta1", "", 100, -1, 1);
  h.push_back(h_costheta1);
  h_costheta1->GetXaxis()->SetTitle("cos(#theta_{1})");

  TH1F* h_costheta2 = new TH1F("costheta2", "", 100, -1, 1);
  h.push_back(h_costheta2);
  h_costheta2->GetXaxis()->SetTitle("cos(#theta_{2})");

  TH1F* h_cosPhi = new TH1F("cosPhi", "", 100, -1, 1);
  h.push_back(h_cosPhi);
  h_cosPhi->GetXaxis()->SetTitle("cos(#Phi)");

  TH1F* h_cosPhi1 = new TH1F("cosPhi1", "", 100, -1, 1);
  h.push_back(h_cosPhi1);
  h_cosPhi1->GetXaxis()->SetTitle("cos(#Phi_{1})");

  for(size_t c=0; c < h.size(); c++)
    {
      h[c]->SetLineWidth(2);
      h[c]->SetLineColor(kBlack);
      h[c]->SetFillStyle(3001);
      h[c]->SetFillColor(kBlue);
    }
  
  // -----------------------------------------
  // EVENT LOOP
  // -----------------------------------------
  cout << endl << "\t=> begin event loop" << endl;
    
  long int numberOfEntries = treeReader->GetEntries();
  if ( numberEvents > 0 ) numberOfEntries = numberEvents;
  std::cout << "\t=> analyzing " << numberOfEntries <<" events" << std::endl;

  // -----------------------------------------
  // cross sections @ 13 TeV
  // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CrossSections#
  // Higgs_cross_sections_and_decay_b
  // BR  = 1.32e-4
  // -----------------------------------------
  double eventCount  = xsection * luminosity;
  double eventWeight = eventCount / numberOfEntries;
  
  char preamble[512];
  sprintf(preamble,
	  "========================================\n"
	  "Luminosity:       %10.2f events/fb\n"
	  "Cross section*BR: %10.2f fb\n"
	  "Event count:      %10.2f events\n"
	  "Event weight:     %10.2e fb/event\n"
	  "========================================\n",
	  luminosity, xsection, eventCount, eventWeight);
  cout << endl << preamble << endl;
  
  // -----------------------------------------
  // Set up output tree
  // -----------------------------------------
  string outfile = namen + string(".root");
  Shrub evtTree(outfile, "Analysis", namen);

  // -----------------------------------------
  // -----------------------------------------
  // ===> START
  // -----------------------------------------
  // -----------------------------------------
  bool debug_hardScatter = false;
  HardScatter hardScatter(branchParticle);
  int passed = 0;
  double totalPassed = 0.0;
  double totalPassedUnc = 0.0;
  
  int nd = 0; // number of diboson events
  
  //numberOfEntries = 2;
  
  for(Int_t entry = 0; entry < numberOfEntries; entry++){
    bool printMe = entry % 500 == 0;
    
    if ( printMe )
      std::cout << "\t=> processing event "
		<< entry
		<< "\t=> selected " << passed << std::endl;
    
    //Load branches for event
    treeReader->ReadEntry(entry);
    
    // get event weight
    double smweight = eventWeight;
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
    // copy data to standard Les Houches objects and
    // apply lepton filter
    // -----------------------------------------------------
    LHParticle::s_UID = 0; // reset unique ID (UID) 
    vector<LHParticle> lepton;
    copyLeptons(branchElectron, branchMuon, lepton);
   
    // apply lepton filter
    filterLeptons(lepton);
    h_nleptons->Fill(lepton.size());
    
    // sort in descending order of pT
    sort(lepton.begin(), lepton.end());
    for(size_t c=0; c < min(lepton.size(), (size_t)maxlep); c++)
      {
	h_PT[c]->Fill(lepton[c].Pt(), smweight);
	h_Eta[c]->Fill(lepton[c].Eta(), smweight);
      }

    // get gen particles
    vector<LHParticle> gparticles;
    hardScatter.get(gparticles, debug_hardScatter);
    
    vector<LHParticle*> genZ;
    vector<LHParticle*> genL;
    LHParticle* genH  = getZbosons(gparticles, genZ, genL);

    LHParticle* genZ1 = 0;
    LHParticle* genL1 = 0;
    LHParticle* genL2 = 0;
    if ( genZ.size() > 0 )
      {
	genZ1 = genZ[0];
	genL1 = genL[0];
	genL2 = genL[1];
      }

    LHParticle* genZ2 = 0;
    LHParticle* genL3 = 0;
    LHParticle* genL4 = 0;
    if ( genZ.size() > 1 )
      {
	genZ2 = genZ[1];
	genL3 = genL[2];
	genL4 = genL[3];    
      }
    vector<LHParticle> gL;
    maxlep = min(genL.size(), (size_t)maxlep);
    for(int c=0; c < maxlep; c++) gL.push_back(*genL[c]);
    sort(gL.begin(), gL.end());
    
    for(int c=0; c < maxlep; c++)
      {
	h_genPT[c]->Fill(gL[c].Pt(), smweight);
	h_genEta[c]->Fill(gL[c].Eta(), smweight);
      }

    if ( genZ1 )
      h_genZ1mass->Fill(genZ1->M(), smweight);
    if ( genZ2 )
      h_genZ2mass->Fill(genZ2->M(), smweight);
   
    if ( genH )
      h_genHmass->Fill(genH->M(), smweight);

    if ( printMe )
      if ( DEBUG > 1 )
	cout << "\t=> number of selected leptons: " << lepton.size() << endl;
    
    // ----------------------------------------------------------
    // CUT 1: number of leptons > 2
    // ----------------------------------------------------------
    if ( !(lepton.size() > 2) )  continue;
    h_nEvent->Fill(1.1, smweight);

    
    // get isolated leptons
    if ( branchRho )
      {
	double rho = 0.0;
	rho = static_cast<Rho*>(branchRho->At(0))->Rho;		
	h_rho->Fill(rho, smweight);
    
	if ( branchTrack && branchTower )
	  {
	    for(size_t k = 0; k < lepton.size(); k++)
	      {
		LHParticle& p = lepton[k]; // make a reference not a copy

		// for now use same isolation for muons and electrons
		int isoType = CMS_MU_ISOL;
		double detisol = nic::leptonIsolation(p.Pt(),
						      p.Eta(),
						      p.Phi(),
						      branchTrack,
						      branchTower,
						      branchRho,
						      isoType);
		h_detisol->Fill(detisol, smweight);

		// need to optimize cut
		if ( !(detisol < 2.0) )  continue;
      
		// non-isolated, so flag as bad
		p.Bad = true;
	      }
	  }
	// remove leptons flagged as bad
	purgeParticles(lepton);	
      }


    // Match gen leptons to reco leptons
    nic::Match match;
    if ( genL.size() > 0 )
      matchObjects(lepton, genL, match);

    // histogram dR
    for(size_t c=0; c < match.order.size(); c++)
      {
	float dR = match.order[c].first;
	h_dRlepton->Fill(dR);
      }
    
    // if no match, ID < 0
    int maxii=-1;
    int nmatch=0;
    for(size_t c=0; c < lepton.size(); c++)
      {
	if ( lepton[c].ID < 0 ) continue;
	maxii = c;
	nmatch++;
      }
    h_maxmatch->Fill(maxii);
    h_nummatch->Fill(nmatch);
    
    // ----------------------------------------------------------
    // CUT 2: number of isolated leptons > 2
    // ----------------------------------------------------------
    if ( !(lepton.size() > 2) )  continue;
    h_nEvent->Fill(2.1, smweight);
    
    // ----------------------------------------------------------
    // CUT 3: number of dileptons > 0
    // ----------------------------------------------------------
    // find all opposite sign same flavor di-leptons
    // that do not share leptons
    vector<DiLepton> dilepton;
    findDileptons(lepton, dilepton, finalState);
    
    if ( !(dilepton.size() > 0) )  continue;
    h_nEvent->Fill(3.1, smweight);

    // set the diBoson flag
    bool diBosonEvent = dilepton.size() > 1;
    
    // number of events that pass selection criteria
    passed++;
    totalPassed += eventWeight;
    totalPassedUnc += eventWeight*eventWeight;
    
    // ----------------------------------------------------------
    // END OF SELECTION
    // ----------------------------------------------------------
    
    if ( printMe )
      if ( DEBUG > 1 )
	{
	  cout << "\t=> number of dileptons: " << dilepton.size() << endl;
	  for(size_t c=0; c < dilepton.size(); c++)
	    cout << "\t" << dilepton[c] << endl;
	}
    
    LHParticle Z1 = dilepton[0];
    LHParticle L1 = dilepton[0].L1;
    LHParticle L2 = dilepton[0].L2;
    assert(L1.PID+L2.PID==0);
    int numberLeptons = 2;
    
    h_Z1mass->Fill(Z1.M(), smweight);

    LHParticle Z2;
    LHParticle L3;
    LHParticle L4;
    LHParticle H;
    if ( diBosonEvent )
      {
	Z2 = dilepton[1];
	L3 = dilepton[1].L1;
	L4 = dilepton[1].L2;
	assert(L3.PID+L4.PID==0);
	numberLeptons += 2;

	h_Z2mass->Fill(Z2.M(), smweight);

	// compute 4-lepton 4-vector
	H = Z1 + Z2;
 
	h_Hmass->Fill(H.M(), smweight);
	h_HmassLarge->Fill(H.M(), smweight);	
      }
    else
      {
	// make sure L3 differs from L1 and L2
	int UID1 = L1.UID;
	int UID2 = L2.UID;
	for(size_t k = 0; k < lepton.size(); k++)
	  {
	    if ( lepton[k].UID == UID1 ) continue;
	    if ( lepton[k].UID == UID2 ) continue;
	    L3 = lepton[k];
	    numberLeptons += 1;
	    break;
	  }
      }	

    // if ( printMe )
    //   {
    // 	cout << endl << "====> Entry: " << entry << endl;
    // 	cout << L1 << endl;
    // 	cout << L2 << endl;
    // 	cout << L3 << endl;
    // 	cout << L4 << endl;
    //   }
    
    // Missing transverse energy
    MissingET* met = static_cast<MissingET*>(branchMET->At(0));    
    h_MET->Fill(met->MET, smweight);
    
    // store jet info
    vector<LHParticle> jet;
    for(int i = 0; i < branchJet->GetEntriesFast(); i++)
      {
	Jet* j = static_cast<Jet*>(branchJet->At(i));    
	if ( ! (j->PT > 20) ) continue; // Run I cut is pT > 30 GeV
	if ( ! (abs(j->Eta) < 4.7) ) continue;
	jet.push_back(LHParticle(81,    // ID for jet
				 j->PT,
				 j->Eta,
				 j->Phi,
				 j->Mass));
      }
    sort(jet.begin(), jet.end());
  
    h_njets->Fill(jet.size(), smweight);

    if ( genZ1 )
      {
	// compute (Z1 - genZ1).M()
	// this will be zero for a perfect match
	// between reco Z1 and gen Z1
	TLorentzVector dZ1 = Z1 - *genZ1;
	h_dZ1->Fill(dZ1.M(), smweight);
      }
    
    if ( diBosonEvent && genZ2 )
      {
	// compute (Z2 - genZ2).M()
	TLorentzVector dZ2 = Z2 - *genZ2;
	h_dZ2->Fill(dZ2.M(), smweight);
      }

    if ( DEBUG > 0 )
      {
	cout << endl;
	cout << "Event summary " << entry << endl;
	if ( genH )
	  cout << " gen(H):   " << *genH  << endl;
	if ( genZ1 )
	  cout << " gen(Z1):  " << *genZ1 << endl;
	if ( genZ2 )
	  cout << " gen(Z2):  " << *genZ2 << endl;
	
	cout << " reco(Z1): " << Z1 << endl;
	if ( diBosonEvent )
	  cout << " reco(Z2): " << Z2 << endl;	
      }
    
    // ----------------------------------------------------------
    // write out some event quantities
    // ----------------------------------------------------------
    
    evtTree.Clear();

    evtTree.weight    = smweight;

    evtTree.njets     = (int)jet.size();
    evtTree.nleps     = numberLeptons;

    if ( genZ1 )
      {
	evtTree.genZ1pt   = genZ1->Pt();
	evtTree.genZ1eta  = genZ1->Eta();
	evtTree.genZ1phi  = genZ1->Phi();
	evtTree.genZ1mass = genZ1->M();
      }
    if ( genZ2 )
      {
	evtTree.genZ2pt   = genZ2->Pt();
	evtTree.genZ2eta  = genZ2->Eta();
	evtTree.genZ2phi  = genZ2->Phi();
	evtTree.genZ2mass = genZ2->M();
      }
    if ( genL1 )
      {
	evtTree.genl1match= genL1->ID;
	evtTree.genl1PID  = genL1->PID;
	evtTree.genl1pt   = genL1->Pt();
	evtTree.genl1eta  = genL1->Eta();
	evtTree.genl1phi  = genL1->Phi();
      }
    if ( genL2 )
      {
	evtTree.genl2match= genL2->ID;
	evtTree.genl2PID  = genL2->PID;    
	evtTree.genl2pt   = genL2->Pt();
	evtTree.genl2eta  = genL2->Eta();
	evtTree.genl2phi  = genL2->Phi();
      }
    if ( genL3 )
      {
	evtTree.genl3match= genL3->ID;
	evtTree.genl3PID  = genL3->PID;    
	evtTree.genl3pt   = genL3->Pt();
	evtTree.genl3eta  = genL3->Eta();
	evtTree.genl3phi  = genL3->Phi();
      }
    if ( genL4 )
      {
	evtTree.genl4match= genL4->ID;
	evtTree.genl4PID  = genL4->PID;    
	evtTree.genl4pt   = genL4->Pt();
	evtTree.genl4eta  = genL4->Eta();
	evtTree.genl4phi  = genL4->Phi();    
      }
    evtTree.Z1pt   = Z1.Pt();
    evtTree.Z1eta  = Z1.Eta();
    evtTree.Z1phi  = Z1.Phi();
    evtTree.Z1mass = Z1.M();

    evtTree.Z2pt   = Z2.Pt();
    evtTree.Z2eta  = Z2.Eta();
    evtTree.Z2phi  = Z2.Phi();
    evtTree.Z2mass = Z2.M();

    evtTree.Hpt    = H.Pt();
    evtTree.Heta   = H.Eta();
    evtTree.Hphi   = H.Phi();
    evtTree.Hmass  = H.M();    

    evtTree.l1match= L1.ID;
    evtTree.l1PID  = L1.PID;
    evtTree.l1pt   = L1.Pt();
    evtTree.l1eta  = L1.Eta();
    evtTree.l1phi  = L1.Phi();

    evtTree.l2match= L2.ID;
    evtTree.l2PID  = L2.PID;    
    evtTree.l2pt   = L2.Pt();
    evtTree.l2eta  = L2.Eta();
    evtTree.l2phi  = L2.Phi();

    evtTree.l3match= L3.ID;
    evtTree.l3PID  = L3.PID;    
    evtTree.l3pt   = L3.Pt();
    evtTree.l3eta  = L3.Eta();
    evtTree.l3phi  = L3.Phi();

    evtTree.l4match= L4.ID;
    evtTree.l4PID  = L4.PID;    
    evtTree.l4pt   = L4.Pt();
    evtTree.l4eta  = L4.Eta();
    evtTree.l4phi  = L4.Phi();

    if ( jet.size() > 0 )
      {
    	evtTree.j1pt   = jet[0].Pt();
    	evtTree.j1eta  = jet[0].Eta();
    	evtTree.j1phi  = jet[0].Phi();
    	evtTree.j1mass = jet[0].M();
      }

    if ( jet.size() > 1 )
      {
    	evtTree.j2pt   = jet[1].Pt();
    	evtTree.j2eta  = jet[1].Eta();
    	evtTree.j2phi  = jet[1].Phi();
    	evtTree.j2mass = jet[1].M();

    	LHParticle dijet = jet[0] + jet[1];
    	evtTree.massjj   = dijet.M();
    	evtTree.deltaetajj = abs(jet[0].Eta() - jet[1].Eta());
      }

    evtTree.met = met->MET;
    
    // add MELA variables
    float costhetastar=-2;
    float costheta1=-2;
    float costheta2=-2;
    float cosPhi=-2;
    float cosPhi1=-2;

  
    if ( diBosonEvent )
      {	
    	nic::computeMELAangles(L1, L1.PID,
    			       L2, L2.PID,
    			       L3, L3.PID,
    			       L4, L4.PID,
    			       costhetastar,
    			       costheta1,
    			       costheta2,
    			       cosPhi,
    			       cosPhi1);
      }
    
    evtTree.costhetastar = costhetastar;
    evtTree.costheta1    = costheta1;
    evtTree.costheta2    = costheta2;
    evtTree.cosPhi       = cosPhi;
    evtTree.cosPhi1      = cosPhi1;

    if ( costhetastar > -2 )
      {
    	h_costhetastar->Fill(costhetastar, smweight);
    	h_costheta1->Fill(costheta1, smweight);
    	h_costheta2->Fill(costheta2, smweight);
    	h_cosPhi->Fill(cosPhi, smweight);
    	h_cosPhi1->Fill(cosPhi1, smweight);
      }

    if ( DEBUG > 0 )
      if ( diBosonEvent )
    	{
    	  char record[1024];
    	  nd++;
    	  sprintf(record,
    		  "\ndiboson event: %5d\tnleps = %3d\n"
    		  "L1.PID = %3d\tL2.PID = %3d\n"
    		  "L3.PID = %3d\tL4.PID = %3d",
    		  nd, evtTree.nleps,
    		  evtTree.l1PID, evtTree.l2PID,
    		  evtTree.l3PID, evtTree.l4PID);
    	  cout << record << endl;
    	}
    
    evtTree.Fill();
    
   } // END OF EVENT LOOP

  evtTree.Close();
  
   // -----------------------------------------------------------
   // Finish up
   // -----------------------------------------------------------
   if ( original )
     {	
       TCanvas* cPT = new TCanvas("PT", "PT", 10, 10, 800, 800);
       cPT->Divide(2, 2);

       TCanvas* cgenPT = new TCanvas("genPT", "genPT", 100, 100, 800, 800);
       cgenPT->Divide(2, 2);       

       for(int c=0; c < 4; c++)
	 {
	   cPT->cd(c+1);
	   h_PT[c]->Draw("hist");

	   cgenPT->cd(c+1);
	   h_genPT[c]->Draw("hist");
	 }
       cPT->Update();
       cgenPT->Update();
       
       theFile->cd();
       cPT->Write();
       cgenPT->Write();
     }
      
   if ( masses )
     {
       TCanvas* cmass = new TCanvas("masses","masses", 200, 200, 800, 800);
       cmass->Divide(2, 2);
     
       cmass->cd(1); 
       h_Z1mass->Draw("hist"); 
     
       cmass->cd(2); 
       h_Z2mass->Draw("hist"); 
       
       cmass->cd(3); 
       h_Hmass->Draw("hist");

       cmass->cd(4); 
       h_HmassLarge->Draw("hist");
       
       cmass->Update();
       
       theFile->cd();
       cmass->Write();
     }
   
   // get total weighted event count
   double summedWeight    = h_nEvent->GetBinContent(1);
   double summedWeightUnc = h_nEvent->GetBinError(1);

   char summary[512];
   totalPassedUnc = sqrt(totalPassedUnc);
   sprintf(summary,
	   "======================================================\n"
	   "Event count:      %10.2f +/- %-10.2f events\n"
	   "Selected count:   %10.2f +/- %-10.2f events\n"
	   "======================================================\n",
	   summedWeight, summedWeightUnc,
	   totalPassed, totalPassedUnc);
   
   cout << summary << endl;   

   h_nEvent->Scale(1.0/summedWeight);
   h_nEvent->SetMaximum(1.e-3);
   h_nEvent->SetMaximum(1.1);

   TCanvas* cuts = new TCanvas("cutflow", "cut flow", 500, 10, 500, 500);
   cuts->cd();
   cuts->SetLogy();     
   h_nEvent->Draw("hist");
   cuts->Update();
   sprintf(filename, "plots/%s_PU%d_cutflow.png", namen.c_str(), pileup);
   cuts->Print(filename, "png");
   
   theFile->cd();
   theFile->Write();
   cuts->Write();
   
   gSystem->Sleep(5000);
}
