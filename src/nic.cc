// ---------------------------------------------------------------------------
// A few simple self-contained utilities for monoHiggs analysis.
// A few simple self-contained string utilities
// Created: 16-Jun-2015 HBP & NDP   Les Houches
// ---------------------------------------------------------------------------
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "nic.h"

using namespace std;

void nic::computeMELAangles(TLorentzVector L1_Z1, int L1_Z1_PID,
			    TLorentzVector L2_Z1, int L2_Z1_PID,
			    TLorentzVector L1_Z2, int L1_Z2_PID,
			    TLorentzVector L2_Z2, int L2_Z2_PID,
			    float& costhetastar, 
			    float& costheta1, 
			    float& costheta2, 
			    float& cosPhi, 
			    float& cosPhi1)
{
  // make sure we have opposite sign leptons
  if ( L1_Z1_PID * L2_Z1_PID > 0 ) return;
  if ( L1_Z2_PID * L2_Z2_PID > 0 ) return;
  
  // sort Z1 leptons so that L1_Z1 is the negative lepton
  // Note: positive PDG ID corresponds to negative leptons
  if ( L2_Z1_PID > 0 ) swap(L1_Z1, L2_Z1);

  // sort Z2 leptons so that negative lepton is first
  if ( L2_Z2_PID < 0 ) swap(L1_Z2, L2_Z2);
  
  // build Z and H 4-vectors
  TLorentzVector Z1 = L1_Z1 + L2_Z1;
  TLorentzVector Z2 = L1_Z2 + L2_Z2;
  TLorentzVector H  = Z1 + Z2;

  // get boosts to H, Z1, and Z2 rest frames
  TVector3 boostH  =  -H.BoostVector();
  TVector3 boostZ1 = -Z1.BoostVector();
  TVector3 boostZ2 = -Z2.BoostVector();
  
  //------ costhetastar (polar angle of Z1 in H frame)
  // boost Z1 to H frame
  TLorentzVector Z1inH(Z1); // make copy of Z1
  Z1inH.Boost(boostH);      // boost Z1 to H frame
  costhetastar = (Z1inH.Vect()).CosTheta();
  if ( isnan(costhetastar) )
    {
      cout << "\t*** costhetastar is NAN" << endl;
      costhetastar = 2;
    }
  
  //------ costheta1 (cosine of polar angle of negative lepton in Z1 frame
  // where "z" axis is taken to be the negative Z2 direction in the Z1 frame.
  
  // boost L1 to Z1 frame
  TLorentzVector L1inZ1(L1_Z1);
  L1inZ1.Boost(boostZ1);
  TVector3 UL1inZ1 = (L1inZ1.Vect()).Unit();
  
  // boost Z2 to Z1 frame
  TLorentzVector Z2inZ1(Z2);
  Z2inZ1.Boost(boostZ1);
  TVector3 UZ2inZ1 = (Z2inZ1.Vect()).Unit();

  // dot product between -Unit(Z2) and Unit(L1)
  costheta1 = -UZ2inZ1.Dot(UL1inZ1);
  if ( isnan(costheta1) )
    {
      cout << "\t*** costheta1 is NAN" << endl;
      costheta1 = 2;
    }
  
  //------ costheta2 (cosine of polar angle of negative lepton in Z2 frame
  // where "z" axis is taken to be the negative Z1 direction in the Z2 frame.

  // boost L1 to Z2 frame
  TLorentzVector L1inZ2(L1_Z2);
  L1inZ2.Boost(boostZ2);
  TVector3 UL1inZ2 = (L1inZ2.Vect()).Unit();
  
  // boost Z1 to Z2 frame
  TLorentzVector Z1inZ2(Z1);  
  Z1inZ2.Boost(boostZ2);
  TVector3 UZ1inZ2 = (Z1inZ2.Vect()).Unit();

  costheta2 = -UZ1inZ2.Dot(UL1inZ2);
  if ( isnan(costheta2) )
    {
      cout << "\t*** costheta2 is NAN" << endl;
      costheta2 = 2;
    }
  
  //------ cos(Phi) (cosine of angle between lepton decay planes in H frame)
  TLorentzVector L1_Z1inH(L1_Z1);
  TLorentzVector L2_Z1inH(L2_Z1);
  TLorentzVector L1_Z2inH(L1_Z2);
  TLorentzVector L2_Z2inH(L2_Z2);

  // boost leptons to H frame
  L1_Z1inH.Boost(boostH);
  L2_Z1inH.Boost(boostH);
  L1_Z2inH.Boost(boostH);
  L2_Z2inH.Boost(boostH);

  // Z1: U1 = L1 x L2/|L1 x L2| in H frame
  TVector3 U1 = (L1_Z1inH.Vect().Cross(L2_Z1inH.Vect())).Unit();

  // Z2: U2 = L1 x L2/|L1 x L2| in H frame
  TVector3 U2 = (L1_Z2inH.Vect().Cross(L2_Z2inH.Vect())).Unit();
  
  cosPhi = U1.Dot(U2);
  if ( isnan(cosPhi) )
    {
      cout << "\t*** cosPhi is NAN" << endl;
      cosPhi = 2;
    }

  //------ cosPhi1 (cosine of angle between (Z1,z-axis) plane and (L1,L2) plane
  // of Z1 leptons, in H frame)
  TVector3 beamAxis(0,0,1);
  TVector3 UZ1inH = (Z1inH.Vect()).Unit();
  TVector3 UbeamZ1inH = beamAxis.Cross(UZ1inH);
  cosPhi1 = U1.Dot(UbeamZ1inH);    
  if ( isnan(cosPhi1) )
    {
      cout << "\t*** cosPhi1 is NAN" << endl;
      cosPhi1 = 2;
    }
}


double nic::leptonIsolation(double PT, double Eta, double Phi,
			    TClonesArray* tracks,
			    TClonesArray* towers,
			    TClonesArray* rhos,
			    int isoType)
{
  if ( tracks == 0 ) return 0;
  if ( towers == 0 ) return 0;
  if ( rhos   == 0 ) return 0;

  // isoType
  // 1         CMS muon isolation
  // 2         CMS electron isolation
  // 3         CMS photon isolation
  //
  // 11        ATLAS muon isolation
  // 12        ATLAS electron isolation
  // 13        ATLAS photon isolation

  double isolation = 0.0;
  switch (isoType)
    {
    case 1:
    case 2:
      {
	double dRcut=0.3;
	double dRMax=0.5;
	
	int ntracks = tracks->GetEntriesFast();
	double sumTrackPt = 0.0;
	for(int i=0; i < ntracks; i++)
	  {
	    Track* track = static_cast<Track*>(tracks->At(i));
	    double dR = nic::deltaR(Eta, Phi,
				    track->Eta, track->Phi);
	    if ( !(dR < dRcut) ) continue;
	    sumTrackPt += track->PT;
	  }
	
	int ntowers = towers->GetEntriesFast();
	double sumTowerPt = 0.0;
	for(int i=0; i < ntowers; i++)
	  {
	    Tower* tower = static_cast<Tower*>(towers->At(i));
	    double dR = nic::deltaR(Eta, Phi,
				    tower->Eta, tower->Phi);
	    if ( !(dR < dRcut) ) continue;
	    sumTowerPt += tower->ET;
	  }
	isolation = sumTrackPt + sumTowerPt;
	
	// correct for pile-up by subtracting the event-by-event
	// average transverse momentum due to pileup.
	// use the first rho value, which is for the eta range [-2.5, 2.5]
	double rho = static_cast<Rho*>(rhos->At(0))->Rho;
	double pileupOffset = rho * dRMax*dRMax*M_PI;
	isolation -= pileupOffset;

	// now "normalize" 
	isolation /= PT;
      }
      break;
    default:
      isolation = 0.0;
    } 
  return isolation;
}

void nic::ciao(string message)
{
  cout << "** ciao! " << message << endl;
  exit(0);
}
std::string nic::strip(std::string line)
{
  int l = line.size();
  if ( l == 0 ) return std::string("");
  int n = 0;
  while (((line[n] == 0)    ||
	  (line[n] == ' ' ) ||
	  (line[n] == '\n') ||
	  (line[n] == '\t')) && n < l) n++;
  
  int m = l-1;
  while (((line[m] == 0)    ||
	  (line[m] == ' ')  ||
	  (line[m] == '\n') ||
	  (line[m] == '\t')) && m > 0) m--;
  return line.substr(n,m-n+1);
}

std::vector<std::string> nic::split(std::string str)
{
  vector<string> vstr;
  std::istringstream stream(str);
  while ( stream )
    {
      std::string istr;
      stream >> istr;
      if ( stream ) vstr.push_back(istr);
    }
  return vstr;
}

std::string nic::replace(std::string str,
			 std::string oldstr,
			 std::string newstr)
{
  return std::string(TString(str).ReplaceAll(oldstr, newstr).Data());
}

std::string nic::nameonly(std::string filename)
{
  int i = filename.rfind("/");
  int j = filename.rfind(".");
  if ( j < 0 ) j = filename.size();
  return filename.substr(i+1,j-i-1);
}

std::string nic::shell(std::string cmd)
{
  FILE* f = popen(cmd.c_str(),"r");
  int buffsize=8192;
  char s[8192];
  int n = fread(s,1,buffsize,f);
  pclose(f);
  std::string result = strip(std::string(s).substr(0,n));
  return result;
}

///
double nic::deltaPhi(double phi1, double phi2)
{
  double deltaphi = phi2 - phi1;
  if ( fabs(deltaphi) > M_PI ) deltaphi = 2 * M_PI - fabs(deltaphi);
  return deltaphi;
}

double nic::deltaR(double eta1, double phi1, double eta2, double phi2)
{
  double deltaeta = eta1 - eta2;
  double deltaphi = deltaPhi(phi1, phi2);
  return sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
}

void nic::Match::add(int ii, float eta1, float phi1,
		     int jj, float eta2, float phi2)
{
  double dR = nic::deltaR(eta1, phi1, eta2, phi2);
  order.push_back(pair<double, std::pair<int, int> >());
  order.back().first  = dR;
  order.back().second.first = ii;
  order.back().second.second = jj;
}
void nic::Match::run()
{
  sort(order.begin(), order.end());
}

void nic::setStyle()
{
  int TEXTFONT=42;
  int NDIVX=505;
  
  gStyle->SetPalette(1);

  // For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(500);   //Height of canvas
  gStyle->SetCanvasDefW(500);   //Width of canvas
  gStyle->SetCanvasDefX(0);     //Position on screen
  gStyle->SetCanvasDefY(0);

  // For the Pad:
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(kFALSE);
  gStyle->SetPadGridY(kFALSE);
  gStyle->SetGridColor(kGreen);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
	
  // For the frame:
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);
	
  // For the histo:
  gStyle->SetHistLineColor(1);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(2);
	
  gStyle->SetEndErrorSize(2);
  //gStyle->SetErrorX(0.);
	
  gStyle->SetMarkerSize(0.5);
  gStyle->SetMarkerStyle(20);

  // For the fit/function:
  gStyle->SetOptFit(1);
  gStyle->SetFitFormat("5.4g");
  gStyle->SetFuncColor(2);
  gStyle->SetFuncStyle(1);
  gStyle->SetFuncWidth(1);
	
  // For the date:
  gStyle->SetOptDate(0);

  // For the statistics box:
  gStyle->SetOptFile(0);
  //gStyle->SetOptStat("");
  
  // To display the mean and RMS:
  gStyle->SetOptStat("mr"); 
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(TEXTFONT);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.3);
    
  // Margins:
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.20);
  gStyle->SetPadRightMargin(0.10);
	
  // For the Global title:
  gStyle->SetOptTitle(0); 
  gStyle->SetTitleFont(TEXTFONT);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.05);

  // For the axis titles:
  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(TEXTFONT, "XYZ");
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetTitleXOffset(1.25);      //(1.25);
  gStyle->SetTitleYOffset(1.40);      //(1.60);
					
  // For the axis labels:
  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(TEXTFONT, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");
  
  // For the axis:
  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(NDIVX, "X");
  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickX(1);  
  gStyle->SetPadTickY(1);

  // Change for log plots:
  gStyle->SetOptLogx(0);
  gStyle->SetOptLogy(0);
  gStyle->SetOptLogz(0);

  // Postscript options:
  gStyle->SetPaperSize(20.,20.);
  gStyle->cd();
}
