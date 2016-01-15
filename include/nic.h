#ifndef NIC_H
#define NIC_H
// ---------------------------------------------------------------------------
// A few simple self-contained utilities for monoHiggs analysis.
// A few simple self-contained string utilities
// Created: 16-Jun-2015 HBP & NDP   Les Houches
// ---------------------------------------------------------------------------
#include <string>
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"

///
struct nic
{
  /// Compute decay angles from the lepton four-vectors and lepton PDG IDs.
  static void computeMELAangles(TLorentzVector L1_Z1, int L1_Z1_PID,
				TLorentzVector L2_Z1, int L2_Z1_PID,
				TLorentzVector L1_Z2, int L1_Z2_PID,
				TLorentzVector L2_Z2, int L2_Z2_PID,
				float& costhetastar, 
				float& costheta1, 
				float& costheta2, 
				float& cosPhi, 
				float& cosPhi1);
  ///
  static void setStyle();
  ///
  static void ciao(std::string message);
  ///
  static std::string strip(std::string line);
  ///
  static std::vector<std::string> split(std::string str);
  ///
  static std::string replace(std::string str,
			     std::string oldstr,
			     std::string newstr);
  ///
  static std::string nameonly(std::string filename);
  ///
  static std::string shell(std::string cmd);
  /// Return particle name given PDG id
  static std::string particleName(int pdgid);
  ///
  static double deltaPhi(double phi1, double phi2);
  ///
  static double deltaR(double eta1, double phi1, double eta2, double phi2);

  struct Match
  {
    Match() {}
    ~Match() {}
    void add(int ii, float eta1, float phi1,
	     int jj, float eta2, float phi2);
    void run();
    std::vector<std::pair<float, std::pair<int, int> > > order;
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
  static double leptonIsolation(double pt, double eta, double phi,
				TClonesArray* tracks,
				TClonesArray* towers,
				TClonesArray* rho,
				int isoType=1);
};

#endif
