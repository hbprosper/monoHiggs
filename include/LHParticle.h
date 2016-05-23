#ifndef LHPARTICLE_H
#define LHPARTICLE_H
// ---------------------------------------------------------------------------
// File: LHParticle.h
// Description: prototype of a generic Les Houches particle class for use in
//              Les Houches analysis description.
// created: Les Houches 2015 HBP
// ---------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <map>
#include "TLorentzVector.h"
// ---------------------------------------------------------------------------
struct LHParticle : public TLorentzVector
{
  LHParticle();
  LHParticle(int PID_,
	     double PT, double Eta, double Phi, double Mass=0);
  LHParticle(const LHParticle& p);
  ~LHParticle();

  static std::string name(int pdgid);
  
  bool        operator<(const LHParticle& p) const;
  LHParticle& operator=(const LHParticle& p);
  LHParticle  operator+(const LHParticle& o) const;
  LHParticle  operator-(const LHParticle& o) const;
  LHParticle  operator*(double a) const;

  int  ID;    
  int  PID;   // PDG ID
  int  Status;  
  int  Mother;
  bool Bad;
  std::string Name;
  std::vector<int> Daughters;
  std::map<std::string, double> Value;
  int UID; // event unique identifier
  static int s_UID;
};
std::ostream& operator<<(std::ostream& os, const LHParticle& o);

#endif
