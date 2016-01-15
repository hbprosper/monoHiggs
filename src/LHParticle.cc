// ---------------------------------------------------------------------------
// File: LHParticle.cc
// Description: prototype of a generic particle class for use in Les Houches
//              analysis description.
// created: Les Houches 2015 HBP
// ---------------------------------------------------------------------------
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <map>

#include "nic.h"
#include "LHParticle.h"

using namespace std;

int LHParticle::s_UID=0;

LHParticle::LHParticle()
  : TLorentzVector(),
    ID(-1),
    PID(0),
    Status(0),
    Mother(0),
    Bad(false),
    Name(""),
    Daughters(std::vector<int>()),
    Value(std::map<std::string, double>()),
    UID(s_UID++)
{}

LHParticle::LHParticle(int PID_, 
		       double PT, double Eta, double Phi, double Mass)
  : TLorentzVector(),
    ID(-1),
    PID(PID_),
    Status(0),
    Mother(0),
    Bad(false),
    Name(nic::particleName(PID_)),
    Daughters(std::vector<int>()),
    Value(std::map<std::string, double>()),
    UID(s_UID++)
{
  SetPtEtaPhiM(PT, Eta, Phi, Mass);
}


LHParticle::LHParticle(const LHParticle& p)
  : TLorentzVector(),
    ID(p.ID),
    PID(p.PID),
    Status(p.Status),
    Mother(p.Mother),
    Bad(p.Bad),
    Name(p.Name),
    Daughters(p.Daughters),
    Value(p.Value),
    UID(p.UID)
{
  SetPtEtaPhiM(p.Pt(), p.Eta(), p.Phi(), p.M());
}


LHParticle::~LHParticle() {}

bool LHParticle::operator<(const LHParticle& p) const
{
  return p.Pt() < this->Pt();
}

LHParticle& LHParticle::operator=(const LHParticle& rhs)
{
  if ( this != &rhs )
    {
      LHParticle p(rhs); // call copy constructor
      SetPtEtaPhiM(p.Pt(), p.Eta(), p.Phi(), p.M());
      ID     = p.ID;
      PID    = p.PID;
      Status = p.Status;
      Mother = p.Mother;
      Bad    = p.Bad;
      Name   = p.Name;
      Daughters = p.Daughters;
      Value  = p.Value;
      UID    = p.UID;
    }
  return *this;
}

LHParticle LHParticle::operator+(const LHParticle& o) const
{
  const TLorentzVector* p1 = dynamic_cast<const TLorentzVector*>(this);
  const TLorentzVector* p2 = dynamic_cast<const TLorentzVector*>(&o);
  TLorentzVector  p  = *p1 + *p2;
  return LHParticle(PID, p.Pt(), p.Eta(), p.Phi(), p.M());
}

LHParticle LHParticle::operator-(const LHParticle& o) const
{
  const TLorentzVector* p1 = dynamic_cast<const TLorentzVector*>(this);
  const TLorentzVector* p2 = dynamic_cast<const TLorentzVector*>(&o);
  TLorentzVector  p  = *p1 - *p2;
  return LHParticle(PID, p.Pt(), p.Eta(), p.Phi(), p.M());
}

LHParticle LHParticle::operator*(double a) const
{
  TLorentzVector p = a * (*dynamic_cast<const TLorentzVector*>(this));
  return LHParticle(PID, p.Pt(), p.Eta(), p.Phi(), p.M());
}

std::ostream& operator<<(std::ostream& os, const LHParticle& o)
{
  string str;
  char rec[8];
  sprintf(rec, " %4d", o.Mother);
  str += string(rec);
  for(size_t c=0; c < o.Daughters.size(); c++)
    {
      sprintf(rec, " %4d", o.Daughters[c]);
      str += string(rec);
    }

  char record[80];
  sprintf(record, "%4d %4d %-10s %8.2f %8.3f %8.3f %8.3f %s",
	  o.UID, o.ID, o.Name.c_str(), o.Pt(), o.Eta(), o.Phi(), o.M(),
	  str.c_str());
  os << record;
  return os;
}
