#ifndef HARDSCATTER_H
#define HARDSCATTER_H
// ---------------------------------------------------------------------------
// File: HardScatter.j
// Description: extract the hard scatter events from the Delphes record.
// created: Les Houches 2015 Nic & HBP
// ---------------------------------------------------------------------------
#include <vector>
#include <string>

#include "classes/DelphesClasses.h"
#include "TClonesArray.h"
#include "nic.h"
#include "LHParticle.h"

class HardScatter
{
 public:
  HardScatter(TClonesArray* particles_);
  ~HardScatter();
  
  void get(std::vector<LHParticle>& plist, bool debug_=false);
  
 private:
  TClonesArray* particles;
  bool debug;

  void printAt(int ID);
  
  void findMother(GenParticle* mother, int& pos, int depth=0);
  
  void climbTree(std::vector<LHParticle>& plist,
		 GenParticle* mother,
		 int ID,
		 int motherID,
		 int depth=0);
};
#endif

