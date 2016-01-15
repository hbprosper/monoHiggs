// ---------------------------------------------------------------------------
// File: HardScatter
// Description: extract the hard scatter events from the Delphes record.
// created: Les Houches 2015 Nic & HBP
// ---------------------------------------------------------------------------
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <cassert>

#include "HardScatter.h"

using namespace std;

// ---------------------------------------------------------------------------

HardScatter::HardScatter(TClonesArray* particles_)
  : particles(particles_),
    debug(false)
  {}

HardScatter::~HardScatter(){}

void HardScatter::get(vector<LHParticle>& plist, bool debug_)
{
  debug = debug_;
  // find particles with two mothers,
  // then find their status = 62 versions, then
  // climb the decay tree
  for(int c=0; c < particles->GetEntriesFast(); c++)
    {
      GenParticle* p = static_cast<GenParticle*>(particles->At(c));
      if ( p->Status != 22 ) continue;
      if ( p->M1 == p->M2 )  continue;
      if ( p->M1 < 0 )       continue;
      if ( p->M2 < 0 )       continue;

      // particle has two mothers, so continue loop until
      // we find its status = 62 doppelganger
      int ID=c;
      findMother(p, ID);
      p = static_cast<GenParticle*>(particles->At(ID));
 
      // found mother particle with status = 62, now
      // climb decay tree
      int motherID = -1;
      climbTree(plist, p, ID, motherID);      
    }
  
  // create mother -> daughter linkages
  // but first sort according to generation
  // with grandmothers before mothers
  map<int, int> motherMap;
  for(size_t c=0; c < plist.size(); c++)
    motherMap[plist[c].ID] = c;
  
  for(size_t c=0; c < plist.size(); c++)
    {
      int m = plist[c].Mother;
      if ( m < 0 ) continue;
      int motherID = motherMap[m];
      plist[motherID].Daughters.push_back(c);
      plist[c].Mother = motherID;
      plist[c].ID = c;
    }
}

void
HardScatter::findMother(GenParticle* mother, int& pos, int depth)
{
  depth++;
  if ( depth > 100 ) return;
  if ( mother->Status == 62 ) return;
  
  // loop over daughters
  for(int c = mother->D1; c <= mother->D2; c++)
    {
      pos = c;
      GenParticle* d = static_cast<GenParticle*>(particles->At(pos));
      if ( d->PID == mother->PID ) findMother(d, pos, depth);
    }
}
  
void
HardScatter::climbTree(vector<LHParticle>& plist,
		       GenParticle* mother,
		       int ID,
		       int motherID,
		       int depth)
{  
  string tab("...........................................................");
  depth++;
  if ( depth > 50 ) return;
  
  plist.push_back(LHParticle(mother->PID,
			     mother->PT,
			     mother->Eta,
			     mother->Phi,
			     mother->Mass));
  plist.back().ID = ID;
  plist.back().Mother = motherID;
  plist.back().Status = mother->Status;
  
  // note location of current mother in list
  if ( debug )
    cout << tab.substr(0, 2*depth)
	 << plist.back() << endl;
  if ( mother->D2 < 0 ) return;
  
  // loop over daughters
  for(int c = mother->D1; c <= mother->D2; c++)
    {
      GenParticle* d = static_cast<GenParticle*>(particles->At(c));
      climbTree(plist, d, c, ID, depth);
    }
}

