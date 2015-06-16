// main42.cc + bits of main30.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Steve Mrenna.

// Example how to create a copy of the event record, where the original one
// is translated to another format, to meet various analysis needs.
// In this specific case the idea is to set up the history information
// of the underlying hard process to be close to the PYTHIA 6 structure.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include "HepMC/GenEvent.h"   
#include "HepMC/IO_GenEvent.h"

using namespace Pythia8; 
using namespace std;

void translate(Event&, Event&);

int main(int argc, char* argv[]) {

  // Check that correct number of command-line arguments
  if (argc < 2) {
    cerr << " Unexpected number of command-line arguments." << endl
         << " Usage:" << endl
         << "    pythiaRun <card-file> <output-HepMCfile>" << endl;
    return 1;
  }

  // Check that the provided input name corresponds to an existing file.
  ifstream is(argv[1]);  
  if (!is) {
    cerr << " Command-line file " << argv[1] << " was not found. \n"
         << " Program stopped! " << endl;
    return 1;
  }
  // Confirm that external files will be used for input and output.
  cout << "\n >>> PYTHIA settings will be read from file " << argv[1] << endl;


  string HepMCfile(".HepMCfile");
  if (argc > 2)
    {
      HepMCfile = string(argv[2]);
    cout  << " HepMC events will be written to file " 
	  << HepMCfile << " <<< \n" << endl;
    }

  // Interface for conversion from Pythia8::Event to HepMC event. 
  HepMC::Pythia8ToHepMC ToHepMC;
  // Specify file where HepMC events will be stored.
  HepMC::IO_GenEvent ascii_io(HepMCfile.c_str(), std::ios::out);
 
  // Generator. 
  Pythia pythia;

  // Read in commands from external file.
  pythia.readFile(argv[1]);    

  // Extract settings to be used in the main program.
  int    nEvent    = pythia.mode("Main:numberOfEvents");
  int    nAbort    = pythia.mode("Main:timesAllowErrors");
 
  // Initialization.
  pythia.init();

  // Event record for hard interaction and resonance decays.
  Event hard;
  hard.init("(Pythia 6 conventions)", &pythia.particleData);
  
  // Begin event loop.
  int iAbort = 0; 
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate event. 
    if (!pythia.next()) {

      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) {
        cout << " Aborted since reached end of Les Houches Event File\n"; 
        break; 
      }

      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }
    
    if (argc > 2)
      {
	// Reset record for this event. List first few.
	// translate( event, hard );
	
	// Construct new empty HepMC event and fill it.
	// Units will be as chosen for HepMC build, but can be changed
	// by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)  
	HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
	ToHepMC.fill_next_event( pythia, hepmcevt );
	
	// Write the HepMC event to file. Done with it.
	ascii_io << hepmcevt;
	delete hepmcevt;
      }

  // End of event loop. Statistics. 
  }
  pythia.stat();

  // Done.
  return 0;
}

//--------------------------------------------------------------------------

// Routine to identify the hard process of the event after
// accounting for the affects of ISR and resonance decay.

void translate(Event& event, Event& hard) {

    // Reset translated event record.
    hard.reset();

    // Identify primordial partons with kT smearing.
    int iL1 = event[1].daughter1();
    int iL2 = event[2].daughter1();
    int iLast = min(iL1,iL2);

    // Identify initial hard partons.
    int iP1 = 3;
    int iP2 = 4;
    int iEv=-1;

    // Counters to identify ISR radiator, recoil, and radiation.
    int i41 = 0, i42 = 0, i43 = 0;
    while ( iEv < iLast ) {
      iEv++;

      // Identify hard process by status in the 20s.
      int iStatus = event[iEv].status();
      if ( abs(iStatus) > 20 && abs(iStatus) < 30 ) {
        hard.append(event[iEv]);
        if ( iStatus==-21 ) {
          hard[hard.size()-1].mothers(0,0);
        } else {
          hard[hard.size()-1].mothers(1,2);
        }
        continue;
      }

      // Follow the flow of the hard event.
      int iDa1 = event[iEv].daughter2();
      // Test to determine origin of this shower.
      bool hardTest = (iDa1 == iP1) || (iDa1 == iP2);
      // Status -41 is the ISR radiator
      // Status -42 is the ISR recoil
      // Status -61 is the primordial parton with kT smearing
      // Status -53 is a recoil in initial-final dipole radiation
      if ( iStatus == -41 && hardTest ) {
        i41 = iEv;
      } else if ( iStatus == -42 && hardTest ) {
        i42 = iEv;
      } else if ( iStatus == -61 && hardTest ) {
        if ( i41 == 0 ) {
          i41 = iEv;
        } else {
          i42 = iEv;
        }
      } else if ( iStatus == -53 ) {
        if ( iDa1 == iP1 ) {
          iP1 = iEv;
        } else if ( iDa1 == iP2 ) {
          iP2 = iEv;
        }
      }

      if ( !(i41 > 0 && i42 > 0) ) continue;
      int ik2 = event[i41].daughter2();
      int ik1 = event[i42].daughter2();

      if ( event[ik2].pz() > 0 ) {
        int iTemp=ik1;
        ik1=ik2;
        ik2=iTemp;
      }

      // Boost to CM frame of hard partons.
      RotBstMatrix toCMS;
      toCMS.toCMframe(event[ik1].p(),event[ik2].p());

      // Momentum of off-shell incoming parton.
      i43 = event[i41].daughter1();
      Vec4 pBoost = event[i41].p()-event[i43].p();
      RotBstMatrix toHard;
      if ( event[i41].pz() > 0 ) {
        toHard.fromCMframe(pBoost,event[i42].p());
      } else {
        toHard.fromCMframe(event[i42].p(),pBoost);
      }

      // Boost to CM frame of old initiators,
      // then boost from frame of new initiators.
      for( int i = 0; i< hard.size(); ++i ) {
        hard[i].rotbst(toCMS);
        hard[i].rotbst(toHard);
      }

      // Update counter to location of new parton
      iEv = i43;
      // Update event history
      iDa1 = event[i41].daughter2();
      if ( iDa1 == iP1 ) {
        iP1 = i41;
        iP2 = i42;
      } else {
        iP2 = i41;
        iP1 = i42;
      }
      i41=0; i42=0;

    }

    // Handle kT smearing of initial partons.
    RotBstMatrix ref0, ref1;
    if ( event[iP1].pz() > 0 ) {
      ref0.toCMframe(event[iP1].p(),event[iP2].p());
    } else {
      ref0.toCMframe(event[iP2].p(),event[iP1].p());
    }
    if ( event[iL1].pz() > 0 ) {
      ref1.fromCMframe(event[iL1].p(),event[iL2].p());
    } else {
      ref1.fromCMframe(event[iL2].p(),event[iL1].p());
    }
    for( int i=0; i< hard.size(); ++i ) {
      hard[i].rotbst(ref0);
      hard[i].rotbst(ref1);
    }

    hard[1].daughters(3,hard.size()-1);
    hard[2].daughters(3,hard.size()-1);
    int iMax = hard.size();

    // Add resonance decays; start here.
    for( int i = 3; i < iMax && i < event.size(); ++i ) {
      int ilast = -1;
      int ida1, ida2 = -1;
      if ( hard[i].status()!=-22 ) {
        hard[i].statusPos();
        continue;
      }
      if ( hard[i].mother1()==3 ) {
        ilast = hard[i].daughter1();
        ida1 = event[ilast].daughter1();
        ida2 = event[ilast].daughter2();
      } else {
        ilast = i;
        ida1 = hard[ilast].daughter1();
        ida2 = hard[ilast].daughter2(); 
      }

      // Resonance decays occur when there are multiple daughters and
      // it is NOT FSR.
      while( ilast > 0 && ilast < event.size() ) {
        if ( ida1 != ida2 && event[ida1].status() != -51 ) break;
        ilast = ida1;
        ida1 = event[ilast].daughter1();
        ida2 = event[ilast].daughter2();
      }

      // Add daughters to the event record, boosting to the frame
      // of the mother.
      if ( ilast > 0 ) {
        Vec4 pall = Vec4(0,0,0,0);
        for(int ida = ida1; ida <= ida2; ++ida) {
          pall += event[ida].p();
        }
        RotBstMatrix toResonance;
        toResonance.bst( pall, hard[i].p() );
        int nDau = 0;
        for(int ida = ida1; ida <= ida2; ++ida) {
          Particle tmp = Particle(event[ida]);
          tmp.rotbst(toResonance);
          hard.append( tmp );
          hard[hard.size()-1].mothers(i,i);
          iMax++;
          nDau++;
        }
        int ip1 = hard.size() - nDau;
        hard[i].daughters(ip1,ip1+nDau-1);
      }
    }
}

