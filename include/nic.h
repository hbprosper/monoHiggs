#ifndef NIC_H
#define NIC_H
// ---------------------------------------------------------------------------
// A few simple self-contained utilities for monoHiggs analysis.
// A few simple self-contained string utilities
// Created: 16-Jun-2015 HBP & NDP   Les Houches
// ---------------------------------------------------------------------------
#include <string>
#include "TClonesArray.h"
#include "classes/DelphesClasses.h"

///
struct nic
{
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
  ///
  static double leptonIsolation(SortableObject* p,
				TClonesArray* tracks,
				TClonesArray* towers);
};

#endif
