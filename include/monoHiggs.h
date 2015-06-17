#ifndef MONOHIGGS4MU_H
#define MONOHIGGS4MU_H

#include <string>

struct monoHiggs
{
  static void analysis(std::string inputFile,
		       std::string prefix,
		       int pileup,
		       int numberEvents=-1);
};
#endif
