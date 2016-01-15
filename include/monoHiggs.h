#ifndef MONOHIGGS4MU_H
#define MONOHIGGS4MU_H

#include <string>

struct monoHiggs
{
  static void hzz4l(std::string inputFile,
		    std::string prefix,
		    int pileup,
		    double luminosity=30.0, // 1/fb
		    int numberEvents=-1);   // all
};
#endif
