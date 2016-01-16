#ifndef MONOHIGGS4MU_H
#define MONOHIGGS4MU_H

#include <string>

struct monoHiggs
{
  static void hzz4l(std::string inputFile,
		    std::string sample,
		    int numberEvents=-1,    // all events
		    double luminosity=30.0, // 1/fb
		    double xsection=1.0,
		    int pileup=0);
};
#endif
