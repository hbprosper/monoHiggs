#!/usr/bin/env python
#-----------------------------------------------------------------------------
# File: Make a simple struct that can be used in a C++ program to create a
#       simple (flat) Root n-tuple given a file containing the list of
#       variables and their types, with the following format:
#       type list-of-variables
#       example:
#
#          float weight
#          float pt, eta, phi
#          int njets
#
#       Why do this? Because I'm tired of writing the same boilerplate code to
#       write flat ntuples.
#
# Created: 25-Oct-2015 Harrison B. Prosper
#-----------------------------------------------------------------------------
import os, sys
from string import *
from time import ctime
#-----------------------------------------------------------------------------
TEMPLATE = '''#ifndef %(structname)s_H
#define %(structname)s_H

// Created: %(date)s by mktstruct.py

#include <string>
#include <cassert>
#include "TFile.h"
#include "TTree.h"

struct %(structname)s
{
%(variables)s
  //------------------------------------------------------------------------
  %(structname)s()
    : file(0), tree(0), clearvalue(0)
  {}

  %(structname)s(std::string filename,
  %(tab1)sstd::string treename="Analysis",
  %(tab1)sstd::string title="Analysis",
  %(tab1)sfloat clearvalue_=0)
    : file(0), tree(0), clearvalue(clearvalue_)
  {
    Open(filename, treename, title, clearvalue);
  }

  void Open(std::string filename,
            std::string treename="Analysis",
            std::string title="Analysis",
            float clearvalue_=0)  
  {
    file = 0;
    tree = 0;
    clearvalue = clearvalue_;
    
    file = new TFile(filename.c_str(), "recreate");
    assert(file);
    assert(file->IsOpen());
    
    file->cd();
    tree = new TTree(treename.c_str(), title.c_str());
    assert(tree);
    
%(branches)s
  }
  ~%(structname)s() { delete file; }

  void Clear()
  {
%(clear)s
  }
  
  void Fill()
  {
    file->cd();
    tree->Fill();
  }
  
  void Close()
  {
    file->cd();
    tree->Write();
  }
  
  TFile* file;
  TTree* tree;
  float clearvalue;
};
#endif
'''
#-----------------------------------------------------------------------------
def main():
    argv = sys.argv[1:]
    argc = len(argv)
    if argc < 1:
        sys.exit('''
    Usage:
       mktstruct.py variables-file
        ''')

    # get name of file containing variables
    varfile = argv[0]
    if not os.path.exists(varfile):
        sys.exit("** can't open file %s" % varfile)
        
    name = split(varfile, '.')[0]
    structname = capitalize(name)
    names = {'structname': structname,
             'date': ctime(),
             'tab1': ' '*len(structname)+' '
             }

    # read file containing  variables
    records = filter(lambda x: (x[0] != '#') or (x[0] != "/"),
                     filter(lambda x: x != '',
                     map(strip, open(varfile).readlines())))    
    variables = ''
    branches  = ''
    clear     = ''
    for record in records:
        if record[-1] == ';': record = record[:-1]
        t  = split(record)
        ftype = t[0]
        ft = upper(ftype[0])
        
        # get fields
        t  = joinfields(t[1:], '')
        t  = split(t, ",")
        for field in t:
            variables += '  %s\t%s;\n' % (ftype, field)
            branches  += '    tree->Branch("%s", \t&%s, \t"%s/%s");\n' % \
              (field, field, field, ft)
            clear     += '    %s\t= clearvalue;\n' % field               
            print "\t%s\t%s" % (ftype, field)
        
    # write out header
    names['variables'] = variables
    names['branches']  = rstrip(branches)
    names['clear']     = rstrip(clear)
    
    header = "%(structname)s.h" % names
    print "\nwriting %s" % header
    record = TEMPLATE % names
    open(header, "w").write(record)
# -------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print "ciao!"
    
