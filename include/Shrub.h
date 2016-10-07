#ifndef Shrub_H
#define Shrub_H

// Created: Fri Oct  7 11:05:25 2016 by mktstruct.py

#include <string>
#include <cassert>
#include "TFile.h"
#include "TTree.h"

struct Shrub
{
  float	weight;
  float	genZ1pt;
  float	genZ1eta;
  float	genZ1phi;
  float	genZ1mass;
  float	genZ2pt;
  float	genZ2eta;
  float	genZ2phi;
  float	genZ2mass;
  float	genl1pt;
  float	genl1eta;
  float	genl1phi;
  int	genl1match;
  int	genl1PID;
  float	genl2pt;
  float	genl2eta;
  float	genl2phi;
  int	genl2match;
  int	genl2PID;
  float	genl3pt;
  float	genl3eta;
  float	genl3phi;
  int	genl3match;
  int	genl3PID;
  float	genl4pt;
  float	genl4eta;
  float	genl4phi;
  int	genl4match;
  int	genl4PID;
  float	Z1pt;
  float	Z1eta;
  float	Z1phi;
  float	Z1mass;
  float	Z2pt;
  float	Z2eta;
  float	Z2phi;
  float	Z2mass;
  float	Hpt;
  float	Heta;
  float	Hphi;
  float	Hmass;
  int	nleps;
  float	l1pt;
  float	l1eta;
  float	l1phi;
  int	l1match;
  int	l1PID;
  float	l2pt;
  float	l2eta;
  float	l2phi;
  int	l2match;
  int	l2PID;
  float	l3pt;
  float	l3eta;
  float	l3phi;
  int	l3match;
  int	l3PID;
  float	l4pt;
  float	l4eta;
  float	l4phi;
  int	l4match;
  int	l4PID;
  int	njets;
  float	j1pt;
  float	j1eta;
  float	j1phi;
  float	j1mass;
  float	j2pt;
  float	j2eta;
  float	j2phi;
  float	j2mass;
  float	massjj;
  float	deltaetajj;
  float	met;
  float	costhetastar;
  float	costheta1;
  float	costheta2;
  float	cosPhi;
  float	cosPhi1;

  //------------------------------------------------------------------------
  Shrub()
    : file(0), tree(0), clearvalue(0)
  {}

  Shrub(std::string filename,
        std::string treename="Analysis",
        std::string title="Analysis",
        float clearvalue_=0)
    : file(0), tree(0), clearvalue(clearvalue_)
  {
    Open(filename, treename, title, clearvalue);
  }

  void Open(std::string filename,
            std::string treename,
            std::string title,
            float clearvalue_)  
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
    
    tree->Branch("weight", 	&weight, 	"weight/F");
    tree->Branch("genZ1pt", 	&genZ1pt, 	"genZ1pt/F");
    tree->Branch("genZ1eta", 	&genZ1eta, 	"genZ1eta/F");
    tree->Branch("genZ1phi", 	&genZ1phi, 	"genZ1phi/F");
    tree->Branch("genZ1mass", 	&genZ1mass, 	"genZ1mass/F");
    tree->Branch("genZ2pt", 	&genZ2pt, 	"genZ2pt/F");
    tree->Branch("genZ2eta", 	&genZ2eta, 	"genZ2eta/F");
    tree->Branch("genZ2phi", 	&genZ2phi, 	"genZ2phi/F");
    tree->Branch("genZ2mass", 	&genZ2mass, 	"genZ2mass/F");
    tree->Branch("genl1pt", 	&genl1pt, 	"genl1pt/F");
    tree->Branch("genl1eta", 	&genl1eta, 	"genl1eta/F");
    tree->Branch("genl1phi", 	&genl1phi, 	"genl1phi/F");
    tree->Branch("genl1match", 	&genl1match, 	"genl1match/I");
    tree->Branch("genl1PID", 	&genl1PID, 	"genl1PID/I");
    tree->Branch("genl2pt", 	&genl2pt, 	"genl2pt/F");
    tree->Branch("genl2eta", 	&genl2eta, 	"genl2eta/F");
    tree->Branch("genl2phi", 	&genl2phi, 	"genl2phi/F");
    tree->Branch("genl2match", 	&genl2match, 	"genl2match/I");
    tree->Branch("genl2PID", 	&genl2PID, 	"genl2PID/I");
    tree->Branch("genl3pt", 	&genl3pt, 	"genl3pt/F");
    tree->Branch("genl3eta", 	&genl3eta, 	"genl3eta/F");
    tree->Branch("genl3phi", 	&genl3phi, 	"genl3phi/F");
    tree->Branch("genl3match", 	&genl3match, 	"genl3match/I");
    tree->Branch("genl3PID", 	&genl3PID, 	"genl3PID/I");
    tree->Branch("genl4pt", 	&genl4pt, 	"genl4pt/F");
    tree->Branch("genl4eta", 	&genl4eta, 	"genl4eta/F");
    tree->Branch("genl4phi", 	&genl4phi, 	"genl4phi/F");
    tree->Branch("genl4match", 	&genl4match, 	"genl4match/I");
    tree->Branch("genl4PID", 	&genl4PID, 	"genl4PID/I");
    tree->Branch("Z1pt", 	&Z1pt, 	"Z1pt/F");
    tree->Branch("Z1eta", 	&Z1eta, 	"Z1eta/F");
    tree->Branch("Z1phi", 	&Z1phi, 	"Z1phi/F");
    tree->Branch("Z1mass", 	&Z1mass, 	"Z1mass/F");
    tree->Branch("Z2pt", 	&Z2pt, 	"Z2pt/F");
    tree->Branch("Z2eta", 	&Z2eta, 	"Z2eta/F");
    tree->Branch("Z2phi", 	&Z2phi, 	"Z2phi/F");
    tree->Branch("Z2mass", 	&Z2mass, 	"Z2mass/F");
    tree->Branch("Hpt", 	&Hpt, 	"Hpt/F");
    tree->Branch("Heta", 	&Heta, 	"Heta/F");
    tree->Branch("Hphi", 	&Hphi, 	"Hphi/F");
    tree->Branch("Hmass", 	&Hmass, 	"Hmass/F");
    tree->Branch("nleps", 	&nleps, 	"nleps/I");
    tree->Branch("l1pt", 	&l1pt, 	"l1pt/F");
    tree->Branch("l1eta", 	&l1eta, 	"l1eta/F");
    tree->Branch("l1phi", 	&l1phi, 	"l1phi/F");
    tree->Branch("l1match", 	&l1match, 	"l1match/I");
    tree->Branch("l1PID", 	&l1PID, 	"l1PID/I");
    tree->Branch("l2pt", 	&l2pt, 	"l2pt/F");
    tree->Branch("l2eta", 	&l2eta, 	"l2eta/F");
    tree->Branch("l2phi", 	&l2phi, 	"l2phi/F");
    tree->Branch("l2match", 	&l2match, 	"l2match/I");
    tree->Branch("l2PID", 	&l2PID, 	"l2PID/I");
    tree->Branch("l3pt", 	&l3pt, 	"l3pt/F");
    tree->Branch("l3eta", 	&l3eta, 	"l3eta/F");
    tree->Branch("l3phi", 	&l3phi, 	"l3phi/F");
    tree->Branch("l3match", 	&l3match, 	"l3match/I");
    tree->Branch("l3PID", 	&l3PID, 	"l3PID/I");
    tree->Branch("l4pt", 	&l4pt, 	"l4pt/F");
    tree->Branch("l4eta", 	&l4eta, 	"l4eta/F");
    tree->Branch("l4phi", 	&l4phi, 	"l4phi/F");
    tree->Branch("l4match", 	&l4match, 	"l4match/I");
    tree->Branch("l4PID", 	&l4PID, 	"l4PID/I");
    tree->Branch("njets", 	&njets, 	"njets/I");
    tree->Branch("j1pt", 	&j1pt, 	"j1pt/F");
    tree->Branch("j1eta", 	&j1eta, 	"j1eta/F");
    tree->Branch("j1phi", 	&j1phi, 	"j1phi/F");
    tree->Branch("j1mass", 	&j1mass, 	"j1mass/F");
    tree->Branch("j2pt", 	&j2pt, 	"j2pt/F");
    tree->Branch("j2eta", 	&j2eta, 	"j2eta/F");
    tree->Branch("j2phi", 	&j2phi, 	"j2phi/F");
    tree->Branch("j2mass", 	&j2mass, 	"j2mass/F");
    tree->Branch("massjj", 	&massjj, 	"massjj/F");
    tree->Branch("deltaetajj", 	&deltaetajj, 	"deltaetajj/F");
    tree->Branch("met", 	&met, 	"met/F");
    tree->Branch("costhetastar", 	&costhetastar, 	"costhetastar/F");
    tree->Branch("costheta1", 	&costheta1, 	"costheta1/F");
    tree->Branch("costheta2", 	&costheta2, 	"costheta2/F");
    tree->Branch("cosPhi", 	&cosPhi, 	"cosPhi/F");
    tree->Branch("cosPhi1", 	&cosPhi1, 	"cosPhi1/F");
  }
  ~Shrub() { delete file; }

  void Clear()
  {
    weight	= clearvalue;
    genZ1pt	= clearvalue;
    genZ1eta	= clearvalue;
    genZ1phi	= clearvalue;
    genZ1mass	= clearvalue;
    genZ2pt	= clearvalue;
    genZ2eta	= clearvalue;
    genZ2phi	= clearvalue;
    genZ2mass	= clearvalue;
    genl1pt	= clearvalue;
    genl1eta	= clearvalue;
    genl1phi	= clearvalue;
    genl1match	= clearvalue;
    genl1PID	= clearvalue;
    genl2pt	= clearvalue;
    genl2eta	= clearvalue;
    genl2phi	= clearvalue;
    genl2match	= clearvalue;
    genl2PID	= clearvalue;
    genl3pt	= clearvalue;
    genl3eta	= clearvalue;
    genl3phi	= clearvalue;
    genl3match	= clearvalue;
    genl3PID	= clearvalue;
    genl4pt	= clearvalue;
    genl4eta	= clearvalue;
    genl4phi	= clearvalue;
    genl4match	= clearvalue;
    genl4PID	= clearvalue;
    Z1pt	= clearvalue;
    Z1eta	= clearvalue;
    Z1phi	= clearvalue;
    Z1mass	= clearvalue;
    Z2pt	= clearvalue;
    Z2eta	= clearvalue;
    Z2phi	= clearvalue;
    Z2mass	= clearvalue;
    Hpt	= clearvalue;
    Heta	= clearvalue;
    Hphi	= clearvalue;
    Hmass	= clearvalue;
    nleps	= clearvalue;
    l1pt	= clearvalue;
    l1eta	= clearvalue;
    l1phi	= clearvalue;
    l1match	= clearvalue;
    l1PID	= clearvalue;
    l2pt	= clearvalue;
    l2eta	= clearvalue;
    l2phi	= clearvalue;
    l2match	= clearvalue;
    l2PID	= clearvalue;
    l3pt	= clearvalue;
    l3eta	= clearvalue;
    l3phi	= clearvalue;
    l3match	= clearvalue;
    l3PID	= clearvalue;
    l4pt	= clearvalue;
    l4eta	= clearvalue;
    l4phi	= clearvalue;
    l4match	= clearvalue;
    l4PID	= clearvalue;
    njets	= clearvalue;
    j1pt	= clearvalue;
    j1eta	= clearvalue;
    j1phi	= clearvalue;
    j1mass	= clearvalue;
    j2pt	= clearvalue;
    j2eta	= clearvalue;
    j2phi	= clearvalue;
    j2mass	= clearvalue;
    massjj	= clearvalue;
    deltaetajj	= clearvalue;
    met	= clearvalue;
    costhetastar	= clearvalue;
    costheta1	= clearvalue;
    costheta2	= clearvalue;
    cosPhi	= clearvalue;
    cosPhi1	= clearvalue;
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
