#ifndef pythiaEvent_h
#define pythiaEvent_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "Pythia8/Pythia.h"

#include "extraInfo.hh"
#include "Pythia8/HeavyIons.h"
#include "Pythia8/HIUserHooks.h"

//using namespace std;

//---------------------------------------------------------------
// Description
// This class generates a pythia8 event
// Author: M.Bolding
//---------------------------------------------------------------

class pythiaEvent {

private :
  Pythia8::Pythia pythia;
  double pthat_;
  unsigned int tune_;
  double rapMin_;
  double rapMax_;
  bool   partonLevel_;
  bool   vinciaShower_;  
  
public :
  
  pythiaEvent(double pthat = 120., unsigned int tune = 14, double rapMin = -3., double rapMax = 3., bool partonLevel = false, bool vinciaShower = false);

  std::vector<fastjet::PseudoJet> partons;
  std::vector<fastjet::PseudoJet> kaons;
  std::vector<fastjet::PseudoJet> pions;
  std::vector<fastjet::PseudoJet> D0tot;
  std::vector<fastjet::PseudoJet> D0bartot;
  std::vector<fastjet::PseudoJet> glusplitd0;
  std::vector<fastjet::PseudoJet> glusplitd0bar;
  std::vector<fastjet::PseudoJet> pairCred0;
  std::vector<fastjet::PseudoJet> pairCred0bar;
  std::vector<fastjet::PseudoJet> c_quark;
  std::vector<fastjet::PseudoJet> ctot;
  std::vector<fastjet::PseudoJet> cbartot;
  std::vector<fastjet::PseudoJet> flavexd0;
  std::vector<fastjet::PseudoJet> flavexd0bar;
  
  double nchar;
  double npions;
  double hit;

  std::vector<fastjet::PseudoJet> createPythiaEvent(); 
  std::vector<fastjet::PseudoJet> getPartonList() const { return partons; }
  std::vector<fastjet::PseudoJet> getD0totList() const { return D0tot; }
  std::vector<fastjet::PseudoJet> getD0bartotList() const { return D0bartot; }
  std::vector<fastjet::PseudoJet> getKaonsList() const { return kaons; }
  std::vector<fastjet::PseudoJet> getPionsList() const { return pions; }
  std::vector<fastjet::PseudoJet> getgluonSplitd0List() const { return glusplitd0; }
  std::vector<fastjet::PseudoJet> getgluonSplitd0barList() const { return glusplitd0bar; }
  std::vector<fastjet::PseudoJet> getpairCred0List() const { return pairCred0; }
  std::vector<fastjet::PseudoJet> getpairCred0barList() const { return pairCred0bar; }
  std::vector<fastjet::PseudoJet> getc_QuarkList() const { return c_quark; }
  std::vector<fastjet::PseudoJet> getctot() const { return ctot; }
  std::vector<fastjet::PseudoJet> getcbartot() const { return cbartot; }
  std::vector<fastjet::PseudoJet> getflavexd0() const { return flavexd0; }
  std::vector<fastjet::PseudoJet> getflavexd0bar() const { return flavexd0bar; }

  

  double getNpart() { return (pythia.info.hiInfo->nPartProj() + pythia.info.hiInfo->nPartTarg()); }
  double getNpart2() { return (pythia.info.hiInfo->nAbsProj()+pythia.info.hiInfo->nDiffProj()+pythia.info.hiInfo->nElProj()+pythia.info.hiInfo->nAbsTarg()+pythia.info.hiInfo->nDiffTarg()+pythia.info.hiInfo->nElTarg()); }
  double getSubcoll() { return pythia.info.hiInfo->nCollTot();}
  double getNchar() const { return nchar; }
  double getNpions() const { return npions; }
  double getCrossT() { return pythia.info.hiInfo->sigmaTot(); } 
  double getb() { return pythia.info.hiInfo->b();}

  void getStat() {pythia.stat();}
  double getWeight() {return pythia.info.weight();}



  //Define function that find the initial
  int FindFirst(int i,int id) {
    int firstc = i;
    bool fc = true;
    while(fc){
      for(unsigned int w = 0; w < pythia.event[firstc].motherList().size(); ++w) {
	if (pythia.event[pythia.event[firstc].motherList()[w]].id()==id){
	  firstc  = pythia.event[firstc].motherList()[w];}
	else {
	  fc = false;
	}}}
    return firstc;
  }

  // Function that finds the final
  int FindLast(int i, int id) {
    int lastc = i;
    bool lc = true;
    while(lc){
      for(unsigned int w = 0; w < pythia.event[lastc].daughterList().size(); ++w) {
	if (pythia.event[pythia.event[lastc].daughterList()[w]].id()==id){
	  lastc  = pythia.event[lastc].daughterList()[w];}
	else {
	  lc = false;
	}}}
    return lastc;
  }
 
};
  
pythiaEvent::pythiaEvent(double pthat, unsigned int tune, double rapMin, double rapMax, bool partonLevel, bool vinciaShower) :
  pthat_(pthat), tune_(tune), rapMin_(rapMin), rapMax_(rapMax), partonLevel_(partonLevel), vinciaShower_(vinciaShower)
{
  
  // Generator. LHC process and output selection. Initialization.
  // tunes: http://home.thep.lu.se/~torbjorn/pythia82html/Tunes.html
  pythia.readString("HeavyIon:mode = 1");
  pythia.readString("Beams:eCM = 5020");
  pythia.readString("Beams:idA = 1000822080"); // 1000822080 - Pb
  pythia.readString("Beams:idB = 1000822080"); // 2212 - p
  pythia.readString("HardQCD:all = on");
  pythia.readString("SoftQCD:inelastic = on");
  pythia.readString(Form("PhaseSpace:pTHatMin = %.1f",pthat_));
  pythia.readString("Next:numberShowInfo = 1");
  pythia.readString("Next:numberShowProcess = 1");
  pythia.readString("Next:numberShowEvent = 1");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 1");
  pythia.readString("Angantyr:NucleusModelA = 2");
  pythia.readString("Angantyr:NucleusModelB = 2");
  pythia.readString("HeavyIonA:WSR = 6.62");
  pythia.readString("HeavyIonB:WSR = 6.62");
  pythia.readString("HeavyIonA:WSa = 0.546");
  pythia.readString("HeavyIonB:WSa = 0.546");
  
  if(vinciaShower_)
    pythia.readString("PartonShowers:Model = 2"); //activate the VINCIA parton shower
  else
    pythia.readString(Form("Tune:pp = %d",tune_));
  
  pythia.init();
  
}

std::vector<fastjet::PseudoJet> pythiaEvent::createPythiaEvent() {
  
  
  pythia.next(); //generate next event
  
  std::vector<fastjet::PseudoJet> particles;
  partons.clear(); //empty list before storing new event
  particles.clear();
  kaons.clear();
  pions.clear();
  D0tot.clear();
  D0bartot.clear();
  glusplitd0.clear();
  glusplitd0bar.clear();
  pairCred0.clear();
  pairCred0bar.clear();
  c_quark.clear();
  ctot.clear();
  cbartot.clear();
  flavexd0.clear();
  flavexd0bar.clear();
  nchar = 0; //charged particles
  npions = 0;
  hit = 0;
  int d0tot_paircre = 0;
  
  for (int i = 0; i < pythia.event.size(); ++i) {   

    //CHARGED PARTICLES
    if (pythia.event[i].isFinal() && pythia.event[i].isCharged() && (TMath::Abs(pythia.event[i].eta()) < 0.5) && (pythia.event[i].pT() > 0.1) && TMath::Abs(pythia.event[i].id()) != 13) {
      ++nchar;	 	
      fastjet::PseudoJet p(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e());
      p.set_user_info(new extraInfo(pythia.event[i].id(), 0)); 
      if(p.rap()>rapMin_ && p.rap()<rapMax_){
	particles.push_back(p);
      }
      p.reset(0,0,0,0);
    }

    //D0 selection
    if(pythia.event[i].id()==4 && i == FindLast(i,4)){
	for(unsigned int w = 0; w < pythia.event[i].daughterList().size(); ++w) {
	  if (pythia.event[pythia.event[i].daughterList()[w]].id()==421){
	    int d0c1 = pythia.event[i].daughterList()[w]; // we have a d0
	    int ctop = FindFirst(i,4);
	    for(unsigned int v = 0; v < pythia.event[FindFirst(i,4)].sisterList().size(); ++v) {
	      if (pythia.event[pythia.event[FindFirst(i,4)].sisterList()[v]].id()==-4) {
		int c2 = pythia.event[FindFirst(i,4)].sisterList()[v];
		int c2top = c2;
		// the c bar
		for(unsigned int u = 0; u < pythia.event[FindLast(c2top,-4)].daughterList().size(); ++u) {
		  if (pythia.event[pythia.event[FindLast(c2top,-4)].daughterList()[u]].id()==-421){
		    // we have a d0bar
		    int d0c2 = pythia.event[FindLast(c2top,-4)].daughterList()[u];

		    fastjet::PseudoJet p3(pythia.event[d0c1].px(),pythia.event[d0c1].py(),pythia.event[d0c1].pz(),pythia.event[d0c1].e()); // save this d0/d0bar
		    p3.set_user_info(new extraInfo(pythia.event[d0c1].id(), -1));
		    D0tot.push_back(p3);
		    p3.reset(0,0,0,0);
		    fastjet::PseudoJet p4(pythia.event[d0c2].px(),pythia.event[d0c2].py(),pythia.event[d0c2].pz(),pythia.event[d0c2].e()); // save this d0/d0bar
		    p4.set_user_info(new extraInfo(pythia.event[d0c2].id(), -1));
		    D0bartot.push_back(p4);
		    p4.reset(0,0,0,0);
		    fastjet::PseudoJet p9(pythia.event[ctop].px(),pythia.event[ctop].py(),pythia.event[ctop].pz(),pythia.event[ctop].e()); // save this d0/d0bar
		    p9.set_user_info(new extraInfo(pythia.event[ctop].id(), -1));
		    ctot.push_back(p9);
		    p9.reset(0,0,0,0);
		    fastjet::PseudoJet p10(pythia.event[c2top].px(),pythia.event[c2top].py(),pythia.event[c2top].pz(),pythia.event[c2top].e()); // save this d0/d0bar
		    p10.set_user_info(new extraInfo(pythia.event[c2top].id(), -1));
		    cbartot.push_back(p10);
		    p10.reset(0,0,0,0);

		    //pair creation
		    if ((pythia.event[ctop].status() == -23 && pythia.event[c2top].status() == -23) || (pythia.event[ctop].status() == -33 && pythia.event[c2top].status() == -33)) {
		      fastjet::PseudoJet p5(pythia.event[d0c1].px(),pythia.event[d0c1].py(),pythia.event[d0c1].pz(),pythia.event[d0c1].e());
		      p5.set_user_info(new extraInfo(pythia.event[d0c1].id(), -1));
		      pairCred0.push_back(p5);
		      p5.reset(0,0,0,0);
		      fastjet::PseudoJet p6(pythia.event[d0c2].px(),pythia.event[d0c2].py(),pythia.event[d0c2].pz(),pythia.event[d0c2].e());
		      p6.set_user_info(new extraInfo(pythia.event[d0c2].id(), -1));
		      pairCred0bar.push_back(p6);
		      p6.reset(0,0,0,0);}

		    //flavour excitation
		    if ((pythia.event[ctop].status() == -21 && pythia.event[c2top].status() < -39) || (pythia.event[ctop].status() == -31 && pythia.event[c2top].status() < -39)
			|| (pythia.event[c2top].status() == -21 && pythia.event[ctop].status() < -39) || (pythia.event[c2top].status() == -31 && pythia.event[ctop].status() < -39)) {
		      fastjet::PseudoJet p11(pythia.event[d0c1].px(),pythia.event[d0c1].py(),pythia.event[d0c1].pz(),pythia.event[d0c1].e());
		      p11.set_user_info(new extraInfo(pythia.event[d0c1].id(), -1));
		      flavexd0.push_back(p11);
		      p11.reset(0,0,0,0);
		      fastjet::PseudoJet p12(pythia.event[d0c2].px(),pythia.event[d0c2].py(),pythia.event[d0c2].pz(),pythia.event[d0c2].e());
		      p12.set_user_info(new extraInfo(pythia.event[d0c2].id(), -1));
		      flavexd0bar.push_back(p12);
		      p12.reset(0,0,0,0);}

		    //gluon splitting
		    if (pythia.event[ctop].status() < -39 && pythia.event[c2top].status() < -39) {
		      fastjet::PseudoJet p7(pythia.event[d0c1].px(),pythia.event[d0c1].py(),pythia.event[d0c1].pz(),pythia.event[d0c1].e());
		      p7.set_user_info(new extraInfo(pythia.event[d0c1].id(), -1));
		      glusplitd0.push_back(p7);
		      p7.reset(0,0,0,0);
		      fastjet::PseudoJet p8(pythia.event[d0c2].px(),pythia.event[d0c2].py(),pythia.event[d0c2].pz(),pythia.event[d0c2].e());
		      p8.set_user_info(new extraInfo(pythia.event[d0c2].id(), -1));
		      glusplitd0bar.push_back(p8);
		      p8.reset(0,0,0,0);}
		    
		  }}}}}}}}
  return particles;
}
