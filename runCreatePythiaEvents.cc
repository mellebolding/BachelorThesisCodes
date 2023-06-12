#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "TFile.h"
#include "TTree.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "include/ProgressBar.h"

#include "include/pythiaEvent.hh"
#include "include/extraInfo.hh"

#include "PU14/CmdLine.hh"
#include "PU14/EventMixer.hh"
#include "PU14/PU14.hh"

#include "include/treeWriter.hh"

using namespace std;
using namespace fastjet;

std::vector<fastjet::PseudoJet> particlesSig;
std::vector<fastjet::PseudoJet> kaonparticles;
std::vector<fastjet::PseudoJet> D0totparticles;
std::vector<fastjet::PseudoJet> D0bartotparticles;
std::vector<fastjet::PseudoJet> pairCred0;
std::vector<fastjet::PseudoJet> pairCred0bar;
std::vector<fastjet::PseudoJet> pionparticles;
std::vector<fastjet::PseudoJet> partons;
std::vector<fastjet::PseudoJet> glusplitd0;
std::vector<fastjet::PseudoJet> glusplitd0bar;
std::vector<fastjet::PseudoJet> cquarks;
std::vector<fastjet::PseudoJet> ctot;
std::vector<fastjet::PseudoJet> cbartot;
std::vector<fastjet::PseudoJet> flavexd0;
std::vector<fastjet::PseudoJet> flavexd0bar;
std::vector<double> nCharvec;
std::vector<double> parId;

int main (int argc, char ** argv)
{
  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  CmdLine cmdline(argc,argv);
  // inputs read from command line
  unsigned int nEvent = cmdline.value<unsigned int>("-nev",1);  // first argument: command line option; second argument: default value 
  
  // Number of events, generated and listed ones.
  //unsigned int nEvent    = 10000;

  //event generator settings
  
  double       ptHat = cmdline.value<double>("-pthat",0);//120.;
  unsigned int tune  = cmdline.value<int>("-tune",14);

  
  std::cout << "generating " << nEvent << " events with pthat = " << ptHat << " and tune = " << tune << std::endl;  

  pythiaEvent pyt(ptHat, tune, -3.0, 3.0);

  
  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);
  

  //output text file
  ofstream fout;
  const char *dir = getenv("PWD");//"/eos/user/m/mverweij/JetWorkshop2017/samples/";
  TString outFileName = Form("%s/PythiaEventsTune%dPtHat%.0f.pu14",dir,tune,ptHat);
 
  fout.open(outFileName.Data());

 
  TFile *file = new TFile(cmdline.value<string>("-output", "results1.root").c_str(), "RECREATE"); //ROOT file
  
  treeWriter trw("pytTree"); // making a TTree file

  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;

  double cross = 0;
  double d0tot_pair = 0;
  
  for(unsigned int ie = 0; ie < nEvent; ie++) {
    
    Bar.Update(ie);
    
    Bar.PrintWithMod(entryDiv);
    
    //---------------------------------------------------------------------------
    //   produce event
    //---------------------------------------------------------------------------
 
    //fout << "# event " << ie << "\n";
    
    //create pythia event
    partons.clear(); //empty list before storing new event
    D0totparticles.clear();
    D0bartotparticles.clear();
    particlesSig.clear();
    kaonparticles.clear();
    pionparticles.clear();
    glusplitd0.clear();
    glusplitd0bar.clear();
    pairCred0.clear();
    pairCred0bar.clear();
    nCharvec.clear();
    parId.clear();
    cquarks.clear();
    ctot.clear();
    cbartot.clear();
    flavexd0.clear();
    flavexd0bar.clear(); 

    std::vector<fastjet::PseudoJet> particlesSig = pyt.createPythiaEvent();
    std::vector<fastjet::PseudoJet> kaonparticles = pyt.getKaonsList();
    std::vector<fastjet::PseudoJet> D0totparticles = pyt.getD0totList();
    std::vector<fastjet::PseudoJet> D0bartotparticles = pyt.getD0bartotList();
    std::vector<fastjet::PseudoJet> pairCred0 = pyt.getpairCred0List();
    std::vector<fastjet::PseudoJet> pairCred0bar = pyt.getpairCred0barList();
    std::vector<fastjet::PseudoJet> pionparticles = pyt.getPionsList();
    std::vector<fastjet::PseudoJet> partons = pyt.getPartonList();
    std::vector<fastjet::PseudoJet> glusplitd0 = pyt.getgluonSplitd0List();
    std::vector<fastjet::PseudoJet> glusplitd0bar = pyt.getgluonSplitd0barList();
    std::vector<fastjet::PseudoJet> cquarks = pyt.getc_QuarkList();
    std::vector<fastjet::PseudoJet> ctot = pyt.getctot();
    std::vector<fastjet::PseudoJet> cbartot = pyt.getcbartot();
    std::vector<fastjet::PseudoJet> flavexd0 = pyt.getflavexd0();
    std::vector<fastjet::PseudoJet> flavexd0bar = pyt.getflavexd0bar();
    std::vector<double> nCharvec;
    std::vector<double> parId;
    nCharvec.push_back(pyt.getNchar());
    nCharvec.push_back(pyt.getNpart());
    nCharvec.push_back(pyt.getb());
    nCharvec.push_back(pyt.getNpions());
    cross += pyt.getCrossT();

    for(fastjet::PseudoJet p : particlesSig) {
      const int & pdgid = p.user_info<extraInfo>().pdg_id();
      const int & vtx   = p.user_info<extraInfo>().vertex_number();
      parId.push_back(p.user_info<extraInfo>().pdg_id());
    }

    //Fill TTree
    trw.addCollection("particlesSig", particlesSig);
    trw.addCollection("partons", partons);
    trw.addCollection("kaons", kaonparticles);
    trw.addCollection("D0tots", D0totparticles);
    trw.addCollection("D0bartots", D0bartotparticles);
    trw.addCollection("pairCred0", pairCred0);
    trw.addCollection("pairCred0bar", pairCred0bar);
    trw.addCollection("pions", pionparticles);
    trw.addCollection("glusplitd0", glusplitd0);
    trw.addCollection("glusplitd0bar", glusplitd0bar);
    trw.addCollection("nChar", nCharvec);
    trw.addCollection("parId", parId);
    trw.addCollection("cQuarks", cquarks);
    trw.addCollection("ctot", ctot);
    trw.addCollection("cbartot", cbartot);
    trw.addCollection("flavexd0", flavexd0);
    trw.addCollection("flavexd0bar", flavexd0bar);
    
    trw.fillTree();
  }

 
  std::cout << " Cross: " << setprecision (6) << cross/(2000*nEvent) << std::endl;

  file = trw.getTree()->GetCurrentFile();
  file->Write();
  file->Close();
  delete file;
 
  fout.close();
  std::cout << "\n Finished generating PYTHIA events" << std::endl;
}
