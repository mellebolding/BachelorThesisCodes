#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

 // Quadratic background function

Double_t baseline(Double_t *par) {
  return par[0];}

Double_t gaussianLine1(Double_t *x, Double_t *par) {
  return ((par[0]/par[2])*exp((-0.5*pow((x[0]-par[1])/par[2],2))));
}
Double_t gaussianLine2(Double_t *x, Double_t *par) {
  return ((par[0]/par[2])*exp((-0.5*pow((x[0]-par[1])/par[2],2))));
}

// Lorentzian Peak function
Double_t lorentzianPeak(Double_t *x, Double_t *par) {
  return (0.5*par[0]*par[1]/TMath::Pi()) / TMath::Max(1.e-10,(x[0]-par[2])*(x[0]-par[2])+ .25*par[1]*par[1]);
}

// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
  return baseline(par) + gaussianLine1(x,&par[1]) + gaussianLine2(x,&par[4]);
}
Double_t fitFunction2(Double_t *x, Double_t *par) {
  return baseline(par) + gaussianLine1(x,&par[1]);
}


void treeMacro() {
  double binN = 15;
  double binN2 = 24;
   // Create a histogram for the values we read.
   auto myHist = new TH1F("Total","hist1",binN2,-1.29,5);
   auto myHistc = new TH1F("Totalc","hist1c",1,-1.29,5);
   auto myHistc2 = new TH1F("Totalc2","hist1c2",1,-1.29,5);
   auto myHistc3 = new TH1F("Totalc3","hist1c3",1,-1.29,5);
   auto myHiste = new TH1F("Totale","histe1",binN2,-1.29,5);
   auto myHisteb = new TH1F("heb","histeb",binN2,-1.29,5);
   auto myHist2 = new TH1F("h2","hist2",binN2,-1.29,5);
   auto myHiste2 = new TH1F("he2","histe2",binN2,-1.29,5);
   auto myHiste2b = new TH1F("he2b","histe2b",binN2,-1.29,5);
   auto myHiste2v2 = new TH1F("he2v2","histe2v2",binN2,-1.29,5);
   auto myHist3 = new TH1F("h3","hist3",binN2,-1.29,5);
   auto myHiste3 = new TH1F("he3","histe3",binN2+2,-1.29,5);
   auto myHiste3b = new TH1F("he3b","histe3b",binN2+2,-1.29,5);
   auto myHiste3v2 = new TH1F("he3v2","histe3v2",binN2,-1.29,5);
   auto myHistPt = new TH1F("D^{0}-mesons","D^{0}-mesons",35,0,12);
   auto myHist4 = new TH1F("h5","D0Pt",30,-1.29,5);
   auto myHiste4 = new TH1F("he5","histe4",binN2,-1.29,5);
   auto myHiste4b = new TH1F("he4b","histe4b",binN2,-1.29,5);
   auto myHiste4v2 = new TH1F("he4v2","histe4v2",binN2,-1.29,5);
   auto myHist5 = new TH1F("Charm quarks","Charm quarks",35,0,12);
   auto myHist6 = new TH1F("h6","hist6",binN,-1.29,5);
   auto myHist6p = new TH1F("h6p","hist6p",binN,-1.29,5);
   auto myHist6er = new TH1F("h6er","hist6er",15,-1.29,5);
   auto myHist7 = new TH1F("h7","hist7",binN,-1.29,5);
   auto myHist7p = new TH1F("h7p","hist7p",binN,-1.29,5);
   auto myHist8 = new TH1F("h8","hist8",binN,-1.29,5);
   auto myHist8p = new TH1F("h8p","hist8p",binN,-1.29,5);
   auto myHist9 = new TH1F("h9","hist9",binN,-1.29,5);
   auto myHist9p = new TH1F("h9p","hist9p",binN,-1.29,5);
   auto myHist10 = new TH1F("h10","hist10",binN,-1.29,5);
   auto myHist10p = new TH1F("h10p","hist10p",binN,-1.29,5);
   auto myHist11 = new TH1F("h11","hist11",14,-1.29,5);
   auto myHist11p = new TH1F("h11p","hist11p",binN,-1.29,5);
   auto myHist12 = new TH1F("h12","hist12",17,-1.29,5);
   auto myHist12p = new TH1F("h12p","hist12p",binN,-1.29,5);
   auto myHist13 = new TH1F("h13","hist13",binN,-1.29,5);
   auto myHist13p = new TH1F("h13p","hist13p",binN,-1.29,5);
   auto myHist14 = new TH1F("h14","hist14",binN,-1.29,5);
   auto myHist14p = new TH1F("h14p","hist14p",binN,-1.29,5);
   auto myHist15 = new TH1F("h15","hist15",10,-1.29,5);
   auto myHist15p = new TH1F("h15p","hist15p",binN,-1.29,5);
   auto myHist16 = new TH1F("h16","hist16",binN,-1.29,5);
   auto myHist16p = new TH1F("h16p","hist16p",binN,-1.29,5);
   auto myHist17 = new TH1F("h17","hist17",binN,-1.29,5);
   auto myHist17p = new TH1F("h17p","hist17p",binN,-1.29,5);
   auto myHist18 = new TH1F("h18","hist18",binN,-1.29,5);
   auto myHist18p = new TH1F("h18p","hist18p",binN,-1.29,5);
   auto myHist19 = new TH1F("h19","hist19",binN,-1.29,5);
   auto myHist19p = new TH1F("h19p","hist19p",binN,-1.29,5);
   auto myHist20 = new TH1F("h20","hist20",binN,-1.29,5);
   auto myHist20p = new TH1F("h20p","hist20p",binN,-1.29,5);
   auto myHist21 = new TH1F("h21","hist21",17,-1.29,5);
   auto myHist21p = new TH1F("h21p","hist21p",binN,-1.29,5);
   // Open the file containing the tree.
   auto myFile = TFile::Open("results12.root");
   if (!myFile || myFile->IsZombie()) {
      return;
   }
   // Create a TTreeReader for the tree, for instance by passing the
   // TTree's name and the TDirectory / TFile it is in.
   TTreeReader myReader("pytTree", myFile);
   TTreeReaderArray<double> D0totPhi(myReader, "D0totsPhi");
   TTreeReaderArray<double> D0totPt(myReader, "D0totsPt");
   TTreeReaderArray<double> D0bartotPt(myReader, "D0bartotsPt");
   TTreeReaderArray<double> D0totEta(myReader, "D0totsEta");
   TTreeReaderArray<double> D0bartotEta(myReader, "D0bartotsEta");
   TTreeReaderArray<double> nChar(myReader, "nChar");
   TTreeReaderArray<double> ctotPt(myReader, "ctotPt");
   TTreeReaderArray<double> cbartotPt(myReader, "cbartotPt");
   TTreeReaderArray<double> pairCred0Phi(myReader, "pairCred0Phi");
   TTreeReaderArray<double> pairCred0Eta(myReader, "pairCred0Eta");
   TTreeReaderArray<double> pairCred0barEta(myReader, "pairCred0barEta");
   TTreeReaderArray<double> flavexd0barPhi(myReader, "flavexd0barPhi");
   TTreeReaderArray<double> flavexd0barEta(myReader, "flavexd0barEta");
   TTreeReaderArray<double> flavexd0Eta(myReader, "flavexd0Eta");
   TTreeReaderArray<double> flavexd0Phi(myReader, "flavexd0Phi");
   TTreeReaderArray<double> D0bartotPhi(myReader, "D0bartotsPhi");
   TTreeReaderArray<double> glusplitd0Phi(myReader,"glusplitd0Phi");
   TTreeReaderArray<double> glusplitd0barPhi(myReader,"glusplitd0barPhi");
   TTreeReaderArray<double> glusplitd0barEta(myReader,"glusplitd0barEta");
   TTreeReaderArray<double> glusplitd0Eta(myReader,"glusplitd0Eta");
   TTreeReaderArray<double> glusplitd0Pt(myReader,"glusplitd0Pt");
   TTreeReaderArray<double> glusplitd0barPt(myReader,"glusplitd0barPt");
   TTreeReaderArray<double> pairCred0barPhi(myReader, "pairCred0barPhi");
   TTreeReaderArray<double> pairCred0Pt(myReader, "pairCred0Pt");
   TTreeReaderArray<double> pairCred0barPt(myReader, "pairCred0barPt");
   TTreeReaderArray<double> flavexd0barPt(myReader, "flavexd0barPt");
   TTreeReaderArray<double> flavexd0Pt(myReader, "flavexd0Pt");
   
   while (myReader.Next()) {

     /*if(D0totPhi.GetSize() != 0) {
       
       for(int j = 0; j < D0totPhi.GetSize(); ++j) {
       	 ++d0;
	 std::cout << " Pt d0 tot " << D0totPt[j] << std::endl;}}
	 
     if(pairCred0Pt.GetSize() != 0) {
       for(int j = 0; j < pairCred0Pt.GetSize(); ++j) {
	 std::cout << " Pt d0 paircre " << pairCred0Pt[j] << std::endl;}}
     if(glusplitd0Pt.GetSize() != 0) {
       for(int j = 0; j < glusplitd0Pt.GetSize(); ++j) {
	 std::cout << " Pt d0 glu " << glusplitd0Pt[j] << std::endl;}}
     if(flavexd0Pt.GetSize() != 0) {
       for(int j = 0; j < flavexd0Pt.GetSize(); ++j) {
       std::cout << " Pt d0 fla " << flavexd0Pt[j] << std::endl;}}*/
       
     
     if(D0totPhi.GetSize() != 0) {
       if (nChar[2] < 14) {
	 for(int j = 0; j < D0totPhi.GetSize(); ++j) {
	   
	   //DELTA PHI FOR ZERO-PI
	   //double delta_phi = (TMath::Abs(D0totPhi[j] - D0bartotPhi[j]) < 3.1416) * TMath::Abs(D0bartotPhi[j] - D0totPhi[j]) +
	   //(TMath::Abs(D0bartotPhi[j] - D0totPhi[j]) >= 3.1416) * (6.29 - TMath::Abs(D0bartotPhi[j] - D0totPhi[j]));
	   
	   //DELTA PHI FOR -1 - 5
	   double delta_phi6 = (D0bartotPhi[j] - D0totPhi[j] > 0) * (D0bartotPhi[j] - D0totPhi[j]) + (D0bartotPhi[j] - D0totPhi[j] < 0) * (D0bartotPhi[j] - D0totPhi[j] + 6.29);
	   double delta_phi = (delta_phi6 > 5) * (delta_phi6 - 6.29) + (delta_phi6 < 5) * (delta_phi6);
	   
	   // if (D0totPt[j]>0.5 && D0bartotPt[j]>05){
	   myHist->Fill(delta_phi);
	   if ((D0totPt[j] > 8 && D0totPt[j] < 16 && D0bartotPt[j] > 0.3) || (D0totPt[j] > 0.3 && D0bartotPt[j] > 8 && D0bartotPt[j] < 16)) {
	     myHistc->Fill(delta_phi);
	     if (TMath::Abs(D0totEta[j])<0.9 && TMath::Abs(D0bartotEta[j])<0.9){
	       myHiste->Fill(delta_phi);
	       myHistc2->Fill(delta_phi);
	     }
	     if (TMath::Abs(D0totEta[j])<0.5 && TMath::Abs(D0bartotEta[j])<0.5){
	       myHisteb->Fill(delta_phi);
	       myHistc3->Fill(delta_phi);}}
	 }}}
     
     if(pairCred0Phi.GetSize() != 0) {
       if (nChar[2] < 14.0) {
	 for(int j = 0; j < pairCred0Phi.GetSize(); ++j) {
	   
	   //DELTA PHI FOR ZERO-PI
	   //double delta_phi = (TMath::Abs(pairCred0Phi[j] - pairCred0barPhi[j]) < 3.1416) * TMath::Abs(pairCred0barPhi[j] - pairCred0Phi[j]) +
	   //(TMath::Abs(pairCred0barPhi[j] - pairCred0Phi[j]) >= 3.1416) * (6.29 - TMath::Abs(pairCred0barPhi[j] - pairCred0Phi[j]));
	   
	   //DELTA PHI FOR -1.29 - 5
	   double delta_phi6 = (pairCred0Phi[j] - pairCred0barPhi[j] > 0) * (pairCred0Phi[j] - pairCred0barPhi[j]) + (pairCred0Phi[j] - pairCred0barPhi[j] < 0) * (pairCred0Phi[j] - pairCred0barPhi[j] + 6.29);
	   double delta_phi = (delta_phi6 > 5) * (delta_phi6 - 6.29) + (delta_phi6 < 5) * (delta_phi6);
	   
	   // if (pairCred0Pt[j] > 6 && pairCred0barPt[j] > 6) {
	   myHist2->Fill(delta_phi);
	   if (TMath::Abs(pairCred0Eta[j])<0.9 && TMath::Abs(pairCred0barEta[j])<0.9){
	     myHiste2->Fill(delta_phi);
	     myHiste2v2->Fill(delta_phi);}
	   if (TMath::Abs(pairCred0Eta[j])<0.5 && TMath::Abs(pairCred0barEta[j])<0.5){
	     myHiste2b->Fill(delta_phi);}
	 }
       }}
     
     if(glusplitd0barPhi.GetSize() != 0) {
       if (nChar[2] < 14.0) {
	 for(int j = 0; j < glusplitd0barPhi.GetSize(); ++j) {
	   
	   //DELTA PHI FOR ZERO-PI
	   //double delta_phi = (TMath::Abs(glusplitd0barPhi[j] - glusplitd0Phi[j]) < 3.1416) * TMath::Abs(glusplitd0barPhi[j] - glusplitd0Phi[j]) +
	   // (TMath::Abs(glusplitd0barPhi[j] - glusplitd0Phi[j]) >= 3.1416) * (6.29 - TMath::Abs(glusplitd0barPhi[j] - glusplitd0Phi[j]));

	   //DELTA PHI FOR -1.29 - 5
	   double delta_phi6 = (glusplitd0Phi[j] - glusplitd0barPhi[j] > 0) * (glusplitd0Phi[j] - glusplitd0barPhi[j]) + (glusplitd0Phi[j] - glusplitd0barPhi[j] < 0) * (glusplitd0Phi[j] - glusplitd0barPhi[j] + 6.29);
	   double delta_phi = (delta_phi6 > 5) * (delta_phi6 - 6.29) + (delta_phi6 < 5) * (delta_phi6);
	   
	   // if (glusplitd0barPt[j] > 6 && glusplitd0Pt[j] > 6) {
	   myHist3->Fill(delta_phi);
	   if (TMath::Abs(glusplitd0Eta[j])<0.9 && TMath::Abs(glusplitd0barEta[j])<0.9){
	     myHiste3->Fill(delta_phi);
	     myHiste3v2->Fill(delta_phi);}
	   if (TMath::Abs(glusplitd0Eta[j])<0.5 && TMath::Abs(glusplitd0barEta[j])<0.5){
	     myHiste3b->Fill(delta_phi);}
	 }}}
       
     if(flavexd0Phi.GetSize() != 0) {
       if (nChar[2] < 14.0) {
	 for(int j = 0; j < flavexd0Phi.GetSize(); ++j) {
	   
	   //DELTA PHI FOR ZERO-PI
	   // double delta_phi = (TMath::Abs(flavexd0Phi[j] - flavexd0barPhi[j]) < 3.1416) * TMath::Abs(flavexd0barPhi[j] - flavexd0Phi[j]) +
	   //(TMath::Abs(flavexd0barPhi[j] - flavexd0Phi[j]) >= 3.1416) * (6.29 - TMath::Abs(flavexd0barPhi[j] - flavexd0Phi[j]));
	   
	   //DELTA PHI FOR -1.29 - 5
	   double delta_phi6 = (flavexd0Phi[j] - flavexd0barPhi[j] > 0) * (flavexd0Phi[j] - flavexd0barPhi[j]) + (flavexd0Phi[j] - flavexd0barPhi[j] < 0) * (flavexd0Phi[j] - flavexd0barPhi[j] + 6.29);
	   double delta_phi = (delta_phi6 > 5) * (delta_phi6 - 6.29) + (delta_phi6 < 5) * (delta_phi6);
	   
	   myHist4->Fill(delta_phi);
	   if (TMath::Abs(flavexd0Eta[j])<0.9 && TMath::Abs(flavexd0barEta[j])<0.9){
	     myHiste4->Fill(delta_phi);
	     myHiste4v2->Fill(delta_phi);}
	   
	   if (TMath::Abs(flavexd0Eta[j])<0.5 && TMath::Abs(flavexd0barEta[j])<0.5){
	     myHiste4b->Fill(delta_phi);}
	   
	 }}}
       
     if(ctotPt.GetSize() != 0) {
       for(int j = 0; j < ctotPt.GetSize(); ++j) {
	 double ctotpt = ctotPt[j];
	 double cbartotpt = cbartotPt[j];
	 if (ctotPt[j] > 0.1){ // cut for Pt > 0.1 GeV
	   myHist5->Fill(ctotpt);
	   myHistPt->Fill(D0totPt[j]);
	 }
	 if (cbartotPt[j] > 0.1){
	   myHist5->Fill(cbartotpt);
	   myHistPt->Fill(D0bartotPt[j]);}}}


     if(D0totPt.GetSize() != 0 && D0totPhi.GetSize() != 0) {
       if (nChar[2] < 9.92) {
	 for(int j = 0; j < D0totPt.GetSize(); ++j) {
	   if (TMath::Abs(D0totEta[j])<1.0 && TMath::Abs(D0bartotEta[j])<1.0){
	     double delta_phi6 = (D0bartotPhi[j] - D0totPhi[j] > 0) * (D0bartotPhi[j] - D0totPhi[j]) + (D0bartotPhi[j] - D0totPhi[j] < 0) * (D0bartotPhi[j] - D0totPhi[j] + 6.29);
	     double delta_phi = (delta_phi6 > 5) * (delta_phi6 - 6.29) + (delta_phi6 < 5) * (delta_phi6);
	     // double delta_phi = (TMath::Abs(D0totPhi[j] - D0bartotPhi[j]) < 3.1416) * TMath::Abs(D0bartotPhi[j] - D0totPhi[j]) +
	     //(TMath::Abs(D0bartotPhi[j] - D0totPhi[j]) >= 3.1416) * (6.29 - TMath::Abs(D0bartotPhi[j] - D0totPhi[j]));
	     /*if (D0totPt[j] > 0.1) {
	       myHistPt->Fill(D0totPt[j]);}
	       if (D0bartotPt[j] > 0.1) {
	       myHistPt->Fill(D0bartotPt[j]);}*/
	     if ((D0totPt[j] > 2 && D0totPt[j] < 4 && D0bartotPt[j] > 0.1) || (D0totPt[j] > 0.1 && D0bartotPt[j] > 2 && D0bartotPt[j] < 4)) {
	       myHist6->Fill(delta_phi);
	       myHist6er->Fill(delta_phi);}
	     if ((D0totPt[j] > 4 && D0totPt[j] < 6 && D0bartotPt[j] > 0.1) || (D0totPt[j] > 0.1 && D0bartotPt[j] > 4 && D0bartotPt[j] < 6)) {
	       myHist7->Fill(delta_phi);}
	     if ((D0totPt[j] > 6 && D0totPt[j] < 8 && D0bartotPt[j] > 0.1) || (D0totPt[j] > 0.1 && D0bartotPt[j] > 6 && D0bartotPt[j] < 8)) {
	       myHist8->Fill(delta_phi);}
	     if ((D0totPt[j] > 8 && D0totPt[j] < 16 && D0bartotPt[j] > 0.1) || (D0totPt[j] > 0.1 && D0bartotPt[j] > 8 && D0bartotPt[j] < 16)) {
	       myHist9->Fill(delta_phi);}
	     if ((D0totPt[j] > 2 && D0totPt[j] < 4 && D0bartotPt[j] > 0.1 && D0bartotPt[j] < 3) || (D0totPt[j] > 0.1 && D0totPt[j] < 3 && D0bartotPt[j] > 2 && D0bartotPt[j] < 4)) {
	       myHist10->Fill(delta_phi);}
	     if ((D0totPt[j] > 4 && D0totPt[j] < 6 && D0bartotPt[j] > 0.1 && D0bartotPt[j] < 3) || (D0totPt[j] > 0.1 && D0totPt[j] < 3 && D0bartotPt[j] > 4 && D0bartotPt[j] < 6)) {
	       myHist11->Fill(delta_phi);}
	     if ((D0totPt[j] > 6 && D0totPt[j] < 8 && D0bartotPt[j] > 0.1 && D0bartotPt[j] < 3) || (D0totPt[j] > 0.1 && D0totPt[j] < 3 && D0bartotPt[j] > 6 && D0bartotPt[j] < 8)) {
	       myHist12->Fill(delta_phi);}
	     if ((D0totPt[j] > 8 && D0totPt[j] < 16 && D0bartotPt[j] > 0.1 && D0bartotPt[j] < 3) || (D0totPt[j] > 0.1 && D0totPt[j] < 3 && D0bartotPt[j] > 8 && D0bartotPt[j] < 16)) {
	       myHist13->Fill(delta_phi);}
	     if ((D0totPt[j] > 2 && D0totPt[j] < 4 && D0bartotPt[j] > 1 && D0bartotPt[j] < 4) || (D0totPt[j] > 1 && D0totPt[j] < 4 && D0bartotPt[j] > 2 && D0bartotPt[j] < 4)) {
	       myHist14->Fill(delta_phi);}
	     if ((D0totPt[j] > 4 && D0totPt[j] < 6 && D0bartotPt[j] > 1 && D0bartotPt[j] < 4) || (D0totPt[j] > 1 && D0totPt[j] < 4 && D0bartotPt[j] > 4 && D0bartotPt[j] < 6)) {
	       myHist15->Fill(delta_phi);}
	     if ((D0totPt[j] > 6 && D0totPt[j] < 8 && D0bartotPt[j] > 1 && D0bartotPt[j] < 4) || (D0totPt[j] > 1 && D0totPt[j] < 4 && D0bartotPt[j] > 6 && D0bartotPt[j] < 8)) {
	       myHist16->Fill(delta_phi);}
	     if ((D0totPt[j] > 8 && D0totPt[j] < 16 && D0bartotPt[j] > 1 && D0bartotPt[j] < 4) || (D0totPt[j] > 1 && D0totPt[j] < 4 && D0bartotPt[j] > 8 && D0bartotPt[j] < 16)) {
	       myHist17->Fill(delta_phi);}
	     if ((D0totPt[j] > 2 && D0totPt[j] < 4 && D0bartotPt[j] > 2) || (D0totPt[j] > 2 && D0bartotPt[j] > 2 && D0bartotPt[j] < 4)) {
	       myHist18->Fill(delta_phi);}
	     if ((D0totPt[j] > 4 && D0totPt[j] < 6 && D0bartotPt[j] > 2) || (D0totPt[j] > 2 && D0bartotPt[j] > 4 && D0bartotPt[j] < 6)) {
	       myHist19->Fill(delta_phi);}
	     if ((D0totPt[j] > 6 && D0totPt[j] < 8 && D0bartotPt[j] > 2) || (D0totPt[j] > 2 && D0bartotPt[j] > 6 && D0bartotPt[j] < 8)) {
	       myHist20->Fill(delta_phi);}
	     if ((D0totPt[j] > 8 && D0totPt[j] < 16 && D0bartotPt[j] > 2) || (D0totPt[j] > 2 && D0bartotPt[j] > 8 && D0bartotPt[j] < 16)) {
	       myHist21->Fill(delta_phi);}
	   }}}}
   if(pairCred0Pt.GetSize() != 0 && pairCred0Phi.GetSize() != 0) {
       if (nChar[2] < 14.00) {
	 for(int j = 0; j < pairCred0Pt.GetSize(); ++j) {
	   if (TMath::Abs(pairCred0Eta[j])<1.0 && TMath::Abs(pairCred0barEta[j])<1.0){
	      double delta_phi6 = (pairCred0Phi[j] - pairCred0barPhi[j] > 0) * (pairCred0Phi[j] - pairCred0barPhi[j]) + (pairCred0Phi[j] - pairCred0barPhi[j] < 0) * (pairCred0Phi[j] - pairCred0barPhi[j] + 6.29);
	     double delta_phi = (delta_phi6 > 5) * (delta_phi6 - 6.29) + (delta_phi6 < 5) * (delta_phi6);
	     // double delta_phi = (TMath::Abs(D0totPhi[j] - D0bartotPhi[j]) < 3.1416) * TMath::Abs(D0bartotPhi[j] - D0totPhi[j]) +
	     //(TMath::Abs(D0bartotPhi[j] - D0totPhi[j]) >= 3.1416) * (6.29 - TMath::Abs(D0bartotPhi[j] - D0totPhi[j]));
	     /*if (D0totPt[j] > 0.1) {
	       myHistPt->Fill(D0totPt[j]);}
	       if (D0bartotPt[j] > 0.1) {
	       myHistPt->Fill(D0bartotPt[j]);}*/
	     if ((pairCred0Pt[j] > 3 && pairCred0Pt[j] < 6 && pairCred0barPt[j] > 0.3 && pairCred0barPt[j] < 40) || (pairCred0barPt[j] > 3 && pairCred0barPt[j] < 6 && pairCred0Pt[j] > 0.3 && pairCred0Pt[j] < 40)) {
	       myHist6p->Fill(delta_phi);}
	     if ((pairCred0Pt[j] > 4 && pairCred0Pt[j] < 10 && pairCred0barPt[j] > 4 && pairCred0barPt[j] < 10)) {
	       myHist7p->Fill(delta_phi);}
	     if ((pairCred0Pt[j] > 10 && pairCred0Pt[j] < 30 && pairCred0barPt[j] > 10 && pairCred0barPt[j] < 30)) {
	       myHist8p->Fill(delta_phi);}
	     if ((pairCred0Pt[j] > 2 && pairCred0barPt[j] > 2)) {
	       myHist9p->Fill(delta_phi);}
	     if (nChar[2] < 4.28){
	       myHist10p->Fill(delta_phi);}
	     if ((D0totPt[j] > 4 && D0totPt[j] < 6 && D0bartotPt[j] > 0.1 && D0bartotPt[j] < 3) || (D0totPt[j] > 0.1 && D0totPt[j] < 3 && D0bartotPt[j] > 4 && D0bartotPt[j] < 6)) {
	       myHist11p->Fill(delta_phi);}
	     if ((D0totPt[j] > 6 && D0totPt[j] < 8 && D0bartotPt[j] > 0.1 && D0bartotPt[j] < 3) || (D0totPt[j] > 0.1 && D0totPt[j] < 3 && D0bartotPt[j] > 6 && D0bartotPt[j] < 8)) {
	       myHist12p->Fill(delta_phi);}
	     if ((D0totPt[j] > 8 && D0totPt[j] < 16 && D0bartotPt[j] > 0.1 && D0bartotPt[j] < 3) || (D0totPt[j] > 0.1 && D0totPt[j] < 3 && D0bartotPt[j] > 8 && D0bartotPt[j] < 16)) {
	       myHist13p->Fill(delta_phi);}
	     if ((D0totPt[j] > 2 && D0totPt[j] < 4 && D0bartotPt[j] > 1 && D0bartotPt[j] < 4) || (D0totPt[j] > 1 && D0totPt[j] < 4 && D0bartotPt[j] > 2 && D0bartotPt[j] < 4)) {
	       myHist14p->Fill(delta_phi);}
	     if ((D0totPt[j] > 4 && D0totPt[j] < 6 && D0bartotPt[j] > 1 && D0bartotPt[j] < 4) || (D0totPt[j] > 1 && D0totPt[j] < 4 && D0bartotPt[j] > 4 && D0bartotPt[j] < 6)) {
	       myHist15p->Fill(delta_phi);}
	     if ((D0totPt[j] > 6 && D0totPt[j] < 8 && D0bartotPt[j] > 1 && D0bartotPt[j] < 4) || (D0totPt[j] > 1 && D0totPt[j] < 4 && D0bartotPt[j] > 6 && D0bartotPt[j] < 8)) {
	       myHist16p->Fill(delta_phi);}
	     if ((D0totPt[j] > 8 && D0totPt[j] < 16 && D0bartotPt[j] > 1 && D0bartotPt[j] < 4) || (D0totPt[j] > 1 && D0totPt[j] < 4 && D0bartotPt[j] > 8 && D0bartotPt[j] < 16)) {
	       myHist17p->Fill(delta_phi);}
	     if ((D0totPt[j] > 2 && D0totPt[j] < 4 && D0bartotPt[j] > 2) || (D0totPt[j] > 2 && D0bartotPt[j] > 2 && D0bartotPt[j] < 4)) {
	       myHist18p->Fill(delta_phi);}
	     if ((D0totPt[j] > 4 && D0totPt[j] < 6 && D0bartotPt[j] > 2) || (D0totPt[j] > 2 && D0bartotPt[j] > 4 && D0bartotPt[j] < 6)) {
	       myHist19p->Fill(delta_phi);}
	     if ((D0totPt[j] > 6 && D0totPt[j] < 8 && D0bartotPt[j] > 2) || (D0totPt[j] > 2 && D0bartotPt[j] > 6 && D0bartotPt[j] < 8)) {
	       myHist20p->Fill(delta_phi);}
	     if ((D0totPt[j] > 8 && D0totPt[j] < 16 && D0bartotPt[j] > 2) || (D0totPt[j] > 2 && D0bartotPt[j] > 8 && D0bartotPt[j] < 16)) {
	       myHist21p->Fill(delta_phi);}
	   }}}}}
   
   //ERROR ESTIMATION
   double err_corr = .45;
   for (int i = 0; i < 20; ++i) {
     myHist6->SetBinError(i,err_corr*myHist6->GetBinError(i));
     myHist7->SetBinError(i,err_corr*myHist7->GetBinError(i));
     myHist8->SetBinError(i,err_corr*myHist8->GetBinError(i));
     myHist9->SetBinError(i,err_corr*myHist9->GetBinError(i));
     myHist10->SetBinError(i,err_corr*myHist10->GetBinError(i));
     myHist11->SetBinError(i,err_corr*myHist11->GetBinError(i));
     myHist12->SetBinError(i,err_corr*myHist12->GetBinError(i));
     myHist13->SetBinError(i,err_corr*myHist13->GetBinError(i));
     myHist14->SetBinError(i,err_corr*myHist14->GetBinError(i));
     myHist15->SetBinError(i,err_corr*myHist15->GetBinError(i));
     myHist16->SetBinError(i,err_corr*myHist16->GetBinError(i));
     myHist17->SetBinError(i,err_corr*myHist17->GetBinError(i));
     myHist18->SetBinError(i,err_corr*myHist18->GetBinError(i));
     myHist19->SetBinError(i,err_corr*myHist19->GetBinError(i));
     myHist20->SetBinError(i,err_corr*myHist20->GetBinError(i));
     myHist21->SetBinError(i,err_corr*myHist21->GetBinError(i));
     /*std::cout << " old err " << myHist6->GetBinError(i) << std::endl;
       std::cout << " new err " << myHist6er->GetBinError(i) << std::endl;*/
   }

   
   //FIT FOR PT BIN 1
   
   TF1 *fitFcn6 = new TF1("fitFcn6",fitFunction2,-1.29,5,7);
   fitFcn6->SetNpx(300);
   fitFcn6->SetLineColor(kMagenta);
   fitFcn6->SetParameter(3,0.4); // width 1
   fitFcn6->SetParameter(2,0);   // peak 1
   fitFcn6->SetParameter(6,0.6); // width 2
   fitFcn6->SetParameter(5,3.1415);   // peak 2
   fitFcn6->SetParameter(0,0.13); // line height
   double int6 = myHist6->Integral();
   myHist6->Scale(1/int6,"width");
   myHist6->Fit("fitFcn6","V+","ep");
   TF1 *backFcn6 = new TF1("backFcn6",gaussianLine1,-1.29,5,3);
   backFcn6->SetLineColor(kGreen);
   backFcn6->SetLineWidth(1);
   TF1 *signalFcn6 = new TF1("signalFcn6",gaussianLine2,-1.29,5,3);
   signalFcn6->SetLineColor(kBlue);
   signalFcn6->SetNpx(300);
   signalFcn6->SetLineWidth(1);
   Double_t par6[7];
   fitFcn6->GetParameters(par6);
   TLine *line6 = new TLine(-1.29,par6[0],5,par6[0]);
   line6->SetLineColor(kRed);
   line6->SetLineStyle(2);
   line6->SetLineWidth(1);
   backFcn6->SetParameters(&par6[1]);
   signalFcn6->SetParameters(&par6[4]);
   
   
   //FIT FOR PT BIN 2
   
   TF1 *fitFcn7 = new TF1("fitFcn7",fitFunction2,-1.29,5,7);
   fitFcn7->SetNpx(500);
   fitFcn7->SetLineWidth(4);
   fitFcn7->SetLineColor(kMagenta);
   fitFcn7->SetParameter(3,0.8); // width 1
   fitFcn7->SetParameter(2,0);   // peak 1
   fitFcn7->SetParameter(6,0.9); // width 2
   fitFcn7->SetParameter(5,3.1415);   // peak 2
   fitFcn7->SetParameter(0,0.1); // line height
   myHist7->Scale(1/myHist7->Integral(),"width");
   myHist7->Fit("fitFcn7","ME+","ep");
   TF1 *backFcn7 = new TF1("backFcn7",gaussianLine1,-1.29,5,3);
   backFcn7->SetLineColor(kGreen);
   TF1 *signalFcn7 = new TF1("signalFcn7",gaussianLine2,-1.29,5,3);
   signalFcn7->SetLineColor(kBlue);
   signalFcn7->SetNpx(500);
   Double_t par7[7];
   fitFcn7->GetParameters(par7);
   TLine *line7 = new TLine(-1.29,par7[0],5,par7[0]);
   line7->SetLineColor(kRed);
   line7->SetLineStyle(2);
   line7->SetLineWidth(3);
   backFcn7->SetParameters(&par7[1]);
   signalFcn7->SetParameters(&par7[4]);
   
   
   //FIT FOR PT BIN 3
   
   TF1 *fitFcn8 = new TF1("fitFcn8",fitFunction,-1.29,5,7);
   fitFcn8->SetNpx(500);
   fitFcn8->SetLineWidth(4);
   fitFcn8->SetLineColor(kMagenta);
   fitFcn8->SetParameter(3,0.8); // width 1
   fitFcn8->SetParameter(2,0);   // peak 1
   fitFcn8->SetParameter(6,0.2); // width 2
   fitFcn8->SetParameter(5,3.1415);   // peak 2
   fitFcn8->SetParameter(0,0.08); // line height
   myHist8->Scale(1/myHist8->Integral(),"width");
   myHist8->Fit("fitFcn8","ME+","ep");
   TF1 *backFcn8 = new TF1("backFcn8",gaussianLine1,-1.29,5,3);
   backFcn8->SetLineColor(kGreen);
   TF1 *signalFcn8 = new TF1("signalFcn8",gaussianLine2,-1.29,5,3);
   signalFcn8->SetLineColor(kBlue);
   signalFcn8->SetNpx(500);
   Double_t par8[7];
   fitFcn8->GetParameters(par8);
   TLine *line8 = new TLine(-1.29,par8[0],5,par8[0]);
   line8->SetLineColor(kRed);
   line8->SetLineStyle(2);
   line8->SetLineWidth(3);
   backFcn8->SetParameters(&par8[1]);
   signalFcn8->SetParameters(&par8[4]);

   //FIT FOR PT BIN 4
   
   TF1 *fitFcn9 = new TF1("fitFcn9",fitFunction,-1.29,5,7);
   fitFcn9->SetNpx(500);
   fitFcn9->SetLineWidth(4);
   fitFcn9->SetLineColor(kMagenta);
   fitFcn9->SetParameter(3,0.6); // width 1
   fitFcn9->SetParameter(2,0);   // peak 1
   fitFcn9->SetParameter(6,0.4); // width 2
   fitFcn9->SetParameter(5,3.1415);   // peak 2
   fitFcn9->SetParameter(0,0.05); // line height
   myHist9->Scale(1/myHist9->Integral(),"width");
   myHist9->Fit("fitFcn9","ME+","ep");
   TF1 *backFcn9 = new TF1("backFcn9",gaussianLine1,-1.29,5,3);
   backFcn9->SetLineColor(kGreen);
   TF1 *signalFcn9 = new TF1("signalFcn9",gaussianLine2,-1.29,5,3);
   signalFcn9->SetLineColor(kBlue);
   signalFcn9->SetNpx(500);
   Double_t par9[7];
   fitFcn9->GetParameters(par9);
   TLine *line9 = new TLine(-1.29,par9[0],5,par9[0]);
   line9->SetLineColor(kRed);
   line9->SetLineStyle(2);
   line9->SetLineWidth(3);
   backFcn9->SetParameters(&par9[1]);
   signalFcn9->SetParameters(&par9[4]);
   
   //FIT FOR PT BIN 5
   
   TF1 *fitFcn10 = new TF1("fitFcn10",fitFunction2,-1.29,5,7);
   fitFcn10->SetNpx(500);
   fitFcn10->SetLineWidth(4);
   fitFcn10->SetLineColor(kMagenta);
   fitFcn10->SetParameter(3,0.6); // width 1
   fitFcn10->SetParameter(2,0);   // peak 1
   fitFcn10->SetParameter(6,1.8); // width 2
   fitFcn10->SetParameter(5,5);   // peak 2
   fitFcn10->SetParameter(0,0.09); // line height
   myHist10->Scale(1/myHist10->Integral(),"width");
   myHist10->Fit("fitFcn10","V+","ep");
   TF1 *backFcn10 = new TF1("backFcn10",gaussianLine1,-1.29,5,3);
   backFcn10->SetLineColor(kGreen);
   TF1 *signalFcn10 = new TF1("signalFcn10",gaussianLine2,-1.29,5,3);
   signalFcn10->SetLineColor(kBlue);
   signalFcn10->SetNpx(500);
   Double_t par10[7];
   fitFcn10->GetParameters(par10);
   TLine *line10 = new TLine(-1.29,par10[0],5,par10[0]);
   line10->SetLineColor(kRed);
   line10->SetLineStyle(2);
   line10->SetLineWidth(3);
   backFcn10->SetParameters(&par10[1]);
   signalFcn10->SetParameters(&par10[4]);

   //FIT FOR PT BIN 6
   
   TF1 *fitFcn11 = new TF1("fitFcn11",fitFunction2,-1.29,5,7);
   fitFcn11->SetNpx(500);
   fitFcn11->SetLineWidth(4);
   fitFcn11->SetLineColor(kMagenta);
   fitFcn11->SetParameter(3,0.8); // width 1
   fitFcn11->SetParameter(2,0);   // peak 1
   fitFcn11->SetParameter(6,0.8); // width 2
   fitFcn11->SetParameter(5,3.1415);   // peak 2
   fitFcn11->SetParameter(0,0.08); // line height
   myHist11->Scale(1/myHist11->Integral(),"width");
   myHist11->Fit("fitFcn11","ME+","ep");
   TF1 *backFcn11 = new TF1("backFcn11",gaussianLine1,-1.29,5,3);
   backFcn11->SetLineColor(kGreen);
   TF1 *signalFcn11 = new TF1("signalFcn11",gaussianLine2,-1.29,5,3);
   signalFcn11->SetLineColor(kBlue);
   signalFcn11->SetNpx(500);
   Double_t par11[7];
   fitFcn11->GetParameters(par11);
   TLine *line11 = new TLine(-1.29,par11[0],5,par11[0]);
   line11->SetLineColor(kRed);
   line11->SetLineStyle(2);
   line11->SetLineWidth(3);
   backFcn11->SetParameters(&par11[1]);
   signalFcn11->SetParameters(&par11[4]);

   //FIT FOR PT BIN 7
   
   TF1 *fitFcn12 = new TF1("fitFcn12",fitFunction,-1.29,5,7);
   fitFcn12->SetNpx(500);
   fitFcn12->SetLineWidth(4);
   fitFcn12->SetLineColor(kMagenta);
   fitFcn12->SetParameter(3,0.3); // width 1
   fitFcn12->SetParameter(2,0);   // peak 1
   fitFcn12->SetParameter(6,0.25); // width 2
   fitFcn12->SetParameter(5,3.1415);   // peak 2
   fitFcn12->SetParameter(0,0.15); // line height
   myHist12->Scale(1/myHist12->Integral(),"width");
   myHist12->Fit("fitFcn12","ME+","ep");
   TF1 *backFcn12 = new TF1("backFcn12",gaussianLine1,-1.29,5,3);
   backFcn12->SetLineColor(kGreen);
   TF1 *signalFcn12 = new TF1("signalFcn12",gaussianLine2,-1.29,5,3);
   signalFcn12->SetLineColor(kBlue);
   signalFcn12->SetNpx(500);
   Double_t par12[7];
   fitFcn12->GetParameters(par12);
   TLine *line12 = new TLine(-1.29,par12[0],5,par12[0]);
   line12->SetLineColor(kRed);
   line12->SetLineStyle(2);
   line12->SetLineWidth(3);
   backFcn12->SetParameters(&par12[1]);
   signalFcn12->SetParameters(&par12[4]);

   //FIT FOR PT BIN 8

   TF1 *fitFcn13 = new TF1("fitFcn13",fitFunction,-1.29,5,7);
   fitFcn13->SetNpx(500);
   fitFcn13->SetLineWidth(4);
   fitFcn13->SetLineColor(kMagenta);
   fitFcn13->SetParameter(3,0.5); // width 1
   fitFcn13->SetParameter(2,0);   // peak 1
   fitFcn13->SetParameter(6,0.35); // width 2
   fitFcn13->SetParameter(5,3.1415);   // peak 2
   fitFcn13->SetParameter(4,0.025);   // peak 2
   fitFcn13->SetParameter(0,0.06); // line height
   myHist13->Scale(1/myHist13->Integral(),"width");
   myHist13->Fit("fitFcn13","ME+","ep");
   TF1 *backFcn13 = new TF1("backFcn13",gaussianLine1,-1.29,5,3);
   backFcn13->SetLineColor(kGreen);
   TF1 *signalFcn13 = new TF1("signalFcn13",gaussianLine2,-1.29,5,3);
   signalFcn13->SetLineColor(kBlue);
   signalFcn13->SetNpx(500);
   Double_t par13[7];
   fitFcn13->GetParameters(par13);
   TLine *line13 = new TLine(-1.29,par13[0],5,par13[0]);
   line13->SetLineColor(kRed);
   line13->SetLineStyle(2);
   line13->SetLineWidth(3);
   backFcn13->SetParameters(&par13[1]);
   signalFcn13->SetParameters(&par13[4]);

   //FIT FOR PT BIN 9

   TF1 *fitFcn14 = new TF1("fitFcn14",fitFunction2,-1.29,5,7);
   fitFcn14->SetNpx(500);
   fitFcn14->SetLineWidth(4);
   fitFcn14->SetLineColor(kMagenta);
   fitFcn14->SetParameter(3,0.7); // width 1
   fitFcn14->SetParameter(2,0);   // peak 1
   fitFcn14->SetParameter(6,1.2); // width 2
   fitFcn14->SetParameter(5,3.1415);   // peak 2
   fitFcn14->SetParameter(0,0.1); // line height
   myHist14->Scale(1/myHist14->Integral(),"width");
   myHist14->Fit("fitFcn14","V+","ep");
   TF1 *backFcn14 = new TF1("backFcn14",gaussianLine1,-1.29,5,3);
   backFcn14->SetLineColor(kGreen);
   TF1 *signalFcn14 = new TF1("signalFcn14",gaussianLine2,-1.29,5,3);
   signalFcn14->SetLineColor(kBlue);
   signalFcn14->SetNpx(500);
   Double_t par14[7];
   fitFcn14->GetParameters(par14);
   TLine *line14 = new TLine(-1.29,par14[0],5,par14[0]);
   line14->SetLineColor(kRed);
   line14->SetLineStyle(2);
   line14->SetLineWidth(3);
   backFcn14->SetParameters(&par14[1]);
   signalFcn14->SetParameters(&par14[4]);

   //FIT FOR PT BIN 10

   TF1 *fitFcn15 = new TF1("fitFcn15",fitFunction2,-1.29,5,7);
   fitFcn15->SetNpx(500);
   fitFcn15->SetLineWidth(4);
   fitFcn15->SetLineColor(kMagenta);
   fitFcn15->SetParameter(3,0.8); // width 1
   fitFcn15->SetParameter(2,0);   // peak 1
   fitFcn15->SetParameter(6,1.6); // width 2
   fitFcn15->SetParameter(5,3.1415);   // peak 2
   fitFcn15->SetParameter(0,0.11); // line height
   myHist15->Scale(1/myHist15->Integral(),"width");
   myHist15->Fit("fitFcn15","ME+","ep");
   TF1 *backFcn15 = new TF1("backFcn15",gaussianLine1,-1.29,5,3);
   backFcn15->SetLineColor(kGreen);
   TF1 *signalFcn15 = new TF1("signalFcn15",gaussianLine2,-1.29,5,3);
   signalFcn15->SetLineColor(kBlue);
   signalFcn15->SetNpx(500);
   Double_t par15[7];
   fitFcn15->GetParameters(par15);
   TLine *line15 = new TLine(-1.29,par15[0],5,par15[0]);
   line15->SetLineColor(kRed);
   line15->SetLineStyle(2);
   line15->SetLineWidth(3);
   backFcn15->SetParameters(&par15[1]);
   signalFcn15->SetParameters(&par15[4]);

   //FIT FOR PT BIN 11

   TF1 *fitFcn16 = new TF1("fitFcn16",fitFunction,-1.29,5,7);
   fitFcn16->SetNpx(500);
   fitFcn16->SetLineWidth(4);
   fitFcn16->SetLineColor(kMagenta);
   fitFcn16->SetParameter(3,0.3); // width 1
   fitFcn16->SetParameter(2,0);   // peak 1
   fitFcn16->SetParameter(6,0.6); // width 2
   fitFcn16->SetParameter(5,3.1415);   // peak 2
   fitFcn16->SetParameter(0,0.1); // line height
   myHist16->Scale(1/myHist16->Integral(),"width");
   myHist16->Fit("fitFcn16","ME+","ep");
   TF1 *backFcn16 = new TF1("backFcn16",gaussianLine1,-1.29,5,3);
   backFcn16->SetLineColor(kGreen);
   TF1 *signalFcn16 = new TF1("signalFcn16",gaussianLine2,-1.29,5,3);
   signalFcn16->SetLineColor(kBlue);
   signalFcn16->SetNpx(500);
   Double_t par16[7];
   fitFcn16->GetParameters(par16);
   TLine *line16 = new TLine(-1.29,par16[0],5,par16[0]);
   line16->SetLineColor(kRed);
   line16->SetLineStyle(2);
   line16->SetLineWidth(3);
   backFcn16->SetParameters(&par16[1]);
   signalFcn16->SetParameters(&par16[4]);
   
   //FIT FOR PT BIN 12

   TF1 *fitFcn17 = new TF1("fitFcn17",fitFunction,-1.29,5,7);
   fitFcn17->SetNpx(500);
   fitFcn17->SetLineWidth(4);
   fitFcn17->SetLineColor(kMagenta);
   fitFcn17->SetParameter(3,0.3); // width 1
   fitFcn17->SetParameter(2,0);   // peak 1
   fitFcn17->SetParameter(6,0.3); // width 2
   fitFcn17->SetParameter(5,3.1415);   // peak 2
   fitFcn17->SetParameter(0,0.04); // line height
   myHist17->Scale(1/myHist17->Integral(),"width");
   myHist17->Fit("fitFcn17","ME+","ep");
   TF1 *backFcn17 = new TF1("backFcn17",gaussianLine1,-1.29,5,3);
   backFcn17->SetLineColor(kGreen);
   TF1 *signalFcn17 = new TF1("signalFcn17",gaussianLine2,-1.29,5,3);
   signalFcn17->SetLineColor(kBlue);
   signalFcn17->SetNpx(500);
   Double_t par17[7];
   fitFcn17->GetParameters(par17);
   TLine *line17 = new TLine(-1.29,par17[0],5,par17[0]);
   line17->SetLineColor(kRed);
   line17->SetLineStyle(2);
   line17->SetLineWidth(3);
   backFcn17->SetParameters(&par17[1]);
   signalFcn17->SetParameters(&par17[4]);

   //FIT FOR PT BIN 13

   TF1 *fitFcn18 = new TF1("fitFcn18",fitFunction2,-1.29,5,7);
   fitFcn18->SetNpx(500);
   fitFcn18->SetLineWidth(4);
   fitFcn18->SetLineColor(kMagenta);
   fitFcn18->SetParameter(3,0.6); // width 1
   fitFcn18->SetParameter(2,0);   // peak 1
   fitFcn18->SetParameter(6,1.6); // width 2
   fitFcn18->SetParameter(5,3.1415);   // peak 2
   fitFcn18->SetParameter(0,0.1); // line height
   myHist18->Scale(1/myHist18->Integral(),"width");
   myHist18->Fit("fitFcn18","V+","ep");
   TF1 *backFcn18 = new TF1("backFcn18",gaussianLine1,-1.29,5,3);
   backFcn18->SetLineColor(kGreen);
   TF1 *signalFcn18 = new TF1("signalFcn18",gaussianLine2,-1.29,5,3);
   signalFcn18->SetLineColor(kBlue);
   signalFcn18->SetNpx(500);
   Double_t par18[7];
   fitFcn18->GetParameters(par18);
   TLine *line18 = new TLine(-1.29,par18[0],5,par18[0]);
   line18->SetLineColor(kRed);
   line18->SetLineStyle(2);
   line18->SetLineWidth(3);
   backFcn18->SetParameters(&par18[1]);
   signalFcn18->SetParameters(&par18[4]);

   //FIT FOR PT BIN 14

   TF1 *fitFcn19 = new TF1("fitFcn19",fitFunction2,-1.29,5,7);
   fitFcn19->SetNpx(500);
   fitFcn19->SetLineWidth(4);
   fitFcn19->SetLineColor(kMagenta);
   fitFcn19->SetParameter(3,0.4); // width 1
   fitFcn19->SetParameter(2,0);   // peak 1
   fitFcn19->SetParameter(6,0.9); // width 2
   fitFcn19->SetParameter(5,3.1415);   // peak 2
   fitFcn19->SetParameter(0,0.08); // line height
   myHist19->Scale(1/myHist19->Integral(),"width");
   myHist19->Fit("fitFcn19","ME+","ep");
   TF1 *backFcn19 = new TF1("backFcn19",gaussianLine1,-1.29,5,3);
   backFcn19->SetLineColor(kGreen);
   TF1 *signalFcn19 = new TF1("signalFcn19",gaussianLine2,-1.29,5,3);
   signalFcn19->SetLineColor(kBlue);
   signalFcn19->SetNpx(500);
   Double_t par19[7];
   fitFcn19->GetParameters(par19);
   TLine *line19 = new TLine(-1.29,par19[0],5,par19[0]);
   line19->SetLineColor(kRed);
   line19->SetLineStyle(2);
   line19->SetLineWidth(3);
   backFcn19->SetParameters(&par19[1]);
   signalFcn19->SetParameters(&par19[4]);

   //FIT FOR PT BIN 15

   TF1 *fitFcn20 = new TF1("fitFcn20",fitFunction,-1.29,5,7);
   fitFcn20->SetNpx(500);
   fitFcn20->SetLineWidth(4);
   fitFcn20->SetLineColor(kMagenta);
   fitFcn20->SetParameter(3,0.3); // width 1
   fitFcn20->SetParameter(2,0);   // peak 1
   fitFcn20->SetParameter(6,0.3); // width 2
   fitFcn20->SetParameter(5,3.1415);   // peak 2
   fitFcn20->SetParameter(0,0.14); // line height
   myHist20->Scale(1/myHist20->Integral(),"width");
   myHist20->Fit("fitFcn20","ME+","ep");
   TF1 *backFcn20 = new TF1("backFcn20",gaussianLine1,-1.29,5,3);
   backFcn20->SetLineColor(kGreen);
   TF1 *signalFcn20 = new TF1("signalFcn20",gaussianLine2,-1.29,5,3);
   signalFcn20->SetLineColor(kBlue);
   signalFcn20->SetNpx(500);
   Double_t par20[7];
   fitFcn20->GetParameters(par20);
   TLine *line20 = new TLine(-1.29,par20[0],5,par20[0]);
   line20->SetLineColor(kRed);
   line20->SetLineStyle(2);
   line20->SetLineWidth(3);
   backFcn20->SetParameters(&par20[1]);
   signalFcn20->SetParameters(&par20[4]);

   //FIT FOR PT BIN 16

   TF1 *fitFcn21 = new TF1("fitFcn21",fitFunction,-1.29,5,7);
   fitFcn21->SetNpx(500);
   fitFcn21->SetLineWidth(4);
   fitFcn21->SetLineColor(kMagenta);
   fitFcn21->SetParameter(3,0.5); // width 1
   fitFcn21->SetParameter(2,0);   // peak 1
   fitFcn21->SetParameter(6,0.4); // width 2
   fitFcn21->SetParameter(5,3.14);   // peak 2
   fitFcn21->SetParameter(0,0.1); // line height
   myHist21->Scale(1/myHist21->Integral(),"width");
   myHist21->Fit("fitFcn21","ME+","ep");
   TF1 *backFcn21 = new TF1("backFcn21",gaussianLine1,-1.29,5,3);
   backFcn21->SetLineColor(kGreen);
   TF1 *signalFcn21 = new TF1("signalFcn21",gaussianLine2,-1.29,5,3);
   signalFcn21->SetLineColor(kBlue);
   signalFcn21->SetNpx(500);
   Double_t par21[7];
   fitFcn21->GetParameters(par21);
   TLine *line21 = new TLine(-1.29,par21[0],5,par21[0]);
   line21->SetLineColor(kRed);
   line21->SetLineStyle(2);
   line21->SetLineWidth(3);
   backFcn21->SetParameters(&par21[1]);
   signalFcn21->SetParameters(&par21[4]);

   double sig_errns1 =  fitFcn6->GetParError(3)/fitFcn6->GetParameter(3);
   std::cout << "NS sig error 2 " << fitFcn7->GetParError(3)/fitFcn7->GetParameter(3) << std::endl;
   std::cout << "NS sig error 3 " << fitFcn8->GetParError(3)/fitFcn8->GetParameter(3) << std::endl;
   std::cout << "NS sig error 4 " << fitFcn9->GetParError(3)/fitFcn9->GetParameter(3) << std::endl;
   std::cout << "NS sig error 5 " << fitFcn10->GetParError(3)/fitFcn10->GetParameter(3) << std::endl;
   std::cout << "NS sig error 6 " << fitFcn11->GetParError(3)/fitFcn11->GetParameter(3) << std::endl;
   std::cout << "NS sig error 7 " << fitFcn12->GetParError(3)/fitFcn12->GetParameter(3) << std::endl;
   std::cout << "NS sig error 8 " << fitFcn13->GetParError(3)/fitFcn13->GetParameter(3) << std::endl;
   std::cout << "NS sig error 9 " << fitFcn14->GetParError(3)/fitFcn14->GetParameter(3) << std::endl;
   std::cout << "NS sig error 10 " << fitFcn15->GetParError(3)/fitFcn15->GetParameter(3) << std::endl;
   std::cout << "NS sig error 11 " << fitFcn16->GetParError(3)/fitFcn16->GetParameter(3) << std::endl;
   std::cout << "NS sig error 12 " << fitFcn17->GetParError(3)/fitFcn17->GetParameter(3) << std::endl;
   std::cout << "NS sig error 13 " << fitFcn18->GetParError(3)/fitFcn18->GetParameter(3) << std::endl;
   std::cout << "NS sig error 14 " << fitFcn19->GetParError(3)/fitFcn19->GetParameter(3) << std::endl;
   std::cout << "NS sig error 15 " << fitFcn20->GetParError(3)/fitFcn20->GetParameter(3) << std::endl;
   std::cout << "NS sig error 16 " << fitFcn21->GetParError(3)/fitFcn21->GetParameter(3) << std::endl;
   

   std::cout << "AS sig error 1 " << fitFcn6->GetParError(6)/fitFcn6->GetParameter(6) << std::endl;
   std::cout << "AS sig error 2 " << fitFcn7->GetParError(6)/fitFcn7->GetParameter(6) << std::endl;
   std::cout << "AS sig error 3 " << fitFcn8->GetParError(6)/fitFcn8->GetParameter(6) << std::endl;
   std::cout << "AS sig error 4 " << fitFcn9->GetParError(6)/fitFcn9->GetParameter(6) << std::endl;
   std::cout << "AS sig error 5 " << fitFcn10->GetParError(6)/fitFcn10->GetParameter(6) << std::endl;
   std::cout << "AS sig error 6 " << fitFcn11->GetParError(6)/fitFcn11->GetParameter(6) << std::endl;
   std::cout << "AS sig error 7 " << fitFcn12->GetParError(6)/fitFcn12->GetParameter(6) << std::endl;
   std::cout << "AS sig error 8 " << fitFcn13->GetParError(6)/fitFcn13->GetParameter(6) << std::endl;
   std::cout << "AS sig error 9 " << fitFcn14->GetParError(6)/fitFcn14->GetParameter(6) << std::endl;
   std::cout << "AS sig error 10 " << fitFcn15->GetParError(6)/fitFcn15->GetParameter(6) << std::endl;
   std::cout << "AS sig error 11 " << fitFcn16->GetParError(6)/fitFcn16->GetParameter(6) << std::endl;
   std::cout << "AS sig error 12 " << fitFcn17->GetParError(6)/fitFcn17->GetParameter(6) << std::endl;
   std::cout << "AS sig error 13 " << fitFcn18->GetParError(6)/fitFcn18->GetParameter(6) << std::endl;
   std::cout << "AS sig error 14 " << fitFcn19->GetParError(6)/fitFcn19->GetParameter(6) << std::endl;
   std::cout << "AS sig error 15 " << fitFcn20->GetParError(6)/fitFcn20->GetParameter(6) << std::endl;
   std::cout << "AS sig error 16 " << fitFcn21->GetParError(6)/fitFcn21->GetParameter(6) << std::endl;

   std::cout << " chi sq 1 " << fitFcn6->GetChisquare() << std::endl;
   std::cout << " chi sq 2 " << fitFcn7->GetChisquare() << std::endl;
   std::cout << " chi sq 3 " << fitFcn8->GetChisquare() << std::endl;
   std::cout << " chi sq 4 " << fitFcn9->GetChisquare() << std::endl;

   //Pt VS SIGMA PLOT
   double sigerrNS1[] = {TMath::Abs(fitFcn6->GetParError(3)/fitFcn6->GetParameter(3)),TMath::Abs(fitFcn7->GetParError(3)/fitFcn7->GetParameter(3)),TMath::Abs(fitFcn8->GetParError(3)/fitFcn8->GetParameter(3)),TMath::Abs(fitFcn9->GetParError(3)/fitFcn9->GetParameter(3))};
   double sigerrNS2[] = {TMath::Abs(fitFcn10->GetParError(3)/fitFcn10->GetParameter(3)),TMath::Abs(fitFcn11->GetParError(3)/fitFcn11->GetParameter(3)),TMath::Abs(fitFcn12->GetParError(3)/fitFcn12->GetParameter(3)),TMath::Abs(fitFcn13->GetParError(3)/fitFcn13->GetParameter(3))};
   double sigerrNS3[] = {TMath::Abs(fitFcn14->GetParError(3)/fitFcn14->GetParameter(3)),TMath::Abs(fitFcn15->GetParError(3)/fitFcn15->GetParameter(3)),TMath::Abs(fitFcn16->GetParError(3)/fitFcn16->GetParameter(3)),TMath::Abs(fitFcn17->GetParError(3)/fitFcn17->GetParameter(3))};
   double sigerrNS4[] = {TMath::Abs(fitFcn18->GetParError(3)/fitFcn6->GetParameter(3)),TMath::Abs(fitFcn19->GetParError(3)/fitFcn19->GetParameter(3)),TMath::Abs(fitFcn20->GetParError(3)/fitFcn20->GetParameter(3)),TMath::Abs(fitFcn21->GetParError(3)/fitFcn21->GetParameter(3))};
   double sigerrAS1[] = {TMath::Abs(fitFcn6->GetParError(6)/fitFcn6->GetParameter(6)),TMath::Abs(fitFcn7->GetParError(6)/fitFcn7->GetParameter(6)),TMath::Abs(fitFcn8->GetParError(6)/fitFcn8->GetParameter(6)),TMath::Abs(fitFcn9->GetParError(6)/fitFcn9->GetParameter(6))};
   double sigerrAS2[] = {TMath::Abs(fitFcn10->GetParError(6)/fitFcn10->GetParameter(6)),TMath::Abs(fitFcn11->GetParError(6)/fitFcn11->GetParameter(6)),TMath::Abs(fitFcn12->GetParError(6)/fitFcn12->GetParameter(6)),TMath::Abs(fitFcn13->GetParError(6)/fitFcn13->GetParameter(6))};
   double sigerrAS3[] = {TMath::Abs(fitFcn14->GetParError(6)/fitFcn14->GetParameter(6)),TMath::Abs(fitFcn15->GetParError(6)/fitFcn15->GetParameter(6)),TMath::Abs(fitFcn16->GetParError(6)/fitFcn16->GetParameter(6)),TMath::Abs(fitFcn17->GetParError(6)/fitFcn17->GetParameter(6))};
   double sigerrAS4[] = {TMath::Abs(fitFcn18->GetParError(6)/fitFcn18->GetParameter(6)),TMath::Abs(fitFcn19->GetParError(6)/fitFcn19->GetParameter(6)),TMath::Abs(fitFcn20->GetParError(6)/fitFcn20->GetParameter(6)),TMath::Abs(fitFcn21->GetParError(6)/fitFcn21->GetParameter(6))};
   double ptbins[] = {3,5,7,12};
   double interx[] = {1,1,1,4};
   double intery[4];

   auto NSplot1  = new TGraphErrors(4,ptbins,sigerrNS1,interx,intery);
   auto NSplot2  = new TGraphErrors(4,ptbins,sigerrNS2,interx,intery);
   auto NSplot3  = new TGraphErrors(4,ptbins,sigerrNS3,interx,intery);
   auto NSplot4  = new TGraphErrors(4,ptbins,sigerrNS4,interx,intery);
   auto ASplot1  = new TGraphErrors(4,ptbins,sigerrAS1,interx,intery);
   auto ASplot2  = new TGraphErrors(4,ptbins,sigerrAS2,interx,intery);
   auto ASplot3  = new TGraphErrors(4,ptbins,sigerrAS3,interx,intery);
   auto ASplot4  = new TGraphErrors(4,ptbins,sigerrAS4,interx,intery);

   
   //CANVASES
   
   //TCanvas *c1 = new TCanvas("c1","c1");
   //TCanvas *c2 = new TCanvas("c2","c2");
   /*
   TCanvas *c3 = new TCanvas("c3","c3");
   TCanvas *c4 = new TCanvas("c4","c4");
   
   
   TCanvas *c7 = new TCanvas("c7","c7");
   TCanvas *c8 = new TCanvas("c8","c8");
   TCanvas *c9 = new TCanvas("c9","c9");*/

   TCanvas *c6 = new TCanvas("c6","c6");
   gStyle->SetOptStat();
   myHistc->Draw("E1");
   c6->Update();
   TCanvas *c7 = new TCanvas("c7","c7");
   gStyle->SetOptStat();
   myHistc2->Draw("E1");
   c7->Update();
   TCanvas *c8 = new TCanvas("c8","c8");
   gStyle->SetOptStat();
   myHistc3->Draw("E1");
   c8->Update();
   
   auto *c5 = new TCanvas("c5","c5");
   TCanvas *c10 = new TCanvas("c10","c10");
   Double_t w = 900;
   Double_t h = 900;
   TCanvas *c2 = new TCanvas("c2", "c2", w, h);


   c10->SetLeftMargin(0.25);
   c10->SetRightMargin(0.1);
   c10->SetTopMargin(0.1);
   c10->SetBottomMargin(0.25);
   c10->Divide(2,1,0,0);
   c10->cd(1);
   TMultiGraph *mg = new TMultiGraph("mg","; p_{T}^{D} (GeV/c); #epsilon_{#sigma}/#sigma");
   NSplot1->SetTitle("p_{T}^{assoc} > 0.3");
   NSplot1->SetLineColor(kBlack);
   NSplot1->SetLineStyle(1);
   NSplot1->SetMarkerColor(kBlack);
   NSplot1->SetMarkerStyle(2);
   mg->Add(NSplot1);
   NSplot2->SetTitle("0.3 < p_{T}^{assoc} < 3");
   NSplot2->SetLineColor(kRed);
   NSplot2->SetLineStyle(1);
   NSplot2->SetMarkerColor(kRed);
   NSplot2->SetMarkerStyle(2);
   mg->Add(NSplot2); 
   NSplot3->SetTitle("1 < p_{T}^{assoc} < 4");
   NSplot3->SetLineColor(kBlue);
   NSplot3->SetLineStyle(1);
   NSplot3->SetMarkerColor(kBlue);
   NSplot3->SetMarkerStyle(2);
   mg->Add(NSplot3);
   NSplot4->SetTitle("p_{T}^{assoc} > 2");
   NSplot4->SetLineColor(kGreen);
   NSplot4->SetMarkerColor(kGreen);
   NSplot4->SetLineStyle(1);
   NSplot4->SetMarkerStyle(2);
   mg->Add(NSplot4);
   mg->GetYaxis()->SetRangeUser(0,.499);
   mg->GetXaxis()->SetRangeUser(0.01,15.9);
   mg->GetXaxis()->CenterTitle(true);
   mg->GetYaxis()->CenterTitle(true);
   mg->SetTitle("; p_{T} (GeV/c); #sigma/#epsilon_{sigma}");
   mg->Draw("AP");
   c10->cd(1)->BuildLegend();
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();

   c10->cd(2);
   TMultiGraph *ng = new TMultiGraph("ng","; p_{T}^{D} (GeV/c); #epsilon_{#sigma}/#sigma");
   ASplot1->SetLineWidth(1);
   ASplot1->SetMarkerStyle(2);
   ASplot1->SetMarkerColor(kBlack);
   ASplot1->SetLineColor(kBlack);
   ASplot1->SetTitle("Away Side, p_{T}^{assoc} > 0.3");
   ng->Add(ASplot1);
   ASplot2->SetLineWidth(1);
   ASplot2->SetMarkerStyle(2);
   ASplot2->SetMarkerColor(kRed);
   ASplot2->SetLineColor(kRed);
   ASplot2->SetTitle("Away Side, 0.3 < p_{T}^{assoc} < 3");
   ng->Add(ASplot2);
   ASplot3->SetLineWidth(1);
   ASplot3->SetMarkerStyle(2);
   ASplot3->SetMarkerColor(kBlue);
   ASplot3->SetLineColor(kBlue);
   ASplot3->SetTitle("Away Side, 1 < p_{T}^{assoc} < 4");
   ng->Add(ASplot3);
   ASplot4->SetLineWidth(1);
   ASplot4->SetMarkerStyle(2);
   ASplot4->SetMarkerColor(kGreen);
   ASplot4->SetLineColor(kGreen);
   ASplot4->SetTitle("Away Side, p_{T}^{assoc} > 2");
   ng->Add(ASplot4);
   ng->GetYaxis()->SetRangeUser(0,.499);
   ng->GetXaxis()->SetRangeUser(6.001,15.9);
   ng->GetXaxis()->CenterTitle(true);
   ng->GetYaxis()->CenterTitle(true);
   ng->SetTitle("hello ; p_{T} (GeV/c); #sigma/#epsilon_{sigma}");
   ng->Draw("AP");
   
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();

   gStyle->SetOptStat(0);
   TCanvas *c = new TCanvas("c", "multipads", 900,450);
   c->SetLeftMargin(0.25);
   c->SetRightMargin(0.2);
   c->SetTopMargin(0.1);
   c->SetBottomMargin(0.);
   TPad *p1 = new TPad("p1","p1",0,0,0.65,1.,0,0,0);
   p1->SetRightMargin(0);
   p1->SetLeftMargin(0.15);
   p1->SetBottomMargin(0.1);
   p1->Draw(); 
   TPad *p2 = new TPad("p2","p2",0.65,0,1.,1.,0,0,0);
   p2->SetLeftMargin(0);
   p2->SetBottomMargin(0.1);
   p2->Draw();
   gPad->SetTickx();
   gPad->SetTicky();
   p1->cd();
   TMultiGraph *mg2 = new TMultiGraph("mg2","; p_{T}^{D} (GeV/c); #epsilon_{#sigma}/#sigma");
   NSplot1->SetTitle("p_{T}^{assoc} > 0.3");
   NSplot1->SetLineColor(kBlack);
   NSplot1->SetLineStyle(1);
   NSplot1->SetMarkerColor(kBlack);
   NSplot1->SetMarkerStyle(2);
   mg2->Add(NSplot1);
   NSplot2->SetTitle("0.3 < p_{T}^{assoc} < 3");
   NSplot2->SetLineColor(kRed);
   NSplot2->SetLineStyle(1);
   NSplot2->SetMarkerColor(kRed);
   NSplot2->SetMarkerStyle(2);
   mg2->Add(NSplot2); 
   NSplot3->SetTitle("1 < p_{T}^{assoc} < 4");
   NSplot3->SetLineColor(kBlue);
   NSplot3->SetLineStyle(1);
   NSplot3->SetMarkerColor(kBlue);
   NSplot3->SetMarkerStyle(2);
   mg2->Add(NSplot3);
   NSplot4->SetTitle("p_{T}^{assoc} > 2");
   NSplot4->SetLineColor(kGreen);
   NSplot4->SetMarkerColor(kGreen);
   NSplot4->SetLineStyle(1);
   NSplot4->SetMarkerStyle(2);
   mg2->Add(NSplot4);
   mg2->GetYaxis()->SetRangeUser(0,.5);
   mg2->GetXaxis()->SetRangeUser(0.01,15.9);
   mg2->GetXaxis()->CenterTitle(true);
   mg2->GetYaxis()->CenterTitle(true);
   mg2->GetXaxis()->SetTitleSize(0.04);
   mg2->GetXaxis()->SetLabelSize(0.04);
   mg2->GetYaxis()->SetTitleSize(0.05);
   mg2->GetYaxis()->SetLabelSize(0.04);
   mg2->GetXaxis()->SetNdivisions(8,5,0,kTRUE);
   mg2->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   mg2->SetTitle("; p_{T} (GeV/c); #sigma/#epsilon_{sigma}");
   mg2->Draw("AP");
   p1->BuildLegend();
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   auto text10 = new TLatex(0,0.015,"#splitline{#splitline{PYTHIA8, Pb-Pb}{#sqrt{s_{NN}} = 5.02 TeV}}{0-40%, |#eta| < 0.9}");
   auto text11 = new TLatex(0,0.035,"Near side");
   text10->Draw();
   text11->Draw();
   
   p2->cd();
   TMultiGraph *ng2 = new TMultiGraph("ng2","; p_{T}^{D} (GeV/c); #epsilon_{#sigma}/#sigma");
   ASplot1->SetLineWidth(1);
   ASplot1->SetMarkerStyle(2);
   ASplot1->SetMarkerColor(kBlack);
   ASplot1->SetLineColor(kBlack);
   ASplot1->SetTitle("Away Side, p_{T}^{assoc} > 0.3");
   ng2->Add(ASplot1);
   ASplot2->SetLineWidth(1);
   ASplot2->SetMarkerStyle(2);
   ASplot2->SetMarkerColor(kRed);
   ASplot2->SetLineColor(kRed);
   ASplot2->SetTitle("Away Side, 0.3 < p_{T}^{assoc} < 3");
   ng2->Add(ASplot2);
   ASplot3->SetLineWidth(1);
   ASplot3->SetMarkerStyle(2);
   ASplot3->SetMarkerColor(kBlue);
   ASplot3->SetLineColor(kBlue);
   ASplot3->SetTitle("Away Side, 1 < p_{T}^{assoc} < 4");
   ng2->Add(ASplot3);
   ASplot4->SetLineWidth(1);
   ASplot4->SetMarkerStyle(2);
   ASplot4->SetMarkerColor(kGreen);
   ASplot4->SetLineColor(kGreen);
   ASplot4->SetTitle("Away Side, p_{T}^{assoc} > 2");
   ng2->Add(ASplot4);
   ng2->GetYaxis()->SetRangeUser(0,.5);
   ng2->GetXaxis()->SetRangeUser(5.9,16);
   ng2->GetXaxis()->CenterTitle(true);
   ng2->GetYaxis()->CenterTitle(true);
   ng2->GetXaxis()->SetTitleSize(0.055);
   ng2->GetXaxis()->SetLabelSize(0.055);
   ng2->GetYaxis()->SetTitleSize(0.05);
   ng2->GetYaxis()->SetLabelSize(0.06);
   ng2->GetXaxis()->SetLabelOffset(-0.005);
   ng2->GetXaxis()->SetTitleOffset(0.8);
   ng2->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   ng2->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   ng2->SetTitle("hello ; p_{T} (GeV/c); #sigma/#epsilon_{sigma}");
   ng2->Draw("AP");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   auto text12 = new TLatex(7,0.035,"Away side");
   text12->Draw();
  
   TCanvas *c2p = new TCanvas("c2p", "c2p", w, h);
   c2p->SetLeftMargin(0.25);
   c2p->SetRightMargin(0.1);
   c2p->SetTopMargin(0.1);
   c2p->SetBottomMargin(0.25);
   c2p->Divide(4,4,0,0);
   
   c2p->cd(1);
   myHist6p->Draw("hist P E1");
   c2p->cd(2);
   myHist7p->Draw("hist P E1");
   c2p->cd(3);
   myHist8p->Draw("hist P E1");
   c2p->cd(4);
   myHist9p->Draw("hist P E1");
   c2p->cd(5);
   myHist10p->Draw("hist P E1");
   c2p->cd(6);
   myHist11p->Draw("hist P E1");
   c2p->cd(7);
   myHist12p->Draw("hist P E1");
   c2p->cd(8);
   myHist13p->Draw("hist P E1");
   c2p->cd(9);
   myHist14p->Draw("hist P E1");
   c2p->cd(10);
   myHist15p->Draw("hist P E1");

   c2->SetLeftMargin(0.25);
   c2->SetRightMargin(0.1);
   c2->SetTopMargin(0.1);
   c2->SetBottomMargin(0.25);
   c2->Divide(4,4,0,0);
   
   c2->cd(1);
   myHist6->SetTitle("#splitline{2 < p_{T}^{D} < 4}{p_{T}^{assoc} > 0.3};#Delta#phi (rad); (1/N)dN/#Delta#phi");
   myHist6->GetXaxis()->CenterTitle(true);
   myHist6->GetYaxis()->CenterTitle(true);
   myHist6->GetXaxis()->SetTitleSize(0.1);
   myHist6->GetXaxis()->SetLabelSize(0.1);
   myHist6->GetYaxis()->SetTitleSize(0.1);
   myHist6->GetYaxis()->SetLabelSize(0.1);
   myHist6->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHist6->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHist6->SetStats(kTRUE);
   myHist6->GetYaxis()->SetRangeUser(0.0001,0.4999);
   myHist6->SetLineColor(kBlack);
   myHist6->SetLineStyle(1);
   myHist6->SetLineWidth(1);
   myHist6->SetMarkerStyle(2);
   myHist6->SetMarkerSize(1);
   myHist6->SetMarkerColor(kBlack);
   myHist6->Draw("hist P E1");
   fitFcn6->SetLineWidth(3);
   fitFcn6->SetLineColorAlpha(kMagenta, 0.35);
   fitFcn6->Draw("same");
   line6->Draw("same");
   backFcn6->Draw("same");
   signalFcn6->Draw("same");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   c2->cd(1)->BuildLegend();
   c2->cd(1)->BuildLegend();
   c2->cd(2);
   myHist7->SetTitle("#splitline{4 < p_{T}^{D} < 6}{p_{T}^{assoc} > 0.3};#Delta#phi (rad); (1/N)dN/#Delta#phi");
   myHist7->GetXaxis()->CenterTitle(true);
   myHist7->GetYaxis()->CenterTitle(true);
   myHist7->GetXaxis()->SetTitleSize(0.1);
   myHist7->GetXaxis()->SetLabelSize(0.1);
   myHist7->GetYaxis()->SetTitleSize(0.09);
   myHist7->GetYaxis()->SetLabelSize(0.1);
   myHist7->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHist7->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHist7->SetStats(kFALSE);
   myHist7->GetYaxis()->SetRangeUser(0.0001,0.4999);
   myHist7->SetLineColor(kBlack);
   myHist7->SetLineStyle(1);
   myHist7->SetLineWidth(1);
   myHist7->SetMarkerStyle(2);
   myHist7->SetMarkerSize(1);
   myHist7->SetMarkerColor(kBlack);
   fitFcn7->SetLineWidth(3);
   fitFcn7->SetLineColorAlpha(kMagenta, 0.35);
   backFcn7->SetLineWidth(1);
   signalFcn7->SetLineWidth(1);
   line7->SetLineWidth(1);
   myHist7->Draw("hist P E1");
   fitFcn7->Draw("same");
   backFcn7->Draw("same");
   signalFcn7->Draw("same");
   line7->Draw("same");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   auto text = new TLatex(0,0.08,"#splitline{#splitline{PYTHIA8, Pb-Pb}{#sqrt{s_{NN}} = 5.02 TeV}}{#splitline{|#eta| < 0.9}{0-40%}}");
   text->Draw();
   //c2->cd(2)->BuildLegend();
   c2->cd(3);
   myHist8->SetTitle("6 < p_{T}^{D} < 8, p_{T}^{assoc} > 0.3;#Delta#phi (rad); (1/N)dN/#Delta#phi");
   myHist8->GetXaxis()->CenterTitle(true);
   myHist8->GetYaxis()->CenterTitle(true);
   myHist8->GetXaxis()->SetTitleSize(0.1);
   myHist8->GetXaxis()->SetLabelSize(0.1);
   myHist8->GetYaxis()->SetTitleSize(0.08);
   myHist8->GetYaxis()->SetLabelSize(0.1);
   myHist8->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHist8->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHist8->SetStats(kFALSE);
   myHist8->GetYaxis()->SetRangeUser(0.0001,0.4999);
   myHist8->SetLineStyle(1);
   myHist8->SetLineWidth(1);
   myHist8->SetMarkerStyle(2);
   myHist8->SetMarkerSize(1);
   myHist8->SetMarkerColor(kBlack);
   fitFcn8->SetLineWidth(3);
   fitFcn8->SetLineColorAlpha(kMagenta, 0.35);
   backFcn8->SetLineWidth(1);
   signalFcn8->SetLineWidth(1);
   line8->SetLineWidth(1);
   myHist8->Draw("hist P E1");
   fitFcn8->Draw("same");
   backFcn8->Draw("same");
   signalFcn8->Draw("same");
   line8->Draw();
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   c2->cd(4);
   myHist9->SetTitle("8 < p_{T}^{D} < 16, p_{T}^{assoc} > 0.3;#Delta#phi (rad); (1/N)dN/#Delta#phi");
   myHist9->GetXaxis()->CenterTitle(true);
   myHist9->GetYaxis()->CenterTitle(true);
   myHist9->GetXaxis()->SetTitleSize(0.1);
   myHist9->GetXaxis()->SetLabelSize(0.1);
   myHist9->GetYaxis()->SetTitleSize(0.08);
   myHist9->GetYaxis()->SetLabelSize(0.1);
   myHist9->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHist9->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHist9->SetStats(kFALSE);
   myHist9->GetYaxis()->SetRangeUser(0.0001,0.4999);
   myHist9->SetLineStyle(1);
   myHist9->SetLineWidth(1);
   myHist9->SetMarkerStyle(2);
   myHist9->SetMarkerSize(1);
   myHist9->SetMarkerColor(kBlack);
   fitFcn9->SetLineWidth(3);
   fitFcn9->SetLineColorAlpha(kMagenta, 0.35);
   backFcn9->SetLineWidth(1);
   signalFcn9->SetLineWidth(1);
   line9->SetLineWidth(1);
   myHist9->Draw("hist P E1");
   fitFcn9->Draw("same");
   backFcn9->Draw("same");
   signalFcn9->Draw("same");
   line9->Draw();
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   c2->cd(5);
   myHist10->SetTitle("2 < p_{T}^{D} < 4, 0.3 < p_{T}^{assoc} < 3;#Delta#phi (rad); (1/N)dN/#Delta#phi");
   myHist10->GetXaxis()->CenterTitle(true);
   myHist10->GetYaxis()->CenterTitle(true);
   myHist10->GetXaxis()->SetTitleSize(0.1);
   myHist10->GetXaxis()->SetLabelSize(0.1);
   myHist10->GetYaxis()->SetTitleSize(0.1);
   myHist10->GetYaxis()->SetLabelSize(0.1);
   myHist10->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHist10->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHist10->SetStats(kFALSE);
   myHist10->GetYaxis()->SetRangeUser(0.0001,0.4999);
   myHist10->SetLineColor(kBlack);
   myHist10->SetLineStyle(1);
   myHist10->SetLineWidth(1);
   myHist10->SetMarkerStyle(2);
   myHist10->SetMarkerSize(1);
   myHist10->SetMarkerColor(kBlack);
   fitFcn10->SetLineWidth(3);
   fitFcn10->SetLineColorAlpha(kMagenta, 0.35);
   backFcn10->SetLineWidth(1);
   signalFcn10->SetLineWidth(1);
   line10->SetLineWidth(1);
   myHist10->Draw("hist P E1");
   fitFcn10->Draw("same");
   backFcn10->Draw("same");
   signalFcn10->Draw("same");
   line10->Draw("same");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   c2->cd(6);
   myHist11->SetTitle("4 < p_{T}^{D} < 6, 0.3 < p_{T}^{assoc} < 3;#Delta#phi (rad); (1/N)dN/#Delta#phi");
   myHist11->GetXaxis()->CenterTitle(true);
   myHist11->GetYaxis()->CenterTitle(true);
   myHist11->GetXaxis()->SetTitleSize(0.1);
   myHist11->GetXaxis()->SetLabelSize(0.1);
   myHist11->GetYaxis()->SetTitleSize(0.08);
   myHist11->GetYaxis()->SetLabelSize(0.1);
   myHist11->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHist11->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHist11->SetStats(kFALSE);
   myHist11->GetYaxis()->SetRangeUser(0.0001,0.4999);
   myHist11->SetLineColor(kBlack);
   myHist11->SetLineStyle(1);
   myHist11->SetLineWidth(1);
   myHist11->SetMarkerStyle(2);
   myHist11->SetMarkerSize(1);
   myHist11->SetMarkerColor(kBlack);
   fitFcn11->SetLineWidth(3);
   fitFcn11->SetLineColorAlpha(kMagenta, 0.35);
   backFcn11->SetLineWidth(1);
   signalFcn11->SetLineWidth(1);
   line11->SetLineWidth(1);
   myHist11->Draw("hist P E1");
   fitFcn11->Draw("same");
   backFcn11->Draw("same");
   signalFcn11->Draw("same");
   line11->Draw("same");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   c2->cd(7);
   myHist12->SetTitle("6 < p_{T}^{D} < 8, 0.3 < p_{T}^{assoc} < 3;#Delta#phi (rad); (1/N)dN/#Delta#phi");
   myHist12->GetXaxis()->CenterTitle(true);
   myHist12->GetYaxis()->CenterTitle(true);
   myHist12->GetXaxis()->SetTitleSize(0.1);
   myHist12->GetXaxis()->SetLabelSize(0.1);
   myHist12->GetYaxis()->SetTitleSize(0.08);
   myHist12->GetYaxis()->SetLabelSize(0.1);
   myHist12->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHist12->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHist12->SetStats(kFALSE);
   myHist12->GetYaxis()->SetRangeUser(0.0001,0.4999);
   myHist12->SetLineColor(kBlack);
   myHist12->SetLineStyle(1);
   myHist12->SetLineWidth(1);
   myHist12->SetMarkerStyle(2);
   myHist12->SetMarkerSize(1);
   myHist12->SetMarkerColor(kBlack);
   fitFcn12->SetLineWidth(3);
   fitFcn12->SetLineColorAlpha(kMagenta, 0.35);
   backFcn12->SetLineWidth(1);
   signalFcn12->SetLineWidth(1);
   line12->SetLineWidth(1);
   myHist12->Draw("hist P E1");
   fitFcn12->Draw("same");
   backFcn12->Draw("same");
   signalFcn12->Draw("same");
   line12->Draw("same");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   c2->cd(8);
   myHist13->SetTitle("8 < p_{T}^{D} < 16, 0.3 < p_{T}^{assoc} < 3;#Delta#phi (rad); (1/N)dN/#Delta#phi");
   myHist13->GetXaxis()->CenterTitle(true);
   myHist13->GetYaxis()->CenterTitle(true);
   myHist13->GetXaxis()->SetTitleSize(0.1);
   myHist13->GetXaxis()->SetLabelSize(0.1);
   myHist13->GetYaxis()->SetTitleSize(0.08);
   myHist13->GetYaxis()->SetLabelSize(0.1);
   myHist13->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHist13->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHist13->SetStats(kFALSE);
   myHist13->GetYaxis()->SetRangeUser(0.0001,0.4999);
   myHist13->SetLineColor(kBlack);
   myHist13->SetLineStyle(1);
   myHist13->SetLineWidth(1);
   myHist13->SetMarkerStyle(2);
   myHist13->SetMarkerSize(1);
   myHist13->SetMarkerColor(kBlack);
   fitFcn13->SetLineWidth(3);
   fitFcn13->SetLineColorAlpha(kMagenta, 0.35);
   backFcn13->SetLineWidth(1);
   signalFcn13->SetLineWidth(1);
   line13->SetLineWidth(1);
   myHist13->Draw("hist P E1");
   fitFcn13->Draw("same");
   backFcn13->Draw("same");
   signalFcn13->Draw("same");
   line13->Draw("same");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   c2->cd(9);
   myHist14->SetTitle("2 < p_{T}^{D} < 4, 1 < p_{T}^{assoc} < 4;#Delta#phi (rad); (1/N)dN/#Delta#phi");

   myHist14->GetXaxis()->SetTitleSize(0.1);
   myHist14->GetXaxis()->SetLabelSize(0.1);
   myHist14->GetYaxis()->SetTitleSize(0.1);
   myHist14->GetYaxis()->SetLabelSize(0.1);
   myHist14->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHist14->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHist14->SetStats(kFALSE);
   myHist14->GetYaxis()->SetRangeUser(0.0001,0.4999);
   myHist14->SetLineColor(kBlack);
   myHist14->SetLineStyle(1);
   myHist14->SetLineWidth(1);
   myHist14->SetMarkerStyle(2);
   myHist14->SetMarkerSize(1);
   myHist14->SetMarkerColor(kBlack);
   fitFcn14->SetLineWidth(3);
   fitFcn14->SetLineColorAlpha(kMagenta, 0.35);
   backFcn14->SetLineWidth(1);
   signalFcn14->SetLineWidth(1);
   line14->SetLineWidth(1);
   myHist14->Draw("hist P E1");
   fitFcn14->Draw("same");
   backFcn14->Draw("same");
   signalFcn14->Draw("same");
   line14->Draw("same");
   gPad->Modified();
   myHist14->GetXaxis()->CenterTitle(true);
   myHist14->GetYaxis()->CenterTitle(true);
   gPad->SetTickx();
   gPad->SetTicky();
   c2->cd(10);
   myHist15->SetTitle("4 < p_{T}^{D} < 6, 1 < p_{T}^{assoc} < 4;#Delta#phi (rad); (1/N)dN/#Delta#phi");
   myHist15->GetXaxis()->CenterTitle(true);
   myHist15->GetYaxis()->CenterTitle(true);
   myHist15->GetXaxis()->SetTitleSize(0.1);
   myHist15->GetXaxis()->SetLabelSize(0.1);
   myHist15->GetYaxis()->SetTitleSize(0.08);
   myHist15->GetYaxis()->SetLabelSize(0.1);
   myHist15->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHist15->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHist15->SetStats(kFALSE);
   myHist15->GetYaxis()->SetRangeUser(0.0001,0.4999);
   myHist15->SetLineColor(kBlack);
   myHist15->SetLineStyle(1);
   myHist15->SetLineWidth(1);
   myHist15->SetMarkerStyle(2);
   myHist15->SetMarkerSize(1);
   myHist15->SetMarkerColor(kBlack);
   fitFcn15->SetLineWidth(3);
   fitFcn15->SetLineColorAlpha(kMagenta, 0.35);
   backFcn15->SetLineWidth(1);
   signalFcn15->SetLineWidth(1);
   line15->SetLineWidth(1);
   myHist15->Draw("hist P E1");
   fitFcn15->Draw("same");
   backFcn15->Draw("same");
   signalFcn15->Draw("same");
   line15->Draw("same");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   c2->cd(11);
   myHist16->SetTitle("6 < p_{T}^{D} < 8, 1 < p_{T}^{assoc} < 4;#Delta#phi (rad); (1/N)dN/#Delta#phi");
   myHist16->GetXaxis()->CenterTitle(true);
   myHist16->GetYaxis()->CenterTitle(true);
   myHist16->GetXaxis()->SetTitleSize(0.1);
   myHist16->GetXaxis()->SetLabelSize(0.1);
   myHist16->GetYaxis()->SetTitleSize(0.08);
   myHist16->GetYaxis()->SetLabelSize(0.1);
   myHist16->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHist16->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHist16->SetStats(kFALSE);
   myHist16->GetYaxis()->SetRangeUser(0.0001,0.4999);
   myHist16->SetLineColor(kBlack);
   myHist16->SetLineStyle(1);
   myHist16->SetLineWidth(1);
   myHist16->SetMarkerStyle(2);
   myHist16->SetMarkerSize(1);
   myHist16->SetMarkerColor(kBlack);
   fitFcn16->SetLineWidth(3);
   fitFcn16->SetLineColorAlpha(kMagenta, 0.35);
   backFcn16->SetLineWidth(1);
   signalFcn16->SetLineWidth(1);
   line16->SetLineWidth(1);
   myHist16->Draw("hist P E1");
   fitFcn16->Draw("same");
   backFcn16->Draw("same");
   signalFcn16->Draw("same");
   line16->Draw("same");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   c2->cd(12);
   myHist17->SetTitle("8 < p_{T}^{D} < 16, 1 < p_{T}^{assoc} < 4;#Delta#phi (rad); (1/N)dN/#Delta#phi");
   myHist17->GetXaxis()->CenterTitle(true);
   myHist17->GetYaxis()->CenterTitle(true);
   myHist17->GetXaxis()->SetTitleSize(0.1);
   myHist17->GetXaxis()->SetLabelSize(0.1);
   myHist17->GetYaxis()->SetTitleSize(0.08);
   myHist17->GetYaxis()->SetLabelSize(0.08);
   myHist17->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHist17->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHist17->SetStats(kFALSE);
   myHist17->GetYaxis()->SetRangeUser(0.0001,0.4999);
   myHist17->SetLineColor(kBlack);
   myHist17->SetLineStyle(1);
   myHist17->SetLineWidth(1);
   myHist17->SetMarkerStyle(2);
   myHist17->SetMarkerSize(1);
   myHist17->SetMarkerColor(kBlack);
   fitFcn17->SetLineWidth(3);
   fitFcn17->SetLineColorAlpha(kMagenta, 0.35);
   backFcn17->SetLineWidth(1);
   signalFcn17->SetLineWidth(1);
   line17->SetLineWidth(1);
   myHist17->Draw("hist P E1");
   fitFcn17->Draw("same");
   backFcn17->Draw("same");
   signalFcn17->Draw("same");
   line17->Draw("same");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   c2->cd(13);
   myHist18->SetTitle("2 < p_{T}^{D} < 4, p_{T}^{assoc} > 2;#Delta#phi (rad); (1/N)dN/#Delta#phi");
   myHist18->GetXaxis()->CenterTitle(true);
   myHist18->GetYaxis()->CenterTitle(true);
   myHist18->GetXaxis()->SetTitleSize(0.08);
   myHist18->GetXaxis()->SetLabelSize(0.07);
   myHist18->GetYaxis()->SetTitleSize(0.08);
   myHist18->GetYaxis()->SetLabelSize(0.07);
   myHist18->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHist18->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHist18->SetStats(kFALSE);
   myHist18->GetYaxis()->SetRangeUser(0.0,0.4999);
   myHist18->SetLineColor(kBlack);
   myHist18->SetLineStyle(1);
   myHist18->SetLineWidth(1);
   myHist18->SetMarkerStyle(2);
   myHist18->SetMarkerSize(1);
   myHist18->SetMarkerColor(kBlack);
   fitFcn18->SetLineWidth(3);
   fitFcn18->SetLineColorAlpha(kMagenta, 0.35);
   backFcn18->SetLineWidth(1);
   signalFcn18->SetLineWidth(1);
   line18->SetLineWidth(1);
   myHist18->Draw("hist P E1");
   fitFcn18->Draw("same");
   backFcn18->Draw("same");
   signalFcn18->Draw("same");
   line18->Draw("same");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();

   c2->cd(14);
   myHist19->SetTitle("4 < p_{T}^{D} < 6, p_{T}^{assoc} > 2;#Delta#phi (rad); (1/N)dN/#Delta#phi");
   myHist19->GetXaxis()->CenterTitle(true);
   myHist19->GetYaxis()->CenterTitle(true);
   myHist19->GetXaxis()->SetTitleSize(0.1);
   myHist19->GetXaxis()->SetLabelSize(0.1);
   myHist19->GetYaxis()->SetTitleSize(0.08);
   myHist19->GetYaxis()->SetLabelSize(0.08);
   myHist19->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHist19->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHist19->SetStats(kFALSE);
   myHist19->GetYaxis()->SetRangeUser(0.0,0.4999);
   myHist19->SetLineColor(kBlack);
   myHist19->SetLineStyle(1);
   myHist19->SetLineWidth(1);
   myHist19->SetMarkerStyle(2);
   myHist19->SetMarkerSize(1);
   myHist19->SetMarkerColor(kBlack);
   fitFcn19->SetLineWidth(3);
   fitFcn19->SetLineColorAlpha(kMagenta, 0.35);
   backFcn19->SetLineWidth(1);
   signalFcn19->SetLineWidth(1);
   line19->SetLineWidth(1);
   myHist19->Draw("hist P E1");
   fitFcn19->Draw("same");
   backFcn19->Draw("same");
   signalFcn19->Draw("same");
   line19->Draw("same");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   c2->cd(15);
   myHist20->SetTitle("#splitline{6 < p_{T}^{D} < 8}{p_{T}^{assoc} > 2};#Delta#phi (rad); (1/N)dN/#Delta#phi");
   myHist20->GetXaxis()->CenterTitle(true);
   myHist20->GetYaxis()->CenterTitle(true);
   myHist20->GetXaxis()->SetTitleSize(0.1);
   myHist20->GetXaxis()->SetLabelSize(0.1);
   myHist20->GetYaxis()->SetTitleSize(0.08);
   myHist20->GetYaxis()->SetLabelSize(0.08);
   myHist20->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHist20->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHist20->SetStats(kFALSE);
   myHist20->GetYaxis()->SetRangeUser(0.0,0.4999);
   myHist20->SetLineColor(kBlack);
   myHist20->SetLineStyle(1);
   myHist20->SetLineWidth(1);
   myHist20->SetMarkerStyle(2);
   myHist20->SetMarkerSize(1);
   myHist20->SetMarkerColor(kBlack);
   fitFcn20->SetLineWidth(3);
   fitFcn20->SetLineColorAlpha(kMagenta, 0.35);
   backFcn20->SetLineWidth(1);
   signalFcn20->SetLineWidth(1);
   line20->SetLineWidth(1);
   myHist20->Draw("hist P E1");
   fitFcn20->Draw("same");
   backFcn20->Draw("same");
   signalFcn20->Draw("same");
   line20->Draw("same");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   c2->cd(16);
   myHist21->SetTitle("#splitline{8 < p_{T}^{D} < 16}{p_{T}^{assoc} > 2};#Delta#phi (rad); (1/N)dN/#Delta#phi");
   myHist21->GetXaxis()->CenterTitle(true);
   myHist21->GetYaxis()->CenterTitle(true);
   myHist21->GetXaxis()->SetTitleSize(0.1);
   myHist21->GetXaxis()->SetLabelSize(0.1);
   myHist21->GetYaxis()->SetTitleSize(0.08);
   myHist21->GetYaxis()->SetLabelSize(0.08);
   myHist21->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHist21->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHist21->SetStats(kFALSE);
   myHist21->GetYaxis()->SetRangeUser(0.0,0.4999);
   myHist21->SetLineColor(kBlack);
   myHist21->SetLineStyle(1);
   myHist21->SetLineWidth(1);
   myHist21->SetMarkerStyle(2);
   myHist21->SetMarkerSize(1);
   myHist21->SetMarkerColor(kBlack);
   fitFcn21->SetLineWidth(3);
   fitFcn21->SetLineColorAlpha(kMagenta, 0.35);
   backFcn21->SetLineWidth(1);
   signalFcn21->SetLineWidth(1);
   line21->SetLineWidth(1);
   myHist21->Draw("hist P E1");
   fitFcn21->Draw("same");
   backFcn21->Draw("same");
   signalFcn21->Draw("same");
   line21->Draw("same");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   

   
   c5->cd();
   double inte = myHist5->Integral();
   myHist5->SetLineColor(kBlack);
   myHist5->SetLineWidth(2);
   myHist5->SetLineStyle(1);
   myHist5->SetMarkerStyle(21);
   myHist5->SetMarkerSize(1.5);
   myHist5->SetMarkerColor(kBlack);
   myHist5->Scale(1./ (inte),"width");
   myHist5->SetTitle(";#it{p_{T}} (GeV); (1/N_{D})dN_{D}/d#it{p_{T}}");
   myHist5->SetStats(kFALSE);
   myHist5->Draw("hist PC E0");
   myHistPt->Scale(1. / (inte), "width");
   myHistPt->SetLineWidth(2);
   myHistPt->SetLineColor(kRed);
   myHistPt->SetLineStyle(1);
   myHistPt->SetMarkerStyle(21);
   myHistPt->SetMarkerSize(1.5);
   myHistPt->SetMarkerColor(kRed);
   
   myHistPt->Draw("same hist PC E0");
   
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   myHist5->GetXaxis()->SetLimits(0.1,12);
   
   c5->BuildLegend();

   TLatex* ltx1 = new TLatex();
   ltx1->DrawLatex(8,0.03,"#it{p_{T,charm}} > 0.1 GeV/c");
   c5->Update();
   
   
   
   
   /* c6->cd();
   myHist6->SetTitle("2 < p_{T}^{D} < 4, p_{T}^{assoc} > 0.3 ;#Delta#phi (rad); (1/N)dN/#Delta#phi");
   myHist6->SetStats(kFALSE);
   myHist6->GetYaxis()->SetRangeUser(0,0.25);
   myHist6->SetLineColor(kBlack);
   myHist6->SetLineStyle(2);
   myHist6->SetLineWidth(2);
   myHist6->SetMarkerStyle(50);
   myHist6->SetMarkerSize(2);
   myHist6->SetMarkerColor(kBlack);
   myHist6->Draw("hist P");
   fitFcn2->Draw("same");
   backFcn2->Draw("same");
   signalFcn2->Draw("same");
   line2->Draw("same");  
   c7->cd();
   myHist7->SetTitle("D0-D0bar 4 < Pt,Pt' < 6;#Delta#phi (rad); (1/N)dN/#Delta#phi");
   myHist7->SetStats(kFALSE);
   myHist7->GetYaxis()->SetRangeUser(0,0.4);
   myHist7->SetLineColor(kBlack);
   myHist7->SetLineStyle(2);
   myHist7->SetLineWidth(2);
   myHist7->SetMarkerStyle(50);
   myHist7->SetMarkerSize(2);
   myHist7->SetMarkerColor(kBlack);
   myHist7->Draw("hist P");
   fitFcn->Draw("same");
   backFcn->Draw("same");
   signalFcn->Draw("same");
   line->Draw("same");
   c7->BuildLegend();
   c8->cd();
   myHist8->SetTitle("D0-D0bar 6 < Pt,Pt' < 12;#Delta#phi (rad); (1/N)dN/#Delta#phi");
   myHist8->SetStats(kFALSE);
   myHist8->GetYaxis()->SetRangeUser(0,0.1);
   myHist8->SetMarkerStyle(50);
   myHist8->SetMarkerSize(2);
   myHist8->SetMarkerColor(kBlack);
   myHist8->Draw("hist P");
   fitFcn3->Draw("same");
   backFcn3->Draw("same");
   signalFcn3->Draw("same");
   line3->Draw();
   c8->BuildLegend();*/
   auto *c1 = new TCanvas("c1","multipads",1100,300);
   c1->SetLeftMargin(0.20);
   c1->SetRightMargin(0.08);
   c1->SetTopMargin(0.08);
   c1->SetBottomMargin(0.17);
   c1->Divide(3,1,0,0);
  
   c1->cd(1);
   double inthist1 = myHist->Integral();
   myHist->GetXaxis()->CenterTitle(true);
   myHist->GetYaxis()->CenterTitle(true);
   myHist->GetXaxis()->SetTitleSize(0.08);
   myHist->GetXaxis()->SetLabelSize(0.08);
   myHist->GetYaxis()->SetTitleSize(0.08);
   myHist->GetYaxis()->SetLabelSize(0.08);
   myHist->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHist->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHist->SetTitle(";#Delta#phi (rad);(1/#it{N_{tot}})d#it{N_{tot}} /d#Delta#phi");
   myHist->Scale(1/inthist1);
   myHist2->SetTitle("Pair creation;dphi (rad); # D0D0bar pairs");
   myHist3->SetTitle("Gluon splitting;dphi (rad); # D0D0bar pairs");
   myHist4->SetTitle("Flavour excitation;dphi (rad); # D0D0bar pairs");
   myHist->GetYaxis()->SetRangeUser(0.0001,0.07399);
   myHist->SetStats(kFALSE);
   myHist->SetLineWidth(1);
   myHist->SetMarkerStyle(1);
   myHist->SetMarkerSize(1);
   myHist->Draw("E1");
   myHist3->Scale(1/inthist1);
   myHist3->SetLineColor(3);
   myHist3->SetLineStyle(1);
   myHist3->SetLineWidth(2);
   myHist3->SetMarkerStyle(1);
   myHist3->SetMarkerSize(1);
   myHist3->SetMarkerColor(3);
   myHist3->Draw("same E1");
   myHist2->Scale(1/inthist1);
   myHist2->SetLineColor(2);
   myHist2->SetLineStyle(1);
   myHist2->SetLineWidth(2);
   myHist2->SetMarkerStyle(1);
   myHist2->SetMarkerSize(1);
   myHist2->SetMarkerColor(2);
   myHist2->Draw("same E1");
   myHist4->Scale(1/inthist1);
   myHist4->SetLineColor(4);
   myHist4->SetLineStyle(1);
   myHist4->SetLineWidth(2);
   myHist4->SetMarkerStyle(1);
   myHist4->SetMarkerSize(1);
   myHist4->SetMarkerColor(4);
   myHist4->Draw("same E1");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   c1->cd(1)->BuildLegend();
   auto text2 = new TLatex(0,0.015,"#splitline{#splitline{PYTHIA8, Pb-Pb}{#sqrt{s_{NN}} = 5.02 TeV}}{0-80%}");
   text2->Draw();
   
   c1->cd(2);
   double inthiste1 = myHiste->Integral();
   myHiste->GetXaxis()->CenterTitle(true);
   myHiste->GetYaxis()->CenterTitle(true);
   myHiste->GetXaxis()->SetTitleSize(0.08);
   myHiste->GetXaxis()->SetLabelSize(0.08);
   myHiste->GetYaxis()->SetTitleSize(0.08);
   myHiste->GetYaxis()->SetLabelSize(0.08);
   myHiste->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHiste->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHiste->SetTitle(";#Delta#phi (rad);(1/#it{N_{tot}})d#it{N_{tot}} /d#Delta#phi");
   myHiste->Scale(1/inthiste1);
   myHiste2->SetTitle("Pair creation;dphi (rad); # D0D0bar pairs");
   myHiste3->SetTitle("Gluon splitting;dphi (rad); # D0D0bar pairs");
   myHiste4->SetTitle("Flavour excitation;dphi (rad); # D0D0bar pairs");
   myHiste->GetYaxis()->SetRangeUser(0.0000001,0.07399);
   myHiste->SetStats(kFALSE);
   myHiste->SetLineWidth(1);
   myHiste->SetMarkerStyle(1);
   myHiste->SetMarkerSize(1);
   myHiste->Draw("E1");
   myHiste3->Scale(1/inthiste1);
   myHiste3->SetLineColor(3);
   myHiste3->SetLineStyle(1);
   myHiste3->SetLineWidth(2);
   myHiste3->SetMarkerStyle(1);
   myHiste3->SetMarkerSize(1);
   myHiste3->SetMarkerColor(3);
   myHiste3->Draw("same E1");
   myHiste2->Scale(1/inthiste1);
   myHiste2->SetLineColor(2);
   myHiste2->SetLineStyle(1);
   myHiste2->SetLineWidth(2);
   myHiste2->SetMarkerStyle(1);
   myHiste2->SetMarkerSize(1);
   myHiste2->SetMarkerColor(2);
   myHiste2->Draw("same E1");
   myHiste4->Scale(1/inthiste1);
   myHiste4->SetLineColor(4);
   myHiste4->SetLineStyle(1);
   myHiste4->SetLineWidth(2);
   myHiste4->SetMarkerStyle(1);
   myHiste4->SetMarkerSize(1);
   myHiste4->SetMarkerColor(4);
   myHiste4->Draw("same E1");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   //c1->cd(2)->BuildLegend();
   auto texte2 = new TLatex(0,0.015,"|#eta| < 0.9");
   texte2->Draw();

   c1->cd(3);
   double inthiste1b = myHisteb->Integral();
   myHisteb->GetXaxis()->CenterTitle(true);
   myHisteb->GetYaxis()->CenterTitle(true);
   myHisteb->GetXaxis()->SetTitleSize(0.08);
   myHisteb->GetXaxis()->SetLabelSize(0.08);
   myHisteb->GetYaxis()->SetTitleSize(0.08);
   myHisteb->GetYaxis()->SetLabelSize(0.08);
   myHisteb->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHisteb->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
   myHisteb->SetTitle(";#Delta#phi (rad);(1/#it{N_{tot}})d#it{N_{tot}} /d#Delta#phi");
   myHisteb->Scale(1/inthiste1b);
   myHiste2b->SetTitle("Pair creation;dphi (rad); # D0D0bar pairs");
   myHiste3b->SetTitle("Gluon splitting;dphi (rad); # D0D0bar pairs");
   myHiste4b->SetTitle("Flavour excitation;dphi (rad); # D0D0bar pairs");
   myHisteb->GetYaxis()->SetRangeUser(0.000001,0.07399);
   myHisteb->SetStats(kFALSE);
   myHisteb->SetLineWidth(1);
   myHisteb->SetMarkerStyle(1);
   myHisteb->SetMarkerSize(1);
   myHisteb->Draw("E1");
   myHiste3b->Scale(1/inthiste1b);
   myHiste3b->SetLineColor(3);
   myHiste3b->SetLineStyle(1);
   myHiste3b->SetLineWidth(2);
   myHiste3b->SetMarkerStyle(1);
   myHiste3b->SetMarkerSize(1);
   myHiste3b->SetMarkerColor(3);
   myHiste3b->Draw("same E1");
   myHiste2b->Scale(1/inthiste1b);
   myHiste2b->SetLineColor(2);
   myHiste2b->SetLineStyle(1);
   myHiste2b->SetLineWidth(2);
   myHiste2b->SetMarkerStyle(1);
   myHiste2b->SetMarkerSize(1);
   myHiste2b->SetMarkerColor(2);
   myHiste2b->Draw("same E1");
   myHiste4b->Scale(1/inthiste1b);
   myHiste4b->SetLineColor(4);
   myHiste4b->SetLineStyle(1);
   myHiste4b->SetLineWidth(2);
   myHiste4b->SetMarkerStyle(1);
   myHiste4b->SetMarkerSize(1);
   myHiste4b->SetMarkerColor(4);
   myHiste4b->Draw("same E1");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   //c1->cd(2)->BuildLegend();
   auto texte2b = new TLatex(0,0.015,"|#eta| < 0.5");
   texte2b->Draw();

   auto *c12 = new TCanvas("c12","multipads",1100,300);
   c12->SetLeftMargin(0.20);
   c12->SetRightMargin(0.08);
   c12->SetTopMargin(0.08);
   c12->SetBottomMargin(0.17);
   c12->Divide(3,1,0,0);
   
   c12->cd(1);
   double inthiste3v2 = myHiste3v2->Integral();
   myHiste3v2->GetXaxis()->CenterTitle(true);
   myHiste3v2->GetYaxis()->CenterTitle(true);
   myHiste3v2->GetXaxis()->SetTitleSize(0.08);
   myHiste3v2->GetXaxis()->SetLabelSize(0.08);
   myHiste3v2->GetYaxis()->SetTitleSize(0.08);
   myHiste3v2->GetYaxis()->SetLabelSize(0.08);
   myHiste3v2->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHiste3v2->GetYaxis()->SetNdivisions(5,10,0,kTRUE);
   myHiste3v2->SetTitle(";#Delta#phi (rad);(1/#it{N_{#alpha}})d#it{N_{#alpha}} /d#Delta#phi");
   myHiste3v2->Scale(1/inthiste3v2);
   myHiste3v2->GetYaxis()->SetRangeUser(0.0000001,0.13995);
   myHiste3v2->SetStats(kFALSE);
   myHiste3v2->SetLineWidth(1);
   myHiste3v2->SetLineColor(3);
   myHiste3v2->SetMarkerStyle(1);
   myHiste3v2->SetMarkerSize(1);
   myHiste3v2->Draw("E1");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   auto legend = new TLegend(0.1,0.7,0.48,0.9);
   //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   legend->AddEntry(myHiste3v2,"Gluon splitting","l");
   legend->AddEntry(myHiste2v2,"Pair creation","l");
   legend->AddEntry(myHiste4v2,"Flavour excitation","l");
   legend->Draw();
   auto texte4 = new TLatex(0,0.015,"#splitline{#splitline{PYTHIA8, Pb-Pb}{#sqrt{s_{NN}} = 5.02 TeV}}{0-80%, |#eta| < 0.9}");
   texte4->Draw();
   
   c12->cd(2);
   double inthiste2v2 = myHiste2v2->Integral();
   myHiste2v2->GetXaxis()->CenterTitle(true);
   myHiste2v2->GetYaxis()->CenterTitle(true);
   myHiste2v2->GetXaxis()->SetTitleSize(0.08);
   myHiste2v2->GetXaxis()->SetLabelSize(0.08);
   myHiste2v2->GetYaxis()->SetTitleSize(0.08);
   myHiste2v2->GetYaxis()->SetLabelSize(0.08);
   myHiste2v2->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHiste2v2->GetYaxis()->SetNdivisions(5,10,0,kTRUE);
   myHiste2v2->SetTitle(";#Delta#phi (rad);(1/#it{N_{tot}})d#it{N_{tot}} /d#Delta#phi");
   myHiste2v2->Scale(1/inthiste2v2);
   myHiste2v2->GetYaxis()->SetRangeUser(0.0000001,0.13995);
   myHiste2v2->SetStats(kFALSE);
   myHiste2v2->SetLineColor(2);
   myHiste2v2->SetLineWidth(1);
   myHiste2v2->SetMarkerStyle(1);
   myHiste2v2->SetMarkerSize(1);
   myHiste2v2->Draw("E1");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   /*auto texte5 = new TLatex(0,0.015,"|#eta| < 0.9");
     texte5->Draw();*/

   c12->cd(3);
   double inthiste4v2 = myHiste4v2->Integral();
   myHiste4v2->GetXaxis()->CenterTitle(true);
   myHiste4v2->GetYaxis()->CenterTitle(true);
   myHiste4v2->GetXaxis()->SetTitleSize(0.08);
   myHiste4v2->GetXaxis()->SetLabelSize(0.08);
   myHiste4v2->GetYaxis()->SetTitleSize(0.08);
   myHiste4v2->GetYaxis()->SetLabelSize(0.08);
   myHiste4v2->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
   myHiste4v2->GetYaxis()->SetNdivisions(5,10,0,kTRUE);
   myHiste4v2->SetTitle(";#Delta#phi (rad);(1/#it{N_{tot}})d#it{N_{tot}} /d#Delta#phi");
   myHiste4v2->Scale(1/inthiste4v2);
   myHiste4v2->GetYaxis()->SetRangeUser(0.0000001,0.13995);
   myHiste4v2->SetStats(kFALSE);
   myHiste4v2->SetLineColor(4);
   myHiste4v2->SetLineWidth(1);
   myHiste4v2->SetMarkerStyle(1);
   myHiste4v2->SetMarkerSize(1);
   myHiste4v2->Draw("E1");
   gPad->Modified();
   gPad->SetTickx();
   gPad->SetTicky();
   /*auto texte6 = new TLatex(0,0.015,"|#eta| < 0.9");
     texte6->Draw();*/
   
}
