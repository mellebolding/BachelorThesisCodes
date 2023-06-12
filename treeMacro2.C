#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TCanvas.h"

void treeMacro2() {
  // Create a histogram for the values we read.
  auto myHist = new TH1F("h1","particles",100,-4000,4000);
  auto myHist2 = new TH1F("h2","particles",100,-4000,4000);
  auto myHist3 = new TH1F("h3","particles",100,-4000,4000);
  // Open the file containing the tree.
  auto myFile = TFile::Open("results12.root");
  if (!myFile || myFile->IsZombie()) {
    return;
  }
  
  // Create a TTreeReader for the tree, for instance by passing the
  // TTree's name and the TDirectory / TFile it is in.
  TTreeReader myReader("pytTree", myFile);
  // The branch "px" contains floats; access them as myPx.
  TTreeReaderArray<double> myCh(myReader, "nChar");
  TTreeReaderArray<double> myPar(myReader, "parId");
  TTreeReaderArray<double> myEta(myReader, "particlesSigEta");
  
  // Loop over all entries of the TTree or TChain.
  int m = 0;
  double cross = 0;
  int ntot = 0;
  double parttot = 0;
  double avparttot = 0;
  double chartot = 0;
  double avchartot = 0;
  double ytot = 0;
  int n11 = 0;
  double part11 = 0;
  double avpart11 = 0;
  double char11 = 0;
  double avchar11 = 0;
  int n2 = 0;
  double part2 = 0;
  double avpart2 = 0;
  double char2 = 0;
  double avchar2 = 0;
  int n3 = 0;
  double part3 = 0;
  double avpart3 = 0;
  double char3 = 0;
  double avchar3 = 0;
  int n4 = 0;
  double part4 = 0;
  double avpart4 = 0;
  double char4 = 0;
  double avchar4 = 0;
  int n5 = 0;
  double part5 = 0;
  double avpart5 = 0;
  double char5 = 0;
  double avchar5 = 0;
  int n6 = 0;
  double part6 = 0;
  double avpart6 = 0;
  double char6 = 0;
  double avchar6 = 0;
  int n7 = 0;
  double part7 = 0;
  double avpart7 = 0;
  double char7 = 0;
  double avchar7 = 0;
  int n8 = 0;
  double part8 = 0;
  double avpart8 = 0;
  double char8 = 0;
  double avchar8 = 0;
  int n12 = 0;
  double part12 = 0;
  double avpart12 = 0;
  double char12 = 0;
  double avchar12 = 0;
  int n13 = 0;
  double part13 = 0;
  double avpart13 = 0;
  double char13 = 0;
  double avchar13 = 0;
  int n14 = 0;
  double part14 = 0;
  double avpart14 = 0;
  double char14 = 0;
  double avchar14 = 0;
  int n9 = 0;
  double part9 = 0;
  double avpart9 = 0;
  double char9 = 0;
  double avchar9 = 0;
  int n10 = 0;
  double part10 = 0;
  double avpart10 = 0;
  double char10 = 0;
  double avchar10 = 0;
  double devsqpart11 = 0;
  double devsqpart12 = 0;
  double devsqpart13 = 0;
  double devsqpart14 = 0;
  double devsqpart2 = 0;
  double devsqpart3 = 0;
  double devsqpart4 = 0;
  double devsqpart5 = 0;
  double devsqpart6 = 0;
  double devsqpart7 = 0;
  double devsqpart8 = 0;
  double devsqchar11 = 0;
  double devsqchar12 = 0;
  double devsqchar13 = 0;
  double devsqchar14 = 0;
  double devsqchar2 = 0;
  double devsqchar3 = 0;
  double devsqchar4 = 0;
  double devsqchar5 = 0;
  double devsqchar6 = 0;
  double devsqchar7 = 0;
  double devsqchar8 = 0;
  double testcherr = 0;
  
  
  while (myReader.Next()) {
    
    // Cross section
    cross += myCh[4];

    // Centrality bins
    if (myCh[2] < 2.47) {
      n11 += 1;
      part11 += myCh[1];
      char11 += myCh[0];
      avpart11 = part11/n11;
      devsqpart11 += pow((myCh[1]-avpart11),2.0);
      avchar11 = char11/n11;
      devsqchar11 += pow((myCh[0]-avchar11),2);
      testcherr += pow(myCh[0],2);
      
      for(int i = 0;i<myPar.GetSize();++i) {
	myHist->Fill(myPar[i]);
	myHist3->Fill(myPar[i]);}
      
    }
    if (myCh[2] > 2.47 && myCh[2] < 3.50) {
      n12 += 1;
      part12 += myCh[1];
      char12 += myCh[0];
      avpart12 = part12/n12;
      devsqpart12 += pow((myCh[1]-avpart12),2);
      avchar12 = char12/n12;
      devsqchar12 += pow((myCh[0]-avchar12),2);
      for(int i = 0;i<myPar.GetSize();++i) {
	myHist->Fill(myPar[i]);
	myHist3->Fill(myPar[i]);}
    }
    if (myCh[2] > 3.50 && myCh[2] < 4.28) {
      n13 += 1;
      part13 += myCh[1];
      char13 += myCh[0];
      avpart13 = part13/n13;
      devsqpart13 += pow((myCh[1]-avpart13),2);
      avchar13 = char13/n13;
      devsqchar13 += pow((myCh[0]-avchar13),2);
      for(int i = 0;i<myPar.GetSize();++i) {
	myHist->Fill(myPar[i]);
	myHist3->Fill(myPar[i]);}
    }
    if (myCh[2] > 4.28 && myCh[2] < 4.96) {
      n14 += 1;
      part14 += myCh[1];
      char14 += myCh[0];
      avpart14 = part14/n14;
      devsqpart14 += pow((myCh[1]-avpart14),2);
      avchar14 = char14/n14;
      devsqchar14 += pow((myCh[0]-avchar14),2);
      for(int i = 0;i<myPar.GetSize();++i) {
	myHist->Fill(myPar[i]);
	myHist3->Fill(myPar[i]);}
    }
    if (myCh[2] > 4.96 && myCh[2] < 7.01) {
      n2 += 1;
      part2 += myCh[1];
      char2 += myCh[0];
      avpart2 = part2/n2;
      devsqpart2 += pow((myCh[1]-avpart2),2);
      avchar2 = char2/n2;
      devsqchar2 += pow((myCh[0]-avchar2),2);}
    if (myCh[2] > 7.01 && myCh[2] < 8.59) {
      n3 += 1;
      part3 += myCh[1];
      char3 += myCh[0];
      avpart3 = part3/n3;
      devsqpart3 += pow((myCh[1]-avpart3),2);
      avchar3 = char3/n3;
      devsqchar3 += pow((myCh[0]-avchar3),2);}
    if (myCh[2] > 8.59 && myCh[2] < 9.92) {
      n4 += 1;
      part4 += myCh[1];
      char4 += myCh[0];
      avpart4 = part4/n4;
      devsqpart4 += pow((myCh[1]-avpart4),2);
      avchar4 = char4/n4;
      devsqchar4 += pow((myCh[0]-avchar4),2);}
    if (myCh[2] > 9.92 && myCh[2] < 11.1) {
      n5 += 1;
      part5 += myCh[1];
      char5 += myCh[0];
      avpart5 = part5/n5;
      devsqpart5 += pow((myCh[1]-avpart5),2);
      avchar5 = char5/n5;
      devsqchar5 += pow((myCh[0]-avchar5),2);}
    if (myCh[2] > 11.1 && myCh[2] < 12.1) {
      n6 += 1;
      part6 += myCh[1];
      char6 += myCh[0];
      avpart6 = part6/n6;
      devsqpart6 += pow((myCh[1]-avpart6),2);
      avchar6 = char6/n6;
      devsqchar6 += pow((myCh[0]-avchar6),2);}
    if (myCh[2] > 12.1 && myCh[2] < 13.1) {
      n7 += 1;
      part7 += myCh[1];
      char7 += myCh[0];
      avpart7 = part7/n7;
      devsqpart7 += pow((myCh[1]-avpart7),2);
      avchar7 = char7/n7;
      devsqchar7 += pow((myCh[0]-avchar7),2);}
    if (myCh[2] > 13.1 && myCh[2] < 14.0) {
      n8 += 1;
      part8 += myCh[1];
      char8 += myCh[0];
      avpart8 = part8/n8;
      devsqpart8 += pow((myCh[1]-avpart8),2);
      avchar8 = char8/n8;
      devsqchar8 += pow((myCh[0]-avchar8),2);}
    if (myCh[2] > 14.0 && myCh[2] < 15.0) {
      n9 += 1;
      part9 += myCh[1];
      char9 += myCh[0];
      for(int i = 0;i<myPar.GetSize();++i) {
	myHist2->Fill(myPar[i]);}}
    if (myCh[2] > 15.0 && myCh[2] < 19.6) {
      n10 += 1;
      part10 += myCh[1];
      char10 += myCh[0];
      for(int i = 0;i<myPar.GetSize();++i) {
	myHist2->Fill(myPar[i]);}}
  
  }

  //Error calculations
  
  double testvar = pow((testcherr/n11 - pow(char11/n11,2)),0.5)/pow(n11,0.5);
  double stderr_part11 = pow((devsqpart11/(n11-1)),0.5)/pow(n11,0.5);
  double stderr_part12 = pow((devsqpart12/(n12-1)),0.5)/pow(n12,0.5);
  double stderr_part13 = pow((devsqpart13/(n13-1)),0.5)/pow(n13,0.5);
  double stderr_part14 = pow((devsqpart14/(n14-1)),0.5)/pow(n14,0.5);
  double stderr_part2 = pow((devsqpart2/(n2-1)),0.5)/pow(n2,0.5);
  double stderr_part3 = pow((devsqpart3/(n3-1)),0.5)/pow(n3,0.5);
  double stderr_part4 = pow((devsqpart4/(n4-1)),0.5)/pow(n4,0.5);
  double stderr_part5 = pow((devsqpart5/(n5-1)),0.5)/pow(n5,0.5);
  double stderr_part6 = pow((devsqpart6/(n6-1)),0.5)/pow(n6,0.5);
  double stderr_part7 = pow((devsqpart7/(n7-1)),0.5)/pow(n7,0.5);
  double stderr_part8 = pow((devsqpart8/(n8-1)),0.5)/pow(n8,0.5);
  double stderr_avpa[] = {stderr_part11,stderr_part12,stderr_part13,stderr_part14,stderr_part2,stderr_part3,stderr_part4,stderr_part5,stderr_part6,stderr_part7,stderr_part8};

  double stderr_char11 = pow((devsqchar11/(n11-1)),0.5)/pow(n11,0.5);
  double stderr_char12 = pow((devsqchar12/(n12-1)),0.5)/pow(n12,0.5);
  double stderr_char13 = pow((devsqchar13/(n13-1)),0.5)/pow(n13,0.5);
  double stderr_char14 = pow((devsqchar14/(n14-1)),0.5)/pow(n14,0.5);
  double stderr_char2 = pow((devsqchar2/(n2-1)),0.5)/pow(n2,0.5);
  double stderr_char3 = pow((devsqchar3/(n3-1)),0.5)/pow(n3,0.5);
  double stderr_char4 = pow((devsqchar4/(n4-1)),0.5)/pow(n4,0.5);
  double stderr_char5 = pow((devsqchar5/(n5-1)),0.5)/pow(n5,0.5);
  double stderr_char6 = pow((devsqchar6/(n6-1)),0.5)/pow(n6,0.5);
  double stderr_char7 = pow((devsqchar7/(n7-1)),0.5)/pow(n7,0.5);
  double stderr_char8 = pow((devsqchar8/(n8-1)),0.5)/pow(n8,0.5);
  double stderr_avchar[] = {stderr_char11,stderr_char12,stderr_char13,stderr_char14,stderr_char2,stderr_char3,stderr_char4,stderr_char5,stderr_char6,stderr_char7,stderr_char8};


  double avch[13] = {avchar11,avchar12,avchar13,avchar14,avchar2,avchar3,avchar4,avchar5,avchar6,avchar7,avchar8,avchar9,avchar10};
  double avpa[13] = {avpart11,avpart12,avpart13,avpart14,avpart2,avpart3,avpart4,avpart5,avpart6,avpart7,avpart8,avpart9,avpart10};
  double y[13] = {(2*avchar11)/avpart11,(2*avchar12)/avpart12,(2*avchar13)/avpart13,(2*avchar14)/avpart14,(2*avchar2)/avpart2,(2*avchar3)/avpart3,(2*avchar4)/avpart4,(2*avchar5)/avpart5,(2*avchar6)/avpart6,(2*avchar7)/avpart7,(2*avchar8)/avpart8,(2*avchar9)/avpart9,(2*avchar10)/avpart10};
  
  double stderr_y[11] = {};
  for (int i=0; i < 11; ++i){
    stderr_y[i] = 2*y[i]*pow((pow(stderr_avchar[i]/avch[i],2)+pow(stderr_avpa[i]/avpa[i],2)),0.5);
    std::cout << " y  "<< y[i] << " pm " << stderr_y[i] << std::endl;
    
    }
  double stderry05 = pow((pow(stderr_y[0],2)+pow(stderr_y[1],2)),0.5);
  std::cout << " stderr05 " << stderry05 << std::endl;
  double bins[] = {0.25,0.5,0.75,1,2,3,4,5,6,7,8,9,10};
  double err_bins[13] = {};
    
  ntot = n11+n12+n13+n14+n2+n3+n3+n4+n5+n6+n7+n8+n9+n10;
  parttot = part11+part12+part13+part14+part2+part3+part4+part5+part6+part7+part8+part9+part10;
  chartot = char11+char12+char13+char14+char2+char3+char4+char5+char6+char7+char8+char9+char10;
  double parttot05 = part11+part12;
  double chartot05 = char11+char12;
  avparttot = parttot/ntot;
  avchartot = chartot/ntot;
  ytot = (2*chartot05)/parttot05;
  cross = cross/ntot;
  std::cout << " ntot " << ntot << std::endl;
  std::cout << " cross-section " << cross << std::endl;
  std::cout<< " ytot " << ytot << std::endl;
  std::cout<< " av11 " << avpart11 << " pm " << stderr_part11 << " pm " << pow(n11,-0.5) << std::endl;
  std::cout<< " av12 " << avpart12 << " pm " << stderr_part12 <<std::endl;
  std::cout<< " av13 " << avpart13 << " pm " << stderr_part13 <<std::endl;
  std::cout<< " av14 " << avpart14 << " pm " << stderr_part14 <<std::endl;
  std::cout<< " av2 " << avpart2 << " pm " << stderr_part2 <<std::endl;
  std::cout<< " av3 " << avpart3 << " pm " << stderr_part3 <<std::endl;
  std::cout<< " av4 " << avpart4 << " pm " << stderr_part4 <<std::endl;
  std::cout<< " av5 " << avpart5 << " pm " << stderr_part5 <<std::endl;
  std::cout<< " av6 " << avpart6 << " pm " << stderr_part6 <<std::endl;
  std::cout<< " av7 " << avpart7 << " pm " << stderr_part7 <<std::endl;
  std::cout<< " av8 " << avpart8 <<" pm " << stderr_part8 << std::endl;
  
  std::cout<< " avC11 " << avchar11 << " pm " << stderr_char11 << std::endl;
  std::cout<< " avC12 " << avchar12 << " pm " << stderr_char12 << std::endl;
  std::cout<< " avC13 " << avchar13 << " pm " << stderr_char13 <<  std::endl;
  std::cout<< " avC14 " << avchar14 << " pm " << stderr_char14 <<  std::endl;
  std::cout<< " avC2 " << avchar2 << " pm " << stderr_char2 <<  std::endl;
  std::cout<< " avC3 " << avchar3 << " pm " << stderr_char3 <<  std::endl;
  std::cout<< " avC4 " << avchar4 << " pm " << stderr_char4 <<  std::endl;
  std::cout<< " avC5 " << avchar5 << " pm " << stderr_char5 <<  std::endl;
  std::cout<< " avC6 " << avchar6 << " pm " << stderr_char6 <<  std::endl;
  std::cout<< " avC7 " << avchar7 << " pm " << stderr_char7 <<  std::endl;
  std::cout<< " avC8 " << avchar8 << " pm " << stderr_char8 <<  std::endl;
  
  
  double totdata[] = {10.23,9.94,9.64,9.4,8.98,8.37,7.81,7.37,6.84,6.33,5.76};
  double err_tot[] = {0.27, 0.31, 0.29, 0.29, 0.27, 0.26, 0.27, 0.32, 0.34, 0.41, 0.47};
  double npartdata[] = {398,372.2,345.6,320.1,263,188,131,86.3,53.6,30.4,15.6};
  double err_npart[] = {2,3,4,4,4,3,2,1.7,1.2,0.8,0.5};
  double chardata[] = {2035,1850,1666,1505,1180,786,512,318,183,96.3,44.9};
  double err_char[] = {52,55,48,44,31,20,15,12,8,5.8,3.4};


  //plot for comparison with ALICE data
  gStyle->SetPadTickX(1);
  TCanvas *c1 = new TCanvas("c1","multipads",700,700);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  Float_t small = 1e-5;
  c1->Divide(2,1,small,small); 
  TMultiGraph *mg = new TMultiGraph("mg","mg");
  auto g = new TGraphErrors(11,avpa,y,stderr_avpa,stderr_y);
  g->SetLineColor(2);
  g->SetFillStyle(0);
  g->SetMarkerStyle(21);
  g->SetMarkerColor(kRed);
  g->SetTitle("PYTHIA 8.3");
  auto tg = new TGraphErrors(11,npartdata,totdata,err_npart,err_tot);
  tg->SetFillStyle(1);
  tg->SetMarkerStyle(21);
  tg->SetMarkerColor(kBlack);
  tg->SetLineColor(kBlack);
  tg->SetTitle("ALICE Data");
  mg->Add(g);
  mg->Add(tg);
  mg->SetTitle(";#LTN_{part}#GT     ;");
  c1->cd(2);
  gPad->SetLeftMargin(small);
  gPad->SetTickx();
  mg->Draw("APE1");
  gPad->Modified();
  mg->GetYaxis()->SetLimits(0,12);
  mg->SetMinimum(0.);
  mg->SetMaximum(12);
  mg->GetXaxis()->SetNdivisions(5, 5, 0, kTRUE);
  mg->GetYaxis()->SetNdivisions(7, 5, 0, kTRUE);
  mg->GetXaxis()->SetTitleSize(0.05);
  mg->GetYaxis()->SetTitleSize(0.05);
   // Draw the axis of the 2nd TMultiGraph
  TAxis *mgYaxis = mg->GetHistogram()->GetYaxis();
  TAxis *mgXaxis = mg->GetHistogram()->GetXaxis();
  Double_t ymin = mgYaxis->GetXmin();
  Double_t ymax = mgYaxis->GetXmax();
  Double_t xmin = mgXaxis->GetXmin();                       
  Double_t xmax = mgXaxis->GetXmax();
  mgYaxis->SetLabelSize(0);
  TGaxis *yaxis = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,510,"+L");     
  yaxis->SetLabelColor(kBlack);
  yaxis->Draw();  
  gStyle->SetPadTickY(1); 
  auto text = new TLatex(8.2,7.5,"#bf{|#eta|} < 0.5"); 
  text->Draw();


  
  c1->cd(1)->SetLogx();
  int n = 4;
  double x[] = {20,200,2760,5020};
  double y2[] = {2.35,4.16,8.28,10.213};
  double ydata[] = {1.9,3.8,8.2,10.1};
  double err_x[] = {0,0,0,0};
  double err_ydata[] = {0,0.1,0.1,0.1};
  double xsim_err[4] = {};
  double ysim_err[] = {0.03,0.03, 0.03, stderry05}; // uncedrtainty on mean: sig/sqrt(n)

  TMultiGraph *ng = new TMultiGraph();
  auto g2 = new TGraphErrors(n,x,y2,xsim_err,ysim_err);
  g2->SetLineColor(2);
  g2->SetMarkerStyle(21);
  g2->SetMarkerSize(1.2);
  g2->SetMarkerColor(kRed);
  g2->SetTitle("PYTHIA 8.3");
  auto tg2 = new TGraphErrors(n,x,ydata,err_x,err_ydata);
  ng->SetTitle(";#sqrt{s_{NN}} (GeV)  ;#frac{2}{#LTN_{part}#GT}#LTdN_{ch}/d#eta#GT");
  tg2->SetFillStyle(1);
  tg2->SetTitle("ALICE Data");
  tg2->SetMarkerStyle(21);
  tg2->SetMarkerSize(1.2);
  tg2->SetMarkerColor(kBlack);
  tg2->SetLineColor(kBlack);
  tg2->SetMinimum(0.001);
  
  ng->Add(g2);
  ng->Add(tg2);
  c1->cd(1);
  gStyle->SetPadTickY(1);
  gPad->SetRightMargin(small);
  ng->Draw("APE1");
  gPad->Modified();
  
  ng->GetXaxis()->SetLimits(9.9,9999);
  ng->GetXaxis()->SetTitleSize(0.05);
  ng->GetYaxis()->SetTitleSize(0.05);
  ng->SetMinimum(0.);
  ng->SetMaximum(12);
  c1->cd(1)->BuildLegend();
  auto text2 = new TLatex(8.2,7.5,"#splitline{#bf{|#eta|} < 0.5}{0-5% centrality}"); 
  text2->Draw();
  
}
