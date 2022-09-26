#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>

#include <TH1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TTree.h>
#include <TTreeCache.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TSystem.h>

#define DEBUG (0)
#include "PlotTools.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"

using namespace std;

void printRunTime(TStopwatch timer_)
{
  Double_t cpuTime = timer_.CpuTime();
  Double_t realTime = timer_.RealTime();

  cout << endl;
  cout << "************************************************" << endl;
  cout << "Total real time: " << realTime << " (seconds)" << endl;
  cout << "Total CPU time:  " << cpuTime << " (seconds)" << endl;
  cout << "  CPU time / real time = " << cpuTime / realTime << endl;
  cout << "************************************************" << endl;
}

static inline void printMemory( TString tab = "" )
{
  ifstream proc_status("/proc/self/status");
  string buffer;
  while (proc_status.peek() != EOF) {
    getline(proc_status, buffer);
    TString str = buffer;
    if(str.Contains("RSS")) {
      cout << tab << str << endl;
      break;
    }
  }
}

static inline void loadBar(int x, int num, int r, int w)
{
  // Only update r times.
  if( x == num )
    cout << endl;

  if ( x % (num/r +1) != 0 ) return;

  // Calculuate the ratio of complete-to-incomplete.
  float ratio = x/(float)num;
  int   c     = ratio * w;

  // Show the percentage complete.
  printf("%3d%% [", (int)(ratio*100) );

  // Show the load bar.
  for (int x=0; x<c; x++) cout << "=";

  for (int x=c; x<w; x++) cout << " ";

  // ANSI Control codes to go back to the
  // previous line and clear it.
  cout << "]\r" << flush;
}


// echo 'gROOT->LoadMacro("drawEfficiencyIsoMu24.C+"); gSystem->Exit(0);' | root -b -l
// rootbq 'drawEfficiencyIsoMu24.C("v03", "DY PU200", "PU200-DYToLL_M50", "L1Tk")'
// rootbq 'drawEfficiencyIsoMu24.C("v03", "DY PU200", "PU200-DYToLL_M50", "")'

void drawEfficiencyIsoMu24(
  // TString ver = "v90", TString SAMPLE = "t#bar{t} PU200", TString tag = "PU200-TTToSemiLep",
  TString ver = "v03", TString SAMPLE = "t#bar{t} PU200", TString tag = "PU200-TTTo2L2Nu",
  TString eff_tag = "Iso / L3", bool isLogy = false  // HERE
) {
  TStopwatch timer_total;
  timer_total.Start();

  gStyle->SetPalette(kRainBow);
  TH1::SetDefaultSumw2(kTRUE);

  TString Dir = "./plots_EffPath_"+ver+"/"+tag+"/";
  if (gSystem->mkdir(Dir,kTRUE) != -1)
    gSystem->mkdir(Dir,kTRUE);


  TString L3_pt_min_str = "p_{T}^{HLT} > 24 GeV";
  TString gen_pt_min_str = "p_{T}^{gen} > 26 GeV";

  // TString fileName = TString::Format("/Users/min/Desktop/HLTUpgrade/Local/test_20210504/Outputs_%s/hist-%s-%s-Eff.root", ver.Data(), ver.Data(), tag.Data());
  TString fileName = TString::Format("./hist-%s-%s-Eff.root", ver.Data(), tag.Data());

  vector<Color_t> v_color = {
    kGray+2,
    kBlack,
    kBlue,
    kRed,
    kMagenta,
    kGreen+2,
  };

  vector<int> v_marker = {
    // 21,
    20,
    22,
    23,
    21,
    24,
    25
  };


  vector<TString> v_var = {"pt", "eta"};  // , "pu"};
  vector< vector<double> > range = {
    {1, 10, 100},  // pt
    {1, -2.4, 2.4},  // eta
    {1, 200, 201}  // PU
  };

  if(tag.Contains("PU140")) {
    range.at(2).at(1) = 140;
    range.at(2).at(2) = 141;
  }

  int n_pt_bins = 13-1;
  double pt_bins[13] = {
    0, 10, 20, 22, 23, 24, 26, 30, 40, 60, 100, 200, 1000
  };

  int n_pt_bins_L122 = 25-1;
  double pt_bins_L122[25] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 18, 20, 21, 22, 23, 26, 30, 40, 60, 100, 200, 1000
  };

  int n_pt_bins_50 = 14-1;
  double pt_bins_50[14] = {
    0, 40, 45, 47, 48, 49, 50, 51, 52, 55, 60, 100, 200, 1000
  };

  // int n_eta_bins = 15-1;
  // double eta_bins[15] = {
  //   -2.4, -2.1, -1.6, -1.2, -0.9,
  //   -0.3, -0.2,  0.0,  0.2,  0.3,
  //    0.9,  1.2,  1.6,  2.1,  2.4
  // };

  int n_eta_bins = 21-1;
  double eta_bins[21] = {
    -2.4, -2.1, -1.9, -1.7, -1.6, -1.5, -1.2, -0.9,
    -0.3, -0.2,  0.0,  0.2,  0.3,
     0.9,  1.2,  1.5, 1.6, 1.7, 1.9, 2.1,  2.4
  };

  vector<TString> types = {
    "Eff/num_Eff_hltL1TkSingleMuFiltered22_genpt26",
    "Eff/num_Eff_hltL3fL1TkSingleMu22L3Filtered24Q_genpt26",
    "Eff/num_Eff_hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41_genpt26",
    "Eff/num_Eff_hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40_genpt26",
    "Eff/num_Eff_hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70_genpt26",
    "Eff/num_Eff_hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk_genpt26"
  };

  vector<TString> types_den = {
    "Eff/den_Eff_hltL1TkSingleMuFiltered22_genpt26",
    "Eff/den_Eff_hltL3fL1TkSingleMu22L3Filtered24Q_genpt26",
    "Eff/den_Eff_hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41_genpt26",
    "Eff/den_Eff_hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40_genpt26",
    "Eff/den_Eff_hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70_genpt26",
    "Eff/den_Eff_hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk_genpt26"
  };

  vector<TString> types_str = {
    "Single L1TkMuon with p_{T} > 22 GeV",  // + quality cuts",
    "Nonisolated single muon trigger with p_{T} > 24 GeV",
    "ECAL isolation filter",
    "HCAL isolation filter",
    "HGCAL isolation filter",
    "Tracker isolation filter",
    // "Isolated single muon trigger with p_{T} > 24 GeV"
  };

  for(int ivar=0; ivar<(int)v_var.size(); ++ivar) {

    double xmin = range[ivar][1];
    double xmax = range[ivar][2];
    double ymin = 0.0;
    double ymax = 1.5;

    if(!v_var[ivar].Contains("pt")) {
      ymin = 0.5;
      ymax = 1.35;
    }
    ymin = 0.0;
    ymax = 1.55;

    TString canvasName = TString::Format("EfficiencyIsoMu24_%s_%s", tag.Data(), v_var[ivar].Data() );
    canvasName.ReplaceAll(".","p").ReplaceAll("-","_");
    TCanvas *c;
    SetCanvas_Square( c, canvasName, kFALSE, kFALSE, 900, 900 );
    c->cd();
    if(isLogy)
      c->SetLogy();

    TLegend *legend;
    SetLegend( legend, 0.15, 0.66, 0.90, 0.83, -1);

    bool isFirst = true;
    for(int i = 0; i<(int)types.size(); ++i) {

      TString the_type_num = types[i];
      TString the_type_den = types_den[i];
      TString the_type_str = types_str[i];

      if(v_var[ivar] == "eta" && the_type_num.Contains("L3pt50"))
        continue;

      TString titleX = GetTitleX(v_var[ivar]+"_gen");  // .ReplaceAll("{gen}", "{HLT}");
      TString titleY = "L1 + HLT efficiency";

      TString den_name = TString::Format("%s_%s", the_type_den.Data(), v_var[ivar].Data() );
      TString num_name = TString::Format("%s_%s", the_type_num.Data(), v_var[ivar].Data() );

      if(v_var[ivar] == "pt") {
        den_name = den_name.ReplaceAll("genpt26", "genpt0");
        num_name = num_name.ReplaceAll("genpt26", "genpt0");
      }

      TH1F* den = Get_Hist( fileName, den_name );
      TH1F* num = Get_Hist( fileName, num_name );

      if(v_var[ivar] == "pt") {
        if(the_type_num.Contains("hltL1TkSingleMuFiltered22")) {
          den = (TH1F*)den->Rebin(n_pt_bins_L122, den_name+"_rb", pt_bins_L122);
          num = (TH1F*)num->Rebin(n_pt_bins_L122, num_name+"_rb", pt_bins_L122);
        }
        else {
          den = (TH1F*)den->Rebin(n_pt_bins, den_name+"_rb", pt_bins);
          num = (TH1F*)num->Rebin(n_pt_bins, num_name+"_rb", pt_bins);
        }
      }
      else if(v_var[ivar] == "eta") {
        den = (TH1F*)den->Rebin(n_eta_bins, den_name+"_rb", eta_bins);
        num = (TH1F*)num->Rebin(n_eta_bins, num_name+"_rb", eta_bins);
      }
      else{
        den = (TH1F*)den->Rebin(range[ivar][0]);
        num = (TH1F*)num->Rebin(range[ivar][0]);
      }

      int nbins = den->GetNbinsX();

      c->cd();

      TGraphAsymmErrors* g = new TGraphAsymmErrors(nbins);
      g->Divide(num, den, "n e0");

      for(int ip=0; ip<nbins; ++ip) {
        if(g->GetPointY(ip) == 0.)  g->SetPointEYhigh(ip, 0.0);
      }

      g->SetTitle("");
      g->SetMarkerSize(1.5);
      g->SetMarkerStyle(v_marker[i]);
      g->SetMarkerColor(v_color[i]);
      g->SetLineColor(  v_color[i]);
      g->SetLineWidth(1);

      g->GetXaxis()->SetLimits( xmin, xmax );
      g->GetXaxis()->SetRangeUser( xmin, xmax );
      g->GetYaxis()->SetRangeUser( ymin, ymax );

      SetAxis_SinglePad( g->GetXaxis(), g->GetYaxis(), titleX, titleY );

      if(isFirst) {
        g->Draw("APE");
        isFirst = false;
      }
      else {
        g->Draw("PE same");
      }

      if( v_var[ivar] == "pu" ) {
        int thep = tag.Contains("PU200") ? 200 : 140;
        legend->AddEntry( g, TString::Format("%s  %.3f", the_type_str.Data(), g->GetPointY(thep)), "lep" );
      }
      else
        legend->AddEntry( g, TString::Format("%s", the_type_str.Data()), "lep" );
    }

    legend->Draw();

    TLatex latex;
    // Latex_Simulation_14TeV( latex );
    latex.DrawLatexNDC( 0.45,0.96, "#scale[0.8]{#font[42]{"+SAMPLE+"}}");
    if(v_var[ivar] != "pt" )
      latex.DrawLatexNDC(0.65, 0.87, "#font[42]{#scale[0.8]{"+gen_pt_min_str+"}}");
    // latex.DrawLatexNDC(0.50, 0.89, "#font[42]{#scale[0.8]{"+eff_tag+"}}");

    TString logy_tag = isLogy ? "_log" : "";
    CMS_lumi(c, 98, 11);
    c->Modified();  c->Update();  c->RedrawAxis();
    gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");
    c->SaveAs(Dir+canvasName+logy_tag+".pdf","pdf");
    c->SaveAs(Dir+canvasName+logy_tag+".C","C");
    c->SaveAs(Dir+canvasName+logy_tag+".root","root");
    gROOT->ProcessLine( "gErrorIgnoreLevel = kPrint;");

    c->Close();
  }


  printRunTime(timer_total);
}
