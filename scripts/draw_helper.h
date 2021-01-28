#pragma once

#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TH1.h>

static int HISTO_INDEX = 0; //current canvas number
static int CANVAS_INDEX = 0; //current canvas number

#include "PhysConst.h"
#include "ScanPoint.h"

//inline TCanvas * get_new_tailed_canvas(std::string title)
inline TCanvas * make_canvas(std::string title)
{
  auto c = new TCanvas;
  const int Nx = 4; //numbe of canvases on x size;
  const int Ny = 3; //number of canvases on y size;
  const int window_title_bar_ysize = 45;
  const int window_border_xsize = 5;
  const int panel_ysize = 40;
  const int display_size_x = 1920*2;
  const int display_size_y = 1080*2;
  const int canvas_width_x = display_size_x/Nx;
  const int canvas_width_y = (display_size_y-panel_ysize)/Ny;
  int canvas_pos_x = (CANVAS_INDEX % Nx)* canvas_width_x;
  int canvas_pos_y = ((CANVAS_INDEX/Nx)%Ny)* canvas_width_y;
  c->SetWindowPosition(canvas_pos_x,canvas_pos_y);
  c->SetWindowSize(canvas_width_x-window_border_xsize,canvas_width_y-window_title_bar_ysize);
  c->SetTitle(title.c_str());
  c->cd();
  CANVAS_INDEX++;
  return c;
}


template<typename ProjX, typename ProjY>
inline TGraph *  draw(const Scan_t & SPL, ProjX projx, ProjY projy) {
  make_canvas("draw");
  TGraphErrors * g = new TGraphErrors;
  for(const auto & sp : SPL) {
    double x = std::invoke(projx, sp).value;
    double y = std::invoke(projy, sp).value;
    double dx = std::invoke(projx, sp).error;
    double dy = std::invoke(projy, sp).error;
    int n = g->GetN();
    g->SetPoint(n, x,y);
    g->SetPointError(n,dx,dy);
  }
  g->SetMarkerStyle(21);
  g->Draw("ap");
  return g;
}


inline TGraph * draw_efficiency(const Scan_t & SR) {
  return draw(SR, [](const auto & sp) { return sp.energy; }, [](const auto & sp) { return sp.tt.efficiency; });
}

inline void draw_lum_per_run2(const std::vector<ScanPoint_t> & SPL)
{
  TMultiGraph * mg = new TMultiGraph;
  int i=0;
  std::vector<TGraphErrors*> g(SPL.size());
  int point=0;
  for (const auto & p : SPL)
  {
    g[point] = new TGraphErrors;
    g[point]->SetMarkerStyle(21);
    g[point]->SetMarkerColor(point+1);
    for(auto & ri : p.runs)
    {
      auto run = ri.first;
      auto lum = ri.second;
      auto chain = get_chain("gg", "gg","g", run, run);
      double Ngg = chain->GetEntries();
      double dNgg = sqrt(Ngg);
      int n = g[point]->GetN();
      g[point]->SetPoint(n,run, Ngg/lum);
      g[point]->SetPointError(n, 0, dNgg/lum);
      i++;
    }
    mg->Add(g[point],"p");
    point++;
  }
  mg->Draw("a");
  mg->GetXaxis()->SetTitle("run");
  mg->GetYaxis()->SetTitle("N_{e^{+}e^{-1} #rightarrow #gamma #gamma} / L_{online}, nb");
}

inline void  draw_lum_per_run(const char * dir, int begin_run, int end_run)
{
  auto * g = new TGraphErrors;
  for ( int run = begin_run; run<=end_run; ++run)
  {
    auto chain = get_chain("gg", "gg","g", run, run);
    double Ngg = chain->GetEntries();
    double dNgg = sqrt(Ngg);
    int n = g->GetN();
    g->SetPoint(n,run, Ngg);
    g->SetPointError(n, 0, dNgg);
  }
  g->Draw("a*");
}

inline void draw_tau_per_run(const std::vector<ScanPoint_t> & SPL, const char * selection = "")
{
  TMultiGraph * mg = new TMultiGraph;
  int i=0;
  std::vector<TGraphErrors*> g(SPL.size());
  int point=0;
  for (const auto & p : SPL)
  {
    g[point] = new TGraphErrors;
    g[point]->SetMarkerStyle(21);
    g[point]->SetMarkerColor(point+1);
    for(auto & ri : p.runs)
    {
      auto run = ri.first;
      auto lum = ri.second;
      auto chain = get_chain("gg", "gg","g", run, run);
      auto tt = get_chain("tt","tt","t",run,run);
      double Ngg = chain->GetEntries();
      double dNgg = sqrt(Ngg);
      double L = Ngg/10*pow(p.energy/MTAU*0.5,2);
      double dL = dNgg/10*pow(p.energy/MTAU*0.5,2);
      double Ntt = tt->GetEntries(selection);
      double dNtt = sqrt(Ntt);
      double sigma = Ntt/L;
      double dsigma = Ntt/L* sqrt(1./Ntt + 1.0/Ngg); 
      int n = g[point]->GetN();
      g[point]->SetPoint(n,run, sigma);
      g[point]->SetPointError(n, 0, dsigma);
      i++;
    }
    mg->Add(g[point],"p");
    point++;
  }
  mg->Draw("a");
  mg->GetXaxis()->SetTitle("run");
  mg->GetYaxis()->SetTitle("N_{#tau #tau} / L_{online}, nb");
}

inline TGraphErrors * draw_result(const char * selection, const std::vector<ScanPoint_t> & Points)
{
  TGraphErrors * g = new TGraphErrors;
  long totalNtt=0;
  for(int i=0; i<Points.size();++i)
  {
    double Ntt = Points[i].tt.tree->GetEntries(selection);
    double Ngg = Points[i].gg.tree->GetEntries(selection);
    double xs  = 0;
    double dxs = 0;
    totalNtt += Ntt;
    if(Ngg != 0 )
    {
      xs = Ntt/Ngg;
      dxs = sqrt( Ntt/(Ngg*Ngg) + pow(Ntt/(Ngg*Ngg), 2.0)*Ngg );
    }
    g->SetPoint(i, 0.5*Points[i].energy-MTAU, xs);
    g->SetPointError(i, 0.5*Points[i].energy.error, dxs);
    std::cout << i << " " << 0.5*Points[i].energy-MTAU << "  " << Ngg << "  " << Ntt << "   " <<  xs << std::endl;
  }
  std::cout << "Total number of tau-tau candidates:" << totalNtt << std::endl;
  g->SetMarkerStyle(21);
  g->Draw("ap");

  TCanvas * cacop = new TCanvas("acop","acop");
  cacop->Divide(2,3);
  for( int i=0;i<Points.size();++i) 
  { 
    cacop->cd(i+1);
    Points[i].tt.Draw("acop",selection);
  }
  TCanvas * cptem = new TCanvas("ptem","ptem");
  cptem->Divide(2,3);
  for( int i=0;i<Points.size();++i) 
  { 
    cptem->cd(i+1);
    Points[i].tt.Draw("ptem",selection);
  }
  TCanvas * cp = new TCanvas("p","p");
  cp->Divide(2,3);
  for( int i=0;i<Points.size();++i) 
  { 
    cp->cd(i+1);
    Points[i].tt.tree->Draw("p",selection);
  }
  return g;
}
