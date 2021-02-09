#ifndef IBN_TAU_DRAW_H
#define IBN_TAU_DRAW_H
#include <type_traits>
#include <tuple>

#include <TLegend.h>

#include "draw_helper.h"
#include "utils.h"
#include "fold.h"

//draw graph extracted with x and y values extracted by proj
template<typename ProjX, typename ProjY>
inline TGraph *  draw(const Scan_t & SPL, ProjX projx, ProjY projy);


//draw registration efficiency from tt data sample
inline TGraph * draw_efficiency(const Scan_t & SR);

//direct or NORM comparison of the scan point. The Nbin, Min and Max is obliged for right normalization
template<class Proj>
inline std::tuple<TH1*,TH1*> cmp(const ScanPoint_t & sp1, const ScanPoint_t & sp2, Proj proj, std::string var, std::string sel, std::string gopt, int Nbin, double Min, double Max);

//direct or NORM comparison of the two scans. This function draws histograms comparison for each point on separate canvas
template<class Proj>
inline void cmp(const Scan_t & SP1, const Scan_t & SP2, Proj proj, std::string var, std::string sel, std::string gopt, int Nbin, double Min, double Max);

inline TH1 * SelectHistogram(TTree * t, std::string var, std::string sel, int Nbin=0, double min=std::numeric_limits<double>::min(), double max=std::numeric_limits<double>::max());

inline TH1 * SelectHistogram(const DataSample_t & d, std::string var, std::string sel, int Nbin=0, double min=std::numeric_limits<double>::min(), double max=std::numeric_limits<double>::max());

//##################### IMPLEMENTATION #####################################################################


//draw graph extracted with x and y values extracted by proj
template<typename ProjX, typename ProjY>
inline TGraph *  draw(const Scan_t & SPL, ProjX projx, ProjY projy)
{
  make_canvas("draw");
  TGraphErrors * g = new TGraphErrors;
  for(const auto & sp : SPL) {
    auto  X = std::invoke(projx, sp);
    auto  Y = std::invoke(projy, sp);
    double x,y,ex{0},ey{0};
    if constexpr (  std::is_same_v< ibn::valer<double>, decltype(X)> )  {
      x=X.value;
      ex=X.error;
    } else if constexpr (  std::is_arithmetic_v<std::remove_cv_t<decltype(X)> > ) {
      x=X;
      ex=0;
    }    

    if constexpr (  std::is_same_v< ibn::valer<double>, decltype(Y)> )  {
      y=Y.value;
      ey=Y.error;
    } else if constexpr (  std::is_arithmetic_v<std::remove_cv_t<decltype(Y)> > ) {
      y=Y;
      ey=0;
    }    
    int n = g->GetN();
    g->SetPoint(n, x,y);
    g->SetPointError(n,ex,ey);
  }
  g->SetMarkerStyle(21);
  g->Draw("ap");
  return g;
}


//draw registration efficiency from tt data sample
inline TGraph * draw_efficiency(const Scan_t & SR) {
  return draw(SR, [](const auto & sp) { return sp.energy; }, [](const auto & sp) { return sp.tt.efficiency; });
}


inline TH1 * SelectHistogram(TTree * t, std::string var, std::string sel, int Nbin, double min, double max)
{
  if(t==nullptr) return nullptr;
  std::string hist_name; // = "h"+std::to_string(++HISTO_INDEX);
  if(Nbin>0)  {
    if(min > std::numeric_limits<double>::min() ||  max<std::numeric_limits<double>::max()) {
      hist_name = ibn::format(">>h%d(%d,%f,%f)", ++HISTO_INDEX,  Nbin, min,max);
    } else {
      hist_name = ibn::format(">>h%d(%d", ++HISTO_INDEX,  Nbin);
    }
  } else {
      hist_name = ibn::format(">>h%d", ++HISTO_INDEX);
  }
  t->Draw((var+hist_name).c_str(),  sel.c_str(), "goff");
  return t->GetHistogram();
}

inline TH1 * SelectHistogram(const DataSample_t & d, std::string var, std::string sel, int Nbin, double min, double max)
{
  return SelectHistogram(d.tree.get(), var, sel, Nbin,min,max );
}




template<class Proj>
inline std::tuple<TH1*,TH1*> cmp(const ScanPoint_t & sp1, const ScanPoint_t & sp2, Proj proj, std::string var, std::string sel, std::string gopt, int Nbin, double Min, double Max)
{
  auto c  = make_canvas(var + " " + sel);
  const auto & d1   = std::invoke(proj,sp1);
  const auto & d2   = std::invoke(proj,sp2);
  auto h1  =  SelectHistogram(d1.tree.get(), var,sel, Nbin, Min,Max);
  auto h2  =  SelectHistogram(d2.tree.get(), var,sel, Nbin, Min,Max);
  h1->SetLineColor(kBlue);
  h2->SetLineColor(kRed);
  TLegend * leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->AddEntry(h1,sp1.title.c_str(), "lp");
  leg->AddEntry(h2,sp2.title.c_str(), "lp");
  if(gopt=="NORM") {
    h1->DrawNormalized();
    h2->DrawNormalized("same");
  } else  {
    h1->Draw();
    h2->Draw("same");
  }
  leg->Draw();
  c->Modified();
  c->Update();
  return {h1,h2};
};

template<class Proj>
inline void cmp(const Scan_t & SP1, const Scan_t & SP2, Proj proj, std::string var, std::string sel, std::string gopt, int Nbin, double Min, double Max)
{
  for(int i=0;i!=std::min(SP1.size(), SP2.size()); ++i ) {
    cmp(SP1[i],SP2[i],proj, var,sel,gopt, Nbin, Min,Max);
  };
};

template<class Proj, class Norm>
inline  std::tuple<TH1*,TH1*> cmp(const ScanPoint_t & sp1, const ScanPoint_t & sp2, Proj proj, Norm norm, std::string var, std::string sel, std::string gopt, int Nbin, double Min, double Max)
{
  auto c  = make_canvas(var + " " + sel);
  const auto & d1   = std::invoke(proj,sp1);
  const auto & d2   = std::invoke(proj,sp2);
  //double norm1      = std::invoke(norm1,sp1);
  //double norm2      = std::invoke(norm2,sp2);
  auto h1  =  SelectHistogram(d1.tree.get(), var,sel, Nbin, Min,Max);
  auto h2  =  SelectHistogram(d2.tree.get(), var,sel, Nbin, Min,Max);
  h1->SetLineColor(kBlue);
  h2->SetLineColor(kRed);
  TLegend * leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->AddEntry(h1,sp1.title.c_str(), "lp");
  leg->AddEntry(h2,sp2.title.c_str(), "lp");
  h1->DrawNormalized("", std::invoke(norm,sp1));
  h2->DrawNormalized("same", std::invoke(norm,sp1));
  leg->Draw();
  c->Modified();
  c->Update();
  return {h1,h2};
};

template<typename HIS>
inline TH1 * fold(const std::vector<HIS> & hs) {
  if(hs.empty()) return nullptr;
  int Nbin = hs[0]->GetNbinsX();
  double min =  hs[0]->GetBinLowEdge(1);
  double max  = hs[0]->GetBinLowEdge(Nbin)+hs[0]->GetBinWidth(Nbin);
  auto hsum = new TH1F(get_unique_histogram_name().c_str(), "", Nbin, min,max);
  if constexpr (  std::is_same_v<std::unique_ptr<TH1>, HIS> ) {
    for(auto & h : hs) hsum->Add(h.get());
  } else {
    for(auto & h : hs) hsum->Add(h);
  }
  return hsum;
}

template<class Proj,class Norm>
inline std::tuple<TH1*,TH1*> cmp(const Scan_t & SP1, const Scan_t & SP2, Proj proj, Norm norm, std::string var, std::string sel, std::string gopt, int Nbin, double Min, double Max)
{
  if(SP1.empty() || SP2.empty() ) return {nullptr,nullptr};
  auto c  = make_canvas(var + " " + sel);
  std::vector<std::unique_ptr<TH1>> hs1,hs2;
  for(int i=0;i!=std::min(SP1.size(), SP2.size()); ++i ) {
    const auto & d1   = std::invoke(proj,SP1[i]);
    const auto & d2   = std::invoke(proj,SP2[i]);
    auto h1  =  SelectHistogram(d1, var,sel, Nbin, Min,Max);
    auto h2  =  SelectHistogram(d2, var,sel, Nbin, Min,Max);
    if(h1==nullptr || h2==nullptr) continue;
    double N1 = h1->GetEntries();
    double N2 = h2->GetEntries();
    //h1->Scale(std::invoke(norm,SP1[i]));
    h2->Scale(N1/N2);
    hs1.push_back(std::unique_ptr<TH1>(h1));
    hs2.push_back(std::unique_ptr<TH1>(h2));
  };

  auto h1 = fold(hs1);
  auto h2 = fold(hs2);
  //for(auto h : hs1) delete h;
  //for(auto h : hs2) delete h;
  h1->SetLineColor(kBlue);
  h2->SetLineColor(kRed);
  h1->Draw();
  h2->Draw("same");

  TLegend * leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->AddEntry(h1,SP1[0].title.c_str(), "lp");
  leg->AddEntry(h2,SP2[0].title.c_str(), "lp");
  leg->Draw();


  return {h1,h2};
};
#endif
