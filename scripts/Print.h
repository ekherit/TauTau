/*
 * =====================================================================================
 *
 *       Filename:  Print.h
 *
 *    Description: Some printer functions 
 *
 *        Version:  1.0
 *        Created:  30.12.2020 11:43:43
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */

#ifndef IBN_TAU_PRINT_H
#define IBN_TAU_PRINT_H
#include <string>
#include <vector>
#include <iostream>
#include <regex>
#include <list>

#include "PhysConst.h"
#include "utils.h"
#include "ScanPoint.h"
#include "Selection.h"

struct PrintConfig_t
{
  int title_width=15;
  int data_width=20;
  int total_width=10;
  int vline_width=6;
  std::string thline=""; //top horizontal line
  std::string bhline=""; //bottom horizontal line
  PrintConfig_t() = default;
  void hline(int N, std::string symb="─", std::string name="")
  {
    int hline_width = title_width + total_width + vline_width + N*data_width+3;
    if(name=="") {
      //int w = width == 0 ? hline_width : width;
      for(int i=0;i<hline_width;++i) std::cout << symb;
    }
    else {
      int f = (hline_width - 2 - name.length())/2;
      for(int  i =0; i< f; ++i) std::cout << symb;
      std::cout << " " << name << " ";
      for(int  i =0; i< f; ++i) std::cout << symb;
    }
    std::cout << std::endl;
  }
};

class printer
{
  struct print_cfg_t
  {
    std::string head;
    std::string fmt;
  };
  std::vector<print_cfg_t> formats;
  std::string result;
  std::string Head;
  char buf[65535];
  int idx;
  public:
  printer(void) {idx=0;}
  printer & add(std::string h, std::string f) 
  {
    formats.push_back({h,f});
    static const std::regex re(R"(%([^\d]?\d+)(?:\.\d+)?[[:alpha:]])");
    std::smatch sm;
    int size;
    std::string head_format;
    if(std::regex_search(f,sm,re)) head_format = "%"+std::string(sm[1])+"s";
    else std::cerr << "ERROR: bad print format for printer" << std::endl;
    sprintf(buf,head_format.c_str(),h.c_str());
    Head+=buf;
    return *this;
  };
  printer & add(std::string f) 
  {
    return add("",f);
  };
  printer & operator()(std::string h, std::string f) { return add(h,f); }
  template <typename Data>
  printer & operator%(const Data  & d) 
  {
    if(idx==0) result="";
    sprintf(buf,formats[idx].fmt.c_str(),d);
    result+=buf;
    ++idx%=formats.size();
    return *this;
  }
  void print(std::ostream & os)
  {
    os << result;
  }
  const std::string  & head(void)  const { return Head; }
};

std::ostream & operator<<(std::ostream & os, printer & p)
{
  p.print(os);
  return os;
};

std::string get_run_formula(const std::list<int> & runs)
{
  if(runs.empty()) return "";
  std::string formula;
  for ( auto r = std::begin(runs); r!=std::end(runs); ++r)
  {
    if( r == std::begin(runs) ) formula = std::to_string(*r); //start sequence
    else
    {
      if ( *r == *std::prev(r) + 1 )  //different in run number is 1
      {
        if( r == std::prev(std::end(runs)) || *std::next(r) != *r + 1 )  //this is last item or next item is not next run and we need to finish sequency
        {
          formula +="-"+std::to_string(*r); 
        }
        continue;
      }
      else
      {
        formula+=" "+std::to_string(*r);
      }
    }
  }
  return formula;
}

//inline void print(const std::vector<ScanPoint_t> & SPL, std::vector<const DataSample_t & (*)(const ScanPoint_t &)> extra_fields={})
//{
//  int point_number = 1;
//  printer pr;
//  //pr.add("#point",     "%-7d");
//  pr.add("#name",       "%-10s");
//  pr.add("W,GeV",      "%15.6f");
//  pr.add("dW,GeV",     "%15.6f");
////  pr.add("E-Mtau,MeV", "%15.3f");
//  pr.add("Sw,GeV",     "%15.6f");
//  pr.add("dSw,GeV",    "%15.6f");
//  pr.add("L,pb^-1",    "%15.6f");
//  pr.add("dL,pb^-1",   "%15.6f");
//  pr("  run list", "  %-10s");
//  std::cout << pr.head() << std::endl;
//  for (const auto & p : SPL) 
//    std::cout <<  pr 
//   //   % point_number++ 
//      % p.title.c_str() 
//      % p.energy.value % p.energy.error 
////      % ((0.5*p.energy.value-MTAU)/MeV) 
//      % p.energy_spread.value
//      % p.energy_spread.error
//      % p.luminosity.value 
//      % p.luminosity.error
//      % get_run_formula(p.run_list).c_str() << std::endl;
//}
//
//

template<typename ... Projs>
inline void print(const std::vector<ScanPoint_t> & SPL, Projs... ps)
{
  std::cout << ibn::mformat("15", "#name","W,GeV", "dW,GeV", "Sw,GeV", "dSw,GeV", "L,pb^-1", "dL,pb^-1");
  for(size_t i=0;i!=sizeof...(ps);++i) {
    std::cout << ibn::format("%15s","");
  }
  std::cout << "\n";
  for (const auto & sp : SPL) {
    std::cout << ibn::format("%15s", sp.title.c_str());
    std::cout << ibn::mformat("15.6", 
        sp.energy.value, 
        sp.energy.error, 
        sp.energy_spread.value, 
        sp.energy_spread.error, 
        sp.luminosity.value, 
        sp.luminosity.error
        );
    sp.print_column(ps...);
    std::cout << ibn::format(" %20s", get_run_formula(sp.run_list).c_str());
    std::cout << "\n";
  }
}

SelectionResult_t  fold(const std::vector<SelectionResult_t>  & SR, std::string name = "all")
{
  SelectionResult_t result;
  result.title = name;
  if(SR.empty()) return result;
  result.resize(SR[0].size());
  result.Ntt=0;
  for(int i=0;i<result.size();++i) 
  {
    auto & rp   = result[i];
    rp.tt.efficiency.value = 0;
    rp.tt.efficiency.error = 0;
    double eps_error2=0;
    for(const auto & r : SR )  
    {
      auto & p    = r[i];
      result.Ntt += p.tt.N;
      rp.tt.efficiency += p.tt.efficiency;
      //rp.eps     += p.eps;
      //eps_error2 += p.eps.error*p.eps.error;
      rp.tt.N      += p.tt.N;
      rp.gg      = p.gg;
      rp.bb      = p.bb;
      rp.energy  = p.energy;
      rp.luminosity = p.luminosity;
    }
    //rp.eps_error = sqrt(eps_error2);
  }
  auto & reference_point = * find_minimum(result, [](const auto & p) { return  std::abs(p.energy*0.5-MTAU+0.5); } );
  //normalize efficiency correction to threshold efficiency
  for(size_t i=0;i<result.size();++i)
  {
    auto & rp   = result[i];
    rp.tt.effcor = rp.tt.efficiency / reference_point.tt.efficiency.value;
    //rp.effcor.error = rp.eps.error / reference_point.eps;
  }
  reference_point.tt.effcor.error = 0;
  return result;
};

SelectionResult_t  fold(std::vector<std::vector<SelectionResult_t>> & SR, std::string name = "all")
{
  std::vector<SelectionResult_t> VSR;
  for(const auto & vsr : SR) VSR.push_back(fold(vsr));
  return fold(VSR);
};

SelectionResult_t  fold(const std::vector<SelectionResult_t>  & s1, const std::vector<SelectionResult_t>  & s2)
{
  std::vector<SelectionResult_t> VSR;
  VSR.push_back(fold(s1));
  VSR.push_back(fold(s2));
  return fold(VSR);
};
void print(const SelectionResult_t & sr, int opt=1 , int first_column_width=10, int last_column_width=5, int  column_width= 8, int vline_width=6)
{
  int N = sr.size(); //number of points
  int hline_width = first_column_width + last_column_width + vline_width + N*column_width+3;
  auto hline = [&hline_width](std::string  symb="─", int width = 0) 
  { 
    int w = width == 0 ? hline_width : width;
    for(int i=0;i<w;++i)
    {
      std::cout << symb;
    }
    std::cout << std::endl;
  };
  auto vline = [&vline_width](char symb='-', int width = 0) { std::cout << std::setw(width == 0? (vline_width-1) : (width-1)) << symb; };

  if(opt<0) hline(); //preparing for total
  switch (opt)
  {
    case 0:
      hline("━");
      std::cout << std::setw(first_column_width) << "CHNNL/PNT";
      for(int i=0;i<N;++i) std::cout << std::setw(column_width)  << sr[i].name;
      std::cout << std::setw(vline_width) << " │ " << std::setw(last_column_width) << "TOTAL";
      std::cout << std::endl;
      hline();
      std::cout << std::setw(first_column_width+count_utf8_extra_byte(std::string("ΔE,MeV"))) << "ΔE,MeV";
      for(int i=0;i<N;++i) std::cout << std::setw(column_width) << (sr[i].energy*0.5-MTAU)*1e3;
      std::cout << std::setw(vline_width) << " │ ";
      std::cout << std::endl;
      hline();
    default:
      std::cout << std::setw(first_column_width+count_utf8_extra_byte(sr.title)) << sr.title;
      for(int i=0; i < N; i++) std::cout << std::setw(column_width) <<  sr[i].tt.N;
      std::cout << std::setw(vline_width) << " │ " << std::setw(last_column_width) <<  sr.Ntt;
      std::cout << std::endl;
  }
  if(opt<0) hline("━");
};



PrintConfig_t PCFG;

template < typename Item, typename Ftitle , typename Fdata, typename Ftotal>
void print_smth(const std::vector<Item> pts, PrintConfig_t cfg, Ftitle ftitle, Fdata fdata, Ftotal ftotal)
{
  int N = pts.size(); //number of points
  int hline_width = cfg.title_width + cfg.total_width + cfg.vline_width + N*cfg.data_width+3;
  auto hline = [&hline_width](std::string  symb="─", int width = 0) 
  { 
    int w = width == 0 ? hline_width : width;
    for(int i=0;i<w;++i)
    {
      std::cout << symb;
    }
    std::cout << std::endl;
  };
  auto vline_width = cfg.vline_width;
  auto vline = [&vline_width](char symb='-', int width = 0) { std::cout << std::setw(width == 0? (vline_width-1) : (width-1)) << symb; };
  print_utf(cfg.title_width,ftitle());
  for(int i=0;i<N;++i) print_utf(cfg.data_width, fdata(pts[i]));
  std::cout << std::setw(vline_width) << " │ ";
  print_utf(cfg.total_width, ftotal(pts));
  std::cout << std::endl;
};

void print_head(const std::vector<PointSelectionResult_t> pts, PrintConfig_t cfg)
{
  print_smth(pts, cfg, [](){ return "CHNL/PNT"; }, [](auto & p) { return p.name; }, [](auto & s) { return "TOTAL"; } );
};

//void print_eps(const std::vector<PointSelectionResult_t> pts, PrintConfig_t cfg=PrintConfig_t())
//{
//  print_smth(pts, cfg, [](){ return "ε"; }, [](auto & p) { return p.tt.efficiency.value; }, [](auto & s) { return ""; } );
//};

//void print_Ntt(const std::vector<PointSelectionResult_t> pts, PrintConfig_t cfg, std::string title="Ntt")
//{
//  print_smth(pts, cfg, 
//      [&title](){ return title; }, 
//      [](auto & p) { return p.tt.N; }, 
//      [](auto & s) { 
//        long N=0; 
//        for(auto & x: s)  N+=x.tt.N;  
//        return N; 
//        } 
//      );
//};
//
//void print_Ngg(const std::vector<PointSelectionResult_t> pts, PrintConfig_t cfg, std::string title="Ngg")
//{
//  print_smth(pts, cfg, 
//      [&title](){ return title; }, 
//      [](auto & p) { return p.gg.N; }, 
//      [](auto & s) { 
//        long N=0; 
//        for(auto & x: s)  N+=x.gg.N;  
//        return N; 
//        } 
//      );
//};
//
//void print_Nbb(const std::vector<PointSelectionResult_t> pts, PrintConfig_t cfg, std::string title="Nbb")
//{
//  print_smth(pts, cfg, 
//      [&title](){ return title; }, 
//      [](auto & p) { return p.bb.N; }, 
//      [](auto & s) { 
//        long N=0; 
//        for(auto & x: s)  N+=x.bb.N;  
//        return N; 
//        } 
//      );
//};

template<typename Fdata, typename Ftotal>
void print_smth(const std::vector<SelectionResult_t> SR, PrintConfig_t cfg, Fdata fd, Ftotal ft)
{
  if(SR.empty()) return;
  auto hline = [&SR, &cfg](std::string s="─") { cfg.hline(SR[0].size(),s); };
  hline("━");
  print_head(SR[0],cfg);
  hline();
  auto prn = [&](std::string title, auto & sr) {
    print_smth(sr, cfg, 
        [&title](){ return title; }, 
        fd, 
        ft
        );
  };
  for( auto &  sr : SR) 
     prn(sr.title,sr);
  hline();
  auto f = fold(SR);
  prn("all", f);
  hline("━");
};

//void print_Ntt(const std::vector<SelectionResult_t> SR, PrintConfig_t cfg = PCFG)
//{
//  print_smth(SR,cfg, 
//      [](auto & p) { 
//        return p.tt.N; 
//      },
//      [](auto & s) { 
//        double sum=0;
//        for(auto & x: s)  sum+=x.tt.N;  
//        return sum; 
//      } 
//  );
//};
//
//void print_Ngg(const std::vector<SelectionResult_t> SR, PrintConfig_t cfg = PCFG)
//{
//  print_smth(SR,cfg, 
//      [](auto & p) { 
//        return p.gg.N; 
//      },
//      [](auto & s) { 
//        double sum=0;
//        for(auto & x: s)  sum+=x.gg.N;  
//        return sum; 
//      } 
//  );
//};
//
//void print_Nbb(const std::vector<SelectionResult_t> SR, PrintConfig_t cfg = PCFG)
//{
//  print_smth(SR,cfg, 
//      [](auto & p) { 
//        return p.bb.N; 
//      },
//      [](auto & s) { 
//        double sum=0;
//        for(auto & x: s)  sum+=x.bb.N;  
//        return sum; 
//      } 
//  );
//};


void print_efficiency(const std::vector<SelectionResult_t> SR, PrintConfig_t cfg = PCFG)
{
  std::cout << "In print_efficiency " << std::endl;
  if(SR.empty()) return;
  std::string format_str = "%4.3f ± %4.3f";
  cfg.total_width = format_str.length();
  auto hline = [&SR, &cfg](std::string s="─",std::string title="") { cfg.hline(SR[0].size(),s,title); };
  hline("━","REGISTRATION EFFICIENCY");
  print_smth(SR[0], cfg, [](){ return "CHNL/PNT"; }, [](auto & p) { return p.name; }, [](auto & s) { return "AVERAGE"; } );
  hline();
  auto prn = [&](std::string title, auto & sr) {
    print_smth(sr, cfg, 
          [&title](){ return title; }, 
          [&format_str](auto & p) { 
            char buf[1024];
            sprintf(buf,format_str.c_str(), p.tt.efficiency.value*100,p.tt.efficiency.error*100);
            return std::string(buf); 
          },
          [&format_str](auto & P) { 
            double sum=0;
            double error2_sum=0;
            for(auto & p: P)  
            {
              sum+=p.tt.efficiency.value;  
              error2_sum+=p.tt.efficiency.error*p.tt.efficiency.error;
            }
            double average = sum/P.size(); 
            double error = sqrt(error2_sum)/P.size();
            char buf[1024];
            sprintf(buf,format_str.c_str(), average*100,error*100);
            return std::string(buf); 
          } 
        );
  };
  for( auto &  sr : SR) 
     prn(sr.title,sr);
  hline();
  auto f = fold(SR);
  prn("all", f);
  hline("━");
};

void print_effcor(const std::vector<SelectionResult_t> SR, PrintConfig_t cfg = PCFG)
{
  if(SR.empty()) return;
  std::string format_str = " %+5.4f ± %5.4f ";
  cfg.total_width = format_str.length();
  auto hline = [&SR, &cfg](std::string s="─",std::string title="") { cfg.hline(SR[0].size(),s,title); };
  hline("━","CORRECTION TO EFFICIENCY");
  print_smth(SR[0], cfg, [](){ return "CHNL/PNT"; }, [](auto & p) { return p.name; }, [](auto & s) { return "AVERAGE"; } );
  hline();
  auto prn = [&](std::string title, auto & sr) {
    print_smth(sr, cfg, 
          [&title](){ return title; }, 
          [&format_str](auto & p) { 
            char buf[1024];
            sprintf(buf,format_str.c_str(), p.tt.effcor.value,p.tt.effcor.error);
            if(p.tt.effcor.error == 0 && p.tt.effcor==1) 
            {
              return std::string("      1      ");
            }
            return std::string(buf); 
          },
          [&format_str](auto & P) { 
            double sum=0;
            double error2_sum=0;
            for(auto & p: P)  
            {
              sum+=p.tt.effcor.value;  
              error2_sum+=p.tt.effcor.error*p.tt.effcor.error;
            }
            double average = sum/P.size(); 
            double error = sqrt(error2_sum)/P.size();
            char buf[1024];
            sprintf(buf,format_str.c_str(), average,error);
            return std::string(buf); 
          } 
        );
  };
  auto prn2 = [&](std::string title, auto & sr) {
    print_smth(sr, cfg, 
          [&title](){ return title; }, 
          [&format_str](auto & p) { 
            char buf[1024];
            sprintf(buf,format_str.c_str(), (p.tt.effcor.value-1)*100,p.tt.effcor.error*100);
            if(p.tt.effcor.error == 0 && p.tt.effcor==1) 
            {
              return std::string("      0      ");
            }
            return std::string(buf); 
          },
          [&format_str](auto & P) { 
            double sum=0;
            double error2_sum=0;
            for(auto & p: P)  
            {
              sum+=p.tt.effcor.value;  
              error2_sum+=p.tt.effcor.error*p.tt.effcor.error;
            }
            double average = sum/P.size(); 
            double error = sqrt(error2_sum)/P.size();
            char buf[1024];
            sprintf(buf,format_str.c_str(), (average-1.0)*100,error*100);
            return std::string(buf); 
          } 
        );
  };
  for( auto &  sr : SR) 
     prn(sr.title,sr);
  hline();
  auto f = fold(SR);
  prn("cor", f);
  prn2("(cor-1)x100", f);
  hline("━");
};

inline void print_luminosity(Scan_t & data) {
  //char buf[1024*16];
  printf("%5s %10.6s %10.6s %10.6s %10.6s %10s %10.6s %10.6s %10s %10s\n",
      "point", "E,GeV", "Lonline,pb^-1","Lgg, pb^-1", "dLgg,pb^-1", "Ngg","Lbb,pb^-1", "dLbb,pb^-1", "Nbb", "Lbb/Lgg"
      );
  ibn::valer<double> Lgg={0,0};
  ibn::valer<double> Lee={0,0};
  double Lonline=0;
  long Ngg=0;
  long Nee=0;
  for(auto & sp : data) {
    Lgg+=sp.gg.luminosity;
    Lee+=sp.bb.luminosity;
    Ngg+=sp.gg.N;
    Nee+=sp.bb.N;
    Lonline+=sp.luminosity.value;
    printf("%5s %10.6f %10.6f %10.6f %10.6f %10ld %10.6f %10.6f %10ld %10.6f \n",
        sp.title.c_str(), sp.energy.value, sp.luminosity.value, 
        sp.gg.luminosity.value*1e-3, sp.gg.luminosity.error*1e-3, sp.gg.N,
        sp.bb.luminosity.value*1e-3, sp.bb.luminosity.error*1e-3, sp.bb.N,
        sp.bb.luminosity.value/sp.gg.luminosity.value
        );
  }
  printf("%5s %10s %10.6f %10.6f %10.6f %10ld %10.6f %10.6f %10ld %10.6f %10.6f\n",
      "total", "", Lonline, Lgg.value*1e-3, Lgg.error*1e-3, Ngg,
      Lee.value*1e-3, Lee.error*1e-3, Nee,
      (Lee/Lgg).value, (Lee/Lgg).error
      );
}

//void print(const std::vector<PointSelectionResult_t> & Points, int opt=1 , int first_column_width=10, int last_column_width=5, int  column_width= 8, int vline_width=6)
//{
//  SelectionResult_t sr;
//  sr = Points;
//  print(sr,0);
//}

//void print(const  std::vector<SelectionResult_t> & SR )
//{
//  PrintConfig_t cfg;
//  print_Ntt(SR,cfg);
//  print_Ngg(SR[0],cfg);
//  print_Nbb(SR[0],cfg);
//};
inline std::string print_tex(const std::vector<SelectionResult_t> & SR,std::string ResultTitle="", std::string fit_file="")
{
  std::ostringstream os;
  os << R"(
  \documentclass[a4paper,12pt]{article}
  \usepackage[T2A]{fontenc}
  \usepackage[utf8]{inputenc} 
  \usepackage[english,russian]{babel}
  \usepackage{graphics}
  \usepackage{epsfig}
  \usepackage{amsmath,amssymb}
  \usepackage{ctable}
  \usepackage[
    a4paper, 
    mag=1000,
    nofoot, 
    nohead,
    left=0.5cm, 
    right=0.5cm, 
    top=0.5cm, 
    bottom=0.5cm
  ]{geometry}
  \renewcommand{\arraystretch}{1.2}
  \newcommand{\myC}[1]{\multicolumn{1}{c}{#1}}
  \newcommand{\myCF}[1]{\multicolumn{1}{c}{\hspace{4ex}#1\hspace{4ex}}}
  \newcommand{\myL}[1]{\multicolumn{1}{l}{#1}}
  \newcommand{\myR}[1]{\multicolumn{1}{r}{#1}}
  \begin{document}
  )" << "\n";
  os << R"(\centering \textbf{ \Large {)" << ResultTitle << R"(}}\\)" << "\n";
  os << R"(\today)" << "\n";



  os << R"(\begin{table}[h!])" << "\n";
  os << R"(\caption{Количество отобранных $\tau$ пар})" << "\n";
  os << R"(\centering)" << "\n";
  os << R"(\begin{tabular}[h!]{)";
  os << "c";
  for(int i=0; i< SR[0].size(); ++i) os << "r";
  os << "|"<<"r}" << R"(\\)";
  os << R"(\specialrule{.1em}{.05em}{.05em})" << "\n";
  int col_width=5;
  os << std::setw(20) << "channel";
  for( auto & p : SR[0]) {
    std::string name = R"(\myCF{)" + p.name + R"(})";
    os << " & " << std::setw(col_width) << name;
  }
  os << std::setw(col_width) << " & " << "total"  << R"(\\ \hline)";
  os << "\n";

  os << std::setw(20) << R"($\int L$, $pb^{-1}$)";
  for( auto & p : SR[0]) {
    os << " & ";
    char buf[1024];
    sprintf(buf,"%5.1f", p.luminosity.value);
    os << std::setw(col_width) << buf;
  }
  os << " & ";
  {
    char buf[1024];
    sprintf(buf,"%5.1f", (SR[0].L));
    os << std::setw(col_width) << buf;
  }
  os<< R"(\\ \hline)" <<  "\n";

  os << std::setw(20) << R"($E-M_{\tau}^{\text{PDG}}$, MeV)";
  for( auto & p : SR[0]) {
    os << " & ";
    char buf[1024];
    sprintf(buf,"%3.2f", (p.energy.value*0.5 - MTAU)*1000);
    os << std::setw(col_width) << R"(\myR{)" <<  buf << "}";
  }
  os << " & " << R"(\\ \hline)" <<  "\n";

  

  os << R"(\renewcommand{\arraystretch}{1.1})" << "\n";
  for (auto & sr : SR) {
    std::string title = sr.title; 
    title = sub(title,"μ",R"(\mu)");
    title = sub(title,"π",R"(\pi)");
    title = sub(title,"ρ",R"(\rho)");
    title = "$ " + title + " $";
    os << std::setw(20) << title << " & ";
    for(auto & p : sr) {
      os << std::setw(col_width) << p.tt.N << " & ";
    }
    os << std::setw(col_width) << sr.Ntt << R"(\\)" << "\n";
  }
  os << R"(\hline)" << "\n";
  auto total = fold(SR);
  os << std::setw(20) << "all" << " & "; 
  for( auto & p : total ) {
    os << std::setw(col_width) << p.tt.N << " & ";
  }
  os << std::setw(col_width) << total.Ntt << R"(\\)" << "\n";
  os << R"(\specialrule{.1em}{.05em}{.05em})" << "\n";
  os << R"(\end{tabular})" << "\n";
  os << R"(\end{table})" << "\n";

  os << R"(\begin{table}[h!])" << "\n";
  os << R"(\caption{Эффективность регистрации $\tau$ пар})" << "\n";
  os << R"(\centering)" << "\n";
  os << R"(\resizebox{\textwidth}{!}{)" << "\n";
  os << R"(\begin{tabular}[h!]{)";
  os << "c";
  for(int i=0; i< SR[0].size(); ++i) os << "r";
  //os << "|"<<"r";
  os << "}" << R"(\\)";
  os << R"(\specialrule{.1em}{.05em}{.05em})" << "\n";
  col_width=20;
  os << std::setw(20) << "channel";
  for( auto & p : SR[0]) {
    std::string name = R"(\myC{)" + p.name + "}";
    os << " & " << std::setw(col_width) << name;
  }
  //os << std::setw(col_width) << "average ";
  os << R"(\\ \hline)";
  os << "\n";
  os << R"(\renewcommand{\arraystretch}{1.1})" << "\n";
  for (auto & sr : SR) {
    std::string title = sr.title; 
    title = sub(title,"μ",R"(\mu)");
    title = sub(title,"π",R"(\pi)");
    title = sub(title,"ρ",R"(\rho)");
    title = "$ " + title + " $";
    os << std::setw(20) << title;
    for(auto & p : sr) {
      os << " & ";
      char buf[1024];
      sprintf(buf, "$%4.3f \\pm %4.3f$", p.tt.efficiency.value*100, p.tt.efficiency.error*100);
      os << std::setw(col_width) << buf;
    }
    //os << std::setw(col_width) << sr.Ntt;
    os << R"(\\)" << "\n";
  }
  os << R"(\hline)" << "\n";
  os << std::setw(20) << "all"; 
  for( auto & p : total ) {
      os << " & ";
      char buf[1024];
      sprintf(buf, "$%4.3f \\pm %4.3f$", p.tt.efficiency.value*100, p.tt.efficiency.error*100);
      os << std::setw(col_width) << buf;
  }
  //os << std::setw(col_width) << total.Ntt;
  os << R"(\\)" << "\n";
  os << R"(\specialrule{.1em}{.05em}{.05em})" << "\n";
  os << R"(\end{tabular})" << "\n";
  os << R"(})" << "\n";
  os << R"(\end{table})" << "\n";

  os << R"(\begin{table}[h!])" << "\n";
  os << R"(\caption{Поправка к эффективности регистрации $\tau$ пар})" << "\n";
  os << R"(\centering)" << "\n";
  os << R"(\resizebox{\textwidth}{!}{)" << "\n";
  os << R"(\begin{tabular}[h!]{)";
  os << "c";
  for(int i=0; i< SR[0].size(); ++i) os << "r";
  //os << "|"<<"r";
  os << "}" << R"(\\)";
  os << R"(\specialrule{.1em}{.05em}{.05em})" << "\n";
  col_width=20;
  os << std::setw(20) << "channel";
  for( auto & p : SR[0]){
    std::string name = R"(\myC{)" + p.name + "}";
    os << " & " << std::setw(col_width) << name;
  }
  //os << std::setw(col_width) << "average ";
  os << R"(\\ \hline)";
  os << "\n";
  os << R"(\renewcommand{\arraystretch}{1.1})" << "\n";
  for (auto & sr : SR) {
    std::string title = sr.title; 
    title = sub(title,"μ",R"(\mu)");
    title = sub(title,"π",R"(\pi)");
    title = sub(title,"ρ",R"(\rho)");
    title = "$ " + title + " $";
    os << std::setw(20) << title;
    for(auto & p : sr) {
      os << " & ";
      char buf[1024];
      sprintf(buf, "$%4.3f \\pm %4.3f$", p.tt.effcor.value, p.tt.effcor.error);
      os << std::setw(col_width) << buf;
    }
    //char buf[1024];
    //sprintf(buf, "$%4.3f \\pm %4.3f$", sr.effcor, sr.effcor_error);
    //os << std::setw(col_width) << sr.Ntt;
    os << R"(\\)" << "\n";
  }
  os << R"(\hline)" << "\n";
  os << std::setw(20) << "all"; 
  for( auto & p : total ) {
      os << " & ";
      char buf[1024];
      sprintf(buf, "$%4.3f \\pm %4.3f$", p.tt.effcor.value, p.tt.effcor.error);
      os << std::setw(col_width) << buf;
  }
  //os << std::setw(col_width) << total.Ntt;
  os << R"(\\)" << "\n";
  os << R"(\specialrule{.1em}{.05em}{.05em})" << "\n";
  os << R"(\end{tabular})" << "\n";
  os << R"(})" << "\n";
  os << R"(\end{table})" << "\n";
  if(fit_file!="")
  {
    auto total = fold(SR);
    os << R"(\begin{table}[h!])" << "\n";
    os << R"(\caption{Количество отобранных $\tau$ пар})" << "\n";
    os << R"(\centering)" << "\n";
    os << R"(\resizebox{\textwidth}{!}{)" << "\n";
    os << R"(\begin{tabular}[h!]{)";
    os << "c";
    for(int i=0; i< SR[0].size(); ++i) os << "r";
    os << "|"<<"r}" << R"(\\)";
    //os << R"(\specialrule{.1em}{.05em}{.05em})" << "\n";
    int col_width=5;

    auto make_row = [&col_width](
        std::ostream & os, 
        const SelectionResult_t & sr, 
        std::string title, 
        std::function<std::string(const PointSelectionResult_t & p)> Data,
        std::function<std::string(void)> LastColumn = [](){ return ""; })
    {
      os << std::setw(20) << title;
      for( auto & p : sr) {
        os << " & ";
        os << std::setw(col_width) << Data(p);
      }
      os << " & ";
      os << LastColumn();
      os << R"(\\)" << "\n";
    };


    make_row(os, total, 
        R"(point)", 
        [](const PointSelectionResult_t & p) { 
        //return ibn::format(R"(\myCF{%s})", p.name.c_str()); },
        return ibn::format(R"(%s)", p.name.c_str()); },
        [&total]() { return "Total"; }
        );
    os<< R"(\hline)" <<  "\n";


    //print luminosity
    make_row(os, total, 
        R"($\int L, pb^{-1}$)", 
        //[](const PointSelectionResult_t & p) { return ibn::format(R"(\myR{%4.1f})", p.luminosity.value); },
        [](const PointSelectionResult_t & p) { return ibn::format(R"(%4.1f)", p.luminosity.value); },
        [&total]() { return ibn::format(R"(%4.1f)", total.L); }
        );

    //print energy for point
    make_row(os, total, 
        R"($E-M_{\tau}$, MeV)", 
        //[](const PointSelectionResult_t & p) { return ibn::format(R"(\myR{%3.3f})", (p.energy.value*0.5 - MTAU)*1000.0); }
        [](const PointSelectionResult_t & p) { return ibn::format(R"(%3.3f)", (p.energy.value*0.5 - MTAU)*1000.0); }
        );
    os<< R"(\hline)" <<  "\n";

    os << R"(\renewcommand{\arraystretch}{1.1})" << "\n";
    for (auto & sr : SR) {
      std::string title = sr.title; 
      title = sub(title,"μ",R"(\mu)");
      title = sub(title,"π",R"(\pi)");
      title = sub(title,"ρ",R"(\rho)");
      title = "$ " + title + " $";
      make_row(os, sr, title,  
        [](const PointSelectionResult_t & p) { return ibn::format(R"(%d)", p.tt.N); },
        [&sr]() { return R"(\textit{)"+std::to_string(sr.Ntt)+"}"; });
    }
    os << R"(\hline)" << "\n";

    //print total event number
    make_row(os, total, 
        "all", 
        [](const PointSelectionResult_t & p) { return ibn::format(R"(\textbf{%d})", p.tt.N); },
        [&total]() { return R"(\textbf{)"+std::to_string(total.Ntt)+"}"; }
        );
    os<< R"(\hline)" <<  "\n";

    //print epsilon
    make_row(os, total, 
        R"($\varepsilon$, \%)", 
        [](const PointSelectionResult_t & p) { return ibn::format(R"($%5.3f\pm%5.3f$)", p.tt.efficiency.value*100, p.tt.efficiency.error*100); });

    //print epsilon correction
    make_row(os, total, 
        R"($\epsilon^{cor}$)", 
        [](const PointSelectionResult_t & p) { return ibn::format(R"($%6.4f\pm%6.4f$)", p.tt.effcor.value, p.tt.effcor.error); });

    make_row(os, total, 
        R"($(\epsilon^{cor}-1)\cdot100\%$)", 
        [](const PointSelectionResult_t & p) { return ibn::format(R"($%5.3f\pm%5.3f$)", (p.tt.effcor.value-1)*100.0, p.tt.effcor.error*100.0); });

    //os << R"(\specialrule{.1em}{.05em}{.05em})" << "\n";
    os << R"(\hline)" << "\n";
    os << R"(\end{tabular})" << "\n";
    os << "}\n";
    os << R"(\end{table})" << "\n";


    os << R"(\begin{figure}[h!])" << "\n";
    os << R"(\centering)" << "\n";
    os << R"(\includegraphics[width=0.8\textwidth]{)" << fit_file << "}\n";
    os << R"(\caption{)" <<sub(fit_file,"_",R"(\_)")  << "}\n";
    os << R"(\end{figure})" << "\n";
  }
  os << R"(\end{document})";
  os << "\n";
  return os.str();
};

void make_tex(std::string latex, std::string filename = "scan.tex")
{
  std::ofstream ofs(filename);
  ofs << latex;
  ofs.close();
  char command[65536];
  std::cout << "Compiling latex file..." << filename << std::endl;
  sprintf(command, "pdflatex %s > /dev/null", filename.c_str());
  system(command);
};

#endif
