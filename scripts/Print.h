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

#include "utils.h"
#include "ScanPoint.h"
#include "Selection.h"

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

void print(const std::vector<ScanPoint_t> & SPL)
{
  int point_number = 1;
  printer pr;
  //pr.add("#point",     "%-7d");
  pr.add("#name",       "%-10s");
  pr.add("W,GeV",      "%15.6f");
  pr.add("dW,GeV",     "%15.6f");
//  pr.add("E-Mtau,MeV", "%15.3f");
  pr.add("Sw,GeV",     "%15.6f");
  pr.add("dSw,GeV",    "%15.6f");
  pr.add("L,pb^-1",    "%15.6f");
  pr.add("dL,pb^-1",   "%15.6f");
  pr("  run list", "  %-10s");
  std::cout << pr.head() << std::endl;
  for (const auto & p : SPL) 
    std::cout <<  pr 
   //   % point_number++ 
      % p.title.c_str() 
      % p.energy.value % p.energy.error 
//      % ((0.5*p.energy.value-MTAU)/MeV) 
      % p.energy_spread.value
      % p.energy_spread.error
      % p.luminosity.value 
      % p.luminosity.error
      % get_run_formula(p.run_list).c_str() << std::endl;
}

ChannelSelectionResult_t  fold(const std::vector<ChannelSelectionResult_t>  & SR, std::string name = "all")
{
  ChannelSelectionResult_t result;
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
  auto & reference_point = * find_best(result, [](const auto & p) { return  std::abs(p.energy*0.5-MTAU+0.5); } );
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

ChannelSelectionResult_t  fold(std::vector<std::vector<ChannelSelectionResult_t>> & SR, std::string name = "all")
{
  std::vector<ChannelSelectionResult_t> VSR;
  for(const auto & vsr : SR) VSR.push_back(fold(vsr));
  return fold(VSR);
};

ChannelSelectionResult_t  fold(const std::vector<ChannelSelectionResult_t>  & s1, const std::vector<ChannelSelectionResult_t>  & s2)
{
  std::vector<ChannelSelectionResult_t> VSR;
  VSR.push_back(fold(s1));
  VSR.push_back(fold(s2));
  return fold(VSR);
};
void print(const ChannelSelectionResult_t & sr, int opt=1 , int first_column_width=10, int last_column_width=5, int  column_width= 8, int vline_width=6)
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


struct PrintConfig_t
{
  int title_width=10;
  int data_width=15;
  int total_width=5;
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

void print_eps(const std::vector<PointSelectionResult_t> pts, PrintConfig_t cfg=PrintConfig_t())
{
  print_smth(pts, cfg, [](){ return "ε"; }, [](auto & p) { return p.tt.efficiency.value; }, [](auto & s) { return ""; } );
};

void print_Ntt(const std::vector<PointSelectionResult_t> pts, PrintConfig_t cfg, std::string title="Ntt")
{
  print_smth(pts, cfg, 
      [&title](){ return title; }, 
      [](auto & p) { return p.tt.N; }, 
      [](auto & s) { 
        long N=0; 
        for(auto & x: s)  N+=x.tt.N;  
        return N; 
        } 
      );
};

void print_Ngg(const std::vector<PointSelectionResult_t> pts, PrintConfig_t cfg, std::string title="Ngg")
{
  print_smth(pts, cfg, 
      [&title](){ return title; }, 
      [](auto & p) { return p.gg.N; }, 
      [](auto & s) { 
        long N=0; 
        for(auto & x: s)  N+=x.gg.N;  
        return N; 
        } 
      );
};

void print_Nbb(const std::vector<PointSelectionResult_t> pts, PrintConfig_t cfg, std::string title="Nbb")
{
  print_smth(pts, cfg, 
      [&title](){ return title; }, 
      [](auto & p) { return p.bb.N; }, 
      [](auto & s) { 
        long N=0; 
        for(auto & x: s)  N+=x.bb.N;  
        return N; 
        } 
      );
};

template<typename Fdata, typename Ftotal>
void print_smth(const std::vector<ChannelSelectionResult_t> SR, PrintConfig_t cfg, Fdata fd, Ftotal ft)
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

void print_Ntt(const std::vector<ChannelSelectionResult_t> SR, PrintConfig_t cfg = PCFG)
{
  print_smth(SR,cfg, 
      [](auto & p) { 
        return p.tt.N; 
      },
      [](auto & s) { 
        double sum=0;
        for(auto & x: s)  sum+=x.tt.N;  
        return sum; 
      } 
  );
};

void print_Ngg(const std::vector<ChannelSelectionResult_t> SR, PrintConfig_t cfg = PCFG)
{
  print_smth(SR,cfg, 
      [](auto & p) { 
        return p.gg.N; 
      },
      [](auto & s) { 
        double sum=0;
        for(auto & x: s)  sum+=x.gg.N;  
        return sum; 
      } 
  );
};

void print_Nbb(const std::vector<ChannelSelectionResult_t> SR, PrintConfig_t cfg = PCFG)
{
  print_smth(SR,cfg, 
      [](auto & p) { 
        return p.bb.N; 
      },
      [](auto & s) { 
        double sum=0;
        for(auto & x: s)  sum+=x.bb.N;  
        return sum; 
      } 
  );
};


void print_efficiency(const std::vector<ChannelSelectionResult_t> SR, PrintConfig_t cfg = PCFG)
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

void print_effcor(const std::vector<ChannelSelectionResult_t> SR, PrintConfig_t cfg = PCFG)
{
  if(SR.empty()) return;
  std::string format_str = " %5.4f ± %5.4f ";
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
  prn("all", f);
  prn2("(cor-1)x100", f);
  hline("━");
};


//void print(const std::vector<PointSelectionResult_t> & Points, int opt=1 , int first_column_width=10, int last_column_width=5, int  column_width= 8, int vline_width=6)
//{
//  ChannelSelectionResult_t sr;
//  sr = Points;
//  print(sr,0);
//}

//void print(const  std::vector<ChannelSelectionResult_t> & SR )
//{
//  PrintConfig_t cfg;
//  print_Ntt(SR,cfg);
//  print_Ngg(SR[0],cfg);
//  print_Nbb(SR[0],cfg);
//};

#endif
