#ifndef IBN_TAUSELECTION_UTILS_H
#define IBN_TAUSELECTION_UTILS_H


#include <string>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <type_traits>

/* ==================== WORKING WITH UTF-8 ================================================== */


template <class String > 
inline size_t count_utf8_symbols( String s )
{
  return std::count_if(begin(s),end(s),[](char c) { return (c & 0xc0) != 0x80; } );
}
template <class String > 
inline size_t count_utf8_extra_byte( String s )
{
  return s.size() - std::count_if(begin(s),end(s),[](char c) { return (c & 0xc0) != 0x80; } );
}
template< typename D>
inline void print_utf(int width, D d)
{
  std::cout << std::setw(width) << d;
};

template<>
inline void print_utf(int width, std::string s)
{

  std::cout << std::setw(width+count_utf8_extra_byte(std::string(s))) << s;
};

template<>
inline void print_utf(int width, const char * s)
{
  print_utf(width, std::string(s));
};

namespace ibn {
  template<typename...Ts>
   inline  std::string format(std::string f, Ts...ts)
    {
      static char MY_FORMAT_BUF[65535];
      snprintf(MY_FORMAT_BUF, sizeof(MY_FORMAT_BUF), f.c_str(), ts...);
      return std::string(MY_FORMAT_BUF);
    };

  template<typename...Ts>
    std::string mformat(std::string format, Ts...ts);

  template<>
    inline std::string mformat(std::string f) { return ""; };

  template<typename T, typename...Ts>
    inline std::string mformat(std::string f, T t, Ts...ts)
    {
      std::string lf="%";
      if constexpr (std::is_same_v<double,decltype(t)> ){
        lf+=f+"f";
      }
      if constexpr (std::is_same_v<int,decltype(t)> ){
        lf  += f.substr(0,f.find('.')) + "d";
      }
      if constexpr (std::is_same_v<long,decltype(t)> ){
        lf  += f.substr(0,f.find('.')) + "ld";
      }
      if constexpr (std::is_same_v<std::string,decltype(t)>  || std::is_same_v<const char *,decltype(t)> ){
        lf  += f.substr(0,f.find('.')) + "s";
      }
      //std::cout << f << "  " << lf << "  " << t << std::endl;
      return format(lf,t) + mformat(f, ts...);
    };
}

/* ========================================================================================== */

#include <regex>
#include <list>
#include <TSystemDirectory.h>
#include <TSystem.h>
#include <TChain.h>
/* =========================== WORKING WITH ROOT FILES IN DIRECTORY ================================*/
//Get recursive file list
inline std::list<std::string> get_recursive_file_list(std::string dirname)
{
  std::list<std::string> flst;
  TSystemDirectory dir(dirname.c_str(), dirname.c_str());
  for(auto  file: * dir.GetListOfFiles())
  {
    std::string file_name(file->GetName());
    if(file->IsFolder())
    {
      if(file_name != "." && file_name != "..")
      {
        flst.merge(get_recursive_file_list(dirname+"/"+file_name));
      }
    }
    else
    {
      flst.push_back(dirname+"/"+file_name);
    }
  }
  flst.sort();
  return flst;
}

inline std::list<std::string> filter_file_list(const std::list<std::string> LST, std::string regexpr=R"(.+\.root)")
{
  std::list<std::string> flst;
  std::regex file_re(regexpr); //regular expression to filter files
  std::smatch file_match;
  for( auto & file : LST)
  {
    if(std::regex_match (file,file_match,file_re))
    {
      flst.push_back(file);
    }
  }
  return flst;
}


template<class Discriminator>
inline std::map<std::string, std::list<std::string> > combine(const std::list<std::string> LST, Discriminator D)
{
  std::map<std::string, std::list<std::string> > fmap;
  for(const std::string & file : LST)
  {
    fmap[D(file)].push_back(file);
  }
  return fmap;
}

inline std::map<std::string, std::list<std::string> > combine(const std::list<std::string> LST, std::string regexpr=R"(^\D*(\d+\.?\d+).root$)")
{
  std::regex re(regexpr);
  return combine(LST, [&re](const std::string & file)
                          {  
                            std::string result;
                            std::smatch match;
                            if(std::regex_match(file, match,re))
                            {
                              if(match.size()>1)
                              {
                                result =  match[1];
                              }
                            }
                            return result;
                          }
                 );
};


//template < typename Map>
//void print(const Map & fmap)
inline void print(const std::map<std::string, std::list<std::string> > & fmap)
{
  for(const auto & point :  fmap)
  {
    for(const auto & file : point.second)
    {
      std::cout << point.first << ": " << file << std::endl;
    }
    std::cout << std::endl;
  }
}

inline TChain * get_chain(std::string name, std::string newname, std::string title, std::string filename)
{
  TChain * chain = new TChain(name.c_str(), title.c_str());
  chain->AddFile(filename.c_str());
  chain->SetName(newname.c_str());
  return chain;
}

inline TChain * get_chain(const char * name, const char * newname, std::string title, const char * file_name)
{
  TChain * chain = new TChain(name, title.c_str());
  chain->AddFile(file_name);
  chain->SetName(newname);
  return chain;
}

inline TChain * get_chain(const char * name, const char * newname, std::string title, int run_begin, int run_end, std::string dir="")
{
  TChain * chain = new TChain(name, title.c_str());
  for(int run = run_begin; run<=run_end; ++run)
  {
    char buf[1024];
    if(dir=="") sprintf(buf, "%d.root",run);
    else sprintf(buf, "%s/%d.root", dir.c_str(), run);
    chain->AddFile(buf);
  }
  chain->SetName(newname);
  return chain;
}

//substitute in string s by substring according to regexpr
//inline std::string sub(std::string s, std::string regexpr, std::string substr)
//{
//  std::string result;
//  std::regex re(regexpr);
//  std::regex_replace(std::back_inserter(result), s.begin(), s.end(), re, substr);
//  return result;
//};

inline std::string sub(std::string_view s, std::string regexpr, std::string substr)
{
  std::string result;
  std::regex re(regexpr);
  std::regex_replace(std::back_inserter(result), s.begin(), s.end(), re, substr);
  return result;
};

/*
std::vector<std::string> split(std::string s) {
  std::regex re(R"#(((.+)&&)+())#");



}
*/


//std::tuple<std::string_view, std::string_view> head_tail(std::string_view s) {
//  int open_brase = 0; //number of open brases
//  int damprs = 0; //number of && ampersands
//  auto it = s.begin();
//  while(it!=s.end() && !( open_brase == 0 && damprs>0) ) {
//    damprs=0;
//    switch(*it) {
//      case '(' : ++open_brase; break;
//      case ')' : --open_brase; break;
//      case '&' : 
//                 ++it;
//                 if(*it == '&') ++damprs; 
//                 break;
//      default: break;
//    };
//    ++it;
//  }
//  auto n = std::distance(s.begin(),it);
//  if(damprs!=0) n-=2;
//  return {std::string_view(s.begin(),n), std::string_view(it,std::distance(it,s.end()))};
//}
//
//
//std::vector<std::string_view> split_cut(std::string & str) {
//  std::vector<std::string_view> result;
//  std::string_view tail(str);
//  do {
//    auto  [s,t] = head_tail(tail);
//    result.push_back(s);
//    tail = t;
//  } while(!tail.empty());
//  return result;
//};


inline std::tuple<const std::string_view, const std::string_view> head_tail2(const std::string_view s, const std::string_view delim) {
  size_t open_brase{0}; 
  size_t delim_idx{0};
  const size_t delim_size{delim.size()};
  auto it = s.cbegin();
  for (;
      !( it == s.cend() || (delim_idx == delim_size && open_brase == 0));
        ++it
      )
  {
    switch(*it) {
      case '(' : 
        ++open_brase; 
        break;
      case ')' : 
        --open_brase; 
        break;
      default  : 
        break;
    };
    if( delim[delim_idx] == *it ) ++delim_idx;
    else delim_idx =0;
  }
  return {
    std::string_view(s.cbegin(), std::distance(s.cbegin(),      it)-delim_idx), 
    std::string_view(it,         std::distance(it,        s.cend()))
  };
}

inline std::vector<std::string_view> split(const std::string & str,const char * delim = "&&") {
  std::vector<std::string_view> result;
  std::string_view tail(str);
  do {
    auto  [s,t] = head_tail2(tail, delim);
    result.push_back(s);
    tail = t;
  } while(!tail.empty());
  return result;
};

/* Remove cuts connected with specified var */
inline std::string  remove_single_cuts(std::string var, std::string cut) {
  auto cuts = split(cut);
  //delete cuts which has 
  //std::string var_re = ;
  //var_re = var_re + var + ;
  var.erase(std::remove(std::begin(var), std::end(var), ' '), var.end());
  var.erase(std::remove(std::begin(var), std::end(var), '\t'), var.end());
  std::regex re(R"(\b)"+var+R"(\b)");
  //std::regex re("var");
  //std::smatch sm;
  std::match_results<std::string_view::const_iterator> sm;
  auto ptr = std::remove_if (cuts.begin(), cuts.end(), [&re, &sm](auto s) { 
        return std::regex_search(s.begin(), s.end(), sm, re);
        //return s.find(var) != std::string::npos; 
        }
        );
  //std::cout << "ptr distance = " << std::distance(cuts.begin(), ptr) << std::endl;
  cuts.erase( ptr , cuts.end());
  std::string result;
  //int idx=0;
  if(cuts.size()!=0) 
  {
    result = cuts[0];
    for(size_t i=1; i!=cuts.size(); ++i) {
      //result += "&&" + cuts[i];
      result += "&&";
      result += cuts[i];
      //result +=cuts[i];
      //if(cuts[i].find(var) !=  std::string::npos)  continue;
      //rcuts.push_back(
      //result += cuts[i];
      ////if(i < cuts.size()-1 ) result+="&&";
      //std::cout << i << " / " << cuts.size() << "  " << cuts[i]<< " " << result << std::endl;
    }
  }
  return result;
}
/* Remove cuts connected with specified var */

inline std::string  remove_some_cuts(std::string vars, std::string cut) {
  auto v = split(vars,"&&"); //split what cuts to remove
  std::string result = cut;
  for( auto & var : v ) {
    std::cout << var <<std::endl;
    result = remove_single_cuts(std::string(var),result);
  }
  return result;
}

#include <TGraph.h>
//find point with maximum y. Return corresponding <x,y>
inline std::pair<double, double> get_maximum(const TGraph & g) {
  //double max = std::numeric_limits<double>::lowest();
  auto ptr = std::max_element(g.GetY(), g.GetY()+g.GetN());
  size_t n = std::distance(g.GetY(), ptr);
  return { g.GetX()[n], *ptr };
}

#include "../ibn/valer.h"
#include <TGraphErrors.h>
inline void SetPoint(TGraphErrors * g, size_t i, const ibn::valer<double> & x, const ibn::valer<double> & y ) {
  g->SetPoint(i, x.value, y.value);
  g->SetPointError(i, x.error, y.error);
}

inline void AddPoint(TGraphErrors * g, const ibn::valer<double> & x, const ibn::valer<double> & y ) {
  size_t i = g->GetN();
  g->SetPoint(i, x.value, y.value);
  g->SetPointError(i, x.error, y.error);
}

inline std::string operator&&(const std::string & a, const std::string & b) {
  if (  a.empty() &&  b.empty() ) return "";
  if (  a.empty() && !b.empty() ) return b;
  if ( !a.empty() &&  b.empty() ) return a;
  return "("+a+")&&("+b+")";
}

inline void operator&=(std::string & a, const std::string & b) {
  if (  a.empty() &&  b.empty() ) return;
  if (  a.empty() && !b.empty() ) a=b;
  if ( !a.empty() &&  b.empty() ) return;
  a = a && b;
}

inline std::string operator||(const std::string & a, const std::string & b) {
  if (  a.empty() &&  b.empty() ) return "";
  if (  a.empty() && !b.empty() ) return b;
  if ( !a.empty() &&  b.empty() ) return a;
  //return "(("+a+")||("+b+"))";
  return a+"||"+b;
}

inline void operator|=(std::string & a, const std::string & b) {
  if (  a.empty() &&  b.empty() ) return;
  if (  a.empty() && !b.empty() ) a=b;
  if ( !a.empty() &&  b.empty() ) return;
  a = a || b;
}


inline std::string initial_common_part(const std::string & a, const std::string & b) {
  std::string s;
  size_t nmax = std::min(a.length(),b.length());
  for(size_t i=0; i!=nmax && a[i]==b[i] ;++i)  s+=a[i];
  return s;
}

std::string and_fold(const std::vector<std::string> & v) {
  if(v.size()==0) return "";
  std::string s = "("+v[0]+")";
  for(size_t i=1;i!=v.size(); ++i) {
    s&="("+v[i]+")";
  }
  return s;
};  
inline std::string common_cut(const std::string & a, const std::string & b) {
  auto prepare = [](const auto & v ) -> std::vector<std::string> {
    auto  va = split(v, "&&");
    std::vector<std::string> sa;
    //copy from string view to string and remove spaceses
    for( auto  sv :  va) {
      sa.push_back(sub(sv,R"(\s)",""));
    }
    //std::sort(std::begin(sa),std::end(sa));
    return sa;
  };
  auto sa = prepare(a);
  auto sb = prepare(b);
  //size_t nmax = std::min(sa.size(), sa.size());
  //std::string s;
  //for(size_t i=0; i!=nmax;++i)  {
  //  if( sa[i] == sb[i] ) s&=sa[i];
  //}
  std::vector<std::string> common_cuts(std::max(sa.size(), sb.size()));
  auto it = std::set_intersection(sa.begin(),sa.end(),sb.begin(),sb.end(), common_cuts.begin());
  common_cuts.resize(it-common_cuts.begin());
  //remove empty cuts
  common_cuts.erase(std::remove(common_cuts.begin(), common_cuts.end(), ""), common_cuts.end());
  return and_fold(common_cuts);
}


inline std::string differece_cut(const std::string & a, const std::string & b) {
  auto prepare = [](const auto & v ) -> std::vector<std::string> {
    auto  va = split(v, "&&");
    std::vector<std::string> sa;
    //copy from string view to string and remove spaceses
    for( auto  sv :  va) {
      sa.push_back(sub(sv,R"(\s)",""));
    }
    std::sort(std::begin(sa),std::end(sa));
    return sa;
  };
  auto sa = prepare(a);
  auto sb = prepare(b);
  std::vector<std::string> cuts(std::max(sa.size(), sb.size()));
  auto it = std::set_difference(sa.begin(),sa.end(),sb.begin(),sb.end(), cuts.begin());
  cuts.resize(it-cuts.begin());
  //remove empty cuts
  cuts.erase(std::remove(cuts.begin(), cuts.end(), ""), cuts.end());
  return and_fold(cuts);
}

template<typename Iterator, typename Func>
Iterator find_best(Iterator it_begin, Iterator it_end, Func F)
{
  using ItemType = typename Iterator::value_type;
  using Type = typename std::result_of<Func(ItemType)>::type;
  Type chi2  = std::numeric_limits<Type>::max();
  Iterator result=it_end;
  for(auto it = it_begin; it!=it_end; ++it)
  {
    auto x = F(*it);
    if( x < chi2 ) 
    {
      result = it;
      chi2 = x;
    }
  }
  return result;
};

template<typename Container, typename Func>
typename Container::iterator find_best(Container & c, Func F)
{
  return find_best(std::begin(c), std::end(c), F);
}

#include <sys/ioctl.h>
#include <termios.h>

//inline bool kbhit(void)
//{
//    termios term;
//    tcgetattr(0, &term);
//
//    termios term2 = term;
//    term2.c_lflag &= ~ICANON;
//    tcsetattr(0, TCSANOW, &term2);
//
//    int byteswaiting;
//    ioctl(0, FIONREAD, &byteswaiting);
//
//    tcsetattr(0, TCSANOW, &term);
//
//    return byteswaiting > 0;
//}

//inline bool kbhit() {
//  static const int STDIN = 0;
//  static bool initialized = false;
//
//  if (! initialized) {
//    // Use termios to turn off line buffering
//    termios term;
//    tcgetattr(STDIN, &term);
//    term.c_lflag &= ~ICANON;
//    tcsetattr(STDIN, TCSANOW, &term);
//    setbuf(stdin, NULL);
//    initialized = true;
//  }
//
//  int bytesWaiting;
//  ioctl(STDIN, FIONREAD, &bytesWaiting);
//  return bytesWaiting>0;
//}

inline bool kbhit(void)
{
  struct termios oldt, newt;
  int ch;
  int oldf;

  tcgetattr(STDIN_FILENO, &oldt);
  newt = oldt;
  newt.c_lflag &= ~(ICANON | ECHO);
  tcsetattr(STDIN_FILENO, TCSANOW, &newt);
  oldf = fcntl(STDIN_FILENO, F_GETFL, 0);
  fcntl(STDIN_FILENO, F_SETFL, oldf | O_NONBLOCK);

  ch = getchar();

  tcsetattr(STDIN_FILENO, TCSANOW, &oldt);
  fcntl(STDIN_FILENO, F_SETFL, oldf);

  if(ch != EOF)
  {
    //ungetc(ch, stdin);
    return true;
  }

  return false;
}

#endif
