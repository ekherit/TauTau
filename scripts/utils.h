#pragma once
#ifndef IBN_TAUSELECTION_UTILS_H

#include <string>
#include <algorithm>
#include <iostream>
#include <iomanip>

/* ==================== WORKING WITH UTF-8 ================================================== */


template <class String > 
size_t count_utf8_symbols( String s )
{
  return std::count_if(begin(s),end(s),[](char c) { return (c & 0xc0) != 0x80; } );
}
template <class String > 
size_t count_utf8_extra_byte( String s )
{
  return s.size() - std::count_if(begin(s),end(s),[](char c) { return (c & 0xc0) != 0x80; } );
}
template< typename D>
void print_utf(int width, D d)
{
  std::cout << std::setw(width) << d;
};

template<>
void print_utf(int width, std::string s)
{

  std::cout << std::setw(width+count_utf8_extra_byte(std::string(s))) << s;
};

template<>
void print_utf(int width, const char * s)
{
  print_utf(width, std::string(s));
};

template<typename...Ts>
std::string myfmt(std::string format, Ts...ts)
{
  static char MY_FORMAT_BUF[65535];
  sprintf(MY_FORMAT_BUF, format.c_str(), ts...);
  return std::string(MY_FORMAT_BUF);
};

/* ========================================================================================== */

#include <regex>
#include <list>
#include <TSystemDirectory.h>
#include <TSystem.h>
#include <TChain.h>
/* =========================== WORKING WITH ROOT FILES IN DIRECTORY ================================*/
//Get recursive file list
std::list<std::string> get_recursive_file_list(std::string dirname)
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

std::list<std::string> filter_file_list(const std::list<std::string> LST, std::string regexpr=R"(.+\.root)")
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
std::map<std::string, std::list<std::string> > combine(const std::list<std::string> LST, Discriminator D)
{
  std::map<std::string, std::list<std::string> > fmap;
  for(const std::string & file : LST)
  {
    fmap[D(file)].push_back(file);
  }
  return fmap;
}

std::map<std::string, std::list<std::string> > combine(const std::list<std::string> LST, std::string regexpr=R"(^\D*(\d+\.?\d+).root$)")
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
void print(const std::map<std::string, std::list<std::string> > & fmap)
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

TChain * get_chain(std::string name, std::string newname, std::string title, std::string filename)
{
  TChain * chain = new TChain(name.c_str(), title.c_str());
  chain->AddFile(filename.c_str());
  chain->SetName(newname.c_str());
  return chain;
}

TChain * get_chain(const char * name, const char * newname, std::string title, const char * file_name)
{
  TChain * chain = new TChain(name, title.c_str());
  chain->AddFile(file_name);
  chain->SetName(newname);
  return chain;
}

TChain * get_chain(const char * name, const char * newname, std::string title, int run_begin, int run_end, std::string dir="")
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
std::string sub(std::string s, std::string regexpr, std::string substr)
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


std::tuple<const std::string_view, const std::string_view> head_tail2(const std::string_view s, const std::string_view delim) {
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

std::vector<std::string_view> split(std::string & str,const char * delim = "&&") {
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
std::string  remove_some_cuts(std::string var, std::string cut) {
  auto cuts = split(cut);
  //delete cuts which has 
  //std::string var_re = ;
  //var_re = var_re + var + ;
  std::regex re(R"(\b)"+var+R"(\b)");
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
      result += "&&" + cuts[i];
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
#endif
