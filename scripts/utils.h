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


#endif
