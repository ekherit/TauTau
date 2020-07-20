#pragma once
#include <string>
#include <vector>
#include "utils.h"

struct Restriction_t
{
  std::string name;
  double min;
  double max;
};

Restriction_t  restrict(std::string value, double min, double max) {
  return {value, min, max};
};

struct ParticleID_t {
  std::string name;
  std::vector<std::string> cuts;
  std::vector<Restriction_t> restrictions;
};

struct ChannelSelection_t
{
  std::string title; //channel name  (temrinal)
  std::string cut;
  std::vector<ParticleID_t> pid;
  //std::string root_title; //title for root to print on canvas
  std::string root_title() const { //remove greek letters
    std::string result = title;
    result = sub(result,"μ", "#mu");
    result = sub(result,"π", "#pi");
    result = sub(result,"ρ", "#rho");
    //std::cout << result << std::endl;
    return result;
  };

  std::string root_title2() const { //remove greek letters
    std::string result = title;
    result = sub(result,"μ", "mu");
    result = sub(result,"π", "pi");
    result = sub(result,"ρ", "rho");
    //std::cout << result << std::endl;
    return result;
  };
  std::string common_cut(void) const {return cut; }
};

struct Selection_t :  public ChannelSelection_t,  public std::vector<Selection_t>
{

  void apply_common_cut(void) {
    std::regex empty_re(R"(\s*)");
    std::smatch sm;
    for(auto & s: *this) {
      if(std::regex_match(s.cut, sm, empty_re)) s.cut = cut;
      else s.cut = cut + "&&" + s.cut;
      s.pid=pid;
    }
  }


  Selection_t & operator()(std::string channel_name) {
    std::string name = channel_name;
    name = sub(name,"u","μ");
    name = sub(name,"pi","π");
    name = sub(name,"rho","p");
    std::cout << name << std::endl;
    auto it =std::find_if(std::begin(*this),std::end(*this),[&name](auto & sel) { return name==sel.title; }); 
    if(it==end()) {
      std::cerr << "ERROR: Unable to find " << channel_name << " aka " << name << " returning first one" << std::endl;
      return front();
    }
    return *it;
  };

  const Selection_t & operator()(std::string channel_name) const {
    return (*this)(channel_name);
  };

  void print(std::ostream & os = std::cout, std::string opts="") const {
    os << "Title: " << title << '\n';
    os << "Common cuts: " << cut << '\n';
    std::vector<std::string_view> cc = split(cut);
    os << "Splited by && cuts:\n";
    char buf[1024];
    for(auto & c : cc) {
      std::string tmp(c);
      sprintf(buf, "%30s\n", tmp.c_str());
      os << buf;
    }
    int idx=0;
    for(const auto & s : *this) {
      sprintf(buf, "%5d %10s %20s\n", idx, s.title.c_str(), s.cut.c_str());
      std::cout << buf;
      ++idx;
    }
  }
};

typedef std::reference_wrapper<Selection_t> SelRef_t;
