/*
 * =====================================================================================
 *
 *       Filename:  RunTable.h
 *
 *    Description: Helper function to read runtables 
 *
 *        Version:  1.0
 *        Created:  30.12.2020 10:16:32
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */
#ifndef IBN_TAU_RUNTABLE_H
#define IBN_TAU_RUNTABLE_H 
#include <fstream>
#include "ScanPoint.h"

//create run list from single line, where runs in format [run1] [run2] [run3-run4] [run5]
inline std::list<int> get_run_list(std::string line);

//read privalov runtable 
Scan_t read_privalov_runtable(std::string filename);

enum {
  DATA_ALL,
  DATA_TAU,
  DATA_JPSI,
  DATA_PSIP,
  DATA_RES
};
//read my runtalbe
Scan_t read_my_runtable(std::string filename, int type=DATA_ALL);

//create run list from single line, where runs in format [run1] [run2] [run3-run4] [run5]
inline std::list<int> get_run_list(std::string line)
{
  std::list<int> run_list;
  if(line.empty()) return run_list;
  std::string str = line;
  bool is_match=true;
  static std::regex run_range_re(R"(^\s*(\d+)\s*-\s*(\d+))"); //look for the run range
  static std::regex run_re(R"(^\s*\d+)"); //look for the single run
  std::smatch match;
  while(is_match)
  {
    is_match = false;
    if(std::regex_search(str,match,run_range_re)) 
    {
      is_match = true;
      //match[0] -- <begin_run>-<end_run>
      //match[1] -- <begin_run>
      //match[2] -- <end_run>
      //std::cout << match[1] << '-'<<match[2] << std::endl;
      int begin_run = stod(match[1]);
      int end_run = stod(match[2]);
      for(int run = begin_run; run<=end_run; ++run) run_list.push_back(run);
      str = match.suffix();
    }
    else if( std::regex_search(str,match,run_re))
    {
      is_match = true;
      run_list.push_back(stod(match[0]));
      str = match.suffix();
    }
  }
  return run_list;
}

Scan_t read_privalov_runtable(std::string filename)
{
  Scan_t theScan;
  std::fstream ifs(filename);
  if(!ifs) { std::cerr << "Unable to open file: " << filename << std::endl; return theScan; }
  std::string point_name;
  std::string point_energy; //beam energy in GeV
  std::string point_runs;
  //std::regex re(R"((\d+-\d+)");
  //std::regex run_range_re(R"(^\s*(\d+)\s*-\s*(\d+))"); //look for the run range
  //std::regex run_re(R"(^\s*\d+)"); //look for the single run
  std::regex energy_re(R"(\d+\.?\d+)"); //for extract energy
  std::smatch match;
  while( ifs >> point_name >> point_energy) 
  { 
    std::getline(ifs,point_runs); 
    ScanPoint_t sp;
    sp.run_list = get_run_list(point_runs);
    if(sp.run_list.empty()) continue;
    sp.title = point_name;
    if(std::regex_search(point_energy, match, energy_re))
    {
      sp.energy = stod(match[0])*2.0;
    }
    else {
      sp.energy = 0;
    }
    theScan.emplace_back(std::move(sp));
  };
  return theScan;
}


inline Scan_t read_my_runtable(std::string filename, int type)
{
  Scan_t theScan;
  std::fstream ifs(filename);
  if(!ifs) { std::cerr << "Unable to open file: " << filename << std::endl; return theScan; }
  std::string point_name;
  std::string point_runs;
  ibn::valer<double>point_energy; //GeV
  ibn::valer<double>point_lum; 
  ibn::valer<double>point_spread;
  //double point_energy; //beam energy in GeV
  //double point_energy_error;
  //double point_spread;
  //double point_spread_error;
  std::regex comment_re(R"(^\s*#.*)");
  std::regex empty_re(R"(^\s*)");
  std::string line;
  std::smatch sm;
  while( std::getline(ifs,line)) 
  { 
    //if ( line.find_first_of('#') != std::string::npos ) continue;
    if( std::regex_match(line,sm,comment_re)  ) continue;
    if( std::regex_match(line,sm,empty_re)  ) continue;
    std::istringstream iss(line);
    iss >> point_name >> point_energy.value >> point_energy.error >> point_spread.value >> point_spread.error >> point_lum.value >> point_lum.error;
    std::getline(iss,point_runs); 
    ScanPoint_t sp;
    sp.run_list = get_run_list(point_runs);
    if(sp.run_list.empty()) 
    {
      std::istringstream iss2(point_runs);
      std::string run_combine_regexpr;
      while(iss2>>run_combine_regexpr) sp.regexprs.push_back(run_combine_regexpr);
    }
    else for(auto & r : sp.run_list) sp.regexprs.push_back(std::to_string(r));
    sp.title = point_name;
    sp.energy = point_energy;
    sp.energy_spread = point_spread; 
    sp.luminosity = point_lum;
    bool is_add=false;
    switch (type) {
      case DATA_JPSI:
        if(3000*MeV < sp.energy.value && sp.energy.value < 3200*MeV) is_add=true;
        break;
      case DATA_TAU:
        if(3500*MeV < sp.energy.value && sp.energy.value < 3630*MeV) is_add=true;
        break;
      case DATA_PSIP:
        if(3650*MeV < sp.energy.value && sp.energy.value < 3800*MeV) is_add=true;
        break;
      case DATA_RES:
        if(3000*MeV < sp.energy.value && sp.energy.value < 3200*MeV) is_add=true;
        if(3650*MeV < sp.energy.value && sp.energy.value < 3800*MeV) is_add=true;
        break;
      default:
        is_add=true;
    };
    if(is_add) theScan.emplace_back(std::move(sp));
    //sp.energy_spread.error = point_spread_error;
  //    std::cout << sp.energy.value << " " << M << std::endl;
    //if(fabs(sp.energy.value - M) < 60*MeV || type == DATA_ALL) {
    //  std::cout << sp.energy.value << " " << M << std::endl;
    //}
  };
  return theScan;
}

Scan_t read_simple_runtable(std::string filename)
{
  Scan_t theScan;
  std::fstream ifs(filename);
  if(!ifs) { std::cerr << "Unable to open file: " << filename << std::endl; return theScan; }
  std::string point_name;
  std::string point_runs;
  ibn::valer<double> point_energy; //GeV
  double point_spread;
  double point_spread_error;
  double point_lum;
  std::smatch sm;
  std::regex sharp_re(R"(^\s*#.*)"); //for extract energy
  std::regex files_re(R"(\w+\.root)");
  std::string line;
  while(std::getline(ifs,line)) 
  { 
    if(regex_match(line,sm,sharp_re)) continue;
    std::stringstream iss(line);
    //std::cout << point_name << std::endl;
    iss >> point_name >> point_energy.value >> point_energy.error >> point_spread >> point_spread_error >> point_lum;
    std::getline(iss,point_runs); 
    ScanPoint_t sp;
    for (std::sregex_iterator it(point_runs.begin(), point_runs.end(), files_re); it != std::sregex_iterator(); ++it)
    {
      sp.file_list.push_back(it->str());
    }
    //sp.begin_run=-1;
    //sp.end_run=-1;
    //std::cout << point_name << " " << point_energy << " " << std::endl;
    //sp.run_list = get_run_list(point_runs);
    //if(sp.run_list.empty()) continue;
    sp.title = point_name;
    sp.energy = point_energy;
    //sp.dW = point_energy_error;
    //sp.Sw.value = point_spread;
    //sp.Sw.error = point_spread_error;
    sp.energy_spread = point_spread;
    sp.luminosity = point_lum;
    theScan.emplace_back(std::move(sp));
  };
  return theScan;
}

#endif
