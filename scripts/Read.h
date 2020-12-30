/*
 * =====================================================================================
 *
 *       Filename:  Read.h
 *
 *    Description: Read data from files 
 *
 *        Version:  1.0
 *        Created:  30.12.2020 11:40:29
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */
#ifndef IBN_TAU_READ_H
#define IBN_TAU_READ_H
#include <string>

#include "utils.h"

#include "ScanPoint.h"
#include "Print.h"

struct PointDiscriminatorByRunList
{
  const Scan_t * scan;
  PointDiscriminatorByRunList(const Scan_t * scan) : scan(scan) {}
  std::string operator()(std::string file)
  {
    std::regex re(R"(.*\D+(\d+)\.root)");
    std::smatch match;
    std::string result;
    int run=0;
    if(std::regex_match(file, match,re))
    {
      if(match.size()>1)
      {
        run = std::stoi(match[1]);
        for( auto & s : *scan)
        {
          if(std::count(std::begin(s.run_list),std::end(s.run_list), run) > 0)
          {
            result=s.title;
          }
        }
      }
    }
    return result;
  }
};

void set_alias(TTree * tt, double W, double L=1.0)
{
  tt->SetAlias("good_emc_time", "0<=ntemc[0]&&ntemc[0]<=14&&0<=ntemc[1]&&ntemc[1]<=14");
  tt->SetAlias("cgood","abs(cos(theta[0]))<0.93 && abs(cos(theta[1]))<0.93");
  tt->SetAlias("barrel","abs(cos(theta[0]))<0.8 && abs(cos(theta[1]))<0.8");
  tt->SetAlias("missed_photon_angle2","( abs(cos_theta_mis2) < 0.8 || ( 0.92 > abs(cos_theta_mis2) && abs(cos_theta_mis2) > 0.86) )");
  tt->SetAlias("missed_photon_angle", "( abs(cos_theta_mis) < 0.8 || ( 0.92 > abs(cos_theta_mis) && abs(cos_theta_mis) > 0.86) )");
  tt->SetAlias("missed_photon_angle3", "( abs(cos(d4p3[3])) < 0.8 || ( 0.92 > abs(cos(d4p3[3])) && abs(cos(d4p3[3])) > 0.86))");
  tt->SetAlias("missed_photon_angle4", "( 0.92 > abs(cos_theta_mis) )");
  tt->SetAlias("missed_photon_angle5", "( ( 0.93 > abs(cos_theta_mis) && abs(cos_theta_mis) > 0.84) || (0.83 > abs(cos_theta_mis)) )");
  tt->SetAlias("MM","(p[0]+p[1])**2 - (px[0]+px[1])**2 - (py[0]+py[1])**2 - (pz[0]+pz[1])**2");
  tt->SetAlias("Epi0","sqrt(p[0]**2 + 0.1396**2)");
  tt->SetAlias("Epi1","sqrt(p[1]**2 + 0.1396**2)");
  tt->SetAlias("Emu0","sqrt(p[0]**2 + 0.1056**2)");
  tt->SetAlias("Emu1","sqrt(p[1]**2 + 0.1056**2)");
  tt->SetAlias("EK0","sqrt(p[0]**2 + 0.4937**2)");
  tt->SetAlias("EK1","sqrt(p[1]**2 + 0.4937**2)");
  tt->SetAlias("psum","sqrt((px[0]+px[1])**2 + (py[0]+py[1])**2 + (pz[0]+pz[1])**2)");
  tt->SetAlias("M2mumu","(Emu0+Emu1)**2-psum*psum");
  tt->SetAlias("M2pipi","(Epi0+Epi1)**2-psum*psum");
  tt->SetAlias("M2emu","(p[0]+Emu1)**2-psum*psum");
  tt->SetAlias("M2mue","(p[1]+Emu0)**2-psum*psum");
  tt->SetAlias("M2KK","(EK0+EK1)**2-psum*psum");
  //tt->SetAlias("NnE100","Sum$(nE>0.1)");
  //tt->SetAlias("NnE50","Sum$(nE>0.050)");
  //tt->SetAlias("NnE25","Sum$(nE>0.025)");
  //tt->SetAlias("Ngbarrel","Sum$( abs(cos(ntheta))<0.8 && nE>0.025)");
  //tt->SetAlias("Ngendcup","Sum$( abs(cos(ntheta))<0.92 && abs(cos(ntheta)) >0.86  && nE>0.05)");
  //tt->SetAlias("Ng","Ngbarrel+Ngendcup");



  tt->SetAlias("lum",myfmt("%6.2f*1",L).c_str());
  char Eb[1024];
  sprintf(Eb,"%5.3f*1",W*0.5);
  tt->SetAlias("Eb",Eb);
  tt->SetAlias("MPI",(std::to_string(MPION)+"*1").c_str());
  tt->SetAlias("Emis2","(2*Eb-p[0]-p[1])");
  //tt->SetAlias("cos_theta_mis","(pz[0]+pz[1])/Emis");
  tt->SetAlias("cos_theta_mis2","(pz[0]+pz[1])/hypot(hypot(px[0]+px[1], py[0]+py[1]), pz[0]+pz[1])");
  tt->SetAlias("acol2","(px[0]*px[1]+py[0]*py[1]+pz[0]*pz[1])/(p[0]*p[1])");
  tt->SetAlias("MM2","Emis**2 - (px[0]+px[1])**2 - (py[0]+py[1])**2 - (pz[0]+pz[1])**2");
  tt->SetAlias("MM2pi","(2*Eb-sqrt(p[0]**2+MPI**2)-sqrt(p[1]**2+MPI**2))**2 - (px[0]+px[1])**2 - (py[0]+py[1])**2 - (pz[0]+pz[1])**2");
  tt->SetAlias("M2pi","(sqrt(p[0]**2+MPI**2)+sqrt(p[1]**2+MPI**2))**2 - (px[0]+px[1])**2 - (py[0]+py[1])**2 - (pz[0]+pz[1])**2");
  tt->SetAlias("k_ptem","(1.0*1.0)");

  tt->SetAlias("ptem2","sqrt((px[0]+px[1]+npx[0]+npx[1])^2 + (py[0]+py[1]+npy[0]+npy[1])^2)/Emis");

  auto set_alias = [](TTree * tt, std::string track_name_prefix, int track, std::string selection_template) -> void
  {
    std::string track_name = track_name_prefix + std::to_string(track);
    tt->SetAlias(track_name.c_str(), sub(selection_template,R"(#)",std::to_string(track)).c_str());
  };

  for(int i=0;i<4;i++)
  {
    set_alias(tt,"Mrho_cut", i, "0.5 < Mrho[#] && Mrho[#] < 1.0");
  }

  //define channels
  auto make_xy  = [](TTree * tt,  std::string x, std::string y)
  {
    auto s1 = sub("((x0&&y1)||(y0&&x1))", "x", x);
    auto alias = sub(s1, "y", y);
    tt->SetAlias( (x+y).c_str(), alias.c_str());
  };
  for (auto & x : {"e","u", "pi", "K","rho","X"} )
    for (auto & y : {"e","u", "pi", "K","rho","X"} )
      make_xy(tt, x,y);

  tt->SetAlias("pil", "(pi0 && (u1 || e1 ))");
  tt->SetAlias("lpi", "(pi1 && (u0 || e0 ))");

  //tt->SetAlias("eX",   "((!e0 && e1)||(!e1 && e0))");
  tt->SetAlias("Xrho",   "((!rho0 && rho1)||(!rho1 && rho0))");
  tt->SetAlias("lpipi0_cut","(Nc==2 && Npi0==1 && ( (pil && Mrho[0]>0.5 && Mrho[0]<1.0 && ) || ((lpi && (Mrho[1]>0.5 && Mrho[1]<1.0))) ) && acop>1.4)");
  tt->SetAlias("all_cut",   "lpipi0_cut || eu_cut || epi_cut");
}

Scan_t read_data(std::string data_dir, const Scan_t & cfg, std::string filter=R"(\.root$)")
{
  std::cout << "Reading data directory \"" << data_dir << "\":\n";
  auto fl  = filter_file_list(get_recursive_file_list(data_dir)); //get file list *.root
  //for(auto & f : fl) std::cout << f << std::endl;
  auto ml  = combine(fl,PointDiscriminatorByRunList(&cfg)); //combine files into points
  //for(auto & m : ml) 
  //{
  //  std::cout << m.first << ": ";
  //  for(auto & f: m.second) std::cout << f<<" "; 
  //  std::cout << std::endl;
  //}
  Scan_t scan;
  scan.reserve(cfg.size());
  for(auto & point : cfg ) 
  {
    ScanPoint_t sp;
    sp = point;
    sp.file_list = ml[point.title];
    if(sp.file_list.empty()) 
    {
      std::cout << "WARNING: no data for point: " << point.title << std::endl;
      continue;
    }
    scan.push_back(sp);
  }
  //read TChain
  int index = 0;
  for(auto & sp : scan)
  {
    //create chain
    //sp.tt= new TChain((, sp.title.c_str());
    //read tau tau events
    auto tt = new TChain("tt", sp.title.c_str());
    for(auto & file : sp.file_list) tt->AddFile(file.c_str());
    //read gamma gamma events
    auto gg = new TChain("gg", sp.title.c_str());
    for(auto & file : sp.file_list) gg->AddFile(file.c_str());

    auto bb = new TChain("bb", sp.title.c_str());
    for(auto & file : sp.file_list) bb->AddFile(file.c_str());
    //change chain names
    tt->SetName(("tt"+std::to_string(index)).c_str());
    bb->SetName(("bb"+std::to_string(index)).c_str());
    set_alias(tt,sp.energy,sp.luminosity);
    sp.tt.tree.reset(tt);
    sp.gg.tree.reset(gg);
    sp.bb.tree.reset(bb);
    sp.tt.energy = sp.energy;
    sp.gg.energy = sp.energy;
    sp.bb.energy = sp.energy;
    sp.type = data_dir;
  }
  print(scan);
  return scan;
}

Scan_t read_mh(std::string data_dir, const Scan_t & cfg, std::string filter=R"(\.root$)")
{
  std::cout << "Reading data directory \"" << data_dir << "\":\n";
  auto fl  = filter_file_list(get_recursive_file_list(data_dir)); //get file list *.root
  //for(auto & f : fl) std::cout << f << std::endl;
  auto ml  = combine(fl,PointDiscriminatorByRunList(&cfg)); //combine files into points
  //for(auto & m : ml) 
  //{
  //  std::cout << m.first << ": ";
  //  for(auto & f: m.second) std::cout << f<<" "; 
  //  std::cout << std::endl;
  //}
  Scan_t scan;
  scan.reserve(cfg.size());
  for(auto & point : cfg ) 
  {
    ScanPoint_t sp;
    sp = point;
    sp.file_list = ml[point.title];
    if(sp.file_list.empty()) 
    {
      std::cout << "WARNING: no data for point: " << point.title << std::endl;
      continue;
    }
    scan.push_back(sp);
  }
  //read TChain
  int index = 0;
  for(auto & sp : scan)
  {
    //create chain
    //sp.tt= new TChain((, sp.title.c_str());
    //read tau tau events
    auto tt = new TChain("data", sp.title.c_str());
    for(auto & file : sp.file_list) { 
      //std::cout << "Adding file " << file << std::endl;
      tt->AddFile(file.c_str());
    }
    //read gamma gamma events
    auto gg = new TChain("gg", sp.title.c_str());
    for(auto & file : sp.file_list) gg->AddFile(file.c_str());

    auto bb = new TChain("bb", sp.title.c_str());
    for(auto & file : sp.file_list) bb->AddFile(file.c_str());
    //change chain names
    tt->SetName(("tt"+std::to_string(index)).c_str());
    bb->SetName(("bb"+std::to_string(index)).c_str());
    set_alias(tt,sp.energy,sp.luminosity);
    sp.tt.tree.reset(tt);
    sp.gg.tree.reset(gg);
    sp.bb.tree.reset(bb);
    sp.tt.energy = sp.energy;
    sp.gg.energy = sp.energy;
    sp.bb.energy = sp.energy;
    sp.type = data_dir;
  }
  print(scan);
  return scan;
}

Scan_t read_privalov_lum(std::string data_dir, const Scan_t & cfg, std::string filter=R"(\.root$)")
{
  std::cout << "Reading data directory \"" << data_dir << "\":\n";
  auto fl  = filter_file_list(get_recursive_file_list(data_dir)); //get file list *.root
  //for(auto & f : fl) std::cout << f << std::endl;
  auto ml  = combine(fl,PointDiscriminatorByRunList(&cfg)); //combine files into points
  //for(auto & m : ml) 
  //{
  //  std::cout << m.first << ": ";
  //  for(auto & f: m.second) std::cout << f<<" "; 
  //  std::cout << std::endl;
  //}
  Scan_t scan;
  scan.reserve(cfg.size());
  for(auto & point : cfg ) 
  {
    ScanPoint_t sp;
    sp = point;
    sp.file_list = ml[point.title];
    if(sp.file_list.empty()) 
    {
      std::cout << "WARNING: no data for point: " << point.title << std::endl;
      continue;
    }
    scan.push_back(sp);
  }
  //read TChain
  int index = 0;
  for(auto & sp : scan)
  {
    //create chain
    //sp.tt= new TChain((, sp.title.c_str());
    //read tau tau events
    auto tt = new TChain("tt", sp.title.c_str());
    for(auto & file : sp.file_list) { 
      //std::cout << "Adding file " << file << std::endl;
      tt->AddFile(file.c_str());
    }
    //read gamma gamma events
    auto gg = new TChain("event", sp.title.c_str());
    for(auto & file : sp.file_list) gg->AddFile(file.c_str());

    auto bb = new TChain("bb", sp.title.c_str());
    for(auto & file : sp.file_list) bb->AddFile(file.c_str());
    //change chain names
    tt->SetName(("tt"+std::to_string(index)).c_str());
    tt->SetName(("gg"+std::to_string(index)).c_str());
    bb->SetName(("bb"+std::to_string(index)).c_str());
    set_alias(tt,sp.energy,sp.luminosity);
    sp.tt.tree.reset(tt);
    sp.gg.tree.reset(gg);
    sp.bb.tree.reset(bb);
    sp.tt.energy = sp.energy;
    sp.gg.energy = sp.energy;
    sp.bb.energy = sp.energy;
    sp.type = data_dir;
  }
  print(scan);
  return scan;
}


std::vector<ScanPoint_t> read_mc(std::string  dirname=".", Scan_t cfg={}, long N0mc=1e6, std::string regexpr=R"(.+\.root)")
{
  std::vector<ScanPoint_t> P;
  TSystemDirectory dir(dirname.c_str(), dirname.c_str());
  std::regex file_re(regexpr); //regular expression to filter files
  std::regex energy_re(R"(^\D*(\d+\.?\d+).root$)");
  std::smatch file_match;
  std::smatch energy_match; 
  int point=0;
  for(auto  file: * dir.GetListOfFiles())
  {
    std::string file_name(file->GetName());
    if(std::regex_match (file_name,file_match,file_re)) {
      //extracting the energy
      if(std::regex_match(file_name,energy_match,energy_re)) {
        double W = std::stod(energy_match[1]);
        if(cfg.empty()){
        }
        else {
          //find closest energy
          auto & sp  = *std::min_element(cfg.begin(), cfg.end(), [W](auto a, auto b){ return fabs(a.energy-W)<fabs(b.energy-W); } );
          if( abs(sp.energy - W) < 0.001 ) {
            P.push_back(sp);
            auto & p = P.back();
            P.back().type = dirname;
            p.title = dirname+p.title;
            p.title = sub(p.title,R"(mc/)","");
            p.title = sub(p.title,R"(T)","");
            p.energy = W;
            p.luminosity = -p.luminosity;
            p.tt.tree.reset(get_chain("tt",("tt"+p.title).c_str(), "signal", (dirname+"/"+file_name).c_str()));
            p.gg.tree.reset(get_chain("gg",("gg"+p.title).c_str(), "gg lum", (dirname+"/"+file_name).c_str()));
            p.bb.tree.reset(get_chain("bb",("bb"+p.title).c_str(), "bb lum", (dirname+"/"+file_name).c_str()));
            p.tt.energy = W;
            p.gg.energy = W;
            p.bb.energy = W;
            p.tt.N0mc = N0mc;
            p.gg.N0mc = N0mc;
            p.bb.N0mc = N0mc;
            set_alias(p.tt.tree.get(), p.energy,p.luminosity);
          }
        }
      }
      else
      {
        std::cout << "Unable to extract energy from filename: " << file_name << std::endl;
      }
    }
  }
  std::sort(P.begin(),P.end(),
      [](auto & sp1, auto &sp2) {
        return sp1.energy < sp2.energy;
        });

  return P;
}


std::vector<ScanPoint_t> read_mc_mh(std::string  dirname=".", Scan_t cfg={}, long N0mc=1e6, std::string regexpr=R"(.+\.root)")
{
  std::vector<ScanPoint_t> P;
  if(cfg.empty()) return P;
  TSystemDirectory dir(dirname.c_str(), dirname.c_str());
  std::regex file_re(regexpr); //regular expression to filter files
  std::regex take_re(R"(^(\D\d+).+.root$)");
  std::smatch file_match;
  std::smatch take_match; 
  int point=0;
  for(auto  file: * dir.GetListOfFiles()) { //looop over all files in derectory
    std::string file_name(file->GetName());
    if(std::regex_match (file_name,file_match,file_re)) { //take only *.root files
      if(std::regex_match(file_name,take_match,take_re)) { //take special file
        auto it = std::find_if(cfg.begin(), cfg.end(), [&take_match](const auto & a) { return take_match[1] == a.title; });
        if ( it != cfg.end() )  {
          auto & sp = *it;
          P.push_back(sp);
          auto & p = P.back();
          P.back().type = dirname;
          p.title = dirname+p.title;
          p.title = sub(p.title,R"(mc/)","");
          p.title = sub(p.title,R"(T)","");
          p.luminosity = -p.luminosity; //this is simulation
          p.tt.tree.reset(get_chain("data",("data"+p.title).c_str(), "signal", (dirname+"/"+file_name).c_str()));
          //p.gg.tree.reset(get_chain("gg",("gg"+p.title).c_str(), "gg lum", (dirname+"/"+file_name).c_str()));
          //p.bb.tree.reset(get_chain("bb",("bb"+p.title).c_str(), "bb lum", (dirname+"/"+file_name).c_str()));
          p.tt.energy = p.energy;
          p.gg.energy = p.energy;
          p.bb.energy = p.energy;
          p.tt.N0mc = N0mc;
          p.gg.N0mc = N0mc;
          p.bb.N0mc = N0mc;
        }
      }
    }
  }
  //sort by energy
  std::sort(P.begin(),P.end(),
      [](auto & sp1, auto &sp2) {
        return sp1.energy < sp2.energy;
        });

  return P;
}

std::vector<ScanPoint_t> read_mc_mhlum(std::string  dirname=".", Scan_t cfg={}, long N0mc=1e6, std::string regexpr=R"(.+\.root)")
{
  std::vector<ScanPoint_t> P;
  if(cfg.empty()) return P;
  TSystemDirectory dir(dirname.c_str(), dirname.c_str());
  std::regex file_re(regexpr); //regular expression to filter files
  std::regex take_re(R"(^(\D\d+).+.root$)");
  std::smatch file_match;
  std::smatch take_match; 
  int point=0;
  for(auto  file: * dir.GetListOfFiles()) { //looop over all files in derectory
    std::string file_name(file->GetName());
    if(std::regex_match (file_name,file_match,file_re)) { //take only *.root files
      if(std::regex_match(file_name,take_match,take_re)) { //take special file
        auto it = std::find_if(cfg.begin(), cfg.end(), [&take_match](const auto & a) { return take_match[1] == a.title; });
        if ( it != cfg.end() )  {
          auto & sp = *it;
          P.push_back(sp);
          auto & p = P.back();
          P.back().type = dirname;
          //p.title = dirname+p.title;
          //p.title = sub(p.title,R"(mc/)","");
          //p.title = sub(p.title,R"(T)","");
          p.luminosity = -p.luminosity; //this is simulation
          p.tt.tree.reset(get_chain("data",("data"+p.title).c_str(), "signal", (dirname+"/"+file_name).c_str()));
          p.gg.tree.reset(get_chain("event",("gg"+p.title).c_str(), "gg lum", (dirname+"/"+file_name).c_str()));
          p.bb.tree.reset(get_chain("bb",("bb"+p.title).c_str(), "bb lum", (dirname+"/"+file_name).c_str()));
          p.tt.energy = p.energy;
          p.gg.energy = p.energy;
          p.bb.energy = p.energy;
          p.tt.N0mc = N0mc;
          p.gg.N0mc = N0mc;
          p.bb.N0mc = N0mc;
        }
      }
    }
  }
  //sort by energy
  std::sort(P.begin(),P.end(),
      [](auto & sp1, auto &sp2) {
        return sp1.energy < sp2.energy;
        });

  return P;
}

#endif
