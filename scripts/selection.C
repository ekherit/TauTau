/*
 * =====================================================================================
 *
 *       Filename:  selection.C
 *
 *    Description:  Tau tau pair selection
 *
 *        Version:  1.0
 *        Created:  15.06.2018 09:26:36
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */

#include <regex>
#include <string>
#include <list>
#include <fmt/printf.h>
#include <algorithm>

const double GeV=1.0;
const double MeV=1e-3*GeV;

const double MTAU=1776.86*MeV;
const double MPI=0.13957061*GeV; 
const double MPI0=0.134977*GeV;

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

struct ScanPoint_t
{
  std::string title;
  int begin_run;
  int end_run;
  double W;
  double dW;
  TTree * tt;
  TTree * gg;
  double L=0;
  std::list<std::pair<int,double> > runs;
  int Ntt;
  int Ngg;
  std::string selection;
  std::map<std::string, int> NttMap;
  double ppi;
  double ppi_max;
  double ppi_min;
  double eps=1.0; //registration efficiency
  double eps_error=1.0; //registration efficiency error
  double effcor=1.0; //registration efficiency correction
  double Sw, dSw; //energy spread
  std::list<int> run_list;
  std::list<std::string> file_list;
  std::string scan_title;
};
void print(const std::vector<ScanPoint_t> & SPL)
{
  int point_number = 1;
  std::cout << setw(5) << "n/n" << setw(10) << "begin run" << setw(10) <<  "end_run" << setw(15) << "IL, nb^-1" << setw(20) << "W/2-MTAU, MeV" << setw(20) << "dE, MeV" << std::endl;
  for (const auto & p : SPL)
  {
    std::cout << setw(5) << point_number << setw(10) << p.begin_run << setw(10) <<  p.end_run << setw(15) << p.L << setw(20) << (p.W/2-MTAU)/MeV << setw(20) << p.dW/2/MeV << std::endl;
    point_number++;
  }
}

typedef std::vector<ScanPoint_t>  Scan_t;
std::map<std::string, std::string > SelMap;

struct PointDiscriminatorByRunList
{
  const Scan_t * scan;
  PointDiscriminatorByRunList(const Scan_t * scan) : scan(scan) {}
  std::string operator()(std::string file)
  {
    std::regex re(R"(\D+(\d+)\.root)");
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


std::string get_run_formula(std::list<int> & runs)
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


void print_runs(const Scan_t  & scan)
{
  std::cout << setw(10) << "name" << setw(10) << "W,MeV" << setw(15) << "E-MTAU,MeV" << "       runs" << std::endl;
  for(auto s : scan)
  {
    std::cout << setw(10) << s.title << setw(10) << s.W/MeV << setw(15) << (s.W*0.5-MTAU)/MeV << "       " << get_run_formula(s.run_list) << std::endl;
  }
}

//create run list from single line, where runs in format [run1] [run2] [run3-run4] [run5]
std::list<int> get_run_list(std::string line)
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
      std::cout << match[1] << '-'<<match[2] << std::endl;
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
      sp.W = stod(match[0])*2.0;
    }
    else
    {
      sp.W = 0;
    }
    theScan.emplace_back(std::move(sp));
  };
  return theScan;
}

Scan_t read_my_runtable(std::string filename)
{
  Scan_t theScan;
  std::fstream ifs(filename);
  if(!ifs) { std::cerr << "Unable to open file: " << filename << std::endl; return theScan; }
  std::string point_name;
  std::string point_runs;
  double point_energy; //beam energy in GeV
  double point_energy_error;
  double point_spread;
  double point_spread_error;
  double point_lum;
  while( ifs >> point_name >> point_energy >> point_energy_error >> point_spread >> point_spread_error >> point_lum) 
  { 
    std::getline(ifs,point_runs); 
    ScanPoint_t sp;
    sp.run_list = get_run_list(point_runs);
    if(sp.run_list.empty()) continue;
    sp.title = point_name;
    sp.W = point_energy;
    sp.dW = point_energy_error;
    sp.Sw = point_spread;
    sp.dSw = point_spread_error;
    sp.L = point_lum;
    theScan.emplace_back(std::move(sp));
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
  double point_energy; //beam energy in GeV
  double point_energy_error;
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
    iss >> point_name >> point_energy >> point_energy_error >> point_spread >> point_spread_error >> point_lum;
    std::getline(iss,point_runs); 
    ScanPoint_t sp;
    for (sregex_iterator it(point_runs.begin(), point_runs.end(), files_re); it != sregex_iterator(); ++it)
    {
      sp.file_list.push_back(it->str());
    }
    sp.begin_run=-1;
    sp.end_run=-1;
    //std::cout << point_name << " " << point_energy << " " << std::endl;
    //sp.run_list = get_run_list(point_runs);
    //if(sp.run_list.empty()) continue;
    sp.title = point_name;
    sp.W = point_energy;
    sp.dW = point_energy_error;
    sp.Sw = point_spread;
    sp.dSw = point_spread_error;
    sp.L = point_lum;
    theScan.emplace_back(std::move(sp));
  };
  return theScan;
}

const char * make_alias(int channel, const char * templ )
{
  return templ;
}


void set_alias(TTree * tt, double W)
{
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

  char Eb[1024];
  sprintf(Eb,"%5.3f*1",W*0.5);
  tt->SetAlias("Eb",Eb);
  tt->SetAlias("Emis","(2*Eb-p[0]-p[1])");
  tt->SetAlias("cos_theta_mis","(pz[0]+pz[1])/Emis");
  tt->SetAlias("cos_theta_mis2","(pz[0]+pz[1])/hypot(hypot(px[0]+px[1], py[0]+py[1]), pz[0]+pz[1])");
  tt->SetAlias("acol2","(px[0]*px[1]+py[0]*py[1]+pz[0]*pz[1])/(p[0]*p[1])");
  tt->SetAlias("MM2","Emis**2 - (px[0]+px[1])**2 - (py[0]+py[1])**2 - (pz[0]+pz[1])**2");

  //particles id

  //define electrons
  for(int i=0;i<2;i++)
  {
    //electron id
    char alias[65535];
    sprintf(alias,"0.8 < Ep[%1$d] && Ep[%1$d]<1.05 && chi2_dedx_e[%1$d]<5 && abs(delta_tof_e[%1$d])<0.3",i);
    char particle[16];
    sprintf(particle,"e%d",i);
    std::cout << particle << "=" << alias << std::endl;
    tt->SetAlias(particle,alias);
    //muon id
  }
  //tt->SetAlias("e0","0.8 < Ep[0] && Ep[0]<1.05 && chi2_dedx_e[0]<5 && abs(delta_tof_e[0])<0.3");
  //tt->SetAlias("e1","0.8 < Ep[1] && Ep[1]<1.05 && chi2_dedx_e[1]<5 && abs(delta_tof_e[0])<0.3");
  tt->SetAlias("ee", "e0 && e1");

  //define muons
  //tt->SetAlias("u0","0.1 < E[0] && E[0] < 0.3 && (depth[0]>(80*p[0]-50) || depth[0]>40) && Nmuhit[0]>1 && chi2_dedx_mu[0] < 5  && abs(delta_tof_mu[0]) < 0.3 && Ep[0]<0.8");
  tt->SetAlias("u0","0.1 < E[0] && E[0] < 0.3 && depth[0]>0 && chi2_dedx_mu[0] < 5  && abs(delta_tof_mu[0]) < 0.3 && Ep[0]<0.8");
  ///tt->SetAlias("u1","0.1 < E[1] && E[1] < 0.3 && (depth[1]>(80*p[1]-50) || depth[1]>40) && Nmuhit[1]>1 && chi2_dedx_mu[1] < 5  && abs(delta_tof_mu[1]) < 0.3 && Ep[1]<0.8");
  tt->SetAlias("u1","0.1 < E[1] && E[1] < 0.3 && depth[1]>0 && chi2_dedx_mu[1] < 5  && abs(delta_tof_mu[1]) < 0.3 && Ep[1]<0.8");
  tt->SetAlias("uu", "u0 && u1");
  tt->SetAlias("eu", "(e0 && u1) || (u0 && e1)");

  //define pions
  tt->SetAlias("pi0","!u0 && 0.8 < p[0] && p[0] < 1.05 && Ep[0]<0.6 && chi2_dedx_pi[0] < 5 && abs(delta_tof_pi[0])<0.3");
  tt->SetAlias("pi1","!u1 && 0.8 < p[1] && p[1] < 1.05 && Ep[1]<0.6 && chi2_dedx_pi[1] < 5 && abs(delta_tof_pi[1])<0.3");
  //tt->SetAlias("pi0","!u0 &&   Ep[0]<0.6 && chi2_dedx_pi[0] < 5 && abs(delta_tof_pi[0])<0.3");
  //tt->SetAlias("pi1","!u1 &&   Ep[1]<0.6 && chi2_dedx_pi[1] < 5 && abs(delta_tof_pi[1])<0.3");


  tt->SetAlias("pipi","pi0 && pi1");
  tt->SetAlias("epi", "(e0 && pi1) || (pi0 && e1)");
  tt->SetAlias("upi", "(u0 && pi1) || (pi0 && u1)");

  tt->SetAlias("pil", "(pi0 && (u1 || e1 ))");
  tt->SetAlias("lpi", "(pi1 && (u0 || e0 ))");

  //define kaons
  tt->SetAlias("K0","!u0 && !pi0 && 0.8 < p[0] && p[0] < 1.05 && Ep[0]<0.6 && chi2_dedx_K[0] < 3 && abs(delta_tof_K[0])<0.3");
  tt->SetAlias("K1","!u1 && !pi1 && 0.8 < p[1] && p[1] < 1.05 && Ep[1]<0.6 && chi2_dedx_K[1] < 3 && abs(delta_tof_K[1])<0.3");

  tt->SetAlias("eK", "(e0 && K1) || (K0 && e1)");
  tt->SetAlias("uK", "(u0 && K1) || (K0 && u1)");
  tt->SetAlias("piK","(pi0 && K1) || (K0 && pi1)");
  tt->SetAlias("KK","K0 && K1");

  //final selection
  tt->SetAlias("eu_cut",    "Nc==2 && Nn==0 && eu   && ptem>0.25");
  tt->SetAlias("epi_cut",   "Nc==2 && Nn==0 && epi  && ptem>0.3");
  tt->SetAlias("upi_cut",   "Nc==2 && Nn==0 && upi  && ptem>0.62");

  tt->SetAlias("ee_cut",    "Nc==2 && Nn==0 && ee   && ptem>0.5 && acop<2.7 && p[0]<0.9 && p[1]<0.9");
  tt->SetAlias("uu_cut",    "Nc==2 && Nn==0 && uu   && ptem>0.5 && acop<2.7 && p[0]<0.9 && p[1]<0.9");
  tt->SetAlias("pipi_cut",  "Nc==2 && Nn==0 && pipi && ptem>0.55");

  tt->SetAlias("eK_cut","Nc==2 && Nn==0 && eK");
  tt->SetAlias("piK_cut","Nc==2 && Nn==0 && piK && ptem >0.4");
  tt->SetAlias("uK_cut","Nc==2 && Nn==0 && uK && ptem >0.15");
  tt->SetAlias("KK_cut","Nc==2 && Nn==0 && KK && ptem >0.15");

  for(int i=0;i<4;i++)
  {
    const double MRHO_MIN = 0.5;
    const double MRHO_MAX = 1.0;
    char name[128];
    char cut[65535];
    sprintf(name, "Mrho%d_cut",i);
    sprintf(cut,"%.2f < Mrho[%d] && Mrho[%d] < %.2f", MRHO_MIN, i, i, MRHO_MAX);
    tt->SetAlias(name,cut);
  }

  tt->SetAlias("rho0", "pi0 && Mrho0_cut && Mpi0[0]>0.1128 && Mpi0[0]<0.1464");
  tt->SetAlias("rho1", "pi1 && Mrho1_cut && Mpi0[0]>0.1128 && Mpi0[0]<0.1464");

  //tt->SetAlias("lpipi0_cut","(Nc==2 && Npi0==1 && ( (pil && Mrho[0]>0.5 && Mrho[0]<1.0 && ) || ((lpi && (Mrho[1]>0.5 && Mrho[1]<1.0))) ) && acop>1.4)");
  tt->SetAlias("Xrho", "!rho0 && (e0 || u0 || pi0 || K0) && rho1 && !(e1 || u1 || K0)");
  tt->SetAlias("rhoX", "!rho1 && (e1 || u1 || pi1 || K1) && rho0 && !(e0 || u0 || K0)");
  tt->SetAlias("rhorho", "!(e0 || u0 || K0) &&  !(e1 || u1 || K1) && Mpi0[0]>0.1128 && Mpi0[0]<0.1464 && Mpi0[1]>0.1128 && Mpi0[1]<0.1464  && pi0 && pi1  && ((Mrho0_cut && Mrho2_cut) || (Mrho1_cut && Mrho3_cut))");

  tt->SetAlias("Xrho_cut","Nc==2 && Npi0==1 && acop>1.4 && (Xrho || rhoX)");
  tt->SetAlias("rhorho_cut","Nc==2 && Npi0 ==2 && rhorho && ptem>0.5");

  //tt->SetAlias("pipi0pipi0_cut","Nc==2 && Npi0==2 && pipi && ptem>0.6");

  tt->SetAlias("all_cut",   "lpipi0_cut || eu_cut || epi_cut");
}

std::vector<ScanPoint_t> read_data(std::string data_dir, const Scan_t & cfg)
{
  //auto cfg = read_privalov_runtable(privalov_runtable); //read scan configuration
  auto fl  = filter_file_list(get_recursive_file_list(data_dir)); //get file list *.root
  auto ml  = combine(fl,PointDiscriminatorByRunList(&cfg)); //combine files into points
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
    auto gg = new TChain("tt", sp.title.c_str());
    for(auto & file : sp.file_list) gg->AddFile(file.c_str());
    //change chain names
    tt->SetName(("tt"+to_string(index)).c_str());
    gg->SetName(("gg"+to_string(index)).c_str());
    set_alias(tt,sp.W);
    sp.tt = tt;
    sp.gg = gg;
  }
  return scan;
}

std::vector<ScanPoint_t> read_data2(std::string data_dir, std::string privalov_runtable)
{
  auto cfg = read_privalov_runtable(privalov_runtable); //read scan configuration
  auto fl  = filter_file_list(get_recursive_file_list(data_dir)); //get file list *.root
  auto ml  = combine(fl,PointDiscriminatorByRunList(&cfg)); //combine files into points
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
    auto gg = new TChain("tt", sp.title.c_str());
    for(auto & file : sp.file_list) gg->AddFile(file.c_str());
    //change chain names
    tt->SetName(("tt"+to_string(index)).c_str());
    gg->SetName(("gg"+to_string(index)).c_str());
    set_alias(tt,sp.W);
    sp.tt = tt;
    sp.gg = gg;
  }
  return scan;
}

std::vector<ScanPoint_t> read_data3(std::string data_dir, std::string cfg_file)
{
  Scan_t scan = read_simple_runtable(cfg_file); //read scan configuration
  print(scan);
  int index = 0;
  for(auto & sp : scan)
  {
    //create chain
    //sp.tt= new TChain((, sp.title.c_str());
    //read tau tau events
    auto tt = new TChain("tt", sp.title.c_str());
    for(auto & file : sp.file_list) tt->AddFile((data_dir+"/"+file).c_str());
    //read gamma gamma events
    auto gg = new TChain("tt", sp.title.c_str());
    for(auto & file : sp.file_list) gg->AddFile((data_dir+"/"+file).c_str());
    //change chain names
    tt->SetName(("tt"+to_string(index)).c_str());
    gg->SetName(("gg"+to_string(index)).c_str());
    set_alias(tt,sp.W);
    sp.tt = tt;
    sp.gg = gg;
  }
  return scan;
}


std::vector<ScanPoint_t> read_data(const char * dir)
{
  std::vector<ScanPoint_t> Points;
  Points.push_back({"Point1", 55116, 55155,  3539.482,  0.110,0,0});

  Points.push_back({"Point1p", 55157,55161,3550.872,0.182,0,0});
  //remove points above this is tune EMS

  Points.push_back({"Point2", 55162,55199,3552.849,0.093,0,0});

  Points.push_back({"Point3", 55200,55231,3553.934,0.08,0,0});
  Points.push_back({"Point4", 55232,55239,3560.356,0.157,0,0});
  Points.push_back({"Point5", 55240,55257, 3599.572,0.117,0,0});
  //Points energy are in MeV now I need to convert it to GeV
  for(int i=0;i<Points.size(); i++)
  {
    Points[i].W*=MeV;
    Points[i].dW*=MeV;
  }
  for(int i=0;i<Points.size(); i++)
  {
    char name_tt[1024];
    char name_gg[1024];
    sprintf(name_tt,"tt%d", i+1);
    sprintf(name_gg,"gg%d", i+1);
    Points[i].tt = get_chain("tt", name_tt, Points[i].title, Points[i].begin_run, Points[i].end_run, dir);
    Points[i].gg = get_chain("gg", name_gg, Points[i].title, Points[i].begin_run, Points[i].end_run, dir);
    set_alias(Points[i].tt, Points[i].W);
  }
  ifstream runinfo("tauscan2018_runinfo.txt");
  if(!runinfo) 
  { 
    std::cout << "Unable to open runinfo file" << std::endl;
  }
  int run;
  double lum;
  while(runinfo >> run >> lum)
  {
    char run_root_file[1024];
    sprintf(run_root_file, "%s/%d.root",dir,run);
    std::cout << run_root_file << std::endl;
    if(!gSystem->AccessPathName(run_root_file))
    {
      //std::cout << "File " << run_root_file << " exists" << " lum = " << lum << " nb^-1" << std::endl;
      for(int i=0;i<Points.size();++i)
      {
        if( run >= Points[i].begin_run && run<= Points[i].end_run )
        {
          Points[i].L += lum;
          Points[i].runs.push_back({run,lum});
          cout << "run = " << run <<  " br=" << Points[i].begin_run << " er=" << Points[i].end_run << " point = " << i+1 << "  IL = " << Points[i].L << endl;
        }
      }
    }
    runinfo.ignore(65535,'\n');
  }
  ////add monte carlo
  //int pmc = Points.size();
  //Points.push_back({"MC", 55300,55300, 1777*2,0.1,0,0});
  //Points[pmc].tt = new TChain("tt","monte carlo tau");
  //((TChain*)Points[pmc].tt)->AddFile("test-sim.root");
  //set_alias(Points[pmc].tt,Points[pmc].W);
  //Points[pmc].gg = new TNtupleD("ggmc","ggmc","i");
  //for(int i=0;i<1000;i++) Points[pmc].gg->Fill();
  return Points;
}

#include <regex>
//read monte carlo Galuga
std::vector<ScanPoint_t> read_galuga(std::string  dirname=".", std::string regexpr=R"(.+\.root)")
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
    if(std::regex_match (file_name,file_match,file_re))
    {
      //extracting the energy
      if(std::regex_match(file_name,energy_match,energy_re))
      {
        double W = std::stod(energy_match[1]);
        std::cout << file_name <<  "  W = " << W << "  " << W-MTAU << std::endl;
        if(W*0.5-MTAU > -0.010 &&  W*0.5-MTAU < 0.030) 
        {
          P.push_back({"P"+to_string(point), 55116, 55155,  W, 1e-5,0,0});
          P.back().tt = get_chain("tt",("tt"+to_string(point)).c_str(), "MC GALUGA", (dirname+"/"+file_name).c_str());
          P.back().gg = get_chain("gg",("gg"+to_string(point)).c_str(), "gg lum", (dirname+"/"+file_name).c_str());
          set_alias(P.back().tt, P.back().W);
        }
      }
      else
      {
        std::cout << "Unable to extract energy from filename: " << file_name << std::endl;
      }
    }
  }
  std::cout << "THE END" << std::endl;
  return P;
}


void draw_lum_per_run2(const std::vector<ScanPoint_t> & SPL)
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

void  draw_lum_per_run(const char * dir, int begin_run, int end_run)
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

void draw_tau_per_run(const std::vector<ScanPoint_t> & SPL, const char * selection = "")
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
      double L = Ngg/10*pow(p.W/MTAU*0.5,2);
      double dL = dNgg/10*pow(p.W/MTAU*0.5,2);
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

TGraphErrors * draw_result(const char * selection, const std::vector<ScanPoint_t> & Points)
{
  TGraphErrors * g = new TGraphErrors;
  long totalNtt=0;
  for(int i=0; i<Points.size();++i)
  {
    double Ntt = Points[i].tt->GetEntries(selection);
    double Ngg = Points[i].gg->GetEntries(selection);
    double xs  = 0;
    double dxs = 0;
    totalNtt += Ntt;
    if(Ngg != 0 )
    {
      xs = Ntt/Ngg;
      dxs = sqrt( Ntt/(Ngg*Ngg) + pow(Ntt/(Ngg*Ngg), 2.0)*Ngg );
    }
    g->SetPoint(i, Points[i].W/2.0-MTAU, xs);
    g->SetPointError(i, Points[i].dW/2.0, dxs);
    std::cout << i << " " << Points[i].W/2.0-MTAU << "  " << Ngg << "  " << Ntt << "   " <<  xs << std::endl;
  }
  std::cout << "Total number of tau-tau candidates:" << totalNtt << std::endl;
  g->SetMarkerStyle(21);
  g->Draw("ap");

  TCanvas * cacop = new TCanvas("acop","acop");
  cacop->Divide(2,3);
  for( int i=0;i<Points.size();++i) 
  { 
    cacop->cd(i+1);
    Points[i].tt->Draw("acop",selection);
  }
  TCanvas * cptem = new TCanvas("ptem","ptem");
  cptem->Divide(2,3);
  for( int i=0;i<Points.size();++i) 
  { 
    cptem->cd(i+1);
    Points[i].tt->Draw("ptem",selection);
  }
  TCanvas * cp = new TCanvas("p","p");
  cp->Divide(2,3);
  for( int i=0;i<Points.size();++i) 
  { 
    cp->cd(i+1);
    Points[i].tt->Draw("p",selection);
  }
  return g;
}

void select(std::vector<ScanPoint_t> & P, const char * varexp, const char * selection="", const char * gopt="col", std::string opt="")
{
  bool goff  = !(bool(strcmp("goff",gopt)));
  static int canvas_counter=0;
  char canvas_name[65535];
  sprintf(canvas_name, "tau_sel_%d",canvas_counter++);
  TCanvas * c = nullptr;
  if(!goff)
  {
    c = new TCanvas(canvas_name,"selection",0,0,1920*2/5.0*P.size(),1080*0.5);
    //auto c = new TCanvas();
    c->SetTitle(selection);
    //int ny = floor(sqrt(P.size()+1));
    //int nx = ny+1;
    c->Divide(P.size(),1);
  }
  int i = 1;
  bool add_to_prev_selection=false;
  if(opt=="add") add_to_prev_selection=true;
  for (auto & p : P)
  {
    if(!goff) c->cd(i++);
    int Ntt; 
    if(!goff) Ntt=p.tt->Draw(varexp,selection,gopt);
    else Ntt = p.tt->GetEntries(selection);
    p.Ntt=Ntt;
    if(opt=="clear") p.NttMap.clear();
    if(opt=="add") p.NttMap[selection]=Ntt; 
    char title[1024];
    sprintf(title, "%d events - %s", p.Ntt, p.title.c_str());
    if(!goff) p.tt->GetHistogram()->SetTitle(title);
    p.selection = selection;
    if(opt=="save")
    {
      p.Ngg = p.gg->GetEntries();
      std::cout << i-2 << " " << p.Ngg << endl;
    }
  }
  if(opt=="save")
  {
    TGraphErrors * g = new TGraphErrors;
    long totalNtt=0;
    long totalNgg = 0;
    std::ofstream ofs("scan.txt");
    double Ltot=0;
    for(int i=0; i<P.size();++i)
    {
      totalNtt += P[i].Ntt;
      totalNgg+=P[i].Ngg;
      Ltot += P[i].L;
    }
    double sigma_gg = totalNgg/Ltot;
    for(int i=0; i<P.size();++i)
    {
      double xs  = 0;
      double dxs = 0;
      int Ntt = P[i].Ntt;
      int Ngg = P[i].Ngg;
      double L = Ngg/(sigma_gg*pow(P[i].W/(2*MTAU),2.0));
      if(Ngg != 0 )
      {
        xs = Ntt/L;
        dxs = xs*sqrt( 1./Ntt + 1./Ngg);
      }
      g->SetPoint(i, P[i].W/2.0-MTAU, xs);
      g->SetPointError(i, P[i].dW/2.0, dxs);
      ofs << setw(5) << i <<  setw(15) << P[i].L << "  " << 10 << setw(15) << P[i].W  << setw(15) << P[i].dW;
      ofs << setw(10) << 1.256 << " " << setw(10) << 0.019;
      ofs << setw(10) << Ntt << setw(10) << " " << 1 << "  " << Ngg <<  " " << 1 << std::endl;
    }
    if(!goff)
    {
      new TCanvas;
      g->SetMarkerStyle(21);
      g->Draw("ap");
    }
  }
}


void add_last(std::vector<ScanPoint_t> & P)
{
  for ( auto & p : P)
  {
    p.NttMap[p.selection]=p.Ntt;
  }
}
void set_last(std::vector<ScanPoint_t> & P)
{
  for ( auto & p : P)
  {
    p.NttMap.clear();
    p.NttMap[p.selection]=p.Ntt;
  }
}

void fit(std::vector<ScanPoint_t> & P, const char * filename="scan.txt", bool nofit=false, std::string title="")
{
  long totalNtt=0;
  long totalNgg = 0;
  double totalL=0;
  std::string total_title="";
  for(int i=0; i<P.size();++i)
  {
    for(auto & item: P[i].NttMap)
    {
      totalNtt += item.second;
    }
    if(P[i].Ngg==0) P[i].Ngg = P[i].gg->GetEntries();
    totalNgg+=P[i].Ngg;
    totalL += P[i].L;
  }
  //total_title = total_title + " " + P[i].selection;
  double sigma_gg;
  if( totalNgg == 0 || totalL == 0) sigma_gg = 1;
  else sigma_gg = totalNgg == 0 ? 1 : totalNgg/totalL;
  std::ofstream ofs(filename);
  ofs << "#Selection: " << P[0].selection << std::endl;
  ofs << "#Aliases: " << std::endl;
  ofs << "#    e0 = " << P[0].tt->GetAlias("e0") << std::endl;
  ofs << "#    e1 = " << P[0].tt->GetAlias("e1") << std::endl;
  ofs << "#    u0 = " << P[0].tt->GetAlias("u0") << std::endl;
  ofs << "#    u1 = " << P[0].tt->GetAlias("u1") << std::endl;
  ofs << "#    pi0 = " << P[0].tt->GetAlias("pi0") << std::endl;
  ofs << "#    pi1 = " << P[0].tt->GetAlias("pi1") << std::endl;
  ofs << "#" << setw(4) << " " << setw(15) << "L, nb-1" << setw(10) << "dL, nb-1" << setw(15) << "W, MeV" << setw(15) << "dW, MeV" << setw(10) << "SW, MeV" << setw(10) << "dSW, MeV";
  ofs << setw(10) << "Ntt" << setw(10) << "Nee" << setw(10) << "Ngg" << setw(10) << "effcor" << std::endl;
  for(int i=0; i<P.size();++i)
  {
    int Ntt=0;
    for(auto & item: P[i].NttMap) Ntt+=item.second;
    int Ngg = P[i].Ngg;
    double L = Ngg==0 ? 1 : Ngg/(sigma_gg*pow(P[i].W/(2*MTAU),2.0));
    ofs << setw(5) << i <<  setw(15) << L << setw(10) << 10 << setw(15) << P[i].W/MeV  << setw(15) << P[i].dW/MeV;
    ofs << setw(10) << 1.256 << " " << setw(10) << 0.019;
    ofs << setw(10) << Ntt << setw(10) <<  1 << setw(10) << Ngg <<  setw(10) << P[i].effcor << std::endl;
  }
  if(!nofit)
  {
    char command[65536];
    //if(title=="") title=total_title;
    sprintf(command, "taufit --tau-spread=1.256 --title='sigma: %s' '%s' --output '%s.txt' &", title.c_str(), filename,filename);
    //sprintf(command, "taufit --tau-spread=1.256  '%s' --output '%s.txt' &",  filename,filename);
    system(command);
  }
}

void fit_last(std::vector<ScanPoint_t> & P, const char * filename="scan.txt", bool nofit=false, std::string title="")
{
  set_last(P);
  fit(P,filename,nofit,title);
}



void mccmp(const std::vector<ScanPoint_t> & P, const char * selection, const char * cut, const char * his="")
{
  std::string title = std::string("mc_data: ") + selection + " cut= " + cut;
  TCanvas * c = new TCanvas("data_mc_diff", title.c_str());
  TChain * chain = new TChain("data_mc_diff","data_mc_diff");
  set_alias((TTree*)chain, MTAU*2);
  for(int i=1;i<P.size();++i) chain->Add((TChain*)P[i].tt);
  std::string sel(selection);
  sel+=">>h1(";
  sel+=his;
  sel+=")";
  chain->Draw(sel.c_str(),cut,"goff");
  auto hdata  = chain->GetHistogram();
  P[P.size()-1].tt->SetLineColor(kRed);
  sel=selection;
  sel+=">>h2(";
  sel+=his;
  sel+=")";
  P[P.size()-1].tt->Draw(sel.c_str(), cut, "goff");
  auto hmc = P[P.size()-1].tt->GetHistogram();
  hmc->Draw();
  hmc->GetXaxis()->SetTitle(selection);
  hmc->SetLineWidth(2);
  hdata->Scale(double(hmc->GetEntries())/hdata->Integral());
  hdata->SetLineWidth(2);
  hdata->Draw("SAME");
}

TGraphErrors * draw_result2(const std::vector<ScanPoint_t> & Points, const char * selection)
{
  TGraphErrors * g = new TGraphErrors;
  long totalNtt=0;
  std::ofstream ofs("scan.txt");
  for(int i=0; i<Points.size();++i)
  {
    double Ntt = Points[i].tt->GetEntries(selection);
    double Ngg = Points[i].gg->GetEntries();
    double xs  = 0;
    double dxs = 0;
    totalNtt += Ntt;
    if(Ngg != 0 )
    {
      xs = Ntt/Ngg;
      dxs = sqrt( Ntt/(Ngg*Ngg) + pow(Ntt/(Ngg*Ngg), 2.0)*Ngg );
    }
    g->SetPoint(i, Points[i].W/2.0-MTAU, xs);
    g->SetPointError(i, Points[i].dW/2.0, dxs);
    std::cout << i << " " << Points[i].W/2.0-MTAU << "  " << Ngg << "  " << Ntt << "   " <<  xs << std::endl;
    auto & P = Points[i];
    ofs << setw(5) << i <<  setw(15) << 20000 << "  " << 10 << setw(15) << P.W/MeV  << setw(15) << P.dW/MeV;
    ofs << setw(10) << 1.256 << " " << setw(10) << 0.019;
    ofs << setw(10) << Ntt << setw(10) << " " << 1 << "  " << Ngg <<  " " << 1 << std::endl;
  }
  std::cout << "Total number of tau-tau candidates:" << totalNtt << std::endl;
  TCanvas * c = new TCanvas("cross","cross");
  c->cd();
  g->SetMarkerStyle(21);
  g->Draw("ap");

  std::vector<const char *> var = { "acop", "ptem", "p" ,"acol"};

  TCanvas * can = new TCanvas("param","param");
  can->Divide(var.size(),Points.size());

  for( int i=0;i<Points.size();++i) 
  { 
    for(int v=0;v<var.size();v++)
    {
      can->cd(1+i*var.size()+v);
      Points[i].tt->Draw(var[v],selection);
    }
  }
  return g;
}

auto DATA = read_data("data");
auto DATA11 = read_data3("tau2011","tau2011/runtable.txt");
auto MC = read_galuga("mc/signal");
auto EEee = read_galuga("mc/EEee");
auto EEuu = read_galuga("mc/EEuu");
auto EEkk = read_galuga("mc/EEkk");
auto EEpipi = read_galuga("mc/EEpipi");
auto BB = read_galuga("mc/BB");
auto UU = read_galuga("mc/uu");
auto HADR = read_galuga("mc/hadrons");
auto S = read_data("data", read_privalov_runtable("../runtable.txt"));


struct Selection_t
{
  std::string title; //channel name  (temrinal)
  std::string root_title; //title for root to print on canvas
  std::string cut;
};

//std::vector<std::string> CHANNELS =
//{
//  "eu_cut"     ,
//  "epi_cut"    ,
//  "upi_cut"    ,
//  "ee_cut"     ,
//  "uu_cut"     ,
//  "pipi_cut"   ,
//  "eK_cut"     ,
//  "uK_cut"     ,
//  "piK_cut"    ,
//  "KK_cut"     ,
//  "Xrho_cut"   ,
//  "rhorho_cut"
//};


std::vector<Selection_t> SELECTION
{
  {"eμ", "e#mu"     , "eu_cut"}     ,
  {"eπ", "e#pi"     , "epi_cut"}    ,
  {"μπ", "#mu#pi"   , "upi_cut"}    ,
  {"ee", "ee"       , "ee_cut"}     ,
  {"μμ", "#mu#mu"   , "uu_cut"}     ,
  {"ππ", "#pi#pi"   , "pipi_cut"}   ,
  {"eK", "eK"       , "eK_cut"}     ,
  {"πK", "#piK"     , "piK_cut"}    ,
  {"KK", "KK"       , "KK_cut"}     ,
  {"Xρ", "X#rho"    , "Xrho_cut"}   ,
  {"ρρ", "#rho#rho" , "rhorho_cut"}
};

std::vector<Selection_t> & DEFAULT_SELECTION=SELECTION;

void measure_efficiency(std::vector<ScanPoint_t> & PNTS, std::vector<Selection_t> & SEL=DEFAULT_SELECTION, long N0=100000)
{
  const char * varexp = "ptem:acop";
  std::vector<long> Ntt(PNTS.size(),0);
  std::vector<TGraphErrors *> graphs(SEL.size());
  for(auto & g : graphs) g = new TGraphErrors;
  TGraphErrors * total_g = new TGraphErrors;
  const int whead = 15;
  const int wcol = 20;
  std::cout << setw(whead) << "W,MeV";
  for(int p = 0; p < PNTS.size(); ++p)
  {
    std::cout << setw(wcol) << PNTS[p].W*1e3;
  }
  std::cout << std::endl;
  std::cout << setw(whead) << "E-MPDG,MeV";
  for(int p = 0; p < PNTS.size(); ++p)
  {
    std::cout << setw(wcol) << (PNTS[p].W*0.5-MTAU)*1e3;
  }
  std::cout << std::endl;
  for(int i = 0; i < SEL.size(); i++)
  {
    select(PNTS, varexp, SEL[i].cut.c_str(), "col" , "add");
    std::cout << setw(whead+count_utf8_extra_byte(SEL[i].title)) << SEL[i].title;
    //std::cout << setw(whead) << CHANNELS[i];
    for(int p = 0; p < PNTS.size(); ++p)
    {
      Ntt[p]+=PNTS[p].Ntt;
      PNTS[p].eps = double(PNTS[p].Ntt)/double(N0);
      graphs[i]->SetPoint(p,(PNTS[p].W*0.5-MTAU)*1e3, PNTS[p].eps);
      graphs[i]->SetPointError(p,PNTS[p].dW*1e3*0.5, sqrt(PNTS[p].eps*(1.0 - PNTS[p].eps)/N0));
      std::cout << setw(wcol) << PNTS[p].eps;
    }
    std::cout << std::endl;
  }
  std::cout << setw(whead) << "TOTAL";
  for(int p = 0;p<PNTS.size();++p) 
  {
    double eps = double(Ntt[p])/double(N0);
//    std::cout << setw(10) << PNTS[p].W << setw(10) << PNTS[p].dW;
    total_g->SetPoint(p, (PNTS[p].W*0.5-MTAU)*1e3, eps);
    double deps = sqrt(eps*(1.0-eps)/N0);
    total_g->SetPointError(p, PNTS[p].dW*1e3*0.5, deps);
    char buf[1024];
    sprintf(buf, "%7.5f±%7.5f", eps,deps);
    std::cout << setw(wcol)<<buf;
  }
  std::cout << std::endl;
  auto * mg = new TMultiGraph;
  mg->Add(total_g, "lp");
  total_g->SetLineColor(1);
  total_g->SetMarkerColor(1);
  total_g->SetMarkerStyle(20);
  total_g->SetMarkerSize(2);
  total_g->SetLineWidth(2);
  for(int c = 0;c<graphs.size();++c) 
  {
    graphs[c]->SetLineColor(c+2);
    graphs[c]->SetMarkerColor(c+2);
    graphs[c]->SetMarkerStyle(21+c);
    graphs[c]->SetMarkerSize(2);
    graphs[c]->SetLineWidth(2);
    mg->Add(graphs[c],"lp");
  }
  auto eff_c = new TCanvas;
  eff_c->Divide(1,2);
  eff_c->cd(1);
  mg->Draw("apl");
  mg->GetXaxis()->SetTitle("E-M_{#tau}, MeV");
  mg->GetYaxis()->SetTitle("efficiency ( #sigma_{vis} / #sigma )");
  //mg->SetMarkerStyle(21);
  //mg->SetMarkerSize(2);
  //mg->SetLineWidth(2);
  TGraphErrors * eps_cor_g = new TGraphErrors;
  std::cout << setw(15) << "COR";
  for(int p = 0; p<PNTS.size();++p)
  {
    eps_cor_g->SetPoint(p, total_g->GetX()[p], total_g->GetY()[p]/total_g->GetY()[0]);
    eps_cor_g->SetPointError(p, total_g->GetErrorX(p), total_g->GetErrorY(p)/total_g->GetY()[0]);
    PNTS[p].eps = total_g->GetY()[p];
    PNTS[p].eps_error = total_g->GetErrorY(p);
    PNTS[p].effcor = eps_cor_g->GetY()[p];
    std::cout << setw(10) << eps_cor_g->GetY()[p];
  }
  std::cout << std::endl;
  eff_c->cd(2);
  eps_cor_g->Draw("ap");
  eps_cor_g->Draw("apl");
  eps_cor_g->SetMarkerStyle(21);
  eps_cor_g->SetMarkerSize(2);
  eps_cor_g->SetLineWidth(2);
  eps_cor_g->GetXaxis()->SetTitle("E-M_{#tau}, MeV");
  eps_cor_g->GetYaxis()->SetTitle("efficiency correction");
}

void select_all(std::vector<ScanPoint_t> & P=DATA, std::vector<Selection_t> & SEL=SELECTION)
{
  //some printing configuration
  int first_column_width = 15;
  int last_column_width = 10;
  int column_width = 10;
  int vline_width = 10;
  int hline_width = first_column_width + last_column_width + vline_width + P.size()*column_width+3;
  auto hline = [&hline_width](std::string  symb="─", int width = 0) 
  { 
    int w = width == 0 ? hline_width : width;
    for(int i=0;i<w;++i)
    {
      std::cout << symb;
    }
    std::cout << std::endl;
    //auto c = cout.fill(); std::cout << setfill(symb) << setw() << symb << setfill(c) << std::endl; 
  };
  auto vline = [&vline_width](char symb='-', int width = 0) { std::cout << setw(width == 0? (vline_width-1) : (width-1)) << symb; };

  //default expression to show
  const char * varexp = "ptem:acop";
  std::vector<long> Ntt;
  std::vector<long> totalEventInChannel(SEL.size());
  Ntt.resize(P.size(),0);
  hline("━");
  std::cout << setw(first_column_width) << "CHANNEL/POINT";
  for(int i=0;i<Ntt.size();++i) std::cout << setw(column_width)  << P[i].title;
  std::cout << setw(vline_width) << " │ " << setw(last_column_width) << "TOTAL";
  std::cout << std::endl;
  hline();
  std::cout << setw(first_column_width) << "E-MTAU,MeV";
  for(int i=0;i<Ntt.size();++i) std::cout << setw(column_width) << (P[i].W*0.5-MTAU)*1e3;
  std::cout << setw(vline_width) << " │ ";
  std::cout << std::endl;
  hline();
  //for(int c = 0; c < CHANNELS.size(); c++)
  int channel=0;
  for(auto & sel : SEL)
  {
    select(P, varexp, sel.cut.c_str(), "goff" , "add");
    //std::cout << setw(first_column_width) << CHANNELS[c];
    std::cout << setw(first_column_width+count_utf8_extra_byte(sel.title)) << sel.title;
    totalEventInChannel[channel]=0;
    for(int i=0; i < Ntt.size(); i++)
    {
      totalEventInChannel[channel]+=P[i].Ntt;
      Ntt[i]+=P[i].Ntt;
      std::cout << setw(column_width) <<  P[i].Ntt;
    }
    std::cout << setw(vline_width) << " │ " << setw(last_column_width) << totalEventInChannel[channel];
    std::cout << std::endl;
    channel++;
  }
  hline();
  std::cout << setw(first_column_width) << "total";
  for(int i=0;i<Ntt.size();i++) std::cout  << setw(column_width) << Ntt[i];
  long totalEvent=0;
  for(auto n :totalEventInChannel) totalEvent+=n;
  std::cout << setw(vline_width) << " │ " << setw(last_column_width) << totalEvent;
  std::cout << std::endl;
  hline("━");
  fit(P, "scan.txt", 1);
}


void draw_hist_vs_energies(Scan_t & scan, const char * var, int Nbin, double xmin, double xmax, const char *selection)
{
  std::vector<int> color={kBlack, kRed, kBlue,kGreen+2};
  std::vector<int> line = {1,2,3};
  TLegend * l = new TLegend(0.8,0.8,1.0,1.0);
  for(int i=0;i<scan.size();++i)
  {
    std::string s = var;
    int c = i % color.size();
    int l = i / color.size();
    std::cout << i << " " << l << " " << c << std::endl;
    scan[i].tt->SetLineColor(color[c]);
    scan[i].tt->SetLineStyle(line[l]);
    scan[i].tt->SetLineWidth(3);
    s+=" >> h"+to_string(i)+"("+to_string(Nbin)+","+to_string(xmin)+","+to_string(xmax)+")";
    if(i==0) scan[i].tt->Draw(s.c_str(), selection);
    else scan[i].tt->Draw(s.c_str(),selection,"same");
  }
}






void compare(const ScanPoint_t & p1, const ScanPoint_t &p2, const char * varexp, const char * selection, const char * gopt)
{
  auto  c = new TCanvas;
  c->SetWindowSize(1200, 1800);
  c->Divide(1,2);
  c->cd(1);
  p1.tt->Draw(varexp,selection,gopt);
  c->cd(2);
  p2.tt->Draw(varexp,selection,gopt);
}


void correct_efficiency(Scan_t & D /*  data  */, Scan_t & MC /*  MC signal */)
{
  for(auto & P : D)
  {
    //if( (P.W*0.5-MTAU)/MeV < 1 ) P.effcor = 1.0;
    //double effcor = 1.0;
    double energy_distance = 1e10;
    for ( auto  & mc : MC)
    {
      if ( fabs(P.W-mc.W) < energy_distance )
      {
        std::cout << "mc.effcor = " << mc.effcor << " " << P.W << " " << mc.W << std::endl;
        P.effcor = mc.effcor;
        energy_distance = fabs(P.W-mc.W);
      }
    }
  }

  std::cout << "#Correct efficiency: " << std::endl;
  std::cout << "#" << setw(4) << " " << setw(10) << "W,MeV" << setw(10) << "E, MeV" << setw(15) << "E-Mtau, MeV" << setw(10) << "Ntt" << setw(10) << "effcor" << std::endl;
  for( int i = 0;i<D.size(); i++)
  {
    std::cout << setw(5) << i+1 << setw(10) << D[i].W/MeV << setw(10) << D[i].W*0.5/MeV << setw(15) << (D[i].W*0.5-MTAU)/MeV <<  setw(10) << D[i].Ntt << setw(10) << D[i].effcor << std::endl;
  }
}

void set_color_and_line(TAttLine *obj, int index)
{
  std::vector<int> color={kBlack, kRed, kBlue,kGreen+2};
  std::vector<int> line = {1,2,3};
  int c = index % color.size();
  int l = index / color.size();
  obj->SetLineColor(color[c]);
  obj->SetLineStyle(line[l]);
}

void compare(std::vector<ScanPoint_t*> SP, std::vector<std::string> varexp, const char * selection, std::string gopt)
{
  auto  c = new TCanvas;
  int win_size=600;
  c->SetWindowSize(win_size*(SP.size()+1), win_size*varexp.size());
  int index =0;
  c->Divide(SP.size()+1,varexp.size());
  for(auto & ve : varexp)
  {
    TLegend * l = new TLegend(0.8,0.9,1.0,1.0);
    auto hs = new THStack(ve.c_str(),ve.c_str());
    int draw_index=0;
    for(auto sp :  SP)
    {
      std::string  gopt_tmp = gopt;
      sp->tt->SetLineWidth(1);
      set_color_and_line(sp->tt, draw_index);
      c->cd(index+1);
      sp->tt->Draw(ve.c_str(),selection, gopt.c_str());
      auto h = sp->tt->GetHistogram();
      h->SetName(("h"+sp->title+to_string(index)).c_str());
      h->SetTitle((sp->title).c_str());
      hs->Add(h);
      index++;
      l->AddEntry(sp->tt,sp->title.c_str());
      draw_index++;
    }
    c->cd(index+1);
    hs->Draw();
    l->Draw();
    index++;
  }
}

void compare(std::vector<ScanPoint_t*> SP, const char * varexp_str, const char * selection, std::string gopt)
{
  compare(SP,std::vector<std::string>{varexp_str}, selection, gopt);
}

