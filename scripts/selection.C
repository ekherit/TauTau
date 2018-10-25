/*
 * =====================================================================================
 *
 *       Filename:  draw.C
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

TChain * get_chain(const char * name, const char * newname, std::string title, int run_begin, int run_end)
{
  TChain * chain = new TChain(name, title.c_str());
  for(int run = run_begin; run<=run_end; ++run)
  {
    char buf[1024];
    sprintf(buf, "%d.root",run);
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
};

std::map<std::string, std::string > SelMap;

#include <regex>
const char * make_alias(int channel, const char * templ )
{
  return templ;
}

const double GeV=1.0;
const double MeV=1e-3*GeV;

const double MTAU=1776.86*MeV;
const double MPI=0.13957061*GeV; 
const double MPI0=0.134977*GeV;

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

std::vector<ScanPoint_t> read_data(void)
{
  std::vector<ScanPoint_t> Points;
  
  //Points.push_back({"Point1", 55115,55119,3538.945,0.3, 0, 0});
  //Points.push_back({"Point1", 55120,55126,3539.004,0.116, 0, 0});
  //Points.push_back({"Point1", 55127,55139,3538.646,0.09,0,0});
  //Points.push_back({"Point1", 55143,55155,3539.482,0.110,0,0});
  ///combine points above
  Points.push_back({"Point1", 55116,55155,3539.482,0.110,0,0});

  //Points.push_back({"Point test", 55157,55161,3550.872,0.182,0,0});
  //remove points above this is tune EMS

  //Points.push_back({"Point2", 55162,55177,3552.885,0.093,0,0});
  //Points.push_back({"Point2", 55179,55199,3552.849,0.093,0,0});
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
    Points[i].tt = get_chain("tt", name_tt, Points[i].title, Points[i].begin_run, Points[i].end_run);
    Points[i].gg = get_chain("gg", name_gg, Points[i].title, Points[i].begin_run, Points[i].end_run);
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
    sprintf(run_root_file, "%d.root",run);
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
  //add monte carlo
  int pmc = Points.size();
  Points.push_back({"MC", 55300,55300, 1777*2,0.1,0,0});
  Points[pmc].tt = new TChain("tt","monte carlo tau");
  ((TChain*)Points[pmc].tt)->AddFile("test-sim.root");
  set_alias(Points[pmc].tt,Points[pmc].W);
  Points[pmc].gg = new TNtupleD("ggmc","ggmc","i");
  for(int i=0;i<1000;i++) Points[pmc].gg->Fill();
  return Points;
}

void print(const std::vector<ScanPoint_t> & SPL)
{
  int point_number = 1;
  std::cout << setw(5) << "n/n" << setw(10) << "begin run" << setw(10) <<  "end_run" << setw(15) << "IL, nb^-1" << setw(20) << "W/2-MTAU, MeV" << setw(20) << "dE, MeV" << std::endl;
  for (const auto & p : SPL)
  {
    std::cout << setw(5) << point_number << setw(10) << p.begin_run << setw(10) <<  p.end_run << setw(15) << p.L << setw(20) << p.W/2-MTAU << setw(20) << p.dW/2 << std::endl;
    point_number++;
  }
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
  static int canvas_counter=0;
  char canvas_name[65535];
  sprintf(canvas_name, "tau_sel_%d",canvas_counter++);
  auto c = new TCanvas(canvas_name,"selection",0,0,1920,1080);
  //auto c = new TCanvas();
  c->SetTitle(selection);
  int ny = floor(sqrt(P.size()+1));
  int nx = ny+1;
  c->Divide(nx,ny);
  int i = 1;
  bool add_to_prev_selection=false;
  if(opt=="add") add_to_prev_selection=true;
  for (auto & p : P)
  {
    c->cd(i++);
    int Ntt=p.tt->Draw(varexp,selection,gopt);
    p.Ntt=Ntt;
    if(opt=="clear") p.NttMap.clear();
    if(opt=="add") p.NttMap[selection]=Ntt; 
    char title[1024];
    sprintf(title, "%d events - %s", p.Ntt, p.title.c_str());
    p.tt->GetHistogram()->SetTitle(title);
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
    for(int i=0; i<P.size()-1;++i)
    {
      totalNtt += P[i].Ntt;
      totalNgg+=P[i].Ngg;
      Ltot += P[i].L;
    }
    double sigma_gg = totalNgg/Ltot;
    for(int i=0; i<P.size()-1;++i)
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
    new TCanvas;
    g->SetMarkerStyle(21);
    g->Draw("ap");
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
  for(int i=0; i<P.size()-1;++i)
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
  double sigma_gg = totalNgg/totalL;
  std::ofstream ofs(filename);
  ofs << "#Selection: " << P[0].selection << std::endl;
  ofs << "#Aliases: " << std::endl;
  ofs << "#    e0 = " << P[0].tt->GetAlias("e0") << std::endl;
  ofs << "#    e1 = " << P[0].tt->GetAlias("e1") << std::endl;
  ofs << "#    u0 = " << P[0].tt->GetAlias("u0") << std::endl;
  ofs << "#    u1 = " << P[0].tt->GetAlias("u1") << std::endl;
  ofs << "#    pi0 = " << P[0].tt->GetAlias("pi0") << std::endl;
  ofs << "#    pi1 = " << P[0].tt->GetAlias("pi1") << std::endl;
  for(int i=0; i<P.size()-1;++i)
  {
    int Ntt=0;
    for(auto & item: P[i].NttMap) Ntt+=item.second;
    int Ngg = P[i].Ngg;
    double L = Ngg/(sigma_gg*pow(P[i].W/(2*MTAU),2.0));
    ofs << setw(5) << i <<  setw(15) << P[i].L << "  " << 10 << setw(15) << P[i].W/MeV  << setw(15) << P[i].dW/MeV;
    ofs << setw(10) << 1.256 << " " << setw(10) << 0.019;
    ofs << setw(10) << Ntt << setw(10) << " " << 1 << "  " << Ngg <<  " " << 1 << std::endl;
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
  for(int i=1;i<P.size()-1;++i) chain->Add((TChain*)P[i].tt);
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

auto P = read_data();

void do_all(void)
{
  const char * varexp = "ptem:acop";
  std::vector<std::string> channel =
  {
    "eu_cut"     ,
    "epi_cut"    ,
    "upi_cut"    ,
    "ee_cut"     ,
    "uu_cut"     ,
    "pipi_cut"   ,
    "eK_cut"     ,
    "uK_cut"     ,
    "piK_cut"    ,
    "KK_cut"     ,
    "Xrho_cut"   ,
    "rhorho_cut"
  };
  std::vector<long> Ntt(P.size()-1,0);
  for(int i = 0; i < channel.size(); i++)
  {
    select(P, varexp, channel[i].c_str(), "col" , "add");
    std::cout << setw(10) << channel[i];
    for(int i=0; i < P.size()-1; i++)
    {
      Ntt[i]+=P[i].Ntt;
      std::cout << setw(5) <<  P[i].Ntt;
    }
    std::cout << std::endl;
  }
  std::cout << setw(10) << "total";
  for(int i=0;i<Ntt.size();i++) std::cout  << setw(5) << Ntt[i];
  std::cout << std::endl;
  fit(P);
}
