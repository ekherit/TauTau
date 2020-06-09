#include <iostream>
#include <string>
//void sys(std::string filename) {
//  vector<double> p;
//  vector<double> N;
//  vector<double> M;
//  ifstream ifs(filename);
//  vector<double> tmps(3);
//  while(ifs >> tmps[0]){
//    for(size_t i=1; i<tmps.size(); ++i) ifs>>tmps[i];
//    for(auto & x : tmps) cout << x << " ";
//    std::cout << "\n";
//    p.push_back(tmps[0]);
//    N.push_back(tmps[1]);
//    M.push_back(tmps[2]);
//  };
//
//  TGraph * pg = new TGraph;
//  TGraph * Ng = new TGraph;
//  for(int i=0;i<p.size();++i) {
//    double dM = M[i]-M[0];
//    pg->SetPoint(i, p[i], dM);
//    Ng->SetPoint(i, (N[i]-N[0])/N[0]*100.0, dM);
//  }
//  pg->Sort();
//  new TCanvas;
//  Ng->Draw("a*");
//  new TCanvas;
//  pg->Draw("al*");
//}

void sys(void) {
}

void sys(std::string filename, std::string xaxis) {
  vector<double> p;
  vector<double> N;
  vector<double> M;
  ifstream ifs(filename);
  vector<double> tmps(3);
  while(ifs >> tmps[0]){
    for(size_t i=1; i<tmps.size(); ++i) ifs>>tmps[i];
    for(auto & x : tmps) cout << x << " ";
    std::cout << "\n";
    p.push_back(tmps[0]);
    N.push_back(tmps[1]);
    M.push_back(tmps[2]);
  };

  TGraph * pg = new TGraph;
  TGraph * Ng = new TGraph;
  for(int i=0;i<p.size();++i) {
    double dM = M[i]-M[0];
    pg->SetPoint(i, p[i], dM);
    Ng->SetPoint(i, (N[i]-N[0])/N[0]*100.0, dM);
  }
  pg->Sort();
  new TCanvas;
  Ng->Sort();
  Ng->Draw("al*");
  //Ng->Fit("pol1");
  auto fun = new TF1("fun_bias","[0]+[1]*x",-3*100./sqrt(N[0]), +3*100./sqrt(N[0]));
  Ng->Fit("fun_bias","R");
  gStyle->SetOptFit();
  double p0 = fun->GetParameter(0);
  double p1 = fun->GetParameter(1);
  double systematics = hypot(p0,  p1*100./sqrt(N[0]));
  std::cout << "Systematics = " << systematics*1000.0 << " keV " << std::endl;
  new TCanvas;
  pg->Draw("al*");
  pg->GetXaxis()->SetTitle(xaxis.c_str());
  pg->GetYaxis()->SetTitle("#Delta M, keV");
}
