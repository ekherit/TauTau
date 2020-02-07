#pragma once
/* Fold something */
#include <memory>
#include <TH1.h>
#include "draw_helper.h"

TH1 * fold(std::vector<TH1*> Hs, int Nbin,double xmin, double xmax) {
  //if(Hs.empty()) return nullptr;
  //int Nbin = Hs[0]->GetNbinsX();
  //double xmin,xmax;
  //Hs[0]->GetMinimumAndMaximum(xmin,xmax);
  char buf[1024];
  sprintf(buf,"hfold%d", ++HISTO_INDEX);
  auto his = new TH1F(buf,"his fold", Nbin,xmin,xmax);
  for(auto h : Hs) his->Add(h); 
  return his;
}


/*
template <typename T> 
struct raw_pointer {
  using type = T*;
  raw_pointer() = default;
  type operator()(T * ptr) const { return ptr; }
};

template<typename T>
struct raw_pointer< std::shared_ptr<T> > {
  using type = T*;
  raw_pointer() = default;
  type operator()(std::shared_ptr<T> ptr) const { return ptr.get(); }
};

template<typename T>
using raw_pointer_t = typename raw_pointer<T>::type;

//template <typename T>
//auto get_raw_pointer( T * t ) -> T* { return t; }
//
//template <typename T>
//auto get_raw_pointer( std::shared_ptr<T>  t ) -> T* { return t.get(); }

*/

#include "../ibn/pointer.h"

/* This code can fold any ROOT hystogram 2 */
template<typename Container >
auto fold(const Container & V) -> typename Container::value_type {
  if(V.empty()) return nullptr;
  using T = typename Container::value_type;
  std::string name = "hadd"+std::to_string(++HISTO_INDEX);

  auto it = std::begin(V);
  auto hsum = static_cast<ibn::raw_pointer_t<T>>((*it)->Clone(name.c_str()));
  while(++it != std::end(V)) {
    hsum->Add(ibn::raw_pointer<T>(*it));
  }
  return T(hsum);
};
