/*
 * =====================================================================================
 *
 *       Filename:  comb_pair.C
 *
 *    Description:  Create pair combinations 
 *
 *        Version:  1.0
 *        Created:  13.07.2018 17:09:58
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */
#pragma once

#include <list>
#include <vector>
//template<class T> using combination_t =  std::vector< std::pair<T, T> >;
//template<class T> using combination_t =  std::list< std::pair<T, T> >;
//template<class T> using combination_list_t =  std::vector< std::list< std::pair<T, T> >  >;

template< class T> 
inline void make_combination_list(
    std::vector< std::list< std::pair<T, T> >  > & R, 
    std::list< std::pair<T, T> > & comb, 
    size_t N, 
    T * v)
{
  if(N==1) 
  {
    comb.push_back({v[0],v[1]}); //last pair
    R.push_back(comb); //store this combination to list
    comb.pop_back(); //resotre the combination state
  }
  else
  {
    for(size_t i=1;i<2*N; i++)
    {
      if(1!=i) swap(v[1], v[i]);
      comb.push_back({v[0],v[1]});
      make_combination_list(R,comb,N-1,&v[2]);
      comb.pop_back(); //restore combination state
      if(1!=i) swap(v[1], v[i]);
    }
  }
}


//make combination list from vector of object 
//the combination list will containt pointer to initial objects
template<class T> 
inline std::vector< std::list< std::pair<T*, T*> >  > make_combination_list(const std::vector<T> & V)
{
  std::vector<T*> v(V.size());
  for(size_t i = 0;i<v.size();++i) v[i] = const_cast<T*>(&V[i]);
  std::list< std::pair<T*, T*> >  comb;
  std::vector< std::list< std::pair<T*, T*> >  > comb_list;
  make_combination_list(comb_list, comb, V.size()/2, &v[0]);
  return comb_list;
}

void test_make_combination_list(int N=3)
{
  std::list< std::pair<int, int> >C;
  std::vector< std::list< std::pair<int, int> >  > R;
  std::vector<int> v(2*N);
  for(int i=0;i<2*N;i++) v[i]=i;
  make_combination_list(R,C,N,&v[0]);
  std::map<std::string, int> check;
  for(int i=0;i<R.size();++i)
  {
    std::string s;
    for(std::list<int,int>::iterator p=R[i].begin(); p!=R[i].end(); ++p)
    {
      s+=to_string(p->first)+" "+to_string(p->second)+" ";
    }
    check[s]++;
  }
  
  int count = 0;
  std::cout << setw(10) << "id" << setw(10) << "#" << "           " << "combination" << endl;
  for(auto & m : check)
  {
    std::cout << setw(10) << ++count << setw(10) << m.second << "             " << m.first << endl;
    if(m.second!=1) std::cout << "ERROR wrong combination entry" << std::endl;
  }
}

void test_make_combination_list2(int N=3)
{
  std::vector<int> v(2*N);
  for(int i=0;i<2*N;i++) v[i]=i;
  std::vector< std::list<std::pari<int*, int*> > >  R = make_combination_list(v);
  std::map<std::string, int> check;
  //for (auto & c : R)
  for(int i=0;i<R.size();++i)
  {
    std::string s;
    for(std::list<int*,int*>::iterator p=R[i].begin(); p!=R[i].end(); ++p)
    //for(auto & p : c)
    {
      s+=to_string(*p->first)+" "+to_string(*p->second)+" ";
    }
    check[s]++;
  }
  
  int count = 0;
  std::cout << setw(10) << "id" << setw(10) << "#" << "           " << "combination" << endl;
  for(auto & m : check)
  {
    std::cout << setw(10) << ++count << setw(10) << m.second << "             " << m.first << endl;
    if(m.second!=1) std::cout << "ERROR wrong combination entry" << std::endl;
  }
}
