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
#include <ostream>

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
    for(std::list<std::pair<int,int> >::iterator p=R[i].begin(); p!=R[i].end(); ++p)
    {
      std::ostringstream os;
      os << p->first << " " << p->second << " ";
      //s+=to_string(p->first)+" "+to_string(p->second)+" ";
      s+=os.str();
    }
    check[s]++;
  }
  
  int count = 0;
  std::cout << setw(10) << "id" << setw(10) << "#" << "           " << "combination" << endl;
  for(std::map<std::string, int>::iterator it=check.begin(); it!=check.end(); ++it)
  {
    std::cout << setw(10) << ++count << setw(10) << it->second << "             " << it->first << endl;
    if(it->second !=1 ) std::cout << "ERROR wrong combination entry" << std::endl;
  }
}

void test_make_combination_list2(int N=3)
{
  std::vector<int> v(2*N);
  for(int i=0;i<2*N;i++) v[i]=i;
  std::vector< std::list<std::pair<int*, int*> > >  R = make_combination_list(v);
  std::map<std::string, int> check;
  //for (auto & c : R)
  for(int i=0;i<R.size();++i)
  {
    std::string s;
    for(std::list< std::pair< int*,int*> >::iterator p=R[i].begin(); p!=R[i].end(); ++p)
    //for(auto & p : c)
    {
      std::ostringstream os;
      os << *(p->first) << " " << *(p->second) << " ";
      //s+=to_string(p->first)+" "+to_string(p->second)+" ";
      s+=os.str();
      //s+=to_string(*p->first)+" "+to_string(*p->second)+" ";
    }
    check[s]++;
  }
  
  int count = 0;
  std::cout << setw(10) << "id" << setw(10) << "#" << "           " << "combination" << endl;
  //for(auto & m : check)
  for(std::map<std::string, int>::iterator m=check.begin(); m!=check.end(); ++m)
  {
    std::cout << setw(10) << ++count << setw(10) << m->second << "             " << m->first << endl;
    if(m->second!=1) std::cout << "ERROR wrong combination entry" << std::endl;
  }
}
#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>
#include <iomanip>

template< class Ps>
void print_pairs( const Ps & P)
{
  for( const auto & p: P)
  {
    std::cout << p.first<<p.second<< " ";
  }
};
template< class Cs>
void print( const Cs & C)
{
  int idx=1;
  for(const auto & P: C)
  {
    std::cout << std::setw(20) << idx << "  : ";
    print_pairs(P);
    std::cout << std::endl;
    idx++;
  }
};

template
< 
    class ItA,   //iterator in the first container of objects type A
    class ItB,   //iterator in the second container of object type B
    class R      //result contains all combinations of pairs:  
    //   { (a1,b1) (a2,b2) (a3,b3) }
    //   { (a1,b2) (a2,b1) (a3,b3) }
    //   ...
>
void  make_unique_pairs
(
   ItA a_begin, ItA a_end, 
   ItB b_begin, ItB b_end, 
   R & result
)
{
  //this is the type of pair's items (A,B)
  using A = typename std::iterator_traits<ItA>::value_type; 
  using B = typename std::iterator_traits<ItB>::value_type;
  //the end of recursion
  if( a_begin==a_end || b_begin==b_end) return;
  if( result.empty() ) result.push_back( typename R::value_type() ); //initial start
  typename R::value_type begin_combination; 
  std::copy(result.back().begin(),result.back().end(), std::back_inserter(begin_combination));
  for(auto a = a_begin; a!=a_end; ++a)
  {
    if(a != a_begin) std::swap(*a, *a_begin);
    for(auto b = b_begin; b!=b_end; ++b)
    {
      if(b!= b_begin) std::swap(*b,*b_begin);
      R tmp{typename R::value_type({std::pair<A,B> (*a_begin, *b_begin)})};
      make_unique_pairs(next(a_begin), a_end, next(b_begin), b_end, tmp);
      for(auto it = tmp.begin(); it!=tmp.end(); ++it)
      {
        std::copy(it->begin(), it->end(),std::back_inserter(result.back()));
        if(std::next(a) != a_end || std::next(b) != b_end || std::next(it)!=tmp.end()) result.push_back(begin_combination); 
      }
      if(b!= b_begin) std::swap(*b,*b_begin);//swap back
    }
    if(a != a_begin) std::swap(*a, *a_begin); //swap back
  }
  for(auto & r : result) std::sort(r.begin(),r.end());
  std::sort(result.begin(),result.end());
  result.erase(std::unique(result.begin(),result.end()), result.end());
}

inline void test_make_pairs(int NA=1, int NB=1)
{
  std::vector<int>  A={1};
  std::vector<char> B={'a'};
  for( int i=1;i<NA;++i) A.push_back(i+1);
  for( int i=1;i<NB;++i) B.push_back('a'+i);
  auto print_array = [](std::string name="", const auto & A)
  {
    std::cout << name << " = {";
    for(auto it = A.begin(); it!=A.end(); ++it)
    {
      std::cout << *it;
      if( std::next(it)!=A.end() ) std::cout << ',';
    }
    std::cout << "}\n";
  };
  print_array("A",A);
  print_array("B",B);
  std::vector< std::vector < std::pair<int, char> > > R;
  make_unique_pairs(A.begin(), A.end(), B.begin(), B.end(), R);
  std::cout << "The result of algorithm: " << R.size() << " combinations:" << std::endl;
  print(R);
}



