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

#include "combinator.h"

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


template<class It>
It ibn_next(const It it)
{
  It new_it = it;
  new_it++;
  return new_it;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  make_unique_pairs
 *  Description:  Make all combination of unique pairs from two arrays of size NA and NB any Sortable type 
 *  The number of pairs in one combination  is min(NA,NB)
 *  The total number of all such combinations is NA!/(NA-NB)! if NA<=NB
 *  My algorithm use swap of objects, do a lot of memory allocations for the temporary result
 *  and need memory for the final results.
 *  I could not find the algorithm wich generate unique combinations
 *  So I need to sort pairs in each combinations then sort all combinations
 *  and remove duplicated combinations. If you know better algorithm - your wellcome
 * =====================================================================================
 */
/*
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
   ItA a_begin, ItA a_end,  //range in first container A
   ItB b_begin, ItB b_end,  //range in second container B
   R & result
)
{
  //this is the type of pair's items (A,B)
  //using A = typename std::iterator_traits<ItA>::value_type; 
  //using B = typename std::iterator_traits<ItB>::value_type;
  typedef typename std::iterator_traits<ItA>::value_type A;
  typedef typename std::iterator_traits<ItB>::value_type B;
  typedef typename R::iterator ItR;
  typedef typename R::value_type Comb_t; //combination of pairs
  //The end of recursion.
  if( a_begin==a_end || b_begin==b_end) return;
  //Start of the algorithm
  if( result.empty() ) result.push_back( Comb_t() );
  //Save state of total result. 
  //The reason is:
  //Current state of the result (stars is note one pair):
  //n-1: ***************
  //n  : ******          <- save the state of result
  //then new combination added
  //n:   ******xxxxxxxxx
  //n+1: ******xxxxxxxxx  
  //where x is the new pairs 
  Comb_t begin_combination; 
  std::copy(result.back().begin(),result.back().end(), std::back_inserter(begin_combination));
  //loop over all pairs
  for(ItA a = a_begin; a!=a_end; ++a)
  {
    //I need swap current pair to separate items from
    //rest of array which recursively goes to make_unique_pairs
    if(a != a_begin) std::swap(*a, *a_begin);
    for(ItB b = b_begin; b!=b_end; ++b)
    {
      //this swap is only need to work with
      //any NA>=NB or NB<=NA
      if(b!= b_begin) std::swap(*b,*b_begin);
      R tmp;
      Comb_t cmb;
      //add current pair to the combination
      cmb.push_back(std::pair<A,B> (*a_begin, *b_begin));
      tmp.push_back(cmb);
      //and recursively goes to the remains of the arrays
      make_unique_pairs(ibn_next(a_begin), a_end, ibn_next(b_begin), b_end, tmp);
      //save the result to final result
      for(ItR it = tmp.begin(); it!=tmp.end(); ++it)
      {
        std::copy(it->begin(), it->end(),std::back_inserter(result.back()));
        //suppress new combinations at the most end of loops over pairs
        if(ibn_next(a) != a_end || ibn_next(b) != b_end || ibn_next(it)!=tmp.end()) result.push_back(begin_combination); 
      }
      if(b!= b_begin) std::swap(*b,*b_begin);//swap back item of B array
    }
    if(a != a_begin) std::swap(*a, *a_begin); //swap back items of A array
  }
  //sort each combination
  for(ItR r=result.begin(); r!=result.end(); ++r) std::sort(r->begin(),r->end());
  //sort final result
  std::sort(result.begin(),result.end());
  //remove duplicated combinations
  result.erase(std::unique(result.begin(),result.end()), result.end());
}


//  Some helper function to print the test result 
template< class Array>
void print_array(std::string name, const Array & A)
{
  std::cout << name << " = {";
  for(typename Array::const_iterator it = A.begin(); it!=A.end(); ++it)
  {
    std::cout << *it;
    if( ibn_next(it)!=A.end() ) std::cout << ',';
  }
  std::cout << "}\n";
};

template< class Ps>
void print_pairs( const Ps & P)
{
  //for( const auto & p: P)
  for( typename Ps::const_iterator p = P.begin(); p!=P.end(); ++p)
  {
    std::cout << p->first<<p->second<< " ";
  }
};

template< class Cs>
void print( const Cs & C)
{
  int idx=1;
  for(typename Cs::const_iterator P=C.begin();P!=C.end();++P)
  {
    std::cout << std::setw(20) << idx << "  : ";
    print_pairs(*P);
    std::cout << std::endl;
    idx++;
  }
};
*/
//test for algorithm
inline void test_make_unique_pairs(int NA=1, int NB=1)
{
  std::vector<int>  A;
  std::vector<char> B;
  for( int i=0;i<NA;++i) A.push_back(i+1);
  for( int i=0;i<NB;++i) B.push_back('a'+i);
  print_array("A",A);
  print_array("B",B);
  std::vector< std::vector < std::pair<int, char> > > R;
  make_unique_pairs(A.begin(), A.end(), B.begin(), B.end(), R);
  std::cout << "The result of algorithm: " << R.size() << " combinations:" << std::endl;
  print(R);
}


