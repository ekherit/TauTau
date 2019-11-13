/*
 * =====================================================================================
 *
 *       Filename:  combinator.h
 *
 *    Description: Make all compinations of unique pairs from list of itemes
 *
 *        Version:  1.0
 *        Created:  13.12.2018 15:28:30
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */
#pragma once
#ifndef TAUTAU_COMBINATOR_H
#define TAUTAU_COMBINATOR_H

#define DEBUG_COMBINATIONS
#undef DEBUG_COMBINATIONS

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  make_unique_paris
 *  Description:  Create all compinations of unique pairs from list of items presented
 *  by begin and and iterators. This function allocates memory into result, presented
 *  by type Container1< Container2 < std::pairs<A,A> > >
 *  if list has 2*N items:
 *              Number of pairs is N
 *              All combiantions is  (2*N)!/(2^N *N!)
 *  if list has 2*N-1 items:
 *              Number of pairs is N-1
 *              All combinations is  (2*N)!/(2^N*N!)
 * =====================================================================================
 */
template< class It, class CombinationList> 
inline void make_unique_pairs(It it_begin, It it_end, CombinationList & result);

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  make_unique_pairs
 *  Description:  Make all combination of unique pairs from two arrays of size NA and NB any type
 *  The number of pairs in one combination  is min(NA,NB)
 *  The total number of all such combinations is NA!/(NA-NB)! if NA<=NB
 *  My algorithm use swap of objects, do a lot of memory allocations for the  result
 * =====================================================================================
 */
template < class ItA,   class ItB,   class R>     
inline void  make_unique_pairs
(
 ItA a_begin, ItA a_end, 
 ItB b_begin, ItB b_end, 
 R & result
 );


/* some helper to debug algorithms  */
template< class Ps> inline void print_pairs( const Ps & P);
template< class Cs> inline void  print( const Cs & C);

/**************************************************************************/
/*************** IMPLEMENTATION ******************************************/

#include <vector> //need for realization

//for debug
#include <iostream>
#include <iomanip>

template<class Cont> inline void print_array(std::string name, const Cont & A);
template<class It>   inline void print_array(std::string name, const It b, const It e);

inline void shift(int d){for(int i=0;i<d;++i) std::cout << "     ";};

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  make_unique_pairs_impl
 *  Description:  Make all combination of unique pairs from two arrays of size NA and NB any type
 *  The number of pairs in one combination  is min(NA,NB)
 *  The total number of all such combinations is NA!/(NA-NB)! if NA<=NB
 *  My algorithm use swap of objects, do a lot of memory allocations for the  result
 * =====================================================================================
 */
template < 
class ItA,   //iterator in the first container of objects type A
      class ItB,   //iterator in the second container of object type B
      class CombinationList>     //result contains all combinations of pairs:  
      //   { (a1,b1) (a2,b2) (a3,b3) }
      //   { (a1,b2) (a2,b1) (a3,b3) }
      //   ...
  inline CombinationList  make_unique_pairs_impl
(
 ItA a_begin, ItA a_end, 
 ItB b_begin, ItB b_end 
 ) 
{
  //this is the type of pair's items (A,B)
  typedef typename std::iterator_traits<ItA>::value_type A; 
  typedef typename std::iterator_traits<ItB>::value_type B;

  CombinationList result;

  //the end of recursion
  if( b_begin==b_end || a_begin==a_end) return result;

#ifdef DEBUG_COMBINATIONS

     static int depth=-1;
     depth++;
     static ItA BEGIN_A = a_begin;
     shift(depth);
     std::cout << "depth = " << depth;
     print_array(".  Total array: ", BEGIN_A, a_end);
     std::cout << ",  b_begin == b_end = " << (b_begin == b_end) << std::endl;
     std::cout << std::endl;
     int a_idx=0;
#endif 
  ItA  first   = a_begin;
  ItB  second  = b_begin;
  for(ItA a = a_begin; a!=a_end; ++a)
  {
#ifdef DEBUG_COMBINATIONS
       shift(depth);
       std::cout << "   a_idx = " << a_idx << " ";
       print_array("swap ", a_begin, a_end);
#endif

    if(a != a_begin) std::swap(*a, *a_begin);
#ifdef DEBUG_COMBINATIONS
       print_array(" -> ", a_begin, a_end);
       std::cout << std::endl;
    //R tmp{typename R::value_type({std::pair<A,B> (*a_begin, *b_begin)})};
    shift(depth);
    std::cout << "   depth = " << depth <<" a_idx=" << a_idx << " pair = " << *first <<  *second<< std::endl;
    shift(depth);
    std::cout << "   a_begin = " << &(*a_begin) << " next(a_begin) = " << &(*mynext(a_begin)) << "  b_begin = " << &(*b_begin) <<"  next(b_begin) = " << &(*mynext(b_begin)) << std::endl;
    std::cout << "   a_end = " << &(*a_end) << " b_end = " << &(*b_end) << std::endl;
#endif
    CombinationList local_result = make_unique_pairs_impl<ItA,ItB,CombinationList>(mynext(a_begin), a_end, mynext(b_begin), b_end);
    if(local_result.empty())  local_result.push_back(typename CombinationList::value_type()); 
    //now save the local_result to total result
    //for(const CombinationList::value_type() & pair_list: local_result)
    for(typename CombinationList::const_iterator pair_list_itr=local_result.begin(); pair_list_itr!=local_result.end();++pair_list_itr)
    {
      const typename CombinationList::value_type & pair_list = * pair_list_itr;
      //start new combination
      result.push_back(typename CombinationList::value_type());
      result.back().push_back(std::pair<A,B>(*first,*second));
      std::copy(pair_list.begin(), pair_list.end(),std::back_inserter(result.back()));
    }
#ifdef DEBUG_COMBINATIONS
       shift(depth);
       print_array("   back swap ", a_begin, a_end);
#endif
    if(a!=a_begin) std::swap(*a, *a_begin);//swap back
#ifdef DEBUG_COMBINATIONS
    //if(a!=a_begin) std::swap(*a, *a_begin); //swap back
    print_array(" -> ", a_begin, a_end);
    std::cout << std::endl;
    a_idx++;
#endif
  }
#ifdef DEBUG_COMBINATIONS
     depth--;
#endif
  return result;
} 

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  make_unique_pairs
 *  Description:  Make all combination of unique pairs from two arrays of size NA and NB any type
 *  The number of pairs in one combination  is min(NA,NB)
 *  The total number of all such combinations is NA!/(NA-NB)! if NA<=NB
 *  My algorithm use swap of objects, do a lot of memory allocations for the  result
 * =====================================================================================
 */
template < class ItA,   class ItB,   class CombinationList>     
  inline void  make_unique_pairs
(
 ItA a_begin, ItA a_end, 
 ItB b_begin, ItB b_end, 
 CombinationList & result
 )
{
  result = make_unique_pairs_impl<ItA,ItB,CombinationList>(a_begin,a_end,b_begin,b_end);
}

template< class It, class CombinationList> 
inline void make_unique_pairs(It it_begin, It it_end, CombinationList & result);


  template<class It>
It mynext(const It it)
{
  It new_it = it;
  new_it++;
  return new_it;
}

  template< class It, class CombinationList> 
inline CombinationList make_unique_pairs_impl(It it_begin, It it_end)
{
  //some code for debug
  /*
     static int depth=-1;
     depth++;
     auto shift = [](int d){for(int i=0;i<d;++i) std::cout << "    ";};
     static auto BEGIN_A = it_begin;
     shift(depth);
     print_array("  Total array: ", BEGIN_A, it_end);
     std::cout << std::endl;
     */

  //definitions
  typedef typename std::iterator_traits<It>::value_type A; 
  //      typedef typename  CombinationList::value_type PairList;
  // using PairList = typename  CombinationList::value_type();

  //result.clear();

  CombinationList result; 

  //nothing to do with this case
  if( mynext(it_begin) == it_end ||  it_begin == it_end ) return result; 


  //pointers to pair candidates
  It  first = it_begin;
  It  second = mynext(it_begin);

  /* int idx=0; */
  for(It it = second; it!=it_end; ++it)
  {
    /*
       shift(depth);
       std::cout <<"  depth = " << depth << " idx = " << idx << " *it=" <<  *it << " ";
       std::cout << std::endl;
       */
    //prepare for next iteration
    std::swap(*second, *it);
    CombinationList local_result =  make_unique_pairs_impl<It,CombinationList>(mynext(second), it_end);
    /*
       std::cout << "local_result = " << std::endl;;
       print(local_result);
       */
    if(local_result.empty())  local_result.push_back(typename CombinationList::value_type()); 
    //now save the local_result to total result
    //for(const auto & pair_list: local_result)
    //{
    for(typename CombinationList::const_iterator pair_list_itr=local_result.begin(); pair_list_itr!=local_result.end();++pair_list_itr)
    {
      const typename CombinationList::value_type & pair_list = * pair_list_itr;
      //start new combination
      result.push_back(typename CombinationList::value_type());
      result.back().push_back(std::pair<A,A>(*first,*second));
      std::copy(pair_list.begin(), pair_list.end(),std::back_inserter(result.back()));
    }
    std::swap(*second, *it);//swap back
    /*
       depth--;
       idx++;
       */
  }
  /*
     shift(depth);
     std::cout << " total loop: " << idx << std::endl;
  //print(result);
  */
  return result;
  }


  template< class It, class CombinationList> 
    inline CombinationList make_unique_pairs(It it_begin, It it_end) 
    {
      //definitions
      typedef typename std::iterator_traits<It>::value_type A;
      //using PairList = typename  CombinationList::value_type();

      //separate algorithm for odd and even array size
      size_t distance = std::distance(it_begin,it_end);
      if( distance % 2 == 0) { //even size: just call algorithm
        return make_unique_pairs_impl<It,CombinationList>(it_begin,it_end);
      }
      else { //odd size:  convert array to vector of pointers and add nullptr pointer to
        //get even number and apply common algorithm
        std::vector<A*> V;
        V.resize(distance+1); //allocate a
        int idx=0;
        for(It it = it_begin; it!= it_end; ++it,++idx) V[idx] = &(*it);
        V[distance]=nullptr;
        typedef  typename std::vector<A*>::iterator ItPtr;
        typedef  typename std::vector< std::vector< std::pair<A*, A*> > > CombPtr;
        CombPtr tmp_result_with_nullptr_pairs = make_unique_pairs_impl<ItPtr, CombPtr> (V.begin(),V.end());
        //copy to result and remove nullptr pair
        CombinationList result;
        //for(auto & ptr_pair_list: tmp_result_with_nullptr_pairs)
        for(typename CombPtr::iterator it = tmp_result_with_nullptr_pairs.begin(); it!=tmp_result_with_nullptr_pairs.end(); ++it)
        {
          typedef typename CombPtr::value_type Pairs_t; 
          Pairs_t & ptr_pair_list = *it;
          //PairList res_lst;
          result.push_back(typename CombinationList::value_type());
          //for(const auto & p : ptr_pair_list) 
          for(typename Pairs_t::const_iterator p = ptr_pair_list.begin();  p!=ptr_pair_list.end(); ++p) 
            if(p->first!=nullptr && p->second!=nullptr) 
              result.back().push_back({*(p->first),*(p->second)}); 
        }
        return result;
      }
    }

  /* This function is adoption for deduce the result container type from argumets */
  template< class It, class CombinationList> 
    inline void make_unique_pairs(It it_begin, It it_end, CombinationList & result)
    {
      result = make_unique_pairs<It, CombinationList>(it_begin,it_end);
    }



  template< class Ps>
    inline void print_pairs( const Ps & P)
    {
      for(typename Ps::const_iterator pitr = P.begin();pitr!=P.end(); ++pitr)
        //for( const auto & p: P)
      {
        const typename Ps::value_type & p = *pitr;
        std::cout << p.first<<p.second<< " ";
      }
    };
  template< class Cs>
    inline void print( const Cs & C)
    {
      int idx=1;
      for(typename  Cs::const_iterator it = C.begin();it!=C.end(); ++it)// for(const auto & P: C)
      {
        const typename Cs::value_type & P = *it;
        std::cout << std::setw(20) << idx << "  : ";
        print_pairs(P);
        //for( const auto & p: P)
        //{
        //  std::cout << p.first<<p.second<< " ";
        //}
        std::cout << std::endl;
        idx++;
      }
    };


  template<class Cont>
    void print_array(std::string name, const Cont & A)
    {
      std::cout << name << " = {";
      for(typename Cont::const_iterator it = A.begin(); it!=A.end(); ++it)
      {
        std::cout << *it;
        if( mynext(it)!=A.end() ) std::cout << ',';
      }
      std::cout << "}\n";
    };

  template<class It>
    void print_array(std::string name, const It b, const It e)
    {
      std::cout << name << "{";
      for(It it = b; it!=e; ++it)
      {
        std::cout << *it;
        if( mynext(it)!=e ) std::cout << ',';
      }
      std::cout << "}";
    };

#endif //TAUTAU_COMBINATOR_H
