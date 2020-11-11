#ifndef INCLUDED_defs
#define INCLUDED_defs

#ifndef INCLUDED_std_algorithm
#include <algorithm>
#endif

#ifndef INCLUDED_std_iostream
#include <iostream>
#endif

#ifndef INCLUDED_std_iomanip
#include <iomanip>
#endif

#ifndef INCLUDED_std_string
#include <string>
#endif

#ifndef INCLUDED_std_map
#include <map>
#endif

#ifndef INCLUDED_std_sstream
#include <sstream>
#endif

#ifndef INCLUDED_std_vector
#include <vector>
#endif

#ifndef INCLUDED_std_list
#include <list>
#endif

#ifndef INCLUDED_std_fstream
#include <fstream>
#endif

#ifndef INCLUDED_std_set
#include <set>
#endif

#ifndef INCLUDED_std_getopt_h
#include <getopt.h>
#endif

#ifndef INCLUDED_std_assert_h
#include <assert.h>
#endif

#ifndef INCLUDED_std_ext_hash_map
//#include <ext/hash_map>
#include <unordered_map>
#endif


#define MAX_LINE_WIDTH 200
#define EPS 0.00000001
#define BLAST_ZERO_EVAL_LOG 500
#define BLAST_EVAL_CUTOFF 10
#define BLAST_ZERO_BIT_SCORE 500
#define BLAST_BIT_SCORE_CUTOFF 60
#define MAX_SCORE 1000000

/*

class dipidhasher : public __gnu_cxx::hash_compare <std::string>
{
public:
  size_t operator() (const std::string& s) const
  {
    size_t h = 0;
    std::string::const_iterator p, p_end;
    int i=0;
    for(p = s.begin()+4, p_end = s.end(),i=0; p != p_end && i < 4; ++p, ++i)
    {
      h = 31 * h + (*p);
    }
    return h;
  }

  bool operator() (const std::string& s1, const std::string& s2) const
  {
    return s1 < s2;
  }
};

*/
struct Node
{
  int index;
  int left;
  int right;
  int parent;
  std::list<int> elements;
};

typedef std::list<Node> ListNode;
typedef std::vector<Node> VecNode;
typedef std::vector<std::string> VecStr;
typedef std::vector<bool> VecBool;
typedef std::vector<int> VecInt;
typedef std::vector<float> VecFloat;
typedef std::vector<VecInt> VecVecInt;
typedef std::vector<VecBool> VecVecBool;
typedef std::vector<std::pair<int,int> > VecIntPair;
typedef std::list<std::pair<int,int> > ListIntPair;
typedef std::vector<double> VecDbl;
typedef std::vector<VecDbl> VecVecDbl;
typedef std::map<std::pair<std::string,std::string>, double> StrPair2Dbl;
typedef std::map<std::pair<int,int>, std::vector<int> > IntPair2VecInt;
typedef std::map<std::string, std::string> Str2Str;
typedef std::pair<std::string,std::string> StrPair;
typedef std::map<StrPair,int> StrPair2Int;
typedef std::map<StrPair, StrPair2Int> StrPair2StrPair2Int;
typedef std::pair<int,int> IntPair;
typedef std::pair<int,double> IntDblPair;
typedef std::map<int, IntDblPair> Int2IntDblPair;
typedef std::pair<double,int> DblIntPair;
typedef std::vector<DblIntPair> VecDblIntPair;
typedef std::map<int,IntPair> Int2IntPair;
typedef std::map<std::string, int> Str2Int;
typedef std::map<std::string, double> Str2Dbl;
typedef std::map<std::string, bool> Str2Bool;
typedef std::map<int, int> Int2Int;
typedef std::map<int, bool> Int2Bool;
typedef std::map<int, double> Int2Dbl;
//typedef __gnu_cxx::hash_map<int, double> Hash_Int2Dbl;
typedef std::unordered_map<int, double> Hash_Int2Dbl;
typedef std::vector<Int2Dbl> VecInt2Dbl;
typedef std::map<std::string,std::vector<double> > Str2VecDbl;
typedef std::map<std::string,std::map<std::string,double> > Str2Str2Dbl;
typedef std::map<std::string,std::map<std::string,int> > Str2Str2Int;
typedef std::vector<std::map<std::string,int> > VecStr2Int;
typedef std::vector<std::vector<std::string> > VecVecStr;
typedef std::vector<std::vector<std::pair<int,double> > > VecVecIntDblPair;
typedef std::map<int, std::vector<std::pair<int,double> > > Int2VecIntDblPair;
typedef std::vector<std::pair<int,double> > VecIntDblPair;
typedef std::list<std::pair<int,double> > ListIntDblPair;
typedef std::vector<std::list<std::pair<int,double> > > VecListIntDblPair;
typedef std::map<std::pair<std::string,std::string>, int> StrPair2Int;
typedef std::map<std::string,std::vector<std::string> > Str2VecStr;
typedef std::map<std::string,std::vector<int> > Str2VecInt;
typedef std::list<std::string> ListStr;
typedef std::map<int,ListStr> Int2ListStr;
typedef std::list<int> ListInt;
typedef std::vector<std::list<int> > VecListInt;
typedef std::list<std::list<int> > ListListInt;
typedef std::map<IntPair,double> IntPair2Dbl;
typedef std::map<IntPair,bool> IntPair2Bool;
typedef std::map<std::string,bool> Str2Bool;
typedef std::map<IntPair,int> IntPair2Int;
typedef std::map<std::string, ListStr> Str2ListStr;
typedef std::map<int, Str2ListStr> Int2Str2ListStr;
typedef std::vector<Str2ListStr> VecStr2ListStr;
typedef std::map<std::string, VecVecInt> Str2VecVecInt;
typedef std::map<int, std::string> Int2Str;
typedef std::set<int> IntSet;
typedef std::set<std::string> StrSet;
#endif
