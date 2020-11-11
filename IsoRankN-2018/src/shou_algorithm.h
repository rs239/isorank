#ifndef INCLUDED_shou_algorithm
#define INCLUDED_shou_algorithm
#endif

#ifndef INCLUDED_defs
#include "defs.h"
#endif

#ifndef INCLUDED_utils
#include "Utils.h"
#endif

#define GAMMA 10 //candidate set
#define EXAMINE_SUBSET 2 //number of branches to be examine
#define BETTA 1e-1
#define SUBSET_THRESHOLD 1e-6 //after the weight of the star becomes small
                             //stop the star modification process
#define FIRST_SPECIES "ecoli"
#define MAX_SPECIES 3

//#define EXPORT_FILENAME "stars_ecoli.txt"
#define READ_CACHED_STARS true

//typedef int Examine_Subset[EXAMINE_SUBSET];
typedef std::vector<int> Examine_Subset;
typedef std::string Star_Branches[MAX_SPECIES][GAMMA];
struct Star
{
  double weight;
  std::string top_protein;
  std::string** branches;
  //Star_Branches branches;
  bool operator <(const Star& other){
    return (weight > other.weight);
  }
};

void doKpartiteMatchingStartWithStar_re(const StrPair2Dbl & idpair2MatchScore, 
            const VecStr & idList,
            const Str2Int & id2idx,
            const Str2Str & id2sp,
            const VecStr & spList,
            double min_primary_fraction,
            double min_secondary_fraction,
            int max_per_species,
            Str2Int & id2clstr, bool starflag, std::string outputClusterFile);

void doKpartiteMatchingStartWithStar_nore(const StrPair2Dbl & idpair2MatchScore, 
            const VecStr & idList,
            const Str2Int & id2idx,
            const Str2Str & id2sp,
            const VecStr & spList,
            double min_primary_fraction,
            double min_secondary_fraction,
            int max_per_species,
            int shouthreadid,
            Str2Int & id2clstr, bool starflag, std::string outputClusterFile);

void export_stars_and_ranks (const std::vector<Star> & stars,
                 const int nsp,
                 const std::vector<int> & ranks,
                 const std::string filename);

void import_stars_and_ranks (std::vector<Star> & stars,
                 const int nsp,
                 std::vector<int> & ranks,
                 const std::string filename);

bool sort_star (const Star& star1, const Star& star2);

double getScore (const std::string & protein1, const std::string & protein2,
         const StrPair2Dbl & idpair2MatchScore);

VecStr this_star_is_ok(Star & star,
               const std::vector<bool> & occupied,
               const std::vector<int> & ranks,
               const StrPair2Dbl & idpair2MatchScore,
               const Str2Int & id2idx,
               int nsp,
               const std::vector<VecStr> & reverse_candidates);
               
VecStr this_star_is_ok_re(Star & star,
               const std::vector<bool> & occupied,
               const std::vector<int> & ranks,
               const StrPair2Dbl & idpair2MatchScore,
               const Str2Int & id2idx,
               int nsp,
               const std::vector<VecStr> & reverse_candidates);               

void try_candidates(Star & star,
            const StrPair2Dbl & idpair2MatchScore,
            const std::vector<bool> & occupied,
            const std::vector<int> & ranks,
            const Str2Int & id2idx,
            int sp1, int sp2);
