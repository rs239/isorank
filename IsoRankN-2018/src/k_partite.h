#ifndef INCLUDED_k_partite
#define INCLUDED_k_partite

#ifndef INCLUDED_defs
#include "defs.h"
#endif

#ifndef INCLUDED_utils
#include "Utils.h"
#endif


#include <math.h>

struct ValIdxPair {
  double val;
  int first;
  int second;
};

typedef std::vector<ValIdxPair> VecValIdxPair;



void doKpartiteMatching(const StrPair2Dbl & idpair2MatchScore, 
			const VecStr & idList,
			const Str2Int & id2idx,
			const Str2Str & id2sp,
			const VecStr & spList,
			double min_primary_fraction,
			double min_secondary_fraction,
			int max_per_species,
			Str2Int & id2clstr);


bool isOKasSecondaryNodeInClstr(int idx,
				int numSpInPrimaryClstr,
				const IntSet & currNodesInClstr,
				const Int2VecIntDblPair & idx2neighs,
				const Int2Str & idx2sp,
				const StrPair2Dbl & sppair2clstrTopScore,
				double min_secondary_fraction,
				int max_per_species);


void addSecondaryNodes(const VecValIdxPair & scores, 
		       IntSet & currNodesInClstr,
		       const Int2VecIntDblPair & idx2neighs, 
		       VecBool & haveProcessedIdx, 
		       VecBool & ignoreScore,
		       const Int2Str & idx2sp, 
		       const VecStr & spList,
		       Int2Int & idx2clstr,
		       int currClstrIdx,
		       double min_secondary_fraction,
		       int max_per_species);


void  processNewScore(int index,
		      const VecValIdxPair & scores, 
		      const Int2VecIntDblPair & idx2neighs, 
		      VecBool & haveProcessedIdx, 
		      VecBool & ignoreScore,
		      const Int2Str & idx2sp, 
		      const VecStr & spList,
		      Int2Int & idx2clstr,
		      int & currClstrIdx,
		      double min_primary_fraction,
		      double min_secondary_fraction,
		      int max_per_species);


void updateNeighFreqAndScore(int idx, 
			     const Int2Str & idx2sp, 
			     const Str2Bool & sp2seen, 
			     const Int2VecIntDblPair & idx2neighs, 
			     Int2IntDblPair & neighIdx2freq_and_score);


void  removeNeighbors(int idx, 
		      const Int2Str & idx2sp, 
		      Str2Bool & sp2seen, 
		      Int2IntDblPair & neighIdx2freq_and_score);



void ignoreAllEdgesInvolvingIndex(int idx, 
				  const VecValIdxPair & scores, 
				  VecBool & ignoreScore);


Int2IntDblPair::const_iterator findBestNeighbor(
			   const Int2IntDblPair & neighIdx2freq_and_score, 
			   int clstrSize,
			   const Str2Bool & sp2seen, 
			   const Int2Str & idx2sp);

#endif 
