#ifndef INCLUDED_utils
#define INCLUDED_utils

#ifndef INCLUDED_defs
#include "defs.h"
#endif

#ifndef INCLUDED_std_iostream
#include <iostream>
#endif


#ifndef INCLUDED_std_iterator
#include <iterator>
#endif


std::string strstrip(const char * str);




template <class T>
void printVector(const std::vector<T> & v,
                 bool doAppendNewLine=true,
                 std::ostream & os= std::cout);



template <class T>
std::ostream & operator <<(std::ostream& os, const std::vector<T> & v);



template <class T1, typename T2>
  std::ostream & operator <<(std::ostream& os, const std::pair<T1,T2>  & p);



template <class T1, typename T2>
std::ostream & operator <<(std::ostream& os, const std::map<T1,T2>  & m);



void DEBUGP(double d);



int readWBMscoreFile(const std::string & fname, 
		     const VecStr & idList, 
		     StrPair2Dbl & idpair2wbmScore);




int  readIntxnFile(const std::string & fname, 
		   const std::string & sp,
		   VecStr & idList,
		   Str2Str & id2sp,
		   StrPair2Dbl & idpair2IntxnScore);



int  readBlastScoresFile(const std::string & fname, 
			 const std::string & sp1,
			 const std::string & sp2,
			 VecStr & idList, 
			 Str2Str & id2sp, 
			 StrPair2Dbl & idpair2BlastScore);


void grabSequenceData(const std::string & dataDescFile, 
		      VecStr & spList, 
		      VecStr & idList, 
		      Str2Str & id2sp,
		      StrPair2Dbl & idpair2BlastScore);



void grabSequenceAndInteractionData(const std::string & dataDescFile, 
				    VecStr & spList, 
				    VecStr & idList, 
				    Str2Str & id2sp,
				    StrPair2Dbl & idpair2BlastScore, 
				    StrPair2Dbl & idpair2IntxnScore);



void grabSequenceAndInteractionData(const std::string & dataDescFile, 
				    VecStr & spList, 
				    VecStr & idList, 
				    Str2Str & id2sp,
				    StrPair2Dbl & idpair2BlastScore, 
				    StrPair2Dbl & idpair2IntxnScore,
				    StrPair2Dbl & idpair2wbmScore);



void grabInteractionData(const std::string & dataDescFile, 
			 VecStr & spList, 
			 VecStr & idList, 
			 Str2Str & id2sp,
			 StrPair2Dbl & idpair2IntxnScore);


void grabClusteringData(Str2Int & id2clstr,
			const std::string & clstrFile,
			const VecStr & idList,
			const Str2Str & id2sp);


double getSpeciesDistance(const std::string & sp1, const std::string &sp2);



//normalizeBlastScores follows the approach of Li, Stoecker and Roos
// from Genome Research 13:2178-2189 (2003), (see Fig 2 in that paper)
//
// doWithinSpecies is useful if you don't want to do within-species
//  normalization (e.g. with network-based scores)
//
//  doRemoveLowEntries is not a Li,Stoecker,Roos idea: i put it in to remove
//    entries that are lower than some fraction of the best score. 
void normalizeBlastLikeScores(StrPair2Dbl & idpair2BlastScore,
			      const VecStr & idList,
			      const VecStr & spList,
			      const Str2Str & id2sp,
			      bool doWithinSpecies = true,
			      bool doRescaling=true,
			      bool doRemoveLowEntries=true,
			      double blastLowThreshold = 0.75);



// normalize scores so they are between 0 and 1
//  for non-diagonal entries, i.e. score(i,j) with i!=j, 
//  we do,   score <- score / sqrt(max_score(i) * max_score(j))
//  the diagonal entries are set to 1
void normalizeToZeroOneScores(VecVecIntDblPair & scoreAdjList);


int countNumClusters(const Str2Int & id2clstr);


//clusters which have >1 size
int countNumNontrivialClusters(const Str2Int & id2clstr);


//returns fraction of edges within clusters
double countWithinClusterEdges(const Str2Int & id2clstr,
			       const StrPair2Dbl & idpair2IntxnScore,
			       Int2Int & clstr2withinEdgeCount);


void grabGO_DIPdata(Str2Int & idSetInGO,
		    StrPair2Dbl & gopairSet);



void subsetDataBySetOfIds(StrPair2Dbl & idpair2value,
			  const Str2Int smallIdSet);



VecStr getIdListFromPairs(const StrPair2Dbl & idpairSet);

int getConnectedComponents(const VecVecIntDblPair & scoreAdjList,
			   ListListInt & connectedSets);



void writeAllScoresToFile(const char * outfile,
			  const VecVecIntDblPair & scoreAdjList,
			  const VecStr & idList);


#endif
