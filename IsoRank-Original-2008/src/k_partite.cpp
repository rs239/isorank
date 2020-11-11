#include "k_partite.h"

struct GreaterScore
{
  inline bool operator()(const ValIdxPair & p1, const ValIdxPair & p2)
  {
    return p1.val > p2.val;
  }
};




void doKpartiteMatching(const StrPair2Dbl & idpair2MatchScore, 
			const VecStr & idList,
			const Str2Int & id2idx,
			const Str2Str & id2sp,
			const VecStr & spList,
			double min_primary_fraction,
			double min_secondary_fraction,
			int max_per_species,
			Str2Int & id2clstr)
{
  typedef StrPair2Dbl::const_iterator It1;
  typedef Str2Int::const_iterator It2;
  typedef Int2Int::const_iterator It3;

  DEBUGP(601);

  Int2Int idx2clstr;
  int currClstrIdx = 0;

  Int2Str idx2sp;
  int maxIdx = 0;
  for (It2 it2=id2idx.begin(); it2!=id2idx.end(); it2++)
    {
      idx2sp[it2->second] = (id2sp.find(it2->first))->second;
      if (maxIdx < it2->second)
	maxIdx = it2->second;
    }

  DEBUGP(602);

  assert(idx2sp.size() == maxIdx+1);
  VecBool haveProcessedIdx(maxIdx+1, false);


  VecValIdxPair scores;
  Int2VecIntDblPair idx2neighs;
  
  VecIntDblPair V0;
  IntDblPair P0;
  for (It1 it1=idpair2MatchScore.begin(); it1 != idpair2MatchScore.end();it1++)
    {
      //DEBUGP(603);

      ValIdxPair p;
      p.val = it1->second;
      p.first = id2idx.find(it1->first.first)->second;
      p.second = id2idx.find(it1->first.second)->second;

      scores.push_back(p);

      if (idx2neighs.find(p.first)==idx2neighs.end())
	{
	  idx2neighs[p.first] = V0;
	}
      P0.first = p.second;
      P0.second = p.val;
      idx2neighs[p.first].push_back(P0);

      if (idx2neighs.find(p.second)==idx2neighs.end())
	{
	  idx2neighs[p.second] = V0;
	}
      P0.first = p.first;
      P0.second = p.val;
      idx2neighs[p.second].push_back(P0);
    }

  DEBUGP(604);

  VecBool ignoreScore(scores.size(), false);
  
  GreaterScore GS;
  std::sort(scores.begin(), scores.end(), GS);

  for (int i=0; i < scores.size(); i++)
    {
      //DEBUGP(605);
      if (ignoreScore[i] ||
	  haveProcessedIdx[scores[i].first] || 
	  haveProcessedIdx[scores[i].second])
	continue;

      processNewScore(i,
		      scores, idx2neighs, 
		      haveProcessedIdx, ignoreScore,
		      idx2sp, spList,
		      idx2clstr,
		      currClstrIdx,
		      min_primary_fraction,
		      min_secondary_fraction,
		      max_per_species);
		      
    }

  DEBUGP(606);
  id2clstr.clear();
  for (It2 it2=id2idx.begin(); it2!=id2idx.end(); it2++)
    {
      It3 it3 = idx2clstr.find(it2->second);
      if (it3==idx2clstr.end())
	{
	  currClstrIdx++;
	  id2clstr[it2->first] = currClstrIdx;
	}
      else
	{
	  id2clstr[it2->first] = it3->second;
	}
    }
  DEBUGP(607);

}




double getScore(int idx1, int idx2, const Int2VecIntDblPair & idx2neighs)
{
  double v = -1;
  const VecIntDblPair & V = idx2neighs.find(idx1)->second;
  for (int i=0; i < V.size(); i++)
    {
      if (V.at(i).first == idx2)
	{
	  v = V.at(i).second;
	  break;
	}
    }
  return v;
}




bool isOKasSecondaryNodeInClstr(int idx,
				int numSpInPrimaryClstr,
				const IntSet & currNodesInClstr,
				const Int2VecIntDblPair & idx2neighs,
				const Int2Str & idx2sp,
				const StrPair2Dbl & sppair2clstrTopScore,
				double min_secondary_fraction,
				int max_per_species)
{
  typedef IntSet::const_iterator It1;
  typedef Int2Str::const_iterator It2;

  StrSet spSeenSet;

  //make this lower to increase the number of secondary nodes. max val = 1
  double min_fraction = min_secondary_fraction;

  const VecIntDblPair & V = idx2neighs.find(idx)->second;
  std::string currSp = idx2sp.find(idx)->second;

  int numIdxsWithCurrSp = 0;
  for (It1 it1=currNodesInClstr.begin();it1 != currNodesInClstr.end(); it1++)
    {
      if (idx2sp.find(*it1)->second == currSp)
	numIdxsWithCurrSp++;
    }

  if (numIdxsWithCurrSp >= max_per_species)
    {
      std::cerr << "sp: " << currSp 
		<< " numIdxsWithCurrSp: " << numIdxsWithCurrSp << std::endl;
      return false;
    }

  for (int i=0; i < V.size(); i++)
    {
      int neighidx = V.at(i).first;
      double v =  V.at(i).second;
      std::string neighSp = idx2sp.find(neighidx)->second;

      if (currNodesInClstr.find(neighidx) == currNodesInClstr.end())
	continue;
      
      StrPair s(currSp, neighSp);
      if (min_fraction*(sppair2clstrTopScore.find(s)->second) > v)
	continue;
      
      spSeenSet.insert(neighSp);
    }

  // -1 means that 'idx' has high-scoring connections to all other sp 
  if (spSeenSet.size() >0 && spSeenSet.size() >= (numSpInPrimaryClstr-1))
    return true;
  else
    return false;
}
  




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
		       int max_per_species)
{
  typedef Int2IntDblPair::iterator It4;
  typedef IntSet::const_iterator It5;

  //primary cluster (i.e. primary nodes) have at most one node per sp
  int numSpInPrimaryClstr = currNodesInClstr.size();

  //make a map of all nodes that are neighbors of one or more nodes in the set
  //make a node2freq map so that nodes are mapped to the number of 
  //   neigbors they have 
  Str2Bool sp2seen;
  for (int i=0; i < spList.size(); i++)
    {
      sp2seen[spList.at(i)] = false; 
    }

  Int2IntDblPair neighs2freq_and_score;
  for (It5 it5=currNodesInClstr.begin(); it5 != currNodesInClstr.end(); it5++)
    {
      updateNeighFreqAndScore(*it5, idx2sp, sp2seen, 
			      idx2neighs, neighs2freq_and_score);
    }

  int initialNodesInCluster = currNodesInClstr.size();

  //filter the node2freq to remove some
  for (It4 it4= neighs2freq_and_score.begin();
       it4 != neighs2freq_and_score.end();
       )
    {
      //each duplicate node should be in at least 2 species or missing 
      // from at most 1 species (whichever is max)
      int minFreq = std::max(2, initialNodesInCluster-2);
      if (it4->second.first < minFreq)
	{
	  It4 it_a = it4;
	  it4++;
	  neighs2freq_and_score.erase(it_a);
	}
      else
	{
	  it4++;
	}

    }

  //for each node in node2freq, evaluate it individually to see if it should 
  //   be added and if so, add it
  StrPair2Dbl sppair2clstrTopScore;
  for (It5 it5=currNodesInClstr.begin(); it5 != currNodesInClstr.end(); it5++)
    {
      It5 it5b =it5;
      it5b++;
      while (it5b != currNodesInClstr.end())
	{
	  std::string sp1,sp2;
	  sp1 = idx2sp.find(*it5)->second;
	  sp2 = idx2sp.find(*it5b)->second;

	  double v = getScore(*it5,*it5b, idx2neighs);
	  if (v==-1)
	    {
	      std::cerr << "this sp-pair score didn't exist initially: " 
			<< sp1 << " " << sp2 << std::endl;
	    }

	  StrPair s1(sp1,sp2);
	  StrPair s2(sp2,sp1);
	  sppair2clstrTopScore[s1] = v;
	  sppair2clstrTopScore[s2] = v;
	  
	  it5b++;
	}
    }
  
  for (It4 it4=neighs2freq_and_score.begin();
       it4 != neighs2freq_and_score.end();
       it4++)
    {
      int idx = it4->first;
      if ( isOKasSecondaryNodeInClstr(idx,
				      numSpInPrimaryClstr,
				      currNodesInClstr,
				      idx2neighs,
				      idx2sp,
				      sppair2clstrTopScore, 
				      min_secondary_fraction,
				      max_per_species))
	{
	  idx2clstr[idx] = currClstrIdx;
	  ignoreAllEdgesInvolvingIndex(idx, scores, ignoreScore);
	  haveProcessedIdx[idx] = true;	  
	  currNodesInClstr.insert(idx);
	}
    }
}





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
		      int max_per_species)
{
  typedef Int2Str::const_iterator It1;
  typedef Int2VecIntDblPair::const_iterator It2;
  typedef VecIntDblPair::const_iterator It3;
  typedef Int2IntDblPair::const_iterator It4;
  typedef IntSet::const_iterator It5;

  int idx1= scores.at(index).first;
  int idx2 = scores.at(index).second;

  currClstrIdx ++;

  idx2clstr[idx1] = currClstrIdx;
  idx2clstr[idx2] = currClstrIdx;

  double thisClstrMaxScore = scores.at(index).val;
  //DEBUGP(6051);

  ignoreAllEdgesInvolvingIndex(idx1, scores, ignoreScore);
  ignoreAllEdgesInvolvingIndex(idx2, scores, ignoreScore);
  haveProcessedIdx[idx1] = true;
  haveProcessedIdx[idx2] = true;
  ignoreScore[index] = true;

  Str2Bool sp2seen;
  for (int i=0; i < spList.size(); i++)
    {
      sp2seen[spList.at(i)] = false;
    }
  
  std::string sp1,sp2;
  It1 sp1it = idx2sp.find(idx1);
  It1 sp2it = idx2sp.find(idx2);
  assert (sp1it != idx2sp.end() && sp2it != idx2sp.end());

  sp1 = sp1it->second;
  sp2 = sp2it->second;

  if (sp1==sp2)
    return;

  sp2seen[sp1] = true;
  sp2seen[sp2] = true;

  Int2IntDblPair neighIdx2freq_and_score;
  IntSet thisClstrIdxs;

  //DEBUGP(6052);
  updateNeighFreqAndScore(idx1, idx2sp, sp2seen, idx2neighs, 
			  neighIdx2freq_and_score);

  thisClstrIdxs.insert(idx1);

  updateNeighFreqAndScore(idx2, idx2sp, sp2seen, idx2neighs, 
			  neighIdx2freq_and_score);

  thisClstrIdxs.insert(idx2);

  //DEBUGP(6053);
  while (neighIdx2freq_and_score.size()!=0)
    {

      std::cerr << "num neighs: " << neighIdx2freq_and_score.size() 
                << std::endl;
      //DEBUGP(6054);

      It4 best_neigh_it = findBestNeighbor(neighIdx2freq_and_score,
					   thisClstrIdxs.size(),
					   sp2seen, idx2sp);
      //DEBUGP(6055);

      int idx3 = -1;
      if (best_neigh_it != neighIdx2freq_and_score.end())
	{
	  //DEBUGP(6056);

	  //it might be worth modifying the removeNeighbors() code
	  // so that it can be moved outside the then{} block and can
	  // be made common to both the then{} and else{} blocks

          idx3 = best_neigh_it->first;
	  if (best_neigh_it->second.second >= (min_primary_fraction * 
					       thisClstrMaxScore))
	    {
	      if (thisClstrMaxScore < best_neigh_it->second.second)
		thisClstrMaxScore =  best_neigh_it->second.second;

	      idx2clstr[idx3] = currClstrIdx;
	      sp2seen[idx2sp.find(idx3)->second] = true;
	      updateNeighFreqAndScore(idx3, idx2sp, sp2seen, idx2neighs, 
				      neighIdx2freq_and_score);
	      thisClstrIdxs.insert(idx3);

	      ignoreAllEdgesInvolvingIndex(idx3, scores, ignoreScore);
	      haveProcessedIdx[idx3] = true;
	      removeNeighbors(idx3, idx2sp, sp2seen, neighIdx2freq_and_score);
	    }
	  else
	    {
	      break;
	    }
	}
    }
  if (max_per_species > 1)
    {
      //add duplicate nodes 
      addSecondaryNodes(scores, 
			thisClstrIdxs,
			idx2neighs, 
			haveProcessedIdx, 
			ignoreScore,
			idx2sp, 
			spList,
			idx2clstr,
			currClstrIdx,
			min_secondary_fraction,
			max_per_species);
    }
}

		      



void updateNeighFreqAndScore(int idx, 
			     const Int2Str & idx2sp, 
			     const Str2Bool & sp2seen, 
			     const Int2VecIntDblPair & idx2neighs, 
			     Int2IntDblPair & neighIdx2freq_and_score)
{
  typedef Str2Bool::const_iterator It1;
  typedef Int2VecIntDblPair::const_iterator It2;
  typedef VecIntDblPair::const_iterator It3;
  typedef Int2IntDblPair::const_iterator It4;
  
  It2 it2 = idx2neighs.find(idx);
  if (it2 == idx2neighs.end())
    return;

  for (It3 it3= it2->second.begin(); it3!= it2->second.end(); it3++)
    {
      int neighidx = it3->first;
      double val = it3->second;
	  
      It1 it1 = sp2seen.find(idx2sp.find(neighidx)->second);
      if (it1 == sp2seen.end()  || it1->second == true)
	continue;
      
      if (idx2sp.find(idx)->second == idx2sp.find(neighidx)->second)
	continue;
	  
      if (neighIdx2freq_and_score.find(neighidx) != 
	  neighIdx2freq_and_score.end())
	{
	  neighIdx2freq_and_score[neighidx].first += 1;
	  neighIdx2freq_and_score[neighidx].second = 
	    std::min(neighIdx2freq_and_score[neighidx].second, val);
	}
      else
	{
	  IntDblPair p(1,val);
	  neighIdx2freq_and_score[neighidx] = p;
	}
    }      
  

}
  


void  removeNeighbors(int idx, 
		      const Int2Str & idx2sp, 
		      Str2Bool & sp2seen, 
		      Int2IntDblPair & neighIdx2freq_and_score)
{
  typedef Int2IntDblPair::iterator It1;
  std::string sp = idx2sp.find(idx)->second;
  sp2seen[sp] = true;

  for (It1 it1=neighIdx2freq_and_score.begin(); 
       it1 != neighIdx2freq_and_score.end() ; )
    {
      //DEBUGP(1101);
      if (sp2seen[idx2sp.find(it1->first)->second])
	{
	  It1 it1a = it1;
	  it1a++;
	  neighIdx2freq_and_score.erase(it1);
	  it1 = it1a;
	}
      else
	{
	  it1++;
	}
    }  
  //DEBUGP(1102);
}




void ignoreAllEdgesInvolvingIndex(int idx, 
				  const VecValIdxPair & scores, 
				  VecBool & ignoreScore)
{
  for (int i=0; i < scores.size(); i++)
    {
      if (scores.at(i).first == idx || scores.at(i).second == idx)
	{
	  ignoreScore[i] = true;
	}
    }
}



Int2IntDblPair::const_iterator findBestNeighbor(
			   const Int2IntDblPair & neighIdx2freq_and_score, 
			   int clstrSize,
			   const Str2Bool & sp2seen, 
			   const Int2Str & idx2sp)
{
  typedef Int2IntDblPair::const_iterator It1;

  int min_freq = static_cast<int>(floor(EPS + (clstrSize/2.0)) + 1);

  It1 result_it = neighIdx2freq_and_score.end();
  It1 end_it = neighIdx2freq_and_score.end();
  for (It1 it1=neighIdx2freq_and_score.begin();
       it1 != end_it;
       it1++)
    {
      //DEBUGP(901);
      int freq = it1->second.first;
      double val = it1->second.second;
      
      /*if (result_it == end_it ||
	  (freq > result_it->second.first) ||
	  (freq == result_it->second.first && val > result_it->second.second))
	result_it = it1;
      */

      if (result_it == end_it ||
	  (freq >= min_freq && val > result_it->second.second))
	result_it = it1;

    }

  return result_it;
}
