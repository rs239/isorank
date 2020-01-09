#ifndef  INCLUDED_utils
#include "Utils.h"


#include <math.h>
#include <cstring>

std::string strstrip(const char * str)
{
  int len = std::strlen(str);
  if (len==0)
    {
      return std::string();
    }

  int a1 = 0;
  while (!isspace(str[a1]))
    {
      a1++;
    }

  int a2 = len-1;
  while (!isspace(str[a2]))
    {
      a2--;
    }
  if (a1 <= a2)
    {
      return std::string(str+a1, str+a2+1);
    }
  else
    {
      return std::string();
    }
}





template <class T>
void printVector(const std::vector<T> & v,
                 bool doAppendNewLine,
                 std::ostream & os)
{
  os << "[";
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os," "));
  os << "]";
  if (doAppendNewLine)
    os << "\n";
}





template <class T>
std::ostream & operator <<(std::ostream& os, const std::vector<T> & v)
{
  os << "[";
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os," "));
  os << "] ";
  return os;
}





template <class T1, typename T2>
std::ostream & operator <<(std::ostream& os, const std::pair<T1,T2>  & p)
{
  os << "(" << p.first << "," << p.second << ")";
  return os;
}




/* the new g++ compiler is bitching about this one... dunno why
template <class T1, typename T2>
std::ostream & operator << (std::ostream& os, const std::map<T1,T2>  & m)
{
  std::map<T1, T2>::const_iterator it1 = m.begin();

  os << "{";
  for (it1 = m.begin(); it1 != m.end(); it1++)
    {
      os << it1->first << ":" << it1->second << ", ";
    }
  os << "} ";
  return os;
}
*/



void DEBUGP(double d)
{
  std::cerr << "Flag: " << std::setprecision(4) << d << std::endl;
}





int readWBMscoreFile(const std::string & fname, 
		     const VecStr & idList, 
		     StrPair2Dbl & idpair2wbmScore)
{
  typedef VecStr::const_iterator It1;
  typedef Str2Int::const_iterator It2;

  idpair2wbmScore.clear();

  Str2Int idSet;
  for (It1 it1 = idList.begin(); it1 != idList.end(); it1++)
    {
      idSet[*it1] = 1;
    }

  int numVals = 0;
  std::ifstream ifs(fname.c_str());
  
  char line[200];
  std::string id1, id2;
  double val = 0;

  while (ifs.getline(line,200-1))
    {      
      std::istringstream istrstr(line);
      istrstr >> id1 >> id2 >> val;
      assert( idSet.find(id1) != idSet.end() && 
	      idSet.find(id2) != idSet.end());
      if (id1 < id2)
	idpair2wbmScore[StrPair(id1,id2)] = val;
      else
	idpair2wbmScore[StrPair(id2,id1)] = val;

      numVals++;
    }
  return numVals;
}






int  readIntxnFile(const std::string & fname, 
		   const std::string & sp,
		   VecStr & idList,
		   Str2Str & id2sp,
		   StrPair2Dbl & idpair2IntxnScore)
{
  
  typedef Str2Str::const_iterator It;
  
  int numVals = 0;

  std::cerr << "intxn file: " << fname << std::endl;

  std::ifstream ifs(fname.c_str());

  std::string id1, id2;
  bool haveConfVal = false;
  double confVal = 0;
  
  std::string s1,s2,s3;
  ifs >> s1 >> s2 >> s3;
  assert(!ifs.eof());

  assert(s1=="INTERACTOR_A" && s2=="INTERACTOR_B");

  if (s3=="CONF_VAL")
    {
      haveConfVal = true;
    }
  else
    {
      id1 = s3;
    }

  bool firstTime = true;
  std::string x;
  while (ifs >> x)
    {
      if (firstTime && !haveConfVal)
	{
	  id2 = x;//ifs >> id2;
	}
      else
	{
	  id1 = x;//ifs >> id1;
	  ifs >> id2;
	  if (haveConfVal)
	    {
	      ifs >> confVal;
	      assert( 0<= confVal && confVal <= 1);
	    }	  
	}

      It id1_it= id2sp.find(id1);
      if (id1_it == id2sp.end())
	{
	  idList.push_back(id1);
	  id2sp[id1] = sp;
	}
      else
	{
	  assert( id1_it->second == sp);
	}

      It id2_it= id2sp.find(id2);
      if (id2_it == id2sp.end())
	{
	  idList.push_back(id2);
	  id2sp[id2] = sp;
	}
      else
	{
	  assert( id2_it->second == sp);
	}
      
      StrPair idpair = (id1 < id2) ? StrPair(id1,id2) : StrPair(id2, id1);
      if (haveConfVal)
	{
	  idpair2IntxnScore[idpair] = confVal;
	}
      else
	{
	  idpair2IntxnScore[idpair] = 1;
	}
      numVals ++;

      firstTime = false;
    } 
  return numVals;
}






int  readBlastScoresFile(const std::string & fname, 
			 const std::string & sp1,
			 const std::string & sp2,
			 VecStr & idList, 
			 Str2Str & id2sp, 
			 StrPair2Dbl & idpair2BlastScore)
{
  
  typedef Str2Str::const_iterator It;
  
  int numVals = 0;
  std::cerr << "blast file: " << fname << std::endl;
  std::ifstream ifs(fname.c_str());
  
  char line[200];
  std::string id1, id2;
  double bitScore = 0;
  double eVal = 0;

  while (ifs.getline(line,200-1))
    {      
      std::istringstream istrstr(line);
      //istrstr >> id1 >> id2 >> bitScore >> eVal;
      istrstr >> id1 >> id2 >> bitScore;

      /*
	ifs >> id1;
      ifs >> id2;
      ifs >> bitScore;
      if (ifs.eof())
	break;
      ifs >> eVal;
      */
      //std::cerr << id1 << " "<< id2 << " " << bitScore << " eVal: " << eVal << std::endl;
      assert( 0 <= eVal && eVal < 1); 

      It id1_it= id2sp.find(id1);
      if (id1_it == id2sp.end())
	{
	  idList.push_back(id1);
	  id2sp[id1] = sp1;
	}
      else
	{
	  assert( id1_it->second == sp1);
	}

      It id2_it= id2sp.find(id2);
      if (id2_it == id2sp.end())
	{
	  idList.push_back(id2);
	  id2sp[id2] = sp2;
	}
      else
	{
	  assert( id2_it->second == sp2);
	}

      StrPair idpair = (id1 < id2) ? StrPair(id1,id2) : StrPair(id2, id1);
      idpair2BlastScore[idpair] = bitScore; //eVal;
      numVals ++;
    }
  //std::cerr << "Blast score file: "<< fname << " with # entries: " << numVals
  //<< std::endl;
  return numVals;
}






void removeNonParalogWithinSpeciesScores(StrPair2Dbl & idpair2BlastScore,
					 StrPair2Int & inParalogSet,
					 const VecStr & idList,
					 const VecStr & spList,
					 const Str2Str & id2sp)
{
  typedef StrPair2Dbl::iterator It1;
  typedef StrPair2StrPair2Int::const_iterator It2;
  typedef StrPair2Int::const_iterator It3;
  
  
  Str2Dbl I0;
  for (int i=0; i < idList.size(); i++)
    {
      I0[idList.at(i)] = 0;
    }

  Str2Dbl id2highestBlastScore = I0;
  inParalogSet.clear();
  
  //for each id, compute the highest blast score it has with somebody
  //  not in its species
  for (It1 it1=idpair2BlastScore.begin();
	   it1 != idpair2BlastScore.end();
	   it1++)
    {
      const std::string & id1 = it1->first.first;
      const std::string & id2 = it1->first.second;
      if (id1 != id2 && id2sp.find(id1)->second != id2sp.find(id2)->second)
	{
	  if (id2highestBlastScore[id1] < it1->second)
	    id2highestBlastScore[id1] = it1->second;
	  if (id2highestBlastScore[id2] < it1->second)
	    id2highestBlastScore[id2] = it1->second;
	}
    }

  
  //remove all within-species entries for which the blast-score
  // is less than the highest outside-species blast score for either id
  // add this entry to IN-PARALOG list
  int count = 0;
  for (It1 it1=idpair2BlastScore.begin();
       it1 != idpair2BlastScore.end();
	   )
    {
      const std::string & id1 = it1->first.first;
      const std::string & id2 = it1->first.second;
      if (id2sp.find(id1)->second == id2sp.find(id2)->second)
	{
	  if ((id2highestBlastScore[id1] > it1->second ||
	       id2highestBlastScore[id2] > it1->second) &&
	      (id1 != id2))
	    {
	      It1 delIt = it1;
	      it1++;
	      idpair2BlastScore.erase(delIt);
	      count ++;
	    }
	  else
	    {
	      inParalogSet[it1->first] = 1;
	      it1++;
	    }
	}
      else
	it1++;
    }
  std::cerr << "Removed " << count 
	    << " within-species entries" << std::endl;
}

void  keepOnlySomeSpeciesPairScalingFactors(StrPair2Dbl & Wij, 
					    double W)
{
  typedef StrPair2Dbl::iterator It1;

  StrPair2Int validSpPairs;

  validSpPairs[StrPair("dmela","hsapi")] = 1;
  validSpPairs[StrPair("hsapi","dmela")] = 1;
  validSpPairs[StrPair("dmela","mmusc")] = 1;
  validSpPairs[StrPair("mmusc","dmela")] = 1;
  validSpPairs[StrPair("scere","hsapi")] = 1;
  validSpPairs[StrPair("hsapi","scere")] = 1;
  validSpPairs[StrPair("scere","mmusc")] = 1;
  validSpPairs[StrPair("mmusc","scere")] = 1;
  validSpPairs[StrPair("scere","celeg")] = 1;  
  validSpPairs[StrPair("celeg","scere")] = 1;
  
  for (It1 it1=Wij.begin(); it1 != Wij.end(); it1++)
    {
      if (validSpPairs.find(it1->first) == validSpPairs.end())
	{
	  it1->second = W;
	}
    }
}



void rescaleScoresBySpDistance(StrPair2Dbl & idpair2BlastScore,
			       const StrPair2Int & inParalogSet,
			       const VecStr & idList,
			       const VecStr & spList,
			       const Str2Str & id2sp,
			       bool doWithinSpecies)
{
  typedef StrPair2Dbl::iterator It1;
  typedef StrPair2StrPair2Int::const_iterator It2;
  typedef StrPair2Int::const_iterator It3;
  
  
  Str2Dbl I0;
  for (int i=0; i < idList.size(); i++)
    {
      I0[idList.at(i)] = 0;
    }
  
  //for each pair of species, 
  // for each id in either of the species, find the highest blast score it 
  // is involved in. add this entry to ORTHOLOG list
  StrPair2StrPair2Int sppair2orthlogSet;
  StrPair2Int S0;

  for (int i =0; i < spList.size(); i++)
    for (int j=i+1; j < spList.size(); j++)
      {
	std::string sp1 = spList.at(i);
	std::string sp2 = spList.at(j);
	StrPair spair(sp1,sp2);
	sppair2orthlogSet[spair] = S0;

	//first, find the highest scores
	Str2Dbl id2thisSpPairHighScore = I0;
	for (It1 it1 = idpair2BlastScore.begin();
	     it1 != idpair2BlastScore.end();
	     it1++)
	  {
	    const std::string & id1 = it1->first.first;
	    const std::string & id2 = it1->first.second;
	    const std::string & spA = id2sp.find(id1)->second;
	    const std::string & spB = id2sp.find(id2)->second;
	    if ((spA==sp1 && spB==sp2) || (spB==sp1 && spA==sp2))
	      {
		if (id2thisSpPairHighScore[id1] < it1->second)
		  {
		    id2thisSpPairHighScore[id1] = it1->second;
		  }
		if (id2thisSpPairHighScore[id2] < it1->second)
		  {
		    id2thisSpPairHighScore[id2] = it1->second;
		  }
	      }
	  }

	//now, get the ortholog  list
	for (It1 it1 = idpair2BlastScore.begin();
	     it1 != idpair2BlastScore.end();
	     it1++)
	  {
	    const std::string & id1 = it1->first.first;
	    const std::string & id2 = it1->first.second;
	    const std::string & spA = id2sp.find(id1)->second;
	    const std::string & spB = id2sp.find(id2)->second;
	    if ((spA==sp1 && spB==sp2) || (spB==sp1 && spA==sp2))
	      {
		if (fabs(id2thisSpPairHighScore[id1]-(it1->second)) <EPS)
		  {
		    //DEBUGP(9.1);
		    sppair2orthlogSet[spair][it1->first] = 1;
		  }
		if (fabs(id2thisSpPairHighScore[id2]-(it1->second)) <EPS)
		  {
		    //DEBUGP(9.2);
		    sppair2orthlogSet[spair][it1->first] = 1;
		  }		  
	      }    
	  }
      }
  

  //compute W, add the ortholog set
  double W =0;
  int nW = 0;
  for (It2 it2 = sppair2orthlogSet.begin();
       it2 != sppair2orthlogSet.end();
       it2++)
    {
      for (It3 it3 = it2->second.begin();
	   it3 != it2->second.end();
	   it3++)
	{
	  W += idpair2BlastScore[it3->first];
	  nW++;
	}
    }


  if (doWithinSpecies)
    {
      // add the in-paralog set to W
      for (It3 it3 = inParalogSet.begin();
	   it3 != inParalogSet.end();
	   it3++)
	{
	  W += idpair2BlastScore[it3->first];
	  nW++;
	}
    }

  // take the average
  if (nW > 0)
    W = W/nW;
  

  //for each pair of species, compute Wij
  StrPair2Dbl Wij;
  for (int i=0; i < spList.size(); i++)
    for (int j=i+1; j < spList.size(); j++)
      {
	std::string sp1 = spList.at(i);
	std::string sp2 = spList.at(j);
	StrPair spair(sp1,sp2);
	double w = 0;
	int nw = 0;
	Wij[spair] = 0;
	
	//get the scores for orthologs corresponding to this pair
	for (It3 it3 = sppair2orthlogSet[spair].begin();
	     it3 != sppair2orthlogSet[spair].end();
	     it3++)
	  {
	    //DEBUGP(9.5);
	    w += idpair2BlastScore[it3->first];
	    nw++;
	  }
	std::cerr << spair << std::endl;
	if (nw >0)
	  Wij[spair] = w/nw;
      }
  
  //copy entries so that (sp1,sp2) and (sp2,sp1) both exist in Wij
  for (int i=0; i < spList.size(); i++)
    for (int j=0; j < i; j++)
      {
	StrPair spair1(spList.at(i), spList.at(j));
	StrPair spair2(spList.at(j), spList.at(i));
	Wij[spair1] = Wij[spair2];
      }


  //compute Wii for in-paralogs
  for (int i=0; i < spList.size(); i++)
    {
      std::string sp = spList.at(i);
      StrPair spair(sp,sp);

      if (!doWithinSpecies)
	{
	  Wij[spair] = W;
	}

      else
	{
	  double w =0;
	  int nw=0;
	  Wij[spair] = 0;
	  for (It3 it3= inParalogSet.begin();
	       it3 != inParalogSet.end();
	       it3++)
	    {
	      const std::string & id1 = it3->first.first;
	      const std::string & id2 = it3->first.second;
	      const std::string & spA = id2sp.find(id1)->second;
	      const std::string & spB = id2sp.find(id2)->second;
	      assert( spA == spB);
	      if (spA == sp)
		{
		  w += idpair2BlastScore[it3->first];
		  nw++;
		}
	    }
	  if (nw > 0)
	    Wij[spair] = w/nw;
	}
    }
  

  /* Sep 3, 2006: coded this to see if it makes a difference when we only
                  rescale scores between some species pairs.
		  Didn't seem to work, hence commenting it out
  */
  //keepOnlySomeSpeciesPairScalingFactors(Wij, W);


  //print out the normalization weights
  for (int i=0; i < spList.size(); i++)
    for (int j=i; j < spList.size(); j++)
      std::cerr << "Scaling factor for " << spList.at(i) << "-" 
		<< spList.at(j) <<" : " 
		<< Wij[StrPair(spList.at(i),spList.at(j))]/W << std::endl;

  //for each score, divide by Wij/W
  for (It1 it1 = idpair2BlastScore.begin();
       it1 != idpair2BlastScore.end();
       it1++)
    {
      const std::string & id1 = it1->first.first;
      const std::string & id2 = it1->first.second;
      const std::string & spA = id2sp.find(id1)->second;
      const std::string & spB = id2sp.find(id2)->second;      
      StrPair spair(spA,spB);
      it1->second = it1->second / ( Wij[spair]/ W);
    }
}  





void removeLowEntries(StrPair2Dbl & idpair2BlastScore,
		      const VecStr & idList,
		      const VecStr & spList,
		      const Str2Str & id2sp,
		      double lowThreshold = 0.5)
{
  typedef StrPair2Dbl::iterator It1;
  typedef StrPair2StrPair2Int::const_iterator It2;
  typedef StrPair2Int::const_iterator It3;
  
  std::cerr << "in removeLowEntries, threshold: " 
	    << lowThreshold << std::endl;
  
  Str2Dbl I0;
  for (int i=0; i < idList.size(); i++)
    {
      I0[idList.at(i)] = 0;
    }

  int count =0;
  for (int i =0; i < spList.size(); i++)
    for (int j=i+1; j < spList.size(); j++)
      {
	std::string sp1 = spList.at(i);
	std::string sp2 = spList.at(j);
	StrPair spair(sp1,sp2);

	DEBUGP(901);

	//first, find the highest scores
	Str2Dbl id2thisSpPairHighScore = I0;
	for (It1 it1 = idpair2BlastScore.begin();
	     it1 != idpair2BlastScore.end();
	     it1++)
	  {
	    const std::string & id1 = it1->first.first;
	    const std::string & id2 = it1->first.second;
	    const std::string & spA = id2sp.find(id1)->second;
	    const std::string & spB = id2sp.find(id2)->second;
	    if ((spA==sp1 && spB==sp2) || (spB==sp1 && spA==sp2))
	      {
		if (id2thisSpPairHighScore[id1] < it1->second)
		  {
		    id2thisSpPairHighScore[id1] = it1->second;
		  }
		if (id2thisSpPairHighScore[id2] < it1->second)
		  {
		    id2thisSpPairHighScore[id2] = it1->second;
		  }
	      }
	  }

	DEBUGP(902);	

	//now, get the ortholog  list
	for (It1 it1 = idpair2BlastScore.begin();
	     it1 != idpair2BlastScore.end();
	     )
	  {
	    const std::string & id1 = it1->first.first;
	    const std::string & id2 = it1->first.second;
	    const std::string & spA = id2sp.find(id1)->second;
	    const std::string & spB = id2sp.find(id2)->second;
	    if ((spA==sp1 && spB==sp2) || (spB==sp1 && spA==sp2))
	      {
		if (id2thisSpPairHighScore[id1]*lowThreshold > it1->second ||
		    id2thisSpPairHighScore[id2]*lowThreshold > it1->second)
		  {
		    It1 delIt = it1;
		    it1++;
		    idpair2BlastScore.erase(delIt);
		    count ++;
		  }
		else
		  {
		    it1++;
		  }
	      }
	    else
	      {
		it1++;
	      }
	  }

	DEBUGP(903);
      }

  std::cerr << "Removed " << count 
	    << " low entries" << std::endl;
}





//normalizeBlastScores follows the approach of Li, Stoecker and Roos
// from Genome Research 13:2178-2189 (2003), (see Fig 2 in that paper)
//
// doWithinSpecies is useful if you don't want to do within-species
//  normalization (e.g. with network-based scores)
//
//  removeLowEntries is not a LSR thing. I added it to get rid of 
//    spurious BLAST edges
void normalizeBlastLikeScores(StrPair2Dbl & idpair2BlastScore,
			      const VecStr & idList,
			      const VecStr & spList,
			      const Str2Str & id2sp,
			      bool doWithinSpecies,
			      bool doRescaling,
			      bool doRemoveLowEntries,
			      double blastLowThreshold)
{
  typedef StrPair2Dbl::iterator It1;
  typedef StrPair2StrPair2Int::const_iterator It2;
  typedef StrPair2Int::const_iterator It3;
  
  Str2Dbl I0;
  for (int i=0; i < idList.size(); i++)
    {
      I0[idList.at(i)] = 0;
    }

  StrPair2Int inParalogSet;  

  if (doWithinSpecies)
    {
      removeNonParalogWithinSpeciesScores(idpair2BlastScore,
					  inParalogSet,
					  idList, spList, id2sp);
    }

  if (doRescaling)
    {
      rescaleScoresBySpDistance(idpair2BlastScore, inParalogSet,
				idList, spList, id2sp,
				doWithinSpecies);
    }  
  
  if (doRemoveLowEntries)
    {
      removeLowEntries(idpair2BlastScore,
		       idList, spList, id2sp, blastLowThreshold);
    }

}


void grabSequenceData(const std::string & dataDescFile, 
		      VecStr & spList, 
		      VecStr & idList, 
		      Str2Str & id2sp,
		      StrPair2Dbl & idpair2BlastScore)
{
  /* file format:
     line 1: <directory path>
     line 2: <suffix to append>
     line 3: K (=number of species)
     line 4: species_1
     line 3+K: species_K
     
     the blast scores are in <dir>/sp1_sp2<suffix>.evals
  */
  std::ifstream ifs(dataDescFile.c_str()); 
  std::string dirPath, suffix;
  int K = 0;
  
  assert(ifs);
  ifs >> dirPath >> suffix >> K;
  if (suffix == "-")
    suffix = "";
  assert(ifs);
    
  spList.clear();
  idList.clear();
  idpair2BlastScore.clear();
  id2sp.clear();
  for (int i=0; i < K; i++)
    {
      std::string s;
      ifs >> s;
      spList.push_back(s);
    }

  std::string fname;
  std::sort(spList.begin(), spList.end());

  for (int i=0; i < K; i++)
    {
      for (int j=i; j < K; j++)
	{
	  fname = dirPath + '/' + spList[i] + "-" + 
	          spList[j] + suffix + ".evals";
	  readBlastScoresFile(fname, 
			      spList[i], spList[j],
			      idList, id2sp, 
			      idpair2BlastScore);
	}
    }
  
  std::cerr << "Number of blast scores read: " << idpair2BlastScore.size()
	    << std::endl;
  
}




void grabSequenceAndInteractionData(const std::string & dataDescFile, 
				    VecStr & spList, 
				    VecStr & idList, 
				    Str2Str & id2sp,
				    StrPair2Dbl & idpair2BlastScore, 
				    StrPair2Dbl & idpair2IntxnScore)
{
  /* file format:
     line 1: <directory path>
     line 2: <suffix to append>
     line 3: K (=number of species)
     line 4: species_1
     line 3+K: species_K
     
     the blast scores are in <dir>/sp1_sp2<suffix>.evals
     the intxns are in <dir>/sp1<suffix>.tab     
  */
  std::ifstream ifs(dataDescFile.c_str());
  
  std::string dirPath, suffix;
  int K = 0;
  
  assert(ifs);
  ifs >> dirPath >> suffix >> K;
  if (suffix == "-")
    suffix = "";
  assert(ifs);
    
  spList.clear();
  idList.clear();
  idpair2IntxnScore.clear();
  idpair2BlastScore.clear();
  id2sp.clear();
  for (int i=0; i < K; i++)
    {
      std::string s;
      ifs >> s;
      spList.push_back(s);
    }

  std::string fname;
  std::sort(spList.begin(), spList.end());

  for (int i=0; i < K; i++)
    {
      for (int j=i; j < K; j++)
	{
	  fname = dirPath + '/' + spList[i] + "-" + 
	          spList[j] + suffix + ".evals";
	  readBlastScoresFile(fname, 
			      spList[i], spList[j],
			      idList, id2sp, 
			      idpair2BlastScore);
	}
    }
  
  std::cerr << "Number of blast scores read: " << idpair2BlastScore.size()
	    << std::endl;


  // //normalize the Blast Scores as per Li, Stoeckert, Roos (Genome Res 2003)
  // normalizeBlastLikeScores(idpair2BlastScore, idList, spList, id2sp);
  //std::cerr << "Number of blast scores after normalization: " 
  //	    << idpair2BlastScore.size() << std::endl;


  for (int i=0; i < K; i++)
    {
      fname = dirPath + '/' + spList[i] + suffix + ".tab";
      readIntxnFile(fname, spList[i], idList, id2sp, idpair2IntxnScore);
    }
}






void grabInteractionData(const std::string & dataDescFile, 
			 VecStr & spList, 
			 VecStr & idList, 
			 Str2Str & id2sp,
			 StrPair2Dbl & idpair2IntxnScore)
{
  /* file format:
     line 1: <directory path>
     line 2: <suffix to append>
     line 3: K (=number of species)
     line 4: species_1
     line 3+K: species_K
     
     the blast scores are in <dir>/sp1_sp2<suffix>.evals
     the intxns are in <dir>/sp1<suffix>.tab     
  */
  std::ifstream ifs(dataDescFile.c_str());
  
  std::string dirPath, suffix;
  int K = 0;
  
  assert(ifs);
  ifs >> dirPath >> suffix >> K;
  if (suffix == "-")
    suffix = "";
  assert(ifs);
    
  spList.clear();
  idList.clear();
  idpair2IntxnScore.clear();
  id2sp.clear();
  for (int i=0; i < K; i++)
    {
      std::string s;
      ifs >> s;
      spList.push_back(s);
    }

  std::string fname;
  std::sort(spList.begin(), spList.end());


  for (int i=0; i < K; i++)
    {
      fname = dirPath + '/' + spList[i] + suffix + ".tab";
      readIntxnFile(fname, spList[i], idList, id2sp, idpair2IntxnScore);
    }
}






void grabSequenceAndInteractionData(const std::string & dataDescFile, 
				    VecStr & spList, 
				    VecStr & idList, 
				    Str2Str & id2sp,
				    StrPair2Dbl & idpair2BlastScore, 
				    StrPair2Dbl & idpair2IntxnScore,
				    StrPair2Dbl & idpair2wbmScore)
{
   /* file format:
     line 1: <directory path>
     line 2: <suffix to append>
     line 3: K (=number of species)
     line 4: species_1
     line 3+K: species_K
     
     the blast scores are in <dir>/sp1_sp2<suffix>.evals
     the intxns are in <dir>/sp1<suffix>.tab     
  */
  std::ifstream ifs(dataDescFile.c_str());
  
  std::string dirPath, suffix;
  assert(ifs);
  ifs >> dirPath >> suffix;

  grabSequenceAndInteractionData(dataDescFile, 
				 spList, 
				 idList, 
				 id2sp,
				 idpair2BlastScore, 
				 idpair2IntxnScore);

  std::string dummy;
  ifs >> dummy;
  for (int i=0; i < spList.size(); i++)
    {
      ifs >> dummy;
    }
  std::string wbmfname;
  ifs >> wbmfname;
  if (wbmfname.size()==0)
    {
#ifdef USE_OLD_WBMSCORES
      wbmfname = dirPath + "/" + "wbm_l1_fromSeq.dat";
#else
      wbmfname = dirPath + "/" + "wbm_scores_9dec.dat";
#endif
    }
  else
    {
      wbmfname = dirPath + "/" + wbmfname;
    }
  std::cerr << "WBM file name: " << wbmfname << std::endl;

  int numWbmScores = readWBMscoreFile(wbmfname, idList, idpair2wbmScore);
  std::cerr << "Number of wbm scores read: " << numWbmScores << std::endl;

  return;
}







double getSpeciesDistance(const std::string & sp1, const std::string &sp2)
{

#ifndef USE_SPECIES_DISTANCE

  //these species are so far apart that the species distance might well be 1
  return 1.0;

#else

  static bool haveDistMatrix = false;
  static Str2Str2Dbl distMatrix;
  if (!haveDistMatrix)
    {
      Str2Dbl dmela, hpylo, scere, celeg, ecoli;

      dmela["hpylo"] = 1.340; dmela["scere"] = 0.411;  dmela["celeg"] = 0.139; 
      dmela["ecoli"] = 1.738; 

      hpylo["dmela"] = 1.340; hpylo["scere"] = 1.254;  hpylo["celeg"] = 1.384; 
      hpylo["ecoli"] = 0.510; 

      scere["dmela"] = 0.411; scere["hpylo"] = 1.254;  scere["celeg"] = 0.414; 
      scere["ecoli"] = 1.579; 

      celeg["dmela"] = 0.139; celeg["hpylo"] = 1.384;  celeg["scere"] = 0.414; 
      celeg["ecoli"] = 1.740; 

      ecoli["hpylo"] = 0.510; ecoli["scere"] = 1.579;  ecoli["celeg"] = 1.740; 
      ecoli["dmela"] = 1.738; 

      distMatrix["scere"] = scere;
      distMatrix["dmela"] = dmela;
      distMatrix["hpylo"] = hpylo;
      distMatrix["celeg"] = celeg;
      distMatrix["ecoli"] = ecoli;     
      haveDistMatrix = true;
    }

  return ((distMatrix.find(sp1)->second).find(sp2)->second);
#endif
}





int countNumClusters(const Str2Int & id2clstr)
{
  Int2Int clstrSet;

  typedef Str2Int::const_iterator It;
  for (It it=id2clstr.begin(); it != id2clstr.end(); it++)
    {
      clstrSet[it->second] = 1;
    }
  return clstrSet.size();
}




int countNumNontrivialClusters(const Str2Int & id2clstr)
{
  //compute cluster sizes
  Int2Int clstrSizes;
  
  typedef Str2Int::const_iterator It;
  typedef Int2Int::const_iterator It2;

  for (It it1 = id2clstr.begin(); it1 != id2clstr.end(); it1++)
    {
      if (clstrSizes.find(it1->second) == clstrSizes.end())
	clstrSizes[it1->second] = 1;
      else
	clstrSizes[it1->second] += 1;
    }
  

  int numNonTrivialClusters = 0;
  for (It2 it2 = clstrSizes.begin(); it2 != clstrSizes.end(); it2++)
    {
      if (it2->second >1) 
	numNonTrivialClusters++;
    }
  return numNonTrivialClusters;
}





//returns fraction of edges within clusters
double countWithinClusterEdges(const Str2Int & id2clstr,
			       const StrPair2Dbl & idpair2IntxnScore,
			       Int2Int & clstr2withinEdgeCount)
{
  typedef StrPair2Dbl::const_iterator It1;
  typedef Int2Int::iterator It2;
  
  clstr2withinEdgeCount.clear();

  int totalNumEdges = 0;
  for (It1 it1=idpair2IntxnScore.begin();
       it1 != idpair2IntxnScore.end();
       it1++)
    {
      const std::string & id1 = it1->first.first;
      const std::string & id2 = it1->first.second;
      if (id2clstr.find(id1) != id2clstr.end() &&
	  id2clstr.find(id2) != id2clstr.end())
	totalNumEdges ++;
    }


  double withinClstrEdges = 0;
  for (It1 it1=idpair2IntxnScore.begin();
       it1 != idpair2IntxnScore.end();
       it1++)
    {
      const std::string & id1 = it1->first.first;
      const std::string & id2 = it1->first.second;
      if (id2clstr.find(id1) != id2clstr.end() &&
	  id2clstr.find(id2) != id2clstr.end() &&
	  id2clstr.find(id1)->second == id2clstr.find(id2)->second)
	{
	  withinClstrEdges++;
	  int clid = id2clstr.find(id1)->second ;
	  if (clstr2withinEdgeCount.find(clid)==clstr2withinEdgeCount.end())
	    clstr2withinEdgeCount[clid]=0;
	  
	  clstr2withinEdgeCount[clid]++;
	}
    }
  return withinClstrEdges/totalNumEdges;
}



void grabGO_DIPdata(Str2Int & idSetInGO,
		    StrPair2Dbl & gopairSet)
{
  idSetInGO.clear();
  gopairSet.clear();

#ifndef BIOCLUSTER
  std::string fname="/net/snunit/scratch/rsingh/work/phylo-networks/data/go-seq/DIP_GO_assoc.tab";
#else
  std::string fname="/r100/berger/rsingh/phylo-networks/data/go-seq/DIP_GO_assoc.tab";
#endif
  std::ifstream ifs(fname.c_str());
  
  char line[1000];
  std::string id;
  while (ifs.getline(line,1000-1))
    {
      //replace ";" by " "
      for (int i=0; line[i]!=0; i++)
	{
	  if (line[i]==';')
	    line[i] = ' ';
	}
      
      std::istringstream istrstr(line);
      istrstr >> id;
      VecStr goTerms;
      std::copy( std::istream_iterator<std::string>(istrstr),
		 std::istream_iterator<std::string>(),
		 std::back_inserter(goTerms));
      
      idSetInGO[id] = 1;
      for (int i=0; i < goTerms.size(); i++)
	{
	  for (int j=i+1; j < goTerms.size(); j++)
	    {
	      if (goTerms[i] == goTerms[j])
		continue;
	      StrPair spair;
	      if (goTerms[i] < goTerms[j])
		{
		  spair = StrPair(goTerms[i],goTerms[j]);
		}
	      else
		{
		  spair = StrPair(goTerms[j],goTerms[i]);
		}
	      if (gopairSet.find(spair) == gopairSet.end())
		{
		  gopairSet[spair] =1;
		}
	      else
		{
		  gopairSet[spair] += 1;
		}
	    }
	}
    }
}




void subsetDataBySetOfIds(StrPair2Dbl & idpair2value,
			  const Str2Int smallIdSet)
{
  typedef StrPair2Dbl::iterator It1;
  for (It1 it1=idpair2value.begin();
       it1 != idpair2value.end();
       )
    {
      const std::string & id1 = it1->first.first;
      const std::string & id2 = it1->first.second;
      if (smallIdSet.find(id1) == smallIdSet.end() ||
	  smallIdSet.find(id2) == smallIdSet.end())
	{
	  It1 delIt = it1;
	  it1++;
	  idpair2value.erase(delIt);
	}
      else
	{
	  it1++;
	}
    }
}



VecStr getIdListFromPairs(const StrPair2Dbl & idpairSet)
{
  typedef StrPair2Dbl::const_iterator It1;
  typedef Str2Int::const_iterator It2;
  
  Str2Int idSet;
  for (It1 it1=idpairSet.begin(); it1 != idpairSet.end(); it1++)
    {
      idSet[it1->first.first] = 1;
      idSet[it1->first.second] = 1;
    }

  VecStr retVal;
  for (It2 it2 = idSet.begin(); it2 != idSet.end(); it2++)
    {
      retVal.push_back(it2->first);
    }
  return retVal;
}




void grabClusteringData(Str2Int & id2clstr,
			const std::string & clstrFile,
			const VecStr & idList,
			const Str2Str & id2sp)
{
  typedef Str2Str::const_iterator It1;
  
  id2clstr.clear();

  std::ifstream ifs(clstrFile.c_str());

  char line[200];
  std::string id1;
  int clid;

  while (ifs.getline(line,200-1))
    {      
      std::istringstream istrstr(line);
      istrstr >> id1 >> clid;
      if ( id2sp.find(id1) != id2sp.end())
	id2clstr[id1] = clid;
    }
}




int firstNonVisited(const VecInt & haveVisited)
{
  int N = haveVisited.size();
  int i=0;
  while (i<N && haveVisited.at(i))
    i++;
    
  return i;
}


int getConnectedComponents(const VecVecIntDblPair & scoreAdjList,
			   ListListInt & connectedSets)
{
  typedef VecIntDblPair::const_iterator It1;

  int N = scoreAdjList.size();

  VecInt haveVisited(N, 0);

  connectedSets.clear();

  int currIdx = 0;

  for (int i = firstNonVisited(haveVisited); 
       i < N; 
       i = firstNonVisited(haveVisited))
    {
      VecInt seenThisTime(N,0);
      ListInt queue;
      ListInt thisSet;
      int j=-1;

      queue.push_back(i);
      thisSet.push_back(i);
      haveVisited[i] = 1;
      
      while (!queue.empty())
	{
	  j = *(queue.begin());
	  queue.pop_front();

	  for (It1 it1=scoreAdjList.at(j).begin();
	       it1 != scoreAdjList.at(j).end();
	       it1++)
	    {
	      if (!haveVisited[it1->first])
		{
		  queue.push_back(it1->first);
		  thisSet.push_back(it1->first);
		  haveVisited[it1->first] = 1;
		}
	    }
	}
      connectedSets.push_back(thisSet);
    }
  return connectedSets.size();
}



void normalizeToZeroOneScores(VecVecIntDblPair & scoreAdjList)
{
  VecDbl maxScores(scoreAdjList.size(),1);
  for (int i=0; i < scoreAdjList.size(); i++)
    for (int jidx=0; jidx < scoreAdjList.at(i).size(); jidx++)
      {
	int j= scoreAdjList.at(i).at(jidx).first;
	double v = scoreAdjList.at(i).at(jidx).second;
	if (i==j)
	  {
	    if (v > maxScores[i])
	      maxScores[i] = v;
	  }
	if (v > maxScores[i])
	  maxScores[i] = v;
	if (v > maxScores[j])
	  maxScores[j] = v;
      }

  for (int i=0; i < scoreAdjList.size(); i++)
    { 
      bool seenDiagElem = false;
    for (int jidx=0; jidx < scoreAdjList.at(i).size(); jidx++)
      {
	int j= scoreAdjList.at(i).at(jidx).first;
	double v = scoreAdjList.at(i).at(jidx).second;
	if (i==j)
	  {
	    seenDiagElem = true;
	    scoreAdjList.at(i).at(jidx).second = 1;
	  }
	else
	  {
	    scoreAdjList.at(i).at(jidx).second = v / sqrt( maxScores[i]*
							   maxScores[j]);
	  }
      }
    if (!seenDiagElem)
      {
	scoreAdjList[i].push_back(IntDblPair(i,1.0));
      }

    for (int jidx=0; jidx < scoreAdjList.at(i).size(); jidx++)
      {
	assert (scoreAdjList.at(i).at(jidx).second <= 1);
      }
    }
}


void writeAllScoresToFile(const char * outfile,
			  const VecVecIntDblPair & scoreAdjList,
			  const VecStr & idList)
{
  std::ofstream ofs(outfile); 
  
  for (int i=0; i < scoreAdjList.size(); i++)
    {
      for (int jidx =0; jidx < scoreAdjList.at(i).size(); jidx++)
	{
	  int j = scoreAdjList.at(i).at(jidx).first;
	  double v = scoreAdjList.at(i).at(jidx).second;
	  if (i <= j) 
	    {
	      ofs << idList.at(i) << " " << idList.at(j) << " " << v
		  << std::endl;
	    }
	}
    }
  ofs.close();
}
#endif
