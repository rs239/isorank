#ifndef INCLUDED_defs
#include "defs.h"
#endif

#ifndef INCLUDED_utils
#include "Utils.h"
#endif


#ifndef INCLUDED_k_partite
#include "k_partite.h"
#endif


#include <math.h>



/* TO IMPLEMENT:
 
*/

void usage(const std::string & errMsg="")
{
  std::string usageStr="Usage: multiway_kpartite [--prefix <file-prefix>] [--scorefile <match-score-file>] [--K <max_iter>] [--thresh <in_percent>] [--alpha <weight_of_ntwk_term>] [--maxveclen <max_len_of_eigenvector>] [--cl_min_pri_thresh <thresh1>] [--cl_min_sec_thresh <thresh2>] [--cl_max_per_sp <max1>] <data-description.txt> \n"
    "\n"
    "--prefix : prefix of the output .txt and .dat files (Default = '')\n\n"
    "--scorefile <match-score-file> : the file containing pairwise scores\n"
    " after doing network alignment. Forms input to k-partite clustering\n\n"
    "--K <maxiter>: max number of iterations of the power method (>3, <100)\n"
    "--thresh <pct-thresh>: the algorithm stops when the average percentage \n"
    "\tchange in value of Rij is less than <pct-thresh>\n"
    "--cl_min_pri_thresh <thresh1> : in the k_partite step, the fractional \n"
    "\tthreshold to use during constructing the primary cluster\n"
    "--cl_min_sec_thresh <thresh2> : in the k_partite step, the fractional \n"
    "\tthreshold to use when adding the secondary nodes\n"
    "--cl_max_per_sp <max1> : in the k_partite step, the max number \n"
    "\t of species per cluster\n"
    "The program stops when either of the above limits are violated\n\n"

    "--alpha: weight of network data (in [0,1])\n"

    "--maxveclen: for really large graphs, limits the number of non-zero \n"
    "\tentries in R\n\n"
    "The argument (data-description.txt) should specify the list of species\n"
    "\tto be used as well as the path to the various files. Look at src/data.in\n"
    "X.tab should have a first line as 'INTERACTOR_A INTERACTOR_B'\n"
    "\tEach following line should be of type 'id1 id2 wt', where 3rd column\n" 
    "\tis optional\n"
    "\n"
    "X-Y.evals should contain the non-network information. Each line is of \n"
    "\tthe form 'id_X id_Y score'. The scores can be BLAST Bit-scores or \n"
    "\tsomething else (even binary values!). This is also used to initialize\n"
    "\tthe power method. You can set alpha=1 and specify random values \n"
    "\there-- the file will then only be used to set the initial value \n"
    "\tof R in the power method\n"
    "\n"
    "X1, X2, X3... MUST NOT SHARE ANY NODE LABELS! Thus, if you are \n"
    "\tcomparing two (or more) conformations of a single protein, rename all nodes \n"
    "\tin one of the conformations (e.g., by prefixing/suffixing with \n"
    "\tsomething)\n";



  if (errMsg.size() > 0) 
    {
      std::cerr << "Error: " << errMsg << "\n\n";
    }
  std::cerr << usageStr;
}



void getoptions(int argc, char ** argv,
		int & K,
		double & convg_thresh,
		double & alpha,
		int & prunedVecSize,
		double & min_primary_fraction,
		double & min_secondary_fraction,
		int & max_per_species,
		std::string & strategy,
		std::string & filePrefix,
		std::string & dataDescFile,
		std::string & matchScoreFile)
{
  struct option longopts[] = {
    {"prefix",1,0,0},
    {"K",1,0,0},
    {"thresh",1,0,0},
    {"alpha",1,0,0},
    {"maxveclen",1,0,0},
    {"scorefile",1,0,0},
    {"cl_min_pri_thresh",1,0,0},
    {"cl_min_sec_thresh",1,0,0},
    {"cl_max_per_sp",1,0,0},
    {"help",0,0,0},
    {0,0,0,0}
  };
  
  int optindex=0;
  optind = 0;
  K = 20;
  convg_thresh = 1;
  alpha = 0.9;
  prunedVecSize = -1;
  strategy = "hsp";
  matchScoreFile = "";
  filePrefix = "tmp";
  min_primary_fraction = 0.1;
  min_secondary_fraction = 0.8;
  max_per_species = 5;
  while (1)
    {
      int c = getopt_long_only(argc, argv, "", longopts, &optindex);
      if (c == -1)
	break;
      if (optindex == 0) //"prefix"
	{
	  filePrefix = optarg;
	}
      /* "splist"
      if (optindex == 1) 
	{
	  std::string splist_str(optarg);
	  for (int i=0;i<splist_str.size();i++)
	    if (splist_str[i]==',')
	      splist_str[i] = ' ';

	  std::istringstream str_buf(splist_str.c_str());
	  splist.clear();
	  std::copy(std::istream_iterator<std::string>(str_buf),
		    std::istream_iterator<std::string>(),
		    std::back_inserter(splist));
	}
      */
      if (optindex == 1) //"K"
	{
	  K = atoi(optarg);
	  if (K<3 || K >= 100)
	    {
	      usage("K should be between 3 and 100B");
	      exit(1);
	    }
	}
      if (optindex == 2) //"thresh"
	{
	  convg_thresh = double(atof(optarg));
	  if (convg_thresh <= 0 || convg_thresh > 10)
	    {
	      usage("Need convg_thresh in (0,10)");
	      exit(1);
	    }	  
	}
      if (optindex == 3) //"alpha"
	{
	  alpha = double(atof(optarg));
	  if (alpha < 0 || alpha > 1)
	    {
	      usage("Bad 'alpha'");
	      exit(1);
	    }	 
	}
      if (optindex == 4) //"prunedVecSize"
	{
	  prunedVecSize = atoi(optarg);
	  if (prunedVecSize <= 0)
	    {
	      usage("Bad 'prunedvecsize' option");
	      exit(1);
	    }
	}
      if (optindex == 5) //"scorefile"
	{
	  matchScoreFile = optarg;
	}
      if (optindex == 6) //"cl_min_pri_thresh"
	{
	  min_primary_fraction = atof(optarg);
	  if (min_primary_fraction <= 0 || min_primary_fraction > 1)
	    {
	      usage("Bad 'cl_min_pri_thresh' option");
	      exit(1);
	    }
	}
      if (optindex == 7) //"cl_min_sec_thresh"
	{
	  min_secondary_fraction = atof(optarg);
	  if (min_secondary_fraction <= 0 || min_secondary_fraction > 1)
	    {
	      usage("Bad 'cl_min_sec_thresh' option");
	      exit(1);
	    }
	}
      if (optindex == 8) //"cl_max_per_sp"
	{
	  max_per_species = atoi(optarg);
	  if (max_per_species < 1 )
	    {
	      std::cerr << "blah: '" << optarg << "' " << max_per_species << std::endl;
	      usage("Bad 'cl_max_per_sp' option");
	      exit(1);
	    }
	}
      if (optindex == 9) //"help"
	{
	  usage("");
	  exit(1);
	}
    }

  assert(argc > optind);
  dataDescFile = argv[optind];
  return;
}




void countEdgesNodesBySp(const StrPair2Dbl & idpair2IntxnScore, 
			 const Str2Str & id2sp, 
			 const std::string & sp1, 
			 int & numNodes, 
			 int & numEdges)
{
  typedef StrPair2Dbl::const_iterator It1;
  typedef Str2Str::const_iterator It2;
  Str2Int seenId;
  numNodes =0;
  numEdges =0;
  for (It1 it1=idpair2IntxnScore.begin();
       it1 != idpair2IntxnScore.end();
       it1++)
    {
      std::string sp = id2sp.find(it1->first.first)->second;
      if (sp==sp1)
	{
	  numEdges++;

	  if (seenId.find(it1->first.first) == seenId.end())
	    {
	      seenId[it1->first.first] = 1;
	      numNodes++;
	    }

	  if (seenId.find(it1->first.second) == seenId.end())
	    {
	      seenId[it1->first.second] = 1;
	      numNodes++;
	    }
	}
    }
  
}



void getAdjListBySp_withAllGenes(const StrPair2Dbl & idpair2IntxnScore, 
		       const Str2Str & id2sp,
		       const Str2Int & id2idx,
		       const std::string & sp, 
		       VecInt & sp1_idxList, 
		       VecVecInt & sp1_adjList)
{
  typedef StrPair2Dbl::const_iterator It1;
  typedef Str2Int::const_iterator It2;
  typedef Str2Str::const_iterator It3;

  sp1_idxList.clear();
  sp1_adjList.clear();

  VecInt V0;

  Str2Int sp1_id2idx;
  for (It3 it3=id2sp.begin(); it3!=id2sp.end(); it3++)
    {
      if (it3->second == sp)
	{
	  sp1_idxList.push_back(id2idx.find(it3->first)->second);
	  sp1_adjList.push_back(V0);
	  sp1_id2idx[it3->first] = sp1_idxList.size()-1;
	}
    }


  int currIdx = 0;
  for (It1 it1=idpair2IntxnScore.begin(); 
       it1!=idpair2IntxnScore.end();
       it1++)
    {
      std::string sp1 = id2sp.find(it1->first.first)->second;
      if (sp1!=sp)
	continue;
      assert (sp1_id2idx.find(it1->first.first) != sp1_id2idx.end() &&
	      sp1_id2idx.find(it1->first.second) != sp1_id2idx.end());
      
      int a = sp1_id2idx.find(it1->first.first)->second;
      int b = sp1_id2idx.find(it1->first.second)->second;
      
      sp1_adjList[a].push_back(b);
      sp1_adjList[b].push_back(a);
      
    }
  
}


void  getAdjListBySp(const StrPair2Dbl & idpair2IntxnScore, 
		       const Str2Str & id2sp,
		       const Str2Int & id2idx,
		       const std::string & sp, 
		       int numNodes,
		       VecInt & sp1_idxList, 
		       VecVecInt & sp1_adjList)
{
  typedef StrPair2Dbl::const_iterator It1;
  typedef Str2Int::const_iterator It2;

  sp1_idxList.clear();
  sp1_adjList.clear();
  //DEBUGP(301);
  VecInt v0;
  sp1_adjList = VecVecInt(numNodes,v0);
  sp1_idxList = VecInt(numNodes,-1);

  Str2Int haveSeenId;
  int currIdx = 0;
  for (It1 it1=idpair2IntxnScore.begin(); 
       it1!=idpair2IntxnScore.end();
       it1++)
    {
      assert (haveSeenId.size() == currIdx);
      std::string sp1 = id2sp.find(it1->first.first)->second;
      if (sp1!=sp)
	continue;
      It2 it2a = haveSeenId.find(it1->first.first);
      int a,b;
      a=b=-1;
      if (it2a == haveSeenId.end())
	{
	  assert (haveSeenId.size() == currIdx);
	  haveSeenId[it1->first.first] = currIdx;
	  sp1_idxList[currIdx] = id2idx.find(it1->first.first)->second;
	  a = currIdx;
	  currIdx++;
	  if (currIdx != haveSeenId.size())
	    {
	      std::cerr << "error: " << currIdx 
			<< " " << haveSeenId.size() << std::endl;
	    }
	}
      else
	a = it2a->second;

      //DEBUGP(302);
      It2 it2b = haveSeenId.find(it1->first.second);
      if (it2b == haveSeenId.end())
	{
	  assert (haveSeenId.size() == currIdx);
	  haveSeenId[it1->first.second] = currIdx;
	  sp1_idxList[currIdx] = id2idx.find(it1->first.second)->second;
	  b = currIdx;
	  currIdx++;
	  /*if (currIdx != haveSeenId.size())
	    {
	      std::cerr << "error: " << currIdx 
			<< " " << haveSeenId.size() << std::endl;
	    }
	  */
	}
      else
	b = it2b->second;

      //DEBUGP(303);
      //std::cerr << currIdx << " " << numNodes << " " << haveSeenId.size() 
      //	<< " " << a << " " << b << std::endl;

      sp1_adjList[a].push_back(b);
      sp1_adjList[b].push_back(a);
      
    }
  //DEBUGP(304);
  assert (currIdx == haveSeenId.size() && currIdx == numNodes);
}





void getMaxBlastScoreById(const StrPair2Dbl & idpair2BlastScore, 
			  const VecStr & idList,
			  Str2Dbl & id2maxBlastScore)
{
  typedef StrPair2Dbl::const_iterator It1;

  id2maxBlastScore.clear();
  for (int i=0; i < idList.size(); i++)
    {
      id2maxBlastScore[idList.at(i)] = 1;
    }
  
  for (It1 it1=idpair2BlastScore.begin();
       it1 != idpair2BlastScore.end();
       it1++)
    {
      std::string id1 = it1->first.first;
      std::string id2 = it1->first.second;
      double v= it1->second;
      if (id2maxBlastScore[id1] < v)
	{
	  id2maxBlastScore[id1] = v;
	}
      if (id2maxBlastScore[id2] < v)
	{
	  id2maxBlastScore[id2] = v;
	}
    }
}




inline 	double getBlastScoreOrZero(const StrPair2Dbl & idpair2BlastScore, 
				   const std::string & i,
				   const std::string & j)
{
  typedef StrPair2Dbl::const_iterator It1;
  It1 it1;
  if (i<j)
    {
      it1 = idpair2BlastScore.find(StrPair(i,j));
    }
  else
    {
      it1 = idpair2BlastScore.find(StrPair(i,j));
    }
  if (it1==idpair2BlastScore.end())
    return 0;
  else
    return it1->second;
}



inline int getIdx(int i,int j, int imax, int jmax)
{
  return (i*jmax) + j;
}







inline void getSubscripts(int idx, int sp1nodes, int sp2nodes,
			  int & i, int & j)
{
  i = idx / sp2nodes;
  j = idx % sp2nodes;
}




double normalizeR_byL1norm(Hash_Int2Dbl & vec, VecBool & vecMembership)
{
  typedef Hash_Int2Dbl::iterator It1;
  double L1_sum = 0;
  for (It1 it1=vec.begin(); it1!=vec.end(); it1++)
    {
      L1_sum += fabs(it1->second);
    }
  std::cerr << "size before: " << vec.size() << " ";
  for (It1 it1=vec.begin(); it1!=vec.end(); )
    {
      it1->second /= L1_sum;
      if (fabs(it1->second) < EPS)
	{
	  It1 it1a = it1;
	  it1a++;
	  vecMembership[it1->first] = false;
	  vec.erase(it1);
	  it1 = it1a;
	}
      else
	{
	  it1++;
	}
    }
  std::cerr << "size after: " << vec.size() << std::endl;
  return L1_sum;
}





double pruneVecTopK(Hash_Int2Dbl & vec, 
		    VecBool & vecMembership,
		    int K=1000000)
{
  typedef Hash_Int2Dbl::iterator It1;
  std::string pruneStyle("top_K");
  int oldSize = vec.size();

  VecFloat vals(vec.size(),0);
  int i=0;
  for (It1 it1=vec.begin(); it1!=vec.end(); it1++)
    {
      vals[i] = float(it1->second);
      i++;
    }
  std::sort(vals.begin(), vals.end());
  
  if (vals.size() <= K)
    return vals[0];

  double threshVal = vals[vec.size()-K-1];
      
  It1 it1=vec.begin();
  while (it1 != vec.end())
    {
      if (it1->second < threshVal)
	{
	  It1 oldIt  = it1;
	  oldIt++;
	  vecMembership[it1->first] = false;
	  vec.erase(it1->first);
	  it1 = oldIt;
	}
      else
	{
	  it1++;
	}
    }
  std::cerr << "After pruning in style=" << pruneStyle
	    << " pruned " << oldSize << " entries to " << vec.size() 
	    << std::endl;
  return threshVal;
}




double change_L1norm(const VecBool & finalVecMembership, 
		      const Hash_Int2Dbl & finalVec,
		      const VecBool & oldFinalVecMembership, 
		      const Hash_Int2Dbl & oldFinalVec)
{
  double sumChange = 0;
  int Nchange = 1;
  assert (oldFinalVecMembership.size() == finalVecMembership.size());
  for (int i=0; i < oldFinalVecMembership.size(); i++)
    {
      if (!oldFinalVecMembership.at(i) || !finalVecMembership.at(i))
	continue;
      double d1 = finalVec.find(i)->second;
      double d2 = oldFinalVec.find(i)->second;
      sumChange += fabs(d1-d2);//*2/(d1+d2);
      //Nchange ++;
    }
  return sumChange;
  //return 100*(sumChange/Nchange);
}





void iterateRecursion(int K,
		      const VecStr & idList,
		      const Str2Int & id2idx,
		      const Str2Str & id2sp,
		      const std::string & sp1,
		      const std::string & sp2,
		      const StrPair2Dbl & idpair2IntxnScore,
		      const StrPair2Dbl & idpair2BlastScore,
		      VecBool & finalVecMembership,
		      Hash_Int2Dbl & finalVec,
		      VecBool & initVecMembership,
		      Hash_Int2Dbl & initVec,
		      VecInt & sp1_idxList, 
		      VecInt & sp2_idxList,
		      double convg_thresh,
		      double alpha,
		      bool doPruning,
		      int prunedVecSize)
{
  DEBUGP(201);
  int sp1edges, sp1nodes, sp2edges, sp2nodes, sp1graph_nodes, sp2graph_nodes;
  sp1edges = sp1graph_nodes = sp2edges = sp2graph_nodes = -1;
  sp1nodes = sp2nodes = -1;

  countEdgesNodesBySp(idpair2IntxnScore, id2sp, sp1, sp1graph_nodes, sp1edges);
  countEdgesNodesBySp(idpair2IntxnScore, id2sp, sp2, sp2graph_nodes, sp2edges);
  
  std::cerr << sp1edges << " " <<  sp1graph_nodes << " " 
	    << sp2edges << " " << sp2graph_nodes << std::endl;
  
  sp1_idxList.clear();
  sp2_idxList.clear();
  VecVecInt sp1_adjList, sp2_adjList;
  
  DEBUGP(202);

  /* 7/10/7 replaced these by getAdjListBySp_withAllGenes() which creates an 
   adj list  with all genes in the species, not just those with a PPI hit

  getAdjListBySp(idpair2IntxnScore, id2sp, id2idx,
		   sp1, sp1nodes, sp1_idxList, sp1_adjList);
  getAdjListBySp(idpair2IntxnScore, id2sp, id2idx,
		   sp2, sp2nodes, sp2_idxList, sp2_adjList);
  assert( sp1_idxList.size() == sp1nodes);
  assert( sp2_idxList.size() == sp2nodes);


  */

  getAdjListBySp_withAllGenes(idpair2IntxnScore, id2sp, id2idx,
			      sp1, sp1_idxList, sp1_adjList);
  getAdjListBySp_withAllGenes(idpair2IntxnScore, id2sp, id2idx,
			      sp2, sp2_idxList, sp2_adjList);



  sp1nodes = sp1_adjList.size();
  sp2nodes = sp2_adjList.size();

  std::cerr << "nodes vs graph_nodes: " 
	    << sp1nodes << " " <<  sp1graph_nodes << " " 
	    << sp2nodes << " " << sp2graph_nodes << std::endl;

  DEBUGP(202.5);
  Str2Dbl id2maxBlastScore;
  getMaxBlastScoreById(idpair2BlastScore, idList, id2maxBlastScore);

  DEBUGP(203);
  
  int R = sp1nodes * sp2nodes;
  
  finalVecMembership.clear();
  finalVec.clear();

  VecBool oldFinalVecMembership;
  Hash_Int2Dbl oldFinalVec;

  initVecMembership.clear();
  initVec.clear();
  
  
  IntPair P0(-1,-1);
  initVecMembership = VecBool(R, false);

  for (int i=0; i < sp1nodes; i++)
    for (int j=0; j < sp2nodes; j++)
      {
	int idx1, idx2, a, b, c ,d;
	idx1 = idx2 = a = b= c = d= -1;
	idx1 = getIdx(i, j, sp1nodes, sp2nodes);

	//DEBUGP(204);
	
	if (idx1 % 1000000 == 0) std::cerr << idx1 << std::endl;
	
	a = sp1_idxList[i];
	b = sp2_idxList[j];	
	double blst = getBlastScoreOrZero(idpair2BlastScore, idList.at(a),
					  idList.at(b));
	//DEBUGP(205);

	if (fabs(blst) > EPS)
	  {
	    initVecMembership[idx1] = true;
	    initVec[idx1] = blst / sqrt(id2maxBlastScore[idList.at(a)]*
					id2maxBlastScore[idList.at(b)]);
	  }
      }

  normalizeR_byL1norm(initVec, initVecMembership);
  
  finalVecMembership = initVecMembership;
  finalVec = initVec;
  
  //double sumChange = 0;
  //int Nchange = 1;

  double avgChange = 100;

  //BROKEN: limit on avgChange should be 2 or 0.5
  for (int iter=0; iter < K  && avgChange > convg_thresh; iter++)
    {
      double smallestNewValAllowed = 0;
      std::cerr << "K: " << iter << std::endl;

      if (iter>2)
	{
	  avgChange = change_L1norm(finalVecMembership, finalVec,
				    oldFinalVecMembership, oldFinalVec);
	  
	  //avgChange = sumChange/Nchange;
	  std::cerr << "curr size: " << finalVec.size() 
		    << " avg change in last iter: " << avgChange
		    << std::endl;
	}

      oldFinalVec = finalVec;
      oldFinalVecMembership = finalVecMembership;

      //sumChange = 0;
      //Nchange = 1;
      
      for (int i=0; i < sp1nodes; i++)
	for (int j=0; j < sp2nodes; j++)
	  {
	    int idx1,idx2,a,b,c,d;
	    idx1=idx2=a=b=c=d=-1;
	    idx1 = getIdx(i, j, sp1nodes, sp2nodes);
	    
	    double newVal = 0;

	    if (idx1 % 1000000 == 0) std::cerr << idx1 << std::endl;
	    //DEBUGP(205.5);

	    for (int u=0; u < sp1_adjList[i].size(); u++)
	      for (int v=0; v < sp2_adjList[j].size(); v++)
		{
		  c = sp1_adjList[i][u];
		  d = sp2_adjList[j][v];
		  
		  idx2 = getIdx(c, d, sp1nodes, sp2nodes);
		  double wt = 1.0 / (sp1_adjList[c].size() *
				     sp2_adjList[d].size());
		  //scoreAdjList[idx1].push_back(IntDblPair(idx2,wt));
		  if (finalVecMembership[idx2])
		    newVal += wt*finalVec[idx2];
		}


	    if (initVecMembership[idx1])
	      {
		newVal = alpha*newVal + (1-alpha)*initVec[idx1];
	      }
	    else
	      {
		newVal = alpha*newVal;
	      }

	    if (newVal > EPS && 
		(!doPruning || finalVecMembership[idx1]==true ||
		 newVal > smallestNewValAllowed))
	      {
		if (!finalVecMembership[idx1])
		  {
		    finalVecMembership[idx1] = true;
		    a = sp1_idxList[i];
		    b = sp2_idxList[j];		
		  }
		else if (finalVec[idx1] > EPS) 
		  {
		    ;
		    /*
		      double change = 
		      fabs(200*((newVal - finalVec[idx1])/
				(newVal + finalVec[idx1])));
		    sumChange += change;
		    Nchange ++;
		    */
		  }
		finalVec[idx1] = newVal;
	      }
	    else
	      {
		finalVecMembership[idx1] = false;
	      }
	  }


      if (doPruning)
	{
	  int MAX_VEC_SIZE = prunedVecSize;
	  double vecSmallestVal =  pruneVecTopK(finalVec, finalVecMembership, 
						MAX_VEC_SIZE);
	  if (iter <= 1)
	    {
	      smallestNewValAllowed = 0.1* vecSmallestVal;
	    }
	  else
	    {
	      smallestNewValAllowed = vecSmallestVal;
	    }
	}
      
      normalizeR_byL1norm(finalVec, finalVecMembership);
  
    }
  normalizeR_byL1norm(finalVec, finalVecMembership);
  DEBUGP(210);
  
}


/* convert values in vec to StrPair2Dbl entries, by matching with idList
   and sp1_idxList & sp2_idxList
*/
int readMatchScores(const VecStr & idList, 
		    const VecBool & vecMembership,
		    const Hash_Int2Dbl &  vec, 
		    const VecInt & sp1_idxList, 
		    const VecInt & sp2_idxList, 
		    StrPair2Dbl & idpair2MatchScore)
{
  typedef Hash_Int2Dbl::const_iterator It1;
  
  int i,j;
  i=j=-1;

  int sp1nodes = sp1_idxList.size();
  int sp2nodes = sp2_idxList.size();
  for(It1 it1=vec.begin(); it1!=vec.end(); it1++)
    {
      getSubscripts(it1->first, sp1nodes, sp2nodes, i,j);
      StrPair spair(idList[sp1_idxList[i]], idList[sp2_idxList[j]]);
      idpair2MatchScore[spair] = it1->second;
    }
  return 1;
}




void writeOutScores(const StrPair2Dbl & idpair2Score,
		    const std::string & outScoreFile)
{
  typedef StrPair2Dbl::const_iterator It;

  std::ofstream os(outScoreFile.c_str());

  for (It it=idpair2Score.begin(); it!=idpair2Score.end(); it++)
    {
      os << it->first.first << " " << it->first.second << " " 
	 << it->second << "\n";
    }
  os.close();

}




void  writeClusterData(const Str2Int & id2clstr, 
		       const std::string & outClstrFile)
{
  typedef Str2Int::const_iterator It;

  std::ofstream os(outClstrFile.c_str());

  //check if id2clstr is not specified
  if (id2clstr.size()==0 || id2clstr.begin()->second == -1)
    {
      std::cerr << "No output cluster\n";
      os.close();
      return;
    }

  for (It it=id2clstr.begin(); it!=id2clstr.end(); it++)
    {
      os << it->first << " " << it->second << "\n";
    }
  os.close();
}



int readMatchScoreFile(const std::string & matchScoreFile, 
		       const VecStr & idList, 
		       StrPair2Dbl & idpair2MatchScore)
{
  int numVals = 0;
  std::ifstream ifs(matchScoreFile.c_str());
  
  StrSet idSet;
  for (int i=0; i < idList.size(); i++)
    {
      idSet.insert(idList.at(i));
    }

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
	idpair2MatchScore[StrPair(id1,id2)] = val;
      else
	idpair2MatchScore[StrPair(id2,id1)] = val;

      numVals++;
    }
  return numVals;
}




void  computeSumOfCommonScores(const StrPair2Dbl & idpair2BlastScore, 
			       const StrPair2Dbl & idpair2MatchScore,
			       double & sumCommonScoresInMatch, 
			       double & sumCommonScoresInBlast)
{
  typedef StrPair2Dbl::const_iterator  It1;
  typedef StrPair2Dbl::iterator It2;
  sumCommonScoresInBlast = 0;
  sumCommonScoresInMatch = 0;

  for (It1 it1=idpair2BlastScore.begin(); it1!=idpair2BlastScore.end(); it1++)
    {
      It1 it_a = idpair2MatchScore.find(it1->first);
      if (it_a != idpair2MatchScore.end())
	{
	  sumCommonScoresInBlast += it1->second;
	  sumCommonScoresInMatch += it_a->second;
	}
    }
  std::cerr << "total of common scores: blast: " << sumCommonScoresInBlast
	    << " match: " << sumCommonScoresInMatch << std::endl;
}






void fillInMissingBlastScores_L1norm_A(StrPair2Dbl & idpair2MatchScore, 
				     double alpha,
				     const StrPair2Dbl & idpair2IntxnScore,
				     const StrPair2Dbl & idpair2BlastScore, 
				     const VecStr & idList,
				     const Str2Str & id2sp,
				     const Str2Int & id2idx,
				     const VecStr & spList)
{
  typedef StrPair2Dbl::const_iterator It1;
  typedef StrPair2Dbl::iterator It2;

  Str2Dbl id2maxBlastScore;
  getMaxBlastScoreById(idpair2BlastScore, idList, id2maxBlastScore);
  
  DEBUGP(501);

  StrPair2Dbl idpair2normalizedBlastScore;
  for (It1 it1=idpair2BlastScore.begin(); it1!=idpair2BlastScore.end(); it1++)
    {
      idpair2normalizedBlastScore[it1->first] = it1->second / 
	sqrt( id2maxBlastScore[it1->first.first] *
	      id2maxBlastScore[it1->first.second]);
    }
  
  DEBUGP(502);
  
  double sumCommonScoresInMatch = 0;
  double sumCommonScoresInBlast = 0;
  computeSumOfCommonScores(idpair2normalizedBlastScore, idpair2MatchScore,
			   sumCommonScoresInMatch, sumCommonScoresInBlast);

  double scale = sumCommonScoresInBlast / sumCommonScoresInMatch;

  double totalMatchScores = 0;
  for (It2 it2=idpair2MatchScore.begin(); it2 !=idpair2MatchScore.end(); it2++)
    {
      it2->second =  it2->second * scale;
      totalMatchScores += it2->second;
    }

  std::cerr << "Total Match Scores before filling in Blast: " 
	    << totalMatchScores << std::endl;
  
  int numScoresFilledIn = 0;

  for (It1 it1=idpair2normalizedBlastScore.begin();
       it1 != idpair2normalizedBlastScore.end();
       it1++)
    {
      std::string spA = id2sp.find(it1->first.first)->second;
      std::string spB = id2sp.find(it1->first.second)->second;
      
      if (spA==spB)
	continue;
      
      if (idpair2MatchScore.find(it1->first) != idpair2MatchScore.end())
	{
	  continue;
	}
      
      idpair2MatchScore[it1->first] = it1->second;
      
      totalMatchScores += it1->second;
      numScoresFilledIn += 1;
    }

  DEBUGP(506);
  std::cerr << "# of Blast scores added in : " << numScoresFilledIn
	    << " total score now: " << totalMatchScores
	    << std::endl;
  
}


void fillInMissingBlastScores_L1norm(StrPair2Dbl & idpair2MatchScore, 
				     double alpha,
				     const StrPair2Dbl & idpair2IntxnScore,
				     const StrPair2Dbl & idpair2BlastScore, 
				     const VecStr & idList,
				     const Str2Str & id2sp,
				     const Str2Int & id2idx,
				     const VecStr & spList)
{
  typedef StrPair2Dbl::const_iterator It1;
  

  Str2Dbl id2maxBlastScore;
  getMaxBlastScoreById(idpair2BlastScore, idList, id2maxBlastScore);
  
  DEBUGP(501);

  StrPair2Dbl idpair2normalizedBlastScore;
  for (It1 it1=idpair2BlastScore.begin(); it1!=idpair2BlastScore.end(); it1++)
    {
      idpair2normalizedBlastScore[it1->first] = it1->second / 
	sqrt( id2maxBlastScore[it1->first.first] *
	      id2maxBlastScore[it1->first.second]);
    }
  
  DEBUGP(502);

  for (int i_sp=0; i_sp < spList.size(); i_sp++)
    for (int j_sp=i_sp+1; j_sp < spList.size(); j_sp++)
      {
	  std::string sp1=spList.at(i_sp);
	  std::string sp2=spList.at(j_sp);
	  
	  int sp1edges, sp1nodes, sp2edges, sp2nodes;
	  sp1edges = sp1nodes = sp2edges = sp2nodes = -1;
	  countEdgesNodesBySp(idpair2IntxnScore, id2sp, sp1, 
			      sp1nodes, sp1edges);
	  countEdgesNodesBySp(idpair2IntxnScore, id2sp, sp2, 
			      sp2nodes, sp2edges);

	  VecInt sp1_idxList, sp2_idxList;
	  VecVecInt sp1_adjList, sp2_adjList;

	  DEBUGP(503);

	  getAdjListBySp(idpair2IntxnScore, id2sp, id2idx,
			 sp1, sp1nodes, sp1_idxList, sp1_adjList);
	  getAdjListBySp(idpair2IntxnScore, id2sp, id2idx,
			 sp2, sp2nodes, sp2_idxList, sp2_adjList);
	  DEBUGP(5031);

	  assert( sp1_idxList.size() == sp1nodes);
	  assert( sp2_idxList.size() == sp2nodes);
	  
	  DEBUGP(504);

	  int R = sp1nodes * sp2nodes;

	  double added_scores =0;
	  
	  VecDbl initVec;
	  for (int i=0; i < sp1nodes; i++)
	    for (int j=0; j < sp2nodes; j++)
	      {
		int idx1, idx2, a, b, c ,d;
		idx1 = idx2 = a = b= c = d= -1;
		idx1 = getIdx(i, j, sp1nodes, sp2nodes);
		
	
		if (idx1 % 1000000 == 0) std::cerr << "fillBlast:" 
						   << idx1 << std::endl;
	
		a = sp1_idxList[i];
		b = sp2_idxList[j];	
		double blst = getBlastScoreOrZero(idpair2BlastScore, 
						  idList.at(a),
						  idList.at(b));

		if (fabs(blst) > EPS)
		  {
		    initVec.push_back(blst / sqrt(id2maxBlastScore[idList.at(a)]*
						  id2maxBlastScore[idList.at(b)]));

		  }
	      }

	  DEBUGP(505);
	  
	  /* need to compute initVec, separately of idpair2normalizedBlastScore
	     because the L1 normalization should only be done for those scores
	     which were present in the network alignment phase
	  */

	  double sum_of_initVec=0;
	  for (int r=0; r<initVec.size(); r++)
	    {
	      sum_of_initVec += initVec[r];
	    }
  
	  double scale_factor = 1.0/sum_of_initVec;
	  
  
	  for (It1 it1=idpair2normalizedBlastScore.begin();
	       it1 != idpair2normalizedBlastScore.end();
	       it1++)
	    {
	      std::string spA = id2sp.find(it1->first.first)->second;
	      std::string spB = id2sp.find(it1->first.second)->second;
	      
	      if (spA==spB || 
		  (spA != sp1 && spA != sp2) || 
		  (spB != sp1 && spB != sp2))
		continue;

	      if (idpair2MatchScore.find(it1->first) != 
		  idpair2MatchScore.end())
		{
		  std::cerr << "CHECK! this score shouldn't be here: " 
			    << it1->first.first << " " << it1->first.second
			    << std::endl;
		  continue;
		}
	      
	      idpair2MatchScore[it1->first] = (1-alpha)*scale_factor*(it1->second);
	      added_scores +=  (1-alpha)*scale_factor*(it1->second);
	    }
	  DEBUGP(506);
	  std::cerr << "filling in blast scores for " << sp1 << " " 
		    << sp2 << ": " << added_scores << std::endl;
      }
}





int main(int argc, char ** argv)
{
  typedef Str2Int::const_iterator It1;
  typedef StrPair2Dbl::const_iterator It2;
  typedef Str2Str::const_iterator It3;
  typedef StrPair2Dbl::iterator It4;

  std::string dataDescFile;
  std::string filePrefix = "";
  std::string sp1in, sp2in;
  int K = 20;
  double convg_thresh = 1;
  double alpha = 0.9;
  int prunedVecSize = -1;
  double min_primary_fraction = 0.1; 
  double min_secondary_fraction = 0.8;
  int max_per_species = 5;

  std::string strategy("hsp");
  std::string matchScoreFile;
  Str2Int id2idx;

  VecStr spList, idList;
  Str2Str id2sp;
  StrPair2Dbl idpair2BlastScore, idpair2IntxnScore, idpair2MatchScore;
  
  //get cmdline options

  getoptions(argc, argv, 
	     K, convg_thresh, alpha, prunedVecSize, 
	     min_primary_fraction, min_secondary_fraction, max_per_species,
	     strategy, filePrefix, dataDescFile, matchScoreFile);

  std::cerr << "Command line summary: "
	    << " K: " << K << " convergence_thresh: " << convg_thresh << "% "
	    << " alpha: " << alpha 
	    << " doPruning: " << (prunedVecSize<=0?"false":"true")
	    << " prune_vec_size: " << prunedVecSize
	    << " min_secondary_fraction: " << min_secondary_fraction
	    << " min_primary_fraction: " << min_primary_fraction
	    << " max_per_species: " << max_per_species
	    << " strategy: " << strategy
	    << " dataDescFile: " << dataDescFile
	    << " match-score-file: '" << matchScoreFile << "'" 
	    << std::endl;

  DEBUGP(1);
  grabSequenceAndInteractionData(dataDescFile, 
				 spList, 
				 idList, 
				 id2sp,
				 idpair2BlastScore, 
				 idpair2IntxnScore);

  std::cerr << "size of sp-list: " << spList.size() << std::endl;

  int numBetweenSpBlastScores = 0;
  for (It2 it2=idpair2BlastScore.begin(); it2!=idpair2BlastScore.end(); it2++)
    {
      if (id2sp.find(it2->first.first)->second != id2sp.find(it2->first.second)->second)
	numBetweenSpBlastScores ++;
    }

  std::cerr << "number of blast scores: " << idpair2BlastScore.size()
	    << " # scores between species: " << numBetweenSpBlastScores
	    << std::endl;

  bool zeroOneBlastScores = false;
  if (zeroOneBlastScores)
    {
      for (It4 it4=idpair2BlastScore.begin(); it4!=idpair2BlastScore.end();
	   it4++)
	{
	  it4->second = 1.0;
	}
    }
  
  
  //create a id2int map
  std::sort(idList.begin(), idList.end());
  for (int i=0; i < idList.size(); i++)
    {
      id2idx[idList[i]] = i;
    }



  DEBUGP(101);
  idpair2MatchScore.clear();


  if (matchScoreFile.size()>0)
    {
      DEBUGP(102);
      readMatchScoreFile(matchScoreFile, idList, idpair2MatchScore);
      DEBUGP(103);
    }
  else
    {
      for (int spIdx1=0; spIdx1 < spList.size(); spIdx1++)
	for (int spIdx2 = spIdx1+1; spIdx2 < spList.size(); spIdx2++)
	  {
	    if (spList[spIdx1] == spList[spIdx2])
	      continue;
  
	    VecBool finalVecMembership, initVecMembership;
	    Hash_Int2Dbl finalVec, initVec;
	    VecInt sp1_idxList, sp2_idxList;
	    //Int2IntPair finalCoords, initCoords;

	    DEBUGP(2);
	
	    double doPruning = false;
	    if (prunedVecSize > 0)
	      doPruning = true;

	    iterateRecursion(K,
			     idList,
			     id2idx,
			     id2sp,
			     spList[spIdx1], spList[spIdx2],
			     idpair2IntxnScore,
			     idpair2BlastScore,
			     finalVecMembership,
			     finalVec,
			     initVecMembership,
			     initVec,
			     sp1_idxList,
			     sp2_idxList,
			     convg_thresh,
			     alpha,
			     doPruning,
			     prunedVecSize);
	
	    //save the match scores
	    DEBUGP(4);
	    int numScoresInFinalNotInInit = 0;
	    for (int i=0; i < initVecMembership.size(); i++)
	      {
		if (!initVecMembership[i] && finalVecMembership[i])
		  numScoresInFinalNotInInit++;
	      }

	    readMatchScores(idList, finalVecMembership, finalVec, 
			    sp1_idxList, sp2_idxList, idpair2MatchScore);
	    
	    std::cerr << spList[spIdx1] << " " << spList[spIdx2]
		      << " # scores in finalVec but not in initVec: "
		      << numScoresInFinalNotInInit
		      <<  " init: " << initVec.size() 
		      << " final: " << finalVec.size() 
		      << std::endl;
	    std::cerr << spList[spIdx1] << " " << spList[spIdx2]
		      << " # cumulative scores matchscores: "
		      << idpair2MatchScore.size() << std::endl;
	     

	  }

      DEBUGP(5);


      /* 7/10/7: changed iterateRecursion to directly includes genes with
	 no PPI hits, so blast scores shouldn't be missing 

      fillInMissingBlastScores_L1norm_A(idpair2MatchScore, 
					alpha,
					idpair2IntxnScore,
					idpair2BlastScore, 
					idList,
					id2sp,
					id2idx,
					spList);
      */

      char matchScoreFile[200];
      sprintf(matchScoreFile,"%s_match-score.txt",filePrefix.c_str());
      writeOutScores(idpair2MatchScore, matchScoreFile);
    }
  DEBUGP(6);

  //do k-partite matching
  Str2Int id2clstr;
  doKpartiteMatching(idpair2MatchScore, idList, id2idx, 
		     id2sp, spList, 
		     min_primary_fraction, 
		     min_secondary_fraction, 
		     max_per_species,
		     id2clstr);
  
  DEBUGP(7);
  //write out stuff;
  char outClstrFile[200];
  sprintf(outClstrFile, "%s_final_cluster.txt", filePrefix.c_str());
  writeClusterData(id2clstr, outClstrFile);
}
