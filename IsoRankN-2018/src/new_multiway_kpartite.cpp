// Date : 2017.10.14
/* memo :
  we use NSD algorithm and the skill of Parallel calculation in order to 
  accelerate the computation of similarity score. Beside,with the help of
  'armadillo' library, we can deal with matrix multiplication in an efficient 
  way
*/
#ifndef INCLUDED_defs
#include "defs.h"
#endif
#ifndef INCLUDED_utils
#include "Utils.h"
#endif
#ifndef INCLUDED_k_partite
#include "k_partite.h"
#endif
#ifndef INCLUDED_shou_algorithm
#include "shou_algorithm.h"
#endif

#include <armadillo>  
#include <math.h>
#include <algorithm>
#include <iostream>
#include <time.h>
#include <omp.h>
/* TO IMPLEMENT:
 
*/
using namespace std;
void usage(const std::string & errMsg="")
{
  std::string usageStr="Usage: multiway_kpartite [--prefix <file-prefix>] [--scorefile <match-score-file>] [--K <max_iter>] [--thresh <in_percent>] [--alpha <weight_of_ntwk_term>] [--maxveclen <max_len_of_eigenvector>] [--cl_min_pri_thresh <thresh1>] [--cl_min_sec_thresh <thresh2>] [--cl_max_per_sp <max1>] [--readstar] [--reverse] [--o <final-cluster-file>] <data-description.txt> \n"
    "\n"
    "--prefix : prefix of the output .txt and .dat files (Default = '').\n\n"
    "--scorefile <match-score-file> : the file containing pairwise scores\n"
    " after doing network alignment. Forms input to k-partite clustering.\n\n"
    "--K <maxiter>: max number of iterations of the power method (>3, <100).\n"
    "--thresh <pct-thresh>: the algorithm stops when the average percentage \n"
    "\tchange in value of Rij is less than <pct-thresh>\n"
    "--cl_min_pri_thresh <thresh1> : in the k_partite step, the fractional \n"
    "\tthreshold to use during constructing the primary cluster.\n"
    "--cl_min_sec_thresh <thresh2> : in the k_partite step, the fractional \n"
    "\tthreshold to use when adding the secondary nodes.\n"
    "--cl_max_per_sp <max1> : in the k_partite step, the max number \n"
    "\t of species per cluster.\n"
    "The program stops when either of the above limits are violated.\n\n"

    "--alpha: weight of network data (in [0,1])\n"

    "--maxveclen: for really large graphs, limits the number of non-zero \n"
    "\tentries in R.\n"
    "--readstar: loading exist star files.\n"
    "--reverse: use reverse mode.\n"
    "--o <final-cluster-file>: the final output cluster file, the default file name is\n"
    "\t\"tmp_output_cluster.txt\".\n\n"
    "The argument (data-description.txt) should specify the list of species\n"
    "\tto be used as well as the path to the various files. Look at src/data.in\n"
    "X.tab should have a first line as 'INTERACTOR_A INTERACTOR_B'\n"
    "\tEach following line should be of type 'id1 id2 wt', where 3rd column\n" 
    "\tis optional.\n"
    "\n"
    "X-Y.evals should contain the non-network information. Each line is of \n"
    "\tthe form 'id_X id_Y score'. The scores can be BLAST Bit-scores or \n"
    "\tsomething else (even binary values!). This is also used to initialize\n"
    "\tthe power method. You can set alpha=1 and specify random values \n"
    "\there-- the file will then only be used to set the initial value \n"
    "\tof R in the power method.\n"
    "\n"
    "X1, X2, X3... MUST NOT SHARE ANY NODE LABELS! Thus, if you are \n"
    "\tcomparing two (or more) conformations of a single protein, rename all nodes \n"
    "\tin one of the conformations (e.g., by prefixing/suffixing with \n"
    "\tsomething).\n";
  if (errMsg.size() > 0) 
    {
      std::cerr << "Error: " << errMsg << "\n\n";
    }
  std::cerr << usageStr;
}


//get the information from cmd
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
		std::string & matchScoreFile,
		bool & starflag,
		bool & reverseflag, 
		std::string & outputClusterFile,
		int & threadid)
{
  struct option longopts[] = {
    {"prefix",1,0,0},
    {"K",1,0,0},
    {"threadid",1,0,0},
    {"thresh",1,0,0},
    {"alpha",1,0,0},
    {"maxveclen",1,0,0},
    {"scorefile",1,0,0},
    {"cl_min_pri_thresh",1,0,0},
    {"cl_min_sec_thresh",1,0,0},
    {"cl_max_per_sp",1,0,0},
    {"help",0,0,0},
    {"readstar",0,0,0},
    {"reverse",0,0,0},
    {"o",1,0,0},
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
  starflag = false;
  reverseflag = false;
  threadid=4;
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
	   if (optindex == 2) //"K"
	  {
	   threadid = atoi(optarg);
	 }
      if (optindex == 3) //"thresh"
	{
	  convg_thresh = double(atof(optarg));
	  if (convg_thresh <= 0 || convg_thresh > 10)
	    {
	      usage("Need convg_thresh in (0,10)");
	      exit(1);
	    }	  
	}
      if (optindex == 4) //"alpha"
	{
	  alpha = double(atof(optarg));
	  if (alpha < 0 || alpha > 1)
	    {
	      usage("Bad 'alpha'");
	      exit(1);
	    }	 
	}
      if (optindex == 5) //"prunedVecSize"
	{
	  prunedVecSize = atoi(optarg);
	  if (prunedVecSize <= 0)
	    {
	      usage("Bad 'prunedvecsize' option");
	      exit(1);
	    }
	}
      if (optindex == 6) //"scorefile"
	{
	  matchScoreFile = optarg;
	}
      if (optindex == 7) //"cl_min_pri_thresh"
	{
	  min_primary_fraction = atof(optarg);
	  if (min_primary_fraction <= 0 || min_primary_fraction > 1)
	    {
	      usage("Bad 'cl_min_pri_thresh' option");
	      exit(1);
	    }
	}
      if (optindex == 8) //"cl_min_sec_thresh"
	{
	  min_secondary_fraction = atof(optarg);
	  if (min_secondary_fraction <= 0 || min_secondary_fraction > 1)
	    {
	      usage("Bad 'cl_min_sec_thresh' option");
	      exit(1);
	    }
	}
      if (optindex == 9) //"cl_max_per_sp"
	{
	  max_per_species = atoi(optarg);
	  if (max_per_species < 1 )
	    {
	      std::cerr << "blah: '" << optarg << "' " << max_per_species << std::endl;
	      usage("Bad 'cl_max_per_sp' option");
	      exit(1);
	    }
	}
      if (optindex == 10) //"help"
	{
	  usage("");
	  exit(1);
	}
      if(optindex == 11)  
    {
      starflag = true;
    }
     if(optindex == 12)
    { 
      reverseflag = true;
    }
     if (optindex == 13) //"outputClusterFile"
	{
	  outputClusterFile = optarg;
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

      //DEBUGP(302);#include <time.h>
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
	vecMembership[it1->first] = true;	
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
struct matrix
{
int row;
int col;
double value;	

};

struct blastscore
{
std::string name1;
std::string name2;
int no1;
int no2;
double blast;
};
struct worklist
{
int number;
int graphA;
int graphB;
int sp1edges;
int sp1nodes;
int sp2edges;
int sp2nodes; 
int sp1graph_nodes;
int sp2graph_nodes;
VecInt sp1_idxList, sp2_idxList;
VecVecInt sp1_adjList; 
VecVecInt sp2_adjList;
};

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
void writematrix(std::vector<matrix>Matrix,char name[])
{
	std::ofstream os(name);
	for(int i=0;i<Matrix.size();i++)
	{
		os<<Matrix[i].row<<" "<<Matrix[i].col<<" "<<Matrix[i].value<<"\n";
	}
	os.close();
	
}
void WriteDynamicData(VecStr &spList, VecStr &idList,std::vector<worklist> work,VecInt numberList,std::vector<VecVecInt> adjancymatrix)
{
	  typedef VecStr::const_iterator It;
      std::ofstream os("Temp/NameTable.txt");
      for (It it=spList.begin(); it!=spList.end(); it++)
       {
         os <<*it<< "\n";
       }
       os.close();
	   std::ofstream os1("Temp/NoTable.txt");
	   for (It it=idList.begin(); it!=idList.end(); it++)
       {
         os1<<*it<< "\n";
       }
       os1.close();
       std::ofstream os2("Temp/Worklist.txt");
       for (int i=0;i<work.size();i++)
       {
         os2<<i<<" "<<work[i].graphA<<" "<<work[i].graphB<<" "<<work[i].sp1nodes<<" "<<work[i].sp2nodes<<" "<<work[i].sp1graph_nodes<<" "<<work[i].sp2graph_nodes<< "\n";
       }
	   os2.close();
	   std::ofstream os3("Temp/NumberTable.txt");
	   os3<<numberList[0]<<"\n";
	   for (int i=1;i<numberList.size();i+=2)
       {
		   os3<<numberList[i]<<"\n";
		}
       os3.close();	
       std::ofstream os4("Temp/Matrix0.txt");
	   for(int i=0;i<adjancymatrix[0].size();i++)
	   {
           os4<<i<<" ";
          for(int j=0;j<adjancymatrix[0][i].size();j++)
          {
           os4<<adjancymatrix[0][i][j]<<" ";
	      }
	       os4<<"\n";
	   }
	   char jonsnow[300]="";
	   for(int v=1;v<adjancymatrix.size();v+=2)
	   {
	     sprintf(jonsnow,"Temp/Matrix%d.txt",v);
	     std::ofstream os5(jonsnow);
	     for(int i=0;i<adjancymatrix[v].size();i++)
	      {
            os5<<i<<" ";
            for(int j=0;j<adjancymatrix[v][i].size();j++)
            {
              os5<<adjancymatrix[v][i][j]<<" ";
	        }
	        os5<<"\n";
 	      }  
	     os5.close();	
	   }
	  

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
  if(ifs.fail())
  {
   std::cerr<< "The match score file is not exist~!!!\n";
   assert(0);
  }
  if(ifs.peek() == 32 || ifs.peek() == -1 || ifs.peek() == 9)
  {
   std::cerr<< "The match score file is empty~!!!\n";
   assert(0);	
  }
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
      //if ( idSet.find(id1) == idSet.end() ||
      //   idSet.find(id2) == idSet.end())
      //	std::cerr << id1 << " "<< id2 << " not in idSet" << std::endl;
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
    {      continue;
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
   printf("Press ANY key to RUN");
   fgetc(stdin);
   // Parameter declaration
   std::string dataDescFile;
   std::string filePrefix = "";
   std::string sp1in, sp2in;
   int K = 20;
   int threadid=4;
   double convg_thresh = 1;
   double alpha = 0.9;
   int prunedVecSize = -1;
   double min_primary_fraction = 0.1; 
   double min_secondary_fraction = 0.8;
   int max_per_species = 5;
   bool starflag = false;
   bool reverseflag = false;
   std::string strategy("hsp");
   std::string matchScoreFile;
   std::string outputClusterFile("tmp_output_cluster.txt");
   std::string readstarstate("No");
   std::string reversestate("No");
   Str2Int id2idx;
   //std::cout<<outputClusterFile<<" "<<strategy<<std::endl; 
   VecStr spList, idList;
   VecInt numberList;
   Str2Str id2sp;
   StrPair2Dbl idpair2BlastScore, idpair2IntxnScore, idpair2MatchScore,idpair2MatchScoreNew;
   //Get Cmdline options
   getoptions(argc,argv,K,convg_thresh,alpha,prunedVecSize,min_primary_fraction,min_secondary_fraction,max_per_species,strategy,filePrefix,dataDescFile,matchScoreFile,starflag,reverseflag,outputClusterFile,threadid);
   if(starflag){readstarstate = "Yes";}
   if(reverseflag){reversestate = "Yes";}
   std::cerr << "\033[1;36mCommand line summary\033[0m: "<<std::endl
	    << "\n\033[1;36mK \033[0m: " << K << "\n\033[1;36mconvergence_thresh\033[0m: " << convg_thresh << "% "
	    << "\n\033[1;36malpha\033[0m: " << alpha 
	    << "\n\033[1;36mdoPruning\033[0m: " << (prunedVecSize<=0?"false":"true")
	    << "\n\033[1;36mprune_vec_size\033[0m: " << prunedVecSize
	    << "\n\033[1;36mmin_secondary_fraction\033[0m: " << min_secondary_fraction
	    << "\n\033[1;36mmin_primary_fraction\033[0m: " << min_primary_fraction
	    << "\n\033[1;36mmax_per_species\033[0m: " << max_per_species
	    << "\n\033[1;36mstrategy\033[0m: " << strategy
	    << "\n\033[1;36mdataDescFile\033[0m: " << dataDescFile
	    << "\n\033[1;36mmatch-score-file\033[0m: "<< matchScoreFile 
	    << "\n\033[1;36mread_star_files\033[0m: "<< readstarstate 
	    << "\n\033[1;36mreverse_mode\033[0m: "<< reversestate
	    << "\n\033[1;36mfinal-cluster-file\033[0m: "<< outputClusterFile
	    << "\n\033[1;36mThread number\033[0m: "<< threadid
	    << std::endl;
   DEBUGP(1);
   if (outputClusterFile.size() == 0)
   {
  	 std::cerr<<"Output cluster file name is null\n";
     assert(0);
   }
   // Check whether Similarity data is given

   if (matchScoreFile.size()> 0) //we have Rij, we don't need idpair2BlastScore and idpair2IntxnScore
   { 
    grabSequenceData(dataDescFile,spList,idList,id2sp);
    std::cerr << "size of sp-list: " << spList.size() << std::endl;
    std::cerr << "size of id_list: " << idList.size() << std::endl;
   }
   else
   {
     grabSequenceAndInteractionData(dataDescFile,spList,idList,id2sp,idpair2BlastScore,idpair2IntxnScore);
     std::cerr << "size of sp-list: " << spList.size() << std::endl;
     int numBetweenSpBlastScores = 0;
     for (It2 it2=idpair2BlastScore.begin(); it2!=idpair2BlastScore.end(); it2++)
     {
      if (id2sp.find(it2->first.first)->second != id2sp.find(it2->first.second)->second)numBetweenSpBlastScores++;
     }
     std::cerr << "number of blast scores: " << idpair2BlastScore.size()<< " # scores between species: " << numBetweenSpBlastScores << std::endl;
     bool zeroOneBlastScores = false;
     if (zeroOneBlastScores)
     {
      for (It4 it4=idpair2BlastScore.begin(); it4!=idpair2BlastScore.end();it4++){it4->second = 1.0;}
     }
   }
   DEBUGP(101);
   //create a id2int map
   std::sort(idList.begin(), idList.end());
   for (int i=0; i < idList.size(); i++){id2idx[idList[i]] = i;}
   idpair2MatchScore.clear();
   std::cout<<"Creating WorkList...."<<std::endl;
   // Build the parallel task list
   std::vector<worklist> work;
   std::vector<VecVecInt> adjancymatrix;
   adjancymatrix.clear();
   int workcounter;
   for (int spIdx1=0; spIdx1 < spList.size(); spIdx1++)
   {
	for (int spIdx2 = spIdx1; spIdx2 < spList.size(); spIdx2++)
	  {
	    if (spList[spIdx1] == spList[spIdx2]){}	     
	    else
	     {		
			 // Declaration :
			 worklist tempwork;
			 int sp1edges, sp1nodes, sp2edges, sp2nodes, sp1graph_nodes, sp2graph_nodes;
             sp1edges = sp1graph_nodes = sp2edges = sp2graph_nodes = -1;
             sp1nodes = sp2nodes = -1;
             VecInt sp1_idxList, sp2_idxList;
	         VecVecInt sp1_adjList; 
		     VecVecInt sp2_adjList;
			 countEdgesNodesBySp(idpair2IntxnScore, id2sp, spList[spIdx1], sp1graph_nodes, sp1edges);
             countEdgesNodesBySp(idpair2IntxnScore, id2sp, spList[spIdx2], sp2graph_nodes, sp2edges); 
             getAdjListBySp_withAllGenes(idpair2IntxnScore, id2sp, id2idx,spList[spIdx1], sp1_idxList, sp1_adjList);
             getAdjListBySp_withAllGenes(idpair2IntxnScore, id2sp, id2idx,spList[spIdx2], sp2_idxList, sp2_adjList);
			 sp1nodes = sp1_adjList.size();
             sp2nodes = sp2_adjList.size();
			 if(spIdx1==0)
			 {
				numberList.push_back(sp1nodes);
				numberList.push_back(sp2nodes);
				adjancymatrix.push_back(sp1_adjList);
				adjancymatrix.push_back(sp2_adjList);
			 }
			 //PuT Property
			 tempwork.number=workcounter; //Record the number
			 tempwork.graphA=spIdx1;
			 tempwork.graphB=spIdx2;
             tempwork.sp1edges=sp1edges;
             tempwork.sp2edges=sp2edges;
             tempwork.sp1graph_nodes=sp1graph_nodes;
             tempwork.sp2graph_nodes=sp2graph_nodes;
             tempwork.sp1nodes=sp1nodes;
             tempwork.sp2nodes=sp2nodes;
             tempwork.sp1_idxList=sp1_idxList;
             tempwork.sp2_idxList=sp2_idxList;
             tempwork.sp1_adjList=sp1_adjList;
             tempwork.sp2_adjList=sp2_adjList;
             //std::cerr << sp1edges << " " <<  sp1graph_nodes << " " << sp2edges << " " << sp2graph_nodes << std::endl; 	
             //std::cerr << "nodes vs graph_nodes: " << sp1nodes << " " <<  sp1graph_nodes << " " << sp2nodes << " " << sp2graph_nodes << std::endl;
             //Put into workList
             work.push_back(tempwork);		 
		     workcounter=workcounter+1;
		 }
	  }
   }
   DEBUGP(102);
   //Writing Some data to save
   WriteDynamicData(spList,idList,work,numberList,adjancymatrix);
   std::cout<<"Adjusting Sequence Data... "<<std::endl; 
   DEBUGP(103);

    // Building Blast Score
   StrPair2Dbl::const_iterator it5;
   Str2Int::const_iterator it6,it7;
   Str2Dbl id2maxBlastScore; 
   getMaxBlastScoreById(idpair2BlastScore, idList, id2maxBlastScore);
   std::vector<blastscore> theblast;
   blastscore temp;
   for(it5 = idpair2BlastScore.begin(); it5 != idpair2BlastScore.end(); it5++)
   {
	 temp.name1=it5->first.first;
	 temp.name2=it5->first.second;
	 temp.blast=(it5->second)/sqrt(id2maxBlastScore[(it5->first.first)]*id2maxBlastScore[(it5->first.second)]);
	 theblast.push_back(temp);
    }
	for(int i=0;i<theblast.size();i++)
    {
     it6=id2idx.find(theblast[i].name1);
     theblast[i].no1=it6->second;
     it7=id2idx.find(theblast[i].name2);
     theblast[i].no2=it7->second;
    }
     
    DEBUGP(104);
    std::cout<<"----- Computing -----"<<std::endl;
    int threadnumber;
    int count,round,iter ;
    double doPruning= true;
    arma::sp_mat A;
    arma::sp_mat B;
    arma::sp_mat H;
    arma::sp_mat spanswer;
    arma::mat fanswer;
    
    int sp1_idxListre[idList.size()]; 
    int sp2_idxListre[idList.size()]; 
    int numScoresInFinalNotInInit = 0;
    //If scores are given ,just reading
    if (matchScoreFile.size()>0)
    {
      DEBUGP(2);
      readMatchScoreFile(matchScoreFile,idList,idpair2MatchScore);
    }
    else
    {	 
	  int i,c,u,d,j,w,idx3;  
	  int MAX_VEC_SIZE;
	  double vecSmallestVal;
	  long double z,q; 
	  char matchScoreFile[200];
	  char RawMatrix[1000] ;
	  double couter;
	  double sum;
  	  VecBool finalVecMembership;
	  Hash_Int2Dbl finalVec;  
	  std::vector<matrix> RawH;
	  std::vector<matrix> RawAnswer;
      #pragma omp parallel num_threads(threadid) private(RawH,RawAnswer,sum,threadnumber,round,RawMatrix,count,finalVecMembership,finalVec,A,B,H,spanswer,fanswer,MAX_VEC_SIZE,vecSmallestVal,i,c,idx3,u,d,j,w,z,q,sp1_idxListre,sp2_idxListre,couter,iter)
       {
         for(round=0;round<work.size();round+=threadid)
	      {		     
		     threadnumber = omp_get_thread_num();
             count = threadnumber+round;
             if(count>=work.size())
		     {
				std::cout<<"Idle Thread :"<<threadnumber<<std::endl;
			 }
		     else
		     {
			  std::cout<<"Computing Task ...No. "<<count<<std::endl;
			  couter=0.0;
	          sum=0.0;
	          A.set_size(work[count].sp1nodes,work[count].sp1nodes);
	          B.set_size(work[count].sp2nodes,work[count].sp2nodes);  
	          for ( i=0; i < work[count].sp1nodes; i++)
                     {
	                    for (u=0; u < work[count].sp1_adjList[i].size(); u++)
	                      {                 
		                   c = work[count].sp1_adjList[i][u];
		                   z = (1.0/work[count].sp1_adjList[i].size());
		                   A(c,i)=z;
	                      }   	  
                     }
             std::cout<<"Setting A Matix "<< " element "<<A.n_nonzero<<" Summation "<<accu(A)<<std::endl;            
             for (j=0; j < work[count].sp2nodes; j++)
                     {
	                     for (w=0; w < work[count].sp2_adjList[j].size(); w++)
	                      {                     
		                    d = work[count].sp2_adjList[j][w];
		                    q = (1.0/work[count].sp2_adjList[j].size());
		                    B(d,j)=q;
	                      }    
                     }
              std::cout<<"Setting B Matix "<< " element "<<B.n_nonzero<<" Summation "<<accu(B)<<std::endl;         
	          std::cout<<"ThreadNo."<<threadnumber<<"--";DEBUGP(201);
              clock_t t1, t2;
			  t1 = clock();
              for(i=0;i<idList.size();i++)
                  {
		             sp1_idxListre[i]=0;
		             sp2_idxListre[i]=0;
                  }
              for(i=0;i<work[count].sp1_idxList.size();i++)
                  {
		             sp1_idxListre[work[count].sp1_idxList[i]]=i+1;
                  }
              for(i=0;i<work[count].sp2_idxList.size();i++)
                  {
		            sp2_idxListre[work[count].sp2_idxList[i]]=i+1;
                  }
	          H.set_size(work[count].sp2nodes,work[count].sp1nodes);
			  // construction and log of H matrix
			  sprintf(RawMatrix,"Temp/H0_%d.txt",count);
			  RawH.clear();
              for(i=0;i<theblast.size();i++)
                  {	  
	                   if(sp1_idxListre[theblast[i].no1]>= 1 && sp2_idxListre[theblast[i].no2]>=1)
	                   {
			              if (fabs(theblast[i].blast) > EPS) // larger than 
	                      {         
                           matrix tempmatrix;
                           H((sp2_idxListre[theblast[i].no2])-1,(sp1_idxListre[theblast[i].no1])-1)=theblast[i].blast;
                           //finn3<<sp2_idxListre[theblast[i].no2]-1<<" "<<sp1_idxListre[theblast[i].no1]-1<<" "<<  theblast[i].blast<<endl;       	 
                           tempmatrix.row =(sp2_idxListre[theblast[i].no2])-1;
                           tempmatrix.col =(sp1_idxListre[theblast[i].no1])-1;
                           tempmatrix.value=theblast[i].blast;
                           RawH.push_back(tempmatrix);
                          }
                          else{}         
		               }
                       //if(i%100000==0){std::cout<<"Thread : "<<threadnumber<<" i :"<< i <<std::endl;}	                   
                  } 
              std::cout<<"ThreadNo."<<threadnumber<<"--";DEBUGP(202);
              writematrix(RawH,RawMatrix);
              std::cout<<"ThreadNo."<<threadnumber<<"--";DEBUGP(203);
              //fgetc(stdin);
			  std::cout<<"Setting H0 Matix "<<" H element "<<H.n_nonzero<<std::endl;    
              t2 = clock();   
              couter=accu(H);
              H=(H/couter);  
              printf("Time of  Building  H : %lf\n", (t2-t1)/(double)(CLOCKS_PER_SEC));   
              spanswer.set_size(work[count].sp2nodes,work[count].sp1nodes);
              fanswer.set_size(work[count].sp2nodes,work[count].sp1nodes);
              spanswer=H; //Initial condition
              std::cout<<" NSD (Tensor Product) : "<<std::endl;
			  for (iter=0;iter<K;iter++)
              {
	                  t1 = clock();
	                  std::cout<<" Iteration : "<<iter+1<<std::endl; 
                      spanswer=((alpha)*(B*spanswer*(A.t())))+((1.0-alpha)*H);
                      std::cout<<" Amount of Nonzero :"<<spanswer.n_nonzero<<std::endl;
                      t2 = clock();
                      printf("Time of Tensor Product in this Iteration : %lf-----Task %d \n", (t2-t1)/(double)(CLOCKS_PER_SEC),count);
              }
              fanswer=spanswer;     
              std::cout<<"ThreadNo."<<threadnumber<<"--";DEBUGP(204);
              finalVecMembership.resize(work[count].sp1nodes*work[count].sp2nodes);
              //finalVec.resize(work[count].sp1nodes*work[count].sp2nodes);
              finalVecMembership = VecBool((work[count].sp1nodes*work[count].sp2nodes), false);
              int ttttt=0;
              RawAnswer.clear();
              for (i=0; i <work[count].sp1nodes; i++)
              {
	                  for (j=0; j <work[count].sp2nodes; j++)
	                    {
		                  idx3 = (i*work[count].sp2nodes)+j;
	                      //if (idx3 % 10000000 == 0) {std::cerr << "Thread : "<<threadnumber<<" i :"<< idx3 << std::endl;}
		                  if (fanswer(j,i)>EPS)
		                   {
			                 matrix tempanswer;
			                 tempanswer.row=j;
			                 tempanswer.col=i;
			                 tempanswer.value=fanswer(j,i);
			                 RawAnswer.push_back(tempanswer);
			                 finalVecMembership[idx3] = true;
			                 finalVec[idx3]=fanswer(j,i);	
			                 ttttt++;	 
		                   }
	                    }
              }
             std::cout<<"ThreadNo."<<threadnumber<<"--";DEBUGP(205); 
			 sprintf(RawMatrix,"Temp/Last_%d.txt",count);
			 writematrix(RawAnswer,RawMatrix);
			 std::cout<<"ThreadNo."<<threadnumber<<"--";DEBUGP(206);
             if (doPruning)
	         {
	           MAX_VEC_SIZE = prunedVecSize;
	           vecSmallestVal =pruneVecTopK(finalVec, finalVecMembership,MAX_VEC_SIZE);
	         }
	         std::cout<<"ThreadNo."<<threadnumber<<"--";DEBUGP(207);
             normalizeR_byL1norm(finalVec,finalVecMembership);  
	         std::cout<<"ThreadNo."<<threadnumber<<"--";DEBUGP(208);
	         readMatchScores(idList, finalVecMembership, finalVec,work[count].sp1_idxList,work[count].sp2_idxList,idpair2MatchScore);                                                          
             std::cout<<"ThreadNo."<<threadnumber<<"--";DEBUGP(209);
             //Reset
             finalVecMembership.clear();
			 finalVec.clear();
			 std::cout<<"Total score number : "<<idpair2MatchScore.size()<<std::endl;
		  }    
          #pragma omp barrier
          }
          
       }
       sprintf(matchScoreFile,"%s_match-score.txt",filePrefix.c_str());
       writeOutScores(idpair2MatchScore, matchScoreFile);   
     
   
     }

     
     // fgetc(stdin);
  //do k-partite matching
  Str2Int id2clstr;
  // Preprocess to find the ten score
  
  
  if(reverseflag == true)
  {
   std::cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Run Kpartite reverse!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
   doKpartiteMatchingStartWithStar_re(idpair2MatchScore, idList, id2idx, 
	 	      id2sp, spList, 
		      min_primary_fraction, 
		      min_secondary_fraction, 
		      max_per_species,
		      id2clstr, starflag, outputClusterFile);
  }
  else
  {
   std::cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Run Kpartite non-reverse!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
   doKpartiteMatchingStartWithStar_nore(idpair2MatchScore, idList, id2idx, 
	 	      id2sp, spList, 
		      min_primary_fraction, 
		      min_secondary_fraction, 
		      max_per_species,threadid,
		      id2clstr, starflag, outputClusterFile);	
  }
  
  
  
  DEBUGP(4);

}
