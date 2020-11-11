#include <sys/stat.h>
#include "shou_algorithm.h"
#include <iostream>
#include <time.h>
#include <omp.h>
using namespace std;
/*
VecStr spList2;
bool sp_is_bigger_than(std::string sp1, std::string sp2){
  int rank1, rank2;
  for (int i=0;i<spList2.size();i++){
    if (spList2[i]==sp1)
      rank1 = i;
    if (spList2[i]==sp2)
      rank2 = i;
  }
  return (rank1 < rank2);
}
*/
struct answer
{
	double matchscore;
	int subindex;
};
struct Shouworklist
{
 std::string WorkId;
 std::string WorkId1; 
 vector <string> Worktempstring;
};
bool struct_cmp_by_freq(answer a, answer b)
{
    return a.matchscore > b.matchscore;
}
struct leaderlist
{
    std::vector<answer> allcandidate;
    int key; 
    answer min;
    leaderlist();
};
leaderlist::leaderlist()
{
 min.matchscore=1000000.0;	
 min.subindex=-1;
};
struct leaderboard
{
std::string maingroup;
std::string subgroup; // dmela or scere 
std::vector<leaderlist> dataset;
};
struct spNet
{
 int node_num;
 std::string spName;
};
void FindTheTopGammaScore(int shouthreadid,const Str2VecStr sp2idList ,const StrPair2Dbl idpair2MatchScore,const VecStr & spList,const Str2Int & id2idx,const Str2Str & id2sp,const int *CorrespondingTable,vector<leaderboard> & datawarehouse)

{
  double max;  
  //std::cout<<"come in "<<spList.size()<<std::endl;
  
  std::vector<string> tempstring;
  std::vector<Shouworklist> shouworklists;
  // Preprocessing Table of paralleling
  for (Str2VecStr::const_iterator Id=sp2idList.begin(); Id!=sp2idList.end();Id++)
  {
	 for (Str2VecStr::const_iterator Id1=sp2idList.begin(); Id1!=sp2idList.end();Id1++)
	 { 
       
       if(Id==Id1)
       {}
       else
       {
	    Shouworklist tempShouworklist;
	    tempShouworklist.WorkId=Id->first;
	    tempShouworklist.WorkId1=Id1->first;
	    tempShouworklist.Worktempstring=Id->second;
	    //std::cout<<Id->first<<" "<<Id1->first<<std::endl;
	    shouworklists.push_back(tempShouworklist);	   		   
	   }
     }
  }
  
  // Paralleling  
    int shouthreadnumber;
    //shouthreadid=4; //Default 4 thread
    int count,round,iter ; // Record which round and which one
    std::vector<leaderlist> alllist;
    StrPair2Dbl::const_iterator IdIt;
    Str2Int::const_iterator it_a; // First total number in idList
    Str2Int::const_iterator it_b; // Second total number in idList   
    Str2Str::const_iterator it_c;  //speciecs 1 
    Str2Str::const_iterator it_d; //speciecs 2  
    answer tempanswer;
    answer tempanswer1;
    leaderboard templeaderboard;
    int i;
     //std::cout<<"Thread"<<shouthreadid<<std::endl;
     //fgetc(stdin);
     std::cout<<"work"<<shouworklists.size()<<std::endl;
     std::cout<<"worklist"<<shouthreadid<<std::endl;
      std::cout<<"Finding Gamma Score...."<<std::endl;
      #pragma omp parallel num_threads(shouthreadid) private(tempanswer,tempanswer1,shouthreadnumber,round,count,alllist,it_a,it_b,it_c,it_d,i,templeaderboard,IdIt)
       {
         for(round=0;round<shouworklists.size();round+=shouthreadid)
	      {
		     std::cout<<"Round"<<(round/shouthreadid)<<std::endl;
		     shouthreadnumber = omp_get_thread_num();
             count = shouthreadnumber+round;  //Kernel of this for Loop
             //std::cout<<"Count"<<(count+1)<<std::endl;
			 if(count>=shouworklists.size())
		     {
			 std::cout<<" Idle thread"<<shouthreadnumber<<std::endl;
			 }
			 else
			 {				 
		         //vector<leaderlist> alllist;
		         alllist.clear();
		         alllist.resize(shouworklists[count].Worktempstring.size());                                 
                 //DEBUGP(1);
                 for(IdIt= idpair2MatchScore.begin(); IdIt!= idpair2MatchScore.end(); IdIt++)
                 {
		            it_a=id2idx.find(IdIt->first.first); // First total number in idList
                    it_b=id2idx.find(IdIt->first.second); // Second total number in idList   
                    it_c=id2sp.find(IdIt->first.first);  //speciecs 1 
                    it_d=id2sp.find(IdIt->first.second); //speciecs 2  			 
				    if (it_c->second==shouworklists[count].WorkId && it_d->second==shouworklists[count].WorkId1)  //Situation 1
	                     {
                          answer tempanswer;
                          tempanswer.matchscore=(IdIt->second -((CorrespondingTable[it_b->second])*(1e-20)));
	                      tempanswer.subindex=CorrespondingTable[it_b->second];
		                  alllist[CorrespondingTable[it_a->second]].allcandidate.push_back(tempanswer);
		                  alllist[CorrespondingTable[it_a->second]].key=CorrespondingTable[it_a->second];
		        		 }
		        	if (it_c->second==shouworklists[count].WorkId1 && it_d->second==shouworklists[count].WorkId)  //Situation 2
	                     {
	                      answer tempanswer1;
	                      tempanswer1.matchscore=(IdIt->second-((CorrespondingTable[it_a->second])*(1e-20)));
	                      tempanswer1.subindex=CorrespondingTable[it_a->second];
		                  alllist[CorrespondingTable[it_b->second]].allcandidate.push_back(tempanswer1);
		                  alllist[CorrespondingTable[it_b->second]].key=CorrespondingTable[it_b->second];
		                 }
		        			
					  
			     }
			     std::cout<<"TaskNo."<<count<<"--";DEBUGP(301);
		         for(i=0;i<alllist.size();i++)
                 {
			       if(alllist[i].allcandidate.size()>0)
			       { 
			       sort(alllist[i].allcandidate.begin(),alllist[i].allcandidate.end(),struct_cmp_by_freq);
			       }		
		         }
		         //DEBUGP(3);	
		       leaderboard templeaderboard;	 	           
               templeaderboard.dataset=alllist;
               templeaderboard.maingroup=shouworklists[count].WorkId;
               templeaderboard.subgroup=shouworklists[count].WorkId1;
               //std::cout<<templeaderboard.maingroup<<" "<<templeaderboard.subgroup<<" "<< shouthreadnumber<<std::endl;
               datawarehouse.at(count)=templeaderboard;		 
			     //DEBUGP(4);
			 }
	        
		  }
	   
	   #pragma omp barrier	
	   } 
 
DEBUGP(302);
}





















bool UDgreater (spNet sp1, spNet sp2)
{
 return sp1.node_num > sp2.node_num;
}


void doKpartiteMatchingStartWithStar_re(const StrPair2Dbl & idpair2MatchScore, 
            const VecStr & idList,
            const Str2Int & id2idx,
            const Str2Str & id2sp,
            const VecStr & spList,
            double min_primary_fraction,
            double min_secondary_fraction,
            int max_per_species,
            Str2Int & id2clstr, bool starflag, std::string outputClusterFile)
{
  int number_of_species = spList.size(); //@should equal to MAX_SPECIES
  std::string firstSpecies = FIRST_SPECIES;
  Str2VecStr sp2idList;
  VecStr empty_vector;
  std::cerr << "idSize:"<< idList.size() << std::endl;
  for(VecStr::const_iterator IdIt=idList.begin(); IdIt!=idList.end(); IdIt++)
  {
    std::string the_species = id2sp.find(*IdIt)->second;
    if (sp2idList.find(the_species)==sp2idList.end())
    {
		sp2idList[the_species] = empty_vector;
    }  
    sp2idList[the_species].push_back(*IdIt);
  }
  std::cerr << "sp2idList :"<< sp2idList.size() << std::endl;
  //1.36GB
  std::cerr << "Protein list for species generated" << std::endl;
  std::cerr << "sp2idList :"<< sp2idList.size() << std::endl;
  //@if we know the number of species here, we don't need a vector
   

  std::vector<bool> occupied(idList.size(),false);
  std::vector<int> ranks(idList.size(),GAMMA); //records the best rank of a protein
  //std::vector<VecStr> reverse_candidates(idList.size()); //reverse candidates;
  std::vector<std::vector<std::pair<std::string, double> > > reverse_candidates(idList.size());
  //initialized
  for (int i=0;i<idList.size();i++){
    for (int j=0;j<GAMMA;j++){
      std::pair<std::string, double> empty_pair("",0);
      reverse_candidates[i].push_back(empty_pair);
    }
  }
  VecStr spList2;
  //Loop on all the species (match 22)
  //char* sps[] = {"hsapi","mmusc","dmela","celeg","ecoli"};
  //for (int i=0;i<5;i++){
//    std::string the_species = sps[i];
//    std::string the_species = spList[i];
//    spList2.push_back(the_species);
//  }
/////////////////////////////////////////edit by BrainMa//////////////////////////////////////////////////////////
std::vector<spNet> netseq;
  spNet tmp;
  int spListIndex=0;
  for (int i=0;i<number_of_species;i++){
    tmp.node_num = sp2idList[spList[i]].size();
    tmp.spName = spList[i];
    netseq.push_back(tmp);
  }
  std::sort(netseq.begin(), netseq.end(), UDgreater);
  for (std::vector <spNet>::iterator spNetIt = netseq.begin();spNetIt != netseq.end();spNetIt++){
    spList2.push_back((*spNetIt).spName);
    //std::cerr <<"seq_order=>>spName: "<<(*spNetIt).spName<<"\tNode number: "<<(*spNetIt).node_num<<"\n";
    std::cerr <<"spList_order=>>spName: "<<spList[spListIndex]<<"\n";
    std::cerr <<"spList2_order=>>spName: "<<spList2[spListIndex]<<"\tNode number: "<<(*spNetIt).node_num<<"\n";
    spListIndex++;
  }

///////////////////////////////////////////////////////////////////////////////////////////////////
  std::ofstream OUTPUT (outputClusterFile.c_str());
  if(!OUTPUT.is_open())
  {
   std::cerr<< "Error opening file: "<<outputClusterFile<<"\n";
   assert(0);	
  }
  for (VecStr::const_iterator spIt=spList2.begin();spIt!=spList2.end();spIt++){
  std::string the_species = *spIt;
  std::string prefix = "stars_";
  std::cerr << the_species << " as center..." << std::endl;
  std::vector<Star> stars(0); //should be empty, but I am not so sure

  std::string strFilename = prefix + the_species + ".txt";
  struct stat stFileInfo;
  int intStat = -1;
  if(starflag == true)
  {
   intStat = stat(strFilename.c_str(), &stFileInfo);
  }
  //intStat == 0 means the file exists and don't run this
  if (intStat!=0 || !READ_CACHED_STARS){
  int star_loaded = 0;
  std::cerr<< "The star file,"<<strFilename<<" is not exist!!\nExport new star file: "<<strFilename<<"\n";
  for(VecStr::const_iterator firstNodeIt=sp2idList[the_species].begin();
      firstNodeIt!=sp2idList[the_species].end(); firstNodeIt++){
    star_loaded++;
   // std::cerr << "Prepare to load the " << star_loaded << " star." << std::endl;

    std::string top_protein = *firstNodeIt; double star_weight=0;
    //Star_Branches the_branches;
    std::string** the_branches = new std::string* [number_of_species];
    for (int sp=0;sp<number_of_species-1;sp++)
      the_branches[sp] = new std::string [GAMMA];
    
    //std::cerr << "Star " << star_loaded <<
    //  ": End allocating member for branches" << std::endl;
    
    bool first_species_has_passed = false;
    std::vector<double> max_weights; double max_weight=0;
    for (int sp=0;sp<number_of_species-1;sp++)
      {  //Species Loop
    
    std::string this_species = first_species_has_passed ? spList[sp+1] :
                                                          spList[sp];
    //std::cerr << "sp "<< sp << " : " << this_species << std::endl;
    if (this_species == the_species){
      first_species_has_passed = true;
      sp--;
      continue; 
    }
    double maxs[GAMMA] = {0.0}; //std::string max_ids[GAMMA];
    int max_nums[GAMMA] = {-1,-1,-1}; //@@macro
    //I thought about using std::list<double>
    //but dynamic allocation might turn out unsuitable/slower
    int num_of_protein_in_this_species = sp2idList[this_species].size();
    for (int protein_num = 0; protein_num<num_of_protein_in_this_species;
         protein_num++) 
       {
      //if (protein_num % 100 == 0)
      //  std::cerr << "100 proteins" << std::endl;
      //std::string branch_protein = *branchIt;
      std::string branch_protein = sp2idList[this_species][protein_num];
      StrPair2Dbl::const_iterator it_a = first_species_has_passed ?
        idpair2MatchScore.find(StrPair(top_protein, branch_protein)):
        idpair2MatchScore.find(StrPair(branch_protein, top_protein));
      //I guess this is very slow
      double score = it_a->second;
      std::cout << top_protein << " " << branch_protein << " " << score << std::endl;
      for (int i=0;i<GAMMA;i++)
      {
        if (score > maxs[i]){
          for (int j=GAMMA-1;j>i;j--)
          {
        maxs[j] = maxs[j-1];
        max_nums[j] = max_nums[j-1];
        //max_ids[j] = max_ids[j-1];
          }
          maxs[i] = score;
          max_nums[i] = protein_num;
          //max_ids[i] = branch_protein;
          break;
        }
      }
    }
    //std::cerr << "debug: " << maxs[0] << maxs[1] 
    //    << maxs[2] << maxs[3] << maxs[4] << std::endl;
    //std::cerr << "Prepare to load star.branches for sp "<< sp << std::endl;
    //star_weight += maxs[0];
    double this_max = maxs[0];
    max_weights.push_back(this_max);
    if (this_max > max_weight) max_weight = this_max;
    for (int candidate=0;candidate<GAMMA;candidate++){
      if (max_nums[candidate] < 0) break;
      if (candidate!=0 && 
          maxs[candidate] < BETTA * maxs[candidate-1]) {
        std::cerr << "Candidate gap threshold met at candidate "
              << candidate << " sp:" << sp << std::endl;
        break;
      }
      //std::cerr << "debug: " << max_nums[candidate] << std::endl;
      std::string candidate_protein = 
        sp2idList[this_species][max_nums[candidate]];
      //DEBUGP(92);
      //the_branches[sp][candidate] = max_ids[candidate];
      the_branches[sp][candidate] = candidate_protein;
      //DEBUGP(95);
      //int id = id2idx.find(max_ids[candidate])->second;
      int id = id2idx.find(candidate_protein)->second;
      //DEBUGP(98);
      if (ranks[id] > candidate)
        ranks[id] = candidate;
      //building reverse candidates

      //@@ filtering might have problem
      //reverse_candidates[id].push_back(top_protein);
      
      for (int i=0;i<GAMMA;i++){
        std::vector<std::pair<std::string, double> > list =
          reverse_candidates[id];
        if (maxs[candidate] > list[i].second){
          for (int j=GAMMA-1;j>i;j--)
        list[j] = list[j-1];
          list[i].first = top_protein;
          list[i].second = maxs[candidate];
          list.pop_back();
          //no more than five
          break;
        }
      }
    }// End of candidate building       
      } // End of Species Loop
    //Star complete_star = {star_weight, top_protein, the_branches};

    //Filtering out bad branches <Filter>
    int accepted = 0; double sum_of_weights = 0;
    for (int sp=0;sp<number_of_species-1;sp++){
      //if (max_weights[sp] < BETTA * max_weight){
      //    std::cerr << "filtering a branch of sp "<< sp << std::endl;
      //    for (int candidate=0;candidate<GAMMA;candidate++)
      //      the_branches[sp][candidate] = "";
      //}else{
    sum_of_weights+= max_weights[sp];
    accepted++;
    //}
    }
    //</Filter>
    Star complete_star;    
    complete_star.weight = sum_of_weights / accepted;
    complete_star.top_protein = top_protein;
    complete_star.branches = the_branches;
    stars.push_back(complete_star);
  }

  std::cerr << "Sorting stars..." << std::endl;
  std::sort(stars.begin(), stars.end(), sort_star); //sorted into decending order
  //not working without sort_star, Error: ....discards quality
  std::cerr << "Stars sorted. There are "
        << stars.size() << " stars." << std::endl;
  std::cerr << "The biggest weight is " << stars[0].weight 
        << ". The middle weight is " << stars[stars.size()/2].weight
        << std::endl;
  
  export_stars_and_ranks(stars, number_of_species-1, ranks, strFilename);
  } else {
    std::cerr << "Loading cached stars (sorted) from " 
          << strFilename << std::endl;
    import_stars_and_ranks(stars, number_of_species-1, ranks, strFilename);

    //Rebuilding the reverse_candidate from exported stars with no weight
    //(match 13)
    std::cerr << "Rebuilding reverse_candidates" << std::endl;
    for(std::vector<Star>::const_iterator starIt=stars.begin();
    starIt!=stars.end();starIt++){
      Star this_star = *starIt;
      for (int sp=0;sp<number_of_species-1;sp++){
    for (int j=0;j<GAMMA;j++){
      std::string the_protein = this_star.branches[sp][j];
      int the_id = id2idx.find(the_protein)->second;
      double the_weight = getScore(the_protein, this_star.top_protein,
                       idpair2MatchScore);

      if (the_weight !=0)
        //std::cerr << "The following should be successful with " 
        //        << "the protein: " << the_protein << std::endl;
      //Made a ridiculous fault here, it has to be a reference
        // std::vector<std::pair<std::string, double> >* the_list =
        // reverse_candidates[the_id];    
      for (int re=0;re<GAMMA;re++){
        //std::cerr << "Debugging, " 
        //        << the_id << " " << the_list[re].second << std::endl;
        //mm16475
        if (the_weight > reverse_candidates[the_id][re].second){
          //std::cerr << "The modification of reverse_candidates occurred at" 
          //        << re << std::endl;
          for (int re2=GAMMA-1;re2>re;re2--)
        reverse_candidates[the_id][re2] = reverse_candidates[the_id][re2-1];
          reverse_candidates[the_id][re].first = this_star.top_protein;
          reverse_candidates[the_id][re].second = the_weight;
          //std::cerr << "Debugging, " << the_list[re].second << std::endl;
          //the_list.pop_back();
          break;
        }
      }
    }
      } 
    }
  }
  
  std::cerr << "Transfering reverse_candidates to reverse_candidates2"
        << std::endl;
  std::vector<VecStr> reverse_candidates2(idList.size());
  for (int i=0;i<idList.size();i++){
    for (int j=0;j<GAMMA;j++)
      reverse_candidates2[i].push_back(reverse_candidates[i][j].first);
  }
  
    
  int star_exported = 0;
  //Examing the stars and possibly substitute with candidates
  //@@Issues: 1. always have a human protein in a group
  //          2. not global in sense of Graemin's definition
  std::vector<VecStr> result_stars(idList.size()); 
  std::vector<VecStr> reverse_candidates_for_stars(idList.size());
  for(std::vector<Star>::const_iterator starIt=stars.begin();
      starIt!=stars.end();starIt++){
    Star this_star = *starIt;
    //the star top is occupied by many-to-many maching
    if (occupied[id2idx.find(this_star.top_protein)->second]) continue;
    Star last_star = (starIt==stars.begin()) ? *starIt : *(starIt-1);
    //neglect matching with low weight
    //if (this_star.weight < SUBSET_THRESHOLD) break;
    if (this_star.weight < BETTA * last_star.weight) break;
    VecStr this_matching = this_star_is_ok_re(this_star,
                       occupied,
                       ranks,
                       idpair2MatchScore,
                       id2idx,
                       number_of_species-1,
                       reverse_candidates2);
    //std::cerr << "Returned vector has " << this_matching.size()
    //<< " elements" << std::endl;

    bool space_met = false;
    int star_id = id2idx.find(this_star.top_protein)->second;
    for(VecStr::const_iterator id=this_matching.begin();
    id!=this_matching.end();id++){
      if ((*id).size() == 0) {
    space_met = true;
    continue;
      }
      if (!space_met){
    result_stars[star_id].push_back(*id);
    occupied[id2idx.find(*id)->second] = true;
      }else{
    reverse_candidates_for_stars[star_id].push_back(*id);
      }
      //std::cout << *id << " ";
      //set occupied<> for result matching (match 10 modified)
      //occupied[id2idx.find(*id)->second] = true;  
      //std::cout << "Test vecotr iterator!!!" << std::endl;
    }
    //for (int i=0;i<this_matching.size();i++)
    //  std::cerr << this_matching[i] << " ";
    
    //if (this_matching.size()!=0){
    //  std::cout << std::endl;
    star_exported++;
    std::cerr << "The " << star_exported << "th star prepared to export" 
          << std::endl;
    //}
    
  }//end of result star production
  std::cerr << "Final step: regrouping" << std::endl;
  for (int match=0;match<result_stars.size();match++){
    VecStr this_match = result_stars[match];
    if (this_match.size()==0) continue;
    for(VecStr::const_iterator id=this_match.begin();
    id!=this_match.end();id++){
      std::cout << *id << " ";
      OUTPUT << *id << " ";
    }
    VecStr candidates = reverse_candidates_for_stars[match];
    bool found_it = false;
    for (VecStr::const_iterator candidIt=candidates.begin();
     candidIt!=candidates.end();candidIt++){
      int star_id = id2idx.find(*candidIt)->second;
      if (star_id <= match) continue; // no < to avoid duplications
      VecStr anti_candidates = reverse_candidates_for_stars[star_id];
      for (VecStr::const_iterator anti_candidIt=anti_candidates.begin();
       anti_candidIt!=anti_candidates.end();anti_candidIt++){
    if (*anti_candidIt == this_match.back()){
      found_it = true;
      VecStr this_match = result_stars[star_id];
      for(VecStr::const_iterator id=this_match.begin();
          id!=this_match.end();id++){
        std::cout << *id << " ";
        OUTPUT << *id << " ";
      }
      result_stars[star_id].clear(); //side effect
      break;
    }
      }
      if (found_it) break;
    }
    std::cout << std::endl;
    OUTPUT << std::endl;
  }
  }//end of species loop
  OUTPUT.close();
}

void doKpartiteMatchingStartWithStar_nore(const StrPair2Dbl & idpair2MatchScore, 
		    const VecStr & idList,
			const Str2Int & id2idx,
			const Str2Str & id2sp,
			const VecStr & spList,
			double min_primary_fraction,
			double min_secondary_fraction,
			int max_per_species,
			int shouthreadid,
			Str2Int & id2clstr, bool starflag, std::string outputClusterFile)
{
  int number_of_species = spList.size(); //@should equal to MAX_SPECIES
  //std::string firstSpecies = FIRST_SPECIES;
  Str2VecStr sp2idList;
  std::cerr << idList.size() << std::endl;
  int numberofspecies=spList.size();
  int CorrespondingTable[idList.size()];
  VecStr empty_vector;
  clock_t starttime, finishtime;
  double  duration ;   
  int book=0;
  int page=0;
  for(VecStr::const_iterator IdIt=idList.begin(); IdIt!=idList.end(); IdIt++)
  {
     
      std::string the_species = id2sp.find(*IdIt)->second;
      //cout<<"id2sp.find(*IdIt)->first: " <<id2sp.find(*IdIt)->first <<"\n id2sp.find(*IdIt)->second: "  <<id2sp.find(*IdIt)->second<<"\n";
    if (sp2idList.find(the_species)==sp2idList.end())
      {
		  sp2idList[the_species] = empty_vector;
		  page=0;
	  }
	  CorrespondingTable[book]=page;  
	  book=book+1;
	  page=page+1; 
      sp2idList[the_species].push_back(*IdIt);
  
  }
  //1.36GB
  std::cerr << "Protein list for species generated" << std::endl;

  //@if we know the number of species here, we don't need a vector
  

  std::vector<bool> occupied(idList.size(),false);
  std::vector<int> ranks(idList.size(),GAMMA); //records the best rank of a protein
  //std::vector<VecStr> reverse_candidates(idList.size()); //reverse candidates;
  std::vector<std::vector<std::pair<std::string, double> > > reverse_candidates(idList.size());
  //initialized
  for (int i=0;i<idList.size();i++){
    for (int j=0;j<GAMMA;j++){
      std::pair<std::string, double> empty_pair("",0);
      reverse_candidates[i].push_back(empty_pair);
    }
  }

  //Loop on all the species (match 22)
  VecStr spList2;
  //char* sps[] = {"hsapi","mmusc","dmela","celeg","scere"};
  //char* sps[] = {"scere","celeg","dmela","mmusc","hsapi"};
  //for (int i=0;i<number_of_species;i++){
  //  std::string the_species = sps[i];
  //  spList2.push_back(the_species);
  //}
////////////////////////////////////////////edit by BrainMa///////////////////////////////////////////////////////////////////
  std::vector<spNet> netseq;
  spNet tmp;
  int spListIndex=0;
  for (int i=0;i<number_of_species;i++)
  {
    tmp.node_num = sp2idList[spList[i]].size();
    tmp.spName = spList[i];
    netseq.push_back(tmp);
  }
  //First sort 
  starttime=clock();
  std::sort(netseq.begin(), netseq.end(), UDgreater);
  finishtime=clock();
  duration = (double)(finishtime - starttime) / CLOCKS_PER_SEC;   
  printf( "First sort %f seconds\n", duration);
  for (std::vector <spNet>::iterator spNetIt = netseq.begin();spNetIt != netseq.end();spNetIt++)
  {
    spList2.push_back((*spNetIt).spName);
    //std::cerr <<"seq_order=>>spName: "<<(*spNetIt).spName<<"\tNode number: "<<(*spNetIt).node_num<<"\n";
    std::cerr <<"spList_order=>>spName: "<<spList[spListIndex]<<"\n";
    std::cerr <<"spList2_order=>>spName: "<<spList2[spListIndex]<<"\tNode number: "<<(*spNetIt).node_num<<"\n";
    spListIndex++;
  }
  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  std::ofstream OUTPUT (outputClusterFile.c_str());
  if(!OUTPUT.is_open())
  {
   std::cerr<< "Error opening file: "<<outputClusterFile<<"\n";
   assert(0);	
  }
  std::vector <leaderboard> datawarehouse; // Given the size is better 
  int reallysize=numberofspecies*(numberofspecies-1);
  datawarehouse.resize(reallysize);
  FindTheTopGammaScore(shouthreadid,sp2idList,idpair2MatchScore,spList,id2idx,id2sp,CorrespondingTable,datawarehouse);
  //std::cout<<datawarehouse.size()<<std::endl;
  for (VecStr::const_iterator spIt=spList2.begin();spIt!=spList2.end();spIt++) // Every species do 
  {
     std::string the_species = *spIt;
     std::string prefix = "stars_";
     std::cerr << the_species << " as center..." << std::endl;
     std::vector<Star> stars(0); //should be empty, but I am not so sure
     std::string strFilename = prefix + the_species + ".txt";
     struct stat stFileInfo;
     int intStat = -1;
     if(starflag == true)
     {
       intStat = stat(strFilename.c_str(), &stFileInfo);
     }  //intStat == 0 means the file exists and don't run this
     if (intStat !=0 || !READ_CACHED_STARS)
     {
      //std::cerr<< "The star file,"<<strFilename<<" is not exist!!\nExport new star file: "<<strFilename<<"\n";
      int star_loaded = 0;
      int iterator=0;
      for( VecStr::const_iterator firstNodeIt=sp2idList[the_species].begin();firstNodeIt!=sp2idList[the_species].end(); firstNodeIt++) //
      {
           star_loaded++;
           
           //std::cerr << "Prepare to load the " << star_loaded << " star." << std::endl;
           std::string top_protein = *firstNodeIt; //First name 
           double star_weight=0;
           //Star_Branches the_branches;
           std::string** the_branches = new std::string* [number_of_species];
           for (int sp=0;sp<number_of_species-1;sp++)
           {
			 the_branches[sp] = new std::string [GAMMA];
           }    
           bool first_species_has_passed = false;
           std::vector<double> max_weights; 
           double max_weight=0;
           for (int sp=0;sp <number_of_species-1;sp++)
           {  	
	         std::string this_species = first_species_has_passed ? spList[sp+1] :spList[sp];
	         //std::cerr << "sp "<< sp << " : " << this_species << std::endl;
	         //std::cerr <<this_species<< "  "<< the_species<<" "<<top_protein<< std::endl;
	         //fgetc(stdin);
             if (this_species == the_species)
	         {
	            first_species_has_passed = true;
	            sp--;
	            continue; 
	         }
	         
             int rightkey=0;
             //std::cout<<datawarehouse.size()<<std::endl;
	         for(int j=0;j<datawarehouse.size();j++)
	         {
			 if(the_species==datawarehouse[j].maingroup && this_species==datawarehouse[j].subgroup)
			    {
				  //std::cerr <<this_species<< "  "<< the_species<< std::endl;
				  rightkey=j;
				  //std::cout<<rightkey<<std::endl;
                  break;
				}
			 }
			 //DEBUGP(2);
	         double maxs[GAMMA] = {0.0}; //std::string max_ids[GAMMA];
	         int max_nums[GAMMA] = {-1,-1,-1}; //@@macro
	         if(datawarehouse[rightkey].dataset[iterator].allcandidate.size()<GAMMA)
	         {
				 int remaind=0;
				 for(int j=0;j<datawarehouse[rightkey].dataset[iterator].allcandidate.size();j++)
	             {
				   maxs[j]=datawarehouse[rightkey].dataset[iterator].allcandidate[j].matchscore+((datawarehouse[rightkey].dataset[iterator].allcandidate[j].subindex)*(1e-20));
				   max_nums[j]=datawarehouse[rightkey].dataset[iterator].allcandidate[j].subindex;
			     }
			     // DEBUGP(3);
			      for(int k=datawarehouse[rightkey].dataset[iterator].allcandidate.size();k<GAMMA;k++)
	             {
				    for(int m=0;m<datawarehouse[rightkey].dataset[iterator].allcandidate.size();m++)
				    {
					   if(remaind==max_nums[m])
					   {
						   remaind++;
					   }
					}
					maxs[k]=6.953e-310;
				    max_nums[k]=remaind;
				    remaind++;
			     } 
			 
			 }
	         else
	         {
	             //DEBUGP(40);
	             for(int j=0;j<GAMMA;j++)
	             {
				    maxs[j]=datawarehouse[rightkey].dataset[iterator].allcandidate[j].matchscore+((datawarehouse[rightkey].dataset[iterator].allcandidate[j].subindex)*(1e-20));
				    max_nums[j]=datawarehouse[rightkey].dataset[iterator].allcandidate[j].subindex;
			     }
		     }
	        //std::cerr << "debug: " << maxs1[0] <<maxs1[2]<< maxs1[4] << maxs1[6] << maxs1[8] << std::endl;
	        //std::cerr << "debug: " <<  max_nums1[0] <<" "<< max_nums1[2]<<" "<< max_nums1[4]<<" "<< max_nums1[6] <<" "<< max_nums1[8] << std::endl;
	       double maxs1[GAMMA] = {0.0}; //std::string max_ids[GAMMA];
	       int max_nums1[GAMMA] = {-1,-1,-1}; //@@macro
	         //I thought about using std::list<double>
	         //but dynamic allocation might turn out unsuitable/slower
	       int num_of_protein_in_this_species = sp2idList[this_species].size();
	         //Find corresponding
	         
	       /*
	         for (int protein_num = 0; protein_num<num_of_protein_in_this_species;protein_num++) 
	          {
	           std::string branch_protein = sp2idList[this_species][protein_num];
	           StrPair2Dbl::const_iterator it_a = first_species_has_passed ?idpair2MatchScore.find(StrPair(top_protein, branch_protein)):
	           idpair2MatchScore.find(StrPair(branch_protein, top_protein));
	           //I guess this is very slow
	           double score = it_a->second;
	           for (int i=0;i<GAMMA;i++) // Find the top gamma
	            {
	              if (score > maxs1[i])
	              {
	                for (int j=GAMMA-1;j>i;j--)
	                {
		              maxs1[j] = maxs1[j-1];
		              max_nums1[j] = max_nums1[j-1];
	                }
	              maxs1[i] = score;
	              max_nums1[i] = protein_num;
	              break;
	               }
	            }
	          }*/
	// test
	/*
	for(int i=0;i<GAMMA;i++)
	{
		if(max_nums1[i]!=max_nums[i])
		{
		  std::cerr << "debug: " << max_nums[0]<<" " << max_nums[1]<<" "<< max_nums[2]<<" " << max_nums[3]<<" " << max_nums[4]<<" " << max_nums[5]<<" " << max_nums[6]<<" "<< max_nums[7] <<" "<<max_nums[8] <<" "<< max_nums[9]<< std::endl;
	      std::cerr << "debug: " << max_nums1[0] <<" "<< max_nums1[1]<<" "<< max_nums1[2]<<" " << max_nums1[3]<<" " << max_nums1[4]<<" " << max_nums1[5]<<" " << max_nums1[6]<<" "<< max_nums1[7]<<" " << max_nums1[8] <<" "<< max_nums1[9]<< std::endl;
		  std::cout << i << std::endl;
		  fgetc(stdin);
		  break;
		}
	}*/
	
	

	//std::cerr << "debug: " << sp2idList[this_species][max_nums[0]] << sp2idList[this_species][max_nums[1]]<< sp2idList[this_species][max_nums[2]] << sp2idList[this_species][max_nums[3]] <<sp2idList[this_species][max_nums[4]] << std::endl;
	//std::cerr << "Prepare to load star.branches for sp "<< sp << std::endl;
	//star_weight += maxs[0];
	
	
	
	
	
	
	double this_max = maxs[0];
	max_weights.push_back(this_max);
	if (this_max > max_weight) max_weight = this_max;
	for (int candidate=0;candidate<GAMMA;candidate++){
	  if (max_nums[candidate] < 0) break;
	  if (candidate!=0 && 
	      maxs[candidate] < BETTA * maxs[candidate-1]) 
	     {
	   //std::cerr << "Candidate gap threshold met at candidate "<< candidate << " sp:" << sp << std::endl;
	    break;
	     }
	  //std::cerr << "debug: " << max_nums[candidate] << std::endl;
	  std::string candidate_protein = 
	    sp2idList[this_species][max_nums[candidate]];
	  //DEBUGP(92);
	  //the_branches[sp][candidate] = max_ids[candidate];
	  the_branches[sp][candidate] = candidate_protein;
	  //DEBUGP(95);
	  //int id = id2idx.find(max_ids[candidate])->second;
	  int id = id2idx.find(candidate_protein)->second;
	  //DEBUGP(98);
	  if (ranks[id] > candidate)
	    ranks[id] = candidate;
	  //building reverse candidates

	  //@@ filtering might have problem
	  //reverse_candidates[id].push_back(top_protein);
	  
	  for (int i=0;i<GAMMA;i++){
	    std::vector<std::pair<std::string, double> > list =
	      reverse_candidates[id];
	    if (maxs[candidate] > list[i].second){
	      for (int j=GAMMA-1;j>i;j--)
		list[j] = list[j-1];
	      list[i].first = top_protein;
	      list[i].second = maxs[candidate];
	      list.pop_back();
	      //no more than five
	      break;
	    }
	  }
	}// End of candidate building	    
      } // End of Species Loop
    //Star complete_star = {star_weight, top_protein, the_branches};

    //Filtering out bad branches <Filter>
    int accepted = 0; double sum_of_weights = 0;
    for (int sp=0;sp<number_of_species-1;sp++){
      //if (max_weights[sp] < BETTA * max_weight){
      //	std::cerr << "filtering a branch of sp "<< sp << std::endl;
      //	for (int candidate=0;candidate<GAMMA;candidate++)
      //	  the_branches[sp][candidate] = "";
      //}else{
	sum_of_weights+= max_weights[sp];
	accepted++;
	//}
    }
    //</Filter>
    Star complete_star;    
    complete_star.weight = sum_of_weights / accepted;
    complete_star.top_protein = top_protein;
    complete_star.branches = the_branches;
    stars.push_back(complete_star);
    iterator=iterator+1;
  }

  std::cerr << "Sorting stars..." << std::endl;
  //Second sort
  starttime=clock();
  std::sort(stars.begin(), stars.end(), sort_star); //sorted into decending order
  finishtime=clock(); 
  duration = (double)(finishtime - starttime) / CLOCKS_PER_SEC;   
  printf( " Second sort : %f seconds\n", duration );
  
  
  
  
  
  //not working without sort_star, Error: ....discards quality
  std::cerr << "Stars sorted. There are "
	    << stars.size() << " stars." << std::endl;
  std::cerr << "The biggest weight is " << stars[0].weight 
	    << ". The middle weight is " << stars[stars.size()/2].weight
	    << std::endl;
  
  //export_stars_and_ranks(stars, number_of_species-1, ranks, strFilename);
  } else {
    std::cerr << "Loading cached stars (sorted) from " 
	      <<strFilename << std::endl;
    import_stars_and_ranks(stars, number_of_species-1, ranks, strFilename);
  }
  
  std::cerr << "Transfering reverse_candidates to reverse_candidates2"
	    << std::endl;
  std::vector<VecStr> reverse_candidates2(idList.size());
  for (int i=0;i<idList.size();i++){
    for (int j=0;j<GAMMA;j++)
      reverse_candidates2[i].push_back(reverse_candidates[i][j].first);
  }
  
    
  int star_exported = 0;
  //Examing the stars and possibly substitute with candidates
  //@@Issues: 1. always have a human protein in a group
  //          2. not global in sense of Graemin's definition
  

  for(std::vector<Star>::const_iterator starIt=stars.begin();
      starIt!=stars.end();starIt++)
      {
    Star this_star = *starIt;
    //the star top is occupied by many-to-many maching
    if (occupied[id2idx.find(this_star.top_protein)->second]) continue;
    Star last_star = (starIt==stars.begin()) ? *starIt : *(starIt-1);
    //neglect matching with low weight
    //if (this_star.weight < SUBSET_THRESHOLD) break;
    if (this_star.weight < BETTA * last_star.weight) break;
    VecStr this_matching = this_star_is_ok(this_star,
					   occupied,
					   ranks,
					   idpair2MatchScore,
					   id2idx,
					   number_of_species-1,
					   reverse_candidates2);
    //std::cerr << "Returned vector has " << this_matching.size()
    //<< " elements" << std::endl;

    for(VecStr::const_iterator id=this_matching.begin();
	id!=this_matching.end();id++){
      //std::cout<< *id << " ";
      OUTPUT << *id << " ";
      //set occupied<> for result matching (match 10 modified)
      occupied[id2idx.find(*id)->second] = true;  
      //std::cout << "Test vecotr iterator!!!" << std::endl;
    }
    
    //for (int i=0;i<this_matching.size();i++)
    //  std::cerr << this_matching[i] << " ";
    if (this_matching.size()!=0){
      //std::cout << std::endl;
      OUTPUT << std::endl;
      star_exported++;
     // std::cerr << "The " << star_exported << "th star exported" << std::endl;
    }
  }
  }//end of species loop
  OUTPUT.close();

}

bool sort_star (const Star& star1, const Star& star2){
  return (star1.weight > star2.weight);
}

double getScore (const std::string & protein1, const std::string & protein2,/*,
         int sp1, int sp2*/
         const StrPair2Dbl & idpair2MatchScore){
  //static int first_species_number =
  double candidate1 = idpair2MatchScore.
    find(StrPair(protein1, protein2))->second;
  double candidate2 = idpair2MatchScore.
    find(StrPair(protein2, protein1))->second;
  return (candidate1 > candidate2) ? candidate1 : candidate2;
}

struct Edge
{
  int sp1;
  int candid1;
  int sp2;
  int candid2;
  double score;
};




bool sort_edge (const Edge& edge1, const Edge& edge2){
  return (edge1.score > edge2.score);
}
//////////////////////////random walk ? /////////////////////////////////////////////
VecStr this_star_is_ok(Star & star,
		       const std::vector<bool> & occupied,
		       const std::vector<int> & ranks,
		       const StrPair2Dbl & idpair2MatchScore,
		       const Str2Int & id2idx,
		       int nsp,
		       const std::vector<VecStr> & reverse_candidates){
		   
  //std::cerr << "Begin this_star_is_ok" << std::endl;
  //for every 2 subset of the branches
  //try_all 1 substitution, if fail, start dropping things
  
  //same for every call
  /*
  std::vector<Examine_Subset> subsets;
  //@marcro here?
  for (int i=0;i<nsp;i++){
    for (int j=i+1;j<nsp;j++){
      Examine_Subset subset;
      subset.push_back(i);subset.push_back(j);
      //Examine_Subset subset = {i,j};
      subsets.push_back(subset);
    }
  }
  */

  bool bypass_first_step = false; int no_empty_sp = 0; int last_noempty_sp = -1;
  for (int i=0;i<nsp;i++) {
    if (star.branches[i][0].size()!=0) {
      no_empty_sp++;
      last_noempty_sp = i;
    }
  }
  if (no_empty_sp < 2) bypass_first_step = true;
  if (last_noempty_sp == -1){
    VecStr empty;
   // std::cerr << "There's no branch at all in this star" << std::endl;
    return empty;
  }

  double the_smallest_weight_among_leaves;
  VecStr result_branch_proteins;
  if (!bypass_first_step) {
 // std::cerr << "Starting Part 1 of many to many" << std::endl;
  std::vector<std::string> branch_proteins;
  for (int i=0;i<nsp;i++){
    for (int j=0;j<GAMMA;j++){
      std::string this_protein = star.branches[i][j];
      std::string empty = "";
      if (this_protein.size() == 0 ||
	  occupied[id2idx.find(this_protein)->second] ||
	  j > ranks[id2idx.find(this_protein)->second])
	branch_proteins.push_back(empty);
      else
	branch_proteins.push_back(this_protein);
    }
  }

  //<debug>
  /*
  for (VecStr::const_iterator leafIt=branch_proteins.begin()+1;
       leafIt!=branch_proteins.end();leafIt++)
    std::cerr << *leafIt << " ";
  std::cerr << std::endl;
  */
  //</debug>

  std::vector<Edge> edges;
  for (int i=0;i<branch_proteins.size();i++){
    if (branch_proteins[i].size() == 0) {
      //std::cerr << "Protein " << branch_proteins[i] 
      //<< " is skipped (occupied)" << std::endl;
      continue;
    }	     
    for (int j=i;j<branch_proteins.size();j++){
      if (branch_proteins[j].size() == 0) continue;      
      Edge edge;
      edge.sp1 = i / GAMMA; edge.candid1 = i % GAMMA;
      edge.sp2 = j / GAMMA; edge.candid2 = j % GAMMA;
      //Don't include those that are in the same species
      if (edge.sp1 == edge.sp2) continue;
      edge.score = getScore(branch_proteins[i], branch_proteins[j], idpair2MatchScore);
      edges.push_back(edge);
    }
  }

  //std::cerr << "There are " << edges.size() << " edges"<< std::endl;
  std::sort(edges.begin(), edges.end(), sort_edge);

  VecStr empty;
  if (edges.size() == 0 || edges[0].score == 0) {
    //std::cerr << "No edge or the score of the first edge is 0" << std::endl;
    return empty;
  }
  for (std::vector<Edge>::iterator edgeIt=edges.begin()+1;edgeIt!=edges.end();
       edgeIt++){
    Edge latter = *edgeIt; Edge former = *(edgeIt-1);
    if (latter.score < BETTA*former.score){
      /*std::cerr << "Threshold for Part1(among leaves) met"
		<< " the former score is " << former.score
		<< " the latter score is " << latter.score
		<< std::endl;*/
      the_smallest_weight_among_leaves = former.score;
      edges.erase(edgeIt,edges.end());
      break;
    }
  }

  //std::cerr << "There are " << edges.size() << " edges"<< " after erasing bad many-to-many" << std::endl;  
  int * external_degree = new int[branch_proteins.size()];
  for (int i=0;i<branch_proteins.size();i++) external_degree[i] = 0;

  //Vetex Weight (match 36)
  std::vector<double> vertex_weight(branch_proteins.size(), 0);
  //double * vertex_weight = new double[branch_proteins.size()];
  //for (int i=0;i<branch_proteins.size();i++) vertex_weight[i] = 0;

  for (std::vector<Edge>::iterator edgeIt=edges.begin();edgeIt!=edges.end();
       edgeIt++){
    Edge edge = *edgeIt;
    if (edge.sp1 == edge.sp2) continue;
    external_degree[edge.sp1*GAMMA + edge.candid1]++;
    //std::cerr << "Add degree to the protein " 
    //	      << branch_proteins[edge.sp1*GAMMA + edge.candid1]
    //	      << std::endl;    
    external_degree[edge.sp2*GAMMA + edge.candid2]++;

    vertex_weight[edge.sp1*GAMMA + edge.candid1] += edge.score;
    vertex_weight[edge.sp2*GAMMA + edge.candid2] += edge.score;
  }

  //(vertex weight)/(vertex degree) (match 39)
  for (int leaf=0;leaf<branch_proteins.size();leaf++)
    if (external_degree[leaf] != 0)
      vertex_weight[leaf] = vertex_weight[leaf] / external_degree[leaf];
  
  //Sort the vector to get a threshold
  double threshold = 0;
  std::vector<double> to_sort = vertex_weight;
  std::sort(to_sort.begin(), to_sort.end());
  threshold = *(to_sort.end()-1) / 10;
  /*
  for (std::vector<double>::iterator revIt=to_sort.end()-1;
       revIt!=to_sort.begin();revIt--){
    std::cerr << *revIt << " " << *(revIt-1) << std::endl;
    if (*(revIt-1) < BETTA*(*revIt))
      threshold = *(revIt-1);
  }
  */ //grep '^[0-9\.e-]* [0-9\.e-]*$' _condor_stderr | grep '^.[^ ]' 

  /* skipped (match 36)
  for (int leaf=0;leaf<branch_proteins.size();leaf++){
    int degree = external_degree[leaf];
    std::string this_protein = branch_proteins[leaf];
    if (this_protein.size()==0) continue;
    if (degree > 2 || (degree == 1 || degree == 2) && 
	ranks[id2idx.find(this_protein)->second] == (leaf % GAMMA))
      result_branch_proteins.push_back(this_protein);
    else if (degree == 0)
      std::cerr << "The isolated leaf: "<< this_protein << " is dropped"
		<< std::endl;
    else
      std::cerr << "The protein "<< this_protein << " is dropped"
		<< " because of rank" << std::endl;
  }
  */

  for (int leaf=0;leaf<branch_proteins.size();leaf++){
    std::string this_protein = branch_proteins[leaf];
    if (this_protein.size()!=0 && vertex_weight[leaf] > threshold)
      result_branch_proteins.push_back(this_protein);
  }
    
  /*  
  //Take the former subsequence and check if it is a complete graph. If not, drop the leaf if its "candidate rank" is not the best.
 
  bool complete = true;
  bool occupied[5*GAMMA][5*GAMMA];
  for (int i=0;i<edges.size();i++){
    Edge edge = edges[i];
    int num1 = edge.sp1 * GAMMA + edge.candid1;
    int num2 = edge.sp2 * GAMMA + edge.candid2;
    occupied[num1][num2] = true;
    occupied[num2][num1] = true;
  }

  int num_of_proteins = 0;
  for (int i=0;i<5*GAMMA;i++){
    int num_in_row = 0
    for (int j=i;j<5*GAMMA;j++)
      if (occupied[i][j]) num_in_row++;
    if (num_in_row){
      num_of_proteins = num_in_row;
      break;
    }
  }

  for (int i=0;i<5*GAMMA;i++){
    int num_in_row = 0
    for (int j=i;j<5*GAMMA;j++)
      if (occupied[i][j]) num_in_row++;
    if (num_in_row != 0 && num_in_row != num_of_proteins){
      complete = false;
      break;
    }
  }

  //VecStr result_star;
  std::set<std::string> result_star;
  //result_star.push_back(star.top_protein);
  
  if (complete){
    
  }else{
  }
  */  

  } else {
    //std::cerr << "There's only one branch" << std::endl;
    for(int j=0;j<GAMMA;j++){
      //last_noempty_sp could be 0 here
      std::string this_protein = star.branches[last_noempty_sp][j];
      if (this_protein.size()!=0 && 
	  j == ranks[id2idx.find(this_protein)->second] &&
	  !occupied[id2idx.find(this_protein)->second]){
	result_branch_proteins.push_back(star.branches[last_noempty_sp][j]);
	the_smallest_weight_among_leaves = getScore(star.branches[last_noempty_sp][j], star.top_protein, idpair2MatchScore);
      }
    }
  }

  if (result_branch_proteins.size() == 0){
    //std::cerr << "This is an isolated first species(human) protein "<< "which should not occur"<<std::endl;
    return result_branch_proteins;
  }

  
  /*
  std::cerr << "Starting Part 2 of many to many" << std::endl;
  VecStr common_candidates_for_first;
  VecStr start_seed = 
    reverse_candidates[id2idx.find(result_branch_proteins[0])->second];
  for (VecStr::const_iterator proIt=start_seed.begin();
       proIt!=start_seed.end();proIt++){
    if (result_branch_proteins.size() == 1){
      std::cerr << "Waring: there's only one leaf." << std::endl;
      common_candidates_for_first = start_seed;
      break;
    }
    std::string this_protein = *proIt;
    bool common_intersection = true;
    for (VecStr::const_iterator leafIt=result_branch_proteins.begin()+1;
	 leafIt!=result_branch_proteins.end();leafIt++){
      //DEBUGP(372);
      VecStr reverse_candidate = reverse_candidates
	[id2idx.find(*leafIt)->second];
      //DEBUGP(375);
      bool in_this_subset = false;
      for (VecStr::const_iterator canIt=reverse_candidate.begin();
	   canIt!=reverse_candidate.end();canIt++){
	if (*canIt == this_protein){
	  in_this_subset = true;
	  break;
	}
      }
      if (!in_this_subset){
	common_intersection = false;
	break;
      }
    }
    if (common_intersection && this_protein.size()!=0) 
      common_candidates_for_first.push_back(this_protein);
  }
  
  VecStr final_result = result_branch_proteins;
  for (VecStr::const_iterator firstIt=common_candidates_for_first.begin();
       firstIt!=common_candidates_for_first.end();
       firstIt++){
    bool include = true;
    for (VecStr::const_iterator leafIt=result_branch_proteins.begin()+1;
	 leafIt!=result_branch_proteins.end();leafIt++){
      if (getScore(*leafIt, *firstIt, idpair2MatchScore) 
	  < BETTA * the_smallest_weight_among_leaves){
	include = false;
	break;
      }
    }
    if (include) final_result.push_back(*firstIt);
  }
  */
  
  //One-to-many
  VecStr final_result = result_branch_proteins;
  final_result.push_back(star.top_protein);
  //One-to-many

  /*
  for (int sp=0;sp<nsp;sp++){
    std::string protein=star.branches[sp][0];
    if (protein.size()>0) {
      occupied[id2idx.find(protein)->second] = true;
      result_star.push_back(protein);
    }
  }
  */
  //std::cerr << "This result star has " << result_star.size() 
  //	    << " proteins" << std::endl;
  //return result_star;
  //std::cerr << "This many-to-many result star has " << final_result.size() << " proteins" << std::endl;
  return final_result;
}
//////////////////////////random walk ? end/////////////////////////////////////////////
VecStr this_star_is_ok_re(Star & star,
               const std::vector<bool> & occupied,
               const std::vector<int> & ranks,
               const StrPair2Dbl & idpair2MatchScore,
               const Str2Int & id2idx,
               int nsp,
               const std::vector<VecStr> & reverse_candidates){
  std::cerr << "Begin this_star_is_ok" << std::endl;
  //for every 2 subset of the branches
  //try_all 1 substitution, if fail, start dropping things
  
  //same for every call
  /*
  std::vector<Examine_Subset> subsets;
  //@marcro here?
  for (int i=0;i<nsp;i++){
    for (int j=i+1;j<nsp;j++){
      Examine_Subset subset;
      subset.push_back(i);subset.push_back(j);
      //Examine_Subset subset = {i,j};
      subsets.push_back(subset);
    }
  }
  */

  bool bypass_first_step = false; int no_empty_sp = 0; int last_noempty_sp = -1;
  for (int i=0;i<nsp;i++) {
    if (star.branches[i][0].size()!=0) {
      no_empty_sp++;
      last_noempty_sp = i;
    }
  }
  if (no_empty_sp < 2) bypass_first_step = true;
  if (last_noempty_sp == -1){
    VecStr empty;
    std::cerr << "There's no branch at all in this star" << std::endl;
    return empty;
  }

  double the_smallest_weight_among_leaves;
  VecStr result_branch_proteins;
  if (!bypass_first_step) {
  std::cerr << "Starting Part 1 of many to many" << std::endl;
  std::vector<std::string> branch_proteins;
  for (int i=0;i<nsp;i++){
    for (int j=0;j<GAMMA;j++){
      std::string this_protein = star.branches[i][j];
      std::string empty = "";
      if (this_protein.size() == 0 ||
      occupied[id2idx.find(this_protein)->second] ||
      j > ranks[id2idx.find(this_protein)->second])
    branch_proteins.push_back(empty);
      else
    branch_proteins.push_back(this_protein);
    }
  }

  //<debug>
  /*
  for (VecStr::const_iterator leafIt=branch_proteins.begin()+1;
       leafIt!=branch_proteins.end();leafIt++)
    std::cerr << *leafIt << " ";
  std::cerr << std::endl;
  */
  //</debug>

  std::vector<Edge> edges;
  for (int i=0;i<branch_proteins.size();i++){
    if (branch_proteins[i].size() == 0) {
      //std::cerr << "Protein " << branch_proteins[i] 
      //<< " is skipped (occupied)" << std::endl;
      continue;
    }        
    for (int j=i;j<branch_proteins.size();j++){
      if (branch_proteins[j].size() == 0) continue;      
      Edge edge;
      edge.sp1 = i / GAMMA; edge.candid1 = i % GAMMA;
      edge.sp2 = j / GAMMA; edge.candid2 = j % GAMMA;
      //Don't include those that are in the same species
      if (edge.sp1 == edge.sp2) continue;
      edge.score = getScore(branch_proteins[i], branch_proteins[j], idpair2MatchScore);
      edges.push_back(edge);
    }
  }

  std::cerr << "There are " << edges.size() << " edges"
        << std::endl;
  std::sort(edges.begin(), edges.end(), sort_edge);

  VecStr empty;
  if (edges.size() == 0 || edges[0].score == 0) {
    std::cerr << "No edge or the score of the first edge is 0" << std::endl;
    return empty;
  }
  for (std::vector<Edge>::iterator edgeIt=edges.begin()+1;edgeIt!=edges.end();
       edgeIt++){
    Edge latter = *edgeIt; Edge former = *(edgeIt-1);
    if (latter.score < BETTA*former.score){
      /*std::cerr << "Threshold for Part1(among leaves) met"
        << " the former score is " << former.score
        << " the latter score is " << latter.score
        << std::endl;*/
      the_smallest_weight_among_leaves = former.score;
      edges.erase(edgeIt,edges.end());
      break;
    }
  }

  std::cerr << "There are " << edges.size() << " edges"
        << " after erasing bad many-to-many" << std::endl;  
  int * external_degree = new int[branch_proteins.size()];
  for (int i=0;i<branch_proteins.size();i++) external_degree[i] = 0;

  for (std::vector<Edge>::iterator edgeIt=edges.begin();edgeIt!=edges.end();
       edgeIt++){
    Edge edge = *edgeIt;
    if (edge.sp1 == edge.sp2) continue;
    external_degree[edge.sp1*GAMMA + edge.candid1]++;
    //std::cerr << "Add degree to the protein " 
    //        << branch_proteins[edge.sp1*GAMMA + edge.candid1]
    //        << std::endl;    
    external_degree[edge.sp2*GAMMA + edge.candid2]++;
  }
  for (int leaf=0;leaf<branch_proteins.size();leaf++){
    int degree = external_degree[leaf];
    std::string this_protein = branch_proteins[leaf];
    if (this_protein.size()==0) continue;
    if (degree > 2 || (degree == 1 || degree == 2) && 
    ranks[id2idx.find(this_protein)->second] == (leaf % GAMMA))
      result_branch_proteins.push_back(this_protein);
    else if (degree == 0)
      std::cerr << "The isolated leaf: "<< this_protein << " is dropped"
        << std::endl;
    else
      std::cerr << "The protein "<< this_protein << " is dropped"
        << " because of rank" << std::endl;
  }
    
  /*  
  //Take the former subsequence and check if it is a complete graph. If not, drop the leaf if its "candidate rank" is not the best.
 
  bool complete = true;
  bool occupied[5*GAMMA][5*GAMMA];
  for (int i=0;i<edges.size();i++){
    Edge edge = edges[i];
    int num1 = edge.sp1 * GAMMA + edge.candid1;
    int num2 = edge.sp2 * GAMMA + edge.candid2;
    occupied[num1][num2] = true;
    occupied[num2][num1] = true;
  }

  int num_of_proteins = 0;
  for (int i=0;i<5*GAMMA;i++){
    int num_in_row = 0
    for (int j=i;j<5*GAMMA;j++)
      if (occupied[i][j]) num_in_row++;
    if (num_in_row){
      num_of_proteins = num_in_row;
      break;
    }
  }

  for (int i=0;i<5*GAMMA;i++){
    int num_in_row = 0
    for (int j=i;j<5*GAMMA;j++)
      if (occupied[i][j]) num_in_row++;
    if (num_in_row != 0 && num_in_row != num_of_proteins){
      complete = false;
      break;
    }
  }

  //VecStr result_star;
  std::set<std::string> result_star;
  //result_star.push_back(star.top_protein);
  
  if (complete){
    
  }else{
  }
  */  

  } else {
    std::cerr << "There's only one branch" << std::endl;
    for(int j=0;j<GAMMA;j++){
      //last_noempty_sp could be 0 here
      std::string this_protein = star.branches[last_noempty_sp][j];
      if (this_protein.size()!=0 && 
      j == ranks[id2idx.find(this_protein)->second] &&
      !occupied[id2idx.find(this_protein)->second]){
    result_branch_proteins.push_back(star.branches[last_noempty_sp][j]);
    the_smallest_weight_among_leaves = getScore(star.branches[last_noempty_sp][j], star.top_protein, idpair2MatchScore);
      }
    }

    //If this match is a 1 branch 1-1, than don't ouput it unless this is the
    //last chance (match 32.a) removed in (match 34)
    /*
    if (result_branch_proteins.size() == 1){
      std::cerr << "This is a 1-1 correpondence" << std::endl;
      if (sp_is_bigger_than(star.top_protein, result_branch_proteins[0])){
    VecStr empty;
    return empty;
      }
    }
    */
  }

  if (result_branch_proteins.size() == 0){
    std::cerr << "This is an isolated first species(human) protein "
          << "which should not occur"<<std::endl;
    return result_branch_proteins;
  }

  
  
  std::cerr << "Starting Part 2 of many to many" << std::endl;
  VecStr common_candidates_for_first;
  VecStr start_seed = 
    reverse_candidates[id2idx.find(result_branch_proteins[0])->second];
  for (VecStr::const_iterator proIt=start_seed.begin();
       proIt!=start_seed.end();proIt++){
    if (result_branch_proteins.size() == 1){
      std::cerr << "Waring: there's only one leaf." << std::endl;
      common_candidates_for_first = start_seed;
      break;
    }
    std::string this_protein = *proIt;
    bool common_intersection = true;
    for (VecStr::const_iterator leafIt=result_branch_proteins.begin()+1;
     leafIt!=result_branch_proteins.end();leafIt++){
      //DEBUGP(372);
      VecStr reverse_candidate = reverse_candidates
    [id2idx.find(*leafIt)->second];
      //DEBUGP(375);
      bool in_this_subset = false;
      for (VecStr::const_iterator canIt=reverse_candidate.begin();
       canIt!=reverse_candidate.end();canIt++){
    if (*canIt == this_protein){
      in_this_subset = true;
      break;
    }
      }
      if (!in_this_subset){
    common_intersection = false;
    break;
      }
    }
    if (common_intersection && this_protein.size()!=0) 
      common_candidates_for_first.push_back(this_protein);
  }
  
  VecStr final_result = result_branch_proteins;
  final_result.push_back(star.top_protein);
  //The return format is "ms1112" "hs1113" "" "hs2223" "hs3321" 
  //(match 26)
  std::string empty_string = "";
  final_result.push_back(empty_string);

  //If the common intersection does not include the center,
  //delete the common intersection (match 32)
  bool include_center = false;
  for (VecStr::const_iterator firstIt=common_candidates_for_first.begin();
       firstIt!=common_candidates_for_first.end();
       firstIt++)
    if (star.top_protein == *firstIt)
      include_center = true;
  
  if (include_center)
    for (VecStr::const_iterator firstIt=common_candidates_for_first.begin();
     firstIt!=common_candidates_for_first.end();
     firstIt++)
      final_result.push_back(*firstIt);
  
  
  
  /*
  for (VecStr::const_iterator firstIt=common_candidates_for_first.begin();
       firstIt!=common_candidates_for_first.end();
       firstIt++){
    bool include = true;
    for (VecStr::const_iterator leafIt=result_branch_proteins.begin()+1;
     leafIt!=result_branch_proteins.end();leafIt++){
      if (getScore(*leafIt, *firstIt, idpair2MatchScore) 
      < BETTA * the_smallest_weight_among_leaves){
    include = false;
    break;
      }
    }
    if (include) final_result.push_back(*firstIt);
  }
  */
  
  
  //One-to-many
  //VecStr final_result = result_branch_proteins;
  //final_result.push_back(star.top_protein);
  //One-to-many

  /*
  for (int sp=0;sp<nsp;sp++){
    std::string protein=star.branches[sp][0];
    if (protein.size()>0) {
      occupied[id2idx.find(protein)->second] = true;
      result_star.push_back(protein);
    }
  }
  */
  //std::cerr << "This result star has " << result_star.size() 
  //        << " proteins" << std::endl;
  //return result_star;
  std::cerr << "This many-to-may result star has " << final_result.size() 
        << " proteins" << std::endl;
  return final_result;
}


void export_stars_and_ranks (const std::vector<Star> & stars,
                 const int nsp,
                 const std::vector<int> & ranks,
                 const std::string filename){
  std::ofstream st(filename.c_str());
  //Exporting ranks
  if(!st.is_open())
  {
   std::cerr<< "Error opening file: "<<filename<<"\n";
   assert(0);	
  }
  for (std::vector<int>::const_iterator rankIt=ranks.begin();
       rankIt!=ranks.end();rankIt++)
    st << *rankIt << " ";
  st << std::endl;
    
  for(std::vector<Star>::const_iterator starIt=stars.begin();
      starIt!=stars.end();starIt++){
    Star this_star = *starIt;
    st << this_star.top_protein << " "
       << this_star.weight << std::endl;
    for (int sp=0;sp < nsp;sp++){
      for (int j=0;j < GAMMA;j++)
    st << this_star.branches[sp][j] << " ";
      st << std::endl;
    }
  }
  st.close();
}

void import_stars_and_ranks (std::vector<Star> & stars,
                 const int nsp,
                 std::vector<int> & ranks,
                 const std::string filename){
  std::ifstream st(filename.c_str());
  //char rank_num[1000000];
  //st.getline(rank_num, 1000000-1);
  std::cerr << "Loading ranks" << std::endl;
  if(st.fail())
  {
   std::cerr<< "The star file,"<<filename<<" is not exist~!!!\n";
   assert(0);
  }
  if(st.peek() == 32 || st.peek() == -1 || st.peek() == 9)
  {
   std::cerr<< "The star file,"<<filename<<" is empty~!!!\n";
   assert(0);	
  }
  for (std::vector<int>::iterator rankIt=ranks.begin();
       rankIt!=ranks.end();rankIt++)
    st >> (*rankIt);

  while (st.get() != '\n'); //this moves the "unformatted" pointer
                            //to the second line
  char line[200];
  int counter = 0;
  while (st.getline(line,200-1)){
    counter++;
    std::cerr << "Loading the " << counter << "th star." << std::endl;
    Star this_star;
    std::istringstream istrstr(line);
    istrstr >> this_star.top_protein >> this_star.weight;
    //std::cerr << "top protein and weight loaded" << std::endl;
    
    std::string** the_branches = new std::string* [nsp];
    for (int sp=0;sp<nsp;sp++)
      the_branches[sp] = new std::string [GAMMA];

    for (int sp=0;sp<nsp;sp++){
      st.getline(line,200-1);
      std::istringstream istrstr(line);
      for (int j=0;j<GAMMA;j++)
    istrstr >> the_branches[sp][j];
    }
    this_star.branches = the_branches;
    stars.push_back(this_star);
  }
  st.close();
}

void try_candidates(Star & star,
            const StrPair2Dbl & idpair2MatchScore,
            const std::vector<bool> & occupied,
            const std::vector<int> & ranks,
            const Str2Int & id2idx,
            int sp1, int sp2){
  std::cerr << "Trying candidates..." << std::endl;
  //examnine the candidates for sp1
  for (int i=1;i<GAMMA;i++){
    double s1 = getScore(star.top_protein,star.branches[sp1][i],
             idpair2MatchScore);
    double s2 = getScore(star.top_protein,star.branches[sp2][0],
             idpair2MatchScore);
    double s12 = getScore(star.branches[sp1][i],star.branches[sp2][0],
              idpair2MatchScore);
    //check whether the candidate is already taken
    int id = id2idx.find(star.branches[sp1][i])->second;
    if (occupied[id]) continue;
    if (ranks[id] < i) continue; 
       //the candidate protein has a higher rank in other protein
    if (s12 > BETTA*(s1+s2) / 2){
      //moves the candidate set [0,1,2,3,4] -> [i,i+1,"","",""]
      for (int i2=i;i2<GAMMA;i2++)
    star.branches[sp1][i2-i] = star.branches[sp1][i2];
      for (int i2=GAMMA-i;i2<GAMMA;i2++)
    star.branches[sp1][i2]="";
      return;
    }
  }
  
  //examine the candidates for sp2
  for (int j=1;j<GAMMA;j++){
    double s1 = getScore(star.top_protein,star.branches[sp1][0],
             idpair2MatchScore);
    double s2 = getScore(star.top_protein,star.branches[sp2][j],
             idpair2MatchScore);
    double s12 = getScore(star.branches[sp1][0],star.branches[sp2][j],
              idpair2MatchScore);
    if (occupied[id2idx.find(star.branches[sp2][j])->second]) continue;
    if (s12 > BETTA*(s1+s2) / 2){
      for (int j2=j;j2<GAMMA;j2++)
    star.branches[sp2][j2-j] = star.branches[sp1][j2];
      for (int j2=GAMMA-j;j2<GAMMA;j2++)
    star.branches[sp2][j2]="";
      return;
    }
  }
  //no candidate is satisfied, drop sp1
  for (int k=0;k<GAMMA;k++)
    star.branches[sp1][k]="";
  return;
} 
