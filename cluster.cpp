#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>      // std::setprecision
#include <algorithm>
#include <chrono>

#include "check_info.h"
#include "param.h"
#include "Vector.h"
#include "general_functions.h"
#include "Cluster_Group.h"

int main (int argc, char** argv){

  std::string metric;
  std::ifstream infile;
  std::ifstream confile;
  std::ofstream outfile;

  int clusters_num=-1;
  int functions_num=NUMBER_OF_HASH_FUNCTIONS;
  int hash_tables_num=NUMBER_OF_HASH_TABLES;
  int Ms=NUMBER_OF_VECTORS;
  int probes=PROBES;

  if(check_info_cluster(argc,argv,&infile,&confile,&outfile,&metric)<0)//check parameters
    exit(-1);
    //std::cout << std::setprecision(data.length()) << s << '\n';

  std::transform(metric.begin(), metric.end(), metric.begin(), ::tolower); //set all letters to lowercase
  if(metric!="euclidean"&&metric!="cosine"){
    print_err(-6);
    return -1;
  }

  int check=0;// if -1 free memory and exit

  //get everything from the input file and store it in a vector then check number of dimensions
  std::vector<Vector*> all_vectors;
  read_info(&infile,&all_vectors);
  //all vectors must have the same dimensions
  if(check_dimensions(all_vectors)==-1) check=-1;
  //get everything from the cluster file anc check it
  if(check_parameters(&confile ,&clusters_num ,&functions_num ,&hash_tables_num,&Ms,&probes )==-1) check=-1;
  //check if we have more clusters than nesessary(we can change the maximum at param.h)
  if(check_cluster_num(all_vectors,clusters_num)==-1) check=-1;
  if(check==-1){//if we found sth wrong we stop
    for (std::vector<Vector*>::iterator it=all_vectors.begin(); it<all_vectors.end(); it++)
      delete(it[0]);
    return -1;
  }
  double atime=0;
  Cluster_Group *Clusters;

  Clusters= new Cluster_Group(clusters_num,metric);
  outfile<<"Algorithm: Ι1A1U1 "<<std::endl;
  outfile<<"Metric:  "<<metric<<std::endl;
  auto start = std::chrono::high_resolution_clock::now();
  Clusters->k_unique_rand_init(all_vectors);

  for(int i=0;i<NUM_TIMES;i++){
    Clusters->Lloyd_assignment(all_vectors);
    if(Clusters->k_means_update()==0)
      break;
  }
//  Clusters->Lloyd_assignment(all_vectors);
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);//
  atime=(double)duration.count()/1000000;
  Clusters->print_info(&outfile,atime);
  delete Clusters;
///////////////////////////////////////////
  Clusters= new Cluster_Group(clusters_num,metric);
  outfile<<"Algorithm: Ι1A1U2 "<<std::endl;
  outfile<<"Metric:  "<<metric<<std::endl;
  start = std::chrono::high_resolution_clock::now();
  Clusters->k_unique_rand_init(all_vectors);

  for(int i=0;i<NUM_TIMES;i++){

  Clusters->Lloyd_assignment(all_vectors);
  Clusters->PAM_improvement_update();
  }
  //  Clusters->Lloyd_assignment(all_vectors);
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);//
  atime=(double)duration.count()/1000000;
  Clusters->print_info(&outfile,atime);
  delete Clusters;
//////////////////////////////////////////
Clusters= new Cluster_Group(clusters_num,metric);
outfile<<"Algorithm: Ι1A2U1 "<<std::endl;
outfile<<"Metric:  "<<metric<<std::endl;
start = std::chrono::high_resolution_clock::now();
Clusters->k_unique_rand_init(all_vectors);

for(int i=0;i<NUM_TIMES;i++){

  Clusters->Range_search_assignment_LSH(all_vectors,2,functions_num);
  Clusters->k_means_update();
}
//   Clusters->Range_search_assignment_LSH(all_vectors,2,functions_num);
stop = std::chrono::high_resolution_clock::now();
duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);//
atime=(double)duration.count()/1000000;
Clusters->print_info(&outfile,atime);
delete Clusters;
//////////////////////////////////////////
Clusters= new Cluster_Group(clusters_num,metric);
outfile<<"Algorithm: Ι1A2U2 "<<std::endl;
outfile<<"Metric:  "<<metric<<std::endl;
start = std::chrono::high_resolution_clock::now();
Clusters->k_unique_rand_init(all_vectors);

for(int i=0;i<NUM_TIMES;i++){

  Clusters->Range_search_assignment_LSH(all_vectors,2,functions_num);
  Clusters->PAM_improvement_update();
}
//   Clusters->Range_search_assignment_LSH(all_vectors,2,functions_num);
stop = std::chrono::high_resolution_clock::now();
duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);//
atime=(double)duration.count()/1000000;
Clusters->print_info(&outfile,atime);
delete Clusters;
//////////////////////////////////////////
Clusters= new Cluster_Group(clusters_num,metric);
outfile<<"Algorithm: Ι1A3U1 "<<std::endl;
outfile<<"Metric:  "<<metric<<std::endl;
start = std::chrono::high_resolution_clock::now();
Clusters->k_unique_rand_init(all_vectors);

for(int i=0;i<NUM_TIMES;i++){

  Clusters->Range_search_assignment_Hypercube(all_vectors,functions_num,Ms,probes);
  Clusters->k_means_update();
}
//   Clusters->Range_search_assignment_Hypercube(all_vectors,functions_num,Ms,probes);
stop = std::chrono::high_resolution_clock::now();
duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);//
atime=(double)duration.count()/1000000;
Clusters->print_info(&outfile,atime);
delete Clusters;
//////////////////////////////////////////

Clusters= new Cluster_Group(clusters_num,metric);
outfile<<"Algorithm: Ι1A3U2 "<<std::endl;
outfile<<"Metric:  "<<metric<<std::endl;
start = std::chrono::high_resolution_clock::now();
Clusters->k_unique_rand_init(all_vectors);

for(int i=0;i<NUM_TIMES;i++){

  Clusters->Range_search_assignment_Hypercube(all_vectors,functions_num,Ms,probes);
  Clusters->PAM_improvement_update();
}
//   Clusters->Range_search_assignment_Hypercube(all_vectors,functions_num,Ms,probes);
stop = std::chrono::high_resolution_clock::now();
duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);//
atime=(double)duration.count()/1000000;
Clusters->print_info(&outfile,atime);
delete Clusters;
//////////////////////////////////////////

Clusters= new Cluster_Group(clusters_num,metric);
outfile<<"Algorithm: Ι2A1U1 "<<std::endl;
outfile<<"Metric:  "<<metric<<std::endl;
start = std::chrono::high_resolution_clock::now();
Clusters->k_means_plus_plus_init(all_vectors);

for(int i=0;i<NUM_TIMES;i++){

  Clusters->Lloyd_assignment(all_vectors);
  Clusters->k_means_update();
}
//   Clusters->Lloyd_assignment(all_vectors);
stop = std::chrono::high_resolution_clock::now();
duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);//
atime=(double)duration.count()/1000000;
Clusters->print_info(&outfile,atime);
delete Clusters;
//////////////////////////////////////////


Clusters= new Cluster_Group(clusters_num,metric);
outfile<<"Algorithm: Ι2A1U2 "<<std::endl;
outfile<<"Metric:  "<<metric<<std::endl;
start = std::chrono::high_resolution_clock::now();
Clusters->k_means_plus_plus_init(all_vectors);

for(int i=0;i<NUM_TIMES;i++){

  Clusters->Lloyd_assignment(all_vectors);
  Clusters->PAM_improvement_update();
}
//   Clusters->Lloyd_assignment(all_vectors);
stop = std::chrono::high_resolution_clock::now();
duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);//
atime=(double)duration.count()/1000000;
Clusters->print_info(&outfile,atime);
delete Clusters;
//////////////////////////////////////////

Clusters= new Cluster_Group(clusters_num,metric);
outfile<<"Algorithm: Ι2A2U1 "<<std::endl;
outfile<<"Metric:  "<<metric<<std::endl;
start = std::chrono::high_resolution_clock::now();
Clusters->k_means_plus_plus_init(all_vectors);

for(int i=0;i<NUM_TIMES;i++){

  Clusters->Range_search_assignment_LSH(all_vectors,2,functions_num);
  Clusters->k_means_update();
}
//     Clusters->Range_search_assignment_LSH(all_vectors,2,functions_num);
stop = std::chrono::high_resolution_clock::now();
duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);//
atime=(double)duration.count()/1000000;
Clusters->print_info(&outfile,atime);
delete Clusters;
//////////////////////////////////////////

Clusters= new Cluster_Group(clusters_num,metric);
outfile<<"Algorithm: Ι2A2U2 "<<std::endl;
outfile<<"Metric:  "<<metric<<std::endl;
start = std::chrono::high_resolution_clock::now();
Clusters->k_means_plus_plus_init(all_vectors);

for(int i=0;i<NUM_TIMES;i++){

  Clusters->Range_search_assignment_LSH(all_vectors,2,functions_num);
  Clusters->PAM_improvement_update();
}
//     Clusters->Range_search_assignment_LSH(all_vectors,2,functions_num);
stop = std::chrono::high_resolution_clock::now();
duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);//
atime=(double)duration.count()/1000000;
Clusters->print_info(&outfile,atime);
delete Clusters;
//////////////////////////////////////////

Clusters= new Cluster_Group(clusters_num,metric);
outfile<<"Algorithm: Ι2A3U1 "<<std::endl;
outfile<<"Metric:  "<<metric<<std::endl;
start = std::chrono::high_resolution_clock::now();
Clusters->k_means_plus_plus_init(all_vectors);

for(int i=0;i<NUM_TIMES;i++){

  Clusters->Range_search_assignment_Hypercube(all_vectors,functions_num,Ms,probes);
  Clusters->k_means_update();
}
//   Clusters->Range_search_assignment_Hypercube(all_vectors,functions_num,Ms,probes);
stop = std::chrono::high_resolution_clock::now();
duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);//
atime=(double)duration.count()/1000000;
Clusters->print_info(&outfile,atime);
delete Clusters;
//////////////////////////////////////////

Clusters= new Cluster_Group(clusters_num,metric);
outfile<<"Algorithm: Ι2A3U2 "<<std::endl;
outfile<<"Metric:  "<<metric<<std::endl;
start = std::chrono::high_resolution_clock::now();
Clusters->k_means_plus_plus_init(all_vectors);

for(int i=0;i<NUM_TIMES;i++){

  Clusters->Range_search_assignment_Hypercube(all_vectors,functions_num,Ms,probes);
  Clusters->PAM_improvement_update();
}
//   Clusters->Range_search_assignment_Hypercube(all_vectors,functions_num,Ms,probes);
stop = std::chrono::high_resolution_clock::now();
duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);//
atime=(double)duration.count()/1000000;
Clusters->print_info(&outfile,atime);
delete Clusters;
//////////////////////////////////////////


  for (std::vector<Vector*>::iterator it=all_vectors.begin(); it<all_vectors.end(); it++)
   delete it[0];

  infile.close();
  confile.close();
  outfile.close();
}
