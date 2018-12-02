#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>      // std::setprecision
#include <algorithm>

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


  Cluster_Group *Clusters= new Cluster_Group(clusters_num,metric);
  Clusters->k_means_plus_plus_init(all_vectors);
 //Clusters->Lloyd_assignment(all_vectors);
 // Clusters->Range_search_assignment_LSH(all_vectors,2,functions_num);
  Clusters->Range_search_assignment_Hypercube(all_vectors,functions_num,Ms,probes);

//  for(int i=0;i<clusters_num;i++)
//    std::cout<<Clusters->get_centroid_id_n(i)<<std::endl;

  for (std::vector<Vector*>::iterator it=all_vectors.begin(); it<all_vectors.end(); it++)
   delete(it[0]);

  delete Clusters;

}
/*
for(int i=0;i<clusters_num;i++){
  std::cout<<clusters2.at(i)->get_identity()<<std::endl;
  std::vector<coordinate> *coordinatesss=clusters2.at(i)->get_coordinates();
}*/
  //  for(int j=0;j<203;j++)
  //      std::cout<<coordinatesss->at(j)<<std::endl;
