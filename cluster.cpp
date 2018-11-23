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

int main (int argc, char** argv){

  std::string data;
  std::string metric;

  std::ifstream infile;
  std::ifstream confile;
  std::ofstream outfile;

  int clusters_num;
  int functions_num;
  int hash_tables_num;

  if(check_info_cluster(argc,argv,&infile,&confile,&outfile,&metric)<0)//check parameters
    exit(-1);
    //std::cout << std::setprecision(data.length()) << s << '\n';

  std::transform(metric.begin(), metric.end(), metric.begin(), ::tolower); //set all letters to lowercase

  if(metric!="euclidean"&&metric!="cosine"){
    print_err(-6);
    return -1;
  }
////////////////////////////////////////////
//get everything from the input file and store it in a vector

  coordinate s;
  std::vector<Vector*> all_vectors;
  char *token;

  while(std::getline( infile, data )){

    token = strtok(&data[0u], " ,");
    std::string id=token;

    Vector* V = new Vector (id);
    token = strtok(NULL, " ,");

     while( token != NULL ) {
       coordinate s;
       sscanf(token, "%Lf", &s);
       V->add_coordinate(s);
       token = strtok(NULL, " ,");
     }
     all_vectors.insert (all_vectors.end(), V);
   }
   infile.close();
///////////////////////////////////////
  
  std::vector<Vector*>::iterator V_it;
  V_it=all_vectors.begin();
  std::vector<coordinate> *coordinates=V_it[0]->get_coordinates();
  int dimensions=coordinates->size();
  int vectors=all_vectors.size();
  std::cout<<dimensions<<"  "<<vectors<<std::endl;






  for (std::vector<Vector*>::iterator it=all_vectors.begin(); it<all_vectors.end(); it++)
    delete(it[0]);
}
