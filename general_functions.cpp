#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <iomanip>      // std::setprecision

#include "Vector.h"
#include "general_functions.h"
#include "Distance.h"

void read_info(std::ifstream * infile,std::vector<Vector*> * all_vectors){

  std::string data;
  coordinate s;
  char *token;

  while(std::getline( *infile, data )){

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
     all_vectors->insert (all_vectors->end(), V);
   }
   infile->close();

}

int check_parameters(std::ifstream * confile,int * clusters_num,int * functions_num,int * hash_tables_num,int *Ms,int *probes){

  std::string data;

  *confile>>data;

  while(!confile->eof()){
    std::transform(data.begin(), data.end(), data.begin(), ::tolower);
    if(data.compare("number_of_clusters:")==0){
      *confile>>data;
      *clusters_num=std::stoi(data);
    }
    else if(data.compare("number_of_hash_functions:")==0){
      *confile>>data;
      *functions_num=std::stoi(data);
    }
    else if(data.compare("number_of_hash_tables:")==0){
      *confile>>data;
      *hash_tables_num=std::stoi(data);
    }
    else if(data.compare("number_of_probes:")==0){
      *confile>>data;
      *probes=std::stoi(data);
    }
    else if(data.compare("number_of_vectors:")==0){
      *confile>>data;
      *Ms=std::stoi(data);
    }
    else{
      *confile>>data;

      std::cerr <<"Conf file parameter not recognised.Continue?(Y,N)"<<")"<<std::endl;
      std::string answer;
      std::cin>>answer;
      std::transform(answer.begin(), answer.end(), answer.begin(), ::tolower);

      if(!answer.compare("n")){
        return -1;
      }
    }
    *confile>>data;
  }

  if(*clusters_num==-1){
    std::string answer;
    std::cout<<"Number_of_clusters parameter not found.Do you want to add one now?(Y,N)."<<std::endl;
    std::cin>>answer;
    std::transform(answer.begin(), answer.end(), answer.begin(), ::tolower);
    while(TRUE){
      if(!answer.compare("y")){
        std::cout<<"Value:";
        std::cin>>answer;
        *clusters_num=std::stoi(answer);
        break;
      }else if(!answer.compare("n")){
        std::cout<<"Number_of_clusters has taken a vlue of 3"<<std::endl;
        *clusters_num=3;
        break;
      }
      else{
        std::cout<<"Please enter your answer again  :"<<std::endl;
        std::cin>>answer;

      }
    }
  }
}

int check_dimensions(std::vector<Vector*> & all_vectors){

  std::vector<Vector*>::iterator it;
  it=all_vectors.begin();
  std::vector<coordinate> *coordinates=it[0]->get_coordinates();
  int dimensions=coordinates->size();
  for (std::vector<Vector*>::iterator it=all_vectors.begin(); it<all_vectors.end(); it++){
    std::vector<coordinate> *coordinates=it[0]->get_coordinates();
    if(dimensions!=coordinates->size()){
      return -1;
    }
  }
  return 0;
}

int check_cluster_num(std::vector<Vector*> & all_vectors,int clusters_num){

  int vectors=all_vectors.size();
  if(clusters_num>vectors/CLUSTERS_DIV){
    std::cerr<<"Number of clusters is way too big"<<std::endl;
     return -1;
  }
  return 0;
}

long double absolute(long double x){
  if((long double)x<0)
    return -x;
  return x;
}

long double power(long double x,int n){
  long double r=1.0;

  for(int i=0;i<n;i++)
    r=r*x;
  return r;
}

int mod (long long  int a, unsigned long long  int b)
{
   return (a%b+b)%b;
}
