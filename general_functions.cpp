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

int check_parameters(std::ifstream * confile,int * clusters_num,int * functions_num,int * hash_tables_num){

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
    }else{
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

void k_unique_rand(int n,int m,std::vector<Vector*> & clusters,std::vector<Vector*> & all_vectors){

  std::random_device rd;
  std::default_random_engine generator(rd());
  std::uniform_int_distribution<int> distribution(0,2*n);
  int pos[m];
  int im = 0;

  for (int in = 0; in < n && im < m; ++in) {
    int rn = n - in;
    int rm = m - im;
    if (distribution(generator) % rn < rm)
      pos[im++] = in ;
  }

  for(int i=0;i<m;i++){
    clusters.insert (clusters.end(), all_vectors.at(pos[i]));
  }
}

void k_means_plus_plus(int n,int m,std::vector<Vector*> & clusters,std::vector<Vector*> & all_vectors,std::string & metric){

  Distance **dist= new Distance*[1];
  dist[0] = new Euclidean();

  std::random_device rd;
  std::default_random_engine generator(rd());
  std::uniform_int_distribution<int> distribution_int(0,n);
  int pos[m];

  long double** distance = new long double*[m-1];
  for(int i = 0; i < m-1; ++i)
    distance[i] = new long double[n];

  long double* average_distance = new long double[n];

  pos[0]=distribution_int(generator);
  clusters.insert (clusters.end(), all_vectors.at(pos[0]));

  for(int j=1;j<m;j++){
    long double all_distance=0;

    for(int i=0;i<all_vectors.size();i++){
     average_distance[i]=0;
      int fl=0;
      for(int d=0;d<j;d++){
        if(pos[d]==i){
          fl=1;
          break;
        }
      }

      distance[j-1][i]=(long double)dist[0]->dist(*clusters.at(j-1),*all_vectors.at(i));

      for(int d=0;d<j;d++)
        average_distance[i]=average_distance[i] +distance[d][i];

      average_distance[i]=average_distance[i]/(long double)j;
      all_distance=all_distance+average_distance[i];
  //    std::cout<<std::setprecision(100)<<(long double)distance[j-1][i]<<std::endl;
    }

    std::uniform_real_distribution<long double> distribution_real(0.0,(long double)all_distance);
    long double prob=distribution_real(generator);

    int s;
    for(s=0;s<all_vectors.size();s++){
      prob=prob-average_distance[s];
      if(prob<0.0)
        break;
    }
    pos[j]=s;
    clusters.insert (clusters.end(), all_vectors.at(s));

  }

  for(int i = 0; i < m-1; ++i)
    delete [] distance[i];
  delete [] distance;
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
