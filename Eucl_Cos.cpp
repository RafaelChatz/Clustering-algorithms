#include <iostream>
#include <chrono>
#include <random>
#include <numeric>
#include <climits>
#include <iomanip>      // std::setprecision
#include <algorithm>

#include "Eucl_Cos.h"
#include "param.h"

extern long int bytes_allocated;

//Ls_Hashing class
Ls_Hashing::Ls_Hashing(){
  bytes_allocated+=sizeof(Ls_Hashing);

  //std::cout<<"Ls_Hashing created type :";
}

Ls_Hashing::~Ls_Hashing(){
  //std::cout<<"Ls_Hashing "<<std::endl;
}
////////////////////////////////////////////

//Euclidean class
Euclidean::Euclidean(int w,int dimensions){

  //std::cout<<"Euclidean "<<std::endl;
  bytes_allocated+=sizeof(Euclidean);

  int min=0;
  int max=INT_MAX;

  std::random_device rd;     // only used once to initialise (seed) engine

  std::normal_distribution<double> distribution_normal(0.0,1.0);
  std::uniform_real_distribution<long double> distribution_uniform(0.0,(double)w);

  t=distribution_uniform(rd);

  this->w=w;

  for(int j=0;j<dimensions;j++)
    vector_V.insert (vector_V.end(),distribution_normal(rd));

    bytes_allocated+=sizeof(Vector)+128*sizeof(coordinate);

}

Euclidean::~Euclidean(void){
  //std::cout<<"Euclidean Destoyed:";
}

int Euclidean::h_func(Vector &v){
    std::vector<coordinate> *p=v.get_coordinates();
    double p_v = std::inner_product(p->begin(), p->end(), vector_V.begin(), 0.0);

    return (int)floor((211*p_v + t)/(double)w);
}


//////////////////////////////////////////////

//Cosine class
Cosine::Cosine(int dimensions){
  bytes_allocated+=sizeof(Cosine);

  std::random_device rd;     // only used once to initialise (seed) engine
  srand (time(NULL));
  std::normal_distribution<long double> distribution_normal(0.0,1.0);
  for(int j=0;j<dimensions;j++)
    vector_R.insert (vector_R.end(),distribution_normal(rd));

    bytes_allocated+=sizeof(Vector)+128*sizeof(coordinate);

    //  std::cout<<"Cosine "<<std::endl;
}

Cosine::~Cosine(void){
//  std::cout<<"Cosine Destoyed:";
}

int Cosine::h_func(Vector &x){

  std::vector<coordinate> *p=x.get_coordinates();
  long double x_r = std::inner_product(p->begin(), p->end(), vector_R.begin(), 0.0);
  if(x_r>=0)
    return 1;
  else
    return 0;
}

////////////////////////////////////////////
