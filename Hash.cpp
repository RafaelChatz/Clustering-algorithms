#include <iostream>
#include <chrono>
#include <random>
#include <numeric>
#include <limits>
#include <vector>
#include <unordered_map>

#include "Hash.h"
#include "Eucl_Cos.h"
#include "Vector.h"
#include "general_functions.h"
#include "param.h"

extern long int bytes_allocated;

//Hash class
Hash::Hash(int b_num,int k){
  bucket_num=b_num;
bytes_allocated++;
bytes_allocated+=sizeof(Hash);

  bucket = new std::vector<Vector*>[bucket_num];
//  std::cout<<"buckets of hash created: "<<bucket_num<<std::endl;

}

Hash::~Hash(void){
//  std::cout<<"buckets of hash destroyed: "<<bucket_num<<std::endl;
  delete[] bucket;
}



void Hash::insert_Vector(Vector* v){ return;}

std::vector<Vector*> * Hash::get_similar_Vectors(Vector *v){ return NULL;}

int Hash::get_bucket_num(){
  return bucket_num;
}

std::vector<Vector*> * Hash::get_bucket(int pos){
  std::vector<Vector*> *similar_Vectors=&bucket[pos];

  return similar_Vectors;
}

 int Hash::hashFunction(Vector *){ return -1;}
int Hash::g_fun(Vector *,Vector *){ return -1;}
//////////////////////////////////////////////////

//Euclidean_Hash class
Euclidean_Hash::Euclidean_Hash(int b_num,int k,int dimensions): Hash(b_num,k){
//  std::cout<<"Euclidean_Hash:"<<b_num<<std::endl;
bytes_allocated+=sizeof(Euclidean_Hash);

  Euclidean *H_vec;
  r = new int [k];
  int min=0;
  int max=std::numeric_limits<int >::max();

  std::random_device rd;     // only used once to initialise (seed) engine

  std::uniform_int_distribution<int> random(min,max);
  for(int i=0;i<k;i++)
    r[i]=random(rd);
  for(int i = 0; i < k; i += 1){
    H_vec = new Euclidean(W,dimensions);
    H_vector.push_back(H_vec);
  }

}

Euclidean_Hash::~Euclidean_Hash(void){
  Euclidean *H_vec=NULL;

  for (int i=0; i<H_vector.size(); i++)  {
    delete(H_vector[i]);
  }
  delete[] r;

//  std::cout<<"Euclidean_Hash destroyed"<<std::endl;

}

int Euclidean_Hash::hashFunction(Vector *v){

  unsigned int s=0;
  for (int i=0; i<H_vector.size(); i++){
    s=s+mod(r[i]*H_vector[i]->h_func(*v),(unsigned long int)M);
  }
  s=mod(mod(s,(unsigned long int)M),bucket_num);
  return (unsigned int)s;
}

int Euclidean_Hash::g_fun(Vector *v,Vector *s){
  int d=0;

  for (int i=0; i<H_vector.size(); i++){

    if(H_vector[i]->h_func(*v)==H_vector[i]->h_func(*s)){ //check if they are equal
      d++;
    }
  }
  return d==H_vector.size(); //if all are equal return 1
}

void Euclidean_Hash::insert_Vector(Vector* v){

  int  place =hashFunction(v);
  bucket[place].insert (bucket[place].end(), v);
  bytes_allocated+=sizeof(Vector*);

}

std::vector<Vector*> * Euclidean_Hash::get_similar_Vectors(Vector *v){
  int pos=hashFunction(v);
  std::vector<Vector*> *similar_Vectors=&bucket[pos];

  return similar_Vectors;
}

////////////////////////////////////////////////////

//Cosine_Hash class
Cosine_Hash::Cosine_Hash(int k,int dimensions): Hash(pow(2,k),k){
  //std::cout<<"Cosine_Hash:"<<std::endl;
  Cosine *H_vec;
  bytes_allocated+=sizeof(Cosine_Hash);

  for(int i = 0; i < k; i += 1){
    H_vec = new Cosine(dimensions);
    H_vector.push_back(H_vec);
    }
}

Cosine_Hash::~Cosine_Hash(void){
//  std::cout<<"Cosine_Hash destroyed"<<std::endl;
  Cosine *H_vec=NULL;

  for (int i=0; i<H_vector.size(); i++)  {
    delete(H_vector[i]);
  }
}

int Cosine_Hash::hashFunction(Vector *v){

  int s=0;
  for (int i=0; i<H_vector.size(); i++){
    s=s+H_vector[i]->h_func(*v)*pow(2,H_vector.size()-i-1);
  }
  return s;

}

void Cosine_Hash::insert_Vector(Vector *v){

  int  place =hashFunction(v);
  bucket[place].insert (bucket[place].end(), v);
  bytes_allocated+=sizeof(Vector *);

}

std::vector<Vector*> * Cosine_Hash::get_similar_Vectors(Vector *v){
  int pos=hashFunction(v);
  std::vector<Vector*> *similar_Vectors=&bucket[pos];

  return similar_Vectors;
}
////////////////////////////////////////////////////



Cube_Hash::Cube_Hash(int k,int dimensions): Hash(pow(2,k),k){

  int min=0;
  int max=std::numeric_limits<int >::max();
  Euclidean *H_vec;
  bytes_allocated+=sizeof(Cube_Hash);

  for(int i = 0; i < k; i += 1){
    H_vec = new Euclidean(W,dimensions);
    H_vector.push_back(H_vec);
  }
}

Cube_Hash::~Cube_Hash(void){
//  std::cout<<"buckets of hash destroyed: "<<bucket_num<<std::endl;
 // delete[] r;

  Euclidean *H_vec=NULL;

  for (int i=0; i<H_vector.size(); i++)  {
    delete(H_vector[i]);
  }

}


void Cube_Hash::insert_Vector(Vector* v){
  int  place =hashFunction(v);
  bucket[place].insert (bucket[place].end(), v);
  bytes_allocated+=sizeof(Vector*);

}


int Cube_Hash::hashFunction(Vector *v){
  std::random_device rd;
  std::uniform_int_distribution<int> random(0,1);
  srand (time(NULL));
  int s=0;
  for (int i=0; i<H_vector.size(); i++){

    double  value=H_vector[i]->h_func(*v);
    std::unordered_map<double,int>::const_iterator got = h_map.find (value);

    if ( got == h_map.end() ){
      int ran=random(rd);
      std::pair<double,int> h (value,ran);

      h_map.insert (h); //add it to the map
      s=s+ran*pow(2,H_vector.size()-i-1);
      bytes_allocated+=sizeof(std::pair<double,int>);

    }
    else{
      s=s+got->second*pow(2,H_vector.size()-i-1); //get map value

    }
  }
  return s;
}

std::vector<Vector*> * Cube_Hash::get_similar_Vectors(Vector *v){
  int pos=hashFunction(v);
  std::vector<Vector*> *similar_Vectors=&bucket[pos];

  return similar_Vectors;
}
