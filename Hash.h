#ifndef HASH_H
#define HASH_H

#include<iostream>
#include <vector>
#include <unordered_map>

#include "param.h"
#include "Eucl_Cos.h"
#include "Vector.h"

class Hash
{
  protected:
    int bucket_num;
    std::vector<Vector*> *bucket;

    Hash(int,int);

  public:
    virtual void insert_Vector(Vector*);
    virtual ~Hash(void);
    virtual std::vector<Vector*> * get_similar_Vectors(Vector*);
    int get_bucket_num();
    std::vector<Vector*> * get_bucket(int);
    virtual int hashFunction(Vector *);
    virtual int g_fun(Vector *,Vector *);

};

class Euclidean_Hash : public Hash
{
    int * r;

    std::vector<Euclidean*>  H_vector;

  public:
    Euclidean_Hash(int,int,int);
    ~Euclidean_Hash(void);
    void insert_Vector(Vector*);
    std::vector<Vector*> * get_similar_Vectors(Vector*);
    int hashFunction(Vector *);
    int g_fun(Vector *,Vector *);

};

class Cosine_Hash : public Hash
{
    std::vector<Cosine*>  H_vector;


    public:
      Cosine_Hash(int,int);
      ~Cosine_Hash(void);
      void insert_Vector(Vector*);
      std::vector<Vector*> * get_similar_Vectors(Vector*);
      int hashFunction(Vector*);

};

class Cube_Hash : public Hash
{
    std::vector<Euclidean*>  H_vector;
    std::unordered_map<double,int> h_map;
  public:
    Cube_Hash(int,int);
    ~Cube_Hash(void);
    void insert_Vector(Vector*);
    int hashFunction(Vector *);
    std::vector<Vector*> * get_similar_Vectors(Vector*);

};

#endif
