#ifndef EUCL_COS_H
#define EUCL_COS_H

#include <vector>

#include "param.h"
#include "Vector.h"

class Ls_Hashing
{
  protected:
    Ls_Hashing(void);

  public:
    virtual ~Ls_Hashing(void);
};

class Euclidean : public Ls_Hashing
{
    std::vector<double> vector_V;
    double t;
    int w;

  public:
    Euclidean(int ,int );
    ~Euclidean(void);

    int h_func(Vector & );

};

class Cosine : public Ls_Hashing
{
    std::vector<double> vector_R;

    public:
      Cosine(int);
      ~Cosine(void);

      int h_func(Vector &);

};

#endif
