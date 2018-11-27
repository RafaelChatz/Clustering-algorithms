#ifndef DISTANCE_H
#define DISTANCE_H

#include <iostream>

#include "param.h"
#include "Vector.h"
#include "general_functions.h"

class Distance{

  public:
    virtual  void  Nearest_N(std::vector< std::pair <std::string,double> >  & ,std::vector<Vector*> & ,Vector * ) ;
    virtual  void Range_N(std::vector<std::string > &,std::vector<Vector*> & ,Vector * ,double ) ;
    virtual  long double dist(Vector & ,Vector & );
};


class Euclidean : public Distance{
  public:
    void  Nearest_N(std::vector< std::pair <std::string,double> >  & ,std::vector<Vector*> & ,Vector * ) ;
    void Range_N(std::vector<std::string > &,std::vector<Vector*> & ,Vector * ,double ) ;
    long double dist(Vector & ,Vector & );

};

class Cosine : public Distance{

  private:
    long double cosi(Vector & ,Vector & );

  public:
    void  Nearest_N(std::vector< std::pair <std::string,double> >  & ,std::vector<Vector*> & ,Vector * ) ;
    void Range_N(std::vector<std::string > &,std::vector<Vector*> & ,Vector * ,double ) ;
    long double dist(Vector & ,Vector & );
};

#endif
