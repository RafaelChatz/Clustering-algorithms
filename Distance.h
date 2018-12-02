#ifndef DISTANCE_H
#define DISTANCE_H

#include <iostream>

#include "param.h"
#include "Vector.h"
#include "general_functions.h"

class Distance{

  public:
    virtual  void  Nearest_N(std::vector< std::pair <std::string,double> >  & ,std::vector<Vector*> & ,Vector * ) ;
    virtual  void Range_N(std::vector<Vector*> &,std::vector<Vector*> & ,Vector * ,long double ) ;
    virtual void Range_bet_N(std::vector<Vector*> &,std::vector<Vector*> & ,Vector * ,long double,long double ) ;
    virtual  long double dist(Vector & ,Vector & );
};


class Euclidean_Distance : public Distance{
  public:
    void  Nearest_N(std::vector< std::pair <std::string,double> >  & ,std::vector<Vector*> & ,Vector * ) ;
    void Range_N(std::vector<Vector*> &,std::vector<Vector*> & ,Vector * ,long double ) ;
    void Range_bet_N(std::vector<Vector*> &,std::vector<Vector*> & ,Vector * ,long double,long double ) ;

    long double dist(Vector & ,Vector & );

};

class Cosine_Distance : public Distance{

  private:
    long double cosi(Vector & ,Vector & );

  public:
    void  Nearest_N(std::vector< std::pair <std::string,double> >  & ,std::vector<Vector*> & ,Vector * ) ;
    void Range_N(std::vector<Vector*> &,std::vector<Vector*> & ,Vector * ,long double ) ;
    void Range_bet_N(std::vector<Vector*> &,std::vector<Vector*> & ,Vector * ,long double,long double ) ;
    long double dist(Vector & ,Vector & );
};

#endif
