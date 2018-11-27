#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <limits>
#include <cfenv>
#include <iomanip>      // std::setprecision

#include "Vector.h"
#include "general_functions.h"
#include "param.h"
#include "Distance.h"



  void  Distance::Nearest_N(std::vector< std::pair <std::string,double> >  & ,std::vector<Vector*> & ,Vector * ){}
  void Distance::Range_N(std::vector<std::string > &,std::vector<Vector*> & ,Vector * ,double ) {}
  long double Distance::dist(Vector & ,Vector & ){}

  long double Euclidean::dist(Vector & x,Vector & y){
    long double distance=0.0;
    std::vector<coordinate> *x_coordinates=x.get_coordinates();
    std::vector<coordinate> *y_coordinates=y.get_coordinates();
    coordinate* p_x = x_coordinates->data();
    coordinate* p_y = y_coordinates->data();

    for(int i=0;i<x_coordinates->size();i++){
      distance=distance + (long double)power(absolute(p_x[i]-p_y[i]),2);
    }
    return (long double)std::sqrt((long double)distance);
  }

  void Euclidean::Range_N(std::vector<std::string > &vr,std::vector<Vector*> & l,Vector * v,double r) {

    for (int i=0;i<l.size();i++){
      Vector  *s=l.at(i);
      double distance=dist(*v,*s);
      if(distance<r){
        vr.push_back(s->get_identity()); //add all below r to vector
      }
    }
  }

  void Euclidean::Nearest_N(std::vector< std::pair <std::string,double> >  & vr,std::vector<Vector*> & l,Vector * v) {

    double min=std::numeric_limits<double >::max();
    double distance=0;
    Vector *a;
    for (int i=0;i<l.size();i++){//find the nearest
      Vector  *s=l.at(i);
      distance=dist(*v,*s);
      if(distance<min){
        min=distance;
        a=l.at(i);
      }
      if(i==l.size()-1){
        vr.push_back( make_pair(a->get_identity(),min) ); //add it to the vector
      }
    }

  }



  template<typename Iter_T>
  long double vectorNorm(Iter_T first, Iter_T last) {
    return sqrt(std::inner_product(first, last, first, 0.0L));
  }

  long double Cosine::cosi(Vector & x,Vector & y){
    std::vector<coordinate> *x_coordinates=x.get_coordinates();
    std::vector<coordinate> *y_coordinates=y.get_coordinates();

    long double x_norm=vectorNorm(x_coordinates->begin(),x_coordinates->end());
    long double y_norm=vectorNorm(y_coordinates->begin(),y_coordinates->end());

    return (double)std::inner_product(x_coordinates->begin(), x_coordinates->end(), y_coordinates->begin(), 0.0L)/(x_norm*y_norm);
  }

  long double Cosine::dist(Vector & x,Vector & y){ return 1 - cosi(x,y);}

  void Cosine::Range_N(std::vector<std::string > &vr,std::vector<Vector*> & l,Vector * v,double r) {

    for (int i=0;i<l.size();i++){
      Vector  *s=l.at(i);
      double distance=dist(*v,*s);
      if(distance<r){
        vr.push_back(s->get_identity());
      }
    }
  }

  void  Cosine::Nearest_N(std::vector< std::pair <std::string,double> >  & vr,std::vector<Vector*> & l,Vector * v) {

    double min=std::numeric_limits<double >::max();
    double distance=0;
    Vector *a;
    for (int i=0;i<l.size();i++){
      Vector  *s=l.at(i);
      distance=dist(*v,*s);
      if(distance<min){
        min=distance;
        a=l.at(i);
      }
      if(i==l.size()-1){
        vr.push_back( make_pair(a->get_identity(),min) );
      }
    }

  }
