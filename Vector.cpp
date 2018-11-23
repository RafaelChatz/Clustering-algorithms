#include <iostream>
#include <vector>
#include <string>
#include "Vector.h"
#include "param.h"

long int bytes_allocated=0;

Vector::Vector(std::string id){
  vector_id=id;
}

Vector::~Vector(void){

}

void Vector::add_coordinate(coordinate X){
  vector_X.insert (vector_X.end(), X);

}

std::string Vector::get_identity(void){
  return vector_id;
}

std::vector<coordinate> *Vector::get_coordinates(void){
  return &vector_X;
}
