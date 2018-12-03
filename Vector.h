#ifndef VECTOR_H
#define VECTOR_H
#include <string>
#include <iostream>
#include <vector>
#include "param.h"

class Vector {
 std::string vector_id;
 std::vector<coordinate> vector_X;

 public:
   Vector(std::string );
   ~Vector(void);
   void add_coordinate(coordinate );

   std::string get_identity(void);
   std::vector<coordinate> *get_coordinates(void);
};

#endif
