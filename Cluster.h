#ifndef CLUSTER_H
#define CLUSTER_H

#include <string>
#include <iostream>
#include <vector>
#include "param.h"

class Cluster{

  Vector * Centroid;
  std::vector<Vector *>  Vectors_in_Cluster;

  public:
    Cluster(Vector *);
    ~Cluster();
    void add_Vector(Vector *);
    void update_Centroid();
};

#endif
