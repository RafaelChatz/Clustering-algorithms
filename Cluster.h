#ifndef CLUSTER_H
#define CLUSTER_H

#include <string>
#include <iostream>
#include <vector>

#include "param.h"
#include "Vector.h"

class Cluster{

  Vector * Centroid;
  std::vector<Vector *>  Vectors_in_Cluster;

  public:

    Cluster(Vector *);
    ~Cluster();
    void add_Vector(Vector *);
    void add_Vectors(std::vector<Vector *> &);
    int update_Centroid(Vector * );
    Vector * get_centroid();
    std::string get_centroid_id();
    std::vector<Vector *> * get_cluster_vectors();
};

#endif
