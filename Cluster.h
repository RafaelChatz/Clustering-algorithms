#ifndef CLUSTER_H
#define CLUSTER_H

#include <string>
#include <iostream>
#include <vector>

#include "param.h"
#include "Vector.h"

class Cluster{

  Vector * Centroid;

  public:
    std::vector<Vector *>  Vectors_in_Cluster;

    Cluster(Vector *);
    ~Cluster();
    void add_Vector(Vector *);
    void add_Vectors(std::vector<Vector *> &);

    void update_Centroid();
    Vector * get_centroid();
    std::string get_centroid_id();
};

#endif
