#include <iostream>
#include "Cluster.h"

Cluster::Cluster(Vector * centroid){
  Centroid=centroid;
}

Cluster::~Cluster(){
  Vectors_in_Cluster.erase(Vectors_in_Cluster.begin(),Vectors_in_Cluster.end());
}
