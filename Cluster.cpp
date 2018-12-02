#include <iostream>
#include <vector>

#include "Cluster.h"
#include "Vector.h"

Cluster::Cluster(Vector * centroid){
  Centroid=centroid;
}

Cluster::~Cluster(){

  Vectors_in_Cluster.erase(Vectors_in_Cluster.begin(),Vectors_in_Cluster.end());
}

void Cluster::add_Vector(Vector * vector){

  Vectors_in_Cluster.insert (Vectors_in_Cluster.end(),vector);
}

Vector * Cluster::get_centroid(){

  return Centroid;
}

std::string Cluster::get_centroid_id(){

  return Centroid->get_identity();
}

void Cluster::add_Vectors(std::vector<Vector *>  &vectors){

  for(int i=0;i<vectors.size();i++)
    Vectors_in_Cluster.insert (Vectors_in_Cluster.end(),vectors.at(i));
}
