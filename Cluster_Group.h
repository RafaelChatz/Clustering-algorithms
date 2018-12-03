#ifndef CLUSTER_GROUP_H
#define CLUSTER_GROUP_H

#include <string>
#include <iostream>
#include <vector>
#include "param.h"
#include "Cluster.h"
#include "Vector.h"
#include "Distance.h"

class Cluster_Group{

  std::vector<Cluster *> cluster_group;
  int Cluster_num;
  std::string metric;
  Distance **dist;
  long double Silhouette(int);

public:
  Cluster_Group(int,std::string & );
  ~Cluster_Group();
  void k_unique_rand_init(std::vector<Vector*> & );
  void k_means_plus_plus_init(std::vector<Vector*> & );
  void k_means_average_init(std::vector<Vector*> & );

  void Lloyd_assignment(std::vector<Vector*> &);
  void Range_search_assignment_LSH(std::vector<Vector*> &, int ,int );
  void Range_search_assignment_Hypercube(std::vector<Vector*> &,int,int,int);
  int k_means_update();
  int PAM_improvement_update();
  std::string get_centroid_id_n(int );
  void print_info(std::ofstream *,double);
};


#endif
