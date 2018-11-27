#ifndef CLUSTER_GROUP_H
#define CLUSTER_GROUP_H

#include <string>
#include <iostream>
#include <vector>
#include "param.h"
#include "Cluster.h"

class Cluster_Group{

  Cluster * Cluster_group;

public:
  Cluster_group(int);
  ~Cluster_Group();
};


#endif
