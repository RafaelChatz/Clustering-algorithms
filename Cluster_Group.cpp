#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <iomanip>      // std::setprecision>
#include <limits>
#include <unordered_set>
#include <utility>
#include <functional>

#include "Cluster.h"
#include "Cluster_Group.h"
#include "Vector.h"
#include "Distance.h"
#include "Hash.h"

Cluster_Group::Cluster_Group(int cluster_num,std::string & metric){
  Cluster_num=cluster_num;
  this->metric=metric;
  dist= new Distance*[1];
  if(metric.compare("euclidean")==0)
    dist[0] = new Euclidean_Distance();
  else
    dist[0] = new Cosine_Distance();
}

Cluster_Group::~Cluster_Group(){

  for (std::vector<Cluster *>::iterator it=cluster_group.begin(); it<cluster_group.end(); it++)
   delete(it[0]);

  delete(dist[0]);
  delete[] dist;
}

void Cluster_Group::k_unique_rand_init(std::vector<Vector*> & all_vectors){

  int n=all_vectors.size();

  std::random_device rd;
  std::uniform_int_distribution<int> distribution(0,2*n);
  int pos[Cluster_num];
  int im = 0;

  for (int in = 0; in < n && im < Cluster_num; ++in) { //take unique random numbers that corespond to the location in the vector
    int rn = n - in;
    int rm = Cluster_num - im;
    if (distribution(rd) % rn < rm)
      pos[im++] = in ;
  }

  for(int i=0;i<Cluster_num;i++){
    Cluster * Cl= new Cluster (all_vectors.at(pos[i]));
    cluster_group.insert (cluster_group.end(), Cl); //and add them as centroids in each cluster
  }
}

void Cluster_Group::k_means_average_init(std::vector<Vector*> & all_vectors){
//extra algorithm
  int n=all_vectors.size();

  std::random_device rd;
  std::uniform_int_distribution<int> distribution_int(0,n);
  int pos[Cluster_num];

  long double** distance = new long double*[Cluster_num-1];
  for(int i = 0; i < Cluster_num-1; ++i)
    distance[i] = new long double[n];

  long double* average_distance = new long double[n];

  pos[0]=distribution_int(rd); //take 1 centroid at random

  Cluster * Cl1= new Cluster (all_vectors.at(pos[0]));
  cluster_group.insert (cluster_group.end(), Cl1);
  //find average distance for all vectors from every centroid to find the propability
  for(int j=1;j<Cluster_num;j++){
    long double all_distance=0;

    for(int i=0;i<all_vectors.size();i++){
     average_distance[i]=0;
      int fl=0;
      for(int d=0;d<j;d++){
        if(pos[d]==i){
          fl=1;
          break;
        }
      }
      if(fl==1)
        continue;

      distance[j-1][i]=(long double)dist[0]->dist(*cluster_group.at(j-1)->get_centroid(),*all_vectors.at(i));

      for(int d=0;d<j;d++)
        average_distance[i]=average_distance[i] +distance[d][i];

      average_distance[i]=average_distance[i]/(long double)j;
      all_distance=all_distance+average_distance[i]; //sum of avergae distances
    }

    std::uniform_real_distribution<long double> distribution_real(0.0,(long double)all_distance);
    long double prob=distribution_real(rd);

    int s;
    for(s=0;s<all_vectors.size();s++){
      prob=prob-average_distance[s];//take the one that makes the propability below 0
      if(prob<0.0)
        break;
    }
    pos[j]=s;
    Cluster * Cl= new Cluster (all_vectors.at(pos[j]));
    cluster_group.insert (cluster_group.end(), Cl);
  }

  for(int i = 0; i < Cluster_num-1; ++i)
    delete [] distance[i];
  delete [] distance;
  delete [] average_distance;

}

void Cluster_Group::k_means_plus_plus_init(std::vector<Vector*> & all_vectors){
//the same as above only this time we use the distance that k-means_plus_plus requires for the propability
  int n=all_vectors.size();

  std::random_device rd;
  std::uniform_int_distribution<int> distribution_int(0,n);
  int pos[Cluster_num];

  long double** distance = new long double*[Cluster_num-1];
  for(int i = 0; i < Cluster_num-1; ++i)
    distance[i] = new long double[n];

  long double* nearest_distance = new long double[n];

  pos[0]=distribution_int(rd);

  Cluster * Cl1= new Cluster (all_vectors.at(pos[0]));
  cluster_group.insert (cluster_group.end(), Cl1);

  for(int j=1;j<Cluster_num;j++){
    long double all_distance=0;

    for(int i=0;i<all_vectors.size();i++){
     nearest_distance[i]=0;
      int fl=0;
      for(int d=0;d<j;d++){
        if(pos[d]==i){
          fl=1;
          break;
        }
      }
      if(fl==1)
        continue;

      distance[j-1][i]=(long double)dist[0]->dist(*cluster_group.at(j-1)->get_centroid(),*all_vectors.at(i));

      for(int d=0;d<j;d++){
        nearest_distance[i]=distance[d][i];
        if(d!=0)
          if(distance[d][i]<nearest_distance[i])
            nearest_distance[i]=distance[d][i];
      }
      all_distance=all_distance+nearest_distance[i];
    }

    std::uniform_real_distribution<long double> distribution_real(0.0,(long double)all_distance);
    long double prob=distribution_real(rd);

    int s;
    for(s=0;s<all_vectors.size();s++){
      prob=prob-nearest_distance[s];
      if(prob<0.0)
        break;
    }
    pos[j]=s;
    Cluster * Cl= new Cluster (all_vectors.at(pos[j]));
    cluster_group.insert (cluster_group.end(), Cl);
  }

  for(int i = 0; i < Cluster_num-1; ++i)
    delete [] distance[i];
  delete [] distance;

  delete [] nearest_distance;

}

void Cluster_Group::Lloyd_assignment(std::vector<Vector*> & all_vectors){

  for(int i=0;i<cluster_group.size();i++){ //make sure that cluster has no vectors before assignment
    cluster_group.at(i)->empty_vectors_in_cluster();
  }

  if(cluster_group.size()<1)
    return ;
    for(int i=0;i<all_vectors.size();i++){
      int pos;
      long double min_distance=std::numeric_limits<long double >::max();
      for(int j=0;j<cluster_group.size();j++){//find min distance
        if(!cluster_group.at(j)->get_centroid_id().compare(all_vectors.at(i)->get_identity())){
          pos=-1;
          break;
        }
        long double distance=(long double)dist[0]->dist(*cluster_group.at(j)->get_centroid(),*all_vectors.at(i));
        if(distance<min_distance){
          min_distance=distance;
          pos=j;
        }
      }

      if(pos!=-1) //if vector is not a centroid  assign vector to the  cluster
        cluster_group.at(pos)->add_Vector(all_vectors.at(i));
    }


}

void Cluster_Group::Range_search_assignment_LSH(std::vector<Vector*> & all_vectors, int L,int k){

  for(int i=0;i<cluster_group.size();i++){
    cluster_group.at(i)->empty_vectors_in_cluster();
  }

  if(cluster_group.size()<1)
    return ;

  int dimensions=all_vectors.at(0)->get_coordinates()->size();
  int vectors=all_vectors.size();

  Hash **h_tables= new Hash*[L];

  if(metric.compare("euclidean")==0){ //hashing according to the metric
    for(int j=0 ; j<L ; j++)
       h_tables[j] = new Euclidean_Hash(vectors/TABLE_DIV,k,dimensions);
  }else{
    for(int j=0 ; j<L ; j++)
       h_tables[j] = new Cosine_Hash(k,dimensions);
  }

  for(int i = 0; i < all_vectors.size(); i ++) { //add to the hash table everything but the centroids
    int fl=0;
    for(int j=0;j<cluster_group.size();j++)
      if(all_vectors[i]->get_identity().compare(cluster_group[j]->get_centroid()->get_identity())==0){
        fl=-1;
        break;
      }

    for(int j=0 ; j<L ; j++){
      if(fl==-1)
        break;
      h_tables[j]->insert_Vector(all_vectors[i]);
    }
  }

  long double minimum_range =std::numeric_limits<long double >::max();

  for(int i=0;i<cluster_group.size();i++){
    for(int j=i+1;j<cluster_group.size();j++){
      long double range=(long double)dist[0]->dist(*cluster_group.at(i)->get_centroid(),*cluster_group.at(j)->get_centroid());
      if(range<minimum_range)
        minimum_range=range;
    }
  }//find minimum distance beetwen centroids


  std::vector<Vector*> *Vectors_passed = new std::vector<Vector*>[cluster_group.size()];


  if(metric.compare("euclidean")==0){//take bucket and use g to fine evey vector that we need to check with range assignment
  for(int i=0;i<cluster_group.size();i++)
    for(int j=0 ; j<L ; j++){
      std::vector<Vector*> * similar_Vectors=h_tables[j]->get_similar_Vectors(cluster_group.at(i)->get_centroid()); //get similar vectos from the hash function
      Vector** Vec = similar_Vectors->data();
      for(int s=0;s<similar_Vectors->size();s++)
        if(h_tables[j]->g_fun(Vec[s],cluster_group.at(i)->get_centroid()))
          if(std::find(Vectors_passed[i].begin(), Vectors_passed[i].end(), Vec[s]) == Vectors_passed[i].end())
            Vectors_passed[i].push_back(Vec[s]);
    }
  }
  else{//we dont use g for cosine
    for(int i=0;i<cluster_group.size();i++)
      for(int j=0 ; j<L ; j++){
        std::vector<Vector*> * similar_Vectors=h_tables[j]->get_similar_Vectors(cluster_group.at(i)->get_centroid()); //get similar vectos from the hash function
        Vector** Vec = similar_Vectors->data();//then we take all the VEctors with the same second hash function
        for(int s=0;s<similar_Vectors->size();s++)
          if(std::find(Vectors_passed[i].begin(), Vectors_passed[i].end(), Vec[s]) == Vectors_passed[i].end())
            Vectors_passed[i].push_back(Vec[s]);

      }
  }

    std::vector<Vector*> *vectors_in_range = new std::vector<Vector*>[cluster_group.size()];
    std::unordered_set<Vector*> assigned_Vectors ;

    for(int i=0;i<cluster_group.size();i++){//find every vector that appears in more than one cluster and assign it to the nearest

      for(int s=0;s<Vectors_passed[i].size();s++){

        int flc=0;
        long double cdist=(long double)dist[0]->dist(*cluster_group.at(i)->get_centroid(),*Vectors_passed[i].at(s));
        int cl=i;
        for(int j=i+1;j<cluster_group.size();j++){
          if(std::find(Vectors_passed[j].begin(), Vectors_passed[j].end(), Vectors_passed[i].at(s)) != Vectors_passed[j].end()){
            flc=1;
            long double c2dist=(long double)dist[0]->dist(*cluster_group.at(j)->get_centroid(),*Vectors_passed[i].at(s));
            if(c2dist<cdist)
              cl=j;
          }

        }

        if(flc==1){
            std::unordered_set<Vector*>::const_iterator got = assigned_Vectors.find (Vectors_passed[i].at(s));
            if ( got == assigned_Vectors.end() ){
              cluster_group.at(cl)->add_Vector(Vectors_passed[i].at(s));
              assigned_Vectors.insert(Vectors_passed[i].at(s));
            }
        }

      }

    }

    int *fln= new int [cluster_group.size()];

    for(int i=0;i<cluster_group.size();i++){ //assign vectors to clusters with range check
      dist[0]->Range_bet_N(vectors_in_range[i],Vectors_passed[i] ,cluster_group.at(i)->get_centroid(),0,minimum_range/2 ) ;
      fln[i]=vectors_in_range[i].size();
      for(int j=0;j<vectors_in_range[i].size();j++){
        std::unordered_set<Vector*>::const_iterator got = assigned_Vectors.find (vectors_in_range[i].at(j));
        if ( got == assigned_Vectors.end() ){
          cluster_group.at(i)->add_Vector(vectors_in_range[i].at(j));
          assigned_Vectors.insert(vectors_in_range[i].at(j));
        }
      }
      vectors_in_range[i].clear();
    }

    int n=1;
    while(TRUE){//until there are no more
      int st=0;
      for(int i=0;i<cluster_group.size();i++){ //double the range every time and search for vectors to add to clusters

        if(fln[i]>=Vectors_passed[i].size()){
          st++;
          continue;
        }
        dist[0]->Range_bet_N(vectors_in_range[i],Vectors_passed[i] ,cluster_group.at(i)->get_centroid(),n*2*minimum_range/4 ,(n+1)*2*minimum_range/4 ) ;
        fln[i]=fln[i]+vectors_in_range[i].size();
        for(int j=0;j<vectors_in_range[i].size();j++){
          std::unordered_set<Vector*>::const_iterator got = assigned_Vectors.find (vectors_in_range[i].at(j));
            if ( got == assigned_Vectors.end() ){
              cluster_group.at(i)->add_Vector(vectors_in_range[i].at(j));
              assigned_Vectors.insert(vectors_in_range[i].at(j));
            }
        }
        vectors_in_range[i].clear();

      }
      if(st==cluster_group.size())
        break;
      n++;
    }

    for(int i=0;i<all_vectors.size();i++){ //assign everything outside of the bucket with Lloyds
      std::unordered_set<Vector*>::const_iterator got = assigned_Vectors.find (all_vectors.at(i));
      if ( got != assigned_Vectors.end() )
        continue;

      int pos;
      long double min_distance=std::numeric_limits<long double >::max();
      for(int j=0;j<cluster_group.size();j++){
        if(!cluster_group.at(j)->get_centroid_id().compare(all_vectors.at(i)->get_identity())){
          pos=-1;
          break;
        }
        long double distance=(long double)dist[0]->dist(*cluster_group.at(j)->get_centroid(),*all_vectors.at(i));
        if(distance<min_distance){
          min_distance=distance;
          pos=j;
        }
      }
      if(pos!=-1)
        cluster_group.at(pos)->add_Vector(all_vectors.at(i));

    }

    for(int j=0 ; j<L ; j++)
        delete(h_tables[j]);
    delete[] h_tables;
    delete[] vectors_in_range;
    delete[] Vectors_passed;
    delete[] fln;
}

void Cluster_Group::Range_search_assignment_Hypercube(std::vector<Vector*> & all_vectors,int k,int Ms,int probes){
//everything is the same here (the only difference is the probes )
  for(int i=0;i<cluster_group.size();i++){
    cluster_group.at(i)->empty_vectors_in_cluster();
  }

  if(cluster_group.size()<1)
    return ;

    int dimensions=all_vectors.at(0)->get_coordinates()->size();
    int vectors=all_vectors.size();

    Hash **h_tables= new Hash*[1];

    if(metric.compare("euclidean")==0){
         h_tables[0] = new Cube_Hash(k,dimensions);
    }else{
         h_tables[0] = new Cosine_Hash(k,dimensions);
    }
    int L=1;

    for(int i = 0; i < all_vectors.size(); i ++) {
      int fl=0;
      for(int j=0;j<cluster_group.size();j++)
        if(all_vectors[i]->get_identity().compare(cluster_group[j]->get_centroid()->get_identity())==0){
          fl=-1;
          break;
        }

      for(int j=0 ; j<1 ; j++){
        if(fl==-1)
          break;
          h_tables[0]->insert_Vector(all_vectors[i]);
      }
    }

    long double minimum_range =std::numeric_limits<long double >::max();

    for(int i=0;i<cluster_group.size();i++){
      for(int j=i+1;j<cluster_group.size();j++){
        long double range=(long double)dist[0]->dist(*cluster_group.at(i)->get_centroid(),*cluster_group.at(j)->get_centroid());
        if(range<minimum_range)
          minimum_range=range;
      }
    }

    std::vector<Vector*> *Vectors_passed = new std::vector<Vector*>[cluster_group.size()];
    std::unordered_set<Vector*> assigned_Vectors ;


    for(int i=0;i<cluster_group.size();i++){ //find every vector that we need to check with range search (until we have M or probes are done)
      int bucket=h_tables[0]->hashFunction(cluster_group.at(i)->get_centroid());
      Vectors_passed[i]=*h_tables[0]->get_similar_Vectors(cluster_group.at(i)->get_centroid()); //get similar vectos from the hash function
      int vec_size=Vectors_passed[i].size();
      if(Ms<vec_size)
        continue;
      for(int g=0;g<probes;g++){
        if(Ms>vec_size){
          int neighbor=bucket|(1u << k-g-1);
          if(neighbor==bucket)
            neighbor=bucket&~(1u << k-g-1);
          std::vector<Vector*> Vectors_n=*h_tables[0]->get_bucket(neighbor);
          Vectors_passed[i].insert(Vectors_passed[i].end(), Vectors_n.begin(),Vectors_n.end());
          vec_size=Vectors_passed[i].size();
        }
        else{
          break;
        }
      }
    }

    for(int i=0;i<cluster_group.size();i++){

      for(int s=0;s<Vectors_passed[i].size();s++){

        int flc=0;
        long double cdist=(long double)dist[0]->dist(*cluster_group.at(i)->get_centroid(),*Vectors_passed[i].at(s));
        int cl=i;
        for(int j=i+1;j<cluster_group.size();j++){
          if(std::find(Vectors_passed[j].begin(), Vectors_passed[j].end(), Vectors_passed[i].at(s)) != Vectors_passed[j].end()){
            flc=1;
            long double c2dist=(long double)dist[0]->dist(*cluster_group.at(j)->get_centroid(),*Vectors_passed[i].at(s));
            if(c2dist<cdist)
              cl=j;
          }

        }

        if(flc==1){
            std::unordered_set<Vector*>::const_iterator got = assigned_Vectors.find (Vectors_passed[i].at(s));
            if ( got == assigned_Vectors.end() ){
              cluster_group.at(cl)->add_Vector(Vectors_passed[i].at(s));
              assigned_Vectors.insert(Vectors_passed[i].at(s));
            }
        }

      }

    }

    std::vector<Vector*> *vectors_in_range = new std::vector<Vector*>[cluster_group.size()];

    int *fln= new int [cluster_group.size()];

    for(int i=0;i<cluster_group.size();i++){
      dist[0]->Range_bet_N(vectors_in_range[i],Vectors_passed[i] ,cluster_group.at(i)->get_centroid(),0,minimum_range/2 ) ;
      fln[i]=vectors_in_range[i].size();
      for(int j=0;j<vectors_in_range[i].size();j++){
        std::unordered_set<Vector*>::const_iterator got = assigned_Vectors.find (vectors_in_range[i].at(j));
        if ( got == assigned_Vectors.end() ){
          cluster_group.at(i)->add_Vector(vectors_in_range[i].at(j));
          assigned_Vectors.insert(vectors_in_range[i].at(j));
        }
      }
      vectors_in_range[i].clear();
    }

    int n=1;
    while(TRUE){
      int st=0;
      for(int i=0;i<cluster_group.size();i++){

        if(fln[i]>=Vectors_passed[i].size()){
          st++;
          continue;
        }
        dist[0]->Range_bet_N(vectors_in_range[i],Vectors_passed[i] ,cluster_group.at(i)->get_centroid(),n*2*minimum_range/4 ,(n+1)*2*minimum_range/4 ) ;
        fln[i]=fln[i]+vectors_in_range[i].size();
        for(int j=0;j<vectors_in_range[i].size();j++){
          std::unordered_set<Vector*>::const_iterator got = assigned_Vectors.find (vectors_in_range[i].at(j));
            if ( got == assigned_Vectors.end() ){
              cluster_group.at(i)->add_Vector(vectors_in_range[i].at(j));
              assigned_Vectors.insert(vectors_in_range[i].at(j));
            }
        }
        vectors_in_range[i].clear();

      }
      if(st==cluster_group.size())
        break;
      n++;
    }


    for(int i=0;i<all_vectors.size();i++){
      std::unordered_set<Vector*>::const_iterator got = assigned_Vectors.find (all_vectors.at(i));
      if ( got != assigned_Vectors.end() )
        continue;

      int pos;
      long double min_distance=std::numeric_limits<long double >::max();
      for(int j=0;j<cluster_group.size();j++){
        if(!cluster_group.at(j)->get_centroid_id().compare(all_vectors.at(i)->get_identity())){
          pos=-1;
          break;
        }
        long double distance=(long double)dist[0]->dist(*cluster_group.at(j)->get_centroid(),*all_vectors.at(i));
        if(distance<min_distance){
          min_distance=distance;
          pos=j;
        }
      }
      if(pos!=-1)
        cluster_group.at(pos)->add_Vector(all_vectors.at(i));

    }

    for(int j=0 ; j<1 ; j++)
        delete(h_tables[j]);
    delete[] h_tables;
    delete[] Vectors_passed;
    delete[] vectors_in_range;
    delete[]  fln;
}

int Cluster_Group::k_means_update(){

  int r=0;
  for(int i=0;i<cluster_group.size();i++){ //assign a new centroid to every cluster

    std::string id="CL"+std::to_string(i);
    Vector* Cl = new Vector (id);
    r=r+cluster_group.at(i)->update_Centroid(Cl);
  }

  return r;
}

int Cluster_Group::PAM_improvement_update(){

  int changes=0;
  std::unordered_map<std::string,long double> distances; //we store distances to a map so we wont have to find them each time

  for(int i=0;i<cluster_group.size();i++){//find medoid and use it as centroid

    std::vector<Vector *> * cluster_vectors=cluster_group.at(i)->get_cluster_vectors();
    Vector * cen=cluster_group.at(i)->get_centroid();
    long double min=0;
    int pos=-1;
    for(int j=0;j<cluster_vectors->size();j++){
      long double d=(long double)dist[0]->dist(*cen,*cluster_vectors->at(j));
      std::string idd="C"+std::to_string(j);
      distances.insert (std::make_pair(idd,d));
      min=min+d;
    }
    std::unordered_map<std::string,long double>::const_iterator got ;
    for(int j=0;j<cluster_vectors->size();j++){
      std::string idd="C"+std::to_string(j);
      long double pos_min=0;

      got=distances.find(idd);
      if ( got != distances.end() )
        pos_min=got->second;

      for(int k=0;k<cluster_vectors->size();k++){
        if(j==k)
          continue;

          if(j<k){
            idd=std::to_string(j)+std::to_string(k);
          }
          else{
            idd=std::to_string(k)+std::to_string(j);
          }
          got=distances.find(idd);
          if ( got != distances.end() )
            pos_min=pos_min+got->second;
          else{
            long double d=(long double)dist[0]->dist(*cen,*cluster_vectors->at(j));
            distances.insert (std::make_pair(idd,d));

            pos_min=pos_min+d;
          }
      }

      if(pos_min<min){ //position with min  distance
        min=pos_min;
        pos=j;
      }
    }
    if(pos!=-1){//if not current add the new one
      changes++;
      cluster_group.at(i)->update_Centroid(cluster_vectors->at(pos));
    }


  }

  return changes;
}

long double Cluster_Group::Silhouette(int cl_n){

  long double a=0;
  long double b=0;

  std::vector<Vector *> * cluster_vectors=cluster_group.at(cl_n)->get_cluster_vectors();
  std::unordered_map<std::string,long double> distances;//add only distances found for a
  //we dont need to do that for b because we wont use them again
  std::string idd;
  long double average=0;
  for(int i=0;i<cluster_vectors->size();i++){

    for(int j=0;j<cluster_vectors->size();j++){
      if(i<j){
        idd=std::to_string(i)+std::to_string(j);
      }
      else if(i>j){
        idd=std::to_string(j)+std::to_string(i);
      }else
        continue;

      std::unordered_map<std::string,long double>::const_iterator got ;

      got=distances.find(idd);
      if ( got != distances.end() )
        a=a+got->second;
      else{
        long double d=(long double)dist[0]->dist(*cluster_vectors->at(i),*cluster_vectors->at(j));
        distances.insert (std::make_pair(idd,d));
        a=a+d;
      }

    }

    a=a/cluster_vectors->size();
    int ncl=0;
    long double ndist=std::numeric_limits<long double >::max();

    for(int j=0;j<cluster_group.size();j++){
      if(j==cl_n)
      continue;
      long double di=(long double)dist[0]->dist(*cluster_vectors->at(i),*cluster_group.at(j)->get_centroid());

      if(di<ndist){
        ncl=j;
        di=ndist;
      }

    }

    std::vector<Vector *> * cluster_vectors_b=cluster_group.at(ncl)->get_cluster_vectors();
    for(int j=0;j<cluster_vectors_b->size();j++){
      b=b+(long double)dist[0]->dist(*cluster_vectors->at(i),*cluster_vectors_b->at(j));
    }
    b=b/cluster_vectors_b->size();

    average=average+(b-a)/std::max(a,b);
  }
  return average/cluster_vectors->size();

}

void Cluster_Group::print_info(std::ofstream *outfile,double atime){ //print the info in the output file

  for(int j=0;j<cluster_group.size();j++){
    std::vector<Vector *> *vec=cluster_group.at(j)->get_cluster_vectors();
    *outfile<<"CLUSTER- "<<j+1<<" {size: "<<vec->size()+1<<" centroid: ";

    if(cluster_group.at(j)->get_centroid_id().compare(0,2,"CL")==0)
    {
      //*outfile<<cluster_group.at(j)->get_centroid_id()<<std::endl;
      Vector * vcl=cluster_group.at(j)->get_centroid();
      std::vector<coordinate> * vclc=vcl->get_coordinates();
      for(int i=0;i<vclc->size();i++){
        *outfile<< std::setprecision(PRECISION) <<" "<<vclc->at(i);
      }
      *outfile<<"}"<<std::endl<<std::endl;
    }
    else
      *outfile<<cluster_group.at(j)->get_centroid_id()<<std::endl<<std::endl;

  }
  *outfile<<"clustering_time: "<<atime<<std::endl;
  *outfile<<"Silhouette: [";
  long double Silhouette_av=0;
  for(int j=0;j<cluster_group.size();j++){
    long double Silh=Silhouette(j);
      *outfile<<Silh<<",";
      Silhouette_av=Silhouette_av+Silh;
  }
  *outfile<<Silhouette_av/(long double)cluster_group.size()<<"]";

  *outfile<<std::endl<<std::endl<<std::endl;

  for(int j=0;j<cluster_group.size();j++){
    std::vector<Vector *> *vec=cluster_group.at(j)->get_cluster_vectors();
    *outfile<<"CLUSTER- "<<j+1<<"{";
    for(int i=0;i<vec->size();i++){
      if(i==vec->size()-1)
      {
        *outfile<<vec->at(i)->get_identity();
        continue;
      }
      *outfile<<vec->at(i)->get_identity()<<",";
    }
    *outfile<<"}"<<std::endl;

  }

  *outfile<<"___________________________________________________________________"<<std::endl<<std::endl;

}

std::string Cluster_Group::get_centroid_id_n(int n){

  if(n<Cluster_num)
    return cluster_group.at(n)->get_centroid_id();
  return NULL;
}
