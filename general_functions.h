#ifndef GENERAL_FUNCIONS_H
#define GENERAL_FUNCIONS_H

void read_info(std::ifstream *,std::vector<Vector*> *);
int check_parameters(std::ifstream * ,int * ,int * ,int * );
void k_unique_rand(int ,int ,std::vector<Vector*> &,std::vector<Vector*> &);
void k_means_plus_plus(int ,int ,std::vector<Vector*> &,std::vector<Vector*> &,std::string &);
int check_dimensions(std::vector<Vector*> & );
int check_cluster_num(std::vector<Vector*> & ,int);
long double absolute(long double );
long double power(long double,int );

#endif
