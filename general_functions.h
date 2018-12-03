#ifndef GENERAL_FUNCIONS_H
#define GENERAL_FUNCIONS_H

void read_info(std::ifstream *,std::vector<Vector*> *);

int check_parameters(std::ifstream * ,int * ,int * ,int *,int*,int* );
int check_dimensions(std::vector<Vector*> & );
int check_cluster_num(std::vector<Vector*> & ,int);
long double absolute(long double );
long double power(long double,int );
int mod (long long  int , unsigned long long  int );

#endif
