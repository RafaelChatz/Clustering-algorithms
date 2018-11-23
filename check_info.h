#ifndef CHECK_INFO_H
#define CHECK_INFO_H

#include <iostream>
#include <fstream>
#include <string>

void print_err(int );
int check_info_cluster(int ,char**,std::ifstream *,std::ifstream *,std::ofstream *,std::string *);

#endif
