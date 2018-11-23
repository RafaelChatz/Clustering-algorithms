#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

#include "check_info.h"
#include "param.h"

void print_err(int err){
  if(err==-1)
    std::cerr <<"Not enough arguments-(Value of error: "<<err<<")"<<std::endl<<std::endl;
  else if(err==-2)
    std::cerr <<"Two or more flags have same value-(Value of error: "<<err<<")"<<std::endl<<std::endl;
  else if(err==-4)
    std::cerr <<"Wrong flag-(Value of error: "<<err<<")"<<std::endl<<std::endl;
  else if(err==-5)
    std::cerr <<"A file does not exist or you do not have permission to open it-(Value of error: "<<err<<")"<<std::endl<<std::endl;
  else if(err==-6)
    std::cerr <<"Metric not availiable-(Value of error: "<<err<<")"<<std::endl<<std::endl;

  return;
}

int check_info_cluster(int argc,char** argv,std::ifstream *infile,std::ifstream *confile,std::ofstream *outfile,std::string *metric){

  int err=0;
  int flags=(argc-1)/2;
  char * arg[flags];

  //  $./cluster –i <input file> –c <configuration file> -o <output file> -d
  //  <metric>
  if(argc<9){
    err=-1;
    print_err(err);
  }


    for(int i=0;i<flags;i++){
      arg[i]=(char*)malloc((strlen(argv[1+2*i])+1)*sizeof(char));
      strcpy(arg[i],argv[1+2*i]);
    }

    for(int i=0;i<flags;i++){ //check for mistakes
      for(int j=i+1;j<flags;j++){
        if(!strcmp(arg[i],arg[j])){
          err=-2;
          print_err(err);
          break;
        }
      }
      if(err==-2)
        break;
    }

    if(err<0){
      for(int i=0;i<flags;i++){
        free(arg[i]);
      }
      return err;
    }


    for(int i=0;i<flags;i++){//take the values
      if(!strcmp(arg[i],"-i")){
        infile->open(argv[2+2*i]);
        if (!infile->is_open()) { err=-5; print_err(err); }
      }
      else if(!strcmp(arg[i],"-c")){
        confile->open(argv[2+2*i]);
        if (!confile->is_open()) { err=-5; print_err(err); }
      }
      else if(!strcmp(arg[i],"-d")){
        *metric=argv[2+2*i];
      }
      else if(!strcmp(arg[i],"-o")){
        outfile->open(argv[2+2*i]);
      }
      else{
        err=-4;
        print_err(err);
        break;
      }
    }

  for(int i=0;i<flags;i++){
    free(arg[i]);
  }
  return err;
}
