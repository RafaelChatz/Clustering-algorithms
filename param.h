#ifndef PARAM_H
#define PARAM_H

#define TRUE 1
#define FALSE 0
#define DIMENSIONS 203 //204 if we count the id
#define CLUSTERS_DIV 10
#define NUMBER_OF_HASH_FUNCTIONS 4
#define NUMBER_OF_HASH_TABLES 5
#define PROBES 2
#define NUMBER_OF_VECTORS 100
#define NUM_TIMES 3
#define PRECISION 10

#define W 149//79//149//might be bigger
#define M pow(2,32)-5
#define TABLE_DIV 25

typedef long double coordinate ;
extern long int bytes_allocated;

#endif
