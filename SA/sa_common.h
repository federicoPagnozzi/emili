#ifndef SA_COMMON_H
#define SA_COMMON_H


#define nullptr NULL

typedef struct _sa_status {

  int    counter;
  int    total_counter;
  int    accepted;
  float  rate;
  short *last_accepted;
  int    tenure;
  int    index;

} sa_status;

#endif
