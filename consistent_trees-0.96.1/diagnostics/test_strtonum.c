#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include "strtonum.c"
#define CONV_64BIT
#include "strtonum.c"


inline test_buf(char *buffer) {
  char *ep1;
  char *ep2;
  double d1, d2;
  float f1, f2;
  int32_t i1, i2;
  int64_t j1, j2;

  d1 = strtodouble(buffer, &ep1);
  d2 = strtod(buffer, &ep2);
  assert(ep1==ep2);
  if (d2 && (!(fabs(d1-d2)/fabs(d1)<DBL_EPSILON))) {
    printf("%s %e %e\n", buffer, d1-d2, DBL_EPSILON);
  }
  
  f1 = strtofloat(buffer, &ep1);
  f2 = strtof(buffer, &ep2);
  assert(ep1==ep2);
  if (d2 && (!(fabs(f1-f2)/fabs(d1)<FLT_EPSILON))) {
    printf("%s %e %e %e %e\n", buffer, f1, f1-d1, f2-d2, FLT_EPSILON);
  }  

  i1 = strtoint32(buffer, &ep1);
  i2 = strtol(buffer, &ep2, 10);
  //assert(ep1==ep2);
  if (i1!=i2 && (d1 <= INT32_MAX) && (d1 >= INT32_MIN)) {
    printf("%s %d %d\n", buffer, i1, i2);
  }  

  j1 = strtoint64(buffer, &ep1);
  j2 = strtoll(buffer, &ep2, 10);
  //assert(ep1==ep2);
  if (j1!=j2) {
    printf("%s %lld %lld\n", buffer, j1, j2);
  }  

}

int main(void) {
  int64_t i;
  char buffer[1024];
  double d;
  int64_t j;

  for (i=0; i<1e5; i++) {
    d = drand48()*((rand()%2)*2-1);
    snprintf(buffer, 1024, "%.500f", d);
    test_buf(buffer);
    snprintf(buffer, 1024, "%f", d);
    test_buf(buffer);


    d = pow(FLT_MAX, fabs(d))*((rand()%2)*2-1);
    snprintf(buffer, 1024, "%.500e", d);
    test_buf(buffer);

    snprintf(buffer, 1024, "%e", d);
    test_buf(buffer);

    snprintf(buffer, 1024, "%lld", (int64_t)d);
    test_buf(buffer);

    snprintf(buffer, 1024, "%059dfoo", (int32_t)d);
    test_buf(buffer);

    j = rand();
  }
  return 0;
}
