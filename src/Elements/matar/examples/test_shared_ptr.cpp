#include <stdio.h>

#include "matar.h"

using namespace mtr; // matar namespace

int N = 2;
using ArrayType = CArray <int>;

ArrayType func()
{
  ArrayType A (N);

  for (int i = 0; i < N; i++) {
    A(i) = 2;
  }
  printf("Pointer of A in func = %p\n", A.pointer());
  for (int i = 0; i < N; i++) {
    printf("Value of A(%d) in func = %d\n", i, A(i));
  }

  return A;
}


int main() {

  auto B = func();

  printf("\n");
  printf("Pointer of B in main = %p\n", B.pointer());
  for (int i = 0; i < N; i++) {
    printf("Value of B(%d) in main = %d\n", i, B(i));
  } 
   
  return 0;
}
