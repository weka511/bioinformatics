#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util.hpp"





int main(int argc, char *argv[]) {
  if( argc != 2 ) {
    printf("usage: try './curl [url]' to make a get request.\n");
    return 1;
  }
  CurlWrapper wrapper;
  wrapper.foo(argv);
}
