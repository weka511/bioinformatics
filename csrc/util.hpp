#ifndef _UTIL_HPP
#define _UTIL_HPP

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct MemoryStruct {
  char *memory;
  size_t size;
};

size_t WriteMemoryCallback(void *contents, size_t size, size_t nmemb, void *userp);

#endif //_UTIL_HPP