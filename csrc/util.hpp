#ifndef _UTIL_HPP
#define _UTIL_HPP


#include <curl/curl.h>

class CurlWrapper {
  private:
	  struct MemoryStruct {
	  char *memory;
	  size_t size;
	};
  public:
  static size_t WriteMemoryCallback(void *contents, size_t size, size_t nmemb, void *userp);
	int foo(char *argv[]);
};

#endif //_UTIL_HPP