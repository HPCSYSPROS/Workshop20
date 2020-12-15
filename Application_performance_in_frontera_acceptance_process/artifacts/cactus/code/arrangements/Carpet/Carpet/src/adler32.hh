// ADLER-32 checksum algorithm

#ifndef ADLER32_HH
#define ADLER32_HH

#include <cstdlib>

namespace Carpet {

unsigned long adler32(char const *const data, // location of the data
                      std::size_t const len); // length of the data in bytes
}

#endif // #ifndef ADLER32_HH
