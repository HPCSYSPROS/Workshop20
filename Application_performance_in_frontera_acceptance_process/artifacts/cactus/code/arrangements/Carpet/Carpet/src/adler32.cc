// ADLER-32 checksum algorithm, taken from
// <http://en.wikipedia.org/wiki/Adler-32>

#include "adler32.hh"

namespace Carpet {

unsigned long const mod_adler = 65521;

unsigned long adler32(char const *const data, // location of the data
                      std::size_t const len)  // length of the data in bytes
{
  unsigned long a = 1, b = 0;

  // Process each byte of the data in order
  for (std::size_t i = 0; i < len; ++i) {
    a = (a + (unsigned char)data[i]) % mod_adler;
    b = (b + a) % mod_adler;
  }

  return (b << 16) | a;
}
}
