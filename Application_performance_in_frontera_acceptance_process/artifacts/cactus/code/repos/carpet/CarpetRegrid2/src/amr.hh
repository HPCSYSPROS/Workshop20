#ifndef AMR_HH
#define AMR_HH

#include <carpet.hh>

namespace CarpetRegrid2 {

void evaluate_level_mask(cGH const *restrict cctkGH, vector<ibset> &regions,
                         int rl);

} // namespace CarpetRegrid2

#endif // #ifndef AMR_HH
