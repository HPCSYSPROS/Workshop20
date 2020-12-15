#include "Piraha.hpp"

using namespace cctki_piraha;

bool Start::match(Matcher *m) {
    return m->pos == 0;
}
