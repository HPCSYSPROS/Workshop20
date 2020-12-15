#include <iostream>
#include "Piraha.hpp"

using namespace cctki_piraha;

bool Multi::match(Matcher *m) {
    int chSize;
    for(int i=0;i<maxv;i++) {
        int save = m->pos;
        chSize = m->children->size();
        if(!pattern->match(m) || m->pos == save) {
            m->pos = save;
            m->children->resize(chSize);
            return i >= minv;
        }
    }
    return true;
}
