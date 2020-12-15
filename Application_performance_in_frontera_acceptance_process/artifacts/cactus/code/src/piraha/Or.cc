#include <iostream>
#include <stdarg.h>
#include "Piraha.hpp"

using namespace cctki_piraha;

Or::Or(Pattern *p,...) : patterns() {
    va_list ap;
    va_start(ap,p);
    patterns.push_back(p);
    while(true) {
        Pattern *pat = va_arg(ap,Pattern*);
        if(pat == NULL)
            break;
        //std::cout << "pat=" << pat->fmt() << std::endl;
        patterns.push_back(pat);
        assert(patterns.size()<7);
    }
    va_end(ap);
}

bool Or::match(Matcher *m) {
    typedef vector<smart_ptr<Pattern> >::iterator pattern_iter;
    int save = m->pos;
    int chSave = m->children->size();
    for(pattern_iter p = patterns.begin();p != patterns.end();++p) {
        m->pos = save;
        m->children->resize(chSave);
        if((*p)->match(m))
            return true;
    }
    return false;
}
