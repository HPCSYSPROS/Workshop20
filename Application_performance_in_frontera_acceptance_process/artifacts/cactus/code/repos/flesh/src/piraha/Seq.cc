#include <iostream>
#include <stdarg.h>
#include "Piraha.hpp"

using namespace cctki_piraha;

Seq::Seq(Pattern *p,...) {
    va_list ap;
    va_start(ap,p);
    patterns.push_back(p);
    while(true) {
        Pattern *pat = va_arg(ap,Pattern*);
        if(pat == NULL)
            break;
        assert(patterns.size()<6);
        patterns.push_back(pat);
    }
    va_end(ap);
}
Seq::Seq(vector<smart_ptr<Pattern> > p,bool ign,bool show) : patterns(p) {}

bool Seq::match(Matcher *m) {
    typedef vector<smart_ptr<Pattern> >::iterator pattern_iter;
    for(pattern_iter p = patterns.begin();p != patterns.end();++p) {
        if(!(*p)->match(m))
            return false;
        if(m->max_pos == m->match_to)
            return true;
    }
    return true;
}
