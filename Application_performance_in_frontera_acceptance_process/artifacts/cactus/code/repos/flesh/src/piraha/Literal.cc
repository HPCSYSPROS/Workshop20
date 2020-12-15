#include "Piraha.hpp"

using namespace cctki_piraha;

bool Literal::match(Matcher *m) {
    if(m->pos - m->input_size >= 0) {
        return false;
    }
    if(m->input[m->pos] == c) {
        m->max_pos = std::max(m->pos,m->max_pos);
        m->pos++;
        return true;
    } else {
        if(m->pos == m->max_pos+1) {
        	Bracket bex;
        	bex.addRange(c,c);
            m->fail(&bex);
        }
        return false;
    }
}
