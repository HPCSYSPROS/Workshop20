#include "Piraha.hpp"

using namespace cctki_piraha;

ILiteral::ILiteral(char b) : lc(lc_(b)), uc(uc_(b)) {}

bool ILiteral::match(Matcher *m) {
    if(m->pos >= (int)m->input_size )
        return false;
    char c = m->input[m->pos];
    if(c == uc || c == lc) {
        m->max_pos = std::max(m->pos,m->max_pos);
        m->pos++;
        return true;
    } else {
        if(m->pos == m->max_pos+1) {
            Bracket bex;
            bex.addRange(lc,lc);
            bex.addRange(uc,uc);
            m->fail(&bex);
        }
        return false;
    }
}
