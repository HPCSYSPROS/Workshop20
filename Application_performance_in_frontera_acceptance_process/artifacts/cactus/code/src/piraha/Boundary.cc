#include "Piraha.hpp"

using namespace cctki_piraha;

bool is_c_ident(char c) {
    return ('a' <= c && c <= 'z')
        || ('A' <= c && c <= 'Z')
        || ('0' <= c && c <= '9')
        || c == '_';
}

bool Boundary::match(Matcher *m) {
    if(m->pos == 0 || m->pos == (int)m->input_size)
        return true;
    char c2 =  m->input[m->pos];
    char c1 =  m->input[m->pos-1];
    bool b1 = is_c_ident(c1);
    bool b2 = is_c_ident(c2);
    return !b1 || !b2;
}
