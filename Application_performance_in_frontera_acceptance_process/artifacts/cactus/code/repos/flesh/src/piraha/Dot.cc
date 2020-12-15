#include "Piraha.hpp"

using namespace cctki_piraha;

bool Dot::match(Matcher *m) {
    if(m->pos - m->input_size >= 0)
        return false;
	char c = m->input[m->pos];
	if(c == '\n') {
		return false;
	} else {
		m->pos++;
		return true;
	}
}
