#include "Piraha.hpp"
#include <string.h>

using namespace cctki_piraha;

typedef vector<smart_ptr<Range> >::iterator range_iter;

bool Range::match(Matcher *m) {
    if(m->pos - m->input_size >= 0)
        return false;
    char c = m->input[m->pos];
    if(lo <= c && c <= hi) {
        return true;
    } else {
        return false;
    }
}

Bracket::Bracket(bool b) : neg(b), ranges() {}

Bracket *Bracket::addRange(char lo,char hi) {
	bool done = false;
	while(!done) {
		done = true;
	    for(range_iter ri = ranges.begin();ri != ranges.end(); ++ri) {
			smart_ptr<Range> r = *ri;
			if(hi < r->lo || r->hi < lo) {
				// no intersection
				continue;
			} else {
				lo = std::min(lo,r->lo);
				hi = std::max(hi,r->hi);
				ranges.erase(ri);
				done = false;
				break;
			}
		}
	}
    ranges.push_back(new Range(lo,hi));
    return this;
};

Bracket *Bracket::addRange(char lo,char hi,bool ign) {
    if(ign) {
        char lolc = lc_(lo);
        char hilc = lc_(hi);
        char louc = uc_(lo);
        char hiuc = uc_(hi);
        if(lolc == louc && hilc == hiuc) {
            ranges.push_back(new Range(lo,hi));
        } else {
            ranges.push_back(new Range(lolc,hilc));
            ranges.push_back(new Range(louc,hiuc));
        }
    } else {
        ranges.push_back(new Range(lo,hi));
    }
    return this;
};

static void fail(Bracket *b,Matcher *m) {
    typedef vector<smart_ptr<Range> >::iterator range_iter;
    Bracket bex;
    if(m->pos == m->max_pos+1) {
        for(range_iter r = b->ranges.begin();r != b->ranges.end(); ++r) {
        	bex.addRange((*r)->lo,(*r)->hi);
        }
        m->fail(&bex);
    }
}

bool Bracket::match(Matcher *m) {
    if(m->pos >= (int)m->input_size) {
        return false;
    }
    for(range_iter r = ranges.begin();r != ranges.end(); ++r) {
        if((*r)->match(m)) {
            if(neg) {
                fail(this,m);
                return false;
            } else {
                m->max_pos = std::max(m->pos,m->max_pos);
                m->pos++;
                return true;
            }
        }
    }
    if(!neg) {
        fail(this,m);
        return false;
    } else {
        m->pos++;
        m->max_pos = std::max(m->pos,m->max_pos);
        return true;
    }
}

void cctki_piraha::insertc(std::ostream& o,char c) {
    if(c == '-') {
        o << "\\-";
    } else if(c == '\n') {
        o << "\\n";
    } else if(c == '\r') {
        o << "\\r";
    } else if(c == '\t') {
        o << "\\t";
    } else if(c == '\b') {
        o << "\\b";
    } else if(strchr("\\\"[]-",c) != 0) {
        o << "\\" << c;
    } else {
        o << c;
    }
}

void Bracket::insert(std::ostream& o) {
    o << "[";
    if(neg)
        o << "^";
    for(range_iter r = ranges.begin();r != ranges.end(); ++r) {
    	char lo = (*r)->lo, hi = (*r)->hi;
    	if(lo == hi) {
            insertc(o,lo);
    	} else {
            insertc(o,lo);
            o << '-';
            insertc(o,hi);
    	}
    }
    o << "]";
}
