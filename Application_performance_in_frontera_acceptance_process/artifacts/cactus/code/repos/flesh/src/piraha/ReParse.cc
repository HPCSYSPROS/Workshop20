#include "Piraha.hpp"
#include <stdlib.h>
#include <string.h>

namespace cctki_piraha {

char getChar(smart_ptr<Group> gr) {
    if(gr->groupCount()==1) {
        std::string sub = gr->group(0)->substring();
        int n = 0;
        for(unsigned int i=0;i<sub.size();i++) {
            char c = sub[i];
            if(c >= '0' && c <= '9')
                n = n*16+c-'0';
            else if(c >= 'a' && c <= 'f')
                n = n*16+c-'a'+10;
            else if(c >= 'A' && c <= 'F')
                n = n*16+c-'A'+10;
        }
    }
    std::string gs = gr->substring();
    if(gs.size()==2) {
        char c = gs[1];
        if(c == 'n')
            return '\n';
        else if(c == 'r')
            return '\r';
        else if(c == 't')
            return '\t';
        else if(c == 'b')
            return '\b';
        else
            return c;
    } else {
        return gs[0];
    }
}
smart_ptr<Pattern> mkMulti(smart_ptr<Group> g) {
    if(g->groupCount()==0) {
        std::string s = g->substring();
        if("*" == s) {
            return new Multi(0,max_int);
        } else if("+" == s) {
            return new Multi(1,max_int);
        } else if("?" == s) {
            return new Multi(0,1);
        }
    } else if(g->groupCount()==1) {
        int mn = atol(g->group(0)->substring().c_str());
        return new Multi(mn,mn);
    } else if(g->groupCount()==2) {
        int mn = atol(g->group(0)->substring().c_str());
        if(g->group(1)->groupCount()>0) {
            int mx = atol(g->group(1)->group(0)->substring().c_str());
            return new Multi(mn,mx);
        } else {
            return new Multi(mn,max_int);
        }
    }
    g->dump();
    return NULL;
}


void compileFile(smart_ptr<Grammar> g,const char *buffer,signed long buffersize) {
	if(buffersize < 0)
		buffersize = strlen(buffer);
	smart_ptr<Grammar> grammar = AutoGrammar::fileParserGenerator();
	smart_ptr<Matcher> m = new Matcher(grammar,"file",buffer,buffersize);
	bool b = m->matches();
    if(!b) {
        m->showError();
	    assert(false);
    }

	for(int i=0;i<m->groupCount();i++) {
		smart_ptr<Group> rule = m->group(i);
		smart_ptr<Pattern> ptmp = compile(rule->group(1), false, g);
        std::string nm = rule->group(0)->substring();
		g->patterns.put(nm,ptmp);
        g->default_rule = nm;
	}
}

smart_ptr<Pattern> compile(smart_ptr<Group> g,bool ignCase,smart_ptr<Grammar> gram) {
    std::string pn = g->getPatternName();
    if("literal" == pn) {
        char c = getChar(g);
        if(ignCase)
            return new ILiteral(c);
        else
            return new Literal(c);
    } else if("pattern" == pn) {
        if(g->groupCount()==0)
            return new Nothing();
        return compile(g->group(0),ignCase,gram);
    } else if("pelem" == pn) {
        if(g->groupCount()==2) {
            smart_ptr<Pattern> pm = mkMulti(g->group(1));
            Multi *m = (Multi *)pm.ptr();
            m->pattern = compile(g->group(0),ignCase,gram);
            return pm;
        }
        return compile(g->group(0),ignCase,gram);
    } else if("pelems" == pn||"pelems_top" == pn||"pelems_next" == pn) {
        vector<smart_ptr<Pattern> > li;
        for(int i=0;i<g->groupCount();i++) {
            li.push_back(compile(g->group(i),ignCase,gram));
        }
        if(li.size()==1)
            return li[0];
        return new Seq(li,false,false);
    } else if("group_inside" == pn||"group_top" == pn) {
        if(g->groupCount()==1)
            return compile(g->group(0),ignCase,gram);
        vector<smart_ptr<Pattern> > li;
        for(int i=0;i<g->groupCount();i++) {
            li.push_back(compile(g->group(i),ignCase,gram));
        }
        Or *or_ = new Or(false,false);
        or_->patterns = li;
        smart_ptr<Pattern> orp = or_;
        return orp;
    } else if("group" == pn) {
        Or *or_ = new Or(false,false);
        smart_ptr<Pattern> orp_ = or_;
        bool ignC = ignCase;
        smart_ptr<Group> inside = NULL;
        if(g->groupCount()==2) {
            ignC = or_->igcShow = true;
            std::string ps = g->group(0)->getPatternName();
            if(ps == "ign_on") {
                ignC = or_->ignCase = true;
            } else if(ps == "ign_off") {
                ignC = or_->ignCase = false;
            } else if(ps == "neglookahead") {
                return new NegLookAhead(compile(g->group(1),ignCase,gram));
            } else if(ps == "lookahead") {
                return new LookAhead(compile(g->group(1),ignCase,gram));
            }
            inside = g->group(1);
        } else {
            inside = g->group(0);
        }
        for(int i=0;i<inside->groupCount();i++) {
            or_->patterns.push_back(compile(inside->group(i),ignC,gram));
        }
        if(or_->igcShow == false && or_->patterns.size()==1)
            return or_->patterns[0];
        return orp_;
    } else if("start" == pn) {
        return new Start();
    } else if("end" == pn) {
        return new End();
    } else if("boundary" == pn) {
        return new Boundary();
    } else if("charclass" == pn) {
        Bracket *br = new Bracket();
        smart_ptr<Pattern> brp = br;
        int i=0;
        if(g->groupCount()>0 && g->group(0)->getPatternName() == "neg") {
            i++;
            br->neg = true;
        }
        for(;i < g->groupCount();i++) {
            std::string gn = g->group(i)->getPatternName();
            if("range"==gn) {
                char c0 = getChar(g->group(i)->group(0));
                char c1 = getChar(g->group(i)->group(1));
                br->addRange(c0, c1, ignCase);
            } else {
                char c = getChar(g->group(i));
                br->addRange(c,c, ignCase);
            }
        }
        return brp;
    } else if("named" == pn) {
        std::string lookup = g->group(0)->substring();
        if("brk" == lookup)
            return new Break();
        return new Lookup(lookup, gram);
    } else if("nothing" == pn) {
        return new Nothing();
    } else if("s" == pn||"s0" == pn) {
        return new Lookup("-skipper", gram);
    } else if("dot" == pn) {
        return new Dot();
    } else if("backref" == pn) {
        return new BackRef(g->substring()[1]-'0', ignCase);
    }
    return NULL;
}

}
