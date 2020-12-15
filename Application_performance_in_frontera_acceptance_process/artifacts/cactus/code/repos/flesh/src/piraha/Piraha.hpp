#ifndef PIRAHA_HPP
#define PIRAHA_HPP
#include <assert.h>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <smart_ptr.hpp>
#include <climits>
#include <string.h>

namespace cctki_piraha {

const int max_int = INT_MAX-1;

using std::map;
using std::vector;

inline char uc_(char a) {
    if(a >= 'a' && a <= 'z')
        return a + 'A' - 'a';
    else
        return a;
}

inline char lc_(char a) {
    if(a >= 'A' && a <= 'Z')
        return a + 'a' - 'A';
    else
        return a;
}

class Group {
public:
    std::string pattern;
    const char *input;
    int start_,end_;
    smart_ptr<vector<smart_ptr<Group> > > children;

    Group(const char *p,const char *value)
        : pattern(p), input(value), start_(0), end_(strlen(value)), children(new vector<smart_ptr<Group> >()) {
        }
    Group(std::string p,const char *input_)
        : pattern(p), input(input_), start_(0), end_(0), children(new vector<smart_ptr<Group> >()) {}
    Group(std::string p,const char *input_,int s,int e,
        smart_ptr<vector<smart_ptr<Group> > > ch)
        : pattern(p), input(input_), start_(s), end_(e), children(ch) {}

    virtual ~Group() {}

    int start() { return start_; }
    int end() { return end_; }
    int childCount(), line();
    std::string getPatternName();
    std::string substring();
    smart_ptr<Group> child(int i);
    void dump(std::ostream& o=std::cout);
    void dump(int n,std::ostream& o,int indent=0);
    void dumpPerl(std::ostream&o=std::cout);
    void dumpPerl(std::ostream&o,int indent);
    void dumpPython(std::ostream&o=std::cout);
    void dumpPython(std::ostream&o,int indent);
    int groupCount() { return children->size(); }
    smart_ptr<Group> group(int i) { return (*children)[i]; }
    smart_ptr<Group> group(const char *nm,int ix=0) {
    	for(unsigned int i=0;i<children->size();i++) {
    		if((*children)[i]->getPatternName() == nm) {
    			if(ix == 0) {
    				return (*children)[i];
    			}
    			ix--;
    		}
    	}
    	smart_ptr<Group> ret;
    	return ret;
    }
};

class Grammar;

class Matcher;

class Pattern {
public:
    virtual bool match(Matcher *m)=0;
    Pattern() {}
    virtual ~Pattern() {}
    virtual std::string fmt() { return "blank"; }
    virtual void insert(std::ostream& o) { o << "{?}"; }
};

inline std::ostream& operator<<(std::ostream& o,Pattern& p) {
    p.insert(o);
    return o;
}

class JMap {
    map<std::string,smart_ptr<Pattern> > m;
public:
    JMap() : m() {}
    smart_ptr<Pattern> get(std::string key) {
        typedef map<std::string,smart_ptr<Pattern> >::iterator mit;
        mit it = m.find(key);
        mit me = m.end();
        if(it == me) {
            return NULL;
        }
        smart_ptr<Pattern> res = m[key];
        assert(res.valid());
        return res;
    }
    void put(std::string key,smart_ptr<Pattern> p) {
        assert(p.valid());
        m[key] = p;
    }
    friend std::ostream& operator<<(std::ostream&,JMap&);
};
inline std::ostream& operator<<(std::ostream& o,JMap& jmap) {
    typedef map<std::string,smart_ptr<Pattern> >::iterator mit;
    mit mb = jmap.m.begin();
    mit me = jmap.m.end();
    o << "{";
	for(mit i = mb; i != me;++i) {
		o << "[" << i->first << "]";
	}
	o << "}";
	return o;
}


class Grammar {
public:
    Grammar() {}
    virtual ~Grammar() {}
    JMap patterns;
    std::string default_rule;
};

class Seq : public Pattern {
    vector<smart_ptr<Pattern> > patterns;
public:
    Seq(Pattern *p,...);
    Seq(vector<smart_ptr<Pattern> > patterns,bool ign,bool show);
    virtual ~Seq() {}
    bool match(Matcher *m);
    virtual void insert(std::ostream& o) {
        for(unsigned int i=0;i<patterns.size();i++)
            o << *patterns[i];
    }
};

class Or : public Pattern {
public:
    vector<smart_ptr<Pattern> > patterns;
    bool ignCase, igcShow;
    Or(bool ign,bool show) : ignCase(ign), igcShow(show) {}
    Or(Pattern *p,...);
    virtual ~Or() {}
    bool match(Matcher *m);
    virtual void insert(std::ostream& o) {
        o << "(";
        for(unsigned int i=0;i<patterns.size();i++) {
            if(i > 0) o << "|";
            o << *patterns[i];
        }
        o << ")";
    }
};

class Literal : public Pattern {
public:
    const char c;
    Literal(char b) : c(b) {}
    bool match(Matcher *m);
    std::string fmt() {
        std::string s = "literal(";
        s += c;
        s += ")";
        return s;
    }
    virtual void insert(std::ostream& o) {
        if(c == '\n')
            o << "\\n";
        else if(c == '\r')
            o << "\\r";
        else if(c == '\t')
            o << "\\t";
        else if(c == '\b')
            o << "\\b";
        else if(c >= 'a' && c <= 'z')
            o << c;
        else if(c >= 'A' && c <= 'Z')
            o << c;
        else if(c >= '0' && c <= '9')
            o << c;
        else
            o << "\\" << c;
    }
};

class ILiteral : public Pattern {
public:
    const char lc,uc;
    ILiteral(char b);
    bool match(Matcher *m);
    std::string fmt() {
        std::string s = "Iliteral(";
        s += lc;
        s += ",";
        s += uc;
        s += ")";
        return s;
    }
    virtual void insert(std::ostream& o) {
        char c = lc;
        if(lc != uc)
            o << "[" << lc << uc << "]";
        else if(c == '\n')
            o << "\\n";
        else if(c == '\r')
            o << "\\r";
        else if(c == '\t')
            o << "\\t";
        else if(c == '\b')
            o << "\\b";
        else if(c >= 'a' && c <= 'z')
            o << c;
        else if(c >= 'A' && c <= 'Z')
            o << c;
        else if(c >= '0' && c <= '9')
            o << c;
        else
            o << "\\" << c;
    }
};

class Lookup : public Pattern {
    smart_ptr<Grammar> gram;
    std::string name;
    bool capture;
public:
    Lookup(std::string s,smart_ptr<Grammar> g);
    virtual ~Lookup() {}
    bool match(Matcher *m);
    std::string fmt() {
        return "Literal:"+name;
    }
    virtual void insert(std::ostream& o) {
        o << "{" << name << "}";
    }
};

class Nothing : public Pattern {
public:
    Nothing() {}
    bool match(Matcher *m) { return true; }
};

class Start : public Pattern {
public:
    Start() {}
    bool match(Matcher *m);
    virtual void insert(std::ostream& o) { o << "^"; }
};

class End : public Pattern {
public:
    End() {}
    bool match(Matcher *m);
    virtual void insert(std::ostream& o) { o << "$"; }
};

class Dot : public Pattern {
public:
    Dot() {}
    bool match(Matcher *m);
    virtual void insert(std::ostream& o) { o << "."; }
};

class Multi : public Pattern {
    const int minv,maxv;
public:
    smart_ptr<Pattern> pattern;
    Multi(int min_,int max_) : minv(min_), maxv(max_), pattern(NULL) {}
    Multi(Pattern *p,int min_,int max_) : minv(min_), maxv(max_), pattern(p) {}
    virtual ~Multi() {}
    bool match(Matcher *m);
    virtual void insert(std::ostream& o) {
        o << *pattern << "{" << minv << "," << maxv << "}";
    }
};

class Range : public Pattern {
public:
    char lo,hi;
    bool match(Matcher *m);
    Range(char lo_,char hi_) : lo(lo_), hi(hi_) {}
};

class Bracket : public Pattern {
public:
    bool neg;
    vector<smart_ptr<Range> > ranges;
    Bracket() : neg(false) {}
    virtual ~Bracket() {}
    Bracket(bool b);
    Bracket *addRange(char lo,char hi);
    Bracket *addRange(char lo,char hi,bool ign);
    bool match(Matcher *m);
    virtual void insert(std::ostream& o);
};

class NegLookAhead : public Pattern {
public:
    smart_ptr<Pattern> pattern;
    NegLookAhead(smart_ptr<Pattern> p) : pattern(p) {}
    virtual ~NegLookAhead() {}
    bool match(Matcher *m);
};

class LookAhead : public Pattern {
public:
    smart_ptr<Pattern> pattern;
    LookAhead(smart_ptr<Pattern> p) : pattern(p) {}
    virtual ~LookAhead() {}
    bool match(Matcher *m) { assert(false); }//TODO: Fill in
};

class Boundary : public Pattern {
    virtual bool match(Matcher *m);
};

class Break : public Pattern {
    virtual bool match(Matcher *m) { assert(false); }//TODO: Fill in
};

class BackRef : public Pattern {
public:
    int index;
    bool ignCase;
    BackRef(int in,bool ign) : index(in), ignCase(ign) {}
    virtual bool match(Matcher *m) { assert(false); }//TODO: Fill in
};

class AutoGrammar {
public:
    static smart_ptr<Grammar> reparserGenerator();
    static smart_ptr<Grammar> fileParserGenerator();
};

class Matcher : public Group {
public:
    Matcher(smart_ptr<Grammar> g,const char *pat_,const char *input_,int input_size=-1);
    virtual ~Matcher() {}

    //std::map<std::string,std::vector<smar_ptr<Group> > > packrat;
    const char *input;
    smart_ptr<Grammar> g;
    int input_size;
    int pos;
    int max_pos;
    int match_to;
    const char *pat;
    bool matches();
    bool matchesTo(int mt);
    Bracket expected;
    void showError(std::ostream& out);
    void showError();
    std::string inrule;
    std::string inrule_max;
    int err_pos;
    void fail(Bracket *ex);
    void fail(char lo,char hi);
    char foundChar() { return input[max_pos+1]; }
};

extern smart_ptr<Grammar> pegGrammar;
extern smart_ptr<Pattern> compile(smart_ptr<Group> g,bool ignCase,smart_ptr<Grammar> gram);
extern void compileFile(smart_ptr<Grammar> g,const char *buffer,signed long buffersize=-1);
void compile(smart_ptr<Grammar> thisg,std::string name,std::string pattern);
void compile(smart_ptr<Grammar> thisg,std::string name,smart_ptr<Group> pattern);
void insertc(std::ostream& o,char c);

}

#endif
