#include <iostream>
#include <iomanip>
#include <string>
#include "Piraha.hpp"
#include <cstdlib>
#include <sstream>
#include "cctk_core.h"
#include "cctk_Parameter.h"
#include <map>
#include <cmath>
#include <algorithm>
#include <limits>
#include <fstream>
#include "util_Expression.h"

namespace cctki_piraha {

#define VAR(X) " " #X "=" << X

extern "C" int CCTK_ParameterFilename(int len, char *filename);

smart_ptr<Grammar> create_grammar() {
    smart_ptr<Grammar> grammar = new Grammar();
    const char *par_file_src =
        "skipper = ([ \\t\\r\\n]|\\#.*)*\n"
        "# comment\n"
        "skipeol = ([ \\t\\r]|\\#.*)*($|\\n)\n"
        "any = [^]\n"
        "stringcomment = #.*\n"
        "stringparser = ^({stringcomment}|{var}|{name}|{any})*$\n"

        "# Note that / occurs in some par files. It is my\n"
        "# feeling that this should require quote marks.\n"

        "name = [a-zA-Z][a-zA-Z0-9_]*\n"
        "dname = [0-9][a-zA-Z_]{2,}\n"
        "inquot = ({var}|\\\\.|[^\\\\\"])*\n"
        "fname = \\.?/[-\\./0-9a-zA-Z_]+\n"
        "quot = \"{inquot}\"|{fname}\n"
        "num = ([0-9]+(\\.[0-9]*|)|\\.[0-9]+)([edDE][+-]?[0-9]+|)\n"
        "env = ENV\\{{name}\\}\n"
        "var = \\$({env}|{name}|\\{{name}\\})\n"

        "powexpr = {value}( \\*\\* {value})?\n"
        "mulop = [*/%]\n"
        "mexpr = {powexpr}( {mulop} {powexpr})*\n"
        "addop = [+-]\n"
        "aexpr = {mexpr}( {addop} {mexpr})*\n"
        "compop = [<>]=?\n"
        "compexpr = {aexpr}( {compop} {aexpr})?\n"
        "eqop = [!=]=\n"
        "eqexpr = {compexpr}( {eqop} {eqexpr})?\n"
        "andexpr = {eqexpr}( && {eqexpr})?\n"
        "expr = {andexpr}( \\|\\| {andexpr})?\n"
        "eval = {expr}\n"

        "paren = \\( {expr} \\)\n"
        "par = {name} :: {name}( {parindex})?\n"
        "func = {name} \\( {expr} \\)\n"
        "array = \\[ {expr}( , {expr})* \\]\n"

        "value = {unop}?({par}|{func}|{paren}|{dname}|{num}|{quot}|{name}|{var})\n"
        "unop = [-!]\n"

        "int = [0-9]+\n"
        "index = \\[ {int} \\]\n"
        "parindex = \\[ {expr} \\]\n"
        "active = (?i:ActiveThorns)\n"
        "set = ({active} = ({quot}|{name})|{par}( {index}|) = ({array}|\\+?{expr})){-skipeol}\n"
        "set_var = \\${name} = \\+?{expr}{-skipeol}\n"
        "desc = !DESC {quot}\n"
        "file = ^( ({desc}|{set_var}|{set}|{active}) )*$";
    //std::ofstream peg("/tmp/par.peg");
    //peg << par_file_src;
    //peg.close();

    compileFile(grammar,par_file_src,strlen(par_file_src));
    return grammar;
}

/**
 * This holds the structure required for parsing
 * a Cactus par file.
 */
static smart_ptr<Grammar> par_file_grammar = create_grammar();

static std::string mklower(std::string& in) {
    std::string s = in;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

enum ValueType { PIR_STRING,PIR_INT,PIR_REAL,PIR_BOOL,PIR_VOID };

std::string get_parfile() {
    char path[500];
    CCTK_ParameterFilename(500, path);
    char *value = strrchr (path, '/');
    if (value == NULL) {
        value = path;
    } else {
        value++;
    }
    std::string s(value);
    return s;
}
std::string get_parfilename() {
    std::string ending = ".par";
    std::string s = get_parfile();
    size_t ns = s.length();
    size_t ne = ending.length();
    if(ns > ne && s.substr(ns-ne).compare(ending)==0) {
        s = s.substr(0,ns-ne);
    }
    return s;
}

std::ostream& operator<<(std::ostream& o,const ValueType& vt) {
    if(vt == PIR_STRING)
        o << "STRING";
    else if(vt == PIR_INT)
        o << "INT";
    else if(vt == PIR_REAL)
        o << "REAL";
    else if(vt == PIR_BOOL)
        o << "BOOL";
    else if(vt == PIR_VOID)
        o << "VOID";
    else
        o << "UNDEF";
    return o;
}

std::string current_thorn;

/**
 * The Value class is something like a union. It holds
 * integer, real, and string data and has a type field
 * that identifies which field is currently valid.
 */
struct Value {
    /** This field holds the parse tree element associated with this Value. */
    smart_ptr<Group> hold;
    CCTK_REAL ddata;
    CCTK_INT idata;
    std::string sdata;
    ValueType type;
    Value(smart_ptr<Group> g) : hold(g), ddata(0), idata(0), sdata(), type(PIR_VOID) { smart_ptr<Group> foo(g); }
    ~Value() {}
    /**
     * Create a string representation of the Value.
     */
    std::string copy() {
        assert(type != PIR_VOID);
        if(type == PIR_STRING) {
            return sdata;
        } else if(type == PIR_BOOL) {
            std::string s = idata ? "yes" : "no";
            return s;
        } else {
            std::ostringstream o;
            if(type == PIR_REAL) {
                o << std::setprecision(15) << ddata;
            } else {
                o << idata;
            }
            return o.str();
        }
    }
    /**
     * Check to see if something is a bool and throw an error
     * if it's not.
     */
    void checkBool() {
        if(type != PIR_BOOL) {
            std::ostringstream msg;
            msg << "Does not evaluate to a boolean: " << hold->substring();
            std::string par = get_parfile();
            CCTK_Error(hold->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
        }
    }
    /**
     * Compares with another Value, and throws a CCTK_Error
     * if the comparison doesn't make sense.
     */
    bool equals(smart_ptr<Value> v) {
        if(type == PIR_BOOL && v->type == PIR_BOOL) {
            return idata == v->idata;
        } else if(type == PIR_INT && v->type == PIR_INT) {
            return idata == v->idata;
        } else if(type == PIR_INT && v->type == PIR_REAL) {
            return idata == v->ddata;
        } else if(type == PIR_REAL && v->type == PIR_INT) {
            return ddata == v->idata;
        } else if(type == PIR_REAL && v->type == PIR_REAL) {
            return ddata == v->ddata;
        } else if(type == PIR_STRING && v->type == PIR_STRING) {
            return sdata == v->sdata;
        }
        std::ostringstream msg;
        msg << "Cannot compare " << type << " and " << v->type << std::endl;
        std::string par = get_parfile();
        CCTK_Error(hold->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
    }
    /**
     * Return a real value, whether the underlying
     * quantity is integer or real.
     */
    CCTK_REAL realValue() {
        if(type == PIR_REAL)
            return ddata;
        else if(type == PIR_INT)
            return idata;
        std::ostringstream msg;
        msg << "Cannot convert " << type << " to floating point." << std::endl;
        std::string par = get_parfile();
        CCTK_Error(hold->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
        return 0;
    }
    bool intOrreal() {
        return type == PIR_INT || type == PIR_REAL;
    }
    /**
     * This function converts a real to a real, but
     * only if this can be done without loss of precision.
     */
    void integerize() {
        if(type == PIR_REAL) {
            idata = std::lrint(ddata);
            if(idata == ddata) {
                type = PIR_INT;
            }
        }
    }
    void booleanize(smart_ptr<Group> gr) {
        if(type == PIR_STRING) {
            std::string s = mklower(sdata);
            if(s == "yes" || s == "true") {
                idata = 1;
                type = PIR_BOOL;
            } else if(s == "no" || s == "false") {
                idata = 0;
                type = PIR_BOOL;
            }
        } else if(type == PIR_INT) {
            /// Steven R. Brandt would like to remove this
            /// particular auto-conversion
            if(idata == 1 || idata == 0) {
                type = PIR_BOOL;
                std::ostringstream msg;
                msg << "Boolean variable is set with integer: " << gr->substring() << std::endl;
                std::string par = get_parfile();
                CCTK_Warn(1,gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
            }
        }
    }
};

std::ostream& operator<<(std::ostream& o,const smart_ptr<Value>& val) {
    if(!val.valid()) {
        return o << "NULL";
    }
    o << "(" << val->type << ")";
    if(val->type == PIR_REAL) {
        o << val->ddata;
    } else if(val->type == PIR_INT) {
        o << val->idata;
    } else if(val->type == PIR_BOOL) {
        o << (val->idata ? "true" : "false");
    } else if(val->type == PIR_STRING) {
        o << "\"" << val->sdata << "\"";
    } else {
        o << val->hold->substring();
    }
    return o;
}

typedef std::map<std::string,std::map<std::string,smart_ptr<Value> > >::iterator th_iter;
typedef std::map<std::string,smart_ptr<Value> >::iterator nm_iter;
std::map<std::string,smart_ptr<Value> > variables;

smart_ptr<Value> eval_expr(std::string inp);

enum read_write { init, read, write };

typedef std::map<std::string,read_write> read_write_tracker_type;
read_write_tracker_type read_write_tracker;

read_write& read_write_status(std::string& key) {
  if(read_write_tracker.find(key) == read_write_tracker.end())
    read_write_tracker[key] = init;
  return read_write_tracker[key];
}

read_write& read_write_status(std::string& thorn,std::string& name) {
  std::ostringstream msg;
  msg << thorn << "::" << name;
  std::string key = msg.str();
  return read_write_status(key);
}

// If the value was already defined in this parameter file, look
// it up in the map. Otherwise, get it from Cactus.
smart_ptr<Value> find_val(smart_ptr<Group> gr,std::string thorn,std::string name) {
    smart_ptr<Value> ret = new Value(gr);
    int type;
    const void *result;
    result = CCTK_ParameterGet(name.c_str(),thorn.c_str(),&type);
    if(result == NULL) {
        std::ostringstream msg;
        msg << "Undefined or inaccessible variable: " << thorn << "::" << name << std::endl;
        std::string par = get_parfile();
        CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
    }
    read_write& status = read_write_status(thorn,name);
    if(status == init)
      status = read;

    switch(type) {
        case PARAMETER_REAL:
            ret->type = PIR_REAL;
            ret->ddata = *(const CCTK_REAL*)result;
            break;
        case PARAMETER_INT:
            ret->type = PIR_INT;
            ret->idata = *(const CCTK_INT*)result;
            break;
        case PARAMETER_BOOLEAN:
            ret->type = PIR_BOOL;
            ret->idata = *(const CCTK_INT*)result;
            break;
        case PARAMETER_STRING:
        case PARAMETER_SENTENCE:
        case PARAMETER_KEYWORD:
            {
                ret->type = PIR_STRING;
                const char *s = *(const char *const *)result;
                ret->sdata += s;
            }
            break;
        default:
            std::ostringstream msg;
            msg << "Unexpected type result from ParameterGet=" << type << std::endl;
            std::string par = get_parfile();
            CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
    }
    //std::cout << "Piraha: GET: " << thorn << "::" << name << "=" << ret << std::endl;
    return ret;
}

smart_ptr<Value> lookup_var(smart_ptr<Group> gr) {
    smart_ptr<Value> ret;
    if(gr->group(0)->getPatternName() == "env") {
        const char *env = getenv(gr->group(0)->group(0)->substring().c_str());
        if(env != NULL) {
            ret = new Value(gr);
            ret->type = PIR_STRING;
            ret->sdata = env;
        }
    } else if(gr->group(0)->substring() == "parfile") {
        ret = new Value(gr);
        ret->type = PIR_STRING;
        ret->sdata = get_parfilename();
    } else if(gr->group(0)->substring() == "pi") {
        ret = new Value(gr);
        ret->type = PIR_REAL;
        ret->ddata = 4.0*atan2(1.0,1.0);
    }
    nm_iter iter = variables.find(gr->group(0)->substring());
    if(iter != variables.end()) {
        ret = iter->second;
    }
    if(!ret.valid()) {
        std::ostringstream msg;
        msg << "Unknown variable: " << gr->substring() << std::endl;
        std::string par = get_parfile();
        CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
    }
    return ret;
}


std::string string_reparser(std::string s) {
    smart_ptr<Matcher> m = new Matcher(par_file_grammar,"stringparser",s.c_str(),s.length());
    if(m->matches()) {
        std::string out = "";
        for(int i=0;i < m->groupCount(); i++) {
            std::string pn = m->group(i)->getPatternName();
            if(pn == "any" || pn == "name") {
                out += m->group(i)->substring();
            } else if(pn == "stringcomment") {
                ;
            } else {
                smart_ptr<Value> val = lookup_var(m->group(i));
                out += val->copy();
            }
        }
        return out;
    } else {
        return s;
    }
}

struct ExpressionEvaluationData {
  std::map<std::string,uExpressionValue> values;
  const void *data;
  int (*eval)(int, const char * const *, uExpressionValue *, const void *);
  ExpressionEvaluationData() : values(), data(0), eval(0) {}
  ~ExpressionEvaluationData() {}
};

/**
 * The meval() function takes any node within
 * the parse tree and creates a Value object
 * from it. It is, therefore, designed to be
 * used recursively.
 **/
smart_ptr<Value> meval(smart_ptr<Group> gr,ExpressionEvaluationData *eedata) {
    assert(gr.valid());
    std::string pn = gr->getPatternName();
    smart_ptr<Value> ret = new Value(gr);
    if(pn == "num") {
        std::string s = gr->substring();
        s = mklower(s);
        std::replace(s.begin(),s.end(),'d','e');
        std::istringstream iss(s);
        if(iss >> ret->ddata) {
            ret->ddata = atof(s.c_str());
            ret->idata = std::lrint(ret->ddata);
            if(ret->idata == ret->ddata && (s.find('.') == std::string::npos))
                ret->type = PIR_INT;
            else
                ret->type = PIR_REAL;
        } else {
            std::ostringstream msg;
            msg << "Invalid numerical value: " << s << std::endl;
            std::string par = get_parfile();
            CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
        }
    } else if(pn == "paren" || pn == "parindex") {
        return meval(gr->group(0),eedata);
    } else if(pn == "func") {
        std::string fn = gr->group(0)->substring();
        fn = mklower(fn);
        smart_ptr<Value> val = meval(gr->group(1),eedata);
        if(val->type == PIR_REAL || val->type == PIR_INT) {
            if(fn == "trunc") {
                val->ddata = trunc(val->realValue());
                val->type = PIR_REAL;
                return val;
            } else if(fn == "floor") {
                val->ddata = floor(val->realValue());
                val->type = PIR_REAL;
                return val;
            } else if(fn == "ceil") {
                val->ddata = ceil(val->realValue());
                val->type = PIR_REAL;
                return val;
            } else if(fn == "sqrt") {
                val->ddata = sqrt(val->realValue());
                val->type = PIR_REAL;
                return val;
            } else if(fn == "atan") {
                val->ddata = atan(val->realValue());
                val->type = PIR_REAL;
                return val;
            } else if(fn == "sin") {
                val->ddata = sin(val->realValue());
                val->type = PIR_REAL;
                return val;
            } else if(fn == "cos") {
                val->ddata = cos(val->realValue());
                val->type = PIR_REAL;
                return val;
            } else if(fn == "tan") {
                val->ddata = tan(val->realValue());
                val->type = PIR_REAL;
                return val;
            } else if(fn == "exp") {
                val->ddata = exp(val->realValue());
                val->type = PIR_REAL;
                return val;
            } else if(fn == "log") {
                val->ddata = log(val->realValue());
                val->type = PIR_REAL;
                return val;
            } else if(fn == "abs") {
                val->ddata = fabs(val->realValue());
                val->type = PIR_REAL;
                return val;
            } else if(fn == "acos") {
                val->ddata = acos(val->realValue());
                val->type = PIR_REAL;
                return val;
            } else if(fn == "asin") {
                val->ddata = asin(val->realValue());
                val->type = PIR_REAL;
                return val;
            } else if(fn == "bool") {
                if(val->type == PIR_REAL) {
                    val->idata = std::lrint(val->ddata);
                }
                val->type = PIR_BOOL;
                return val;
            } else if(fn == "int") {
                if(val->type == PIR_REAL) {
                    val->idata = std::lrint(val->ddata);
                    val->type = PIR_INT;
                }
                return val;
            } else if(fn == "real") {
                if(val->type == PIR_INT) {
                    val->ddata = val->idata;
                    val->type = PIR_REAL;
                }
                return val;
            }
        } else if(val->type == PIR_BOOL) {
            if(fn == "int") {
                val->type = PIR_INT;
                return val;
            } else if(fn == "real") {
                val->type = PIR_REAL;
                val->ddata = val->idata;
                return val;
            }
        } else if(val->type == PIR_STRING) {
            if(fn == "int") {
                val->type = PIR_INT;
                std::istringstream buf(val->sdata);
                if(!(buf >> val->idata)) {
                    std::ostringstream msg;
                    msg << "Invalid numerical value: " << val->sdata << std::endl;
                    std::string par = get_parfile();
                    CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
                }
                std::string extra;
                if(buf >> extra) {
                    std::ostringstream msg;
                    msg << "Trailing input in numerical value: " << val->sdata << std::endl;
                    std::string par = get_parfile();
                    CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
                }
                return val;
            } else if(fn == "real") {
                val->type = PIR_REAL;
                std::istringstream buf(val->sdata);
                if(!(buf >> val->ddata)) {
                    std::ostringstream msg;
                    msg << "Invalid numerical value: " << val->sdata << std::endl;
                    std::string par = get_parfile();
                    CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
                }
                std::string extra;
                if(buf >> extra) {
                    std::ostringstream msg;
                    msg << "Trailing input in numerical value: " << val->sdata << std::endl;
                    std::string par = get_parfile();
                    CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
                }
                return val;
            } else if(fn == "bool") {
                val->type = PIR_BOOL;
                std::string s = mklower(val->sdata);
                if(s == "no" || s == "false") {
                    ret->idata = 0;
                } else if(s == "yes" || s == "true") {
                    ret->idata = 1;
                } else {
                    std::ostringstream msg;
                    msg << "Invalid boolean value: " << val->sdata << std::endl;
                    std::string par = get_parfile();
                    CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
                }
                return val;
            }
        }
        std::ostringstream msg;
        msg << "Unknown func: " << fn << "(" << val->type << ")" << std::endl;
        std::string par = get_parfile();
        CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
    } else if(pn == "name"||pn == "dname") {
        std::string s = gr->substring();
        s = mklower(s);
        if(s == "nan") {
            ret->ddata = NAN;
            ret->type = PIR_REAL;
        } else if(s == "inf") {
            ret->ddata = INFINITY;
            ret->type = PIR_REAL;
        } else if(s == "no" || s == "false") {
            ret->type = PIR_BOOL;
            ret->idata = 0;
        } else if(s == "yes" || s == "true") {
            ret->type = PIR_BOOL;
            ret->idata = 1;
        } else {
          ret->type = PIR_STRING;
          ret->sdata = gr->substring();
          bool evaluated = false;
          // vname_string needs to be an actual variable and not just a
          // temporary so that c_str() stays valid
          const std::string vname_string = gr->substring();
          const char *vname = vname_string.c_str();
          if(eedata != 0) {
            uExpressionValue uval;
            std::map<std::string,uExpressionValue>::iterator iter = eedata->values.find(vname);
            if(iter != eedata->values.end()) {
              evaluated = true;
              uval = iter->second;
            } else if(eedata->eval(1,&vname,&uval,eedata->data) == 0) {
              evaluated = true;
            }
            if(evaluated) {
              if(uval.type == uExpressionValue::ival) {
                ret->type = PIR_INT;
                ret->idata = uval.value.ival;
              } else if(uval.type == uExpressionValue::rval) {
                ret->type = PIR_REAL;
                ret->ddata = uval.value.rval;
              } else if(uval.type == uExpressionValue::sval) {
                ret->type = PIR_STRING;
                ret->sdata = uval.value.sval;
                delete uval.value.sval;
              }
            }
          }
        }
        return ret;
    } else if(pn == "par") {
        std::string thorn = gr->group(0)->substring();
        std::string name = gr->group(1)->substring();
        if(gr->groupCount() == 3) {
            std::ostringstream vn;
            smart_ptr<Value> index = meval(gr->group(2),eedata);
            if(index->type == PIR_INT) {
                std::stringstream o;
                o << name << "[" << index->idata << "]";
                o << std::flush;
                std::string keyi = o.str();
                ret = find_val(gr,thorn,keyi);
                return ret;
            }
        } else {
            ret = find_val(gr,thorn,name);
            return ret;
        }
        std::ostringstream msg;
        msg << "Unknown par: " << thorn << "::" << name << std::endl;
        std::string par = get_parfile();
        CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
    } else if(pn == "var") {
        ret = lookup_var(gr);
    } else if(pn == "value") {
        if(gr->groupCount()==2) {
            std::string unop = gr->group(0)->substring();
            ret = meval(gr->group(1),eedata);
            if(unop == "-") {
                if(ret->type == PIR_INT) {
                    ret->idata = -ret->idata;
                } else if(ret->type == PIR_REAL) {
                    ret->ddata = std::copysign(ret->ddata,
                                               std::signbit(ret->ddata) ?
                                               1. : -1.);
                } else {
                    std::ostringstream msg;
                    msg << "Unknown operation: " << unop << ret->type << std::endl;
                    std::string par = get_parfile();
                    CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
                }
            } else if(unop == "!") {
                if(ret->type == PIR_BOOL) {
                    ret->idata = !ret->idata;
                } else {
                    std::ostringstream msg;
                    msg << "Unknown operation: " << unop << ret->type << std::endl;
                    std::string par = get_parfile();
                    CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
                }
            } else {
                std::ostringstream msg;
                msg << "Unknown operation: " << unop << ret->type << std::endl;
                std::string par = get_parfile();
                CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
            }
        } else {
            return meval(gr->group(0),eedata);
        }
    } else if(pn == "quot") {
        ret->type = PIR_STRING;
        ret->sdata = string_reparser(gr->group(0)->substring());
    } else if(pn == "inquot") {
        ret->type = PIR_STRING;
        ret->sdata = gr->substring();
    } else if(pn == "expr") {
        if(gr->groupCount()==1)
            return meval(gr->group(0),eedata);
        smart_ptr<Value> v1 = meval(gr->group(0),eedata);
        smart_ptr<Value> v2 = meval(gr->group(1),eedata);
        v1->checkBool();
        v2->checkBool();
        ret->type = PIR_BOOL;
        ret->idata = v1->idata || v2->idata;
    } else if(pn == "powexpr") {
        if(gr->groupCount()==1)
            return meval(gr->group(0),eedata);
        smart_ptr<Value> v1 = meval(gr->group(0),eedata);
        smart_ptr<Value> v2 = meval(gr->group(1),eedata);
        ret->type = PIR_REAL;
        ret->ddata = pow(v1->realValue(),v2->realValue());
    } else if(pn == "andexpr") {
        if(gr->groupCount()==1)
            return meval(gr->group(0),eedata);
        smart_ptr<Value> v1 = meval(gr->group(0),eedata);
        smart_ptr<Value> v2 = meval(gr->group(1),eedata);
        v1->checkBool();
        v2->checkBool();
        ret->type = PIR_BOOL;
        ret->idata = v1->idata && v2->idata;
    } else if(pn == "eqexpr") {
        if(gr->groupCount()==1)
            return meval(gr->group(0),eedata);
        smart_ptr<Value> v1 = meval(gr->group(0),eedata);
        std::string eqop = gr->group(1)->substring();
        smart_ptr<Value> v2 = meval(gr->group(2),eedata);
        ret->type = PIR_BOOL;
        if(eqop == "==") {
            ret->idata = v1->equals(v2);
        } else if(eqop == "!=") {
            ret->idata = !v1->equals(v2);
        } else {
            std::ostringstream msg;
            msg << "Unknown equality operator: " << eqop << std::endl;
            std::string par = get_parfile();
            CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
        }
    } else if(pn == "compexpr") {
        if(gr->groupCount()==1)
            return meval(gr->group(0),eedata);
        smart_ptr<Value> v1 = meval(gr->group(0),eedata);
        if(gr->groupCount()>0) {
            std::string compop = gr->group(1)->substring();
            smart_ptr<Value> v2 = meval(gr->group(2),eedata);
            CCTK_REAL d1 = v1->realValue();
            CCTK_REAL d2 = v2->realValue();
            ret->type = PIR_BOOL;
            if(compop == "<") {
                ret->idata = (d1 < d2);
                return ret;
            } else if(compop == ">") {
                ret->idata = (d1 > d2);
                return ret;
            } else if(compop == "<=") {
                ret->idata = (d1 <= d2);
                return ret;
            } else if(compop == ">=") {
                ret->idata = (d1 >= d2);
                return ret;
            }
            std::ostringstream msg;
            msg << "Unknown comparison operator: " << compop << std::endl;
            std::string par = get_parfile();
            CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
        }
    } else if(pn == "aexpr") {
        if(gr->groupCount()==1)
            return meval(gr->group(0),eedata);
        smart_ptr<Value> v1 = meval(gr->group(0),eedata);
        for(int i=1;i+1<gr->groupCount();i+=2) {
            std::string addop = gr->group(i)->substring();
            smart_ptr<Value> v2 = meval(gr->group(i+1),eedata);
            assert(v2.valid());
            if(v1->type == PIR_INT && v2->type == PIR_INT) {
                ret->type = PIR_INT;
                if(addop == "+") {
                    ret->idata = v1->idata + v2->idata;
                } else if(addop == "-") {
                    ret->idata = v1->idata - v2->idata;
                } else {
                    std::ostringstream msg;
                    msg << "Unknown add operator: " << addop << std::endl;
                    std::string par = get_parfile();
                    CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
                }
            } else if(v1->type == PIR_STRING && v2->type == PIR_STRING) {
                ret->type = PIR_STRING;
                if(addop == "+") {
                    ret->sdata = v1->sdata + v2->sdata;
                } else {
                    std::ostringstream msg;
                    ret->sdata = v1->sdata + addop + v2->sdata;
                    msg << "Unknown add operator: " << addop << std::endl;
                    msg << "Interpreting as literal string with value '" << ret->sdata << "'" << std::endl;
                    std::string par = get_parfile();
                    CCTK_Warn(1,gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
                }
            } else if(v1->intOrreal() && v2->intOrreal()) {
                ret->type = PIR_REAL;
                if(addop == "+") {
                    ret->ddata = v1->realValue()+v2->realValue();
                } else if(addop == "-") {
                    ret->ddata = v1->realValue()-v2->realValue();
                } else {
                    std::ostringstream msg;
                    msg << "Unknown add operator: " << addop << std::endl;
                    std::string par = get_parfile();
                    CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
                }
            } else {
                std::ostringstream msg;
                msg << "Unknown operation: " << v1->type;
                msg << " " << addop << v2->type << std::endl;
                std::string par = get_parfile();
                CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
            }
            v1 = ret;
        }
    } else if(pn == "mexpr") {
        if(gr->groupCount()==1)
            return meval(gr->group(0),eedata);
        smart_ptr<Value> v1 = meval(gr->group(0),eedata);
        for(int i=1;i+1<gr->groupCount();i+=2) {
            std::string mulop = gr->group(i)->substring();
            smart_ptr<Value> v2 = meval(gr->group(i+1),eedata);
            if(v1->type == PIR_INT && v2->type == PIR_INT) {
                ret->type = PIR_INT;
                if(mulop == "*") {
                    ret->idata = v1->idata * v2->idata;
                } else if(mulop == "/") {
                    ret->idata = v1->idata / v2->idata;
                } else if(mulop == "%") {
                    ret->idata = v1->idata % v2->idata;
                } else {
                    std::ostringstream msg;
                    msg << "Unknown mul operator: " << v1->type << mulop << v2->type << std::endl;
                    std::string par = get_parfile();
                    CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
                }
            } else if(v1->type == PIR_STRING && v2->type == PIR_INT) {
                ret->type = PIR_STRING;
                if(mulop == "*") {
                    ret->sdata = "";
                    for(CCTK_INT i=0;i<v2->idata;i++)
                        ret->sdata += v1->sdata;
                } else {
                    std::ostringstream msg;
                    msg << "Unknown mul operator: " << v1->type << mulop << v2->type << std::endl;
                    std::string par = get_parfile();
                    CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
                }
            } else if(v1->intOrreal() && v2->intOrreal()) {
                ret->type = PIR_REAL;
                if(mulop == "*") {
                    ret->ddata = v1->realValue()*v2->realValue();
                } else if(mulop == "/") {
                    ret->ddata = v1->realValue()/v2->realValue();
                } else {
                    std::ostringstream msg;
                    msg << "Unknown mul operator: " << v1->type << mulop << v2->type << std::endl;
                    std::string par = get_parfile();
                    CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
                }
            } else if(v1->type == PIR_STRING && v2->type == PIR_STRING) {
                std::ostringstream msg;
                ret->type = PIR_STRING;
                ret->sdata = v1->sdata + mulop + v2->sdata;
                msg << "Unknown operation: " << v1->type << " " << mulop << v2->type << std::endl;
                msg << "Interpreting as literal string with value '" << ret->sdata << "'" << std::endl;
                std::string par = get_parfile();
                CCTK_Warn(1,gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
            } else {
                std::ostringstream msg;
                msg << "Unknown operation: " << v1->type << " " << mulop << v2->type << std::endl;
                std::string par = get_parfile();
                CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
            }
            v1 = ret;
        }
    } else {
        std::ostringstream msg;
        std::string par = get_parfile();
        msg << "Pattern not handled[" << gr->getPatternName() << "]=" << gr->substring() <<
            " at " << gr->line() << " in " << par << std::endl;
        CCTK_Error(__LINE__,__FILE__,"piraha",msg.str().c_str());
    }
    return ret;
}


smart_ptr<Value> eval_expr(std::string input,ExpressionEvaluationData *eedata) {
    smart_ptr<Matcher> m = new Matcher(par_file_grammar,"eval",input.c_str());
    smart_ptr<Value> ret;
    if(m->matches()) {
        ret = meval(m->group(0),eedata);
    }
    return ret;
}

void check_types(const char *thorn,int line,smart_ptr<Value> svm,int t) {
    ValueType v = svm->type;
    bool ok = false;
    if(v == PIR_STRING && (t == PARAMETER_STRING || t == PARAMETER_KEYWORD))
        ok = true;
    else if(v == PIR_INT && t == PARAMETER_INT)
        ok = true;
    else if(v == PIR_BOOL && t == PARAMETER_BOOLEAN)
        ok = true;
    else if((v == PIR_INT || v == PIR_REAL) && t == PARAMETER_REAL)
        ok = true;
    if(!ok) {
        std::string par = get_parfile();
        std::ostringstream msg;
        msg << "Invalid assignment: Attempting to set a variable of type ";
        switch(t) {
            case PARAMETER_BOOLEAN:
                msg << "BOOL";
                break;
            case PARAMETER_INT:
                msg << "INT";
                break;
            case PARAMETER_KEYWORD:
                msg << "KEYWORD";
                break;
            case PARAMETER_STRING:
                msg << "STRING";
                break;
            case PARAMETER_REAL:
                msg << "REAL";
                break;
            default:
                msg << "type(" << t << ")";
                break;
        }
        msg << " with " << svm << std::endl;
        CCTK_Error(line,par.c_str(),thorn,msg.str().c_str());
    }
}

extern "C" void *Util_ExpressionParse(const char *expr) {
    int exprsize = strlen(expr);
    Matcher *m2 = new Matcher(par_file_grammar,"eval",expr,exprsize);
    bool b = m2->matches();
    if(b) {
      return m2;
    } else {
      std::ostringstream msg;
      m2->showError(msg);
      CCTK_Error(__LINE__,__FILE__,"Piraha",msg.str().c_str());
      return 0;
    }
}

extern "C" void Util_ExpressionFree(void *m2_) {
  Matcher *m2 = (Matcher *)m2_;
  delete m2;
}

extern "C" int Util_ExpressionEvaluate(void *m2_,
      uExpressionValue *result,
      uExpressionEvaluator eval,
      const void *data) {
  Matcher *m2 = (Matcher *)m2_;

  ExpressionEvaluationData eedata;
  eedata.eval = eval;
  eedata.data = data;

  smart_ptr<Value> value;
  bool set;
  if(m2->group(0)->getPatternName() == "set_var") {
    set = true;
    value = meval(m2->group(0)->group(1),&eedata);
  } else {
    set = false;
    //m2->dump(std::cout);
    value = meval(m2->group(0),&eedata);
  }
  if(value->type == PIR_INT) {
    result->type = uExpressionValue::ival;
    result->value.ival = value->idata;
  } else if(value->type == PIR_REAL) {
    result->type = uExpressionValue::rval;
    result->value.rval = value->ddata;
  } else if(value->type == PIR_BOOL) {
    result->type = uExpressionValue::ival;
    result->value.ival = (value->idata != 0);
  } else if(value->type == PIR_STRING) {
    result->type = uExpressionValue::sval;
    result->value.sval = strdup(value->sdata.c_str());
  } else {
    abort();
    return -1;
  }
  if(set) {
    variables[m2->group(0)->group(0)->substring()] = value;
  }
  return 0;
}

extern "C" uExpressionValue Util_ExpressionParseEvaluate(const char *expr) {
  void *m2 = Util_ExpressionParse(expr);
  uExpressionValue uvalue;
  Util_ExpressionEvaluate(m2,&uvalue,0,0);
  Util_ExpressionFree(m2);
  return uvalue;
}

extern "C" int cctk_PirahaParser(const char *buffer,unsigned long buffersize,int (*set_function)(const char *, const char *, int)) {
    std::string active;
    smart_ptr<Matcher> m2 = new Matcher(par_file_grammar,"file",buffer,buffersize);
    //std::clock_t st = std::clock();
    bool b = m2->matches();
    //std::clock_t en = std::clock();
    //std::cout << "PARSE TIME = " << ((en-st)/CLOCKS_PER_SEC) << std::endl;
    if(b) {
        int line = -1;
        for(int i=0;i<m2->groupCount();i++) {
            smart_ptr<Group> gr = m2->group(i);
            if(gr->group(0)->getPatternName() == "active") {
                smart_ptr<Value> smv = meval(gr->group(1),0);
                std::string val = smv->copy();
                active += val;
                active += ' ';
                line = gr->line();
            }
        }
        set_function("ActiveThorns",active.c_str(),line);
        for(int i=0;i<m2->groupCount();i++) {
            smart_ptr<Group> gr = m2->group(i);
            if(gr->getPatternName() == "set") {
                smart_ptr<Group> par = gr->group("par");
                if(par.valid()) {
                    // add value->quot->inquot
                    std::string key;
                    std::string thorn = par->group("name",0)->substring();
                    std::string name  = par->group("name",1)->substring();
                    key += thorn;
                    key += "::";
                    key += name;
                    const cParamData *data = CCTK_ParameterData(name.c_str(),thorn.c_str());

                    smart_ptr<Group> index = par->group("parindex");
                    if(index.valid()) {
                        key += '[';
                        smart_ptr<Value> vv = meval(index,0);
                        if(vv->type != PIR_INT) {
                            std::ostringstream msg;
                            std::string par = get_parfile();
                            msg << "bad index " << vv << std::endl;
                            CCTK_Error(index->line(),par.c_str(),thorn.c_str(),msg.str().c_str());
                        }
                        std::string vvstr = vv->copy();
                        key += vvstr;
                        key += ']';
                    }

                    std::string val;
                    smart_ptr<Group> aexpr = gr->group("expr");
                    if(aexpr.valid()) {
                        current_thorn = thorn;
                        smart_ptr<Value> smv = meval(aexpr,0);
                        val = smv->copy();
                        assert(smv.valid());
                        smv->integerize();
                        if(data != NULL) {
                            if(data->type == PARAMETER_REAL)
                                smv->integerize();
                            if(data->type == PARAMETER_BOOLEAN)
                                smv->booleanize(gr);
                            check_types(thorn.c_str(),aexpr->line(),smv,data->type);
                        }
                        read_write status = read_write_status(key);
                        if(status == read) {
                          std::ostringstream msg;
                          msg << "Write after read: " << key << std::endl;
                          std::string par = get_parfile();
                          CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
                        } 
                        set_function(
                                key.c_str(),
                                val.c_str(),
                                gr->group(0)->line());
                    } else {
                        smart_ptr<Group> arr = gr->group("array");
                        for(int i=0;i<arr->groupCount();i++) {
                            aexpr = arr->group(i);
                            std::ostringstream keyi;
                            keyi << key << '[' << i << ']';
                            smart_ptr<Value> smv = meval(aexpr,0);
                            val = smv->copy();
                            if(data != NULL) {
                                if(data->type == PARAMETER_REAL)
                                    smv->integerize();
                                if(data->type == PARAMETER_BOOLEAN)
                                    smv->booleanize(gr);
                                check_types(thorn.c_str(),aexpr->line(),smv,data->type);
                            }
                            std::string keyi_str = keyi.str();
                            read_write status = read_write_status(keyi_str);
                            if(status == read) {
                              std::ostringstream msg;
                              msg << "Write after read: " << keyi.str() << std::endl;
                              std::string par = get_parfile();
                              CCTK_Error(gr->line(),par.c_str(),current_thorn.c_str(),msg.str().c_str());
                            } 
                            set_function(
                                    keyi.str().c_str(),
                                    val.c_str(),
                                    aexpr->line());
                        }
                    }
                }
            } else if(gr->getPatternName() == "set_var") {
                variables[gr->group(0)->substring()] = meval(gr->group(1),0);
            }
        }
    } else {
        std::ostringstream msg;
        msg << "ERROR IN PARAMETER FILE:";
        if(m2->inrule_max == "file::set::par" && m2->foundChar() =='=') {
            msg << std::endl;
            msg << "Invalid assignment." << std::endl;
            msg << "Valid assignments are: " << std::endl;
            msg << "ActiveThorns = \"...\"" << std::endl;
            msg << "Arrangement::Thorn = ..." << std::endl;
            // Construct the parse tree for the partial match
            if(m2->matchesTo(m2->max_pos)) {
                // find the last element in the parse tree
                smart_ptr<Group> gr = m2->group(m2->groupCount()-1);
                while(true) {
                    int n = gr->groupCount();
                    if(n > 0)
                        gr = gr->group(n-1);
                    else
                        break;
                }
                msg << "You wrote: " << gr->substring() << " = ..." << std::endl;
            }
        } else {
            m2->showError(msg);
        }
        std::string par = get_parfile();
        CCTK_Warn(0,m2->line(),par.c_str(),"cactus",msg.str().c_str());
        return 1;
    }
    return 0;
}

}
