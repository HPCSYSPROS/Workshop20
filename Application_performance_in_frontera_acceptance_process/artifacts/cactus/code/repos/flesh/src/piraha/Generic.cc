#include "Piraha.hpp"
#include <fstream>
#include <stdlib.h>

using namespace cctki_piraha;

void read_file(const char *file,std::string& buf) {
    std::ifstream in;
    in.open(file);
    while(true) {
        int c = in.get();
        if(c < 0)
            break;
        buf += (char)c;
    }
    in.close();
}

void usage() {
    std::cerr << "usage: generic [--perl|--python] grammar input" << std::endl;
    exit(2);
}

#define VAR(X) " " << #X << "=" << X

bool newEnd(std::string& in,const char *new_end,std::string& out) {
	out = "/dev/null";
	int n = in.size();
	int i;
	for(i=n-1;i>0;i--) {
		if(in[i] == '.') {
			break;
		}
		if(in[i] == '/') {
			std::cout << VAR(in) << VAR(i) << std::endl;
			break;
		}
	}
	if(in[i] != '.') {
		std::cout << VAR(in) << VAR(i) << std::endl;
		return false;
	}
	out.clear();
	for(int j=0;j<i;j++)
		out += in[j];
	out += new_end;
	return true;
}

int main(int argc,char **argv) {
	std::string grammarArg, inputArg;
	bool perlFlag = false;
    bool pythonFlag = false;
	bool oFlag = false;
    std::string outFile;
    int narg = 0;
    for(int n=1;n<argc;n++) {
    	std::string arg = argv[n];
    	if(arg == "--perl") {
    		perlFlag = true;
        } else if(arg == "--python") {
            pythonFlag =  true;
    	} else if(arg == "-o") {
    		outFile = argv[++n];
    		oFlag = true;
    	} else if(arg.size()>2 && arg[0]=='-' && arg[1]=='o') {
    		outFile = arg.substr(2,arg.size());
    		oFlag = true;
    	} else if(narg == 0) {
    		grammarArg = argv[n];
    		narg++;
    	} else if(narg == 1) {
    		inputArg = argv[n];
    		narg++;
    	} else {
    		usage();
    	}
    }
    if(!oFlag) {
    	if(perlFlag) {
    		newEnd(inputArg,".pm",outFile);
        } else if(pythonFlag) {
    		newEnd(inputArg,".py",outFile);
    	} else {
    		newEnd(inputArg,".pegout",outFile);
    	}
    }
    std::cout << "reading file: " << inputArg << std::endl;

    std::string grammar_file, input_file;
    read_file(grammarArg.c_str(),grammar_file);
    read_file(inputArg.c_str(),input_file);

    smart_ptr<Grammar> g = new Grammar();
    compileFile(g,grammar_file.c_str());
    smart_ptr<Matcher> mg =
        new Matcher (g,g->default_rule.c_str(),input_file.c_str());
    if(mg->matches()) {
        //std::vector<char> vec(4096);
        std::ofstream o;
        //o.rdbuf()->pubsetbuf(&vec.front(),vec.size());
        o.open(outFile.c_str());
        std::cout << "writing file: " << outFile << std::endl;
        smart_ptr<Group> src_file =
            new Group("annot:src_file",inputArg.c_str());
        mg->children->push_back(src_file);
    	if(perlFlag) {
    		mg->dumpPerl(o);
        } else if(pythonFlag) {
            mg->dumpPython(o);
    	} else {
    		mg->dump(o);
        }
        o.close();
    } else {
        mg->showError();
        return 1;
    }
    return 0;
}
