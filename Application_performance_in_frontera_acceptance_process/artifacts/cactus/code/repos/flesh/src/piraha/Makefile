# Build the standalone tool
#
BIN ?= .
CXX ?= c++
CXXFLAGS ?= -g -O2 -I . -I ../include
SRCS = AutoGrammar.cc Bracket.cc Dot.cc End.cc \
	Generic.cc Grammar.cc Group.cc ILiteral.cc \
	Literal.cc Lookup.cc Matcher.cc Multi.cc \
	Or.cc ReParse.cc Seq.cc smart_ptr.cc Start.cc \
	Boundary.cc NegLookAhead.cc
OBJS = $(SRCS:%.cc=%.o)

libpiraha.a : $(OBJS)
	ar cr libpiraha.a $(OBJS)

$(BIN)/generic : libpiraha.a ../include/Piraha.h Piraha.hpp smart_ptr.hpp
	mkdir -p $(BIN)
	$(CXX) $(CXXFLAGS) -o $(BIN)/generic $(OBJS)

clean :
	rm -f $(BIN)/generic libpiraha.a $(OBJS)
