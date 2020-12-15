#####
# define version
#####
NAMD_VERSION = 2.12

#####
# specific names of programs used in make commands
#####
CCPP		= CC
CC		= cc
CD		= cd
COPY		= cp
MOVE		= mv -f
DELETE		= rm -f
ECHO		= echo
TAR		= tar
COMPRESS	= compress
AR		= ar rsc
RANLIB		= /bin/true
COCMD		= co
SHELL 		= /bin/sh


DOCDIR		= .
MAIN_UG_FILE	= $(DOCDIR)/ug

SRCFILES = ug.bib namd_begindoc.tex namd_copyright.tex namd_macros.tex \
namd_title.tex namd_version.tex ug.tex ug_accel.tex ug_avail.tex \
ug_dynamics.tex ug_forcefield.tex ug_performance.tex ug_userdef.tex \
ug_intro.tex ug_io.tex ug_macros.tex ug_runit.tex ug_sample.tex ug_start.tex \
ug_xplor.tex ug_fenergy.tex ug_psfgen.tex psfgen_macros.tex ug_alchemy.tex \
ug_colvars-cv.tex colvars-cv.tex ug_colvars.bib ug_colvars_macros.tex \
ug_colvars.tex ug_analysis.tex ug_gbis.tex

PDFFIGS = figures/fmaOn.pdf figures/pairlistdist.pdf \
figures/shifting.pdf figures/switching.pdf figures/dual_top.pdf

### run pdflatex
pdf:	$(DOCDIR)/ug.pdf

### run latex2html
html:	$(DOCDIR)/ug.tar.gz

$(DOCDIR)/ug.tar.gz: $(DOCDIR)/ug.pdf
	latex2html -local_icons -split 5 -toc_depth 4 -toc_stars \
	-t "NAMD $(NAMD_VERSION) User's Guide" -address "http://www.ks.uiuc.edu/Research/namd/" \
	-up_url "../../" -up_title "Main NAMD Page" ug
	chmod -R a+rX ug
	tar czvf ug.tar.gz ug

$(DOCDIR)/ug.pdf: $(SRCFILES) $(PDFFIGS)
	$(ECHO) "Creating NAMD $(NAMD_VERSION) User's Guide ..." ; \
	$(CD) $(DOCDIR); \
	pdflatex $(MAIN_UG_FILE); \
	bibtex $(MAIN_UG_FILE); \
	bibtex ug_alchemy; \
	makeindex $(MAIN_UG_FILE).idx; \
	pdflatex $(MAIN_UG_FILE); \
	makeindex $(MAIN_UG_FILE).idx; \
	pdflatex $(MAIN_UG_FILE); \

DISTILL = distill -compatlevel 2.1

# need to specify page size because distill can't figure it out (yet?)

#figures/fmaOn.pdf: figures/fmaOn.eps
#	$(DISTILL) -pagesize 301 174 pts < $^ > $@

#figures/pairlistdist.pdf: figures/pairlistdist.eps
#	$(DISTILL) -pagesize 240 240 pts < $^ > $@

#figures/shifting.pdf: figures/shifting.eps
#	$(DISTILL) -pagesize 301 165 pts < $^ > $@

#figures/switching.pdf: figures/switching.eps
#	$(DISTILL) -pagesize 310 165 pts < $^ > $@

#figures/dual_top.pdf: figures/dual_top.eps
#	$(DISTILL) -pagesize 584 220 pts < $^ > $@

$(DOCDIR)/ug.dvi: $(SRCFILES)
	$(ECHO) "Creating NAMD $(NAMD_VERSION) User's Guide ..." ; \
	$(CD) $(DOCDIR); \
	latex $(MAIN_UG_FILE); \
	bibtex $(MAIN_UG_FILE); \
	makeindex $(MAIN_UG_FILE).idx; \
	latex $(MAIN_UG_FILE); \
	makeindex $(MAIN_UG_FILE).idx; \
	latex $(MAIN_UG_FILE); \
	# dvips -t letter -o ug.ps ug.dvi

namd_version.tex: Makefile
	$(DELETE) namd_version.tex ; \
	env echo '\newcommand{\NAMDVER}{$(NAMD_VERSION)}' > namd_version.tex

psfgen_macros.tex: ../psfgen/doc/ug_macros.tex
	$(DELETE) psfgen_macros.tex ; \
	grep PSFGEN ../psfgen/doc/ug_macros.tex > psfgen_macros.tex

ug_psfgen.tex: ../psfgen/doc/ug_psfgen.tex
	$(DELETE) ug_psfgen.tex ; \
	sed -e '20,$$ s/\\subsection{/\\subsubsection{/' \
	    -e '20,$$ s/\\section{/\\subsection{/' \
	    -e '20,$$ s/\\subsection\*{/\\subsubsection\*{/' \
	../psfgen/doc/ug_psfgen.tex > ug_psfgen.tex

clean:
	$(DELETE) namd_version.tex ; \
	$(DELETE) psfgen_macros.tex ; \
	$(DELETE) ug_psfgen.tex ; \
	$(DELETE) $(DOCDIR)/ug.ps ; \
	$(DELETE) $(DOCDIR)/ug.aux ; \
	$(DELETE) $(DOCDIR)/ug.bbl ; \
	$(DELETE) $(DOCDIR)/ug.blg ; \
	$(DELETE) $(DOCDIR)/ug.dvi ; \
	$(DELETE) $(DOCDIR)/ug.lof ; \
	$(DELETE) $(DOCDIR)/ug.log ; \
	$(DELETE) $(DOCDIR)/ug.idx ; \
	$(DELETE) $(DOCDIR)/ug.ilg ; \
	$(DELETE) $(DOCDIR)/ug.ind ; \
	$(DELETE) $(DOCDIR)/ug.toc ; \
	$(DELETE) $(DOCDIR)/ug.pdf ; \
	$(DELETE) $(DOCDIR)/ug.tar.gz; \
#	$(DELETE) $(PDFFIGS)
