FIGURES = doc/paper/fully_idealized.pdf doc/paper/fully_idealized_sign.pdf \
	doc/paper/fully_idealized_scaling.pdf doc/paper/concave.pdf

ARXIV_FILES = doc/paper/ms.tex \
	doc/paper/agufull08.bst doc/paper/def.tex \


PAPER_FILES = doc/paper/references.bib doc/paper/body.tex \
doc/paper/abstract.tex doc/paper/agu-james-submission.tex \
doc/paper/intro.tex doc/paper/conclusions.tex doc/paper/james-supplement.tex \
doc/paper/supp-figs.tex

TABLES = doc/paper/param_fixed.tex doc/paper/param_varying.tex

.PHONY: clean figure-all paper arxiv clearn-arxiv clean-bak clean-paper install

.EXPORT_ALL_VARIABLES:
# assumes that you run make from toplevel
PLOTS = $(CURDIR)/etc/plots
DATA = $(CURDIR)/dat

arxiv : doc/paper/arxiv-submission.tar


# installation, data and python dependencies
install : src/FLUXNET_citations dat/changjie/MAT_DATA/US-Ne3.mat
	mkdir -p doc/shared_figs doc/paper/anc

dat/changjie/MAT_DATA/US-Ne3.mat : ${DATA}/changjie/vpd_data.tar.gz
	cd ${DATA}/changjie && tar -xzvf vpd_data.tar.gz

${DATA}/changjie/vpd_data.tar.gz :
# cd ${DATA}/changjie && wget "http://www.columbia.edu/~akm2203/data/vpd_data.tar.gz"
# update May 9, 2022: I moved the data to Google drive, but it is lionmail-restricted.
# email Adam Massmann (akm2203@columbia.edu) if you would like the data.

# Full analysis (Changjie's data -> plot/table data)
dat/changjie/diagnosed_data.pkl : src/analysis.py
	cd ./src && python3 analysis.py

src/FLUXNET_citations : .gitmodules
	git submodule init && git submodule update


# Figures
$(FIGURES) : src/paper_figure_code/fully_idealized_figure.py dat/changjie/diagnosed_data.pkl
	cd src/paper_figure_code && python3 fully_idealized_figure.py
	cd src/paper_figure_code && python3 concave.py


# tables
$(TABLES) : src/paper_figure_code/fully_idealized_figure.py dat/changjie/diagnosed_data.pkl
	cd src/paper_figure_code && python3 fully_idealized_figure.py


# Manuscript targets

clean-arxiv :
	rm doc/paper/arxiv-submission.tar

doc/paper/anc/james-supplement.pdf : doc/paper/james-supplement.pdf
	cp $^ $@

doc/paper/arxiv-submission.tar : $(FIGURES) $(ARXIV_FILES) $(TABLES) doc/paper/ms.bbl doc/paper/anc/james-supplement.pdf
	tar -cvf $@ --transform 's/doc\/paper\///g' $^

doc/paper/ms.bbl : doc/paper/ms.pdf

doc/paper/ms.pdf : $(FIGURES) $(ARXIV_FILES) $(TABLES) doc/paper/references.bib
	cd ./doc/paper && pdflatex ms && \
	bibtex ms && pdflatex ms && \
	bibtex ms && pdflatex ms

doc/paper/james-supplement.pdf : doc/paper/james-supplement.tex doc/paper/references.bib doc/paper/map.pdf  doc/paper/supp-figs.tex doc/paper/supp-figs/0joint_rh_es.pdf
	cd ./doc/paper && pdflatex james-supplement && \
	bibtex james-supplement && pdflatex james-supplement && \
	bibtex james-supplement && pdflatex james-supplement

doc/paper/map.pdf : src/paper_figure_code/map.py
	cd src/paper_figure_code && python3 map.py

doc/paper/supp-figs.tex : src/paper_figure_code/swc_boxplot.py
	cd src/paper_figure_code && python3 swc_boxplot.py

doc/paper/references.bib : doc/paper/paper_references.bib doc/paper/flux_sites.bib
	cat $^ > $@

doc/paper/supp-figs/0joint_rh_es.pdf : src/paper_figure_code/joint_rh_es.py src/codebase/plot_tools.py
	cd src/paper_figure_code && python3 joint_rh_es.py



# Clean targets
clean :
	rm -f $(FIGURES) doc/paper/arxiv-submission.tar && \
	cd doc/paper && rm -f ./*.aux ./*.log ./*.blg ./*.bbl ms.pdf \
	agu-james-submission.pdf james-supplement.pdf  supp-figs/*.pdf supp-figs.tex

# below is if you don't want to regenerate figs
clean-paper :
	cd doc/paper && rm -f ./*.aux ./*.log ./*.blg ./*.bbl \
	data_scatter.png ms.pdf agu-james-submission.pdf

# below is because scatter figure takes so long to make, we don't want
# to always clean it
clean-bak :
	rm doc/paper/data_scatter.bak
