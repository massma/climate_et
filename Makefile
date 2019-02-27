FIGURES = doc/paper/fully_idealized.pdf doc/paper/fully_idealized_sign.pdf \
	doc/paper/fully_idealized_scaling.pdf

ARXIV_FILES = doc/paper/ms.tex doc/paper/abstract.tex \
	doc/paper/agufull08.bst doc/paper/def.tex doc/paper/body.tex \
	doc/paper/intro.tex doc/paper/conclusions.tex

PAPER_FILES = doc/paper/references.bib doc/paper/body.tex \
doc/paper/abstract.tex doc/paper/vpd_et_paper_hess.tex \
doc/paper/intro.tex doc/paper/conclusions.tex

TABLES = doc/paper/param_fixed.tex doc/paper/param_varying.tex

ALL_REQUIRE = dat/changjie/diagnosed_data.pkl

.PHONY: clean figure-all paper arxiv clearn-arxiv clean-bak clean-paper install

.EXPORT_ALL_VARIABLES:
# assumes that you run make from toplevel
PLOTS = $(CURDIR)/etc/plots
DATA = $(CURDIR)/dat


paper : doc/paper/vpd_et_paper_hess.pdf


# installation, data and python dependencies
install : src/FLUXNET_citations dat/changjie/MAT_DATA
	pipenv install && mkdir -p doc/shared_figs

src/FLUXNET_citations : .gitmodules
	git submodule init && git submodule update


# Figures
$(FIGURES) : src/paper_figure_code/fully_idealized_figure.py
	cd src/paper_figure_code && pipenv run python fully_idealized_figure.py


# tables
$(TABLES) : src/paper_figure_code/fully_idealized_figure.py
	cd src/paper_figure_code && pipenv run python fully_idealized_figure.py


# Manuscript targets
arxiv : doc/paper/arxiv-submission.tar

clean-arxiv :
	rm doc/paper/arxiv-submission.tar

doc/paper/arxiv-submission.tar : $(FIGURES) $(ARXIV_FILES) $(TABLES) doc/paper/ms.bbl
	tar -cvf $@ --transform 's?.*/??g' $^

doc/paper/ms.bbl : doc/paper/ms.pdf

doc/paper/ms.pdf : $(FIGURES) $(ARXIV_FILES) $(TABLES) doc/paper/references.bib
	cd ./doc/paper && pdflatex ms && \
	bibtex ms && pdflatex ms && \
	bibtex ms && pdflatex ms

doc/paper/vpd_et_paper_hess.pdf : $(FIGURES) $(PAPER_FILES) $(TABLES)
	cd ./doc/paper && pdflatex vpd_et_paper_hess && \
	bibtex vpd_et_paper_hess && pdflatex vpd_et_paper_hess && \
	bibtex vpd_et_paper_hess && pdflatex vpd_et_paper_hess

doc/paper/references.bib : doc/paper/paper_references.bib doc/paper/flux_sites.bib
	cat $^ > $@


# Clean targets
clean :
	rm -f $(FIGURES) doc/paper/arxiv-submission.tar && \
	cd doc/paper && rm -f ./*.aux ./*.log ./*.blg ./*.bbl ms.pdf \
	vpd_et_paper_hess.pdf

# below is if you don't want to regenerate figs
clean-paper :
	cd doc/paper && rm -f ./*.aux ./*.log ./*.blg ./*.bbl \
	data_scatter.png ms.pdf vpd_et_paper_hess.pdf

# below is because scatter figure takes so long to make, we don't want
# to always clean it
clean-bak :
	rm doc/paper/data_scatter.bak
