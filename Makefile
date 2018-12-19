FIGURES = doc/paper/idealized_scale.pdf doc/paper/idealized_sign.pdf \
	doc/paper/joint_vpd_sigma.pdf doc/paper/map.pdf \
	doc/paper/swc_boxplot.pdf doc/paper/test_sign.pdf doc/paper/data_scatter.png \
doc/paper/concave.pdf

ARXIV_FILES = doc/paper/ms.tex doc/paper/abstract.tex \
	doc/paper/agufull08.bst doc/paper/def.tex doc/paper/body.tex \
	doc/paper/intro.tex doc/paper/conclusions.tex

PAPER_FILES = doc/paper/references.bib doc/paper/body.tex \
doc/paper/abstract.tex doc/paper/vpd_et_paper_hess.tex \
doc/paper/intro.tex doc/paper/conclusions.tex

TABLES = doc/paper/pft_params.tex doc/paper/vpd_crit.tex \
	doc/paper/stats.tex doc/paper/d2_solutions.tex doc/paper/flux_sites.tex

ALL_REQUIRE = dat/changjie/diagnosed_data.pkl

.PHONY: clean figure-all paper arxiv clearn-arxiv clean-bak clean-paper install

.EXPORT_ALL_VARIABLES:
# assumes that you run make from toplevel
PLOTS = $(CURDIR)/etc/plots
DATA = $(CURDIR)/dat


# installation, data and python dependencies
install : src/FLUXNET_citations dat/changjie/MAT_DATA
	pipenv install && mkdir -p doc/shared_figs

src/FLUXNET_citations : .gitmodules
	git submodule init && git submodule update

dat/changjie/MAT_DATA : ${DATA}/changjie/vpd_data.tar.gz
	cd ${DATA}/changjie && tar -xzvf vpd_data.tar.gz

${DATA}/changjie/vpd_data.tar.gz :
	cd ${DATA}/changjie && wget "http://www.columbia.edu/~akm2203/data/vpd_data.tar.gz"



# Full analysis (Changjie's data -> plot/table data)
dat/changjie/diagnosed_data.pkl : src/analysis.py src/codebase/*.py
	cd ./src && pipenv run python analysis.py


# Figures
figure-all :
	cd ./src/paper_figure_code && pipenv run python runall.py

doc/paper/map.pdf : doc/paper/map.pdf.bak
	cp doc/paper/map.pdf.bak doc/paper/map.pdf

doc/paper/idealized_scale.pdf : src/paper_figure_code/idealized_scale.py $(ALL_REQUIRE)
	cd src/paper_figure_code && pipenv run python idealized_scale.py

doc/paper/idealized_sign.pdf : src/paper_figure_code/idealized_sign.py $(ALL_REQUIRE)
	cd src/paper_figure_code && pipenv run python idealized_sign.py

doc/paper/joint_vpd_sigma.pdf : src/paper_figure_code/joint_vpd_sigma.py $(ALL_REQUIRE)
	cd src/paper_figure_code && pipenv run python joint_vpd_sigma.py

doc/paper/swc_boxplot.pdf : src/paper_figure_code/swc_boxplot.py $(ALL_REQUIRE)
	cd src/paper_figure_code && pipenv run python swc_boxplot.py

doc/paper/test_sign.pdf : src/paper_figure_code/test_sign.py $(ALL_REQUIRE)
	cd src/paper_figure_code && pipenv run python test_sign.py

doc/paper/concave.pdf : src/paper_figure_code/concave.py $(ALL_REQUIRE)
	cd src/paper_figure_code && pipenv run python concave.py

doc/paper/data_scatter.png : doc/paper/data_scatter.bak $(ALL_REQUIRE)
	cd doc/paper && cp data_scatter.bak data_scatter.png

doc/paper/data_scatter.bak : src/paper_figure_code/data_scatter.py $(ALL_REQUIRE)
	cd src/paper_figure_code && pipenv run python data_scatter.py
	cd doc/paper && mv data_scatter.png data_scatter.bak


# tables
$(TABLES) : src/paper_figure_code/tables.py
	cd src/paper_figure_code && pipenv run python tables.py
	cd doc/paper && \
	sed -i "s/{AU-Stp}/{AU-DaP}/" flux_sites.tex && \
	sed -i "s/{CA-SF2}/{CA-SF1}/" flux_sites.tex && \
	sed -i "s/{DE-Kli}/{DE-Gri}/" flux_sites.tex && \
	sed -i "s/{US-AR2}/{US-AR1}/" flux_sites.tex && \
	sed -i "s/{US-Ne3}/{US-Ne1}/" flux_sites.tex && \
	sed -i "s/{AU-Rig}/{AU-Gin}/" flux_sites.tex && \
	sed -i "s/{US-Ne2}/{US-Ne1}/" flux_sites.tex


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

paper : doc/paper/vpd_et_paper_hess.pdf

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
