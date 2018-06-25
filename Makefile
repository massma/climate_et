FIGURES = doc/paper/idealized_scale.pdf doc/paper/idealized_sign.pdf \
	doc/paper/joint_vpd_sigma.pdf doc/paper/map.pdf \
	doc/paper/swc_boxplot.pdf doc/paper/test_sign.pdf doc/paper/data_scatter.png \
doc/paper/concave.pdf

ARXIV_FILES = doc/paper/ms.tex doc/paper/abstract.tex \
	doc/paper/agufull08.bst doc/paper/def.tex doc/paper/body.tex

PAPER_FILES = doc/paper/references.bib doc/paper/body.tex \
doc/paper/abstract.tex doc/paper/vpd_et_paper_afm.tex

TABLES = doc/paper/pft_params.tex doc/paper/vpd_crit.tex \
	doc/paper/stats.tex doc/paper/d2_solutions.tex doc/paper/flux_sites.tex

.PHONY: data clean figure-all paper arxiv clearn-arxiv clean-bak clean-paper

# you'll need to rerun data if you make any changes to the analysis script
data :
	cd ./src && python3 analysis.py

figure-all :
	cd ./src/paper_figure_code && python3 runall.py

doc/paper/idealized_scale.pdf : src/paper_figure_code/idealized_scale.py
	cd src/paper_figure_code && python3 idealized_scale.py

doc/paper/idealized_sign.pdf : src/paper_figure_code/idealized_sign.py
	cd src/paper_figure_code && python3 idealized_sign.py

doc/paper/joint_vpd_sigma.pdf : src/paper_figure_code/joint_vpd_sigma.py
	cd src/paper_figure_code && python3 joint_vpd_sigma.py

doc/paper/map.pdf : src/paper_figure_code/map.py
	cd src/paper_figure_code && python3 map.py

doc/paper/swc_boxplot.pdf : src/paper_figure_code/swc_boxplot.py
	cd src/paper_figure_code && python3 swc_boxplot.py

doc/paper/test_sign.pdf : src/paper_figure_code/test_sign.py
	cd src/paper_figure_code && python3 test_sign.py

doc/paper/concave.pdf : src/paper_figure_code/concave.py
	cd src/paper_figure_code && python3 concave.py

doc/paper/data_scatter.png : doc/paper/data_scatter.bak
	cd doc/paper && cp data_scatter.bak data_scatter.png

doc/paper/data_scatter.bak : src/paper_figure_code/data_scatter.py
	cd src/paper_figure_code && python3 data_scatter.py
	cd doc/paper && mv data_scatter.png data_scatter.bak

# note below sed's are to fix duplicate citations
doc/paper/flux_sites.tex : src/paper_figure_code/tables.py \
src/codebase/FLUXNET_citations/F15T1_LaTeX/fluxnet_pycite.py
	cd src/paper_figure_code && python3 tables.py
	cd doc/paper && \
	sed -i "s/{AU-Stp}/{AU-DaP}/" flux_sites.tex && \
	sed -i "s/{CA-SF2}/{CA-SF1}/" flux_sites.tex && \
	sed -i "s/{DE-Kli}/{DE-Gri}/" flux_sites.tex && \
	sed -i "s/{US-AR2}/{US-AR1}/" flux_sites.tex && \
	sed -i "s/{US-Ne3}/{US-Ne1}/" flux_sites.tex && \
	sed -i "s/{AU-Rig}/{AU-Gin}/" flux_sites.tex && \
	sed -i "s/{US-Ne2}/{US-Ne1}/" flux_sites.tex

arxiv : doc/paper/arxiv-submission.tar

clean-arxiv :
	rm doc/paper/arxiv-submission.tar

doc/paper/arxiv-submission.tar : $(FIGURES) $(ARXIV_FILES) $(TABLES) doc/paper/ms.bbl
	tar -cvf $@ --transform 's?.*/??g' $^

doc/paper/ms.bbl : doc/paper/ms.pdf

paper : doc/paper/ms.pdf

doc/paper/ms.pdf : $(FIGURES) $(ARXIV_FILES) $(TABLES) doc/paper/references.bib
	cd ./doc/paper && pdflatex ms && \
	bibtex ms && pdflatex ms && \
	bibtex ms && pdflatex ms

doc/paper/references.bib : doc/paper/paper_references.bib doc/paper/flux_sites.bib
	cat $^ > $@

doc/paper/vpd_et_paper_afm.pdf : $(FIGURES) $(PAPER_FILES) $(TABLES)
	cd ./doc/paper && pdflatex vpd_et_paper_afm && \
	bibtex vpd_et_paper_afm && pdflatex vpd_et_paper_afm && \
	bibtex vpd_et_paper_afm && pdflatex vpd_et_paper_afm

clean :
	rm -f $(FIGURES) doc/paper/arxiv-submission.tar && \
	cd doc/paper && rm -f ./*.aux ./*.log ./*.blg ./*.bbl ms.pdf \
	vpd_et_paper_afm

# below is if you don't want to regenerate figs
clean-paper :
	cd doc/paper && rm -f ./*.aux ./*.log ./*.blg ./*.bbl \
	data_scatter.png ms.pdf vpd_et_paper_afm.pdf

# below is because scatter figure takes so long to make, we don't want
# to always clean it
clean-bak :
	rm doc/paper/data_scatter.bak
