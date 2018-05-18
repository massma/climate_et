FIGURES = doc/paper/idealized_scale.pdf doc/paper/idealized_sign.pdf \
	doc/paper/joint_vpd_sigma.pdf doc/paper/map.pdf \
	doc/paper/swc_boxplot.pdf doc/paper/test_sign.pdf doc/paper/data_scatter.png \
doc/paper/concave.pdf

ARXIV_FILES = doc/paper/ms.tex doc/paper/ms.bbl \
	doc/paper/agufull08.bst doc/paper/def.tex 

TABLES = doc/paper/pft_params.tex doc/paper/vpd_crit.tex \
	doc/paper/stats.tex doc/paper/d2_solutions.tex

.PHONY: data clean figure-all paper arxiv clearn-arxiv clean-bak clean-paper 

# you'll need to rerun data if you make any changes to the analysis script
data :
	cd ./src && python analysis.py

figure-all :
	cd ./src/paper_figure_code && python runall.py

doc/paper/idealized_scale.pdf :
	cd src/paper_figure_code && python idealized_scale.py

doc/paper/idealized_sign.pdf :
	cd src/paper_figure_code && python idealized_sign.py

doc/paper/joint_vpd_sigma.pdf :
	cd src/paper_figure_code && python joint_vpd_sigma.py

doc/paper/map.pdf :
	cd src/paper_figure_code && python map.py

doc/paper/swc_boxplot.pdf :
	cd src/paper_figure_code && python swc_boxplot.py

doc/paper/test_sign.pdf :
	cd src/paper_figure_code && python test_sign.py

doc/paper/concave.pdf :
	cd src/paper_figure_code && python concave.py

doc/paper/data_scatter.png : doc/paper/data_scatter.bak
	cd doc/paper && cp data_scatter.bak data_scatter.png

doc/paper/data_scatter.bak :
	cd src/paper_figure_code && python data_scatter.py
	cd doc/paper && mv data_scatter.png data_scatter.bak

doc/paper/ms.bbl : doc/paper/vpd_et_paper.bbl
	cp $< $@

doc/paper/ms.tex : doc/paper/vpd_et_paper_arxiv.tex
	cp $< $@

arxiv : doc/paper/arxiv-submission.tar

clean-arxiv :
	rm doc/paper/arxiv-submission.tar

doc/paper/arxiv-submission.tar : $(FIGURES) $(ARXIV_FILES) $(TABLES)
	tar -cvf $@ --transform 's?.*/??g' $^

paper : doc/paper/vpd_et_paper.pdf

doc/paper/vpd_et_paper.pdf : $(FIGURES)
	cd ./doc/paper && pdflatex vpd_et_paper && \
	bibtex vpd_et_paper && pdflatex vpd_et_paper && \
	bibtex vpd_et_paper && pdflatex vpd_et_paper

clean :
	rm -f $(FIGURES) doc/paper/arxiv-submission.tar && \
	cd doc/paper && rm -f ./*.aux ./*.log ./*.blg ./*.bbl vpd_et_paper.pdf

# below is if you don't want to regenerate figs
clean-paper :
	cd doc/paper && rm -f ./*.aux ./*.log ./*.blg ./*.bbl \
	data_scatter.png vpd_et_paper.pdf

# below is because scatter figure takes so long to make, we don't want
# to always clean it
clean-bak :
	rm doc/paper/data_scatter.bak
