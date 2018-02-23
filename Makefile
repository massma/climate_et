FIGURES = doc/paper/idealized_scale.pdf doc/paper/idealized_sign.pdf \
	doc/paper/joint_vpd_sigma.pdf doc/paper/map.pdf \
	doc/paper/swc_boxplot.pdf doc/paper/test_sign.pdf doc/paper/data_scatter.png
.PHONY: data clean figure-all paper clean-bak clean-paper

# you'll need to rerun data if you make any changes to the analysis script
data :
	cd ./src && python3 analysis.py

figure-all :
	cd ./src/paper_figure_code && python3 runall.py

doc/paper/idealized_scale.pdf :
	cd src/paper_figure_code && python3 idealized_scale.py

doc/paper/idealized_sign.pdf :
	cd src/paper_figure_code && python3 idealized_sign.py

doc/paper/joint_vpd_sigma.pdf :
	cd src/paper_figure_code && python3 joint_vpd_sigma.py

doc/paper/map.pdf :
	cd src/paper_figure_code && python3 map.py

doc/paper/swc_boxplot.pdf :
	cd src/paper_figure_code && python3 swc_boxplot.py

doc/paper/test_sign.pdf :
	cd src/paper_figure_code && python3 test_sign.py

doc/paper/data_scatter.png : doc/paper/data_scatter.bak
	cd doc/paper && cp data_scatter.bak data_scatter.png

doc/paper/data_scatter.bak :
	cd src/paper_figure_code && python3 data_scatter.py
	cd doc/paper && mv data_scatter.png data_scatter.bak

paper : $(FIGURES) 
	cd ./doc/paper && pdflatex vpd_et_paper && \
	bibtex vpd_et_paper && pdflatex vpd_et_paper && \
	bibtex vpd_et_paper && pdflatex vpd_et_paper

clean :
	rm $(FIGURES) && \
	cd doc/paper && rm ./*.aux ./*.log ./*.blg ./*.bbl vpd_et_paper.pdf

# below is if you don't want to regenerate figs
clean-paper :
	cd doc/paper && rm ./*.aux ./*.log ./*.blg ./*.bbl vpd_et_paper.pdf

# below is because scatter figure takes so long to make, we don't want
# to always clean it
clean-bak :
	rm doc/paper/data_scatter.bak
