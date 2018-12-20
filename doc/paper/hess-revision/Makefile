DISPOSABLES = review-1.aux  review-1.bcf  review-1.log review-1.bbl  review-1.blg  review-1.pdf  review-1.run.xml review-2.aux  review-2.bcf  review-2.log review-2.bbl  review-2.blg  review-2.pdf  review-2.run.xml

review-2.pdf : references.bib review-2.tex
	pdflatex review-2.tex && bibtex review-2 && pdflatex review-2.tex && bibtex review-2 && pdflatex review-2.tex

review-1.pdf : references.bib review-1.tex
	pdflatex review-1.tex && bibtex review-1 && pdflatex review-1.tex && bibtex review-1 && pdflatex review-1.tex


clean :
	rm -rf $(DISPOSABLES)
