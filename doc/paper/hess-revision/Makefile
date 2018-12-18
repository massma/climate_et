DISPOSABLES = review-1.aux  review-1.bcf  review-1.log review-1.bbl  review-1.blg  review-1.pdf  review-1.run.xml

review-1.pdf : references.bib review-1.tex
	pdflatex -shell-escape review-1.tex && bibtex review-1 && pdflatex -shell-escape review-1.tex

clean :
	rm -rf $(DISPOSABLES)
