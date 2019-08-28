all: Manual.pdf

Manual.pdf: Manual.tex intro.tex ch1.tex ch2.tex ap1.tex references.bib
	rubber -d Manual

clean:
	rm Manual.pdf

.PHONY: clean
