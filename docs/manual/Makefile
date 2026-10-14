# Built out of the source tree, into the project's builds/manual, like
# everything else: `make -C docs/manual`. OUT=<dir> to put it elsewhere.
SRC := $(abspath .)
OUT ?= $(abspath ../../../builds/manual)

all: $(OUT)/Manual.pdf

$(OUT)/Manual.pdf: Manual.tex intro.tex ch1.tex ch2.tex ap1.tex references.bib
	mkdir -p $(OUT)
	TEXINPUTS=$(SRC): BIBINPUTS=$(SRC): rubber -d --into $(OUT) $(SRC)/Manual.tex

clean:
	rm -rf $(OUT)

.PHONY: all clean
