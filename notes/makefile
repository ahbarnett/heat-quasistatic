# GNU makefile for latex notes. Thanks to the late D J C MacKay.
# Barnett 4/28/21

.PHONY: clean

.SUFFIXES: .tex .dvi .pdf

.tex.dvi: localrefs.bib
	latex $<
	if grep Citation $*.log > /dev/null; then bibtex $*; latex $*; fi
	if grep Rerun $*.log > /dev/null; then latex $*; fi

.dvi.pdf:
	dvipdf $<

SRC = $(wildcard *.tex)

default: all

all: $(SRC:%.tex=%.pdf)

clean:
	rm -f $(SRC:%.tex=%.pdf)
	rm -f $(SRC:%.tex=%.aux)
	rm -f $(SRC:%.tex=%.blg)
	rm -f $(SRC:%.tex=%.bbl)
	rm -f $(SRC:%.tex=%.out)
	rm -f $(SRC:%.tex=%.log)
	rm -f $(SRC:%.tex=%.dvi)
