all: simple_drifter_dynamics

simple_drifter_dynamics: simple_drifter_dynamics.tex
	pdflatex -halt-on-error simple_drifter_dynamics.tex
	pdflatex -halt-on-error simple_drifter_dynamics.tex
	pdflatex -halt-on-error simple_drifter_dynamics.tex

.PHONY:
clean:
	rm -v *.{dvi,out,aux,toc,bbl,blg,fff,lof,log,lot,ttt,pdf}
