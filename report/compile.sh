PDFLATEX=`which pdflatex`
BIBTEX=`which bibtex`

$PDFLATEX main.tex
$BIBTEX main
$PDFLATEX main.tex
rm *.aux *.log *.toc *.out *.bbl *.blg *.bcf *.xml *.snm

exit 0
