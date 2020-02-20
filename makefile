OUTPUT_NAME=main
INPUT_NAME=main
BUILD_FOLDER=build
BIB_FILE=bibliography.bib

INPUT=${INPUT_NAME}.tex
OUTPUT=${OUTPUT_NAME}.pdf
DEPENDS=$(wildcard *.tex) ${INPUT}

all: ${OUTPUT}

.PHONY: everything
everything: images ${OUTPUT}

.PHONY: images
images:
	$(MAKE) -C $@ all

${OUTPUT}: $(INPUT) $(DEPENDS)
	latexmk -pdf -jobname=$(OUTPUT_NAME) ${INPUT_NAME}

.PHONY: clean
clean:
	rm -f *.aux
	rm -f $(OUTPUT_NAME).dvi $(OUTPUT_NAME).lot $(OUTPUT_NAME).log $(OUTPUT_NAME).toc $(OUTPUT_NAME).fls $(OUTPUT_NAME).fdb_latexmk $(OUTPUT_NAME).out $(OUTPUT_NAME).pdf $(OUTPUT_NAME).lof $(OUTPUT_NAME).bcf $(OUTPUT_NAME).run.xml $(OUTPUT_NAME).blg $(OUTPUT_NAME).bbl
