OUTPUT_NAME=thesis
INPUT_NAME=thesis
BUILD_FOLDER=build
BIB_FILE=bibliography.bib

INPUT=${INPUT_NAME}.tex
OUTPUT=$(BUILD_FOLDER)/${OUTPUT_NAME}.pdf
DEPENDS=$(wildcard *.tex) ${INPUT}

all: ${OUTPUT}

.PHONY: everything
everything: images ${OUTPUT}

.PHONY: images
images:
	$(MAKE) -C $@ all

${OUTPUT}: $(BUILD_FOLDER) $(INPUT) $(DEPENDS)
	latexmk -pdf -jobname=$(BUILD_FOLDER)/$(OUTPUT_NAME) ${INPUT_NAME}

$(BUILD_FOLDER):
	mkdir -p $(BUILD_FOLDER)

.PHONY: clean
clean:
	rm -f *.aux
	rm -rf $(BUILD_FOLDER)
	#rm -f $(OUTPUT_NAME).dvi $(OUTPUT_NAME).lot $(OUTPUT_NAME).log $(OUTPUT_NAME).toc $(OUTPUT_NAME).fls $(OUTPUT_NAME).fdb_latexmk $(OUTPUT_NAME).out $(OUTPUT_NAME).pdf $(OUTPUT_NAME).lof $(OUTPUT_NAME).bcf $(OUTPUT_NAME).run.xml $(OUTPUT_NAME).blg $(OUTPUT_NAME).bbl
