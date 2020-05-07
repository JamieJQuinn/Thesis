OUTPUT_NAME=thesis
INPUT_NAME=thesis
BUILD_FOLDER=build
BIB_FILE=bibliography.bib

INPUT=${INPUT_NAME}.tex
OUTPUT=$(BUILD_FOLDER)/${OUTPUT_NAME}.pdf
DEPENDS=$(wildcard *.tex) ${INPUT}

all: ${OUTPUT}
	cp $(OUTPUT) .

.PHONY: everything
everything: images ${OUTPUT}

.PHONY: images
images:
	$(MAKE) -C $@ all

${OUTPUT}: $(BUILD_FOLDER) $(INPUT) $(DEPENDS)
	latexmk -pdf -outdir=$(BUILD_FOLDER) ${INPUT_NAME}

$(BUILD_FOLDER):
	mkdir -p $(BUILD_FOLDER)

.PHONY: clean
clean:
	rm -rf $(BUILD_FOLDER)
	rm -f $(OUTPUT_NAME).pdf

fetch-images:
	rsync -avz maths-hop:~/thesis/images .
