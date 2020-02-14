NOTEBOOKS=$(wildcard *.Rmd)
SCRIPTS=$(NOTEBOOKS:.Rmd=.py)
RUN_COMMANDS=$(addprefix run_,$(SCRIPTS))
PLOTS=$(wildcard *.pdf)

#%.py: %.ipynb
	#jupyter nbconvert --to script $<

%.py: %.Rmd
	jupytext --to script $<

run_%.py: %.py
	ipython $<

.PHONY: all
all: $(RUN_COMMANDS)

.PHONY: clean
clean:
	rm -f $(SCRIPTS)
	rm -f $(PLOTS)
