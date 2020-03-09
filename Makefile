# This makefile just installs the headers in the main include directory of
# the examples for the PACS course
# For the examples look at the src/examples and the src/pacs_examples directory
#
# Installing to the main include directory is not necessary to test the examples
#
include ../../../Makefile.inc
all: install Description.pdf test_redSVD
install:
	install -d -v $(PACS_INC_DIR)/RedSVD
	install  -p -v -t $(PACS_INC_DIR)/RedSVD/ ./include/RedSVD/*.hpp
Description.pdf:Description.tex
	pdflatex Description.tex
test_redSVD:test_redSVD.cpp 
clean:
	$(RM) *.log *.aux *.dvi *~
distclean: clean
	$(RM) Description.pdf
