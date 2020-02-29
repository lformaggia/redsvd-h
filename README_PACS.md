#Version of RedSVD-h adapted for the course#

This is a software for SVD decomposition (useful for principal
component analysis) To install it and have it available in the pacs
directory hierarchy just type make

It's a template-only library so you just install the header file

**IMPORTANT NOTE** The file in `include/RedSVD/RedSVD-h` is an the
original version, while `include/RedSVD/RedSVD.hpp` is a new version
that uses C++ random number generator classes to create the random
matrix needed by the algorithm.
