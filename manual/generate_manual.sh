#/bin/bash

pandoc -s MANUAL.md INSTALLATION.md MAIN_INPUT.md LATTICE.md XGENESIS.md -V 'mainfont:SansSerif.ttf' -V geometry:margin=1in -V colorlinks=true -V linkcolor=blue -V urlcolor=red -V toccolor=gray -f markdown -t latex -s -o Manual.pdf