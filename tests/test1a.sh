#! /bin/sh

valgrind --leak-check=full --show-leak-kinds=all ../src/pops --pdb 1f3r.pdb --coarse --compositionOut --typeOut --topologyOut --atomOut --residueOut --chainOut || exit 1

