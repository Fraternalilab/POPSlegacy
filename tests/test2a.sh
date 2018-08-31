#! /bin/sh

echo "--------------------------------------------------------------"
echo " test2a                                                       "
echo "--------------------------------------------------------------"

valgrind --leak-check=full --show-leak-kinds=all ../src/pops --pdb 1f3r.pdb --compositionOut --typeOut --topologyOut --atomOut --residueOut --chainOut || exit 1

