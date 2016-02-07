#! /bin/sh

echo "--------------------------------------------------------------"
echo "POPS running on single test structure"
echo "--------------------------------------------------------------"
../src/pops --pdb 1f3r.pdb --compositionOut --typeOut --topologyOut --atomOut --residueOut --chainOut || exit 1

