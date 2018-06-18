#! /bin/sh

echo "--------------------------------------------------------------"
echo "POPS running on single test structure"
echo "--------------------------------------------------------------"
valgrind --leak-check=full --show-leak-kinds=all ../src/pops --pdb 1f3r.pdb --jsonOut || exit 1

