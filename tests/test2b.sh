#! /bin/sh

echo "--------------------------------------------------------------"
echo " test2b                                                       "
echo "--------------------------------------------------------------"

valgrind --leak-check=full --show-leak-kinds=all ../src/pops --pdb 1f3r.pdb --jsonOut || exit 1

