#! /bin/sh

echo "--------------------------------------------------------------"
echo " test2d                                                       "
echo "--------------------------------------------------------------"

valgrind --leak-check=full --show-leak-kinds=all ../src/pops --pdb 1f3r.pdb.gz --zipped --jsonOut || exit 1

