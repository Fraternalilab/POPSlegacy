#! /bin/sh

echo "--------------------------------------------------------------"
echo " test2e                                                       "
echo "--------------------------------------------------------------"

valgrind --leak-check=full --show-leak-kinds=all ../src/pops --pdbml 1f3r.xml.gz --zipped --jsonOut || exit 1

