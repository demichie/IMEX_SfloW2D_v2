#!/bin/bash
# run TEST_2D
cd TEST_2D
echo ""
./cleanFolder.sh
../../bin/SW_VAR_DENS_MODEL
FILE=TEST_2D_0001.p_2d
if test -f "$FILE"; then
    ./cleanFolder.sh
    echo "TEST_2D succesfully run"
else
    ./cleanFolder.sh
    echo "Problems with TEST_2D"
    exit 1
fi
cd ..


echo ""
echo "All tests succesfully passed"
echo ""
exit 0
