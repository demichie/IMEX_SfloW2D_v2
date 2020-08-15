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

# run TEST_BOX
cd TEST_BOX
echo ""
./cleanFolder.sh
../../bin/SW_VAR_DENS_MODEL
FILE=TEST_BOX_400_0001.p_2d
if test -f "$FILE"; then
    ./cleanFolder.sh
    echo "TEST_BOX succesfully run"
else
    ./cleanFolder.sh
    echo "Problems with TEST_BOX"
    exit 1
fi
cd ..

# run TEST_BW
cd TEST_BW
echo ""
./cleanFolder.sh
../../bin/SW_VAR_DENS_MODEL
FILE=TEST_BW_400_0001.p_2d
if test -f "$FILE"; then
    ./cleanFolder.sh
    echo "TEST_BW succesfully run"
else
    ./cleanFolder.sh
    echo "Problems with TEST_BW"
    exit 1
fi
cd ..


# run TEST_CF
cd TEST_CF
echo ""
./cleanFolder.sh
cp ../../EXAMPLES/EXAMPLE_CF/topo.zip .
unzip topo.zip
../../bin/SW_VAR_DENS_MODEL
FILE=TEST_CF_0001.p_2d
if test -f "$FILE"; then
    ./cleanFolder.sh
    echo "TEST_CF succesfully run"
else
    ./cleanFolder.sh
    echo "Problems with TEST_CF"
    exit 1
fi
cd ..

# run TEST_TAAL
cd TEST_TAAL
echo ""
./cleanFolder.sh
ln -s ../../EXAMPLES/EXAMPLE_TAAL/topography_dem.asc .
../../bin/SW_VAR_DENS_MODEL
FILE=TEST_TAAL_0001.p_2d
if test -f "$FILE"; then
    ./cleanFolder.sh
    echo "TEST_TAAL succesfully run"
else
    ./cleanFolder.sh
    echo "Problems with TEST_TAAL"
    exit 1
fi
cd ..


exit 0




