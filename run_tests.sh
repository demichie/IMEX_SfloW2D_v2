cd SW_VAR_DENS_MODEL-master/EXAMPLES

# TEST EXAMPLE_2D
cd EXAMPLE_2D
./cleanfolder.sh
python create_example.py 100 0.5 300 false
../../bin/SW_VAR_DENS_MODEL
FILE=example2D_0030.p_2d
if test -f "$FILE"; then
    ./cleanfolder.sh
    echo "EXAMPLE_2D succesfully run"
else
    echo "Problems with EXAMPLE_2D"
    exit 1
fi

# TEST EXAMPLE_BOX
cd EXAMPLE_BOX
./cleanfolder.sh
python create_example.py 400 0.5 300 false
../../bin/SW_VAR_DENS_MODEL
FILE=exampleBOX_400_0100.p_2d
if test -f "$FILE"; then
    ./cleanfolder.sh
    echo "EXAMPLE_BOX succesfully run"
else
    echo "Problems with EXAMPLE_BOX"
    exit 1
fi

exit 0


