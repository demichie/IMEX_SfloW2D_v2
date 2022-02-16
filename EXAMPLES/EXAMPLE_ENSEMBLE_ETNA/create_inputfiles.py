import pandas as pd
import numpy as np
import os
import shutil as st


def replace_strings(working_dir, df , header , i_row):

    main_dir = os.getcwd()
    os.chdir(working_dir)
    directory = os.getcwd()

    for fname in os.listdir(directory):

        if os.path.isfile(fname) and "template" in fname:
            # print(fname)
            # Full path
            f = open(fname, 'r')
            filedata = f.read()

            for name in header:

                searchstring = 'ENSEMBLE_' + name
                # print(searchstring)
                f.seek(0)

                if searchstring in filedata:

                    print('found string in file %s' % fname)
                    print(searchstring)
                    f.seek(0)

                    for line in f:
                        # replacing the string and write to output file
                        filedata = filedata.replace(searchstring,
                                                    str(df.at[i_row, name]))

            f.close()
            f_out = open(fname.replace('template', 'inp'), 'w')
            f_out.write(filedata)
            f_out.close()

    os.chdir(main_dir)

##########################################

def main():

    df = pd.read_csv('samples.csv')

    templatedir = 'templatedir'

    print(df)

    header = df.columns

    print(header)

    n_rows = len(df.index)
    for i_row in range(n_rows):

        working_dir = "ensemble." + "{:05d}".format(i_row)
        working_dir = os.path.join(os.getcwd(), working_dir)
    
        st.copytree('templatedir', working_dir, symlinks=True)
    
        replace_strings(working_dir, df , header , i_row)
        

if __name__ == '__main__':

    main()
