#%%
module_names = ['numba', 'pandas', 'numpy']
import subprocess
import sys
reqs = subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'])
installed_packages = [r.decode().split('==')[0] for r in reqs.split()]

if not all(elem in installed_packages  for elem in module_names):
    print('Installing required packages')
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'numba', 'pandas', 'numpy'])
else:
    print('All required packages are installed')

import pandas as pd
import numpy as np
from numba import prange, njit, set_num_threads
from numba import vectorize, float64
from pathlib import Path
import argparse

#%%

parser = argparse.ArgumentParser(description='BaCon: Bayesian Correlation Analysis of Networks.')
parser.add_argument('-c','--corr_matrix', help='Correlation Matrix', required=True)
parser.add_argument('-i','--input1', help='Input Matrix 1', required=False)
parser.add_argument('-i2','--input2', help='Input Matrix 2', required=False)
parser.add_argument('-ncpu','--n_cpu', help='Number of CPU', required=False, default=6, type=int)
parser.add_argument('-o','--output', help='Output file', required=False, default='bacon.csv')
args = vars(parser.parse_args())


#%%
if __name__ == "__main__":


    print("""

░▒▓███████▓▒░ ░▒▓██████▓▒░ ░▒▓██████▓▒░ ░▒▓██████▓▒░░▒▓███████▓▒░  
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
░▒▓███████▓▒░░▒▓████████▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
░▒▓███████▓▒░░▒▓█▓▒░░▒▓█▓▒░░▒▓██████▓▒░ ░▒▓██████▓▒░░▒▓█▓▒░░▒▓█▓▒░ 


    """)

    #%%
    # count extreme values
    #####################################################################################
    @njit(parallel=True)
    def count_extreme(mat):
        cols, rows = mat.shape
        z = np.zeros((cols, rows), dtype=np.float64)
        val = 0.0
        for i in prange(cols):
            for j in range(rows):
                cell_old_value = mat[i, j]

                if mat[i, j] < 0:
                    mat[i, j] += 0.05
                    val = -(mat[i, :] < mat[i, j]).sum()
                    z[i, j] =  val
                else:
                    mat[i, j] -= 0.05
                    val = (mat[i, :] > mat[i, j]).sum()
                    z[i, j] =  val

                mat[i, j] = cell_old_value
        return z


    # normalize
    #####################################################################################
    @vectorize([float64(float64)])
    def f(x):
        return -(1+x) if x < 0 else 1-x

    #%%
    # normalize
    #####################################################################################
    @vectorize([float64(float64)])
    def f2(x):
        return x+2 if x>0 else x-2




    # cross correlation
    #####################################################################################

    def cross_correlation(A,B):

        # Get number of rows in either A or B
        N = B.shape[0]

        # Store columnw-wise in A and B, as they would be used at few places
        sA = A.sum(0)
        sB = B.sum(0)

        # Basically there are four parts in the formula. We would compute them one-by-one
        p1 = N*np.dot(B.T,A)
        p2 = sA.values*sB.values[:,None]
        p3 = N*((B**2).sum(0)) - (sB**2)
        p4 = N*((A**2).sum(0)) - (sA**2)

        # Finally compute Pearson Correlation Coefficient as 2D array 
        pcorr = ((p1 - p2)/np.sqrt(p4.values*p3.values[:,None]))

        corr = pd.DataFrame(pcorr)
        corr.columns = A.columns
        corr.index = B.columns

        return corr


    # Start count extreme
    #####################################################################################
    def start_count_extreme(corr, n_core, output_filename=False):
        mat = corr.to_numpy().copy()
        print(f'*** Number of CPU: {n_core}')
        set_num_threads(n_core)
        print(f'*** Counting extreme values')
        matT = np.ascontiguousarray(mat.T)
        print(f'*** Transposing matrix')
        cols  = count_extreme(matT) 
        rows  = count_extreme(mat) 
        print(f'*** Normalizing matrix')
        all = rows + cols.T
        all = f2(all)
        norm = all / (mat.shape[0] + mat.shape[1])
        normed = f(norm)
        print(f'*** Normalized matrix shape: {normed.shape}')
        print(f'*** Creating output dataframe')
        df = pd.DataFrame(normed, index=corr.index, columns=corr.columns)
        if output_filename:
            print(f'*** Saving file to: {output_filename}')
            df.to_csv(output_filename)
        return  df




    #%%

    def empirical_pval(stat_null, stat):
        """Calculate empirical p-value.

        Calculate empirical p-value based on the observed (surrogate maps)
        and expected (reference map) correlation scores.

        Parameters
        ----------
        stat_null: numpy.array
            A vector or matrix (nxm) of simulated or data-resampled null test
            statistics, where m is the number of test statistics and n is the
            number of null values in each distribution.
            The ith column thus corresponds to the null distribution of the ith
            test statistic in 'stat'.
        stat: numpy.array (m,)
            A vector of calculated test statistics.

        Returns
        -------
        p-values : numpy.array
            Calculated empirical pvalues.
        """

        # #%%
        # """Test empirical_pval."""
        # stat_real = np.array([4, 5, 1, 8, 1, 19, 12, 2, 3])
        # stat_null = np.arange(90).reshape(10, 9)
        # pval = h.empirical_pval(stat_null, stat_real)
        # assert pval.shape[0] == 9
        # # %%

        n_null_vals, n_tests = stat_null.shape
        assert n_tests == len(
            stat
        ), "Number of test statistics in 'stat_null' an 'stat' are mismatched"

        check = np.sum(np.abs(stat_null) >= np.abs(stat), axis=0)
        pvalues = (check + 1) / (n_null_vals + 1)
        return pvalues




    def quick_sort(df):
        sorted_indices = np.argsort(df.score.values)[::-1][:len(df.score)]
        #df.iloc[df.score.argsort()[::-1][:len(df.score)]
        return df.iloc[sorted_indices].reset_index(drop=True)



    # check if corr_matrix is a file and exists
    if w['corr_matrix'] != True:
        print(f'*** Correlation matrix: {args["corr_matrix"]}')

        if Path(args['corr_matrix']).is_file():
            corr = pd.read_csv(args['corr_matrix'])
        else:
            print('Correlation matrix file does not exist')
            exit()
    else:
    # check if input1 and input2 are files and exists
        if Path(args['input1']).is_file() and Path(args['input2']).is_file():
            print(f'*** Input matrix 1: {args["input1"]}')
            print(f'*** Input matrix 2: {args["input2"]}')
            f1 = pd.read_csv(args['input1'], index_col=0)
            f2 = pd.read_csv(args['input2'], index_col=0)
            print(f'*** Applying cross-correlation to {f1.shape} and {f2.shape}')
            corr = cross_correlation(f1,f2)
        else:
            print('Input files do not exist')
            exit()

    print(f'*** Correlation matrix shape: {corr.shape}')
    print(f'*** Starting BaCon analysis')
    bacon = start_count_extreme(corr, 6, output_filename=args['output'])




