# BaCon: Bayesian Correlation Analysis of Networks

BaCon is a tool for performing Bayesian Correlation Analysis on network data. This script allows users to analyze correlation matrices with optional input matrices using a specified number of CPUs for parallel processing.

## Requirements

- numba
- pandas
- numpy
- argparse 


## Usage

```bash
python bacon.py -c <correlation_matrix> [-i <input_matrix1>] [-i2 <input_matrix2>] [-ncpu <number_of_cpus>] [-o <output_file>]
```

## Arguments

| command  | command (long)  | description  |
| ------------ | ------------ | ------------ |
|-c   | --corr_matrix (required)  |Path to the correlation matrix file.     |
| -i  |  --input1 (optional) | Path to the first input matrix file.   |
| -i2 |--input2 (optional)   |   Path to the second input matrix file.   |
| -ncpu  | --n_cpu (optional)  | Number of CPUs to use for parallel processing. Default is 6.    |
| -o  |  --output (optional) |  Path to the output file. Default is 'bacon.csv'.    |


## Example
To run the script with a correlation matrix, one input matrix, and default settings for other parameters, use:

```python
python bacon.py -c path/to/correlation_matrix.csv 
```
```python
python bacon.py -i path/to/input_matrix1.csv -i2 path/to/input_matrix2.csv
```

To specify a different number of CPUs and an output file:

```python
python bacon.py -c path/to/correlation_matrix.csv -i path/to/input_matrix1.csv -ncpu 4 -o result.csv
```

## Description
BaCon performs Bayesian Correlation Analysis using the provided correlation matrix and optional input matrices. The results are written to the specified output file, or 'bacon.csv' by default.

## Author
Yasir Demirtaş, Thomas Rohde, Angela Shaw, Maximilian Billmann

## License
This project is licensed under the MIT License - see the LICENSE.md file for details.
