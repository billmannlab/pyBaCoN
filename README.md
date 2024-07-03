# BaCoN

Buffering between genes is fundamental for robust cellular functions. While experimentally testing all possible gene pairs is infeasible, gene buffering can be predicted genome-wide under the assumption that a gene’s buffering capacity depends on its expression level and the absence of this buffering capacity primes a severe fitness phenotype of the buffered gene. We developed BaCoN (Balanced Correlation Network), a post-hoc unsupervised correction method that amplifies specific signals in expression-vs-fitness effect correlation-based networks. We quantified 147 million potential buffering relationships by associating CRISPR-Cas9-screening fitness effects with transcriptomic data across 1019 Cancer Dependency Map (DepMap) cell lines. BaCoN outperformed state-of-the-art methods including multiple linear regression, based on our newly compiled metrics for gene buffering predictions. Combining BaCoN with batch correction or Cholesky data whitening further boosts predictive performance. We characterized a high-confidence list of 899 buffering predictions and found that while buffering genes overall are often syntenic, buffering paralogs are on different chromosomes. BaCoN performance increases with more screens and genes considered, making it a valuable tool for gene buffering predictions from the constantly growing DepMap.


## R package

[https://github.com/billmannlab/BaCoN](https://github.com/billmannlab/BaCoN)


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

## Author
Yasir Demirtaş, Thomas Rohde, Angela Shaw, Maximilian Billmann

## License
This project is licensed under the MIT License - see the LICENSE.md file for details.
