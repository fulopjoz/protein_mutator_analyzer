# Protein Sequence Modification and Analysis

This script allows users to replace a specified percentage of amino acids in a given protein sequence with others that are part of the same chemical group. This is done while ensuring that the resulting sequence's isoelectric point (pI) and molecular weight do not deviate significantly from the original values. The script also provides various analyses of the original and modified sequences, including their molecular weight, isoelectric point, aromaticity, instability index, GRAVY (Grand Average of Hydropathy), and charge at pH 7.

This script utilizes the Biopython library, particularly the `ProteinAnalysis` module.

## Dependencies

This project has the following dependencies:

- Python 3.x
- Biopython
- Matplotlib
- NumPy
- Seaborn

To install the required dependencies, you can run the following commands in your bash shell:

```bash
pip install biopython matplotlib numpy seaborn
```



1. Run the script with:

```
./protein_sequence_analysis.py
```


1. You will be asked to input your protein sequence and the percentage of amino acids you want to replace. After first analysis you will be asked to input the sliding window size.

## Implementation Details

The script maintains a dictionary of amino acid properties, including group (hydrophilic, moderate, hydrophobic), charge (+, -, neutral), pKa values for the amino and carboxyl groups, and solubility.

To replace an amino acid, the script identifies the best replacement within the same group (i.e., with the smallest total difference in pKa values and solubility). If the replacement causes the sequence's pI or molecular weight to deviate beyond the specified thresholds from the original values, the replacement is discarded.

The script provides a color-coded alignment of the original and modified sequences, highlighting the replaced amino acids. It also calculates and displays the differences in the molecular weight, pI, aromaticity, instability index, GRAVY, and charge at pH 7 between the original and modified sequences.

The script also includes functions to generate Kyte-Doolittle and Hopp-Woods hydrophobicity plots for the sequences, providing further insight into their properties.

The main script is contained in `protein_analysis.py`.

## Notes

This script is meant for educational and research purposes and should not be used for critical applications without proper validation.