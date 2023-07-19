# Protein Sequence Modification and Analysis

This Python script is a bioinformatics tool that performs sequence mutation and analysis on a given protein sequence. It uses various parameters related to amino acid properties, such as hydrophilicity, charge, pKa values, and solubility, to perform mutations while trying to maintain important characteristics of the original sequence. Additionally, it analyzes the impact of these mutations on properties such as molecular weight, isoelectric point (pI), aromaticity, instability index, GRAVY (grand average of hydropathy), and charge at pH 7.

This script utilizes the Biopython library, particularly the `ProteinAnalysis` module.

## Requirements

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

## Usage

1. Clone or download the repository.
2. Run the script using Python:

```
./protein_sequence_analysis.py
```


1. The script will prompt you to enter a protein sequence and the percentage of changes you want to introduce to the sequence.
2. You can also enter a random seed for reproducibility, but leaving it empty will generate a random seed based on the current timestamp.
3. The script will perform sequence mutation and analysis. It will display the original and modified properties, differences, and the modified sequence. Additionally, it will show the alignment and the percentage change in amino acids.
4. The script will generate plots for the Kyte-Doolittle and Hopp-Woods scores of the original and modified sequences.
5. All results and plots will be saved to a file named output_<timestamp>.txt, where <timestamp> is the current date and time when the file was created.

## Implementation Details

The script maintains a dictionary of amino acid properties, including group (hydrophilic, moderate, hydrophobic), charge (+, -, neutral), pKa values for the amino and carboxyl groups, and solubility.

To replace an amino acid, the script identifies the best replacement within the same group (i.e., with the smallest total difference in pKa values and solubility). If the replacement causes the sequence's pI or molecular weight to deviate beyond the specified thresholds from the original values, the replacement is discarded.

The script provides a color-coded alignment of the original and modified sequences, highlighting the replaced amino acids. It also calculates and displays the differences in the molecular weight, pI, aromaticity, instability index, GRAVY, and charge at pH 7 between the original and modified sequences.

The script also includes functions to generate Kyte-Doolittle and Hopp-Woods hydrophobicity plots for the sequences, providing further insight into their properties.

The main script is contained in `protein_analysis.py`.

## Output

The script will produce a file named output_<timestamp>.txt, which will contain the following information:

* Original and modified weights
* Weight difference
* Original and modified pI
* pI difference
* Original and modified aromaticity
* Aromaticity difference
* Original and modified instability index
* Instability index difference
* Original and modified GRAVY
* GRAVY difference
* Original and modified charge at pH 7
* Charge at pH 7 difference
* Alignment of the original and modified sequences
* Length of the original and modified sequences
* Number of differing amino acids
* Percent change in amino acids
* Kyte-Doolittle scores of the original and modified sequences
* Hopp-Woods scores of the original and modified sequences

The script will also generate two plots:

* A plot showing the Kyte-Doolittle scores of the original and modified sequences
* A plot showing the Hopp-Woods scores of the original and modified sequences

These plots will be saved as kyte_doolittle_window_<window_size>.png and hopp_woods_window_<window_size>.png, respectively.

Please note that the output will be specific to the input sequence and percentage of changes provided during runtime.

## Disclaimer

This tool is for educational and research purposes only. Use it responsibly and always verify the results. The authors are not responsible for any consequences arising from the use of this tool.

## License

This project is licensed under the MIT License. See the LICENSE file for more details.