#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# created by fulopjoz

from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import numpy as np
import random
import math
import seaborn as sns
import random
from datetime import datetime
import os

random_seed = input("Enter the random seed (leave empty for a random seed): ")
if random_seed.strip() == "":
    random_seed = int(datetime.now().timestamp())
else:
    random_seed = int(random_seed.strip())

random.seed(random_seed)


# Define the amino acid properties based on the table
amino_acid_properties = {
    'R': {'group': 'hydrophilic', 'charge': '+', 'pKa_NH2': 9.09,  'pKa_COOH': 2.18, 'Solubility': 71.8},
    'N': {'group': 'hydrophilic', 'charge': 'N', 'pKa_NH2': 8.8,   'pKa_COOH': 2.02, 'Solubility': 2.4},
    'D': {'group': 'hydrophilic', 'charge': '-', 'pKa_NH2': 9.6,   'pKa_COOH': 1.88, 'Solubility': 0.42},
    'E': {'group': 'hydrophilic', 'charge': '-', 'pKa_NH2': 9.67,  'pKa_COOH': 2.19, 'Solubility': 0.72},
    'Q': {'group': 'hydrophilic', 'charge': 'N', 'pKa_NH2': 9.13,  'pKa_COOH': 2.17, 'Solubility': 2.6},
    'K': {'group': 'hydrophilic', 'charge': '+', 'pKa_NH2': 8.9,   'pKa_COOH': 2.2,  'Solubility': 10.28},
    'S': {'group': 'hydrophilic', 'charge': 'N', 'pKa_NH2': 9.15,  'pKa_COOH': 2.21, 'Solubility': 36.2},
    'T': {'group': 'hydrophilic', 'charge': 'N', 'pKa_NH2': 9.12,  'pKa_COOH': 2.15, 'Solubility': 'freely'},
    'C': {'group': 'moderate',    'charge': 'N', 'pKa_NH2': 10.78, 'pKa_COOH': 1.71, 'Solubility': 'freely'},
    'H': {'group': 'moderate',    'charge': '+', 'pKa_NH2': 8.97,  'pKa_COOH': 1.78, 'Solubility': 4.19},
    'M': {'group': 'moderate',    'charge': 'N', 'pKa_NH2': 9.21,  'pKa_COOH': 2.28, 'Solubility': 5.14},
    'A': {'group': 'hydrophobic', 'charge': 'N', 'pKa_NH2': 9.87,  'pKa_COOH': 2.35, 'Solubility': 15.8},
    'V': {'group': 'hydrophobic', 'charge': 'N', 'pKa_NH2': 9.72,  'pKa_COOH': 2.29, 'Solubility': 5.6},
    'G': {'group': 'hydrophobic', 'charge': 'N', 'pKa_NH2': 9.6,   'pKa_COOH': 2.34, 'Solubility': 22.5},
    'I': {'group': 'hydrophobic', 'charge': 'N', 'pKa_NH2': 9.76,  'pKa_COOH': 2.32, 'Solubility': 3.36},
    'L': {'group': 'hydrophobic', 'charge': 'N', 'pKa_NH2': 9.6,   'pKa_COOH': 2.36, 'Solubility': 2.37},
    'F': {'group': 'hydrophobic', 'charge': 'N', 'pKa_NH2': 9.24,  'pKa_COOH': 2.58, 'Solubility': 2.7},
    'P': {'group': 'hydrophobic', 'charge': 'N', 'pKa_NH2': 10.6,  'pKa_COOH': 1.99, 'Solubility': 1.54},
    'W': {'group': 'hydrophobic', 'charge': 'N', 'pKa_NH2': 9.39,  'pKa_COOH': 2.38, 'Solubility': 1.06},
    'Y': {'group': 'hydrophobic', 'charge': 'N', 'pKa_NH2': 9.11,  'pKa_COOH': 2.2,  'Solubility': 0.038}
}


def replace_amino_acids(sequence, percent_changes, pI_threshold=0.2, weight_threshold=100): # pI_threshold=0.2, weight_threshold=100
    
    # Assign a large solubility value for amino acids with 'freely' solubility
    for aa, properties in amino_acid_properties.items():
        if properties['Solubility'] == 'freely':
            properties['Solubility'] = 1000  # or other large value

    # Create a ProteinAnalysis object for the original sequence
    original_analysis = ProteinAnalysis(sequence)

    # Get the molecular weight and isoelectric point of the original sequence
    original_weight = original_analysis.molecular_weight()
    original_pI = original_analysis.isoelectric_point()
    original_aromaticity = original_analysis.aromaticity()
    original_instability_index = original_analysis.instability_index()
    original_gravy = original_analysis.gravy()
    original_charge_ph_7 = original_analysis.charge_at_pH(7)

    # Calculate the number of changes based on the percentage
    num_changes = math.ceil(len(sequence) * percent_changes / 100)

    # Initialize the modified sequence and its properties
    modified_sequence = sequence
    modified_weight = original_weight
    modified_pI = original_pI

    # Initialize the modified sequence and its properties
    modified_sequence = sequence
    modified_weight = original_weight
    modified_pI = original_pI
    modified_analysis = ProteinAnalysis(modified_sequence)
    modified_aromaticity = modified_analysis.aromaticity()
    modified_instability_index = modified_analysis.instability_index()
    modified_gravy = modified_analysis.gravy()
    modified_charge_ph_7 = modified_analysis.charge_at_pH(7)

    # Count the number of changes made
    changes = 0

    # Iterate over randomly shuffled indices
    indices = list(range(len(sequence)))
    random.shuffle(indices)
    for i in indices:
        # Check if the number of changes limit has been reached
        if changes == num_changes:
            break

        amino_acid = sequence[i]

        # Get the properties of the current amino acid
        properties = amino_acid_properties.get(amino_acid)

        # If the properties are not found, skip to the next amino acid
        if not properties:
            continue

        # Get the group, charge, pKa_NH2, pKa_COOH solubility of the amino acid
        group = properties['group']
        charge = properties['charge']
        pKa_NH2 = properties['pKa_NH2']
        pKa_COOH = properties['pKa_COOH']
        solubility = properties['Solubility']

        # Find the best replacement within the same group
        best_replacement = None
        best_distance = float('inf')  # Initialize with a large value
        for amino_acid, properties in amino_acid_properties.items():
            if properties['group'] == group and amino_acid != sequence[i]:
                distance = abs(properties['pKa_NH2'] - pKa_NH2) + abs(properties['pKa_COOH'] - pKa_COOH) + abs(properties['Solubility'] - solubility)
                if distance < best_distance:
                    best_replacement = amino_acid
                    best_distance = distance


        # If a replacement is found, calculate its properties
        if best_replacement:
            modified = modified_sequence[:i] + best_replacement + modified_sequence[i + 1:]
            modified_analysis = ProteinAnalysis(modified)
            modified_weight = modified_analysis.molecular_weight()
            modified_pI = modified_analysis.isoelectric_point()
            modified_aromaticity = modified_analysis.aromaticity()
            modified_instability_index = modified_analysis.instability_index()
            modified_gravy = modified_analysis.gravy()
            modified_charge_ph_7 = modified_analysis.charge_at_pH(7)

            # Check if the modified sequence satisfies the thresholds
            weight_diff = abs(modified_weight - original_weight)
            pI_diff = abs(modified_pI - original_pI)
            if weight_diff <= weight_threshold and pI_diff <= pI_threshold:
                modified_sequence = modified
                modified_weight = modified_analysis.molecular_weight()
                modified_pI = modified_analysis.isoelectric_point()
                changes += 1
                


    # Calculate the differences
    weight_diff = abs(modified_weight - original_weight)
    pI_diff = abs(modified_pI - original_pI)
    aromaticity_diff = abs(modified_aromaticity - original_aromaticity)
    instability_index_diff = abs(modified_instability_index - original_instability_index)
    gravy_diff = abs(modified_gravy - original_gravy)
    charge_ph_7_diff = abs(modified_charge_ph_7 - original_charge_ph_7)

    # Return the modified sequence and its properties  
    return modified_sequence, modified_weight, modified_pI, modified_aromaticity, modified_instability_index, \
           modified_gravy, modified_charge_ph_7, weight_diff, pI_diff, aromaticity_diff, instability_index_diff, \
           gravy_diff, charge_ph_7_diff, original_weight, original_pI, original_aromaticity, original_instability_index,\
           original_gravy, original_charge_ph_7


# Get user input for the sequence and the percentage of changes
sequence = input("Enter the sequence: ")
percent_changes = float(input("Enter the percentage of changes: "))


modified_sequence, modified_weight, modified_pI, modified_aromaticity, modified_instability_index, \
modified_gravy, modified_charge_ph_7, weight_diff, pI_diff, aromaticity_diff, instability_index_diff, \
gravy_diff, charge_ph_7_diff, original_weight, original_pI, original_aromaticity, original_instability_index,\
original_gravy, original_charge_ph_7 = replace_amino_acids(sequence, percent_changes)



# Print the modified sequence and its properties
print('\n')
# Print the properties and differences
print("Original Weight:   {:.2f}".format(original_weight))
print("Modified Weight:   {:.2f}".format(modified_weight))
print("Weight Difference: {:.2f}".format(weight_diff))
print('\n')

print("Original pI:   {:.2f}".format(original_pI))
print("Modified pI:   {:.2f}".format(modified_pI))
print("pI Difference: {:.2f}".format(pI_diff))
print('\n')

print("Original Aromaticity:   {:.2f}".format(original_aromaticity))
print("Modified Aromaticity:   {:.2f}".format(modified_aromaticity))
print("Aromaticity Difference: {:.2f}".format(aromaticity_diff))
print('\n')

print("Original Instability Index:   {:.2f}".format(original_instability_index))
print("Modified Instability Index:   {:.2f}".format(modified_instability_index))
print("Instability Index Difference: {:.2f}".format(instability_index_diff))
print('\n')

print("Original Gravy:   {:.2f}".format(original_gravy))
print("Modified Gravy:   {:.2f}".format(modified_gravy))
print("Gravy Difference: {:.2f}".format(gravy_diff))
print('\n')

print("Original Charge at pH 7:   {:.2f}".format(original_charge_ph_7))
print("Modified Charge at pH 7:   {:.2f}".format(modified_charge_ph_7))
print("Charge at pH 7 Difference: {:.2f}".format(charge_ph_7_diff))


def colored_diff(seq1, seq2):
    # Prepare alignment strings with mismatches highlighted
    alignment_seq1 = ''
    alignment_seq2 = ''
    diff_count = 0
    for i in range(min(len(seq1), len(seq2))):
        if seq1[i] == seq2[i]:
            alignment_seq1 += seq1[i]
            alignment_seq2 += seq2[i]
        else:
            alignment_seq1 += '\033[91m' + seq1[i] + '\033[0m'
            alignment_seq2 += '\033[91m' + seq2[i] + '\033[0m'
            diff_count += 1

    # Return the alignment strings and diff_count
    return alignment_seq1, alignment_seq2, diff_count



# Call the colored_diff function
alignment_seq1, alignment_seq2, diff_count = colored_diff(sequence, modified_sequence)

# Print the alignment and other information
print('\nAlignment:')
print(alignment_seq1)
print(alignment_seq2)
print('\n')
print(f'Length of sequences: Original: {len(sequence)}, Mutated: {len(modified_sequence)} amino acids')
print(f'Number of differing amino acids: {diff_count}')
print(f'Percent change: {diff_count/len(sequence)*100:.2f}%')
print(f'Mutated sequence: {modified_sequence}')
print('\n')
print('Kyte-Doolittle and Hopp-Woods analysis:')
print()

# Dictionary mapping amino acids to their Kyte-Doolittle scores
kyte_doolittle_scores = {
    'A': 1.80, 'C': 2.50, 'D': -3.50, 'E': -3.50, 'F': 2.80, 'G': -0.40, 'H': -3.20, 'I': 4.50,
    'K': -3.90, 'L': 3.80, 'M': 1.90, 'N': -3.50, 'P': -1.60, 'Q': -3.50, 'R': -4.50, 'S': -0.80,
    'T': -0.70, 'V': 4.20, 'W': -0.90, 'Y': -1.30
}



# Dictionary mapping amino acids to their Hopp-Woods scores
hopp_woods_scores = {
    'A': -0.50, 'C': -1.00, 'D': 3.00, 'E': 3.00, 'F': -2.50, 'G': 0.00, 'H': -0.50, 'I': -1.80,
    'K': 3.00, 'L': -1.80, 'M': -1.30, 'N': 0.20, 'P': 0.00, 'Q': 0.20, 'R': 3.00, 'S': 0.30,
    'T': -0.40, 'V': -1.50, 'W': -3.40, 'Y': -2.30
}



# User input window size
window_size = int(input("Enter the window size for analysis: "))

# Function to apply sliding window
def sliding_window(sequence, scores, window_size):
    half_window = window_size // 2
    padding_sequence = ['-'] * half_window + list(sequence) + ['-'] * half_window
    window_scores = []
    for i in range(half_window, len(padding_sequence) - half_window):
        window = padding_sequence[i - half_window : i + half_window + 1]
        window_score = sum(scores.get(aa, 0) for aa in window)
        window_scores.append(round(window_score, 2))
    return window_scores

# Kyte-Doolittle analysis with sliding window
kyte_doolittle_scores_original = sliding_window(sequence, kyte_doolittle_scores, window_size)
kyte_doolittle_scores_modified = sliding_window(modified_sequence, kyte_doolittle_scores, window_size)

# Hopp-Woods analysis with sliding window
hopp_woods_scores_original = sliding_window(sequence, hopp_woods_scores, window_size)
hopp_woods_scores_modified = sliding_window(modified_sequence, hopp_woods_scores, window_size)


print("Kyte-Doolittle scores:")
print('')

print(f'Original: \n {kyte_doolittle_scores_original}')
print('')
print(f'Modified: \n {kyte_doolittle_scores_modified}')

print('\n')

print("Hopp-Woods scores:")
print(f'Original: \n {hopp_woods_scores_original}')
print(' ')
print(f'Modified: \n {hopp_woods_scores_modified}')

print('\n')
# Highlighting function
def highlight_hydrophobic_regions(sequence, scores):
    highlighted_sequence = ""
    for aa, score in zip(sequence, scores):
        if score > 0:
            highlighted_sequence += '\033[1;31m' + aa + '\033[0m'  # Red color for hydrophobic amino acids
        else:
            highlighted_sequence += aa
    return highlighted_sequence

# Highlight hydrophobic regions
highlighted_sequence_original = highlight_hydrophobic_regions(sequence, kyte_doolittle_scores_original)
highlighted_sequence_modified = highlight_hydrophobic_regions(modified_sequence, kyte_doolittle_scores_modified)
print('Highlighted hydrophobic regions')
print()

# Print highlighted sequences
print("Original sequence:")
print(highlighted_sequence_original)
print()

print("Modified sequence:")
print(highlighted_sequence_modified)
print()
########################################
# Highlighting function
def highlight_antigenic_sites(sequence, scores):
    highlighted_sequence = ""
    for aa, score in zip(sequence, scores):
        if score > 0:
            highlighted_sequence += '\033[1;34m' + aa + '\033[0m'  # Blue color for antigenic amino acids
        else:
            highlighted_sequence += aa
    return highlighted_sequence

# Highlight antigenic sites
highlighted_sequence_original = highlight_antigenic_sites(sequence, hopp_woods_scores_original)
highlighted_sequence_modified = highlight_antigenic_sites(modified_sequence, hopp_woods_scores_modified)
print('Highlighted antigenic sites')
print()

# Print highlighted sequences
print("Original sequence:")
print(highlighted_sequence_original)
print()

print("Modified sequence:")
print(highlighted_sequence_modified)



# create a function to plot scores
def plot_scores(sequence1, scores1, sequence2, scores2, title):
    plt.figure(figsize=(12, 6))
    sns.lineplot(x=range(len(sequence1)), y=scores1, label='Original', alpha=0.7)
    sns.lineplot(x=range(len(sequence2)), y=scores2, label='Modified', alpha=0.7)
    plt.title(title)
    plt.xlabel('Amino Acid Position')
    plt.ylabel('Score')
    plt.grid(True)
    plt.legend()
    plt.show()

# plot Kyte-Doolittle analysis
plot_scores(sequence, kyte_doolittle_scores_original, modified_sequence, kyte_doolittle_scores_modified, 'Kyte-Doolittle scores')

# plot Hopp-Woods analysis
plot_scores(sequence, hopp_woods_scores_original, modified_sequence, hopp_woods_scores_modified, 'Hopp-Woods scores')

# save plot
def save_plot(sequence1, scores1, sequence2, scores2, title, filename):
    plt.figure(figsize=(12, 6))
    sns.lineplot(x=range(len(sequence1)), y=scores1, label='Original', alpha=0.7)
    sns.lineplot(x=range(len(sequence2)), y=scores2, label='Modified', alpha=0.7)
    plt.title(title)
    plt.xlabel('Amino Acid Position')
    plt.ylabel('Score')
    plt.grid(True)
    plt.legend()
    plt.savefig(filename, dpi=300)

def save_output_to_file(output_string):
    # Get the current timestamp in a suitable format for the file name
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_filename = f"output_{timestamp}.txt"

    with open(output_filename, "w") as output_file:
        output_file.write(output_string)

    print(f"Output saved to file: {output_filename}")

output_string = f'''
Original Weight:   {original_weight:.2f}
Modified Weight:   {modified_weight:.2f}
Weight Difference: {weight_diff:.2f}

Original pI:   {original_pI:.2f}
Modified pI:   {modified_pI:.2f}
pI Difference: {pI_diff:.2f}

Original Aromaticity:   {original_aromaticity:.2f}
Modified Aromaticity:   {modified_aromaticity:.2f}
Aromaticity Difference: {aromaticity_diff:.2f}

Original Instability Index:   {original_instability_index:.2f}
Modified Instability Index:   {modified_instability_index:.2f}
Instability Index Difference: {instability_index_diff:.2f}

Original Gravy:   {original_gravy:.2f}
Modified Gravy:   {modified_gravy:.2f}
Gravy Difference: {gravy_diff:.2f}

Original Charge at pH 7:   {original_charge_ph_7:.2f}
Modified Charge at pH 7:   {modified_charge_ph_7:.2f}
Charge at pH 7 Difference: {charge_ph_7_diff:.2f}

Alignment:
{alignment_seq1}
{alignment_seq2}

Length of sequences: Original: {len(sequence)}, Mutated: {len(modified_sequence)} amino acids
Number of differing amino acids: {diff_count}
Percent change: {diff_count/len(sequence)*100:.2f}%
Mutated sequence: {modified_sequence}

Kyte-Doolittle scores:

Original:
{kyte_doolittle_scores_original}

Modified:
{kyte_doolittle_scores_modified}


Hopp-Woods scores:

Original:
{hopp_woods_scores_original}

Modified:
{hopp_woods_scores_modified}

Highlighted hydrophobic regions:
Original sequence:
{highlighted_sequence_original}

Modified sequence:
{highlighted_sequence_modified}

Highlighted antigenic sites:
Original sequence:
{highlighted_sequence_original}

Modified sequence:
{highlighted_sequence_modified}
'''

# Save the output to a file with a timestamp in the name
save_output_to_file(output_string)

# save Kyte-Doolittle analysis plot
save_plot(sequence, kyte_doolittle_scores_original, modified_sequence, kyte_doolittle_scores_modified, 'Hydropath./Kyte-Doolittle scores', f'kyte_doolittle_window_{window_size}.png')

# save Hopp-Woods analysis plot
save_plot(sequence, hopp_woods_scores_original, modified_sequence, hopp_woods_scores_modified, 'Hphob./Hopp-Woods scores', f'hopp_woods_window_{window_size}.png')
