# Dependencies:
# Folder Depend: mutstrings2.csv, fqlist.csv
# Libraries: csv, pandas, gzip, os, multiprocessing, re2
# Purpose: This script takes a list of fastq files and a list of strings to search for in the sequences.
# seqprocess function requires two fastq files, is iterable, compatible with multiprocessing
# re2 may not work on Windows systems, replaced with re in case of error

import csv
import pandas as pd
import gzip
import os
import multiprocessing as mp
try:
    import re2 as re
except ImportError:
    import re


def seqprocess(filea, fileb):
    outputfilemerge = file1.split('-')[0] + '_output.csv'
    sequences1 = []
    sequences2 = []
    filename1 = os.path.basename(filea)
    filename2 = os.path.basename(fileb)
    # Extract sequences from gzip files
    with gzip.open(filename1, 'rt') as fastq_file1, gzip.open(filename2, 'rt') as fastq_file2:
        unique_lines = set()  # Set to store unique lines
        for line_num, line in enumerate(fastq_file1):
            if line_num % 4 == 1 and line.strip() not in unique_lines:
                sequences1.append(line.strip())
                unique_lines.add(line.strip())
        for line_num, line in enumerate(fastq_file2):
            if line_num % 4 == 1 and line.strip() not in unique_lines:
                sequences2.append(line.strip())
                unique_lines.add(line.strip())
    seqall = sequences1 + sequences2
    del unique_lines
    del sequences1
    del sequences2
    print('Sequences extracted')
    # Make folder for output files
    outputfilemerge = filename1.split('-')[0] + '_output.csv'
    # Extract directory name
    output_directory = os.path.splitext(outputfilemerge)[0]
    # Create output directory if it doesn't exist
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    # Change to output directory
    os.chdir(output_directory)
    # Create output DataFrame and write headers
    print('Beginning matching')
    compiled_patterns = [re.compile(pattern) for pattern in left_strings]
    rows = []
    with mp.Pool() as pool:
        results = pool.imap(process_pattern, [(compiled_patterns, seqall)] * len(compiled_patterns))
        for matches in results:
            for sequence, match, char, pattern_str in matches:
                rows.append([pattern_str, sequence])
    # Create output DataFrame directly
    output_df = pd.DataFrame(rows, columns=['Left bases', 'Sequence'])
    # write the dataframe to a csv before processing
    output_df.to_csv(outputfilemerge.split('.')[0] + '_raw.csv', index=False)
    # Perform right matching and other operations on the DataFrame
    output_df['Match index'] = output_df['Left bases'].apply(
        lambda x: [left_strings.index(pattern) for pattern in x.split(',') if pattern in left_strings])
    output_df = output_df.explode('Match index')
    output_df['Right match'] = (output_df['Match index'].astype(bool) &
                                output_df['Sequence'].str.contains('|'.join(right_strings)))
    output_df['Right bases'] = output_df['Match index'].map(lambda x: right_strings[x] if isinstance(x, int) else '')
    # write mutation information to a new column at the same index as left strings
    output_df['Mutation information'] = output_df['Match index'].map(lambda x: string_info[x] if isinstance(x, int) else '')
    output_df = output_df[output_df['Right match']]
    output_df = output_df.drop(columns=['Match index', 'Right match'])
    output_df['Variant'] = output_df['Right bases'].str[0].fillna('')
    output_df = output_df.groupby(['Left bases', 'Variant', 'Right bases','Mutation information']).size().reset_index(name='Count total')
    # Write the final DataFrame to a CSV file
    output_df.to_csv(outputfilemerge.split('.')[0] + '_count.csv', index=False)
    os.chdir('..')
    print('Matching complete')
    print(outputfilemerge.split('.')[0] + ' complete')
    return


def process_pattern(args):
    patterns, seqall = args
    matches = []
    for sequence in seqall:
        for pattern in patterns:
            for match in pattern.finditer(sequence):
                if match.end() < len(sequence):
                    char = sequence[match.end()]
                    matches.append((sequence, match.group(), char, pattern.pattern))
                    break
    return matches

# define search strings
strings = pd.read_csv('Depend/mutstrings3.csv')
left_strings = strings.iloc[:, 0].tolist()
right_strings = strings.iloc[:, 2].tolist()
string_info = strings.iloc[:, 3].tolist()


# open the list of files to process
with open('Depend/fastq_files.csv', 'r') as file_list:
    reader = csv.reader(file_list)
    file_list = list(reader)

# apply function to each row
for row in file_list:
# get the file names from the first two entries in the row
    file1 = row[0]
    file2 = row[1]
    if __name__ == '__main__':
        p = mp.Process(target=seqprocess, args=(file1, file2))
        p.start()
        p.join()
print('All files processed')

