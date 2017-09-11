multifastats
============

Multifastats: Multi-Fasta Sequence Stats. Free python-based program that,
from a set of a set of 'fasta' sequences (as group or individually),
produces basic outputs and statistics (N50, %GC, molecular weight, cut
in kmers and concatenate sequences using defined spacers). And also
allows to apply some filters and select subgroups by length and/or type.

-------------------------------------------------------------------------

>Version 1.4.8, 03-Nov-2014.
Made by David Requena.
Laboratory of Bioinformatics and Molecular Biology.
Universidad Peruana Cayetano Heredia. Lima, Peru.
This program is free and will be updated in the GitHub repository:
https://github.com/lbbm-upch/multifastats

-------------------------------------------------------------------------

*** HELP INFO ***
-------------------------------------------------------------------------
Multifastats: Multi-Fasta Sequence Stats

Python script to produce basic genomic statistics and outputs for a set
of fasta sequences.

- - - - - - -
REQUIREMENTS:
- - - - - - -

* python 2.7 (x86 - 32 bits)
		https://www.python.org/download/releases/2.7.8/
* biopython
		http://biopython.org/wiki/Download
* pyreadline (to use the program in Windows and Mac)
		https://pypi.python.org/pypi/pyreadline/

Preferably, use python v2.7 32-bits because Biopython is only compatible
with x86 but not x64 (a numpy issue). So, you have to use the pyreadline
32-bits version. I made this program using biopython v1.64 and pyreadline
v2.0, but other versions should be compatible.

- - - -
USAGE:
- - - -
The script multifastats.py has to be in the same directory of your file
to be analyzed. If you want to use the program from any location, please
read the 'Optional' section.

There are two modes:

a) 'User-Interactive':
The program will ask for a file or option in each step.
To use this mode, run the script with NO arguments. The program will ask
you the name of the file you want to analyze (autocomplete allowed).

Run the script multifastats.py as follows:

In Windows:
>Just run the script with doble-click! Or from the terminal:
>>~$ multifastats.py

In Linux:
>One of the following two ways:
>>~$ python multifastats.py

>or
>>~$ ./multifastats.py

b) 'Command-Line':
Some options like file name, length cut-off, single analysis and 
pseudosequence options should be given from the beggining.

- -f (or --file): Name of the input file (including the extension). Example: inputfile.fasta
- -l (or --Lmin): Minimum sequence length to be analyzed.
- -L (or --Lmax): Maximum sequence length to be analyzed. Cutf-off value for the length of the sequences to be analyzed: Lmin =< (Sequence Length) =< Lmax. You can provide minimum, maximum or both. Only positive numbers allowed!
- -k (or --kmers):   Cut the sequences in all the possible fragments of length 'k' (k-mers).
- -p (or --pseudo):   Produces the 'pseudo-sequence' output. The sequences will be concatenated with the letter 'N' (for DNA/RNA) or 'X' (for Proteins), repeated the number of times indicated after the option '-p'
- -o (or --outsbg):   Produces an output with the subgroup of sequences analyzed.
- -s (or --single):   Add the single sequence analysis (see INFO).

Examples using all options:
>Windows Terminal:
>>~$ multifastats.py -f inputfile.fasta -l 21 -L 400 -k 9 -p 100 -o -s

>Linux Terminal:
>>~$ python multifastats.py -f inputfile.fasta -l 21 -L 400 -k 9 -p 100 -o -s

- - - - - - -
INFO OPTIONS:
- - - - - - -
You will give some arguments to get some program info:

* -h (or --help)		:	General description and usage of the program.
* -v (or --version)	:	Version and revision of the program.
* -i (or --info)		:	Information about the parameters calculated.
* -n (or --notes)	:	Notes about the current and previous versions
					of the program.

- - - - -
Optional:
- - - - -
To be able to call the program directly from the terminal in any directory
of your computer, just copy the script to '/usr/local/bin' with permisions
(chmod +xr). Example (Linux terminal):

>>~$ sudo chmod +xr multifastats.py
>>~$ sudo mv multifastats.py /usr/local/bin

Then, you could call the program from any location in your user. So, instead of writting '>>$ python multifastats.py' now you just need to write '>>$ multifastats.py'.

-------------------------------------------------------------------------

*** PARAMETERS INFORMATION ***
-------------------------------------------------------------------------

For a set of sequences, the parameters calculated are:

- Type of sequences: Type(s) of sequences (DNA, RNA or protein) in the file.
- Num of Sequences: The total number of sequences in the multifasta file.
- Min, Max and Average length: Gives the minimum, maximum and average
  sequence length in the multifasta file, respectively.
- N50 value: This parameter is calculated using the Broad Institute
  definition: https://www.broad.harvard.edu/crd/wiki/index.php/N50
- Number of residues: The total number of nucleotides (or aminoacids)
  in all the sequences in the multifasta file.
- Total %GC in file: The percentage of G, C or S (G or C) nucleotides in
  the complete set of sequences (no case sensitive).

According to the output required, some files are written in the current
working directory:

- K-mers: kmer_(FILENAME)_(TIME).fasta
- Subgroup: subgroup_(FILENAME)_(TIME).fasta
- Pseudosequence: pseudoseq_(FILENAME)_(TIME).fasta
- Single analysis: single_(FILENAME)_(TIME).csv

The last file contains, for each single sequence, the following stats:

- N: Ordinal number of the sequence in the set of sequences.
- ID: ID (header) of the sequence.
- Length: Number of nucleotides (or aminoacids) in the sequence.
- %GC: Percentage of G, C or S (G or C) nucleotides (no case sensitive)
  in the sequence.
- Mol Weight: Molecular weight of the sequence, calculated using the
  molecular_weight() function from the Biopython SeqUtils (does not works
  with ambiguous residues).

-------------------------------------------------------------------------

*** VERSION NOTES ***
-------------------------------------------------------------------------

Current version:

- Version 1.4.8 (D.R. 03-Nov-2014):
  Functions and messages improved. New option: cut sequences in k-mers.

History:

- Version 1.4.7 (D.R. 30-Oct-2014):
  Functions improved. New options: upper length cut-off. Two new
  outputs: the pseudomolecule and a fasta output of the selected
  subset of sequences. A reduction in the warnings for the sequences
  not considered.
- Version 1.4.2 (D.R. 24-Oct-2014):
  Adding autocomplete using tab in Windows and MAC OS.
- Version 1.4.1 (D.R. 23-Oct-2014):
  Minor bugs fixed.
- Version 1.4 (D.R. 23-Oct-2014):
  Adding type(s) of sequences, length cut-off and Windows compatibility.
  Now publicly available in GitHub and free, under the MIT License.
- Version 1.3 (D.R. 22-Oct-2014):
  Adding the command-line mode with -f and -s options.
  Adding full documentation (info, version, notes and help).
- Version 1.2 (D.R. 21-Oct-2014):
  Adding autocomplete using tab, analysis by single sequence and
  the output (time-csv) in the current directory.
- Version 1.1 (D.R. 21-Oct-2014):
  Adding the user-interactive mode, min, max and average length, and the
  two decimal digits format.
- Version 1.0 (D.R. 20-Oct-2014):
  The basic version, which only calculates the N50, %GC, number of
  sequences and total residues in the multifasta file.

***For the Version 1.5, I will add the following new options:
- Filter by sequence type.
- Filter by any string in sequence name.
- Top Blast-Hit for each single sequence (given certain parameters).

...Available coming soon!!!

Requests for new options are welcome!

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REALLY IMPORTANT TO DECLARE:
- - - - - - - - - - - - - - -
Some parts of the code were adapted from multiple free sources, like
biopython packages, and developed by the scientific community, like
SQUID's seqstat and Peter 'maubp' N50 code (seqanswers.com). This code
is free to use, modify and distribute under the MIT License.

*** Thanks to my lab friends for the exhaustive testing as intended
bad users.

And please, share your doubts, comments, request to add new options or
improves with us! Write to: d.requena.a@gmail.com
Thank you!
