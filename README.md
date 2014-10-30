multifastats
============

Multifastats: Multi-Fasta Sequence Stats. Free program based in a python
script, to calculate basic statistics and outputs and allows to apply some
filters (N50, %GC, molecular weight, pseudosequences and subgroups for
length and type filters) for a set of 'fasta' sequences (and optionally for
each single sequence too).

-------------------------------------------------------------------------------

>Version 1.4.7, 30-Oct-2014.
Made by David Requena. Laboratory of Bioinformatics and Molecular Biology.
Universidad Peruana Cayetano Heredia. Lima, Peru.
This code will be updated and freely available in GitHub:
https://github.com/lbbm-upch/multifastats

-------------------------------------------------------------------------------

*** HELP INFO ***
-------------------------------------------------------------------------------
Multifastats: Multi-Fasta Sequence Stats

Python script to calculate basic genomic statistics for a set of sequences,
and optionally for each single sequence, too.

- - - - - - -
REQUIREMENTS:
- - - - - - -
To run this program, python 2.7 and the following packages are required:
* python 2.7 (x86 - 32 bits) https://www.python.org/download/releases/2.7.8/
* biopython http://biopython.org/wiki/Download
* pyreadline (to use the program in Windows and Mac) https://pypi.python.org/pypi/pyreadline/

Preferably, use python v2.7 32-bits because Biopython is only compatible with
x86 but not x64 (a numpy issue). So, you have to use the pyreadline 32-bits
version. I made this program with biopython v1.64 and pyreadline v2.0, but
other versions will be compatible.

- - - -
USAGE:
- - - -
The script multifastats.py has to be in the same directory of your file to analyze.
If you want to use the program from any location, please read the 'Optional' section.

There are two modes:

a) 'User-Interactive':
The program will ask for a file or option in each step.
To use in this mode, run the script with NO arguments. The program will ask you the name
of the file you want to analyze (autocomplete allowed). To use the program in this mode,
just run the script multifastats.py as follows:

>In Windows:
>>Just run the script with doble-click!
>>Or from the terminal:   ~$ multifastats.py

>In Linux:
>>One of the following two ways:
>>~$ python multifastats.py
>or
>>~$ ./multifastats.py

b) 'Command-Line':
Some options like file name, length cut-off, single analysis and pseudosequence output
would be given from the beggining.

- -f (or --file): Allows to give an input file name (including the extension). As example: inputfile.fasta
- -l (or --Lmin): Minimum sequence length to be analyzed
- -L (or --Lmax): Maximum sequence length to be analyzed. Add a cutf-off value for the length of the sequences to be analyzed: Lmin =< (Sequence Length) =< Lmax. You can provide minimum, maximum or both cut-off values. Only positive numbers allowed!
- -o (or --outsbg): Produces an output with the subgroup of sequences analyzed
- -s (or --single):   Add the single sequence analysis (see INFO)
- -p (or --pseudo):   Produces the 'pseudo-sequence' output. The sequences will be concatenated with the letter 'N' (for DNA/RNA) or 'X' (for Proteins) repeated the number of times indicated after -p

As example, one command line with all options will be:

>In Windows Terminal:
>>~$ multifastats.py -f inputfile.fasta -l 21 -L 400 -p 100 -s

>In Linux Terminal:
>>~$ python multifastats.py -f inputfile.fasta -l 21 -L 400 -p 100 -s

- - - - - - -
INFO OPTIONS:
- - - - - - -
You will give some arguments to get some program info:

* -h (or --help)	 :	General description and usage of the program.
* -v (or --version)	 :	Version and revision of the program.
* -i (or --info)	 :          Information about the parameters calculated and outputs.
* -n (or --notes)	 :	Notes about the current and previous versions program.

- - - - -
Optional:
- - - - -
Maybe you will have more fun if you copy the script to '/usr/local/bin' 
and give some permisions (chmod +xr). This will allow you to call it and use
directly from the terminal in any directory of your computer and allow to
autocomplete the name of any file (not only python scripts).
You will do that with the following commands in the linux terminal:

>>~$ sudo chmod +xr multifastats.py

>>~$ sudo mv multifastats.py /usr/local/bin

Then, you will call the program from any location in your user. So, instead to
write '~$ python multifastats.py' you just need to write '~$ multifastats.py'

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*** PARAMETERS INFORMATION ***
-------------------------------------------------------------------------------

For sets of sequences, the parameters calculated are:

- Type of sequences: Indicates the type(s) of sequences in the file (DNA,
  RNA or protein).
- Num of Sequences: Count the total number of sequences in the multifasta file.
- Min, Max and Average length: Gives the minimum, maximum and average sequence
  length in the multifasta file, respectively.
- N50 value: This parameter is calculated using the Broad Institute definition:
  https://www.broad.harvard.edu/crd/wiki/index.php/N50
- Number of residues: Count the total number of nucleotides (or aminoacids) in
  all the sequences in the multifasta file.
- Total %GC in file: This parameter calculates the percentage of G, C or S (G
  or C) nucleotides (no case sensitive) in the complete set of sequences.

For the output of the subgroup of sequences, a .csv file is write in the
current working directory: subgroup_(FILENAME)_(TIME).csv

For the pseudosequence output, a .csv file is write in the current
working directory: pseudoseq_(FILENAME)_(TIME).csv

For the single sequence analysis, a .csv file is write in the current
working directory: sganl_(FILENAME)_(TIME).csv
This contains stats of each single sequence:

- N: Number of the sequence in the set of sequences.
- ID: ID (header) of the sequence.
- Type: Sequence type (DNA, RNA or protein).
- Length: Number of nucleotides (or aminoacids) in the sequence.
- %GC: Percentage of G, C or S (G or C) nucleotides (no case sensitive) in the
  sequence.
- Mol Weight: Molecular weight of the sequence, calculated according the
  molecular_weight() function from Biopython SeqUtils (does not works with
  ambiguous residues). http://biopython.org/DIST/docs/api/Bio.SeqUtils-module.html#molecular_weight

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*** VERSION NOTES ***
-------------------------------------------------------------------------------
History:
- Version 1.4.7 (D.R. 30-Oct-2014):
  Functions improved. New options: upper length cut-off, a fasta output of the selected
  subset of sequences, the pseudomolecule output and a reduce in the warnings for
  sequences not considered.
- Version 1.4.2 (D.R. 24-Oct-2014):
  Adding autocomplete using tab in Windows and MAC OS.
- Version 1.4.1 (D.R. 23-Oct-2014):
  Including the revisions after testing of many users.
- Version 1.4 (D.R. 23-Oct-2014):
  Adding length cut-off, indicate type(s) of sequences and Windows compatibility.
  Public available in GitHub, under the MIT License.
- Version 1.3 (D.R. 22-Oct-2014):
  Adding command-line mode with -f, -s and info options.
  Adding version, notes, help and full documentation.
- Version 1.2 (D.R. 21-Oct-2014):
  Adding autocomplete using tab, single sequences analysis and time-csv output
  in current directory.
- Version 1.1 (D.R. 21-Oct-2014):
  Adding user-interactive mode, min, max, average length and 0.00 format.
- Version 1.0 (D.R. 20-Oct-2014):
  The basic version, just calculate the N50, %GC, number of sequences and total
  residues in the multifasta by a given file name inside the script.

***For the Version 1.5, I'm implementing the following new options:
- Filter by sequence type.
- Filter by any string in sequence name.
- Top Blast-Hit for each single sequence (given certain parameters).

...Available coming soon!!!

Requests for new options are welcome!

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REALLY IMPORTANT TO DECLARE:
- - - - - - - - - - - - - - -
Some parts of the code was adapted from multiple sources freely available
developed by the scientific community, like SQUID's seqstat, Peter 'maubp' N50
(seqanswers.com) and biopython packages. This code is free to use, modify and
distribute under the MIT License.

***Thanks to my lab friends by the exhaustively testing as intended bad users.

And please, share your doubts, comments, request to add new options or
improves with us! Write to: david.requena.a@upch.pe

Thank you!
