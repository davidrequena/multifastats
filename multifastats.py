#!/usr/bin/python
# -*- coding: utf-8 -*-
import os, sys, platform
#=========================================================================#
#To update in every release:
vers='Version 1.4.7, 30-Oct-2014.'
versnote='''- Version 1.4.7 (D.R. 30-Oct-2014):
  Functions improved. New options: upper length cut-off, a fasta
  output of the selected subset of sequences, the pseudomolecule output and
  a reduce in the warnings for the sequences not considered.'''
#=========================================================================#
dochpini='''
===============================================================================
*** HELP INFO ***
===============================================================================
Multifastats: Multi-Fasta Sequence Stats
-----------------------------------------
Python script to calculate basic genomic statistics for a set of sequences,
and optionally for each single sequence, too.
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REQUIREMENTS:
- - - - - - -
To run this program, python 2.7 and the following packages are required:
* python 2.7 (x86 - 32 bits)
  https://www.python.org/download/releases/2.7.8/
* biopython
  http://biopython.org/wiki/Download
* pyreadline (to use the program in Windows and Mac)
  https://pypi.python.org/pypi/pyreadline/
Preferably, use python v2.7 32-bits because Biopython is only compatible with
x86 but not x64 (a numpy issue). So, you have to use the pyreadline 32-bits
version. I made this program with biopython v1.64 and pyreadline v2.0, but
other versions will be compatible.

USAGE:
- - -
The script multifastats.py has to be in the same directory of your file to
analyze. If you want to use the program from any location, please read the
'Optional' section.

There are two modes:
1) 'User-Interactive': The program will ask for a file or option in each step.
To use in the program in this mode, run the script multifastats.py as follows:
  -In Windows:\tJust run the script with doble-click!
  \t\tOr from the terminal:    >>~$ multifastats.py
  -In Linux:\tOne of the following two ways:
  \t\t>>~$ python multifastats.py   or   >>~$ ./multifastats.py

2) 'Command-Line': Some options like file name, length cut-off, single
analysis and pseudosequence options would be given from the beggining.

-f (or --file): \tAllows to give an input file name (including the
\t\t\textension). As example: inputfile.fasta
-l (or --Lmin): \tMinimum sequence length to be analyzed
-L (or --Lmax): \tMaximum sequence length to be analyzed
\t\t\tAdd a cutf-off value for the length of the sequences
\t\t\tto be analyzed: Lmin =< (Sequence Length) =< Lmax.
\t\t\tYou can provide minimum, maximum or both cut-off
\t\t\tvalues. Only positive numbers allowed!
-o (or --outsbg):\tProduces an output with the subgroup of sequences
\t\t\tanalyzed
-s (or --single):\tAdd the single sequence analysis
-p (or --pseudo):\tProduces the 'pseudo-sequence' output. The sequences will
\t\t\tbe concatenated with the letter 'N' (for DNA/RNA)
\t\t\tor 'X' (for Proteins) repeated the number of times
\t\t\tindicated after -p

As example, one command line with all options will be:
In Windows Terminal:
  >>~$ multifastats.py -f inputfile.fasta -l 21 -L 400 -p 100 -s
In Linux Terminal:
  >>~$ python multifastats.py -f inputfile.fasta -l 21 -L 400 -p 100 -s

INFO OPTIONS:
- - - - - - -
You will give some arguments to get some program info:
-h (or --help)\t :\tGeneral description and usage of the program.
-v (or --version):\tVersion (and revision) of the program.
-i (or --info)\t :\tInformation about the parameters calculated and outputs.
-n (or --notes)\t :\tNotes about the current and previous versions of the program.

Optional:
- - - - -
Maybe you will have more fun if you copy the script to '/usr/local/bin' 
and give some permisions (chmod +xr). This will allow you to call it and use
directly from the terminal in any directory of your computer and allow to
autocomplete the name of any file (not only python scripts).
You will do that with the following commands in the unix terminal:
>>~$ sudo chmod +xr multifastats.py
>>~$ sudo mv multifastats.py /usr/local/bin
Then, you will call the program from any location in your user, directly from
the command line. So, instead to write '>>~$ python multifastats.py' you just
need to write '>>~$ multifastats.py'
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'''
dochpmid=' Made by David Requena. Laboratory of Bioinformatics'
dochpend='''and Molecular Biology, Universidad Peruana Cayetano Heredia. Lima, Peru.
This code will be updated and freely available in GitHub:
https://github.com/lbbm-upch/multifastats

Please, share your doubts, comments, request to add new options or improves
with us! Write to: david.requena.a@upch.pe
Thank you!
===============================================================================
'''
docvrini='''
   Multifastats: Multi-Fasta Sequence Stats
----------------------------------------------'''
docvrend=' Revised by David Requena. LBBM, UPCH. Lima, Peru.\n'
docntini='''
===============================================================================
*** VERSION NOTES ***
===============================================================================
Multifastats: Multi-Fasta Sequence Stats
----------------------------------------
History:
'''
docntend='''- Version 1.4.2 (D.R. 24-Oct-2014):
  Adding autocomplete using tab in Windows and MAC OS.
- Version 1.4.1 (D.R. 23-Oct-2014):
  Including the revisions after testing of many users.
- Version 1.4 (D.R. 23-Oct-2014):
  Adding length cut-off, indicate type(s) of sequences and Windows
  compatibility. Public release in GitHub (under MIT License).
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
- Upper length cut-off.
- Filter by any string in sequence name.
- A fasta output of the selected subset of sequences.
- The pseudomolecule output.
- Top Blast-Hit for each single sequence (given certain parameters).

...Available coming soon!!!
Requests for new options are welcome!

 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Really important to declare:
----------------------------
Some parts of the code was adapted from multiple sources freely available
developed by the scientific community, like SQUID's seqstat, Peter 'maubp' N50
(seqanswers.com) and biopython packages. This code is free to use, modify and
distribute under the MIT License.
***Thanks to my lab friends by the exhaustively testing as intended bad users.

And please, share your doubts, comments, request to add new options or
improves with us! Write to: david.requena.a@upch.pe
Thank you!
===============================================================================
'''
docinfo='''
===============================================================================
*** PARAMETERS INFORMATION ***
===============================================================================
Multifastats: Multi-Fasta Sequence Stats
----------------------------------------
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

For the single sequence analysis, a .csv file is write in the current working
directory: sganl_(FILENAME)_(TIME).csv
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
===============================================================================
'''
#=========================================================================#
usersys=platform.system()
def exitval():
    if usersys=='Windows':
        os.system('pause')
        exit()
    else:
        exit()
#=========================================================================#
maininp=0 #This variable is the state of the manual input requirement for the multifasta analysis (1=YES/0=NO).
mainout=0 #This variable is the state of the type of output for the multifasta analysis (1=YES/0=NO).
testsbg=1 #This variable is the state of produce or not the output of the subgroup of sequences analyzed (1=YES/0=NO).
testsga=1 #This variable is the state of do or not do the single analysis (1=YES/0=NO).
testps=1 #This variable is the state of produce or not the pseudosequence output (1=YES/0=NO).
sbginp=0 #This variable is the state of the manual input requirement for the output of the subgroup of sequences analyzed (1=YES/0=NO).
sgainp=0 #This variable is the state of the manual input requirement for the single analysis (1=YES/0=NO).
psinp=0 #This variable is the state of the manual input requirement for the pseudopsequence (1=YES/0=NO).
lctinp=0 #This variable is the state of the manual input requirement for the length cut-off (1=YES/0=NO).
nrep=0 #This variable is the number of times to repeat the ambiguous residue in the pseudosequence.
minlen=0 #This variable is the minimum length to analyze (0 means no cut-off).
maxlen=0 #This variable is the maximum length to analyze (0 means no cut-off).
#=========================================================================#
args=[x.lower() for x in sys.argv] #Lowcase list of arguments received
hlist=['h','-h','help','-help','--help']
vlist=['v','-v','version','-version','--version']
ilist=['i','-i','info','-info','--info']
nlist=['n','-n','notes','-notes','--notes']
#=========================================================================#
if len(sys.argv)==1: #This is the case when the script NO receive arguments (***THE "USER-INTERACTIVE" SCENARIO***).
    print '- '*29+'-'
    print 'Multifastats: Multi-Fasta Sequence Stats'
    print '-'*41
    print 'Show basic genomic statistics for a given set of sequences\nand optionally for each single sequence, too.'
    print 'Made by D.R. LBBM, UPCH. '+vers
    mainout=1
    maininp=1
    lctinp=1
    psinp=1
    sbginp=1
    sgainp=1
    pass
#=========================================================================#
elif len(sys.argv)==2: #This is the case when the script receive arguments asking for information.
    if args[1] in hlist:
        print dochpini
        print vers+dochpmid
        print dochpend
        exitval()
    elif args[1] in vlist:
        print docvrini
        print vers+docvrend
        exitval()
    elif args[1] in ilist:
        print docinfo
        exitval()
    elif args[1] in nlist:
        print docntini
        print versnote
        print docntend
        exitval()
    else: #We will add more options for info in the future versions.
        print "Incorrect way to give the arguments. See 'help' with '-h' or '--help' option."
        exitval()
#=========================================================================#
elif (len(sys.argv)>=3) and (('-f' or '--file') in args): #This is the case when the script receive arguments (***COMMAND-LINE SCENARIO***) and the file option correctly.
    try:
        indf=args.index('-f')
    except:
        indf=args.index('--file')
    currDir=os.listdir(os.getcwd())
    if sys.argv[indf+1] in currDir: #If we receive the file name correctly.
        filename=sys.argv[indf+1] #Use this file to continue.
    else: #If receive the file name incorrectly.
        print 'Incorrect file name given.'
        exitval()
    if (('-l' in sys.argv) or ('--lmin' in args)):
        try:
            indlmin=sys.argv.index('-l')
        except:
            indlmin=args.index('--lmin')
        try:
            minlen=float(sys.argv[indlmin+1])
        except:
            print "Incorrect value provided for minimum length. See 'help' with '-h' or '--help' option."
            exitval()
    if (('-L' in sys.argv) or ('--lmax' in args)):
        try:
            indlmax=sys.argv.index('-L')
        except:
            indlmax=args.index('--lmax')
        try:
            maxlen=float(sys.argv[indlmax+1])
        except:
            print "Incorrect value provided for maximum length. See 'help' with '-h' or '--help' option."
            exitval()
    if (('-p' or '--pseudo') in args):
        try:
            indps=sys.argv.index('-p')
        except:
            indps=args.index('--pseudo')
        try:
            nrep=int(sys.argv[indps+1])
        except:
            print "Incorrect value provided for the number of repeats. See 'help' with '-h' or '--help' option."
            exitval()
    testsga=0 if (('-s' or '--single') not in args) else 1 #This means that we will do the single analysis if -s option if given.
    testsbg=0 if (('-o' or '--outsbg') not in args) else 1 #This means that we will do the single analysis if -s option if given.
    testps=0 if (('-p' or '--pseudo') not in args) else 1 #This means that we will make the pseudosequence output.
else: #Any other case
    print "Incorrect way to give the arguments. See 'help' with '-h' or '--help' option."
    exitval()
#=========================================================================#
def complete(text, state):
    for cmd in os.listdir(os.getcwd()):
        if cmd.startswith(text):
            if not state:
                return cmd
            else:
                state -= 1
#=========================================================================#
if maininp: #This means that I need a manual input: THE "USER-INTERACTIVE" SCENARIO.
    #Setting the autocomplete capability of the script
    try:
        import readline #If you have the readline package (Linux), autocompletion enabled.
        readline.parse_and_bind("tab: complete")
        readline.set_completer(complete)
    except ImportError: #If you have not the readline package, autocompletion disabled.
        try:
            import pyreadline as readline
            readline.parse_and_bind("tab: complete")
            readline.set_completer(complete)
        except ImportError:
            print "\n"+"*"*79+"\nWARNING!: The readline or pyreadline (32bits x86) packages are necesary to run.\nYou will download it from: https://pypi.python.org/pypi/pyreadline/\nSee the 'help' file below."+"*"*79
            print dochpini
            print vers+dochpmid
            print dochpend
            exitval()
    while maininp: #If we are in the user-interactive mode to input multifasta file, this continue asking for an existing file name in the directory as input file:
        print "-"*79+"\nTO RUN: Just write the name of the file you want to analyze."
        print "See 'help' giving the options '-h' or '--help' below, or from\nthe command-line as follows: >>~$ python multifastats.py -h"
        if usersys in ('Windows','Linux'):
            print "(autocomplete allowed using 'TAB' key)"
        else:
            print "You are in "+usersys+". Please, do not use 'TAB' key. Autocomplete not implemented."
        filename=raw_input("File name or Option: ") #Asking for a filename or option, to store in the variable filename.
        options=hlist+vlist+ilist+nlist
        if filename.lower() in options: #If the input received in filename  is an option instead the name of a file.
            maininp=0 #This means that we get a correct file name or option, so we will use this to continue.
            lctinp=0
            if filename.lower() in hlist:
                print dochpini
                print vers+dochpmid
                print dochpend
            elif filename.lower() in vlist:
                print docvrini
                print vers+docvrend
            elif filename.lower() in ilist:
                print docinfo
            elif filename.lower() in nlist:
                print docntini
                print versnote
                print docntend
            else: #We will add more options for info in the future versions.
                print "\nIncorrect way to give the arguments. See 'help' with '-h' or '--help' option."
            exitval()
        elif filename in os.listdir(os.getcwd()): #If the input received in filename  is the name of a file.
            maininp=0 #This means that we get a correct file name or option, so we will use this to continue.
            pass
        else:
            print "\nBad file name or option, try again or press 'Ctrl+C' to exit."
            pass
#=========================================================================#
if lctinp:
    clinp=0 #This variable allows to verify if a valid value for L is given.
    print '- '*29+'-'
    print "Do you want to analyze all the sequences? Or just a subset\nof them between a minimum and a maximum length?"
    ctlnopt=int(raw_input("Please, choose one option: (1, 2, 3 or 4)\n\n1) Analyze all the sequences\n2) A subset above a minimum sequence length (Lmin)\n3) A subset below a maximum sequence length (Lmax)\n4) A subset between a Lmin and Lmax values\n\nEnter your option:"))
    while (ctlnopt not in (1, 2, 3, 4)):
        ctlnopt=int(raw_input("Bad input given. Please, choose one option: (1, 2, 3 or 4)\n\n1) Analyze all the sequences\n2) A subset above a minimum sequence length (Lmin)\n3) A subset below a maximum sequence length (Lmax)\n4) A subset between a Lmin and Lmax values\n\nEnter your option:"))
    if ctlnopt==1:
        pass
    elif ctlnopt==2:
        while clinp==0:
            tochklc=raw_input("Please, give a positive value for the minimum sequence length to analyze\n(or 0 to not set a minimum value): ")
            try:
                if (float(tochklc)>=0):
                    clinp=1 #Correct value of l given, continue.
                    minlen=float(tochklc)
                else:
                    print "Incorrect value provided for minimum length."
                    pass
            except ValueError:
                print "A positive number is required as value for length."
                pass
    elif ctlnopt==3:
        while clinp==0:
            tochklc=raw_input("Please, give a positive value for the maximum sequence length to analyze\n(or 0 to not set a maximum value): ")
            try:
                if (float(tochklc)>=0):
                    clinp=1 #Correct value of L given, continue.
                    maxlen=float(tochklc)
                else:
                    print "Incorrect value provided for maximum length."
                    pass
            except ValueError:
                print "A positive number is required as a value for length."
                pass
    elif ctlnopt==4:
        while clinp==0:
            tochklc=raw_input("Please, give a positive value for the minimum sequence length to analyze\n(or 0 to not set a minimum value): ")
            try:
                if (float(tochklc)>=0):
                    clinp=1 #Correct value of l given, continue.
                    minlen=float(tochklc)
                else:
                    print "Incorrect value provided for minimum length."
                    pass
            except ValueError:
                print "A positive number is required as a value for length."
                pass
        clinp=0 #Give 0 to clinp to evaluate a new input (now for Lmax).
        while clinp==0:
            tochklc=raw_input("Please, give a positive value for the maximum sequence length to analyze\n(or 0 to not set a maximum value): ")
            try:
                if (float(tochklc)>=0):
                    clinp=1 #Correct value of L given, continue.
                    maxlen=float(tochklc)
                else:
                    print "Incorrect value provided for maximum length."
                    pass
            except ValueError:
                print "A positive number is required as a value for length."
                pass
#=========================================================================#
#Import some extra libraries to use in the next procedures.
try:
    from Bio import SeqIO, Alphabet
    from Bio.Seq import Seq
    from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA, IUPACAmbiguousRNA, ExtendedIUPACProtein
    from Bio.Alphabet import ProteinAlphabet, _verify_alphabet
except ImportError:
    print "The package 'biopython' is necesary to run.\nSee 'help' with '-h' or '--help' option."
    exitval()
#=========================================================================#
def seqcategory(oneseq):
    seqtype=''
    seqDNA=Seq(oneseq,IUPACAmbiguousDNA()) #Produce a sequence using the string received and the DNA alphabet.
    seqRNA=Seq(oneseq,IUPACAmbiguousRNA()) #Produce a sequence using the string received and the RNA alphabet.
    seqProt=Seq(oneseq,ExtendedIUPACProtein()) #Produce a sequence using the string received and the protein alphabet.
    if Alphabet._verify_alphabet(seqDNA): #Verify if is a DNA sequence.
        seqtype='DNA'
    elif Alphabet._verify_alphabet(seqRNA): #Verify if is a RNA sequence.
        seqtype='RNA'
    else:
        if Alphabet._verify_alphabet(seqProt): #Verify if is a protein sequence.
            seqtype='protein'
        else:
            seqtype='noseq' #If any, is not a valid sequence.
    return seqtype
#=========================================================================#
def lenfilter(strseq,strid,lmin,lmax,warn):
    typeseq=seqcategory(strseq)
    if (lmax==0 and len(strseq)>=lmin):
        if typeseq in ('DNA', 'RNA', 'protein'): #If valid format, use this sequence and calculate parameters.
            return 1, typeseq
        elif typeseq=='noseq': #If invalid format, ignore this sequence.
            if warn==1:
                print 'WARNING!: Sequence '+strid+' has invalid format and not taken in account.'
            return 0, ''
        else:
            return 0, ''
    elif len(strseq) in range(int(lmin),int(lmax)+1):
        if typeseq in ('DNA', 'RNA', 'protein'): #If valid format, use this sequence and calculate parameters.
            return 1, typeseq
        elif typeseq=='noseq': #If invalid format, ignore this sequence.
            if warn==1:
                print 'WARNING!: Sequence '+strid+' has invalid format and not taken in account.'
            return 0, ''
        else:
            return 0, ''
    else:
        return 0, ''
#=========================================================================#
def N50(lenlist): # N50 calculation, based on the Broad Institute definition: https://www.broad.harvard.edu/crd/wiki/index.php/N50
    lenlist.sort()
    wghtlist = []
    for i in lenlist:
        wghtlist += [i]*i
    # Take the median of the weighted list of lengths:
    if len(wghtlist) % 2 == 0:
        return float(wghtlist[len(wghtlist)/2] + wghtlist[len(wghtlist)/2-1])/2
    else:
        return wghtlist[len(wghtlist)/2]
#=========================================================================#
def pseudoseq(filenm,numrep,lmin,lmax,ctime,warn):
    from Bio.SeqRecord import SeqRecord
    pseudo_handle = open("pseudoseq_"+filenm.replace(".","")+ctime+".fasta", "w")
    manyseqs = open(filenm,"rU")
    if lmax==0 and lmin>0: #Setting a specific description for each case of length restriction.
        dscrp="FILE: "+filenm+" Sequences from "+str(lmin)+" residues length."
    elif lmax>=lmin>0:
        dscrp="FILE: "+filenm+" Sequences from "+str(lmin)+" to "+str(lmax)+" residues length."
    else:
        dscrp="FILE: "+filenm #The description of the pseudosequence.
    newseq=SeqRecord(Seq(''), id="pseudoseq_"+ctime, description=dscrp) #Creating an empty sequence to append each sequence which pass the length filter.
    for eachseq in SeqIO.parse(manyseqs, "fasta"):
        evalen=lenfilter(str(eachseq.seq),str(eachseq.id),lmin,lmax,warn)
        if evalen[0] and (evalen[1] in ['DNA','RNA']): #If DNA or RNA, concatenate with Ns.
            newseq.seq+=eachseq.seq+numrep*'N'
        elif evalen[0] and evalen[1]=='protein': #If protein, concatenate with Xs.
            newseq.seq+=eachseq.seq+numrep*'X'
        else:
            pass
    newseq.seq=newseq.seq[:-numrep] #This step is to cut the last Ns or Xs added by the for-loop.
    SeqIO.write(newseq, pseudo_handle, "fasta")
    pseudo_handle.close()
#=========================================================================#
def subgroupseq(filenm,lmin,lmax,ctime,warn):
    sbgrp_handle = open("subgroup_"+filenm.replace(".","")+ctime+".fasta", "w") #Creating a new file to write the sequences.
    multifst = open(filenm,"rU")
    for eachseq in SeqIO.parse(multifst, "fasta"):
        evalen=lenfilter(str(eachseq.seq),str(eachseq.id),lmin,lmax,warn)
        if evalen[0]:
            SeqIO.write(eachseq, sbgrp_handle, "fasta") #If the sequence passes the filter is writen in the file.
    multifst.close()
    sbgrp_handle.close()
#=========================================================================#
def singanlys(filenm,lmin,lmax,ctime,warn):
    from Bio.SeqUtils import GC, molecular_weight
    outpbyseqcsv = open("sganl_"+filenm.replace(".","")+"_"+ctime+".csv","w") #Creating the file to write the output
    outpbyseqcsv.write('N,ID,Type,Length,%GC,Mol Weight\n') #Headers line.
    idsnocom=[] #To store the ID of each seq without commas.
    molwght=[] #To store the molecular weight of each sequence.
    multifs = open(filenm,"rU")
    for indsq in SeqIO.parse(multifs, "fasta"):
        evalen=lenfilter(str(indsq.seq),str(indsq.id),lmin,lmax,warn)
        if evalen[0]: #Use the sequence if pass the filter.
            newid = indsq.id.replace(",", " ")
            idsnocom.append(newid)
            try: #Try calculate the molecular weight if possible, else pass.
                molwght.append(molecular_weight(str(indsq.seq),evalen[1]))
            except:
                molwght.append(0)
        else: #If does not pass the filter, ignore.
            pass
    count=0
    for i in range(len(contgslen)): #This is a external resource of the function. If you want to make the function independent from the rest of the code, provide the contgslen list as input.
        count+=1
        newlinecsv=str(count)+','+str(idsnocom[i])+','+typesofseqs[i]+','+str(contgslen[i])+','+"{0:.2f}".format(gcs[i])+','+"{0:.2f}".format(molwght[i])+'\n' #If we want to use the function independently, we need to change this line to: #newlinecsv=str(count)+','+str(idsnocom[i])+','+str(listcntgslen[i])+','+"{0:.2f}".format(listgcs[i])+','+"{0:.2f}".format(molwght[i])+'\n'
        #print newlinecsv #I left this line because maybe someone wants to see the output in the python interpreter.
        outpbyseqcsv.write(newlinecsv)
        multifs.close()
    outpbyseqcsv.close()
#=========================================================================#
sqfs = open(filename,"rU")
numGC=0
currGC=0
contgslen=[] #To store contigs lengths.
gcs=[] #To store %GC values.
typesofseqs=[] #To store the type of sequences in the file.
wrng=1 #Show the warning message for the sequences in the range of lengths not taken in account because they have bad format.
for contg in SeqIO.parse(sqfs, "fasta"):  #If is in fasta format, calculate some parameters for sequences with valid format.
    evalen=lenfilter(str(contg.seq),str(contg.id),minlen,maxlen,wrng)
    if evalen[0]: #Use the sequence if pass the filter.
        typesofseqs.append(evalen[1])
        currGC=sum(contg.seq.count(x) for x in ['G', 'C', 'g', 'c', 'S', 's'])
        numGC+=currGC
        currlen=len(contg.seq)
        contgslen.append(currlen)
        gcs.append(currGC*100.0/currlen)
    else: #If does not pass the filter, ignore.
        pass
#=========================================================================#
if len(contgslen)>0:
    types=sorted([str(x) for x in set(typesofseqs)])
    sqfs.close()
    numseqs=len(contgslen)
    totlen=sum(contgslen)
    avelen=totlen*1.0/numseqs
    totGC=numGC*100.0/totlen
    wrng=0 #Turn off the warnings for ignored sequences, these have to be shown at this point.
else:
    print "There is no sequence to analyze. Please, check your file and\nthe value of 'l' (length) given."
    exitval()
#=========================================================================#
if mainout:
    print '- '*29+'-'
    if (minlen>0 and maxlen==0):
        print '***Only analyzing the sequences of Length>='+str(minlen)
    elif (minlen==0 and maxlen>0):
        print '***Only analyzing the sequences of Length=<'+str(maxlen)
    elif (minlen>0 and maxlen>0):
        print '***Only analyzing the sequences of Length: '+str(minlen)+' =< Length =< '+str(maxlen)
    print 'Filename:\t\t'+filename
    print 'Type(s) of sequences:\t'+', '.join(types)
    print 'Number of sequence:\t'+str(numseqs)
    print 'Min length:\t\t'+str(min(contgslen))
    print 'Max length:\t\t'+str(max(contgslen))
    print 'Average length:\t\t'+"{0:.2f}".format(avelen)
    print 'N50 value:\t\t'+"{0:.2f}".format(N50(contgslen))
    print 'Number of residues:\t'+str(totlen)
    print 'Total %GC in file:\t'+"{0:.2f}".format(totGC)
    print '- '*29+'-'
    pass
else:
    print '- '*34+'-'
    print 'Program:\t\tmultifastats - '+vers+' Rev by D.R.'
    print 'Filename:\t\t'+filename
    if (minlen>0 and maxlen==0):
        print '***Only analyzing the sequences of Length>='+str(minlen)+'.'
    elif (minlen==0 and maxlen>0):
        print '***Only analyzing the sequences of Length=<'+str(maxlen)+'.'
    elif (minlen>0 and maxlen>0):
        print '***Only analyzing the sequences of Length: '+str(minlen)+' =< Length =< '+str(maxlen)+'.'
    print 'Type(s) of sequences:\t'+', '.join(types)
    print 'Number of sequences:\t'+str(numseqs)
    print 'Min length:\t\t'+str(min(contgslen))
    print 'Max length:\t\t'+str(max(contgslen))
    print 'Average length:\t\t'+"{0:.2f}".format(avelen)
    print 'N50 value:\t\t'+"{0:.2f}".format(N50(contgslen))
    print 'Number of residues:\t'+str(totlen)
    print 'Total %GC in file:\t'+"{0:.2f}".format(totGC)
    print '- '*34+'-'
#=========================================================================#
if testps or testsga or testsbg:
    import time
    from _socket import timeout
    currTime = time.strftime("%d")+time.strftime("%m")+time.strftime("%y")+time.strftime("%H")+time.strftime("%M")+time.strftime("%S")
    yesopt=['y', 'yes', 's', 'si']
    noopt=['n', 'no']
#=========================================================================#
if testsbg:
    if sbginp: #If we are in the user-interactive mode, asking for a valid option for the subgroup output to continue.
        while testsbg:
            sbgopt=raw_input("\nDo you want to generate a multifasta file with the sequences analyzed?\nWrite 'Y' to make it or 'N' to continue): ")
            if sbgopt.lower() in yesopt:
                subgroupseq(filename,int(minlen),int(maxlen),currTime,wrng)
                print '\nDONE: A .fasta file with the sequences analyzed has been wrote in\nyour current working directory:'
                print '(subgroup_'+filename.replace(".","")+currTime+'.fasta)'
                print '- '*29+'-'
                testsbg=0
            elif sbgopt.lower() in noopt:
                testsbg=0
            else:
                print "Bad option, try again."
    else: #If we are in the command-line mode and the parameter for subgroup option is given, do the analysis.
        subgroupseq(filename,int(minlen),int(maxlen),currTime,wrng)
        print 'A .fasta file with the sequences analyzed has been wrote in\nyour current working directory:'
        print '(subgroup_'+filename.replace(".","")+currTime+'.fasta)'
        print '- '*34+'-'
        timeout(1)
#=========================================================================#
if testps:
    if psinp: #If we are in the user-interactive mode, asking for a valid option for the pseudosequence ooutput to continue.
        while testps:
            psopt=raw_input("\nDo you want to generate a pseudosequence?\nWrite 'Y' to build it or 'N' to continue): ")
            if psopt.lower() in yesopt:
                pschk=1
                while pschk:
                    nrep=raw_input("Please, give the NUMBER of repeats of the ambiguous\nresidue to concatenate the sequences: ")
                    try:
                        if int(nrep)>=0:
                            pseudoseq(filename,int(nrep),int(minlen),int(maxlen),currTime,wrng) #Correct number of repetitions, continue.
                            print '\nDONE: A .fasta file with the pseudosequence has been wrote in\nyour current working directory:'
                            print '(pseudoseq_'+filename.replace(".","")+currTime+'.fasta)'
                            print '- '*29+'-'
                            pschk=0
                            testps=0
                        else:
                            print "Incorrect value provided for the NUMBER of repeats, try again."
                    except ValueError:
                        print "A positive value is required as a NUMBER of repeats. Try again."
            elif psopt.lower() in noopt:
                testps=0
            else:
                print "Bad option, try again."
    else: #If we are in the command-line mode and the parameter for pseudosequence option is given, do the analysis.
        pseudoseq(filename,int(nrep),int(minlen),int(maxlen),currTime,wrng)
        print 'A .fasta file with the pseudosequence has been wrote in\nyour current working directory:'
        print '(pseudoseq_'+filename.replace(".","")+currTime+'.fasta)'
        print '- '*34+'-'
        timeout(1)
#=========================================================================#
if testsga:
    if sgainp: #If we are in the user-interactive mode, asking for a valid option to make the single analysis to continue.
        while testsga:
            sgaopt=raw_input("\nDo you want statistics for each single sequence?\nWrite 'Y' to calculate or 'N' to finish): ")
            if sgaopt.lower() in yesopt:
                doanlys=singanlys(filename,int(minlen),int(maxlen),currTime,wrng)
                print '\nDONE: A .csv file has been wrote with the single sequence\nstats in your current working directory:'
                print '(sganl_'+filename.replace(".","")+'_'+currTime+'.csv)'
                print '\nThank you!'
                print '- '*29+'-'
                testsga=0
                exitval()
            elif sgaopt.lower() in noopt:
                print '\nThank you!'
                print '- '*29+'-'
                testsga=0
                exitval()
            else:
                print "Bad option, try again."
                pass
    else: #If we are in the command-line mode and the parameter for single analysis is given, do the analysis.
        doanlys=singanlys(filename,int(minlen),int(maxlen),currTime,wrng)
        print 'A .csv file has been wrote with the single sequence stats in\nyour current working directory:'
        print '(sganl_'+filename.replace(".","")+'_'+currTime+'.csv)'
        print '\nThank you!'
        print '- '*34+'-'
        timeout(1)
        exit()
else:
    print 'Thank you!'
    print '- '*34+'-'
    exitval()
