#!/usr/bin/python
# -*- coding: utf-8 -*-
import os, sys, platform
vers='Version 1.4.2, 24-Oct-2014.'
versnote='''- Version 1.4.2 (D.R. 24-Oct-2014):
  Adding autocomplete using tab in Windows and MAC OS.'''
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
* biopython - http://biopython.org/wiki/Download
* pyreadline (for Windows and Mac) - https://pypi.python.org/pypi/pyreadline/
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
To use in this mode, run the script with NO arguments. The program will ask
you the name of the file you want to analyze (autocomplete allowed in LINUX).
To use the program in this mode, just run the script with doble-click or call
it from command-line as follows: >>~$ python multifastats.py

2) 'Command-Line': File name, length cut-off and single analysis options have
to be given from the beggining.
To use the program in this mode, run the script from command-line as follows:
>>~$ python multifastats.py -f (inputfile)
-f (or --file): \tallows to give an input file name. In the above line,
\t\t\treplace (inputfile) with the name of your file
\t\t\t(including the extension), like: myseq.fasta

You will add some options after the (inputfile) name, like:

-l (or --length):\tAdd a cutf-off for the sequences length to be analyzed:
\t\t\tgreater or equal to this value.
\t\t\tOnly positive numbers allowed!
\t\t\tIf you give 0, all sequences will be analyzed.
-s (or --single):\tAdd the single sequence analysis.
\t\t\tIMPORTANT: This option has to be the last given.

As example, one command line with all options will be:
>>~$ python multifastats.py -f inputfile.fasta -l 500 -s

INFO OPTIONS:
- - - - - - -
You will give some arguments to get some program info:
-h (or --help)\t :\tGeneral description and usage of the program.
-v (or --version):\tVersion and revision of the program.
-i (or --info)\t :\tInformation about the parameters calculated.
-n (or --notes)\t :\tNotes about the current and previous versions program.

Optional:
- - - - -
Maybe you will have more fun if you copy the script to '/usr/local/bin' 
and give some permisions (chmod +xr). This will allow you to call it and use
directly from the terminal in any directory of your computer and allow to
autocomplete the name of any file (not only python scripts).
You will do that with the following commands in the linux terminal:
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
docntend='''- Version 1.4.1 (D.R. 23-Oct-2014):
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

***For the Version 1.5, I'm implementing the following options:
- Filter by sequence type.
- Upper length cut-off.
- Filter by any string in sequence name.
- A fasta output of the selected subset of sequences.
- The pseudomolecule output.
- Top Blast-Hit for each single sequence (given certain parameters).
- 'TAB' autocompletion in Windows.
...Available coming soon!!!
Requests for new options are welcome!
===============================================================================
'''
docinfo='''
===============================================================================
*** PARAMETERS INFORMATION ***
===============================================================================
Multifastats: Multi-Fasta Sequence Stats
----------------------------------------
For sets of sequences, the parameters calculated are:

- Num of Sequences: Count the total number of sequences in the multifasta file.
- Min, Max and Average length: Gives the minimum, maximum and average sequence
  length in the multifasta file, respectively.
- N50 value: This parameter is calculated using the Broad Institute definition:
  https://www.broad.harvard.edu/crd/wiki/index.php/N50
- Number of residues: Count the total number of nucleotides (or aminoacids) in
  all the sequences in the multifasta file.
- Total %GC in file: This parameter calculates the percentage of G, C or S (G
  or C) nucleotides (no case sensitive) in the complete set of sequences.
----------------------------------------
For each sequence, the parameters calculated are:

A .csv file is write in the current working directory:
sganl_(FILENAME)_(TIME).csv
This contains each single sequence stats:

- N: Number of the sequence in the set of sequences.
- ID: ID (header) of the sequence.
- Length: Number of nucleotides (or aminoacids) in the sequence.
- %GC: Percentage of G, C or S (G or C) nucleotides (no case sensitive) in the
  sequence.
- Mol Weight: Molecular weight of the sequence, calculated according the
  molecular_weight() function from Biopython SeqUtils:
  http://biopython.org/DIST/docs/api/Bio.SeqUtils-module.html#molecular_weight
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
maininp=0 #This variable is the state of the manual input requirement for the multifasta analysis (1=YES/0=NO)
mainout=0 #This variable is the state of the type of output for the multifasta analysis (1=YES/0=NO)
singinp=0 #This variable is the state of the manual input requirement for the single analysis (1=YES/0=NO)
testval=1 #This variable is the state of do or not do the single analysis (1=YES/0=NO)
lctinp=0 #This variable is the state of the manual input requirement for the length cut-off
lctof=0 #This variable is the length cut-off (0 means no cut-off)
if len(sys.argv)==1: #This is the case when the script NO receive arguments (***THE "USER-INTERACTIVE" SCENARIO***).
    print '- '*29+'-'
    print 'Multifastats: Multi-Fasta Sequence Stats'
    print '-'*41
    print 'Show basic genomic statistics for a given set of sequences\nand optionally for each single sequence, too.'
    print 'Made by D.R. LBBM, UPCH. '+vers
    mainout=1
    singinp=1
    maininp=1
    lctinp=1
    pass
elif len(sys.argv)==2: #This is the case when the script receive arguments asking for information.
    if sys.argv[1].lower() in ('h','-h','help','-help','--help'):
        print dochpini
        print vers+dochpmid
        print dochpend
        exitval()
    elif sys.argv[1].lower() in ('v','-v','version','-version','--version'):
        print docvrini
        print vers+docvrend
        exitval()
    elif sys.argv[1].lower() in ('i','-i','info','-info','--info'):
        print docinfo
        exitval()
    elif sys.argv[1].lower() in ('n','-n','notes','-notes','--notes'):
        print docntini
        print versnote
        print docntend
        exitval()
    else: #We will add more options for info in the future versions.
        print "Incorrect way to give the arguments. See 'help' with '-h' or '--help' option."
        exitval()
elif (len(sys.argv) in range(3,7)) and (('-f' or '--file') in sys.argv[1].lower()): #This is the case when the script receive arguments (***COMMAND-LINE SCENARIO***) and the file option correctly.
    currDir=os.listdir(os.getcwd())
    if (sys.argv[2] in currDir): #If we receive the file name correctly.
        filename=sys.argv[2] #The file name is correct, so we will use this to continue.
        if len(sys.argv)==3:
            testval=0 #This means that we will NOT do the single analysis
            pass
        elif (len(sys.argv)==4) and (('-s' or '--single') in sys.argv[3].lower()): #If receive the single sequence analysis option correctly.
            pass #This means that we will do the single analysis without asking Yes or No.
        elif (len(sys.argv)==5 or len(sys.argv)==6) and (('-l' or '--length') in sys.argv[3].lower()): #If receive the length option correctly.
            if len(sys.argv)==5:
                try:
                    if float(sys.argv[4])>=0:
                        lctof=float(sys.argv[4]) #This means that we will use the length cut-off value given.
                        testval=0 #This means that we will NOT do the single analysis
                        pass
                    else:
                        print "Incorrect cut-off value for length. See 'help' with '-h' or '--help' option."
                        exitval()                    
                except ValueError:
                    print "A number is required as cut-off value for length. See 'help' with '-h' or '--help' option."
                    exitval()
            elif (len(sys.argv)==6) and (('-s' or '--single') in sys.argv[5]): #If receive the single sequence analysis option correctly.
                try:
                    if (float(sys.argv[4])>=0):
                        lctof=float(sys.argv[4]) #This means that we will use the length cut-off value given.
                        pass #This means that we will do the single analysis without asking Yes or No.
                    else:
                        print "Incorrect cut-off value for length. See 'help' with '-h' or '--help' option."
                        exitval()
                except ValueError:
                    print "A number is required as cut-off value for length. See 'help' with '-h' or '--help' option."
                    exitval()
            else:
                print "Incorrect way to give the arguments. See 'help' with '-h' or '--help' option."
                exitval()
        else:
            print "Incorrect way to give the arguments. See 'help' with '-h' or '--help' option."
            exitval()
    else: #If receive the file name incorrectly.
        print 'Incorrect file name'
        exitval()
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
        filename=raw_input("File name or Option: ")
        options=('h','-h','help','-help','--help','v','-v','version','-version','--version','i','-i','info','-info','--info','n','-n','notes','-notes','--notes')
        if filename.lower() in options:
            maininp=0 #This means that we get a correct file name or option, so we will use this to continue.
            lctinp=0
            if filename.lower() in ('h','-h','help','-help','--help'):
                print dochpini
                print vers+dochpmid
                print dochpend
            elif filename.lower() in ('v','-v','version','-version','--version'):
                print docvrini
                print vers+docvrend
            elif filename.lower() in ('i','-i','info','-info','--info'):
                print docinfo
            elif filename.lower() in ('n','-n','notes','-notes','--notes'):
                print docntini
                print versnote
                print docntend
            else: #We will add more options for info in the future versions.
                print "\nIncorrect way to give the arguments. See 'help' with '-h' or '--help' option."
            exitval()
        elif filename in os.listdir(os.getcwd()):
            maininp=0 #This means that we get a correct file name or option, so we will use this to continue.
            pass
        else:
            print "\nBad file name or option, try again or press 'Ctrl+C' to exit."
            pass
#=========================================================================#
if lctinp:
    clinp=0 #This variable allows to verify if a valid value for L is given.
    print '- '*29+'-'
    print "Do you want to analyze all the sequences? Or just a subset\nof the sequences greater or equal of a length 'L'?"
    while clinp==0:
        tochklc=raw_input("Please, give a positive value for L (or 0 to analyze all):")
        try:
            if (float(tochklc)>=0):
                clinp=1
                lctof=float(tochklc)
            else:
                print "Incorrect cut-off value for length."
                pass
        except ValueError:
            print "A number is required as cut-off value for length."
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
    seqDNA=Seq(oneseq,IUPACAmbiguousDNA())
    seqRNA=Seq(oneseq,IUPACAmbiguousRNA())
    seqProt=Seq(oneseq,ExtendedIUPACProtein())
    if Alphabet._verify_alphabet(seqDNA):
        seqtype='DNA'
    elif Alphabet._verify_alphabet(seqRNA):
        seqtype='RNA'
    else:
        if Alphabet._verify_alphabet(seqProt):
            seqtype='Protein'
        else:
            seqtype='invalidseq'
            pass
    return seqtype
#=========================================================================#
def N50(lenlist):
    # N50 calculation, based on the Broad Institute definition:
    # https://www.broad.harvard.edu/crd/wiki/index.php/N50
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
def singanlys(file,lncut): #While the function be nested in this script, just need one argument. But if we want to use independently, we need to change the input to: #def singanlys(file,listcntgslen,listgcs):
    import time
    from _socket import timeout
    from Bio.SeqUtils import GC, molecular_weight
    currTime = time.strftime("%d")+time.strftime("%m")+time.strftime("%y")+time.strftime("%H")+time.strftime("%M")+time.strftime("%S")
    outpbyseqcsv = open("sganl_"+file.replace(".","")+"_"+currTime+".csv","w") #Creating the file to write the output
    outpbyseqcsv.write('N,ID,Length,%GC,Mol Weight\n') #Headers line
    idsnocom=[] #To store the ID of each seq without commas.
    molwght=[] #To store the molecular weight of each sequence.
    multifs = open(file,"rU")
    for indsq in SeqIO.parse(multifs, "fasta"):
        typeindsq=seqcategory(str(indsq.seq))
        if lncut==0:
            if typeindsq=='invalidseq': #If invalid format, ignore this sequence.
                pass
            else:
                newid = indsq.id.replace(",", " ")
                idsnocom.append(newid)
                try:
                    molwght.append(molecular_weight(indsq.seq))
                except:
                    molwght.append(0)
                pass
        else:
            if len(indsq.seq)>=lncut:
                if typeindsq=='invalidseq': #If invalid format, ignore this sequence.
                    pass
                else:
                    newid = indsq.id.replace(",", " ")
                    idsnocom.append(newid)
                    try:
                        molwght.append(molecular_weight(indsq.seq))
                    except:
                        molwght.append(0)
                    pass
            else:
                pass
    count=0
    for i in range(len(contgslen)): #If we want to use the function independently, we need to change this line to:  #for i in range(len(listcntgslen)):
        count+=1
        newlinecsv=str(count)+','+str(idsnocom[i])+','+str(contgslen[i])+','+"{0:.2f}".format(gcs[i])+','+"{0:.2f}".format(molwght[i])+'\n' #If we want to use the function independently, we need to change this line to: #newlinecsv=str(count)+','+str(idsnocom[i])+','+str(listcntgslen[i])+','+"{0:.2f}".format(listgcs[i])+','+"{0:.2f}".format(molwght[i])+'\n'
        #print newlinecsv #I left this line because maybe someone wants to see the output in the python interpreter.
        outpbyseqcsv.write(newlinecsv)
        multifs.close()
        return currTime
#=========================================================================#
sqfs = open(filename,"rU")
numGC=0
currGC=0
contgslen=[] #To store contigs lengths.
gcs=[] #To store %GC values.
typesofseqs=[] #To store the type of sequences in the file.
for contg in SeqIO.parse(sqfs, "fasta"):  #If is in fasta format, calculate some parameters for sequences with valid format.
    typeseq=seqcategory(str(contg.seq))
    if lctof==0: #If we have to analyze all the sequences, use all the valid sequences.
        if typeseq in ('DNA', 'RNA', 'Protein'): #If valid format, use this sequence and calculate parameters.
            typesofseqs.append(typeseq)
            currGC=sum(contg.seq.count(x) for x in ['G', 'C', 'g', 'c', 'S', 's'])
            numGC+=currGC
            currlen=len(contg.seq)
            contgslen.append(currlen)
            gcs.append(currGC*100.0/currlen)
            pass
        elif typeseq=='invalidseq': #If invalid format, ignore this sequence.
            print 'WARNING!: Sequence '+contg.id+' has invalid format and not taken in account.'
            pass
        else:
            pass
    else: #If we have a sequence length as cut-off, use all the valid sequences in the subset.
        if len(contg.seq)>=lctof: #If the sequence length pass the cut-off, use this sequence (if valid).
            if typeseq in ('DNA', 'RNA', 'Protein'): #If valid format, use this sequence and calculate parameters.
                typesofseqs.append(typeseq)
                currGC=sum(contg.seq.count(x) for x in ['G', 'C', 'g', 'c', 'S', 's'])
                numGC+=currGC
                currlen=len(contg.seq)
                contgslen.append(currlen)
                gcs.append(currGC*100.0/currlen)
                pass
            elif typeseq=='invalidseq': #If invalid format, ignore this sequence and continue.
                print 'WARNING!: Sequence '+contg.id+' has invalid format and not taken in account.'
                pass
            else:
                pass
        else: #If the sequence length does not pass the cut-off, ignore this sequence and continue.
            pass
#=========================================================================#
if len(contgslen)>0:
    types=sorted([str(x) for x in set(typesofseqs)])
    sqfs.close()
    numseqs=len(contgslen)
    totlen=sum(contgslen)
    avelen=totlen*1.0/numseqs
    totGC=numGC*100.0/totlen
else:
    print "There is no sequence to analyze. Check your file and the value of 'l' (length) given"
    exit()
#=========================================================================#
if mainout:
    print '- '*29+'-'
    if lctof>0:
        print '***Only analyzing the sequences of Length>='+str(lctof)
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
    print 'Program:\t\tmultifastats.py. '+vers+' Rev by D.R.'
    print 'Filename:\t\t'+filename
    if lctof>0:
        print '***Only analyzing the sequences of Length>='+str(lctof)
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
if testval:
    if singinp: #If we are in the user-interactive mode for single outputs, asking for a valid option to continue:
        yesopt=['y', 'yes', 's', 'si']
        noopt=['n', 'no']
        while testval:
            singopt=raw_input("Do you want statistics for each single sequence?\nWrite 'Y' to calculate or 'N' to finish): ")
            if singopt.lower() in yesopt: 
                doanlys=singanlys(filename,lctof) #While the function be nested in this script, just need one argument. If we want to use the function independently, we need to replace this line by: #doanlys=singanlys(filename,contgslen,gcs)
                print '\nDONE: A .csv file has been wrote with the single sequence\nstats in your current working directory:'
                print '(sganl_'+filename.replace(".","")+'_'+doanlys+'.csv)'
                print '\nThank you!'
                print '- '*29+'-'
                testval=0
                exitval()
            elif singopt.lower() in noopt:
                print '\nThank you!'
                print '- '*29+'-'
                testval=0
                exitval()
            else:
                print "Bad option, try again."
                pass
    else: #If we are in the command-line mode for single outputs, do the analysis:
        doanlys=singanlys(filename,lctof) #While the function be nested in this script, just need one argument. If we want to use the function independently, we need to replace this line by: #doanlys=singanlys(filename,contgslen,gcs)
        testval=0
        print 'A .csv file has been wrote with the single sequence stats in\nyour current working directory:'
        print '(sganl_'+filename.replace(".","")+'_'+doanlys+'.csv)'
        print '\nThank you!'
        print '- '*34+'-'
        exit()
else:
    print 'Thank you!'
    print '- '*34+'-'
    exitval()
