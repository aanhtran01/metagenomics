MetaGenomics 
============

The objective if this Python project is to conduct metagenomic analysis on a given sample.  The project receives multiple genome sequences and reads from the sample as input. Its outcome involves determining the likely source genome for each read.  Each sample contains a specific quantity of genomes.

The program randomly samples 200 reads and aligns it to the all the genomes. Genomes that pass a 10% match threshold will be selected for full alignment and the other genomes will be discarded. 

The input of this project will be a zip file that contains all the genome files and one read file. 

It is important that the genome files contain the word "genome" in the file name and the read file contains the words "reads.fasta" in this exact format because this is how the program distinguishes and extracts the genome files and the read files from the zipfile. 

The output is a list of reads and the genome that it matches to in the following format:
>read_0	Genome_Number_45
>read_1	Genome_Number_65
>read_2	Genome_Number_65
>read_3	Genome_Number_50
>read_4	Genome_Number_32
...


Deliverables:
-------------

metagenomics.py -- code for metagenomic analysis

predictions.csv -- a list of reads and the genome that it matches to 

predictions.zip -- zipped csv of predictions.txt


Usage
-----
The program takes in one input, the reads fasta without the genome positions 

To run the program, navigate to the project directory and run:

> python3 metagenomics.py samples.zip

The program takes the following argument:

* `--samples.zip`: A zip file that contains all the genome files and one read file. It is important that the genome files contain the word "genome" in the file name and the read file contains the words "reads.fasta" in this exact format because this is how the program distinguishes and extracts the genome files and the read files from the zipfile. 

Examples
--------

Here is an example of how to run the program: 

>python3 metagenomics.py project4-sample.zip

Performance
-----------

Runtime is based on laptop computer with the following specifications:
* Processor: 1.1 GHz Quad-Core Intel Core i5
* RAM: 8GB
* Operating System: MacOS Catalina 

For alignment of ~1000 genomes of ~10,000 nucleotides and ~20,000 reads the runtime is: 

real	1m7.535s
user	1m5.817s
sys	0m0.459s

Future Improvements
-------------------
For a small sample of reads like this, it is relatively fast to align sequences. However, for larger amounts of data, more effective data structures would be neccessary for a  reasonable runtime. 