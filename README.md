Counsyl Genome Search
===

This script searches a genome file to find potential instances of the protien sequence "countsyl" and then outputs the
dictionary to output.txt. This file could be then reparsed for further analysis.

### To run:

You must have biopython and genome files. For this example I used the hs_ref_GRCh37.p10_chr*.fa files located at

ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/Assembled_chromosomes/seq

On Windows/Mac/Linux:

  $ python counsyl.py


some simple output will show on the screen and the final dictionary will appear in output.txt

###Results

I only found 2 confirmed instances at chromX: [34354639] and chrom7: [41199145]. I was excited to find one on the
X chromosome since the whole nature of the company is to look at male and female genomes and predict offspring
genetic health. Very exciting! 

There's a potential that there are more instances since there a lot of N bases in this sequence. For simplisity
any codon with an N base is counted as a false hit and the search string combo is broken (if there was one).

Fixing the false hit for codons with N bases would be the first improvement to attend to if this were to be used
for mapping potential probe sites. Since I used a dictionary (hash table) to search for bases I would either have
to add all permutations of codons with N in them and allow for a list of amio acid characters to be associated to the
key OR use a regEx. I personally would think that the dictionary hash lookup would be faster at run time vs a regEx.
Then it would be a matter of seeing if any of the AA in the dictionary list matched the upcoming character in the 
search sequence before preceding to the next codon. The regEx would also need to be generated for every new search
string and could become extremely cumbersome and difficult to maintain.

This was a simple project but fun to do. I hope that you enjoy it and ask me questions on the code quality as well as
some of the choices I made during development.


	
###Enjoy!
