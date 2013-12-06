# Rapsearch To Xml

___________________________________

This program takes the Rapseach output and transforms it in a xml file similar to blast xml format.
This is useful when we want to submit a Rapsearch output file to Blast2Go to perform the annotation step.

____________________________________

##Usage:

In the command line:

<pre><code>python RapsearchToXml.py infile.aln outfile.xml "example@mail.com"
</code></pre>

Notes:

+ This program uses fake values for the fields wich dont appear in the output of Rapsearch.
+ This fields will not interfere with the final result of the Blast2Go annotation.
+ The program uses the Entrez module from [biopython](https://github.com/biopython/biopython) adapted to [NCBI_Mass_Downloader](https://github.com/StuntsPT/NCBI_Mass_Downloader). A big Kudos to the authors;

________________________________________

##Dependencies:

+ Python3

________________________________________

##License:

GPLv2

___________________________________

##Known limitations:

For now:

+ It only accepts the rapsearchoutput.aln file
+ It only accepts one hit per query.

___________________________________

##Found a bug?

Or maybe just wanto to drop some feedback? Just open an issue on github!
   
