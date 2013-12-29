# Rapsearch To Xml



This program takes the Rapseach output and transforms it in a xml file similar to blast xml format.
This is useful when we want to submit a Rapsearch output file to Blast2Go to perform the annotation step.



##Usage:


___To transform the output:___


In the command line:
<pre><code>python3 RapsearchToXml.py infile.aln outfile.xml "example@mail.com"
</code></pre>


*The description field is not working after run the blast2go4pipe.I am still figuring it out, but for now, i developed a work around:*

 
___To get the correct definition on the final annotation file:___


1. Run the RapsearchToXml.py script
2. Run the Blast2Go4Pipe with the output of th previous script
3. Get the output.dat and open it with Blast2Go with graphical interface. You will see the "seq definition" filled with "-NA-" instead of the actual description. 
4. Run Blast-Description-Annotator (BDA): Menu >> Tools >> Run BDA - It will replace the "-NA-" for the acession numbers.
5. Save the annotation to a file
6. **Then, in  the command line**:

<pre><code>python3 setDescription.py annoationfile.annot newannotationfile.annot "example@mail.com"
</code></pre>


Notes:

+ The RapsearchToXml.py uses fake values for the fields wich dont appear in the output of Rapsearch.
+ This fields will not interfere with the final result of the Blast2Go annotation.
+ It takes multiple hits per query
+ The program uses the Entrez module from [biopython](https://github.com/biopython/biopython) adapted to [NCBI_Mass_Downloader](https://github.com/StuntsPT/NCBI_Mass_Downloader). A big Kudos to the authors;

##Dependencies:

+ Python3



##License:

GPLv2



##Known limitations:

For now:
+ It only accepts the rapsearchoutput.aln file
+ Only works in the Blast2Go edition for comand line Blast2Go4pipe (still figuring it out)
+ Only works with a sequence name as "name|size1234", where 1234 is used for the size of the sequence




##Found a bug?

Or maybe just wanto to drop some feedback? Just open an issue on github!
   
