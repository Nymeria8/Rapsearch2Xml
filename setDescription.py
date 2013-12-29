from sys import argv
import Entrez

def get_description(mail,ID):
	"""Used by make description for searching the definition of the hit using the accession number from NCBI
	which appears in the nnotation output
	Adapted from the RapsearchToXml.py file"""
	Entrez.email = mail
	handle = Entrez.efetch(db="protein", id=ID, rettype="gb", retmode="text")
	entry=(handle.read().strip())
	complete=entry.split("\n")
	definition=complete[1][12:] #get the definition camp
	definition2=definition.split("[")#removes the species informaton
	handle.close()
	return (definition2[0])


def make_description(anotation, anotationnew):
	""" makes another anottation file with the same information, but it substitute the acession number
	by the real description"""
	anot= open (anotation)
	anotnew=open(anotationnew, "w")
	for i in anot:
		seq= i.split()
		if (len(seq)==3):
			anotnew.write(seq[0]+"\t"+seq[1]+"\t")
			anotnew.write(get_description(argv[3], seq[2]) + "\n")
		else:
			anotnew.write(i)
	anot.close()
	anotnew.close()
			
make_description(argv[1],argv[2])
