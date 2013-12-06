#!/usr/bin/python3

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.


#Usage: RapsearchToXml.py infile.aln outfile.xml "example@mail.com"

import xml.etree.cElementTree as ET
import xml.dom.minidom
import re
import Entrez
from sys import argv


def tabParse(tabfile):
	""" Not working perfectly yet
	It takes the output tabular file  from rapsearch(output.m8)
	and makes it readable for the xml write function
	It takes into account several hits per query
	"""
	rawscafs=open(tabfile)
	scafs = {}
	for line in rawscafs:
		if line.startswith("#") == False:
			splitvalues=line.split()
			if splitvalues[0] in scafs:
				scafs[splitvalues[0]] += [splitvalues[1:]]
			else:
				scafs[splitvalues[0]] = [splitvalues[1:]]
	rawscafs.close()
	return(scafs)



def alnParse(alnfile):
	"""
	Takes the output list file from rapsearch (outupt.aln)
	and makes it readable for the xml write function
	It returns a dictionary with query and hits information
	"""
	rawscafs=open(alnfile)
	scafs = {}
	for line in rawscafs:
		#Take of the No Hit lines. They wont appear in the final output like blast
		if line.strip().endswith("NO HIT"):
			continue
		elif line.strip().startswith(">"):
			line = re.sub("\S*=","", line)
			idline=line.replace("vs","").split()
			idline[0]=idline[0][1:]
			idline[4]=str(round(float(idline[4][:-1])))
		elif line.strip().startswith("Query"):
			idline += line.split()[1:]
		elif line.startswith("    "):
			idline.append(line.strip())
		elif line.strip().startswith("Sbjct"):
			idline += line.split()[1:]
			if idline[0] in scafs:
				scafs[idline[0]] += [idline[1:]]
			else:
				scafs[idline[0]] = [idline[1:]]
	rawscafs.close()
	return (scafs)
	
	
	
def hsp_checker(scafs, keyvalue):
	"""Not working yet. complementary to the tabParse Function"""
	acessionNumbers=[]
	for value in scafs[keyvalue]:
		acessionNumbers.append(value[0])
	if (len(acessionNumbers)==len(set(acessionNumbers))):
		return True
	else:
		return False

def get_description(mail,ID):
	"""Used by xml write for searching the definition of the hit using the accession number from NCBI
	which appears in the rapsearch output"""
	Entrez.email = mail
	handle = Entrez.efetch(db="protein", id=ID, rettype="gb", retmode="text")
	entry=(handle.read().strip())
	complete=entry.split("\n")
	definition=complete[1][12:]
	handle.close()
	return (definition)
	


def xmlwrite(scafs, outfile):
	"""Parse the information of the rapsearch output to a xml file
	with the same configuration of blast xml output
	This parser uses several default fields that don't appear in rapsearch output.
	Still not working with several hits and tabular file from rapsearch"""
	#makes the xml tree
	#Header of the xmlfile
	root= ET.Element("BlastOutput")
	BOprogram=ET.SubElement(root, "BlastOutput_program")
	BOprogram.text= "blastx"
	BOversion=ET.SubElement(root, "BlastOutput_version")
	BOversion.text="BLASTX 2.2.25+"
	BOreference=ET.SubElement(root, "BlastOutput_reference")
	BOreference.text="Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&amp;auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), &quot;Gapped BLAST and PSI-BLAST: a new generation of protein database search programs&quot;, Nucleic Acids Res. 25:3389-3402."
	BOdb=ET.SubElement(root,"BlastOutput_db" )
	BOdb.text="/home/databases/nr/nr"
	BOqueryID=ET.SubElement(root, "BlastOutput_query-ID")
	BOqueryID.text="Query_1"
	BOquerydef=ET.SubElement(root, "BlastOutput_query-def")
	BOquerydef.text="scaffold" #Default field
	BOquerylen=ET.SubElement(root, "BlastOutput_query-len")
	BOquerylen.text="500" #Default field
	BOparam=ET.SubElement(root, "BlastOutput_param")
	PA=ET.SubElement(BOparam, "Parameters")
	PAmatrix=ET.SubElement(PA, "Parameters_matrix")
	PAmatrix.text="BLOSUM62"
	PAexpect=ET.SubElement(PA, "Parameters_expect")
	PAmatrix.text="0.001"
	PAgapopen=ET.SubElement(PA, "Parameters_gap-open")
	PAgapopen.text="11"
	PAextend=ET.SubElement(PA, "Parameters_gap-extend")
	PAextend.text="1"
	PAfilter=ET.SubElement(PA, "Parameters_filter")
	PAfilter.text="L"
	BOiterations=ET.SubElement(root, "BlastOutput_iterations")
	x=1
	for key in scafs:
		IT=ET.SubElement(BOiterations, "Iteration")
		ITiter=ET.SubElement(IT, "Iteration_iter-num" )
		ITiter.text=str(x)
		ITqueryID=ET.SubElement(IT, "Iteration_query-ID")
		ITqueryID.text="Query_"+str(x)
		ITquerydef=ET.SubElement(IT, "Iteration_query-def")
		ITquerydef.text=key
		ITquerylen=ET.SubElement(IT, "Iteration_query-len")
		size=key.split("e")
		ITquerylen.text=size[1]
		IThits=ET.SubElement(IT, "Iteration_hits")
		x+=1
		y=1
		for value in scafs[key]:
			Hit=ET.SubElement(IThits, "Hit")
			hitnum=ET.SubElement(Hit, "Hit_num")
			hitnum.text=str(y)
			hitid=ET.SubElement(Hit, "Hit_id")
			hitid.text=value[0]
			hitdef=ET.SubElement(Hit,"Iteration_query-def")
			acession=value[0].split("|")
			hitdef.text=get_description(argv[3],acession[3])
			hitacess=ET.SubElement(Hit, "Hit_accession")
			hitacess.text=acession[3]
			hitlen=ET.SubElement(Hit, "Hit_len")
			hitlen.text="1234" #Default field
			hithsp=ET.SubElement(Hit, "Hit_hsps")
			Hsp=ET.SubElement(hithsp, "Hsp")
			hspnum=ET.SubElement(Hsp, "Hsp_num")
			hspnum.text="1"
			hspbitscore=ET.SubElement(Hsp, "Hsp_bit-score")
			hspbitscore.text=value[1]
			hspscore=ET.SubElement(Hsp, "Hsp_score")
			hspscore.text="1234" #Default field
			hspevalue=ET.SubElement(Hsp, "Hsp_evalue")
			hspevalue.text=str(10**(round(float(value[2]))))
			hspqstart=ET.SubElement(Hsp, "Hsp_query-from")
			hspqstart.text=value[8]
			hspqend=ET.SubElement(Hsp, "Hsp_hit-to")
			hspqend.text=value[10]
			hspqframe=ET.SubElement(Hsp, "Hsp_query-frame")
			hspqframe.text=value[7]#Default field
			hsphframe=ET.SubElement(Hsp, "Hsp_hit-frame")
			hsphframe.text=value[7]#Default field
			hspid=ET.SubElement(Hsp, "Hsp_identity")
			ident=round(int(value[3])*0.01*int(value[4]))
			hspid.text= str(ident)
			hsppositive=ET.SubElement(Hsp, "Hsp_positive")
			hsppositive.text=str(ident + value[11].count("+"))
			hspgap=ET.SubElement(Hsp, "Hsp_gaps")
			hspgap.text=value[6]
			hsplen=ET.SubElement(Hsp, "Hsp_align-len")
			hsplen.text=value[4]
			hspqseq=ET.SubElement(Hsp, "Hsp_qseq")
			hspqseq.text=value[9]
			hsphseq=ET.SubElement(Hsp, "Hsp_hseq")
			hsphseq.text=value[13]
			hspmidseq=ET.SubElement(Hsp, "Hsp_midline")
			hspmidseq.text=value[11]
			y+=1
		ITstat=ET.SubElement(IT, "Iteration_stat")
		Stat=ET.SubElement(ITstat, "Statistics")
		Statdbn=ET.SubElement(Stat, "Statistics_db-num")
		Statdbn.text="34201960"#Default field
		Statdbl=ET.SubElement(Stat, "Statistics_db-len")
		Statdbl.text="12001222213"#Default field
		Stathsplen=ET.SubElement(Stat, "Statistics_hsp-len")
		Stathsplen.text="127"#Default field
		Statefspace=ET.SubElement(Stat, "Statistics_eff-space")
		Statefspace.text="299643336524"#Default field
		Statkappa=ET.SubElement(Stat, "Statistics_kappa")
		Statkappa.text="0.041"#Default field
		Statlambda=ET.SubElement(Stat, "Statistics_lambda")
		Statlambda.text="0.267"#Default field
		Statent=ET.SubElement(Stat, "Statistics_entropy")
		Statent.text="0.14"#Default field
	#Write the tree to file	
	tree= ET.ElementTree(root)
	tree = ET.tostring(root, encoding="utf-8", method="xml")
	tree = tree.decode("utf-8")
	o = open(outfile,'w')
	salsa = xml.dom.minidom.parseString(tree)
	#Add header
	o.write("<?xml version=\"1.0\"?>\n<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"NCBI_BlastOutput.dtd\">")
	#Makes pretty identation on the xml text file
	o.write(salsa.toprettyxml(indent="  ")[22:])
	o.close()

xmlwrite(alnParse(argv[1]), argv[2])
