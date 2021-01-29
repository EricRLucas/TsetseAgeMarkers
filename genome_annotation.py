#!/usr/bin/python3
# genome_annotation.py


import os
import string
from re import *
from Bio.Blast import NCBIXML # this import is needed to interpret blast output xml files
from Bio.Blast.Applications import NcbiblastpCommandline # these imports are needed to create command lines that can be run from the system using os.sys. These are used to perform local blasts.
from Bio import SeqIO
import time # required to output the time at which the script was run
from sys import stdout # this import is need to flush python output to the stdout (instead of leaving it
# in the buffer
from sys import argv # this import is needed in order for the script to handle command line arguments
import socket # this import allows us to access the name of the computer that the script is being run on
from collections import OrderedDict as dict

# Create an exception class for the situatin in which the Blast database for the proteome is not found
class NoBlastDatabaseError(Exception):
	def __init__(self, path):
		Exception.__init__(self)
		print("Fail. Couldn't determine the BLAST database " + path + " for the reciprocal blast. Either the database has not been created, or its name does not match the fasta file.")


def write_to_file(string_to_be_written, filename, mode = 'w'):
	
	tofile = open(filename, mode)
	tofile.write(string_to_be_written)
	tofile.close()


def do_blast(tsetse_sequences_file, db):
	# Does a blast from the tsetse sequence to the alien proteome
	
	print('\trunning blast against database ' + db)
	stdout.flush()
	cline = NcbiblastpCommandline(query = tsetse_sequences_file, db=db,  evalue = 0.0001, outfmt=5, out=db + '_blast.xml')
	os.system(str(cline))
	blast_file = open(db + '_blast.xml', 'r')
	blastings = NCBIXML.parse(blast_file)
	return blastings


def find_protein_sequence(gene_name, the_fasta):
	# Given the results of a blast against a proteome, this function extracts the sequence that corresponds to the 
	# name of the top hit
	
	# extract the sequences from the proteome file
	pattern = sub(r'\\', r'\\\\', str(gene_name))
	pattern = sub('\|', '\|', str(pattern))  # The following lines add slashes to some special character so that the following regex parses properly
	pattern = sub('\]', '\]', str(pattern))
	pattern = sub('\[', '\[', str(pattern)) 
	pattern = sub('\(', '\(', str(pattern))
	pattern = sub('\)', '\)', str(pattern))
	pattern = sub('\+', '\+', str(pattern))
	pattern = sub('\*', '\*', str(pattern))
	pattern = sub('\^', '\^', str(pattern))
	pattern = sub('\.', '\.', str(pattern))
	pattern = sub('\?', '\?', str(pattern))
	pattern = sub('\$', '\$', str(pattern))
	sequence_search = findall('(?m)^.*' + pattern + '[^>]*', the_fasta)
	if len(sequence_search) != 1:
		raise Exception('Fail. There should be one and only one hit for a protein name')
	else:
		return sequence_search[0]


# The input Glossina morsitans proteome
morsitans_fastafile = 'Gmorsitans_genome/Glossina-morsitans-Yale_PEPTIDES_GmorY1.9.fa'
# The genes that we want to look for in the proteome
focal_gene_names = ['GMOY000749', 'GMOY001603', 'GMOY002920', 'GMOY003090', 'GMOY003371', 'GMOY003588', 'GMOY003952', 'GMOY005053', 'GMOY005321', 'GMOY009908', 'GMOY010976', 'GMOY011979']
# The filename to which to output the reduced fasta file containing only the genes of interest
reduced_fasta_output_filename = 'Gmorsitans_genome/Glossina-morsitans_genes_of_interest.fa'
# The database that we will blast against
database = 'Drosophila_proteome/Dmelano_reference_proteome'
# The final output filename
output_filename = 'tables/genome_annotation.csv'

# Run the script 

# Gene_identities.py is the file that will hold the dictionary associating contig names with candidate proteins
Gene_identities = dict()

# We bring the query proteome into memory
morsitans_proteome = SeqIO.parse(morsitans_fastafile, 'fasta')
# Pull out the genes of interest and write them to a new fasta file
focal_genes = [x for x in morsitans_proteome if sub('-P.', '', x.name) in focal_gene_names]
SeqIO.write(focal_genes, reduced_fasta_output_filename, 'fasta')

# Do the forward blast
allblasts = do_blast(reduced_fasta_output_filename, database)

# Iterate over all of the blasts (one for each sequence in the proteome)
for blast_record in allblasts:
	query_name = blast_record.query
	query_code = findall('[^ ]+', query_name)[0]
	forward_blast = blast_record.alignments
	# If at least one alien peptide was found, extract its sequence from the proteome file
	if forward_blast:
		Gene_identities[query_name] = forward_blast[0].title
	# If there were no hits for the forward blast, make a note of this in the log file
	else:
		Gene_identities[query_name] = 'XXX'

csvstring = 'Contig name\tAnnotation\n'
for prot, candidate in Gene_identities.items():
	csvstring += prot + '\t' + candidate + '\n'
	
print('Writing csv file')
stdout.flush()
write_to_file(csvstring, output_filename)


