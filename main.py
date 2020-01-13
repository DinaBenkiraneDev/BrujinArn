import gzip
from graph import DeBrujinGraph
import random
import string

# Classe main
#Auteur: Julien Kiang, Dina Benkirane
#Question 3 et 4 incomplete
with gzip.open('GCF_000002985.6_WBcel235_rna.fna.gz', 'rt') as f:
    #fonction pour read le fichier fasta
    def read_fasta(path):
        with gzip.open(path, 'rt') as f:
            accession, description, seq = None, None, None
            for line in f:
                if line[0] == '>':
                    # yield current record
                    if accession is not None:
                        yield accession, description, seq

                    # start a new record
                    accession, description = line[1:].rstrip().split(maxsplit=1)
                    seq = ''
                else:
                    seq += line.rstrip()

next(read_fasta('GCF_000002985.6_WBcel235_rna.fna.gz'))

#Creation fichier contigs
contigs=open("contigs.fa","w+")
#Ecrire dans le fichier contigs
def write_contig(str):
    # Identifiant aléatoire
    random.seed(123)  # recommendé pour des résultats reproduisibles!
    contigs.write('>'.join(random.choices(string.ascii_uppercase, k=10)), str)
    contigs.close()

write_contig(read_fasta("GCF_000002985.6_WBcel235_rna.fna.gz"))

with gzip.open('reads.fastq.gz', 'rt') as f:
    #Fonction pour read les fichier fastq
    def read_fastq(path):
        with gzip.open(path, 'rt') as f:
            for line in f:
                sequence = f.readline().rstrip()
                #Question 2a
                question2=DeBrujinGraph(sequence)
                if sequence[0] is not '~':
                    yield sequence

next(read_fastq('reads.fastq.gz'))

k=21

def kmers():
    for _ in read_fastq('reads.fastq.gz'):
        kmers = [_[i:i + k] for i in range(len(_) - k + 1)]
        yield (kmers)

#Creation fichier occurrence.bed
occurrence=open("occurrences.bed","w+")
def read_bed(path):
    with open(path) as f:
        ref, start, end, name = f.readline().rstrip().split('\t')
        yield ref, int(start) - 1, int(end), name

def write_bed():
    pass
#next(read_bed('example.bed'))
#next(read_bed('occurences.bed'))

def separation():
    for i in kmers():
        for j in i:
            yield j


c = DeBrujinGraph(separation(), 21)
print(c.noPredecessors())


