from Bio import Entrez, SeqIO
import pandas as pd

Entrez.email = "clodagh.mcguire@nhs.net"
chr_accession = pd.read_csv('chr_accession_no.csv')
peaks = pd.read_csv('RPgenePeaks_all_adult_peaks_hg19.csv')
peaks = pd.merge(peaks, chr_accession, on="chr")
peaks["sequence"] = ""

for i in range(len(peaks)):

    handle = Entrez.efetch(db="nucleotide",
                       id=peaks.loc[i, "accession"],
                       rettype="fasta",
                       strand=peaks.loc[i, "strand"],
                       seq_start=peaks.loc[i, "start"],
                       seq_stop=peaks.loc[i, "end"])
    record = SeqIO.read(handle, "fasta")
    handle.close()
    peaks.loc[i, "sequence"] = record.seq

peaks.to_csv('peak_sequences.txt', sep='\t')
