from Bio import Entrez, SeqIO
import pandas as pd

Entrez.email = "clodagh.mcguire@nhs.net"
chr_accession = pd.read_csv('chr_accession_no.csv')
peaks = pd.read_csv('assigned_peaks/RPgenePeaks_all_adult_peaks_hg19.csv')
peaks = pd.merge(peaks, chr_accession, on="chr")
peaks["sequence"] = ""
peaks_list = []

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
    header = "\n>" + str(peaks.loc[i, "chr"]) + ":" + str(peaks.loc[i, "start"]) + "-" + str(peaks.loc[i, "end"])
    peaks_list.append(header)
    peaks_list.append(record.seq)

peaks.to_csv('peak_sequences.txt', sep='\t')
with open('peaks.fasta', 'w') as f:
    for item in peaks_list:
        f.write("%s\n" % item)