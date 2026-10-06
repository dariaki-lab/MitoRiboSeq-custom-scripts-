from Bio import SeqIO
from Bio.Seq import Seq

FASTA = "mtDNA.fa"
GFF3 = "mtDNA1_UTR.gff3"

OUT_CDS = "mtDNA_CDS.fa"
OUT_PROTEIN = "mtDNA_proteins.fa"

# Read genome
genome = SeqIO.read(FASTA, "fasta")
genome_seq = genome.seq

cds_features = []

with open(GFF3) as fh:
    for line in fh:
        if line.startswith("#"):
            continue

        fields = line.strip().split("\t")

        if len(fields) != 9:
            continue

        seqid, source, feature_type, start, end, score, strand, phase, attributes = fields

        if feature_type != "CDS":
            continue

        start = int(start)
        end = int(end)

        feature_id = None

        for item in attributes.split(";"):
            if item.startswith("Name="):
                feature_id = item.replace("Name=", "")
                break
            elif item.startswith("gene="):
                feature_id = item.replace("gene=", "")
                break
            elif item.startswith("ID="):
                feature_id = item.replace("ID=", "")
                break

        if feature_id is None:
            feature_id = f"CDS_{start}_{end}"

        cds_features.append(
            (feature_id, start, end, strand)
        )

with open(OUT_CDS, "w") as cds_out, open(OUT_PROTEIN, "w") as prot_out:

    for gene, start, end, strand in cds_features:

        seq = genome_seq[start-1:end]

        if strand == "-":
            seq = seq.reverse_complement()

        cds_out.write(f">{gene}\n{seq}\n")

        # Vertebrate mitochondrial code
        protein = Seq(str(seq)).translate(table=2)

        prot_out.write(f">{gene}\n{protein}\n")

print("Finished!")
print(f"CDS written to: {OUT_CDS}")
print(f"Proteins written to: {OUT_PROTEIN}")