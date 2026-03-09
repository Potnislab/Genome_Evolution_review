from Bio import SeqIO, SeqFeature
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i',  help = 'Input Genbank format file (include GI)')
parser.add_argument('--output', '-o',  help = 'Output faa format file')
args = parser.parse_args()

faa_out = []
records = list(SeqIO.parse(str(args.input), "genbank"))
for gbk in records:
    # process each gbk record here

#gbk = SeqIO.read(str(args.input), "genbank")
    for gene in gbk.features:
        if (gene.type == "CDS"):
            gi = str(gene.qualifiers.get('locus_tag'))
            if (len(gi.split("'")) > 1):
                gi = gi.split("'")[1]
            else:
                continue
	# Head line
            header = ('>%s' % gi)
            faa_out.append(header + '\n')
        # Extract amino acid sequences
            protein_seq = ''.join(re.findall(r'[A-Za-z]', str(gene.qualifiers.get('translation'))))
        # Insert "\n" per 70 letters into sequences 
            i = 0
            tmp = []
            while True:
                length_protein_seq = len(protein_seq)                
                if ((i + 70) < length_protein_seq):
                    tmp.append(protein_seq[i:(i + 70)])
                    i = i + 70
                else:
                    tmp.append(protein_seq[i:length_protein_seq])
                    break
            protein_seq = "\n".join(tmp) + "\n"
            faa_out.append(protein_seq) 

# Output amino acid sequences and DNA sequences
with open(str(args.output), "w") as f:
    f.writelines(faa_out)   
