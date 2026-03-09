# Prediction of T6SS gene cluster in bacterial genomes

This program is the local version of T6SS prediction tool in SecReT6 v3 database (https://bioinfo-mml.sjtu.edu.cn/SecReT6/index.php)

## 1. Installation

Scripts for data preparation and T6SS prdiction is written in Python. Please ensure that Python >= 3.7.0 is installed in your environment.

Run the following command in shell to install the dependent Python packages listed in `requirements.txt`
`pip install -r requirements.txt`

The homology search tools `BLAST+ 2.11.0` and `HMMER 3.3.2` are already put in the `software` folder, so there is no need to install these tools by the user in general.

## 2. Usage

T6SS prediction tool detects T6SS gene clusters in an annotated or unannotated genome sequence by BLASTp 2.11.0+ or hmmsearch (HMMER v3.3.2) against experimentally validated T6SSs recorded in SecReT6.

(1) If your genome sequence is unannotated, you can input FASTA format nucleotide file and run the prediction. In this case, gene-finding tool Prodigal v2.6.3 is first utilized to annotate the sequence before detection of T6SS gene clusters using BLASTp 2.11.0+ or hmmsearch (HMMER v3.3.2).

(2) If your genome sequence is already annotated, you can input GenBank file containing both annotations and genome sequence. In this case, DNA sequences of annotated genes would be extracted and translated as queries for BLASTp 2.11.0+ or hmmsearch (HMMER v3.3.2).

After the prediction of the homologues of T6SS component proteins (T6CPs), this tool will identify T6SS regions by colocalization of T6CPs in the user-specified distance.

Command:
```
python your_install_path/t6ss_finder.py [-h] -i INPUT -f gbk/fasta [-o OUTDIR] [-t blastp/hmmsearch] [-e EVALUE] [-y IDENTITY] [-d DISTANCE] [-n NUM_PROTEINS]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        (required) Input query file containing genome sequence
  -f gbk/fasta, --format gbk/fasta
                        (required) The format of input file (gbk or fasta)
  -o OUTDIR, --outdir OUTDIR
                        The output directory to store the results (default: current directory)
  -t blastp/hmmsearch, --tool blastp/hmmsearch
                        The tool used for prediction (blastp or hmmsearch). Default is blastp
  -e EVALUE, --evalue EVALUE
                        E-value threshold (default: 0.001)
  -y IDENTITY, --identity IDENTITY
                        Percent identity for BLASTp search against T6SS component proteins (default: 30)
  -d DISTANCE, --distance DISTANCE
                        Maximum interval (bp) between two co-localized T6SS components (default: 15000)
  -n NUM_PROTEINS, --num_proteins NUM_PROTEINS
                        Minimum number of T6SS component proteins for each T6SS region (default: 3)
```

## 3. Output
The file `t6ss_prediction_result.txt` in the output directory contains the predicted T6SS regions and all BLASTp/hmmsearch hits.