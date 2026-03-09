#!/usr/bin/env python3.13
"""
T6SS prediction tool detects T6SS gene clusters in an annotated or unannotated genome sequence by BLASTp 2.11.0+ or hmmsearch (HMMER v3.3.2) against experimentally validated T6SSs recorded in SecReT6.
(1) If your genome sequence is unannotated, you can input FASTA format nucleotide file and run the prediction. In this case, gene-finding tool Prodigal v2.6.3 is first utilized to annotate the sequence before detection of T6SS gene clusters using BLASTp 2.11.0+ or hmmsearch (HMMER v3.3.2).
(2) If your genome sequence is already annotated, you can input GenBank file containing both annotations and genome sequence. In this case, DNA sequences of annotated genes would be extracted and translated as queries for BLASTp 2.11.0+ or hmmsearch (HMMER v3.3.2).
After the prediction of the homologues of T6SS component proteins (T6CPs), this tool will identify T6SS regions by colocalization of T6CPs in the user-specified distance.

The latest edition: 7/3/2022

Copyright 2022 Guan Jiahao <jhguan@sjtu.edu.cn>
 
Licensed under the GNU General Public License, Version 3
"""

from argparse import ArgumentParser
from Bio import SeqIO
import pandas as pd
import os
import re
import sys

def parse_args():
    parser = ArgumentParser(description = "Prediction of T6SS and its related proteins in bacterial genomes by NCBI BLASTp or HMMER")
    parser.add_argument("-i", "--input", type = str, required = True, help = "(required) Input query file containing genome sequence")
    parser.add_argument("-f", "--format", type = str, required = True, metavar = "gbk/fasta", help = "(required) The format of input file (gbk or fasta)")
    parser.add_argument("-o", "--outdir", type = str, required = False, default = os.getcwd(), help = 'The output directory to store the results (default: current directory)')
    parser.add_argument("-t", "--tool", type = str, required = False, default = "blastp", metavar = "blastp/hmmsearch", help = 'The tool used for prediction (blastp or hmmsearch). (default: blastp')
    parser.add_argument("-e", "--evalue", type = float, required = False, default = 0.001, help = 'E-value threshold (default: 0.001)')
    parser.add_argument("-y", "--identity", type = float, required = False, default = 30, help = 'Percent identity for BLASTp search against T6SS component proteins (default: 30)')
    #parser.add_argument("-a", "--havalue", type = float, required = False, default = 0.42, help = 'Ha-value threshold (0-1) for BLASTp (default: 0.42)')
    parser.add_argument("-d", "--distance", type = int, required = False, default = 15000, help = 'Maximum interval (bp) between two co-localized T6SS components (default: 15000)')
    parser.add_argument("-n", "--num_proteins", type = int, required = False, default = 3, help = 'Minimum number of T6SS component proteins for each T6SS region (default: 3)')
    parser.add_argument("-p", "--num_protein_types", type = int, required = False, default = 3, help = 'Minimum number of T6SS component protein types or domains for each T6SS region (default: 3)')
    parser.add_argument("-a", "--all", required = False, action = "store_true", default = False, help = 'Use all T6SS component proteins instead of experimentally validated T6SS component proteins as subject database (default: False)')

    return parser.parse_args()    

    
def main():
    args = parse_args()
    if args.tool != "blastp" and args.tool != "hmmsearch":
        sys.exit("The '-t/--tool' parameter must be either blastp or hmmsearch.")
    if not os.path.exists(args.outdir):
        try:
            os.mkdir(args.outdir)
        except Exception as ex:
            print(ex)

    run_t6ss_prediction(args.input, args.format, args.outdir, args.tool, args.evalue, args.identity, args.distance, args.num_proteins, args.num_protein_types, args.all)


def run_t6ss_prediction(input, format, outdir, tool, evalue, identity, distance, num_proteins, num_protein_types, all):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    blastp_path = os.path.join(script_dir, "software/ncbi-blast-2.11.0+/bin/blastp")
    hmmsearch_path = os.path.join(script_dir, "software/hmmer-3.3.2/src/hmmsearch")
    
    file = input.split('/').pop()
    file_name = os.path.splitext(file)[0]

 #file_name = file.split('.')[0]
    faa_file = os.path.join(outdir, file_name + ".faa")
    ptt_file = os.path.join(outdir, file_name + ".ptt")
    t6ss_result = os.path.join(outdir, file_name + ".txt")

    # Check format
    if format == "fasta":
        if is_fasta(input):
            # Run prodigal
            prodigal_path = os.path.join(script_dir, "software/prodigal_v2.6.3")
            prodigal_out = os.path.join(outdir, "prodigal.gff")
            prodigal_cmd = prodigal_path + " -i " + input + " -o " + prodigal_out + " -f gff"
            os.system(prodigal_cmd)
            # Convert .gff to .ptt and .faa
            print("Preparing for the input...")
            cmd_gff2ptt = "python " + os.path.join(script_dir, "software/gff2ptt.py") + " -i " + prodigal_out + " -o " + ptt_file
            cmd_ptt2faa = "python " + os.path.join(script_dir, "software/ptt2faa.py") + " -i " + ptt_file + " -f " + input + " -o " + faa_file
            os.system(cmd_gff2ptt)
            os.system(cmd_ptt2faa)
        else:
            sys.exit("Your input file is not in standard FASTA format.")
    elif format == "gbk":
        if is_gbk(input):
            # Convert .gbk/.gb to .ptt and .faa
            print("Preparing for the input...")
            cmd_gbk2faa = "python " + os.path.join(script_dir, "software/gbk2faa_locus_tag.py") + " -i " + input + " -o " + faa_file
            cmd_gbk2ptt = "python " + os.path.join(script_dir, "software/gbk2ptt.py") + " -i " + input + " -o " + ptt_file
            os.system(cmd_gbk2faa)
            os.system(cmd_gbk2ptt)
        else:
            sys.exit("Your input file is not in standard GenBank format.")
    else:
        sys.exit("The '-f/--format' parameter must be either gbk or fasta.")
    
    if not os.path.exists(faa_file) or not os.path.exists(ptt_file):
        sys.exit("The PTT file and/or FASTA file containing protein sequences were not existed.")

    if tool == "blastp":
        # Run BLASTp
        if all:
            db_path = os.path.join(script_dir, "db/t6cp_protein_all")
        else:
            db_path = os.path.join(script_dir, "db/t6cp_core_protein_exp")
        blast_out = os.path.join(outdir,file_name + "blast.out")
        blast_log = os.path.join(outdir, file_name + "blast.log")
        blast_cmd = blastp_path + " -query " + faa_file + " -db " + db_path + " -out " +  blast_out  + " -outfmt '6 delim= qaccver saccver pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads 4 -evalue " + str(evalue) + " -max_target_seqs 20 > " + blast_log + " 2>&1"
        print("Running BLASTp...")
        os.system(blast_cmd)
        if not os.path.exists(blast_out):
            sys.exit("The BLAST output file was not found.")
        print("Parsing BLASTp results for T6SS prediction...")
        content_list = parse_blastp_result(blast_out, ptt_file, identity, distance, num_proteins, num_protein_types, all)
        region_content, content = content_list[0], content_list[1]
        
    elif tool == "hmmsearch":
        # Run hmmsearch
        db_path = os.path.join(script_dir, "db/T6CP.hmm") 
        hmmsearch_out = os.path.join(outdir, "hmmsearch.out")
        hmmsearch_domtbl = os.path.join(outdir, "hmmsearch.domtbl")
        hmmsearch_cmd = hmmsearch_path + " --cpu 4 --noali -E " + str(evalue) + " --domE " + str(evalue) + " -o " +  hmmsearch_out  +  " --domtblout " + hmmsearch_domtbl + " " + db_path + " " + faa_file
        print("Running hmmsearch...")
        os.system(hmmsearch_cmd)
        if not os.path.exists(hmmsearch_domtbl):
            sys.exit("The hmmsearch output file was not found.")
        print("Parsing hmmsearch results for T6SS prediction...")
        content_list = parse_hmmsearch_result(hmmsearch_domtbl, ptt_file, evalue, distance, num_proteins, num_protein_types)
        region_content, content = content_list[0], content_list[1]        

    print("Writing T6SS prediction results to file...")
    with open(t6ss_result, 'w') as result:
            result.write(region_content)
            #result.write(content)
    content.to_csv(t6ss_result, index = False, sep = "\t", mode = 'a')
    print("Prediction is done!")

def is_fasta(file):
    with open(file, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file


def is_gbk(file):
    with open(file, "r") as handle:
        gbk = SeqIO.parse(handle, "genbank")
        return any(gbk)  # False when `gbk` is empty, i.e. wasn't a GenBank file


def parse_blastp_result(blast_out, ptt_file, identity, distance, num_proteins, num_protein_types, all):
    left_gi = {} ## locus_tag as key, left end as value
    left_coor = {} ## left end as value, right end as key
    gi_line = {} ## locus_tag as key, line as content;
    script_dir = os.path.dirname(os.path.realpath(__file__))
    if all:
        t6ss_protein_table = os.path.join(script_dir, "db/t6ss_protein_new_subset.csv")
    else:
        t6ss_protein_table = os.path.join(script_dir, "db/t6ss_protein_new.csv")
    df_t6ss_protein = pd.read_csv(t6ss_protein_table, header = None)
    df_t6ss_protein = df_t6ss_protein[[9, 3]]
    df_t6ss_protein.columns = ['T6CP_ID', 'Description']

    # Read PTT info
    with open(ptt_file, 'r') as ptt:        
        lines = ptt.readlines()
        seq_len = int(re.search("([0-9]+)\.\.([0-9]+)", lines[0]).group(2))
        
        seq_lines_started = False
        for line in lines:
            line = line.strip()
            if not line:
                continue  # skip empty lines
            if "Location" in line and "Strand" in line:
                seq_lines_started = True
                continue
            if not seq_lines_started:
                continue  # skip pre-header lines
            line_split = line.split("\t")
        #for line in lines:
         #   line = line.strip()
          #  line_split = line.split("\t")
            coordinates = line_split[0]
            locus_tag = line_split[5]
            coord_search = re.search(r"([0-9]+)\.\.([0-9]+)", coordinates)
            left, right = coord_search.group(1), coord_search.group(2)
            left_gi[locus_tag] = left
            left_coor[left] = right
            gi_line[locus_tag] = line


    # Parse BLASTp result
    hit_left = {} ## left end as key, locus_tag as value
    hit_left_coor = {} ## left end as key, right end as vlue
    gi_tmp = 0
    hit_line = {}
    #content = "###All BLASTp hits\nCoordinates\tStrand\tLength\tPID\tGene\tLocus tag\tCode\tCOG\tProduct\tBest hit\tIdentities (%)\tHa-value\tE-value\tBit score\n"
    content = pd.DataFrame(columns=['Coordinates', 'Strand', 'Length', 'PID', 'Gene', 'Locus_tag', 'Code', 'COG', 'Product', 'T6CP_ID', 'Identities (%)', 'E-value', 'Bit_score'])
    hit_num = 0

    with open(blast_out, 'r') as blast:
        lines = blast.readlines()
        for line in lines:
            line = line.strip()
            if not line:
               continue
            line_split = line.split("\t")

        # Skip lines that don't have enough columns
            if len(line_split) < 14:
                print(f"Skipping malformed line: {line}")
                continue

            identity_predict = float(line_split[2])
            if identity_predict < identity:
                continue

            q_gi = line_split[0]
            s_gi = line_split[1]
            
            #identity_predict = float(line_split[2])
            #q_len = int(line_split[3])
            evalue = line_split[12]
            bitscore = line_split[13]
            #match_length = int(line_split[5])
            #havalue_predict = identity_predict / 100 * match_length / q_len
            #havalue_predict = round(havalue_predict, 3)

            #if havalue_predict < havalue:
            #    continue

            if identity_predict < identity:
                continue

            q_gi = line_split[0]
            if q_gi == gi_tmp:
                continue
            gi_tmp = q_gi
            s_gi = line_split[1]

            #hit_line[q_gi] = gi_line[q_gi] + "\t" + s_gi + "\t" + str(identity_predict) + "\t" + str(havalue_predict) + "\t" + evalue + "\t"  + bitscore
            hit_line[q_gi] = gi_line[q_gi] + "\t" + s_gi + "\t" + str(identity_predict) + "\t" + evalue + "\t"  + bitscore
            #content += hit_line[q_gi]
            content.loc[content.shape[0]] = hit_line[q_gi].split("\t")
            hit_left[left_gi[q_gi]] = q_gi
            hit_left_coor[left_gi[q_gi]] = left_coor[left_gi[q_gi]]
            hit_num += 1

    content = content.merge(df_t6ss_protein, how = 'left', on = 'T6CP_ID')
    content.drop(columns = ['Gene', 'Code', 'COG'], inplace = True)
            
    # Output result
    if hit_num == 0:
        sys.exit("No T6SS was detected.")

    ## co-localization
    region = []
    region_gi = []
    i = 0
    region_num = 0

    hit_left_coor_keys = list(hit_left_coor.keys())
    hit_left_coor_keys = [int(i) for i in (hit_left_coor_keys)]
    hit_left_coor_keys.sort()

    for left in hit_left_coor_keys:
        if i == 0:
            region.append([left, int(hit_left_coor[str(left)])])
            region_gi.append(hit_left[str(left)])
            i += 1
        else:
            if left <= (int(region[region_num][1]) + int(distance)):
                region[region_num][1] = hit_left_coor[str(left)]
                region_gi[region_num] += "\t" + hit_left[str(left)]
                i += 1
            else:
                region_num += 1
                region.append([left, int(hit_left_coor[str(left)])])
                region_gi.append(hit_left[str(left)])
                i = 1

    num = 1
    region_content = ""

    # Deal with situation in which two parts are separated in two ends of a circular genome
    region, region_gi = check_separate_loci(region, region_gi, seq_len)
    region, region_gi = check_num_protein_types(region, region_gi, content, num_proteins, num_protein_types)
    
    for i in range(0, len(region)):
        region_content += "#Region" + str(num) + "\t" + str(region[i][0]) + ".." + str(region[i][1]) + "\t" + region_gi[i] + "\n"
        num += 1

    return region_content, content


def parse_hmmsearch_result(hmmsearch_out, ptt_file, evalue, distance, num_proteins, num_protein_types):
    left_gi = {} ## locus_tag as key, left end as value
    left_coor = {} ## left end as value, right end as key
    gi_line = {} ## locus_tag as key, line as content;

    # Read PTT info
    with open(ptt_file, 'r') as ptt:
        lines = ptt.readlines()
        seq_len = int(re.search("([0-9]+)\.\.([0-9]+)", lines[0]).group(2))
        del(lines[0:3]) # Skip the first 3 rows
        for line in lines:
            line = line.strip()
            line_split = line.split("\t")
            coordinates = line_split[0]
            locus_tag = line_split[5]
            coord_search = re.search(r"([0-9]+)\.\.([0-9]+)", coordinates)
            left, right = coord_search.group(1), coord_search.group(2)
            left_gi[locus_tag] = left
            left_coor[left] = right
            gi_line[locus_tag] = line

    # Parse hmmsearch result
    hit_left = {} ## left end as key, locus_tag as value
    hit_left_coor = {} ## left end as key, right end as vlue
    gi_tmp = 0
    hit_line = {}
    #content = "###All hmmsearch hits\nCoordinates\tStrand\tLength\tPID\tGene\tLocus tag\tCode\tCOG\tProduct\tHMM profile\tPfam accession\tE-value\tScore\t\n"
    content = pd.DataFrame(columns=['Coordinates', 'Strand', 'Length', 'PID', 'Gene', 'Locus_tag', 'Code', 'COG', 'Product', 'HMM_profile', 'Pfam_accession', 'E-value', 'Score'])
    hit_num = 0

    with open(hmmsearch_out, 'r') as hmmsearch:
        lines = hmmsearch.readlines()        
        for line in lines:
            if re.search('^#', line):
                continue
            line = line.strip()
            line_split = line.split()
            hmm_profile = line_split[3]
            accession = line_split[4]
            evalue_predict = line_split[12]
            score = line_split[13]

            if float(evalue_predict) > evalue:
                continue
            q_gi = line_split[0]
            if q_gi == gi_tmp:
                continue
            gi_tmp = q_gi

            hit_line[q_gi]= gi_line[q_gi] + "\t" + hmm_profile + "\t" + accession + "\t" + evalue_predict + "\t"  + score
            #content += hit_line[q_gi]
            content.loc[content.shape[0]] = hit_line[q_gi].split("\t")
            hit_left[left_gi[q_gi]] = q_gi
            hit_left_coor[left_gi[q_gi]] = left_coor[left_gi[q_gi]]
            hit_num += 1

    content.drop(columns = ['Gene', 'Code', 'COG'], inplace = True)
    content.sort_values(by = ['Locus_tag', 'Coordinates'], inplace = True)

    # Output result
    if hit_num == 0:
        sys.exit("No T6SS was detected.")

    ## co-localization
    region = []
    region_gi = []
    i = 0
    region_num = 0

    hit_left_coor_keys = list(hit_left_coor.keys())
    hit_left_coor_keys = [int(i) for i in (hit_left_coor_keys)]
    hit_left_coor_keys.sort()

    for left in hit_left_coor_keys:
        if i == 0:
            region.append([left, int(hit_left_coor[str(left)])])
            region_gi.append(hit_left[str(left)])
            i += 1
        else:
            if left <= (int(region[region_num][1]) + int(distance)):
                region[region_num][1] = hit_left_coor[str(left)]
                region_gi[region_num] += "\t" + hit_left[str(left)]
                i += 1
            else:
                region_num += 1
                region.append([left, int(hit_left_coor[str(left)])])
                region_gi.append(hit_left[str(left)])
                i = 1

    num = 1
    region_content = ""

    # Deal with situation in which two parts are separated in two ends of a circular genome
    region, region_gi = check_separate_loci(region, region_gi, seq_len)
    region, region_gi = check_num_protein_types(region, region_gi, content, num_proteins, num_protein_types)

    for i in range(0, len(region)):
        region_content += "#Region" + str(num) + "\t" + str(region[i][0]) + ".." + str(region[i][1]) + "\t" + region_gi[i] + "\n"
        num += 1

    return region_content, content


def check_separate_loci(region, region_gi, seq_len):
    separate_loci = {}
    for i in range(0, len(region)):
        region[i][0], region[i][1] = int(region[i][0]), int(region[i][1])
        if region[i][0] < 15000:
            separate_loci['left_start'] = region[i][0]
            separate_loci['left_end'] = region[i][1]
            separate_loci['left_region'] = i
        if region[i][1] > 15000 and (seq_len - region[i][1]) < 15000:
            separate_loci['right_start'] = region[i][0]
            separate_loci['right_end'] = region[i][1]
            separate_loci['right_region'] = i
    if len(separate_loci.keys()) == 6:
        if (seq_len - separate_loci['right_end'] + separate_loci['left_start']) < 15000:
            region[separate_loci['left_region']][0] = separate_loci['right_start']
            region[separate_loci['left_region']][1] = separate_loci['left_end']
            region_gi[separate_loci['left_region']] += '\t' + region_gi[separate_loci['right_region']]
            del(region[separate_loci['right_region']])
            del(region_gi[separate_loci['right_region']])
    return region, region_gi


def check_num_protein_types(region, region_gi, content, num_proteins, num_protein_types):
    content_process = content.copy(deep=True)
    content_process['start'] = content_process['Coordinates'].str.split('\.\.').str[0].astype(int)
    content_process['end'] = content_process['Coordinates'].str.split('\.\.').str[1].astype(int)
    
    region_keep = []
    region_gi_keep = []
    for i in range(0, len(region)):
        start_line = content_process[content_process.start == region[i][0]].index[0]
        end_line = content_process[content_process.end == region[i][1]].index[0]

        if region[i][0] < region[i][1]:
            content_tmp = content_process.loc[start_line:end_line]
        else:
            content_tmp = pd.concat([content_process.loc[0:end_line], content_process.loc[start_line:]])      
        num_t6cp = content_tmp.shape[0]

        if "Description" in content_process.columns:
            num_t6cp_type = content_tmp.drop_duplicates('Description').shape[0]
        else:
            num_t6cp_type = content_tmp.drop_duplicates('HMM_profile').shape[0]

        if num_t6cp >= num_proteins and num_t6cp_type >= num_protein_types:
            region_keep.append(region[i])
            region_gi_keep.append(region_gi[i])
        
    return region_keep, region_gi_keep
        

if __name__ == "__main__":
    main()
