module load gcc/14.2.0
module load python/3.13.0-35hwwod
python3 -m venv ~/bioenv
source ~/bioenv/bin/activate
pip install biopython
pip install pandas

for file in ./gbff_files/*.gbff
do

    python3 ./t6ss_finder/t6ss_finder.py -i "$file"  -f gbk -o T6SS -t blastp
done 
