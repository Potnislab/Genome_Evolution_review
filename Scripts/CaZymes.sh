{\rtf1\ansi\ansicpg1252\cocoartf2822
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fnil\fcharset0 HelveticaNeue;}
{\colortbl;\red255\green255\blue255;\red25\green25\blue25;}
{\*\expandedcolortbl;;\cssrgb\c12941\c12941\c12941;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs29\fsmilli14667 \cf2 \expnd0\expndtw0\kerning0
#Code to run Dbcan to identify cazymes\
#!/bin/bash\
# load the module\
source /apps/profiles/modules_asax.sh.dyn\
module load run_dbcan\
\uc0\u8232 # Loop through all .faa files\u8232 for file in *.faa; do\u8232 tag=$\{file%.faa\}\u8232 run_dbcan "$file" protein --out_dir "$\{tag\}_output"\
\uc0\u8232 done\
\'a0\
#Once you get the output run the below code to extract cazymes hit with atleast 2 programs to make a new file\
\'a0\
#!/bin/sh\
source /apps/profiles/modules_asax.sh.dyn\
ls *.faa | sed 's/.faa//g' > file_list.txt\
files=`cat file_list.txt`\
for var in $files\
\'a0\
do\
\'a0\
echo $var\
awk '($5 >= 3) \{print $1 "\\t" $4\}' $\{var\}_output/overview.txt\
awk -F "_" '\{print $1 $2\}' $\{var\}_output/overview.txt\
awk -F "\\t" '\{count[$2]++\} END \{for (word in count) print word "\\t" count[word]\}' $\{var\}_output/overview.txt > output/$\{var\}_count.txt\
sed -i 's/^/'$\{var\}' /g' output/$\{var\}_count.txt\
sed -i '1i 'Genome' 'Cazymes_family' 'Count'' output/$\{var\}_count.txt\
sed -i 's/ /\\t/g' output/$\{var\}_count.txt\
\'a0\
done\
\'a0}