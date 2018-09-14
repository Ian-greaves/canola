#FILE TRANSFER

#sftp -P 1046 gre486@upload.data.csiro.au:/37734/34860v001/data/

#sftp -P 22 gre486@pearcey-login.hpc.csiro.au:/flush1/gre486/Canola_data

#sftp -P 22 gre486$osm-13-cdc.it.csiro.au:OSM/CBR/AF_HETEROSIS/work/2018_dataschool/data/

#put /OSM/CBR/AF_HETEROSIS/work/Anyu/C_21_*.fq.gz


#Pulls out first line of all files
#for file in C_4_1_1_1.fq.gz; do zcat $file | head n1 ; done

#But need to check if a file has multiple lanes or multiple machines
#for file in C_4_*.fq.gz; do zcat $file | cut 


#for file in C_4_1_1_1.fq.gz; do sample=${file%.*} ; zcat $file | grep -E '^@' | sed -E "s/ /:/g" | cut -d ':' -f 1,2,3,4,8,9,11 | sort -u | awk '{print $0 "$sample"}' > x.txt ; done 

#for file in C_4_1_1_1.fq.gz; do zcat $file | grep -E '^@' | sed -E "s/ /:/g" | cut -d ':' -f 1,2,3,4,8,9,11 | sort -u  > x.txt ; done

#for file in C_4_*.fq.gz; do zcat $file | head -n1 | grep -E '^@' | sed "s/ /:/g" | cut -d ':' -f 1,2,3,4,8,9,11 | sort -u > x.txt ; done

#Final script for use
for file in C_4_*.fq.gz; do
	zcat $file | head -n1 | grep -E '^@' | sed "s/ /:/g" | cut -d ':' -f 1,2,3,4,8,9,11 >> x.txt
	echo $file >> y.txt
done
paste x.txt y.txt > metadata.txt




