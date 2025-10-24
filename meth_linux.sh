#!/bin/bash

date

# mkdir {de_dup,methylation,cytosine,clean,repeat,imprinting,graph}

# fastqc -t 5 ./01.RawData/*/*.gz -o ./01.RawData
# multiqc ./01.RawData

fastqc -t 5 ./trim/*.gz -o ./trim
multiqc ./trim -o ./trim

cd ./clean

cat ../config_trim  |while read id;
do echo $id

arr=($id)
fq2=${arr[2]}
fq1=${arr[1]}
sample=${arr[0]}

	bismark -bowtie2 --sam -p 8 /mnt/h/wy/multi-omics/project/genome/mouse/mm9 -1 $fq1 -2 $fq2   >> bowtie2.log 2>&1

done



date

cd ../de_dup

ls ../clean/*.sam | while read id
do
	deduplicate_bismark --sam -p $id --output_dir .
done

date

cd ../methylation

ls ../de_dup/*.sam | while read id
do
        bismark_methylation_extractor --comprehensive --no_overlap -p --bedGraph --CX --cytosine_report --buffer_size 40% --genome_folder /mnt/h/wy/multi-omics/project/genome/mouse/mm9/Bisulfite_Genome/ ${id}
done

date

cd ../cytosine

ls ../methylation/*_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | while read id ;do (coverage2cytosine -CX --gc --genome_folder /mnt/h/wy/multi-omics/project/genome/mouse/mm9  -o $(basename $id "_1_bismark_bt2_pe.deduplicated.bismark.cov.gz") $id) ;done

echo " \n \n \n #  All work down !!! \n \n \n"

date


ls ./*_val.CX_report.txt | while read id

do
awk  '{print $1,$2,$3,$4,$5,$6}' $id | grep CG > $(basename $id "_val.CX_report.txt").CG_report.txt
awk 'BEGIN{ FS=" ";OFS="\t"} {print $1,$2,$2+1,$3,$4,$5,$4+$5,$6}' $(basename $id "_val.CX_report.txt").CG_report.txt > $(basename $id "_val.CX_report.txt").CG_report.bed
rm $(basename $id "_val.CX_report.txt").CG_report.txt
# echo -e " \n \n $(basename $id "_val.CX_report.txt").CG_report.txt was removed!!! \n \n"

# c
sed '/chrM/d' $(basename $id "_val.CX_report.txt").CG_report.bed | awk 'BEGIN{ FS=" ";OFS="\t"} {print $1,$2,$3,$5,$6,$7}' | awk '$6 >= 1 {print $0}' | awk -F "\t" '{OFS = "\t"}{ print $1,$2,$3,$4/$6}' > $(basename $id "_val.CX_report.txt").wig

# bedgraph
# sed '/chrM/d' $(basename $id "_val.CX_report.txt").CG_report.bed | awk 'BEGIN{ FS=" ";OFS="\t"} {print $1,$2,$3,$5,$6,$7}' | awk '$6 >= 1 {print $0}' | awk -F "\t" '{OFS = "\t"}{ print $1,$2,$3,$4/$6}' > $(basename $id "_val.CX_report.txt").bedGraph

# sort wig chr
sort -k1,1V -k2,2n -k3,3n $(basename $id "_val.CX_report.txt").wig | sed -n "/chr/p" | sed 's/chr//g' > $(basename $id "_val.CX_report.txt").sorted.wig

# wig to bw
wigToBigWig $(basename $id "_val.CX_report.txt").sorted.wig /mnt/g/wy/ChIP/fPSC/fPSC_try/align/mm9.chrom.sizes $(basename $id "_val.CX_report.txt").bw


done
