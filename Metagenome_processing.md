
# Metagenome sequence processing

## Raw sequence processing


```R
#Name raw sequences
cat *R1* > samplename_R1.fastq.gz
cat *R2* > samplename_R2.fastq.gz
```


```R
#Filter with Cucurbita pepo reference genome *
#Format reference genome
#!bin/bash
bowtie2-build GCA_002806865.2_ASM280686v2_genomic.fna calabacitagenome

# *The squash genome can be dowanload with the biproject number PRJNA421955:
# Montero-Pau J et al., "De novo assembly of the zucchini genome reveals a whole-genome duplication associated with the origin of the Cucurbita genus.", Plant Biotechnol J, 2017 Dec 4;16(6):1161-1171
# Alverson AJ et al., "Insights into the evolution of mitochondrial genome size from complete sequences of Citrullus lanatus and Cucurbita pepo (Cucurbitaceae).", Mol Biol Evol, 2010 Jan 29;27(6):1436-48
# The file name is GCA_002806865.2_ASM280686v2_genomic.fna, available at:
#ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/806/865/GCF_002806865.1_ASM280686v2
```


```R
#Align against reference genome
#!/bin/bash
#bash nombre_shipt.sh

SEQS=/sequence_location
BIN=/bowtie2_location

COUNT=0
for FAA in `ls *.fastq | perl -pe 's/\_.*//g' | sort | uniq`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.sh
echo "$BIN/bowtie2 -x calabacitagenome -1 $SEQS/$FAA"_R1.fastq" -2 $SEQS/$FAA"_R2.fastq" -S $FAA.sam"  -p 20 >>$*.$COUNT.sh
chmod +x *.sh; done

for N in `ls  bwt.* `; do echo "executable = $N"  >$N.condor; echo "output = $N.out" >>$N.condor; echo "error =$N.out"  >>$N.condor; echo "log = $N.out" >>$N.condor; echo "request_cpus = 20" >>$N.condor; echo "queue 1" >>$N.condor; chmod +x *.condor; done
```


```R
#Convert sam to bam
#!/bin/bash
#bash nombre_shipt.sh

FILES=/file_location
BIN=/samtools-1.3.1_location

COUNT=0
for FAA in `ls *.sam | sed -e 's/.sam//g'`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.sh
echo "$BIN/samtools view -bS $FILES/$FAA".sam"  > $FAA.bam"  >>$*.$COUNT.sh
chmod +x *.sh; done

for N in `ls  bS.* `; do echo "executable = $N"  >$N.condor; echo "output = $N.out" >>$N.condor; echo "error =$N.out"  >>$N.condor; echo "log = $N.out" >>$N.condor; echo "request_cpus = 20" >>$N.condor; echo "queue 1" >>$N.condor; chmod +x *.condor; done
```


```R
#Get unmaped paired sequences
#!/bin/bash
#bash nombre_shipt.sh

FILES=/file_location
BIN=/samtools-1.3.1_location

COUNT=0
for FAA in `ls *.bam | sed -e 's/.bam//g'`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.sh
echo "$BIN/samtools view -b -f 12 -F 256  $FILES/$FAA".bam"  > $FAA_"um.bam  >>$*.$COUNT.sh
chmod +x *.sh; done

#-f 12     Extract only (-f) alignments with both reads unmapped: <read unmapped><mate unmapped>
#-F 256   Do not(-F) extract alignments which are: <not primary alignment>
#see meaning of SAM-flags

for N in `ls  um.* `; do echo "executable = $N"  >$N.condor; echo "output = $N.out" >>$N.condor; echo "error =$N.out"  >>$N.condor; echo "log = $N.out" >>$N.condor; echo "request_cpus = 20" >>$N.condor; echo "queue 1" >>$N.condor; chmod +x *.condor; done
```


```R
#Obtain R1 and R2 paired sequences for host (squash) filtered samples
# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools

#!/bin/bash
#bash nombre_shipt.sh

FILES=/file_location
BIN=/samtools-1.3.1_location

COUNT=0
for FAA in `ls *_um.bam | sed -e 's/_um.bam//g'`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.sh
echo "$BIN/samtools sort -n  $FILES/$FAA"_um.bam"  > $FAA"_sorted.bam  >>$*.$COUNT.sh
chmod +x *.sh; done
```


```R
#Quality trimming
#!/bin/bash
#bash nombre_shipt.sh

SEQS=/sequences_location
BIN=/Trimmomatic-0.36_location

COUNT=0
for FAA in `ls *.fastq | perl -pe 's/\_.*//g' | sort | uniq`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.sh
echo "java -jar $BIN/trimmomatic-0.36.jar PE -threads 2 -phred33 -trimlog trim.log $SEQS/$FAA"_filtered_R1.fastq" $SEQS/$FAA"_filtered_R2.fastq" $SEQS/$FAA"_paired_R1.fastq" $S
EQS/$FAA"_unpaired_R1.fastq" $SEQS/$FAA"_paired_R2.fastq" $SEQS/$FAA"_unpaired_R2.fastq" ILLUMINACLIP:$BIN/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:1
5 MINLEN:36" >>$*.$COUNT.sh
chmod +x *.sh; done

for N in `ls  trim.* `; do echo "executable = $N"  >$N.condor; echo "output = $N.out" >>$N.condor; echo "error =$N.out"  >>$N.condor; echo "log = $N.out" >>$N.condor; echo "request_cpus = 20" >>$N.condor; echo "queue 1" >>$N.condor; chmod +x *.condor; done
```

# Metagenome assembly


```R
#Spades assembly
#!/bin/bash

SEQS=/sequences_location
BIN=/SPAdes-3.12.0_location

COUNT=0
for FAA in `ls *_paired_R*.fastq | perl -pe 's/\_.*//g' | sort | uniq`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo  "$BIN"/spades.py --meta -1 "$SEQS/$FAA"_paired_R1.fastq -2 "$SEQS/$FAA"_paired_R2.fastq -o "$SEQS/contigs_spades/contig_$FAA"/ >>$*.$COUNT.scr
chmod +x *.scr; done
```


```R
#Map unassembled sequences in megahit in order to make a second assembly using velvet. 
#!/bin/bash

CONT=/contigs_location
SEQS=/sequence_location

COUNT=0
for FAA in `ls *paired_R1* | sed -e 's/_paired_R1.fastq//g'`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo /bbwrap_location/bbwrap.sh ref="$CONT/$FAA"_contigs.fa in="$SEQS/
$FAA"_paired_R1.fastq in2="$SEQS/$FAA"_paired_R2.fastq build="$COUNT" out="$SEQS/$FAA".sam kfil
ter=22 subfilter=15 maxindel=80 >> $*.$COUNT.scr
chmod +x *.scr; done
```


```R
#Get unmapped sequences

#!/bin/bash

SEQS=/sequence_location

COUNT=0
for FAA in `ls *.sam | sed -e 's/.sam//g'`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo /home/cristobal/binc/samtools-1.2/samtools view -u -f 12 -F 256 "$SEQS/$FAA".sam "|" /home/cristobal/binc/samtools-1.2/samtools bam2fq - ">" "$FAA"_ump.fastq  >> $*.$COUNT
.scr
chmod +x *.scr; done

#Concatenate unpaired forward and revearse reads
for FAA in `ls *_unpaired_R1.fastq | perl -pe 's/\_.*//g' | sort | uniq`; do cat "$FAA"_unpaired_R1.fastq "$FAA"_unpaired_R2.fastq > "$FAA"_all_unpaired.fastq; done
```


```R
#Use velvet to assembly unassembled sequences
#!/bin/bash

SEQS=/sequence_location
BIN=/velvet_location

COUNT=0
for FAA in `ls *_ump.fastq.gz | perl -pe 's/\_.*//g' | sort | uniq`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo  "$BIN/velveth $SEQS/$FAA"_assemblyVelvet 31 -fastq -short "$SEQS/$FAA"_all_unpaired.fastq.gz -shortPaired "$SEQS/$FAA"_ump.fastq.gz >>$*.$COUNT.scr
echo "$BIN/velvetg $SEQS/$FAA"_assemblyVelvet -exp_cov 2 -ins_length 350 > $*.$COUNT.scr
chmod +x *.scr; done
```


```R
#Concatenate short and long reads
for FAA in `ls *_contigs.fasta | perl -pe 's/\_.*//g' | sort | uniq`
do cat "$FAA"_contigs.fasta /velvet_contigs_location/"$FAA"_short_contigs.fasta > /new_contigs_location/"$FAA"_hibrid_contigs.fasta
done

#Remove contigs shorter than 100pb
for N in `ls *_hybrid_contigs.fasta | sed -e 's/_hybrid_contigs.fasta//g'`; do perl removesmall.pl 100 "$N"_hybrid_contigs.fasta > "$N"_hb100_contigs.fasta; done

#Assembly statistics
python quast-5.0.0/quast.py *contigs.fasta
```


```R
#Map unassembled sequences in hybrid assembly in order to annotate unassembled genes.
#!/bin/bash

CONT=/hybrid_contigs_location/contigs_hb
SEQS=/sequence_location

COUNT=0
for FAA in `ls *_hb100_contigs.fasta | sed -e 's/_hb100_contigs.fasta//g'`
do let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo /home/cristobal/binc/bbmap/bbwrap.sh ref="$CONT/$FAA"_hb100_contigs.fasta in="$SEQS/$FAA"_ump.fastq build="$COUNT" outm="$SEQS/$FAA"_map_hb100.fastq outu="$SEQS/$FAA"_ump_
hb100.fastq kfilter=22 subfilter=15 maxindel=80 >> $*.$COUNT.scr
chmod +x *.scr; done
```


```R
#Convert unmaped sequences to fasta file and remove reads shorter than 100pb
for FAA in `ls *_ump_hb100.fastq | sed -e 's/_ump_hb100.fastq//g'`; do sed -n '1~4s/^@/>/p;2~4p' "$FAA"_ump_hb100.fastq | perl removesmall.pl 100 > "$FAA"_ump_hb100.fasta; done

#Concatenate unmapped sequences with contigs
for FAA in `ls *_ump_hb100.fasta | perl -pe 's/\_.*//g' | sort | uniq`; do cat "$FAA"_ump_hb100.fasta /hybrid_contigs_location/"$FAA"_hb100_contigs.fasta > /hybrid_contigs_location/"$FAA"_hb100seq_contigs.fasta; done


#Concatenate all trimmed sequences and convert them to fasta
for FAA in `ls *_all_unpaired.fastq | perl -pe 's/\_.*//g' | sort | uniq`; do cat "$FAA"_all_unpaired.fastq "$FAA"_paired_R1.fastq "$FAA"_paired_R2.fastq | sed -n '1~4s/^@/>/p;2~4p' > "$FAA"_trimmed.fasta; done
```

# Gene prediction


```R
#Gene prediction in prodigal
#!/bin/bash

SEQS=/hybrid_contigs_location

COUNT=0
for FAA in `ls *_hb100seq_contigs.fasta | perl -pe 's/\_.*//g' | sort | uniq`
do let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo
echo "$SEQS/$FAA"_hb100seq_contigs.fasta | perl -pe 's/\_.*//g' | sort | uniq`;  do prodigal -i "$SEQS/$FAA"_hb100seq_contigs.fast -o  "$SEQS/$FAA"_prodigal.out -d "$SEQS/$FAA".fna -a  "$SEQS/$FAA".faa -p meta
chmod +x *.scr; done

```


```R
#Rename sequences
The sequence name is changed according to the sample name with the perl script header.fasta.numbers.pl.

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
for FAA in `ls *.faa | perl -pe 's/\.faa//g' | sort | uniq`; do perl header.fasta.numbers.pl "$FAA" "$FAA".faa; done

#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
for FAA in `ls *.faa | perl -pe 's/\.faa//g' | sort | uniq`; do perl header.fasta.numbers.pl "$FAA" "$FAA".fna; done
```

## Annotation


```R
#Make diamond blast against m5nr database 
#!/bin/bash
#bash nombre_shipt.sh <nombre-del-trabajo>

SEQS=hybrid_contig_location
DB=m5nr_database_location
BIN=diamond_location

COUNT=0
for FAA in `ls *.faa`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr

echo "$BIN/diamond blastp -d $DB/m5nr -q $SEQS/$FAA" -f 6 -e 1e-10 -k 10 -p 1 -o $SEQS/
"$FAA".bout  >>$*.$COUNT.scr

chmod +x *.scr; done

#Get the best aligments based on bitscores
for FAA in `ls *.faa.bout | sed -e 's/.faa.bout//g'`; do cat "$FAA".faa.bout |  perl -pe ' $name_col=0; $score_col=11; while(<>) { s/\r?\n//; @F=split /\t/, $_; ($n, $s) = @F[$name_col, $score_col]; if (! exists($max{$n})) { push @names, $n }; if (! exists($max{$n}) || $s > $max{$n}) { $max{$n} = $s; $best{$n} = () }; if ($s == $max{$n}) { $best{$n} .= "$_\n" }; } for $n (@names) { print $best{$n} } ' > best;  perl -e ' $column=0; $unique=0; while(<>) { s/\r?\n//; @F=split /\t/, $_; if (! ($save{$F[$column]}++)) { print "$_\n"; $unique++ } } ' best > "$FAA"_best_uniq.tsv; rm best; done

#Simplify output
for FAA in `ls *.faa.bout | sed -e 's/.faa.bout//g'`; do awk '{print $1"\t"$3"\t"$11"\t"$12"\t"$2}' "$FAA"_best_uniq.tsv > "$FAA"_best.simple.tsv; done
```

## Generate tables of predicted proteins counts 


```R
#Map trimmed sequences against hybrid contigs
#Bowtie build
#!/bin/bash

BIN=/bowtie2-build_location
SEQS=/hybrid_contig_location

COUNT=0
for FAA in `ls *.fna.numbered.fas | sed -e 's/.fna.numbered.fas//g'`
do let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo "$BIN/bowtie2-build $SEQS/$FAA".fna.numbered.fas" $SEQS/$FAA" >> $*.$COUNT.scr
chmod +x *.scr; done

#Map trimmed sequences
#!/bin/bash

BIN=/bowtie2_location
SEQS=/hybrid_contig_location

COUNT=0
for FAA in `ls *_trimmed.fasta  | sed -e 's/_trimmed.fasta//g'`
do let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo "$BIN/bowtie2 -f -x $FAA -U $SEQS/$FAA"_trimmed.fasta -S "$SEQS/$FAA"_trimmed.sam --quiet -p 20 --very-sensitive >> $*.$COUNT.scr
chmod +x *.scr; done

for N in `ls *.sam | sed -e 's/.sam//g'`; do grep -v '^@' $N.sam | awk '{if($5 == "42") print $3}' | sort | uniq -c > $N.hits; done

#Make OTU list, where every predicted protein is an OTU
for FAA in `ls *.faa.numbered.fas | sed -e 's/.faa.numbered.fas//g'`; do grep ">" "$FAA".faa.numbered.fas | sed 's/#/\t/g;s/>//g' | cut -f1 | sed 's/ /\t/g' | awk 'BEGIN{i=0} /.*/{printf "%d\t% s\n",i,$0; i++}' | cut -f1,2 > "$FAA".otu; done

```


```R
#Count hits per sample
#!/bin/bash
#bash nombre_shipt.sh <nombre-del-trabajo>

SEQS=/hybrid_contig_location

COUNT=0
for FAA in `ls *_best.simple.tsv | sed -e 's/_best.simple.tsv//g'`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr

echo python3 "$SEQS"/hitter.py "$SEQS/$FAA"_best.simple.tsv "$SEQS/$FAA".hits "$SEQS/$FAA".otu "$SEQS/$FAA"  >>$*.$COUNT.scr

chmod +x *.scr; done

#Create a list of the files previously created
ls *hout > lista

#Join all tables in a list
python3 hitter_table.py lista squash

#Split file containing the checksums
cut -f1 squash.tsv > checkid

#Get subsystems, kos and refseq annotation
/m5nr_tools_location/m5nr-tools.pl --api http://api.metagenomics.anl.gov --option annotation --source Subsystems --md5 checkid > Subsystems_sq

/m5nr_tools_location/m5nr-tools.pl --api http://api.metagenomics.anl.gov --option annotation --source RefSeq --md5 checkid > RefSeq_sq

/m5nr_tools_location/m5nr-tools.pl --api http://api.metagenomics.anl.gov --option annotation --source KO --md5 checkid > KO_sq

```


```R
#Add ontology to table
perl -e ' $col1=0; $col2=0; ($f1,$f2)=@ARGV; open(F2,$f2); while (<F2>) { s/\r?\n//; @F=split /\t/, $_; $line2{$F[$col2]} .= "$_\n" }; $count2 = $.; open(F1,$f1); while (<F1>) { s/\r?\n//; @F=split /\t/, $_; $x = $line2{$F[$col1]}; if ($x) { $num_changes = ($x =~ s/^/$_\t/gm); print $x; $merged += $num_changes } } warn "\nJoining $f1 column $col1 with $f2 column $col2\n$f1: $. lines\n$f2: $count2 lines\nMerged file: $merged lines\n"; '  /subsystems_db_location/SUBSYSTEMS_29_10_18_edit.tsv Subsystems_sq > Ontology_sq

#Take of redundancy and reorder the table
awk '{print $7"\t"$2"\t"$3"\t"$4"\t"$5}' Ontology_sq > Ontology_order_sq

#Get non reduntant md5s
cut -f2 Subsystems_sq | sort | uniq > md5_sq_nr

#Get the non redundant ontology
for i in `cat md5_sq_nr` ; do grep $i Ontology_order_sq | head -n1 >> Ontology_nr_sq.tsv ; done

#There are some repeated lines that we have to remove
cut -f1 Ontology_nr_sq.tsv | sort | uniq -c | grep -w 2 | awk '{print$2}' > duplicates
for i in `cat duplicates`; do grep -n -m1 $i Ontology_nr_sq.tsv | sed -e 's/:/\t/g' | cut -f1; done > remove_from_onto
for i in `cat remove_from_onto`; do sed  - '$id' Ontology_nr_sq.tsv; done | sort | uniq > Ontology_nr_sq.tsv2

#Convert to csv. This step is necessary for an error in R, which made impossible to read the table
perl -lpe 's/"/""/g; s/^|$/"/g; s/\t/","/g' < Ontology_nr_sq.tsv2 > Ontology_nr_sq.csv
perl -lpe 's/"/""/g; s/^|$/"/g; s/\t/","/g' < squash.tsv > squash.csv
```


```R
#The list of proteins annotated in Refseq is cleaned and sorted to analyze diversity along with non annotated (hypothetical proteins)
sed  's/^/Refseq_protein\t/' RefSeq_sq > RefSeq_sq2
sed -i "s/'//g" RefSeq_sq2
sed -i 's/#//g' RefSeq_sq2
awk '{print $3"\t"$1"\t"$4"\t"$5"\t"$4"_from_"$5}' RefSeq_sq2 > Refseq.tsv

#Remove duplicated md5s in R. 
library(tidyverse)
proteins <- read.table("Refseq.tsv")
colnames(proteins) <- c("md5", "Access", "Category", "Protein", "Specie", "Protein_specie")
uproteins <- proteins[!duplicated(proteins$md5), ]
write.table(uproteins, "Refseq_nr.tsv", quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
perl -lpe 's/"/""/g; s/^|$/"/g; s/\t/","/g' < Refseq_nr.tsv > Refseq_nr.csv
#Change uppercase to lower case
tr '[:upper:]' '[:lower:]' < all_proteins.csv > all_nuproteins.csv

```


```R
#KO annotation is processed in the same way that Refseq.
sed -i.bak 's/ /_/g' KO_sq
awk '{print $2"\t"$1"\t"$3}' KO_sq > KO_sq2
perl -lpe 's/"/""/g; s/^|$/"/g; s/\t/","/g' < KO_sq2 > KO_sq.csv

```

## Add non annotated sequences to analysis


```R
#Get list of annotated sequences
for FAA in `ls *_best_uniq.tsv | sed -e 's/_best_uniq.tsv//g'`; do awk '{print $1}' "$FAA"_best_uniq.tsv > "$FAA"_anotados.tsv; done

#Get list of all predicted proteins
for FAA in `ls *.faa.numbered.fas | sed -e 's/.faa.numbered.fas//g'`; do grep '>' "$FAA".faa.numbered.fas | sed 's/>//g' |  awk '{print $1}' > "$FAA"_todos.txt; done

#Get names of non annotated sequences
for FAA in `ls *_anotados.tsv | sed -e 's/_anotados.tsv//g'`; do cat "$FAA"_anotados.tsv "$FAA"_todos.txt | sort | uniq -c | grep '1 ' | awk '{print $2}' > "$FAA"_na.txt; done

#Extract non annotated sequences
for FAA in `ls *.faa.numbered.fas | sed -e 's/.faa.numbered.fas//g'`; do /home/cristobal/binc/seqtk/seqtk subseq "$FAA".faa.numbered.fas "$FAA"_na.txt > "$FAA"_na.faa; done

#Join non annotated sequences from all samples
cat *_na.faa > todos_scna.faa

#Make clusters
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
cd-hit -i todos_scna.faa -o todos70 -c 0.70 -n 4 -aL 0.7 -d 0 -M 3000 -T 2 >  todos70.cdhit.out

#Get list of files with mapped sequences 
ls *.hits > squash_list.txt

#Convert cluster list in to an otu table
perl -pne 's/\t//g;s/^.*,//g;s/\.\.\..*$//g;s/\n/\t/g;s/\>Cluster\ /\n/g;s/\>//g; eof && do{chomp; print "$_ \n"; exit}' todos70.clstr > todos70.otu
sed -i -e "1d" todos70.ot

```


```R
#Recover the representative sequence
sed -i 's/C//' pre_hyp
perl -e ' $col1=0; $col2=0; ($f1,$f2)=@ARGV; open(F2,$f2); while (<F2>) { s/\r?\n//; @F=split /\t/, $_; $line2{$F[$col2]} .= "$_\n" }; $count2 = $.; open(F1,$f1); while (<F1>) { s/\r?\n//; @F=split /\t/, $_; $x = $line2{$F[$col1]}; if ($x) { $num_changes = ($x =~ s/^/$_\t/gm); print $x; $merged += $num_changes } } warn "\nJoining $f1 column $col1 with $f2 column $col2\n$f1: $. lines\n$f2: $count2 lines\nMerged file: $merged lines\n"; '  pre_hyp todos70.utu > representative_seqs
cut -f3 representative_seqs > rep_seqs

#Create object with id, seq and protein names
pr -mts pre_hypo rep_seqs | cat -n > no_match_proteins
awk '{print $3"\t"$2"\t"$4"\t"$2"_"$1"\t"$2"_"$1"\t"$2"_"$1}' no_match_proteins > no_match_proteins.tsv
perl -lpe 's/"/""/g; s/^|$/"/g; s/\t/","/g' < no_match_proteins.tsv > no_match_proteins.csv

#Recover annotated proteins that were not found in Refseq
cut -f1 Refseq_nr.tsv > ref_list
sed '1d' checkid > checkid2
cat checkid2 ref_list | sort | uniq -c | grep -w '1' > non_refseq_list
cat -n non_refseq_list > non_numref_list
sed -i 's/^/Matched_but_not_refseq\t/' non_numref_list
awk '{print $4"\t"$1"\t"$4"\t"$1"_"$2"\t"$1"_"$2"\t"$1"_"$2}' non_numref_list > Matched_but_not_refseq.tsv
perl -lpe 's/"/""/g; s/^|$/"/g; s/\t/","/g' < Matched_but_not_refseq.tsv > Matched_but_not_refseq.csv

#Concatenate refseq protein list, non_refseq protein list and non_annotated (hypothetical proteins).
cat Refseq_nr.csv Matched_but_not_refseq.csv na/no_match_proteins.csv > all_proteins.csv

cat squash.tsv na/squash_nan.tsv > all_squash.tsv
perl -lpe 's/"/""/g; s/^|$/"/g; s/\t/","/g' < all_squash.tsv > all_squash.csv
```

# Align against reference genomes


```R
#Download representative genomes from: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/

#Align using nucmer
for N in `cat /genomes_to_align`; do /nucmer -p "$N"_WGS_3Fr /genome_db/$N /WGS_3Fr_hybrid_contigs.fasta; done
#genomes_to_align is a list of genomes 
#This procedure must be repeated for the contigs of all metagenomes

#Filter alignments with an identity porcentage lower than 80
for N in `.delta | sed 's/.delta//g'`; do delta-filter -i 80 "$N".delta > "$N".filt

#Recover information about percetange of aligments
for N in `cat WGS_3Fr.list | sed 's/.filt//g'`; do dnadiff -p "$N" -d "$N".filt; done
#list is a list of filtered alignments in metagenome WGS_3Fr
#This procedure must be repeated for the contigs of all metagenomes

#The report files contain the information about percentage of aligned bases from each genome

```
