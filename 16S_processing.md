# Amplicon sequence processing


```R
Raw sequences were quality checked with fastx_tools and assembled with Pandaseq.
```


```R
A list of files to assembly is generated.
$ ls *.fastq.gz | perl -pe 's/_R.*.fastq.gz//g' | sort | uniq >lista
```


```R
# assembly.sh

#!/bin/bash

SEQS=$(pwd)
SALIDAS=$(pwd)
BIN=/usr/bin
BIN2=/usr/bin
COUNT=0
for FAA in `cat lista`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo "zcat $SEQS/$FAA"_R1.fastq.gz" | $BIN2/fastx_trimmer -l 250 > $SEQS/$FAA"tr_R1.fastq" &" >>$*.$COUNT.scr mo
echo "zcat $SEQS/$FAA"_R2.fastq.gz" | $BIN2/fastx_trimmer -l 250 > $SEQS/$FAA"tr_R2.fastq"" >>$*.$COUNT.scr
echo "$BIN/pandaseq -B -f $SEQS/$FAA"tr_R1.fastq" -r $SEQS/$FAA"tr_R2.fastq" -t 0.95 -l 250 -L 470 -o 15 -w $SALIDAS/$FAA"_$*.fasta" -G $SALIDAS/$FAA"-$*.log.bz2"" >>$*.$COUNT.scr

done
```


```R
Preceding script utilices “fastx_trimmer” that crop sequence to a size that is defined with -l. The obtained file is assembled with “pandaseq”, -t can take values for 0 to 1 and alignments with lower values are discarded; -l is the minimum sequence length -L maximum sequence length; -o minimum overlap between assembled sequences.
```


```R
Sequences names are changed with “header.fasta.numbers.pl”.

# Luis David Alcaraz 2013-04-11

my $prefix = $ARGV[0]; chomp $prefix;
my $f =  1;

my $fasta_file = $ARGV [1]; chomp $fasta_file;

my $fh;
open($fh, $fasta_file) or die "can't open $fasta_file: $!\n";
open(OUT, ">$fasta_file.numbered.fas") || die "can't open $fasta_file.numbered.f
as\n";

my %sequence_data;
while (read_fasta_sequence($fh, \%sequence_data)) {
   print OUT ">$sequence_data{header}\n$sequence_data{seq}\n";
}

close $fh;
close OUT;

sub read_fasta_sequence {
   my ($fh, $seq_info) = @_;

$seq_info->{seq} = undef; # clear out previous sequence

   # put the header into place
   $seq_info->{header} = $seq_info->{next_header} if $seq_info->{next_header};

   my $file_not_empty = 0;
   while (<$fh>) {
      $file_not_empty = 1;
      next if /^\s*$/;  # skip blank lines
      chomp;    

      if (/^>/) { # fasta header line
         my $h = $_;    
         $h =~ s/>/$prefix\_$f\ /;
     $f++;
         if ($seq_info->{header}) {
            $seq_info->{next_header} = $h;
            return $seq_info;
       
         }              
         else { # first time through only
            $seq_info->{header} = $h;
         }              
      }         
      else {    
         s/\s+//;  # remove any white space
         s/\n\n/\n/;
         $seq_info->{seq} .= $_;
      }         
   }    

   if ($file_not_empty) {
      return $seq_info;
   }    
   else {
      # clean everything up
      $seq_info->{header} = $seq_info->{seq} = $seq_info->{next_header} = undef;

      return;   
   }    
}
```


```R
Each sequence is renamed according to the sample name, followed by a number in consecutive order.
$ perl ./header.fasta.numbers.pl nombre_del_archivo nombre_del_archivo
The previous command must be applied for every assembled sample.
```


```R
Concatenate all samples
$ cat *.numbered.fas > acalabacita.fas
```


```R
Edit sequence name to leave only the first part, which refers to sample name. 
$ perl -i.bak -pe "s/\ .*//g" acalabacita.fas
```


```R
Sequence clustering of OTUs at 97%  of identity was done with cd-hit-est, -c indicates the identity for clustering.
$ cd-hit-est -i calabacita.fas -c 0.97 -o acalabacita.fasout -T 20 -M 0
```


```R
The clustering file is converted in a file that can be read by QIIME.
$ perl -pne 's/\t//g;s/^.*,//g;s/\.\.\..*$//g;s/\n/\t/g;s/\>Cluster\ /\n/g;s/\>//g; eof && do{chomp; print "$_ \n"; exit}' acalabacita.fasout.clstr > acalabacita.otu
$ sed -i '1d' acalabacita.otu
```

# Taxonomic assignment


```R
Taxonomic assignment was done with QIIME using Greengenes database as reference.
```


```R
Extract representative OTUs, -i indicates input file, -f is the fasta file with the sequences to be extracted. 
$ pick_rep_set.py -i acalabacita.otu -f calabacita.fas -o acalabacita.rep.fna
```


```R
Taxonomic assignment was done in parallel, -i is the fasta file with representative sequences, -o is the directory for output, -r the reference sequences for alignments.
$ parallel_assign_taxonomy_blast.py -i acalabacita.rep.fna -o taxonomy -r /qiime/gg_otus-13_8-release/rep_set/97_otus.fasta
```


```R
Make a file with the elements to filter, includes no hits, chloroplast and mitochondria sequences.
$ grep -i "mitochondria" taxonomy/biof_tax_assignments.txt  >> to_remove_from_biom
$ grep -i "Chloroplast" taxonomy/biof_tax_assignments.txt  >> to_remove_from_biom
$ grep -i "no blast hit" taxonomy/biof_tax_assignments.txt  >> to_remove_from_biom
```


```R
Make OTU table with the number of times that every OTU appears in the samples, -i indicates the OTU file, -t taxonomic assignments, -e OTUs to be removed.
$ make_otu_table.py -i acalabacita.otu -t taxonomy/acalabacita.rep_tax_assignments.txt -e to_remove_from_biom 
-o acalabacita.biom
```


```R
Singleton filtering, -n indicates that OTUs with a lower abundance will be removed. 
$ filter_otus_from_otu_table.py -i acalabacita.biom -o nscalabacita.biom -n 2
```


```R
Convert OTU table to tabular format.
$ biom convert --to-tsv -i nscalabacita.biom -o nscalabacita.txt --table-type "Taxon table" --header-key=taxonomy
```


```R
Split taxonomy table from OTU table.
$ perl -pe 's/\; /\;/g' nscalabacita.txt | awk '{print $1,"\t",$NF}' | perl -pe 's/\;/\t/g' > Silva_taxonomy.tsv
```


```R
Split OTU abundance for sample from OTU table: 
$ cut -f 1-48 nscalabacita.txt > otusquash.txt
```


```R
Filter fasta to remove no hits, chloroplast and mitochondria sequences, -f is the fasta file, -b the previously filtered biom file. 
$ filter_fasta.py -f acalabacita.rep.fna -b nscalabacita.biom -o nscalabacita.rep.fna
```


```R
Update taxonomy table with Silva data base
$ R
> calabaza <- readDNAStringSet("calabacita_repset2020.fna") 
> seqs <- getSequences(calabaza)  

> taxa <- assignTaxonomy(seqs, "silva_nr99_v138_train_set.fa.gz", multithread=40); save.image()
> asv_seqs <- colnames(calabaza)
> asv_headers <- names(calabaza)


> asv_tax <- taxa
> row.names(asv_tax) <- sub(">", "", asv_headers)
> write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)
> savehistory()
> q()
```


```R
Remove chimeras
$ identify_chimeric_seqs.py -i nscalabacita.rep_aligned.fasta -t 97_otu_taxonomy.txt -r 97_otus.fasta 
-o cblast.txt -m blast_fragments

Make list of non-chimeric sequences
$ awk '{print $1}' cblast.txt > chim_list

$ cut -f1  Silva_taxonomy.tsv > all_otus

$ cat all_otus chim_list | sort | uniq -c | grep -w '1' | awk '{print$2}' > no_chim_list

$ R
> silva <- read.table("Silva_taxonomy.tsv", header=T, sep="\t", row.names=1)
> no_chim <- read.table("no_chim_list", row.names=1)
> nc_silva <- merge(no_chim, silva, by=0)
> write.table(nc_silva, "nc_silva_taxonomy.tsv", sep="\t", row.names = FALSE, quote = FALSE)
> q()

$ sed -i 's/NA/unclassified/g' nc_silva_taxonomy.tsv

$ sed -i 's/ /_/g' nc_silva_taxonomy.tsv 
```


```R
Remove chimeras from filtered sequences
$ seqtk subseq nscalabacita.rep.fna no_chim_list > nc_nscalabacita.rep.fna

Alignment of filtered sequences

$ parallel_align_seqs_pynast.py -i nc_nscalabacita.rep.fna  -o nc_calabacita.align
```


```R
Phylogenetic tree of aligned representative sequences, -nt indicates nucleotide sequences, -gtr substitution model GTR.

$ fasttree -nt -gtr nc_calabacita.align/nc_nscalabacita.rep_aligned.fasta > nc_squash.tree
```
