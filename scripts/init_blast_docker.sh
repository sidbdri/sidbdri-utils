#!/usr/bin/env bash

#This only need to be run when we want to update the blast database!!!

#https://github.com/ncbi/docker/tree/master/blast
#https://github.com/ncbi/docker/tree/master/blast#install-ncbi-provided-blast-databases
docker pull ncbi/blast
HOST_ROOT=/srv/data/blast
cd ${HOST_ROOT}; mkdir blastdb queries fasta results blastdb_custom


#docker run --rm ncbi/blast update_blastdb.pl --showall --source ncbi
#docker run --rm ncbi/blast update_blastdb.pl --showall pretty --source gcp
#update_blastdb.pl --showall [*]


##This will take a while!! about 70g of data will be downloaded.
docker run --rm \
  -v ${HOST_ROOT}/blastdb:/blast/blastdb:rw \
  -w /blast/blastdb \
  ncbi/blast \
  update_blastdb.pl --blastdb_version 5 --decompress --source ncbi nt_v5 &

wait


#nohup docker run --rm \
#  -v ${HOST_ROOT}/blastdb:/blast/blastdb:rw \
#  -w /blast/blastdb \
#  ncbi/blast \
#  update_blastdb.pl --decompress --source gcp nt_v5 &


docker run --rm \
  -v ${HOST_ROOT}/blastdb:/blast/blastdb:rw \
  -w /blast/blastdb \
  ncbi/blast \
  makembindex -input nt_v5 &
wait

#docker run --rm -v ${HOST_ROOT}/blastdb:/blast/blastdb:ro ncbi/blast blastdbcmd -list /blast/blastdb -remove_redundant_dbs




docker run --rm -dit --name blast \
  -v ${HOST_ROOT}/blastdb:/blast/blastdb:ro -v ${HOST_ROOT}/blastdb_custom:/blast/blastdb_custom:ro \
  -v ${HOST_ROOT}/queries:/blast/queries:ro \
  -v ${HOST_ROOT}/results:/blast/results:rw \
  ncbi/blast \
  sleep infinity

#
#docker run --rm \
#  -v ${HOST_ROOT}/blastdb:/blast/blastdb:rw \
#  -w /blast/blastdb \
#  ncbi/blast \
#  blastdbcmd -db nt_v5 -tax_info
#  		    %a means accession
#   		%g means gi
#   		%o means ordinal id (OID)
#   		%i means sequence id
#   		%t means sequence title
#   		%l means sequence length
#   		%h means sequence hash value
#   		%T means taxid
#   		%X means leaf-node taxids
#   		%e means membership integer
#   		%L means common taxonomic name
#   		%C means common taxonomic names for leaf-node taxids
#   		%S means scientific name
#   		%N means scientific names for leaf-node taxids
#   		%B means BLAST name
#   		%K means taxonomic super kingdom
#   		%P means PIG
#   		%m means sequence masking data.
#docker run --rm \
#  -v ${HOST_ROOT}/blastdb:/blast/blastdb:rw \
#  -w /blast/blastdb \
#  ncbi/blast \
#  blastdbcmd -db nt_v5 -taxids 9606,10090,10116 -outfmt '%a;%g;%o;%i;%t;%l;%h;%T;%X;%e;%L;%C;%S;%N;%B;%K;%P' | head
#
#blastdbcmd -db nt_v5 -taxids 9606 -outfmt '%a;%g;%o;%i;%t;%l;%h;%T;%X;%e;%L;%C;%S;%N;%B;%K;%P' | head
#blastdbcmd -db nt_v5 -taxids 10090 -outfmt '%a;%g;%o;%i;%t;%l;%h;%T;%X;%e;%L;%C;%S;%N;%B;%K;%P' | head
#blastdbcmd -db nt_v5 -taxids 10116 -outfmt '%a;%g;%o;%i;%t;%l;%h;%T;%X;%e;%L;%C;%S;%N;%B;%K;%P' | head





#goinside blast root
#blastn -task blastn -db nt
#
#
#docker exec blast /bin/bash -c "echo AGCATGCTTGAAATGGCTGCGTTGTTGCCGGGAACTCGGTAAATACCTGTATATTCAAGACCCCTTTCTTCAACCAATTTGCAACATATGTCAACGATTA | blastn -db nt"
#docker exec blast /bin/bash -c "echo AGCATGCTTGAAATGGCTGCGTTGTTGCCGGGAACTCGGTAAATACCTGTATATTCAAGACCCCTTTCTTCAACCAATTTGCAACATATGTCAACGATTA | \
#blastn -db nt -outfmt 6 -out /blast/results/test6.out"
#
##         staxid means Subject Taxonomy ID
##   	  ssciname means Subject Scientific Name
##   	  scomname means Subject Common Name
##   	sblastname means Subject Blast Name
#docker exec blast /bin/bash -c "echo AGCATGCTTGAAATGGCTGCGTTGTTGCCGGGAACTCGGTAAATACCTGTATATTCAAGACCCCTTTCTTCAACCAATTTGCAACATATGTCAACGATTA | \
#blastn -db nt -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
#-out /blast/results/test6.out -remote -entrez_query 'rat OR mouse OR human'"
#
#
## Using index make it much faster!
#docker exec blast /bin/bash -c "echo AGCATGCTTGAAATGGCTGCGTTGTTGCCGGGAACTCGGTAAATACCTGTATATTCAAGACCCCTTTCTTCAACCAATTTGCAACATATGTCAACGATTA | \
#blastn -db nt -out /blast/results/test0.out" &
#docker exec blast /bin/bash -c "echo AGCATGCTTGAAATGGCTGCGTTGTTGCCGGGAACTCGGTAAATACCTGTATATTCAAGACCCCTTTCTTCAACCAATTTGCAACATATGTCAACGATTA | \
#blastn -db nt -use_index true -out /blast/results/test1.out" &
#docker exec blast /bin/bash -c "blastn -db nt -query /blast/queries/test20.fasta -out /blast/results/test2.out" &
#docker exec blast /bin/bash -c "blastn -db nt -use_index true -num_threads 20 -query /blast/queries/test20.fasta" &
#
## sed -n '1~4s/^@/>/p;2~4p' test.fastq > test.fasta
#docker exec blast /bin/bash -c "blastn -db nt -use_index true -query /blast/queries/test.fasta -out /blast/results/test.out" &
#
#
#
#
#docker exec blast /bin/bash -c "echo AGCATGCTTGAAATGGCTGCGTTGTTGCCGGGAACTCGGTAAATACCTGTATATTCAAGACCCCTTTCTTCAACCAATTTGCAACATATGTCAACGATTA | \
#blastn -db nt -use_index true -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
#-out /blast/results/test6.out"
#
#docker exec blast /bin/bash -c "echo NCAGGTGGGCCTTCCCGGCCGTCCCGGAGCCGGTCGCGGCGCACCGCCACGGTGGAAGTGCGCCCGGCGGCGGCCGGTCGCCGGCCGGGGGGCGGTCCCC | \
#blastn -db nt -use_index true -max_target_seqs 5 -outfmt " &
#
#docker exec blast /bin/bash -c "blastn -db nt -use_index true -max_target_seqs 10 -outfmt '6 qaccver saccver staxid ssciname scomname pident length mismatch gapopen qstart qend sstart send evalue bitscore' -query /blast/queries/test1.fasta" &
#




#-num_threads 10
#
# -num_descriptions <Integer, >=0>
#   Number of database sequences to show one-line descriptions for
#   Not applicable for outfmt > 4
#   Default = `500'
#    * Incompatible with:  max_target_seqs
# -num_alignments <Integer, >=0>
#   Number of database sequences to show alignments for
#   Default = `250'
#    * Incompatible with:  max_target_seqs
#
# -max_target_seqs <Integer, >=1>
#   Maximum number of aligned sequences to keep
#   (value of 5 or more is recommended)
#   Default = `500'
#    * Incompatible with:  num_descriptions, num_alignments
