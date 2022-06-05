#!/bin/bash

set -x 
set -euo pipefail

# Benchmark
# wget https://dataverse.harvard.edu/api/access/datafile/6176325 -O DNA_2_150_qer_8x.fasta.gz
# wget https://dataverse.harvard.edu/api/access/datafile/6176326 -O DNA_2_150_ref_8x.fasta.gz
# wget https://dataverse.harvard.edu/api/access/datafile/6176327 -O DNA_2_200_qer_128x.fasta.gz
# wget https://dataverse.harvard.edu/api/access/datafile/6176328 -O DNA_2_200_ref_128x.fasta.gz
# wget https://dataverse.harvard.edu/api/access/datafile/6176329 -O DNA_2_250_qer_32x.fasta.gz
# wget https://dataverse.harvard.edu/api/access/datafile/6176330 -O DNA_2_250_ref_32x.fasta.gz
# 
# wget https://dataverse.harvard.edu/api/access/datafile/6176307 -O DNA-big-As.fasta.gz.part-aa
# wget https://dataverse.harvard.edu/api/access/datafile/6176308 -O DNA-big-As.fasta.gz.part-ab
# wget https://dataverse.harvard.edu/api/access/datafile/6176309 -O DNA-big-As.fasta.gz.part-ac
# wget https://dataverse.harvard.edu/api/access/datafile/6176312 -O DNA-big-As.fasta.gz.part-ad
# wget https://dataverse.harvard.edu/api/access/datafile/6176313 -O DNA-big-As.fasta.gz.part-ae
# wget https://dataverse.harvard.edu/api/access/datafile/6176316 -O DNA-big-As.fasta.gz.part-af
# wget https://dataverse.harvard.edu/api/access/datafile/6176317 -O DNA-big-As.fasta.gz.part-ag
# cat DNA-big-As.fasta.gz.part-* > DNA-big-As.fasta.gz
# 
# wget https://dataverse.harvard.edu/api/access/datafile/6176318 -O DNA-big-Bs.fasta.gz.part-aa
# wget https://dataverse.harvard.edu/api/access/datafile/6176319 -O DNA-big-Bs.fasta.gz.part-ab
# wget https://dataverse.harvard.edu/api/access/datafile/6176320 -O DNA-big-Bs.fasta.gz.part-ac
# wget https://dataverse.harvard.edu/api/access/datafile/6176321 -O DNA-big-Bs.fasta.gz.part-ad
# wget https://dataverse.harvard.edu/api/access/datafile/6176322 -O DNA-big-Bs.fasta.gz.part-ae
# wget https://dataverse.harvard.edu/api/access/datafile/6176323 -O DNA-big-Bs.fasta.gz.part-af
# cat DNA-big-Bs.fasta.gz.part-* > DNA-big-Bs.fasta.gz
# 
# wget https://dataverse.harvard.edu/api/access/datafile/6176332 -O PROTEIN_200_que.fasta.gz
# wget https://dataverse.harvard.edu/api/access/datafile/6176333 -O PROTEIN_200_ref.fasta.gz
# wget https://dataverse.harvard.edu/api/access/datafile/6176334 -O PROTEIN_400_que.fasta.gz
# wget https://dataverse.harvard.edu/api/access/datafile/6176335 -O PROTEIN_400_ref.fasta.gz
# wget https://dataverse.harvard.edu/api/access/datafile/6176336 -O PROTEIN_600_que.fasta.gz
# wget https://dataverse.harvard.edu/api/access/datafile/6176337 -O PROTEIN_600_ref.fasta.gz
# 
# wget https://dataverse.harvard.edu/api/access/datafile/6176338 -O PROTEIN-longer_que.fasta.gz
# wget https://dataverse.harvard.edu/api/access/datafile/6176339 -O PROTEIN-longer_ref.fasta.gz
wget https://dataverse.harvard.edu/api/access/datafile/6317573 -O As_new.fasta.gz
wget https://dataverse.harvard.edu/api/access/datafile/6317574 -O Bs_new.fasta.gz



# PASTIS
wget https://dataverse.harvard.edu/api/access/datafile/6176344 -O metaclust50_500K.fasta.gz


# Unzip everything
unpigz *.gz || gunzip *.gz


for d in *.fasta ; do
	cat $d | awk '(NR+1)%2==1' > ${d%.fasta}.txt
done

