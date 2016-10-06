#!/bin/sh

# generate vcf files #-i /net/topmed2/working/khlin/subset_id/subset.sftp-barnes.ids
# for chr in `seq 1 22`;
# do
# 	nohup python /net/topmed2/working/khlin/utilities/lai_2_vcf_hgdp_topmed.py -i /net/topmed2/working/khlin/subset_id/subset.sftp-barnes.ids -c ${chr} > chr${chr}.lai.vcf.txt &
# done

# #bgzip
# for chr in `seq 1 22`;
# do
# 	bgzip  chr${chr}.lai.vcf &
# done

# #tabix
# for chr in `seq 1 22`;
# do
# 	bcftools index -t -f chr${chr}.lai.vcf.gz &
# done

##copy into directory
# cp chr*.lai.vcf.gz barnes/ &
# cp chr*.lai.vcf.gz.tbi barnes/ &

##rename 
# for chr in `seq 1 22`;
# do
# 	mv chr${chr}.lai.vcf.gz topmed_freeze2.chr${chr}.local_ancestry.sftp-barnes.vcf.gz
# 	mv chr${chr}.lai.vcf.gz.tbi topmed_freeze2.chr${chr}.local_ancestry.sftp-barnes.vcf.gz.tbi
# done

## cp to freeze2a   cd /net/topmed2/working/gt-release/sftp-barnes/freeze.2a/local_ancestry/
# cp /net/topmed2/working/khlin/vcf/barnes/*.* ./

##md5sum
# md5sum *.* > 0_md5sum.txt

##chgrp
chgrp topmed *.*

