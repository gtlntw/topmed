#!/bin/sh

bcftools query -l /net/fantasia/home/sayantan/DATABASE/TOPMED/REFERENCE_DATA/VCF/ALL.chr22.freeze3a.pass.gtonly.genotypes.remove.missing.remove.monomorphic.Eagle.Phased.vcf.gz > id_freeze3.txt

comm id_10597.txt id_freeze3.txt -13 > id_freeze3_addittion.txt