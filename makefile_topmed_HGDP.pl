#!/usr/bin/perl -w
#adapted the original perl script from Adrian Tan
 
use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;
 
=head1 NAME
 
generate_simple_stuff_makefile
 
=head1 SYNOPSIS
 
 generate_simple_stuff_makefile [options]
 
  -o     output directory : location of all output files
  -l     local or slurm mode
  -m     output make file
  -j     jobname
  -i     sample id
 
 example: ./generate_simple_stuff_makefile.pl
 
=head1 DESCRIPTION
 
=cut
 
#option variables
my $help;
my $verbose;
my $debug;
my $outputDir = "/net/topmed2/working/khlin";
my $makeFile = "makefile";
my $launchMethod = "local";
my $id = "";
my $jobName = "runmake.time.log";

#initialize options
Getopt::Long::Configure ('bundling');
 
if(!GetOptions ('h'=>\$help, 'v'=>\$verbose, 'd'=>\$debug,
                'o:s'=>\$outputDir,
                'l:s'=>\$launchMethod,
                'm:s'=>\$makeFile,
                'j:s'=>\$jobName,
                'i:s'=>\$id)
  || !defined($outputDir)
  || scalar(@ARGV)!=0)
{
    if ($help)
    {
        pod2usage(-verbose => 2);
    }
    else
    {
        pod2usage(1);
    }
}
 
if ($launchMethod ne "local" && $launchMethod ne "slurm")
{
    print STDERR "Launch method has to be local or slurm\n";
    exit(1);
}
 
##############
#print options
##############
printf("Options\n");
printf("\n");
printf("sample ID        : %s\n", $id);
printf("output directory : %s\n", $outputDir);
printf("launch method    : %s\n", $launchMethod);
printf("output jobName   : %s\n", $jobName);
printf("makefile         : %s\n", $makeFile);
printf("\n");

 
#arrays for storing targets, dependencies and commands
my @tgts = ();
my @deps = ();
my @cmds = ();
 
#temporary variables
my $tgt;
my $dep;
my @cmd;
 
mkpath($outputDir);
mkpath("output");
mkpath("log");
mkpath("LAI/$id");
my $slurmScriptsDir = "$outputDir/slurm_scripts";
mkpath($slurmScriptsDir);
my $slurmScriptNo = 0;
my $toolsDir = '/net/snowwhite/home/khlin/tools';

my $inputFiles = "";
my $inputFilesOK = "";
my $inputFile = "";
my $outputFile = "";
my $param = ""; #parameters for slurm when using shell script


######################
######################
#0.1. log the start time
######################

$tgt = "$outputDir/log/start.runmake.${id}_shapeit.OK";
$dep = "";
@cmd = ("date | awk '{print \"Local ancestry pipeline\\n\\nstart: \"\$\$0}' > $outputDir/log/runmake_${id}_shapeit_time.log");
makeJob("local", $tgt, $dep, @cmd);

######################
#0.2. extract chromosomes of sample with bi-allelic snps
######################
for my $chr (1..22)
{
    $tgt = "${outputDir}/LAI/${id}/${id}_filtered_chr${chr}.vcf.gz.OK";
    $dep = "";
    @cmd = ("bcftools view /net/topmed2/working/gt-release/sftp-barnes/freeze.1a/topmed.freeze1.nhlbi.791.sftp-barnes.keep.chr${chr}.gtonly.vcf.gz \\
        --types snps -M2 --exclude-uncalled -f PASS --regions ${chr} --force-samples -s ${id} \\
        --output-type z --output-file ${outputDir}/LAI/${id}/${id}_filtered_chr${chr}.vcf.gz && \\
    bcftools index -t -f ${outputDir}/LAI/${id}/${id}_filtered_chr${chr}.vcf.gz");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

# ######################
# #0.2. remove duplicate snps and keep only bi-alleic
# ######################
# for my $chr (1..22)
# {
#     $tgt = "${outputDir}/LAI/${id}/${id}_filtered_phased_norm_chr${chr}.vcf.gz.OK";
#     $dep = "${outputDir}/LAI/${id}/${id}_filtered_chr${chr}.vcf.gz.OK";
#     @cmd = ("bcftools norm -m+ ${outputDir}/LAI/${id}/${id}_filtered_phased_chr${chr}.vcf.gz -Ou | bcftools view -M2 -Oz -o ${outputDir}/LAI/${id}/${id}_filtered_phased_norm_chr${chr}.vcf.gz");
#     makeJob($launchMethod, $tgt, $dep, @cmd);
# }

######################
#0.3. run shapeit to phase and exclude missing snps if necessary
######################
$param="--mem=20000"; # make srun sh to specify memomry

for my $chr (1..22)
{
    $tgt = "${outputDir}/LAI/${id}/${id}_phased_chr${chr}.OK";
    $dep = "${outputDir}/LAI/${id}/${id}_filtered_chr${chr}.vcf.gz.OK";
    @cmd = ("/net/snowwhite/home/khlin/tools/shapeit.v2.r837/bin/shapeit -phase \\
        -V ${outputDir}/LAI/${id}/${id}_filtered_chr${chr}.vcf.gz \\
        -M ${outputDir}/genetic_map_GRCh37/genetic_map_chr${chr}_combined_b37.txt \\
        --input-ref ${outputDir}/1000GP_Phase3/1000GP_Phase3_chr${chr}.hap.gz ${outputDir}/1000GP_Phase3/1000GP_Phase3_chr${chr}.legend.gz ${outputDir}/1000GP_Phase3/1000GP_Phase3.sample \\
        -O ${outputDir}/LAI/${id}/${id}_phased_chr${chr} \\
        --output-log ${outputDir}/LAI/${id}/${id}_phased_chr${chr}.Output || \\
        /net/snowwhite/home/khlin/tools/shapeit.v2.r837/bin/shapeit -phase \\
        -V ${outputDir}/LAI/${id}/${id}_filtered_chr${chr}.vcf.gz \\
        -M ${outputDir}/genetic_map_GRCh37/genetic_map_chr${chr}_combined_b37.txt \\
        --input-ref ${outputDir}/1000GP_Phase3/1000GP_Phase3_chr${chr}.hap.gz ${outputDir}/1000GP_Phase3/1000GP_Phase3_chr${chr}.legend.gz ${outputDir}/1000GP_Phase3/1000GP_Phase3.sample \\
        --exclude-snp ${outputDir}/LAI/${id}/${id}_phased_chr${chr}.Output.snp.strand.exclude \\
        -O ${outputDir}/LAI/${id}/${id}_phased_chr${chr} \\
        --output-log ${outputDir}/LAI/${id}/${id}_phased_chr${chr}.Output");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}
$param=""; #clean param

######################
#0.4. convert back to vcf.gz and clean intermediate files
######################
$inputFiles=""; #clean up
$inputFilesOK=""; #clean up

for my $chr (1..22)
{
    $tgt = "${outputDir}/LAI/${id}/${id}_filtered_phased_chr${chr}.vcf.gz.OK";
    $inputFilesOK .= " ${tgt}";
    $dep = "${outputDir}/LAI/${id}/${id}_phased_chr${chr}.OK";
    @cmd = ("/net/snowwhite/home/khlin/tools/shapeit.v2.r837/bin/shapeit -convert \\
        --input-haps ${outputDir}/LAI/${id}/${id}_phased_chr${chr} \\
        --output-vcf ${outputDir}/LAI/${id}/${id}_filtered_phased_chr${chr}.vcf \\
        --output-log ${outputDir}/LAI/${id}/${id}_filtered_phased_chr${chr}.vcf && \\
        bgzip -f ${outputDir}/LAI/${id}/${id}_filtered_phased_chr${chr}.vcf &&\\
        tabix -f -p vcf ${outputDir}/LAI/${id}/${id}_filtered_phased_chr${chr}.vcf.gz");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#0.5. clean intermediate files
######################
for my $chr (1..22)
{
    $tgt = "${outputDir}/LAI/${id}/${id}_rm_intermediate_files_chr${chr}.OK";
    $dep = "$inputFilesOK";
    @cmd = ("rm ${outputDir}/LAI/${id}/${id}_phased_chr${chr}* && \\
        rm ${outputDir}/LAI/${id}/${id}_filtered_chr${chr}*");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}


######################
#0.6. log end time
######################
$tgt = "$outputDir/log/end.runmake.${id}_shapeit.OK";
$dep = "$inputFilesOK";
@cmd = ("date | awk '{print \"\\nend: \"\$\$0}' >> $outputDir/log/runmake_${id}_shapeit_time.log");
makeJob("local", $tgt, $dep, @cmd);



#############################################
#RFMix local ancestry inference using all 7 ancestral populations
#############################################

######################
#1.0. filtered HGDP and keep only bi-allelic
######################
for my $chr (1..22)
{
    $tgt = "$outputDir/HGDP_938/HGDP_938_chr${chr}_filtered_phased.vcf.gz.OK";
    $dep = "";
    @cmd = ("bcftools view $outputDir/HGDP_938/HGDP_938_chr${chr}_phased.vcf.gz \\
        --types snps -M2 --exclude-uncalled -f PASS \\
        --output-type z --output-file $outputDir/HGDP_938/HGDP_938_chr${chr}_filtered_phased.vcf.gz; \\
    bcftools index -t -f $outputDir/HGDP_938/HGDP_938_chr${chr}_filtered_phased.vcf.gz");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

# ######################
# #1.1. ind. common snp between sample and HGDP
# ######################
# for my $chr (1..22)
# {
#     $tgt = "$outputDir/common_site/topmed_chr${chr}_HGDP_common.txt.OK";
#     $dep = "${outputDir}/LAI/${id}/${id}_filtered_phased_chr${chr}.vcf.gz.OK $outputDir/HGDP_938/HGDP_938_chr${chr}_filtered_phased.vcf.gz.OK";
#     @cmd = ("/net/snowwhite/home/khlin/bin/vcftools --gzvcf ${outputDir}/LAI/${id}/${id}_filtered_phased_chr${chr}.vcf.gz --gzdiff $outputDir/HGDP_938/HGDP_938_chr${chr}_filtered_phased.vcf.gz \\
#         --diff-site --stdout | awk '{if(\$4 == \"B\")  print \$1 \"\\t\" \$2}' > $outputDir/common_site/topmed_chr${chr}_HGDP_common.txt");
#     makeJob($launchMethod, $tgt, $dep, @cmd);
# }

######################
#1.2. log the start time
######################

$tgt = "$outputDir/log/start.runmake.${id}_HGDP_RFMix.OK";
$dep = "";
@cmd = ("date | awk '{print \"Local ancestry pipeline\\n\\nstart: \"\$\$0}' > $outputDir/log/runmake_${id}_HGDP_RFMix_time.log");
makeJob("local", $tgt, $dep, @cmd);

######################
#1.3. merge sample and HGDP
######################
for my $chr (1..22)
{
    $tgt = "${outputDir}/LAI/${id}/${id}_HGDP_chr${chr}_filtered_phased.vcf.gz.OK";
    $dep = "${outputDir}/LAI/${id}/${id}_filtered_phased_chr${chr}.vcf.gz.OK $outputDir/HGDP_938/HGDP_938_chr${chr}_filtered_phased.vcf.gz.OK";
    @cmd = ("bcftools merge -O z -o ${outputDir}/LAI/${id}/${id}_HGDP_chr${chr}_filtered_phased.vcf.gz \\
        ${outputDir}/LAI/${id}/${id}_filtered_phased_chr${chr}.vcf.gz $outputDir/HGDP_938/HGDP_938_chr${chr}_filtered_phased.vcf.gz && \\
        bcftools index -t -f ${outputDir}/LAI/${id}/${id}_HGDP_chr${chr}_filtered_phased.vcf.gz");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#1.4. make files needed for RFMix allele, location, classes
######################
for my $chr (1..22)
{
    $tgt = "${outputDir}/LAI/${id}/${id}_HGDP_chr${chr}_RFMix_files_prep.OK";
    $dep = "${outputDir}/LAI/${id}/${id}_HGDP_chr${chr}_filtered_phased.vcf.gz.OK";
    @cmd = ("bcftools view ${outputDir}/LAI/${id}/${id}_HGDP_chr${chr}_filtered_phased.vcf.gz \\
      | bcftools query -f '[%GT]\\n' - | sed 's/|//g' > ${outputDir}/LAI/${id}/${id}_HGDP_chr${chr}.alleles && \\
    bcftools view ${outputDir}/LAI/${id}/${id}_HGDP_chr${chr}_filtered_phased.vcf.gz \\
        | bcftools query -f '%POS\\n' - | Rscript $outputDir/utilities/generate_markerLocations_file.R stdin $outputDir/genetic_map_GRCh37/genetic_map_chr${chr}_combined_b37.txt ${outputDir}/LAI/${id}/${id}_HGDP_chr${chr}.locations && \\
    bcftools query -l ${outputDir}/LAI/${id}/${id}_HGDP_chr${chr}_filtered_phased.vcf.gz | Rscript $outputDir/utilities/make_classes_file.R \\
        stdin ${outputDir}/LAI/${id}/${id}_HGDP_chr${chr}.classes");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#1.5. run RFMix
######################
$inputFiles=""; #clean up
$inputFilesOK=""; #clean up
for my $chr (1..22)
{
    $inputFiles .= " ${outputDir}/LAI/${id}/${id}_HGDP_chr${chr}.0.Viterbi.txt";
    $inputFiles .= " $outputDir/common_site/topmed_chr${chr}_HGDP_common.txt";
    $inputFilesOK .= " ${outputDir}/LAI/${id}/${id}_HGDP_chr${chr}_RFMix_run.OK";
    $tgt = "${outputDir}/LAI/${id}/${id}_HGDP_chr${chr}_RFMix_run.OK";
    $dep = "${outputDir}/LAI/${id}/${id}_HGDP_chr${chr}_RFMix_files_prep.OK";
    @cmd = ("cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && \\
    python ./RunRFMix.py PopPhased ${outputDir}/LAI/${id}/${id}_HGDP_chr${chr}.alleles \\
        ${outputDir}/LAI/${id}/${id}_HGDP_chr${chr}.classes \\
        ${outputDir}/LAI/${id}/${id}_HGDP_chr${chr}.locations \\
        -o ${outputDir}/LAI/${id}/${id}_HGDP_chr${chr} --forward-backward --num-threads 1");
    makeJob($launchMethod, $tgt, $dep, @cmd);
}

######################
#1.6. make local ancestry plot
######################
$outputFile = "${outputDir}/output/${id}_HGDP_local_ancestry_RFMix.png";
$tgt = "${outputDir}/output/${id}_HGDP_local_ancestry_RFMix.png.OK";
$dep = "$inputFilesOK";
@cmd = ("Rscript $outputDir/utilities/plot_local_ancestry_pipeline.R ${id}_HGDP_RFMix $outputFile $inputFiles 40000");
makeJob($launchMethod, $tgt, $dep, @cmd);

######################
#1.7. log end time
######################
$tgt = "$outputDir/log/end.runmake.${id}_HGDP_RFMix.OK";
$dep = "${outputDir}/output/${id}_HGDP_local_ancestry_RFMix.png.OK";
@cmd = ("date | awk '{print \"\\nend: \"\$\$0}' >> $outputDir/log/runmake_${id}_HGDP_RFMix_time.log");
makeJob("local", $tgt, $dep, @cmd);


#############################################
#RFMix local ancestry inference using onlu europe, africa, native american
#############################################


# ######################
# #2.0. Europe, Africa, Native America subset of HGDP
# ######################
# for my $chr (1..22)
# {
#     $tgt = "$outputDir/HGDP_938/HGDP_938_chr${chr}_subset_filtered_phased.vcf.gz.OK";
#     $dep = "";
#     @cmd = ("bcftools view $outputDir/HGDP_938/HGDP_938_chr${chr}_filtered_phased.vcf.gz \\
#         --force-samples -S $outputDir/HGDP_938/HGDP_europe_africa_native_america.txt \\
#         --output-type z --output-file $outputDir/HGDP_938/HGDP_938_chr${chr}_subset_filtered_phased.vcf.gz; \\
#     bcftools index -t -f $outputDir/HGDP_938/HGDP_938_chr${chr}_subset_filtered_phased.vcf.gz");
#     makeJob($launchMethod, $tgt, $dep, @cmd);
# }

# ######################
# #2.1. log the start time
# ######################
# $tgt = "$outputDir/log/start.runmake.${id}_HGDP_subset_RFMix.OK";
# $dep = "";
# @cmd = ("date | awk '{print \"Local ancestry pipeline using HGDP subset \\n\\nstart: \"\$\$0}' > $outputDir/log/runmake_${id}_HGDP_subset_RFMix_time.log");
# makeJob("local", $tgt, $dep, @cmd);

# ######################
# #2.2. merge sample and HGDP subset
# ######################
# for my $chr (1..22)
# {
#     $tgt = "${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}_filtered_phased.vcf.gz.OK";
#     $dep = "$outputDir/common_site/topmed_chr${chr}_HGDP_common.txt.OK $outputDir/HGDP_938/HGDP_938_chr${chr}_subset_filtered_phased.vcf.gz.OK";
#     @cmd = ("bcftools merge -O z -o ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}_filtered_phased.vcf.gz -R $outputDir/common_site/topmed_chr${chr}_HGDP_common.txt \\
#         ${outputDir}/LAI/${id}/${id}_filtered_phased_chr${chr}.vcf.gz $outputDir/HGDP_938/HGDP_938_chr${chr}_subset_filtered_phased.vcf.gz && \\
#         bcftools index -t -f ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}_filtered_phased.vcf.gz");
#     makeJob($launchMethod, $tgt, $dep, @cmd);
# }

# ######################
# #2.3. make files needed for RFMix allele, location, classes
# ######################
# for my $chr (1..22)
# {
#     $tgt = "${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}_RFMix_files_prep.OK";
#     $dep = "${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}_filtered_phased.vcf.gz.OK";
#     @cmd = ("bcftools view ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}_filtered_phased.vcf.gz \\
#       | bcftools query -f '[%GT]\\n' - | sed 's/|//g' > ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}.alleles && \\
#     bcftools view ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}_filtered_phased.vcf.gz \\
#         | bcftools query -f '%POS\\n' - | Rscript $outputDir/utilities/generate_markerLocations_file.R stdin $outputDir/genetic_map_GRCh37/genetic_map_chr${chr}_combined_b37.txt ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}.locations && \\
#     bcftools query -l ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}_filtered_phased.vcf.gz | Rscript $outputDir/utilities/make_classes_file.R \\
#         stdin ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}.classes && \\
#     sed 's/4/2/g' ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}.classes | sed 's/5/3/g' > ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}.classes.tmp && \\
#         mv ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}.classes.tmp ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}.classes");
#     makeJob($launchMethod, $tgt, $dep, @cmd);
# }

# ######################
# #2.4. run RFMix
# ######################
# $inputFiles=""; #clean up
# $inputFilesOK=""; #clean up
# for my $chr (1..22)
# {
#     $inputFiles .= " ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}.0.Viterbi.txt";
#     $inputFiles .= " $outputDir/common_site/topmed_chr${chr}_HGDP_common.txt";
#     $inputFilesOK .= " ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}_RFMix_run.OK";
#     $tgt = "${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}_RFMix_run.OK";
#     $dep = "${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}_RFMix_files_prep.OK";
#     @cmd = ("cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && \\
#     python ./RunRFMix.py PopPhased ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}.alleles \\
#         ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}.classes \\
#         ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}.locations \\
#         -o ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr} --forward-backward --num-threads 1 && \\
#     sed 's/2/4/g' ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}.0.Viterbi.txt | sed 's/3/5/g' > ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}.0.Viterbi.txt.tmp && \\
#         mv ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}.0.Viterbi.txt.tmp ${outputDir}/LAI/${id}/${id}_HGDP_subset_chr${chr}.0.Viterbi.txt");
#     makeJob($launchMethod, $tgt, $dep, @cmd);
# }

# ######################
# #2.5. make local ancestry plot
# ######################
# $outputFile = "${outputDir}/output/${id}_HGDP_subset_local_ancestry_RFMix.png";
# $tgt = "${outputDir}/output/${id}_HGDP_subset_local_ancestry_RFMix.png.OK";
# $dep = "$inputFilesOK";
# @cmd = ("Rscript $outputDir/utilities/plot_local_ancestry_pipeline.R ${id}_HGDP_subset_RFMix $outputFile $inputFiles 40000");
# makeJob($launchMethod, $tgt, $dep, @cmd);

# ######################
# #2.6. log end time
# ######################
# $tgt = "$outputDir/log/end.runmake.${id}_HGDP_subset_RFMix.OK";
# $dep = "${outputDir}/output/${id}_HGDP_subset_local_ancestry_RFMix.png.OK";
# @cmd = ("date | awk '{print \"\\nend: \"\$\$0}' >> $outputDir/log/runmake_${id}_HGDP_subset_RFMix_time.log");
# makeJob("local", $tgt, $dep, @cmd);


# #############################################
# #local ancestry inference using LAMPLD
# #############################################


# ######################
# #3.0. create subset reference panel by ancestral populations
# ######################
# $tgt = "$outputDir/HGDP_938/LAMPLD/create_LAMP_ref_hap.OK";
# $dep = "";
# @cmd = ("./HGDP_938/create_LAMP_ref_hap.sh");
# makeJob("local", $tgt, $dep, @cmd);

# ######################
# #3.1. log the start time
# ######################
# $tgt = "$outputDir/log/start.runmake.${id}_HGDP_subset_LAMPLD.OK";
# $dep = "";
# @cmd = ("date | awk '{print \"Local ancestry pipeline using HGDP subset \\n\\nstart: \"\$\$0}' > $outputDir/log/runmake_${id}_HGDP_subset_LAMPLD_time.log");
# makeJob("local", $tgt, $dep, @cmd);

# ######################
# #3.2. convert sample vcf with common snp to HGDP to genotype dosage
# ######################
# for my $chr (1..22)
# {
#     $tgt = "${outputDir}/LAI/${id}/${id}_chr${chr}_lamp.012.OK";
#     $dep = "";
#     @cmd = ("vcftools --gzvcf ${outputDir}/LAI/${id}/${id}_filtered_phased_chr${chr}.vcf.gz --012 --out ${outputDir}/LAI/${id}/${id}_chr${chr}_lamp --positions $outputDir/common_site/topmed_chr${chr}_HGDP_common.txt && \\
#         rm -f ${outputDir}/LAI/${id}/${id}_chr${chr}_lamp.012.indv ${outputDir}/LAI/${id}/${id}_chr${chr}_lamp.012.pos ${outputDir}/LAI/${id}/${id}_chr${chr}_lamp.log && \\
#         cut -f 2- ${outputDir}/LAI/${id}/${id}_chr${chr}_lamp.012 | sed 's/\t//g' > ${outputDir}/LAI/${id}/${id}_chr${chr}_lamp.012.tmp && mv ${outputDir}/LAI/${id}/${id}_chr${chr}_lamp.012.tmp ${outputDir}/LAI/${id}/${id}_chr${chr}_lamp.012");
#     makeJob($launchMethod, $tgt, $dep, @cmd);
# }

# ######################
# #3.3. run LAMPLD
# ######################
# for my $chr (1..22)
# {
#     $tgt = "${outputDir}/LAI/${id}/${id}_chr${chr}_lampped.out.OK";
#     $dep = "$outputDir/HGDP_938/LAMPLD/create_LAMP_ref_hap.OK ${outputDir}/LAI/${id}/${id}_chr${chr}_lamp.012.OK";
#     @cmd = ("perl /net/snowwhite/home/khlin/tools/LAMPLD-v1.1/run_LAMPLD.pl $outputDir/HGDP_938/LAMPLD/chr${chr}.pos \\
#         $outputDir/HGDP_938/LAMPLD/HGDP_europe_chr${chr}.impute.hap \\
#         $outputDir/HGDP_938/LAMPLD/HGDP_native_america_chr${chr}.impute.hap \\
#         $outputDir/HGDP_938/LAMPLD/HGDP_africa_chr${chr}.impute.hap \\
#         ${outputDir}/LAI/${id}/${id}_chr${chr}_lamp.012 \\
#         ${outputDir}/LAI/${id}/${id}_chr${chr}_lampped.out");
#     makeJob($launchMethod, $tgt, $dep, @cmd);
# }

# ######################
# #3.4. convert plot format 
# ######################
# $inputFiles=""; #clean up
# $inputFilesOK=""; #clean up
# for my $chr (1..22)
# {   $inputFiles .= " ${outputDir}/LAI/${id}/${id}_chr${chr}_lampped_plot.out";
#     $inputFiles .= " $outputDir/common_site/topmed_chr${chr}_HGDP_common.txt";
#     $inputFilesOK .= " ${outputDir}/LAI/${id}/${id}_chr${chr}_lampped_plot.out.OK";
#     $tgt = "${outputDir}/LAI/${id}/${id}_chr${chr}_lampped_plot.out.OK";
#     $dep = "${outputDir}/LAI/${id}/${id}_chr${chr}_lampped.out.OK";
#     @cmd = ("python $outputDir/utilities/lamped_out_2_plot.py -i ${outputDir}/LAI/${id}/${id}_chr${chr}_lampped.out \\
#         -p $outputDir/common_site/topmed_chr${chr}_HGDP_common.txt -o ${outputDir}/LAI/${id}/${id}_chr${chr}_lampped_plot.out");
#     makeJob($launchMethod, $tgt, $dep, @cmd);
# }

# ######################
# #3.5. make local ancestry plot
# ######################
# $outputFile = "${outputDir}/output/${id}_HGDP_subset_local_ancestry_LAMPLD.png";
# $tgt = "${outputDir}/output/${id}_HGDP_subset_local_ancestry_LAMPLD.png.OK";
# $dep = "$inputFilesOK";
# @cmd = ("Rscript $outputDir/utilities/plot_local_ancestry_pipeline.R ${id}_HGDP_subset_LAMPLD $outputFile $inputFiles 40000");
# makeJob($launchMethod, $tgt, $dep, @cmd);

# ######################
# #3.6. log end time
# ######################
# $tgt = "$outputDir/log/end.runmake.${id}_HGDP_subset_LAMPLD.OK";
# $dep = "${outputDir}/output/${id}_HGDP_subset_local_ancestry_LAMPLD.png.OK";
# @cmd = ("date | awk '{print \"\\nend: \"\$\$0}' >> $outputDir/log/runmake_${id}_HGDP_subset_LAMPLD_time.log");
# makeJob("local", $tgt, $dep, @cmd);

#*******************
#Write out make file
#*******************
open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all: @tgts\n\n";

#clean
push(@tgts, "clean");
push(@deps, "");
push(@cmds, "\t-rm -rf *NA* log/*.* output/*.OK ${outputDir}/LAI/NA*/*.* $slurmScriptsDir/*.*");

#cleam_ref
push(@tgts, "clean_ref");
push(@deps, "");
push(@cmds, "\t-rm -rf $outputDir/HGDP_938/*.OK $outputDir/common_site/*.*");

#clean_LAMPLD
push(@tgts, "clean_lampld");
push(@deps, "");
push(@cmds, "\t-rm -rf ${outputDir}/LAI/NA*/*LAMPLD* ${outputDir}/LAI/NA*/*lampped* ${outputDir}/LAI/NA*/*lamp*");

#clean_lanccsv
push(@tgts, "clean_lanccsv");
push(@deps, "");
push(@cmds, "\t-rm -rf ${outputDir}/LAI/NA*/*LancCSV* ");

#clean_jobs
push(@tgts, "clean_job");
push(@deps, "");
push(@cmds, "\tps xu | awk '{if(\$\$11==\"make\") print \$\$2}' | xargs --verbose kill; squeue -ukhlin | awk '{if(\$\$3 ~ /NWD[0-9]/) print \$\$1}' | xargs scancel");

 
for(my $i=0; $i < @tgts; ++$i)
{
    print MAK "$tgts[$i]: $deps[$i]\n";
    print MAK "$cmds[$i]\n";
}
close MAK;
 
##########
#functions
##########
 
#run a job either locally or by slurm
sub makeJob
{
    my ($method, $tgt, $dep, @cmd) = @_;
 
    if ($method eq "local")
    {
        makeLocalStep($tgt, $dep, @cmd);
    }
    elsif ($method eq "slurm")
    {
        makeSlurm($tgt, $dep, @cmd);
    }
}
 
#run slurm jobs
sub makeSlurm
{
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        #contains pipe exactly one time or cd
        if ($c=~/\||^cd/)
        {
            ++$slurmScriptNo;
            my $slurmScriptFile = "$slurmScriptsDir/${slurmScriptNo}_${id}.sh";
            open(IN, ">$slurmScriptFile");
            print IN "#!/bin/bash\n"; 
            print IN "set pipefail; $c"; 
            close(IN);
            chmod(0755, $slurmScriptFile);
            
            # $cmd .= "echo '" . $c . "'\n";
            $cmd .= "\tsrun -p nomosix,main -J $id -D $outputDir $param $slurmScriptFile\n";
        }
        else
        {
            $cmd .= "\tsrun -p nomosix,main -J $id -D $outputDir $param " . $c . "\n";
        }
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}
 
#run a local job
sub makeLocalStep
{
    my ($tgt, $dep, @cmd) = @_;
 
    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        $cmd .= "\t" . $c . "\n";
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}
