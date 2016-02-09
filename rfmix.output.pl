#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

extract_features.pl

=head1 SYNOPSIS

 extract_features.pl [options]

  -s     sample file list giving the location of each sample
         column 1: sample name
         column 2: path of bam file
  -r     reference sequence fasta file
  -b     binaries directory : location of binaries required for this pipeline
  -o     output directory : location of all output files
  -m     output make file

 example:

=head1 DESCRIPTION

=cut

#option variables
my $help;

#
my $outputDir = "/net/fantasia/home/atks/dev/20150812_spectrum";
my $vt = "/net/fantasia/home/atks/dev/vt/vt";
my $clusterDir = "/net/fantasia/home/atks/programs/cluster";
my $makeFile = "spectrum_analysis.mk";
my $intervals = "20";
my $refGenomeFASTAFile = "/net/1000g/atks/ref/genome/hs37d5.fa";
my $launchMethod = "local";

#programs

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help)
  || scalar(@ARGV)!=1)
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

my $inputVCFFile = $ARGV[0];

open(VCF, "bcftools view $inputVCFFile |");

open(VCF1, ">out1.txt");
open(VCF2, ">out2.txt");


while (<VCF>)
{
    chomp;
    
    my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, @gt) = split("\t");
    
    print VCF1 "$pos\t$gt[0]\n";
    print VCF2 "$pos\t$gt[0]\n";
    
    
}

close(VCF1);
close(VCF2);
