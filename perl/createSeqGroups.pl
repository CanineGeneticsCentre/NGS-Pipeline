#! /usr/bin/env perl

use strict;
use warnings;

my $dict = $ARGV[0];
my %seqs = ();
my $longest_seq = 0;

open (DICT, $dict) or die "Unable to find fasta dictionary - $dict\n";
while (<DICT>){
  chomp $_;
  my @tmp = split("\t", $_);
  if($_ =~ /^\@SQ/){
    $tmp[1] =~ s/SN://;
    $tmp[2] =~ s/LN://;
    if ($tmp[2] > $longest_seq){
      $longest_seq = $tmp[2];
    }
    # print $tmp[1]."\t".$tmp[2];
    $seqs{$tmp[1]} = $tmp[2];
  }
}

my $tsv_string = '';
my $temp_size = 0;
foreach my $chr (keys %seqs){
  if ($temp_size == 0){
    $tsv_string = $chr;
    $temp_size = $seqs{$chr};
    next;
  }

  if (($temp_size + $seqs{$chr}) <= $longest_seq){
    $temp_size += $seqs{$chr};
    $tsv_string .= "\t".$chr;
  } else {
    $tsv_string .= "\n".$chr;
    $temp_size = $seqs{$chr};
  }
}

open(OUT, ">sequence_grouping.txt" );
print OUT $tsv_string."\n";
close OUT;

$tsv_string .= "\nunmapped";

open(OUT, ">sequence_grouping_with_unmapped.txt" );
print OUT $tsv_string."\n";
close OUT;
