# /**
#
# Author:      Linzheng
# DateTime:    2016-12-08 17:35:35
# Description: Description
#
#**/

use strict;
use warnings;
use 5.010;
sub code2aa {
	my ($code) = @_ ;
	my(%genetic_code) = (
    'TCA' => 'S',    # Serine 
    'TCC' => 'S',    # Serine 
    'TCG' => 'S',    # Serine 
    'TCT' => 'S',    # Serine 
    'TTC' => 'F',    # Phenylalanine 
    'TTT' => 'F',    # Phenylalanine 
    'TTA' => 'L',    # Leucine 
    'TTG' => 'L',    # Leucine 
    'TAC' => 'Y',    # Tyrosine 
    'TAT' => 'Y',    # Tyrosine 
    'TAA' => '_',    # Stop 
    'TAG' => '_',    # Stop 
    'TGC' => 'C',    # Cysteine 
    'TGT' => 'C',    # Cysteine 
    'TGA' => '_',    # Stop 
    'TGG' => 'W',    # Tryptophan 
    'CTA' => 'L',    # Leucine 
    'CTC' => 'L',    # Leucine 
    'CTG' => 'L',    # Leucine 
    'CTT' => 'L',    # Leucine 
    'CCA' => 'P',    # Proline 
    'CCC' => 'P',    # Proline 
    'CCG' => 'P',    # Proline 
    'CCT' => 'P',    # Proline 
    'CAC' => 'H',    # Histidine 
    'CAT' => 'H',    # Histidine 
    'CAA' => 'Q',    # Glutamine 
    'CAG' => 'Q',    # Glutamine 
    'CGA' => 'R',    # Arginine 
    'CGC' => 'R',    # Arginine 
    'CGG' => 'R',    # Arginine 
    'CGT' => 'R',    # Arginine 
    'ATA' => 'I',    # Isoleucine 
    'ATC' => 'I',    # Isoleucine 
    'ATT' => 'I',    # Isoleucine 
    'ATG' => 'M',    # Methionine 
    'ACA' => 'T',    # Threonine 
    'ACC' => 'T',    # Threonine 
    'ACG' => 'T',    # Threonine 
    'ACT' => 'T',    # Threonine 
    'AAC' => 'N',    # Asparagine 
    'AAT' => 'N',    # Asparagine 
    'AAA' => 'K',    # Lysine 
    'AAG' => 'K',    # Lysine 
    'AGC' => 'S',    # Serine 
    'AGT' => 'S',    # Serine 
    'AGA' => 'R',    # Arginine 
    'AGG' => 'R',    # Arginine 
    'GTA' => 'V',    # Valine 
    'GTC' => 'V',    # Valine 
    'GTG' => 'V',    # Valine 
    'GTT' => 'V',    # Valine 
    'GCA' => 'A',    # Alanine 
    'GCC' => 'A',    # Alanine 
    'GCG' => 'A',    # Alanine 
    'GCT' => 'A',    # Alanine 
    'GAC' => 'D',    # Aspartic Acid 
    'GAT' => 'D',    # Aspartic Acid 
    'GAA' => 'E',    # Glutamic Acid 
    'GAG' => 'E',    # Glutamic Acid 
    'GGA' => 'G',    # Glycine 
    'GGC' => 'G',    # Glycine 
    'GGG' => 'G',    # Glycine 
    'GGT' => 'G',    # Glycine 
    'UGA' => 'U'     # Sec
    ); 
 
    if(exists $genetic_code{$code}) { 
        return $genetic_code{$code}; 
    }else{ 
 
            print STDERR "Bad code \"$code\"!!\n"; 
            exit; 
    } 

}


my $dna = <STDIN>;
chomp($dna);
my $protein = ''; 
my $codon; 
 
# Translate each three-base codon into an amino acid, and append to a protein  
for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
	    $codon = substr($dna,$i,3); 
    $protein .= code2aa($codon); 
} 
 
print "Translate the DNA\n\n$dna\n\n  into the protein\n\n$protein\n\n";