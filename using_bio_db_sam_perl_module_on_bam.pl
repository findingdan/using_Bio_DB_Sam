#!/usr/bin/env perl

# Copyright (c) 2007 Genome Research Ltd.
# Author: Dan King <dk6@sanger.ac.uk>
#
# Any redistribution or derivation in whole or in part including any
# substantial portion of this code must include this copyright and
# permission notice.
#
# THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# This code is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at
# your option) any later version (http://www.gnu.org/copyleft/gpl.txt).

use warnings;
use strict;
use Data::Dumper; 
use Bio::DB::Sam;

Bio::DB::Sam->max_pileup_cnt(250);
my $bam_obj = Bio::DB::Sam->new( -bam => "example.bam");

open (POS, "example.positions") or die $!;
while (<POS>) {
	chomp;
	my ($pos, $ref, $alt) = split;
	my $baf = retrieve_pileup_counts($bam_obj, "chr17", $pos, $ref, $alt);
	print join ("\t", $pos, $ref, $alt, $baf),"\n";
}

sub retrieve_pileup_counts {
    my ($bam_obj, $needed_chr, $needed_pos, $ref, $alt) = @_;
    my %count;
    my $cb = sub {
        my ($seqid, $pos, $pileup_aref, $sam) = @_;
        if ($pos == $needed_pos) {
            for my $pileup (@{$pileup_aref}){
                my $al = $pileup->alignment;
				
				next unless $al->proper_pair == 1; 				# Returns true if mate and pair are both mapped 

                my $baseQual = $al->qscore->[$pileup->qpos];    # Phred qual of base (in numeric format)
                next unless $baseQual >= 0;

				my $alignment_mapping_Qual = $al->qual;
				next unless $alignment_mapping_Qual >= 0;  	# "white" or "yellow" in tview

				(my $cigar_str_letters = $al->cigar_str) =~ s/\d//g;
				next unless $cigar_str_letters =~ m/^M$/;		# Only use alignments with cigar strings of only matches (no soft clipped bases or indels)

                my $qBase = substr($al->qseq, $pileup->qpos, 1);# The base at the needed position
                $count{$qBase}++;
            }
        }
    };
    $bam_obj->fast_pileup("$needed_chr:$needed_pos-$needed_pos", $cb);
	my $altct = $count{$alt} || 0;
	my $refct = $count{$ref} || 0;
	return if ($altct + $refct == 0);
	my $baf = $altct / ($altct + $refct);
	return $baf;
}

