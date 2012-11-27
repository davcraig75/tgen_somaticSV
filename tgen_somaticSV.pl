#!/usr/bin/perl -w
###########################################################################
#
#  Program:  tgen_somaticSV.pl
#  Author:   David W. Craig, Ph.D.
#  Created:  2012
#
#  This script is generic and for any specified database, it will
#  use Bio::Parser::Util::Serialize to parse the text source file into
#  individual record text and object files as requested.
#
#  $Id: tgen_somaticSV.pl,v 0.13 dcraig Exp $
#
###########################################################################

use POSIX;
$|=1;
#######################
# Variables & Main Loop
#######################
MAIN: {

  $mqf=10;
  $flags= "-F 0x0400 -F 0x0200 -q $mqf"; 
  $bin=500;
  $gap_min=$bin;
  $min_hits=3;
  $i_min=50; $i_max=1000;
  &command_arg;
  $skip=1000000;
  ($i_min,$i_max)=calc_insert($tumor_file[0],$i_min,$i_max);
  $i_min=50; $i_max=1500;
  ($i_min,$i_max)=calc_insert($ref_file[0],$i_min,$i_max);
  my $tcnt={};
  my $rcnt={};
  my $tref_count={};
  my $ref_count={};
  for ($i=0;$i<=$#tumor_file;++$i) {
    print "-opening $tumor_file[$i]\n";
    ($tcnt,$tref_count)=main_loop($tumor_file[$i],$tcnt,$tref_count);
  }
  for ($i=0;$i<=$#ref_file;++$i) {
    print "-opening $ref_file[$i]\n";
    ($rcnt,$ref_count)=main_loop($ref_file[$i],$rcnt,$ref_count);
  }
  &print_out($tcnt,$rcnt,$ref_count);
} 
#################
# Main Loop
#################
sub main_loop {
  @time = localtime(time);
  my $in=$_[0];
  my %cnt=%{$_[1]};
  my %ref_cnt=%{$_[2]};
  my $c=0;
  my $last_pos=0;
  open (IN,"$s view $flags $in -L $lo |") || die "Can't open $in\n";
  LOOP: while (<IN>) {
    my @temp=split(/\t/);
    ++$c;
    if ($c % $skip == 0 && $c>1) {
      @prev = @time;
      @time = localtime(time);
      $sec=60*($time[1]-$prev[1])+$time[0]-$prev[0];
      print "$in Read: $c at $temp[2]\t$temp[3] -  $skip reads in $sec seconds\n";
    }
    if ($temp[6] eq "*") {next LOOP;}
    #if (/NH:i:([0-9])/) {if($1>3) {next LOOP}}
    #if (/CM:i:([0-9])/){if ($1>3) {next LOOP}}
    if ($temp[1] & 0x10 ) {$d="-"} else {$d="+";} #direction
    if ($temp[1] & 0x20 ) {$p="-"} else {$p="+";} #direction
    if ($temp[6] eq "=") {$temp[6]=$temp[2]}
    if (abs($temp[8]) > $i_min && abs($temp[8]) < $i_max) {
      $c_pos=$bin*int($temp[3]/$bin);
      if (exists($ref_count{"$temp[2]:$c_pos"})) {
        ++$ref_count{"$temp[2]:$c_pos"};
      } else {
        $ref_count{"$temp[2]:$c_pos"}=1;
      }
      next LOOP
    }
    if ($temp[3] == $last_pos) {next LOOP} else {$last_pos=$temp[3]}
    if ($_=~/MQ:i:(\d*)\t/) {$omq=$1; if ($omq<10) {next LOOP}}
    $do=$bin*int($temp[3]/$bin);
    $ac=$bin*int($temp[7]/$bin);
    $an="$p:$temp[6]:$ac|$d:$temp[2]:$do";
    $ac_be=$ac-1*$bin;
    $ac_af=$ac+1*$bin;
    $do_be=$do-1*$bin;
    $do_af=$do+1*$bin;
    $t0="$p:$temp[6]:$ac|$d:$temp[2]:$do";
    $t1="$p:$temp[6]:$ac_be|$d:$temp[2]:$do";
    $t2="$p:$temp[6]:$ac_af|$d:$temp[2]:$do";
    $t3="$p:$temp[6]:$ac_be|$d:$temp[2]:$do_be";
    $t4="$p:$temp[6]:$ac_af|$d:$temp[2]:$do_be";
    $t5="$p:$temp[6]:$ac_be|$d:$temp[2]:$do_af";
    $t6="$p:$temp[6]:$ac_af|$d:$temp[2]:$do_af";
    $t7="$p:$temp[6]:$ac|$d:$temp[2]:$do_be";
    $t8="$p:$temp[6]:$ac|$d:$temp[2]:$do_af";
    $max=0;
    if(exists($cnt{$t0})){ if ($cnt{$t0} > $max ) {$an=$t0;$max=$cnt{$t0}}}
    if(exists($cnt{$t1})){ if ($cnt{$t1} > $max ) {$an=$t1;$max=$cnt{$t1}}}
    if(exists($cnt{$t2})){ if ($cnt{$t2} > $max ) {$an=$t2;$max=$cnt{$t2}}}
    if(exists($cnt{$t3})){ if ($cnt{$t3} > $max ) {$an=$t3;$max=$cnt{$t3}}}
    if(exists($cnt{$t4})){ if ($cnt{$t4} > $max ) {$an=$t4;$max=$cnt{$t4}}}
    if(exists($cnt{$t5})){ if ($cnt{$t5} > $max ) {$an=$t5;$max=$cnt{$t5}}}
    if(exists($cnt{$t6})){ if ($cnt{$t6} > $max ) {$an=$t6;$max=$cnt{$t6}}}
    if(exists($cnt{$t7})){ if ($cnt{$t7} > $max ) {$an=$t6;$max=$cnt{$t8}}}
    if(exists($cnt{$t8})){ if ($cnt{$t8} > $max ) {$an=$t8;$max=$cnt{$t8}}}
    if (exists($cnt{$an})){ ++$cnt{$an}; } else { $cnt{$an}=1; }
    if ($c % 1000000 == 0) {print "Read: $c at $temp[2]\t$temp[3]\n";}
  }
  close (IN);
  return (\%cnt,\%ref_count);
} 
#################
# Print Out
#################
sub print_out {
  my $tlcnt=$_[0];
  my $rlcnt=$_[1];
  my $href_count=$_[2];
  $href_mean=href_mean_calc($href_count);
  print "-Printing output\n";
  open (OUT,">$out.dat") || die "Can't open $out.dat\n";
  print OUT join("\t","SV_desc","Str\tChr\tPos\tStr\tChr\tPos","Tumor_count","Ref_count","Tumor_FF","Tumor_rr","Tumor_fr","Tumor_rr","rev_Tumor_sff","rev_Tumor_rr","rev_Tumor_fr","rev_Tumor_rr","Distance","Ref_NoSV_pairs_don","Ref_NoSV_pairs_before_don","Ref_NoSV_pairs_after_don","Ref_NoSV_pairs_acc","Ref_NoSV_pairs_before_acc","Ref_NoSV_pairs_after_acc") . "\n";
  LOOP: foreach $key (sort {$tlcnt->{$b} <=> $tlcnt-> {$a}} keys %$tlcnt ) {
    if (!(exists($rlcnt->{$key}))){$rlcnt->{$key}=0}
    $min_con=floor(0.00*$tlcnt->{$key});
    if ($tlcnt->{$key}>$min_hits && $rlcnt->{$key}<=$min_con && $tlcnt->{$key} < 50000) {
       ($acc,$don)=split(/\|/,$key);
       ($dira,$chra,$posa)=split(/\:/,$acc);
       ($dird,$chrd,$posd)=split(/\:/,$don);
       $ff="+:$chra:$posa|+:$chrd:$posd";
       $fr="+:$chra:$posa|-:$chrd:$posd";
       $rr="-:$chra:$posa|-:$chrd:$posd";
       $rf="-:$chra:$posa|+:$chrd:$posd";

       $sff="+:$chrd:$posd|+:$chra:$posa";
       $sfr="+:$chrd:$posd|-:$chra:$posa";
       $srr="-:$chrd:$posd|-:$chra:$posa";
       $srf="-:$chrd:$posd|+:$chra:$posa";

      $tum_count_ff=0;$tum_count_rr=0;$ref_count_rr=0;$ref_count_ff=0;
      $tum_count_sff=0;$tum_count_srr=0;$ref_count_srr=0;$ref_count_sff=0;
      $tum_count_rf=0;$tum_count_rf=0;$ref_count_fr=0;$ref_count_rf=0;
      $tum_count_srf=0;$tum_count_srf=0;$ref_count_sfr=0;$ref_count_srf=0;
      if (exists($tlcnt->{$ff})) {$tum_count_ff=$tlcnt->{$ff};} else {$tum_count_ff=0}
      if (exists($tlcnt->{$fr})) {$tum_count_fr=$tlcnt->{$fr};} else {$tum_count_fr=0}
      if (exists($tlcnt->{$rr})) {$tum_count_rr=$tlcnt->{$rr};} else {$tum_count_rr=0}
      if (exists($tlcnt->{$rf})) {$tum_count_rf=$tlcnt->{$rf};} else {$tum_count_rf=0}
      if (exists($tlcnt->{$sff})) {$tum_count_sff=$tlcnt->{$sff};} else {$tum_count_sff=0}
      if (exists($tlcnt->{$sfr})) {$tum_count_sfr=$tlcnt->{$sfr};} else {$tum_count_sfr=0}
      if (exists($tlcnt->{$srr})) {$tum_count_srr=$tlcnt->{$srr};} else {$tum_count_srr=0}
      if (exists($tlcnt->{$srf})) {$tum_count_srf=$tlcnt->{$srf};} else {$tum_count_srf=0}

      if (exists($rlcnt->{$ff})) {$ref_count_ff=$rlcnt->{$ff};} else {$ref_count_ff=0}
      if (exists($rlcnt->{$fr})) {$ref_count_fr=$rlcnt->{$fr};} else {$ref_count_fr=0}
      if (exists($rlcnt->{$rf})) {$ref_count_rf=$rlcnt->{$rf};} else {$ref_count_rf=0}
      if (exists($rlcnt->{$rr})) {$ref_count_rr=$rlcnt->{$rr};} else {$ref_count_rr=0}
      if (exists($rlcnt->{$sff})) {$ref_count_ff=$rlcnt->{$sff};} else {$ref_count_sff=0}
      if (exists($rlcnt->{$sfr})) {$ref_count_fr=$rlcnt->{$sfr};} else {$ref_count_sfr=0}
      if (exists($rlcnt->{$srf})) {$ref_count_rf=$rlcnt->{$srf};} else {$ref_count_srf=0}
      if (exists($rlcnt->{$srr})) {$ref_count_rr=$rlcnt->{$srr};} else {$ref_count_srr=0}

      $tum_count=$tum_count_ff+$tum_count_rr+$tum_count_fr+$tum_count_rf+$tum_count_sff+$tum_count_srr+$tum_count_sfr+$tum_count_srf;
      $ref_count=$ref_count_ff+$ref_count_rr+$ref_count_fr+$ref_count_rf+$ref_count_sff+$ref_count_srr+$ref_count_sfr+$ref_count_srf;
      if ($chrd eq $chra) {$gap=abs($posa-$posd)} else {$gap = 1000000000;}
      if ($gap>$gap_min) {
        $dons=$posd-1*$bin;
        $done=$posd+1*$bin;
        $accs=$posa-1*$bin;
        $acce=$posa+1*$bin;
        $n=$key;
        $n=~s/\:/\t/g;
        $n=~s/\|/\t/g;
        $r="+:$chra:$accs|+:$chrd:$dons";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="+:$chra:$acce|+:$chrd:$dons";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r} }
        $r="+:$chra:$posa|+:$chrd:$dons";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="+:$chra:$posa|+:$chrd:$done";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="+:$chra:$acce|+:$chrd:$done";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="+:$chra:$accs|+:$chrd:$done";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="+:$chra:$acce|+:$chrd:$posd";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="+:$chra:$accs|+:$chrd:$posd";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}

        $r="-:$chra:$accs|+:$chrd:$dons";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="-:$chra:$acce|+:$chrd:$dons";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="-:$chra:$posa|+:$chrd:$dons";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="-:$chra:$posa|+:$chrd:$done";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="-:$chra:$acce|+:$chrd:$done";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="-:$chra:$accs|+:$chrd:$done";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="-:$chra:$acce|+:$chrd:$posd";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="-:$chra:$accs|+:$chrd:$posd";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}

        $r="+:$chra:$accs|-:$chrd:$dons";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="+:$chra:$acce|-:$chrd:$dons";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="+:$chra:$posa|-:$chrd:$dons";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="+:$chra:$posa|-:$chrd:$done";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="+:$chra:$acce|-:$chrd:$done";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="+:$chra:$accs|-:$chrd:$done";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="+:$chra:$acce|-:$chrd:$posd";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="+:$chra:$accs|-:$chrd:$posd";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="-:$chra:$accs|-:$chrd:$dons";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="-:$chra:$acce|-:$chrd:$dons";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="-:$chra:$posa|-:$chrd:$dons";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="-:$chra:$posa|-:$chrd:$done";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="-:$chra:$accs|-:$chrd:$done";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="-:$chra:$acce|-:$chrd:$done";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="-:$chra:$acce|-:$chrd:$posd";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $r="-:$chra:$accs|-:$chrd:$posd";if (exists($rlcnt->{$r})) {$ref_count=$ref_count+ $rlcnt->{$r}}
        $acnorm_pos="$chra:$posa";
        $acnorm_bef="$chra:$accs";
        $acnorm_aff="$chra:$acce";
        $donorm_pos="$chrd:$posd";
        $donorm_aff="$chrd:$done";
        $donorm_bef="$chrd:$dons";
        if (!(exists($href_count->{$acnorm_pos}))) {$acnorm_pos_m=0;}else{$acnorm_pos_m=sprintf("%1.1f",$href_count->{$acnorm_pos}/$href_mean)}
        if (!(exists($href_count->{$acnorm_bef}))) {$acnorm_bef_m=0;}else{$acnorm_bef_m=sprintf("%1.1f",$href_count->{$acnorm_bef}/$href_mean)}
        if (!(exists($href_count->{$acnorm_aff}))) {$acnorm_aff_m=0;}else{$acnorm_aff_m=sprintf("%1.1f",$href_count->{$acnorm_aff}/$href_mean)}
        if (!(exists($href_count->{$donorm_pos}))) {$donorm_pos_m=0;}else{$donorm_pos_m=sprintf("%1.1f",$href_count->{$donorm_pos}/$href_mean)}
        if (!(exists($href_count->{$donorm_aff}))) {$donorm_aff_m=0;}else{$donorm_aff_m=sprintf("%1.1f",$href_count->{$donorm_aff}/$href_mean)}
        if (!(exists($href_count->{$donorm_bef}))) {$donorm_bef_m=0;}else{$donorm_bef_m=sprintf("%1.1f",$href_count->{$donorm_bef}/$href_mean)}
        if ($ref_count <= $min_con && !(exists($isdone{$key})) ) {
          if ($acnorm_pos_m < 3 && $acnorm_bef_m < 3 && $acnorm_aff_m < 3 && $donorm_pos_m<3 && $donorm_aff_m < 3 && $donorm_bef_m < 3) {
            print OUT join("\t",$key,$n,$tum_count,$ref_count,$tum_count_ff,$tum_count_rr,$tum_count_fr,$tum_count_rf,$tum_count_sff,$tum_count_srr,$tum_count_sfr,$tum_count_srf,$gap,$acnorm_pos_m,$acnorm_bef_m,$acnorm_aff_m,$donorm_pos_m,$donorm_aff_m,$donorm_bef_m,$href_mean) . "\n";
            $isdone{$ff}=1;$isdone{$rr}=1;$isdone{$rf}=1;$isdone{$fr}=1;
            $isdone{$sff}=1;$isdone{$srr}=1;$isdone{$srf}=1;$isdone{$sfr}=1;
          }
        }
      }
    }
  }
  close (OUT);
}

sub href_mean_calc {
  my $href_count=$_[0];
  my $t=0;my $c=0;
  LOOP3: foreach $key (keys %$href_count) {
    $t=$t+$href_count->{$key}; 
    ++$c;
    if ($c>100000) { last LOOP3;} 
  }
  return int($t/$c);
}

sub calc_insert {
 print "-Calculating insert size for $_[0]\n";
 my $in=$_[0];
 open (IN,"$s view -s 0.000001 $in |") || die "Can't open $in\n";
 my @insert_vals=();
 my $i_min=$_[1];
 my $i_max=$_[2];
 LOOP:while (<IN>) {
    my @temp=split(/\t/);
    if ($temp[6] eq "*") {next LOOP;}
    if (abs($temp[8]) > $i_min && abs($temp[8]) < $i_max) {
      push (@insert_vals,abs($temp[8]));
    } #else { print "temp8 $temp[8] - $temp[6]\n"}
    if ($#insert_vals >10000) {last LOOP}
 }
 @sort_insert_vals=sort { $a <=> $b } @insert_vals;
 $i_min=$sort_insert_vals[10];
 $i_max=$sort_insert_vals[9990];
 $i_med=$sort_insert_vals[5000];
 print "-$in max_insert=$i_max\tmin_insert=$i_min\tmed_insert=$i_med\n";
 return ($i_min,$i_max);
}


#################
# Initialize
#################

sub command_arg {
  $out="out.cln";
  $s="samtools";
  $lo="";
  $nz=0;
  if ($#ARGV < 1) {
    print "./tran.9.pl R=mca005.Normal.bwa.i.d.r.bam T=mca005.Tumor.bwa.i.d.r.bam I=1000 O=mca005 M=RG: > MCA.log\n"; 
    die;
  } else {
    foreach $a (@ARGV) {
      if ($a=~/T=/) {
        $a=~s/T=//g;
        push(@tumor_file,$a);
        print "-Tumor_file:\t$tumor_file[$#tumor_file]\n";
       } elsif ($a=~/R=/) {
         $a=~s/R=//g;
         push(@ref_file,$a);
         print "-Reference:\t$ref_file[$#ref_file]\n";
       } elsif ($a=~/L=/) {
        $lo=$a;
        $lo=~s/L=//g;
        print "-Location:\t$lo\n";
       } elsif ($a=~/O=/) {
        $out="$a.trn";
        $out=~s/O=//g;
        print "-Out_head:\t$out\n";
        print "-Max Insert: $i_max\n";
       } elsif ($a eq "--nz") {
        $nz=1;
        print "-Only outputing non-zero coverage\n";
       } elsif ($a=~/B=/) {
        $bin=$a;
        $bin=~s/B=//g;
        print "-Bin: $bin\n";
       }elsif ($a=~/S=/) {
        $s=$a; 
        $s=~s/S=//;
        print "S=$s\n";
       }else {die "DIE: Can't parse $a\n";}
    }
  } 
}
