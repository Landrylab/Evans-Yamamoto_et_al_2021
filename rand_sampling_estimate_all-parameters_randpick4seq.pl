#!/usr/bin/perl -w

use strict;
use constant PI    => 4 * atan2(1, 1);
use constant DEBUG => 0;
use Data::Dumper;

##############################################
# Monte-Carlo simulation of BFG-Y2H screen   #
# Last modified by Nozomu Yachie on 03032014 #
##############################################


# PARAMETER DEFINITIONS ##################################################################################
#
# 01. CV (coefficient of variation) of strain abundances (log-normal distribution) in haploid pools is 30%
# 02. Yeast mating efficiency of Y2H is 1% (underestimation)
# 03. CV of X-Y pair-dependent mating efficiencies (log-normal distribution) is 50%
# 04. CV of X-Y pair-dependent growth amplitudes (log-normal distribution) in liquid media is 50%
# 05. Y2H positive rate is 0.1% (overestimation)
# 06. CV of X-dependent growth amplitudes (DB-X auto-activity, log-normal distribution) in Y2H-selectable condition is 100%
# 07. CV of X-Y pair-dependent growth amplitudes (log-normal distribution) in Y2H-selectable condition is 10000%
# 08. CV of X-Y pair-dependent growth amplitudes (log-normal distribution) in non-selectable condition is 10%
# 09. CV of X-Y pair-dependent PCR amplitudes (log-normal distribution) is 50%
# 10. Read depth of Illumina sequencing in each screen is 10,000,000
#
# 11. Yeast DNA miniprep yield from 3 × 10^7 diploid yeast cells is 30 ng
# 12. Fraction of Y2H plasmids in yeast DNA miniprep product is 6% of the total DNA mass
# 13. Barcode fusion efficiency is 20%
#
# 14. Y2H plasmid sizes are 10 kbp
# 15. 1 mole of double-stranded DNA molecule is 660 g bp^-1
# 16. 1 OD(600 nm) unit for haploid yeast is 3 × 10^7 cells ml^-1
# 17. 1 OD(600 nm) unit for diploid yeast is 1 × 10^7 cells ml^-1
#
##########################################################################################################

my $x_strain_num = shift ;
my $y_strain_num = shift ;
my $header       = 'simulation';

 
my $P1           = shift; # CV (coefficient of variation) of strain abundances (log-normal distribution) in haploid pools is 30%
my $P3           = shift; # CV of X-Y pair-dependent mating efficiencies (log-normal distribution) is 50%
my $P4            = shift; # Yeast mating efficiency of Y2H is 1% (underestimation)
my $P5           = shift; # CV of X-Y pair-dependent growth amplitudes (log-normal distribution) in liquid media is 50%
my $P_plated     = shift; # Amount of cells spread on each condition.
my $P_CFU        = shift; #
my $P10          = shift; # CV of X-Y pair-dependent growth amplitudes (log-normal distribution) in non-selectable condition is 10%
my $P_template   = shift; # Multiplification of the following ;  1.4 * 10**7 as experimental value
                         #  Yeast DNA miniprep yield from 3 × 10^7 diploid yeast cells is 20 ng
                         #  Fraction of Y2H plasmids in yeast DNA miniprep product is 6% of the total DNA mass
                         #  Barcode fusion efficiency is 20%
                         #   4 PCR replicates
my $P16          = shift; # CV of X-Y pair-dependent PCR amplitudes (log-normal distribution) is 50%









###############################################################
#print Dumper "Pooling haploid";

my $x_hap = &haploid($x_strain_num);
my $y_hap = &haploid($y_strain_num);

#for (my $hap_cv = 0; $hap_cv < 0.51; $hap_cv+=0.01){

# CV (coefficient of variation) of strain abundances (log-normal distribution) in haploid pools is 30%
my $hap_cv = sprintf("%.3f", $P1);

#print Dumper $hap_cv1;
#print Dumper $hap_cv;


my $hap_size1 = 1.8*10**10; # 3 × 10^7 /OD x 600 OD
$x_hap = &adjust_num($x_hap,$hap_cv,$hap_size1);
$y_hap = &adjust_num($y_hap,$hap_cv,$hap_size1);







###############################################################
#print Dumper "Yeast mating and query to diploid selection";

my $dip = &mating($x_hap,$y_hap);

my $mating    = $P4;                 # Mating efficiency
my $hap_size2 = 0.75*1.5/2*3 * 10**10;# Amount of haploid/diploid cells to be queried to diploid selection (750mL of OD0.5 = 1500* 0.5 * 3 * 10 ** 7)
my $dip_size1 = $hap_size2 * $mating; # Amount of diploid to be queried to diploid selection

# CV of X-Y pair-dependent mating efficiencies (log-normal distribution) is 50%
$dip = &adjust_num($dip,$P3,$dip_size1); #


###############################################################
#print Dumper "Diploid selection and amount of cells spread per plate";

my $dip_size2 = 5 * 10**8; # Amount of cells to be spread to each selection plate (10 plates with 5 OD unit of cells each = 5 *10 * 10 **7  )

# CV of X-Y pair-dependent growth amplitudes (log-normal distribution) in liquid media is 50%
$dip = &adjust_num($dip,$P5,$dip_size2);




###############################################################
#print Dumper  "Y2H selections";

my $CFU_control = $P_CFU; #Observation = 2 * 10**7;
my $dip_neg  = &adjust_num($dip,$P10,$CFU_control);
#my $pos_rate = 0.001;















###############################################################
#print Dumper "Dox induction and plasmid extraction";

my $plasmids = $P_template; # Number of plasmids before PCR per plated condition per replicate. 

# Input for PCR is 100 pg/µL x 16 µL x 2 reactions = 3.2 ng
# 6 % of the total DNA mass is Y2H products 

# Y2H plasmid sizes are 10 kbp
# 1 mole of double-stranded DNA molecule is 660 g bp^-1

# Molecule number subjected to PCR
# = 3200 pg x 0.06 x 6.02 x 10^23 /(660 x 10^12 pg  x 10000)
# = 192  pg        x 6.02 x 10^7 / 660 pg  x 0.2 (Fused product)
# = 0.35 x 10^7
# = 3.5 * 10**6 x 8 #For each of these, we perfomred 4 PCRs
# = 1.4 * 10 **7



#Final sampling size is PCR input. The growth amplitude is from Dox induction. 
#Do we assume the sampling to the dox induction is representing the pool? We took 2.5 * 10**8 cells from the scraped sample for dox induction. 
my $plasmid_neg = &adjust_num($dip_neg,$P5,$plasmids);  
my $library_input = 1.0* (6.02*10**7) ; #Under estimating 
my $bfg_pcr     = &adjust_num($plasmid_neg,$P16,$library_input); #Illumina input is 2nMx5µL = 10fmol for entire library. 3 condition x 2 rep x 2 barcode fuion x another assay (2) = 24. 10/24 = 0.41 fmol. Here we underestimate each library is 0.1fmol. 

#my $plasmid_pos = &adjust_num($dip_pos,0.5,$plasmids);







###############################################################
#print Dumper "PCR re-amplification and Illumina sequencing";



my $cr = 1.5;
for (my $reads = 10000; $reads < 50000000; $reads*=$cr){#95000000; $reads*=$cr){
#  print "count = $count ¥n";
    my $rounded_reads = int($reads+ 0.5);

#my $reads = 6.0 * 10**4;
    my $read_neg = &random_selection($bfg_pcr,$rounded_reads); #X-Y dependent PCR amplification is CV = 50%
#my $read_pos = &adjust_num($plasmid_pos,0.5,$reads);


#my ($marginal_x_hap,$marginal_y_hap) = &marginal($read_neg);



#  print 'BEGIN: Illumina sequencing'."\n";
#  print 'DESCRIPTION: simulating Illumina library prep and sequencing outputamount'."\n";

  {
#    print 'BEGIN: X-Y coverage in non-selectable condition'."\n";
#    print "\#reads\tcoverage\tTotal_reads\thap_cv\n";
    my $table = &coverage2($read_neg,0,$x_strain_num,$y_strain_num);
    
    &print_table_parameters($table,$rounded_reads,$P1,$P3,$P4,$P5,$P_plated,$P_CFU,$P10,$P_template,$P16);

#my P1           = shift; # CV (coefficient of variation) of strain abundances (log-normal distribution) in haploid pools is 30%
#my P3           = shift; # CV of X-Y pair-dependent mating efficiencies (log-normal distribution) is 50%
#my P5           = shift; # CV of X-Y pair-dependent growth amplitudes (log-normal distribution) in liquid media is 50%
#my P_plated     = shift; # Amount of cells spread on each condition.
#my P_CFU        = shift; #
#my P10          = shift; # CV of X-Y pair-dependent growth amplitudes (log-normal distribution) in non-selectable condition is 10%
#my P_template   = shift; # Multiplification of the following ;  2.8 * 10**7 as experimental value
                         #  Yeast DNA mniprep yield from 3 × 10^7 diploid yeast cells is 20 ng
                         #  Fraction of Y2H plasmids in yeast DNA miniprep product is 6% of the total DNA mass
                         #  Barcode fusion efficiency is 20%
                         #  (DnDn & UpUp) x 4 PCR replicates
#my P16          = shift; # CV of X-Y pair-dependent PCR amplitudes (log-normal distribution) is 50%

    #print 'BEGIN: X-Y coverage in non-selectable condition'."\n";
  
#print Dumper "DONE";




}
}
#}


=C
=cut
















#############
# Functions #
#############

sub marginal(){
  my $data = shift;

  my $x_hap;
  my $y_hap;

  my $sum = 0;

  while(my($label,$amount) = each %$data){
    $label =~ s/\*//g;
    my ($x,$y) = split /\-/,$label;
    $x_hap->{$x} += $amount + 1;
    $y_hap->{$y} += $amount + 1;

    $sum += $amount + 1;
  }

  while(my($x,$amount) = each %$x_hap){
    $x_hap->{$x} = $amount/$sum;
  }
  while(my($y,$amount) = each %$y_hap){
    $y_hap->{$y} = $amount/$sum;
  }

  return ($x_hap,$y_hap);
}



sub print_table(){
  my $table = shift;

  for my $copy (sort{$a <=> $b} keys %$table){
    my $value = $table->{$copy};
    print "$copy\t$value\n";
  }
}

sub print_table_parameters(){
  my $table        = shift;
  my $reads        = shift;
  my $P1           = shift; # CV (coefficient of variation) of strain abundances (log-normal distribution) in haploid pools is 30%
  my $P3           = shift; # CV of X-Y pair-dependent mating efficiencies (log-normal distribution) is 50%
  my $P4           = shift;
  my $P5           = shift; # CV of X-Y pair-dependent growth amplitudes (log-normal distribution) in liquid media is 50%
  my $P_plated     = shift; # Amount of cells spread on each condition.
  my $P_CFU        = shift; #
  my $P10          = shift; # CV of X-Y pair-dependent growth amplitudes (log-normal distribution) in non-selectable condition is 10%
  my $P_template   = shift; # Multiplification of the following ;  2.8 * 10**7 as experimental value
                           #  Yeast DNA miniprep yield from 3 × 10^7 diploid yeast cells is 20 ng
                           #  Fraction of Y2H plasmids in yeast DNA miniprep product is 6% of the total DNA mass
                           #  Barcode fusion efficiency is 20%
                           #  (DnDn & UpUp) x 4 PCR replicates
  my $P16          = shift; # CV of X-Y pair-dependent PCR amplitudes (log-normal distribution) is 50%




  for my $copy (sort{$a <=> $b} keys %$table){
    my $value = $table->{$copy};
    print "$copy\t$value\t$reads\t$P1\t$P3\t$P4\t$P5\t$P_plated\t$P_CFU\t$P10\t$P_template\t$P16\n";
  }
}



sub print_logf_dist(){
  my $dist_f   = shift;
  my $bin_size = shift;
  my $low      = shift;
  my $high     = shift;

  for(my $i=$low;$i<$high;$i+=$bin_size){
    $i = sprintf "%.1f", $i;
    $i *= 1;
    my $value = $dist_f->{$i} || 0;
    print "$i\t$value\n";
  }
}


sub logf_dist(){
  my $data     = shift;
  my $bin_size = shift;
  my $alpha    = shift || 0;

  my $sum = 0;
  while(my($label,$amount) = each %$data){
    $sum += $amount+$alpha;
  }

  my $dist_f;
  my $total;
  while(my($label,$amount) = each %$data){
    my $log_f = sprintf "%.0f",log(($amount+$alpha)/$sum)/log(10)/$bin_size;
    $log_f *= $bin_size;
    $dist_f->{$log_f}++;
    $total++;
  }

  while(my($cell,$amount) = each %$dist_f){
    $dist_f->{$cell}/=$total;
  }

  return $dist_f;
}


sub coverage2(){

  my $data    = shift;
  my $target  = shift || 0;
  
  my $hap_x   = shift;
  my $hap_y   = shift;

  my $filter     = shift;

  my $coverage;

  my @array;


  
  while(my($label,$value) = each %$data){

    if($filter){
      if($label =~ /$filter/){
	push @array,$value;
      }
    }else{
      push @array,$value;
    }
  }

  #print Dumper @array;
  @array = sort{$b<=>$a} @array;
  
  
  my $count = 0;
  for my $i (0..$#array){
      if($array[$i] >=1){
	  $count += 1;
      }
  }


my $table;
my $num = 10**$target;
 
$table->{$num} = $count;
#last;
	  

return $table;
}






















sub auto_activity(){
  my $data      = shift;
  my $CV        = shift;
  my $popl_size = shift;


  my $data2;

  my $mem;
  while(my($label,$amount) = each %$data2){
    my ($x,$y) = split /\-/,$label;
    my $autoactivity = 0;
    if($mem->{$x}){
      $autoactivity = $mem->{$x};
    }else{
      $autoactivity = lograndn($CV);
    }
    my $amount2 = $amount * $autoactivity;
    $data2->{$label} = $amount2;
  }

  my $sum = 0;
  while(my($label,$amount) = each %$data2){
    $sum += $amount;
  }

  while(my($label,$amount) = each %$data2){
    $data2->{$label} = int($amount/$sum*$popl_size);
  }
  return $data2;
}




sub adjust_num(){
  my $data      = shift;
  my $CV        = shift;
  my $popl_size = shift;

 


  my $data2;

  while(my($label,$amount) = each %$data){
    my $amount2 = $amount * lograndn($CV);
    $data2->{$label} = $amount2;
  }

  my $sum = 0;
  while(my($label,$amount) = each %$data2){
    $sum += $amount;
  }

  while(my($label,$amount) = each %$data2){
    $data2->{$label} = int($amount/$sum*$popl_size);
  }
  return $data2;
}


sub random_selection(){
    my $data = shift;
    my $read = shift;

    my $sum = 0;
    while(my($label,$amount) = each %$data){
	$sum += $amount;
    }

    my $data2;
    my $rand_num = 0;
    my $rand_hash;
    while(my($label,$amount) = each %$data){
	$data2->{$label} = 0;

	my $s = $rand_num;
	my $e = $s + $amount;
	#my @array = ($s,$e);
	#print Dumper @a;
     
	$rand_hash->{$label} ->{0} = $s;
	$rand_hash->{$label} ->{1} = $e;

    
	$rand_num += $amount;}
    
    #print Dumper "$_\n" for keys $rand_hash; 

    for (my $i = 0; $i < $read; $i++) {
      
	my $random = rand();
	#print Dumper $random;
	#print Dumper $sum;
	$random *= $sum;
	my $rounded = int($random + 0.5);
	#print Dumper $random,$rounded, $rand_num, $sum;
	my $c = 0;
	foreach my $key (keys %$rand_hash){
	    $c += 1;
  
	    my $s = $rand_hash->{$key}-> {0};
	    my $e = $rand_hash->{$key}->{1};
	    
	    if( ($rounded > $s) && ($rounded <= $e) ){
		$data2 ->{$key} += 1;
	      } 
	    
	}
	#print Dumper $c,$i;



    }
    return $data2;
    
}
    

    #while(my($label,$amount) = each %$data2){
    #$data2->{$label} = int($amount/$sum*$popl_size);
    #}
    #return $data2;









sub mating(){
  my $hap_x      = shift;
  my $hap_y      = shift;

  my $dip;

  while(my($x,$x_value) = each %$hap_x){
    while(my($y,$y_value) = each %$hap_y){
      my $xy_value = $x_value * $y_value;
      $dip->{"$x\-$y"} = $xy_value;
    }
  }

  return $dip;
}



sub haploid(){
  my $size = shift;

  my $hap;

  for my $i (1..$size){
    $hap->{$i} = 1;
  }
  return $hap;
}



sub lograndn(){
  my $CV = shift;

  my $sigma = log(($CV**2)+1);
  $sigma **= 1/2;

  return exp(&randn(0,$sigma));
}

sub mean_sigma(){
  my $data  = shift;

  my @array = values %$data;

  my $mean = 0;
  for my $num (@array){
    $mean += $num;
  }
  $mean /= $#array+1;

  my $sigma = 0;
  for my $num (@array){
    $sigma += ($num-$mean)**2;
  }
  $sigma /= $#array+1;
  $sigma **= 1/2;

  my $CV = $sigma/$mean;

  return ($mean,$sigma,$CV);
}

sub randn(){

  my $m     = shift;
  my $sigma = shift;


  my ($r1, $r2) = (rand(),rand());
  while($r1==0){
    $r1 = rand();
  }
  my $value = ($sigma*sqrt(-2 * log($r1)) * sin(2 * PI * $r2)) + $m;
  return $value;
}
