#!/usr/bin/perl
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);

my %opts;
# 现在只接收 -i 和 -o
GetOptions(\%opts,"i=s","o=s","h"); 

if(!defined($opts{i}) || !defined($opts{o}) || defined($opts{h})) {
    print "Usage: perl $Script -i fragments.tsv.gz -o /path/to/prefix\n";
    print "Example: perl $Script -i input.gz -o /results/M8.filter4\n";
    exit;
}

my $umifile = $opts{i};
my $base_prefix = $opts{o}; # 这里的 -o 是完整路径前缀
my $outfile = "${base_prefix}_result.txt"; 
my $bc_dist = "${base_prefix}.barcode_dist.tmp";

# --- Step 0: 预扫描 ---
&show_log("Step 0: Pre-scanning fragments...");
open (PRE, "gzip -dc $umifile |" ) || die "Error: Cannot open $umifile: $!";
my %bc_all = ();
while(<PRE>){
    next if (m/^\#/);
    my ($chr, $start, $end, $bc, $num) = split;
    next unless ($chr =~ /^chr(?:[1-9]|1[0-9]|X)$/); 
    $bc_all{$bc} += $num; 
}
close PRE;

open (DIST, ">$bc_dist");
print DIST "Barcode\tnFrags\n";
foreach (sort {$bc_all{$b} <=> $bc_all{$a}} keys %bc_all) {
    print DIST "$_\t$bc_all{$_}\n";
}
close DIST;

# --- Step 0.1: 动态寻秩 ---
&show_log("Step 0.1: Detecting Cutoff in R...");
my $cutoff = `Rscript $Bin/draw_ATAC_refined.R $outfile $bc_dist $base_prefix FIND_CUTOFF`;
chomp($cutoff);
$cutoff ||= 1000;
&show_log("  >> Final Cutoff applied: $cutoff fragments.");

# 保存 cutoff 值到文件，供后续步骤使用
my $cutoff_file = "${base_prefix}.cutoff.txt";
open (CUTOFF, ">$cutoff_file") || die "Cannot write cutoff file: $!";
print CUTOFF "$cutoff\n";
close CUTOFF;
&show_log("  >> Cutoff saved to: $cutoff_file");

# --- Step 1: 信号过滤 ---
my %valid_bcs = ();
foreach (keys %bc_all) { 
    $valid_bcs{$_} = 1 if $bc_all{$_} >= $cutoff; 
}
my $signal_filter_count = scalar(keys %valid_bcs);
&show_log("  >> Signal Filter Pass: $signal_filter_count spots.");

# 输出信号过滤后的 barcode list
my $signal_bc_file = "${base_prefix}.barcode_pass_signal.txt";
open (SIGBC, ">$signal_bc_file");
foreach my $bc (sort keys %valid_bcs) {
    print SIGBC "$bc\n";
}
close SIGBC;
&show_log("  >> Signal filter barcodes saved to: $signal_bc_file");

# --- Step 2: 加载数据用于抽样 ---
&show_log("Step 2: Loading valid fragments for sampling...");
open (IN, "gzip -dc $umifile |" ) || die $!;
my ($total_reads, $index, @simu_arr, @idx_to_bc) = (0, 0, (), ());
while(<IN>){
    next if (m/^\#/);
    my ($chr, $start, $end, $bc, $num) = split;
    next unless ($chr =~ /^chr(?:[1-9]|1[0-9]|X)$/); 
    next unless exists $valid_bcs{$bc};
    $index++;
    $total_reads += $num;
    push @idx_to_bc, $bc;
    push @simu_arr, ($index)x$num;
}
close IN;

# --- Step 3: Shuffle & Sampling ---
&show_log("Step 3: Shuffling...");
&shuffle(\@simu_arr);

&show_log("Step 4: Sampling...");
my ($percent, $current_uniq, @is_seen, %hstat, %bc_counts) = (12.5, 0, (), (), ());
for (my $i=0; $i<@simu_arr; $i++){
    my $f_idx = $simu_arr[$i];
    if (!defined $is_seen[$f_idx]){
        $is_seen[$f_idx] = 1;
        $current_uniq++;
        $bc_counts{$idx_to_bc[$f_idx-1]}++;
    }
    if ($i+1 >= $total_reads * $percent / 100 || $i+1 == @simu_arr){
        my @counts = values %bc_counts;
        my $median = &calculate_median_fixed(\@counts, $signal_filter_count);
        $hstat{$percent} = [$current_uniq, $i+1, (1-$current_uniq/($i+1))*100, $median];
        $percent += 12.5;
        last if $percent > 100.1;
    }
}

# --- Step 5: 输出与绘图 ---
open (OUT, ">$outfile");
print OUT "Percent\tUnique_Fragments\tTotal_Fragments\tPercent_Duplicates\tMedian_nFrags\n";
foreach my $p (sort {$a<=>$b} keys %hstat){
    printf OUT "%.1f\t%d\t%d\t%.2f\t%.2f\n", $p, @{$hstat{$p}};
}
close OUT;

&show_log("Step 5: Generating PDF report...");
system("Rscript $Bin/draw_ATAC_refined.R $outfile $bc_dist $base_prefix DRAW $cutoff");

# 保留 barcode_dist.tmp 给后续步骤使用
# unlink $bc_dist;
&show_log("Done. Results: $outfile");
&show_log("  >> Barcode distribution: $bc_dist");
&show_log("  >> Saturation stats: $outfile");

sub calculate_median_fixed {
    my ($arr, $total) = @_;
    my @all = (@$arr);
    push @all, (0) x ($total - @all) if $total > @all;
    my @s = sort {$a<=>$b} @all;
    my $m = int(@s/2);
    return @s % 2 ? $s[$m] : ($s[$m-1]+$s[$m])/2;
}

sub shuffle {
    my $a = shift;
    for (my $i = @$a - 1; $i > 0; $i--){
        my $j = int(rand($i + 1));
        @$a[$i,$j] = @$a[$j,$i];
    }
}

sub show_log {
    my ($y,$m,$d,$h,$mi,$s) = (localtime)[5,4,3,2,1,0];
    printf("[%04d-%02d-%02d %02d:%02d:%02d] %s\n", $y+1900,$m+1,$d,$h,$mi,$s, "@_");
}
