# This script will combine all the replicates to generate a signle combined replicate

# Extracting the frag stats
my @frag_list = @ARGV[0..$#ARGV - 1];
print("Concating the following REfragStats: @frag_list\n");

# extracting the output file
my $out_fn = $ARGV[-1];
chomp ($out_fn);

# merge the split data
$i = 0;
while ($i <= $#frag_list){
        chomp $frag_list[$i];
        open (in,"$frag_list[$i]");
        while (<in>){
                chomp $_;
                @data = split(/\s+/,$_);
                if ($frag_list_check{@data[0]}{@data[1]}{@data[2]} eq ""){
                        $frag_list_check{@data[0]}{@data[1]}{@data[2]} = "yes";
                        push @frag_list_re,"@data[0]\t@data[1]\t@data[2]\n";
                }
		$j = 3;
		while ($j <= 11){
			chomp $data[$j];
			$frag_list_count{$j}{@data[0]}{@data[1]}{@data[2]} += $data[$j];
			$j++;
		}
                undef @data;
        }
        close in;
        $i++;
}

open (out,">$out_fn");
$i = 0;
while ($i <= $#frag_list_re){
        chomp $frag_list_re[$i];
        @frag_list_re_data = split(/\s+/,$frag_list_re[$i]);
        $chr   = @frag_list_re_data[0];
        $start = @frag_list_re_data[1];
        $end   = @frag_list_re_data[2];
        if ($chr ne "chrM" && $chr ne "chrY") {
          print out "$chr\t$start\t$end\t.";
          $j = 3;
          while ($j <= 11){
              $count = int($frag_list_count{$j}{$chr}{$start}{$end}/($#frag_list+1));
              print out "\t$count";
              $j
          }
          print out "\n";
        }
        undef @frag_list_re_data;
        $i++;
}
close out;
undef %frag_list_check;
undef %frag_list_count;
undef @frag_list_re;
undef @frag_list;
