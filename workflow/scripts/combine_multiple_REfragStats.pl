# Super fast perl script to merge REfragStats files
# Usage: perl combine_multiple_REfragStats.pl <txt file with a list of REfragStat paths> <combined output path>
my @list = `cat $ARGV[0]`;
my $combined_file = $ARGV[1];
my %merged;
my %check;
my @bed;
foreach my $l (@list) {
    chomp $l;
    print ("$l\n");
    open(in, "$l");
    while(my $frag = <in>) {
        chomp $frag;
        my $chr = (split(/\s+/,$frag))[0];
        my $start = (split(/\s+/,$frag))[1];
        my $end   = (split(/\s+/,$frag))[2];
        $merged{$chr}{$start}{$end}{1} += (split(/\s+/,$frag))[3];
        $merged{$chr}{$start}{$end}{2} += (split(/\s+/,$frag))[4];
        $merged{$chr}{$start}{$end}{3} += (split(/\s+/,$frag))[5];
        $merged{$chr}{$start}{$end}{4} += (split(/\s+/,$frag))[6];
        $merged{$chr}{$start}{$end}{5} += (split(/\s+/,$frag))[7];
        $merged{$chr}{$start}{$end}{6} += (split(/\s+/,$frag))[8];
        $merged{$chr}{$start}{$end}{7} += (split(/\s+/,$frag))[9];
        $merged{$chr}{$start}{$end}{8} += (split(/\s+/,$frag))[10];
        $merged{$chr}{$start}{$end}{9} += (split(/\s+/,$frag))[11];
        if ($check{$chr}{$start}{$end} eq "" ){
            push(@bed, "$chr\t$start\t$end");
            $check{$chr}{$start}{$end} = 1;
        }
    }
    close in;
}
open(out,">$combined_file");
foreach my $b (@bed) {
    chomp $bed;
    my $chr = (split(/\s+/,$b))[0];
    my $start = (split(/\s+/,$b))[1];
    my $end   = (split(/\s+/,$b))[2];
    print out ("$b");
    for($i = 1; $i <= 9; $i++){
        print out ("\t$merged{$chr}{$start}{$end}{$i}");
    }
    print out ("\n");
}
close out;
