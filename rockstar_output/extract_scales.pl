#!/usr/bin/perl -w

for (@ARGV) {
    next unless /^out_(\d+).list$/;
    my $id = $1;
    my $scale;
    open FILE, "<", $_;
    while ($_ = <FILE>) {
if (/a = (\d\.\d+)/) {
$scale = $1;
}
if (/^\d+/) { $scale{$id} = $scale; last; }
    }
    close FILE;
}

my @ids = sort {$scale{$a} <=> $scale{$b}} keys %scale;
print "$_ $scale{$_}\n" for (@ids);
