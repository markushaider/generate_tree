#!/usr/bin/perl -w
use lib qw(src);
use ReadConfig;
use CachedReadWrite;
use MultiThread;


load_config();

opendir DIR, $TREE_OUTBASE;
my @trees = grep { /^tree_.*\.dat$/ } readdir DIR;
closedir DIR;

foreach (@trees) {
    next unless (MultiThread::ForkChild());
    convert_to_catalog("$TREE_OUTBASE/$_");
    exit;
}
MultiThread::WaitForChildren();

sub convert_to_catalog 
{
    @ARGV = shift;
    while (<>) {
	next unless /\#\s*tree/i; #Skip header
	last;
    }

    while (<>) {
	if (/^\#\s*tree/) {
	    process_tree();
	    %halos = ();
	    @halos = ();
	    next;
	}
	next unless (/^\s*\d+\.?\d*\s+\d+/);
	my $h = new TreeHalo;
	$h->scanline($_);
	push @halos, $h;
	$halos{$h->{scale}}{$h->{id}} = $h;
    }
    process_tree();
    $_->close() foreach (values %TreeHalo::tree_outputs);
}

sub max {
    return (($_[0] > $_[1]) ? $_[0] : $_[1]);
}

sub process_tree {
    for my $h (@halos) {
	$h->{num_prog} = 0;
    }
    for my $h (@halos) {
	my $d = $halos{$h->{desc_scale}}{$h->{descid}};
	next unless defined($d);
	$d->{num_prog}++;
	my $p = $d->{prog};
	if (defined($p) and $p->{mvir} > $h->{mvir}) {
	    next;
	}
	$d->{prog} = $h;
    }
    for my $h (@halos) {
	calc_mass_vmax_acc($h);
	$h->print();
    }
}

sub calc_mass_vmax_acc {
    no warnings 'recursion';
    my $h = shift;
    return if ($h->{seen});
    $h->{seen} = 1;
    if ($h->{prog}) {
	calc_mass_vmax_acc($h->{prog});
	if ($h->{upid}>-1) { # If we are a subhalo
	    $h->{vacc} = $h->{prog}{vacc};
	    $h->{macc} = $h->{prog}{macc};
	} else {
	    $h->{vacc} = $h->{vmax};
	    $h->{macc} = $h->{mvir};
	}
	$h->{vpeak} = max($h->{vmax}, $h->{prog}{vpeak});
	$h->{mpeak} = max($h->{mvir}, $h->{prog}{mpeak});

	#Vpeak / Mpeak *before* accretion
	#$h->{vpeak} = max($h->{vacc}, $h->{prog}{vpeak});
	#$h->{mpeak} = max($h->{macc}, $h->{prog}{mpeak});
    } else {
	$h->{vpeak} = $h->{vmax};
	$h->{vacc} = $h->{vmax};
	$h->{mpeak} = $h->{mvir};
	$h->{macc} = $h->{mvir};
    }
}

sub load_config {
    my $config = new ReadConfig($ARGV[0]);
    $config->set_defaults(
	SCALEFILE => "/Volumes/Peter 1/Bolshoi/DescScales.txt",
	TREE_OUTBASE => "/Volumes/Peter 2/Bolshoi/Trees",
	HLIST_OUTBASE => "/Volumes/Peter 2/Bolshoi/Hlists");
    our $TREE_OUTBASE = $config->{TREE_OUTBASE};
    our $HLIST_OUTBASE = $config->{HLIST_OUTBASE};
    our $OUTLIST = $config->{SCALEFILE};
#    our $SINGLE_THREAD_OUTPUT = $config->{SINGLE_THREAD_OUTPUT};
}

package TreeHalo;

sub _open_scale {
    my $scale = shift;
    my $fn = $main::HLIST_OUTBASE.sprintf("/hlist_%.5f.list", $scale);
    $tree_outputs{$scale} = new CachedReadWrite(file => $fn, cachesize => 1000,
						header => "#Scale(0) Id(1) Desc_scale(2) Descid(3) Num_prog(4) Pid(5) Upid(6) Desc_pid(7) Phantom(8) Mvir(9) Orig_Mvir(10) Rvir(11) Rs(12) Vrms(13) Mmp(14) Last_mm(15) Vmax(16) X(17) Y(18) Z(19) VX(20) VY(21) VZ(22) Macc(23) Mpeak(24) Vacc(25) Vpeak(26) JX(27) JY(28) JZ(29) Spin(30)\n");
    $tree_outputs{$scale}->set_append(1);
}

sub new {
    my $class = shift;
    my $obj = { map { $_ => 0 }
		qw(id descid scale desc_scale
num_prog pid upid desc_pid mtype mmp flags phantom
mvir rvir vmax vrms rs mtrunc orig_mvir rtrunc last_mm np spin) };
    $obj->{pos} = [0,0,0];
    $obj->{vel} = [0,0,0];
    $obj->{J} = [0,0,0];
    bless $obj, $class;
}

sub scanline {
    my $h = shift;
    my $line = shift;
    ($h->{scale}, $h->{id}, $h->{desc_scale}, $h->{descid}, $h->{num_prog},
     $h->{pid}, $h->{upid}, $h->{desc_pid}, $h->{phantom},
     $h->{mvir}, $h->{orig_mvir}, $h->{rvir}, $h->{rs}, $h->{vrms},
     $h->{mmp}, $h->{last_mm}, $h->{vmax},
     $h->{pos}[0], $h->{pos}[1], $h->{pos}[2],
     $h->{vel}[0], $h->{vel}[1], $h->{vel}[2],
     $h->{J}[0], $h->{J}[1], $h->{J}[2], $h->{spin}) = split(" ", $line);
    $h->{mvir} = abs($h->{mvir});
    return "$h->{scale} $h->{id}";
}

sub print {
    my $h = shift;
    return unless (defined $h->{scale} and $h->{scale} > 0);
    _open_scale($h->{scale}) if (!exists($tree_outputs{$h->{scale}}));
    my $file = $tree_outputs{$h->{scale}};
    $file->printf("%.4f %8d %.4f %8d %6d %8d %8d %8d %2d %.5e %.5e %6f %6f %6f %2d %.4f %6f %.5f %.5f %.5f %.3f %.3f %.3f %.5e %.5e %6f %6f %.3e %.3e %.3e %.5f\n",
    $h->{scale}, $h->{id}, $h->{desc_scale}, $h->{descid}, $h->{num_prog},
    $h->{pid}, $h->{upid}, $h->{desc_pid}, $h->{phantom},
    $h->{mvir}, $h->{orig_mvir}, $h->{rvir}, $h->{rs}, $h->{vrms},
    $h->{mmp}, $h->{last_mm}, $h->{vmax},
    $h->{pos}[0], $h->{pos}[1], $h->{pos}[2],
    $h->{vel}[0], $h->{vel}[1], $h->{vel}[2],
    $h->{macc}, $h->{mpeak}, $h->{vacc}, $h->{vpeak},
    $h->{J}[0], $h->{J}[1], $h->{J}[2], $h->{spin});
}

sub print_mass_vmax_acc {
    my $h = shift;
    printf "%.4f %d %.4e %.2f %.4e %.2f %.4e %.2f %.5f %.5f %.5f %.2f %.2f %.2f\n",
    $h->{scale}, $h->{id}, $h->{mvir}, $h->{vmax}, $h->{macc}, $h->{vacc}, $h->{mpeak}, $h->{vpeak}, @{$h->{pos}}, @{$h->{vel}};
}
