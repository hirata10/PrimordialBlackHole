# Script to run several configuration files in parallel.
# Useful for running on a multi-core machine.
#
# It has to be called *from one directory up* (so from pyimcom/ need to do a cd .. before running)

($M, $nstep, $istart, $Nr) = @ARGV;

@omegas = (0.05, 0.06, 0.07, 0.09, 0.11, 0.14, 0.17, 0.2, 0.25, 0.31, 0.37, 0.46,
0.56, 0.68, 0.84, 1.02, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75,
3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75,
6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.5, 9.0, 9.5,
10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0);

print (scalar @omegas);
print "\n";

my $i=0;
for $i (0..$Nr-1) {
  my $pid=fork;
  if (not defined $pid) {
    print STDERR "Error: fork failed\n";
    exit;
  }
  if (not $pid) {
    $j = $i+$istart;
    $command1 = "python correction_factor.py $M $omegas[$j] $nstep > out/correction_$M\_$omegas[$j].dat";
    print "[$i,$pid] $command1\n";
    system "$command1";
    print "[$i,$pid] Finished.\n";
    exit;
  }
}

# Wait for children
my $k;
for $k (1..$Nr) {wait();}
