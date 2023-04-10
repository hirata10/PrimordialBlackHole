$N = 10;


for $i (0..$N-1) {
        print "fork $i: "; system "date";
        my $pid = fork();
        if (not $pid) {
                system("python playground.py $i > bispec_components.$i");
                exit();
        }
}
for $i (0..$N-1) {
        wait();
        print "waited $i/$N: "; system "date";
}

for $i (0..$N-1) {system("cat bispec_components.$i");}

