$N = 100;
$omega = 5;
$total = 2000 - $omega*100;
$j_max = $total/$N;
$l=1;
## i gives h indexing 

for $j (0..$j_max-1){
    print "round $j/$j_max: "; system "date";
    for $i ($j*$N..(($j+1)*$N)-1) {
            print "fork $i: "; system "date";
            my $pid = fork();
            if (not $pid) {
                    system("python trialParallelUpUpUp.py $i UpUpUp/trialParallelresulte$i.fits UpUpUp/trialParallelresulto$i.fits $omega $l");
                    exit();
            }
    }
    for $i (0..$N-1) {
            wait();
            print "waited $i/$N: "; system "date";
    }
}
##for $i (0..$N-1) {system("cat trialParallelresulte$i.fits > testingParallele.fits");}
##for $i (0..$N-1) {system("cat trialParallelresulto$i.fits > testingParallelo.fits");}
