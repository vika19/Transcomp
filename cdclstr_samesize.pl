#!/usr/bin/perl
# evigene/scripts/prot/cdclstr_samesize.pl : cd-hit .clstr parser, from cd-hit scripts

my $pMINLEN = $ENV{minlen} || 0.95;  # as prop 0.90 
   $pMINLEN=$pMINLEN/100 if($pMINLEN>1);
my $pMINID  = $ENV{minid}  || 0; # as pct 90% 
   $pMINID = 100*$pMINID if($pMINID > 0 and $pMINID < 1);
my $doTINY  = $ENV{tiny} || 0;
# my $doSAME  = $ENV{same} || 1;
## add no-same, self only rows ..
my $dosingle= $ENV{nosingle} ? 0 : 1;

my ($no,$nclin,$ncltr,$ngtr,$nitem,$nput,$topid,$topsize) = (0) x 9;
my @trow;

sub putclstr {
	my($topid,$topsize)= @_;
  $nclin++;
  if(@trow) { $ncltr++; $ngtr+= @trow; }
  elsif($dosingle) {
        push(@trow, join("\t","self",$topsize,100));
  }
	
	foreach my $tr (@trow) {
		my($tid,$tsize,$pid)=split"\t",$tr;
		my $pdif= $tsize/$topsize; $pdif= 1.0/$pdif if($pdif>1 and not $doTINY);
		
		my $putrow= ($doTINY) ? $pdif < $pMINLEN : $pdif >= $pMINLEN;
		#old# if($putrow) 
		if($putrow and $pid >= $pMINID)
		{
		 $pdif= int(1000*$pdif)/10;  $nput++;
		 print join("\t",$nclin,$topid,$topsize,$tid,$tsize, $pid, $pdif)."\n";
		}
	}
}

sub putsum {
  #?? ncl2=$ncltr, , nitem2nd=$ngtr
  my $type=($doTINY)?"Tiny":"Samesize";
  print "#Summary: $type nout=$nput, ncluster=$nclin, in_items=$nitem\n";
}

print join("\t",qw(IClstr TopID TopLen QueryID QueryLen pIden pSize))."\n";
while(my $ll=<> ) {
  if ($ll =~ /^>/) {
  	putclstr($topid,$topsize) if($topid or @trow>0); @trow=();
   	$topid = ""; $topsize= $no = 0;
  } else {
    chomp($ll);
    my ($pid,$tid,$tsize) = (0) x 9;
    if ($ll =~ /(\d+)(?:aa|nt), >(.+)\.\.\./) {
       $tsize = $1; $tid = $2;
       ($pid)= $ll =~ /at ([\d\.]+)/; # 0 for \*
    } else {
      die "cd-hit clstr format error: $ll";
    }

    if ($ll =~ /\*$/) {
      $topsize = $tsize; $topid= $tid;
    } else {
    	$pid=~s/100.00/100/; $pid=~s/(\.\d)\d/$1/;
    	push(@trow, join("\t",$tid,$tsize,$pid));
		}
		$nitem++;
  }
}

putclstr($topid,$topsize) if($topid or @trow>0);
putsum();

