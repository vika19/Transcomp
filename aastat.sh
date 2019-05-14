#!/bin/bash

if [ "X" = "X$top" ]; then top=1000; fi
inaa=$*

for atab in $inaa; do {
  if [ $atab = "stdin" ]; then atab="-"; fi
  if [ $atab = "-" ]; then 
    if [ "X" = "X$nam" ]; then nam="stdin"; fi
  else
    nam=`basename $atab .aa.qual | sed 's/\.count//; s/\.aa//;'`
  fi

echo "# aa-quality for $nam : longest $top summary"
cat $atab | egrep -v '^#|^total' | sort -k2,2nr | head -$top | env nam=$nam perl -ne \
'next unless(/^\w/); ($aw,$nn)=(split)[1,2]; $n++; $sw+=$aw; $sn+=$nn; push @aw,$aw;
END{ $aw=int($sw/$n); $an=int(10*$sn/$n)/10; @aw=sort{$b <=> $a}@aw ; ($mx,$md,$mi)=@aw[0,int($n/2),-1];
print "#$ENV{nam}\t  n=$n; average=$aw; median=$md; min,max=$mi,$mx; sum=$sw; gaps=$sn,$an\n"; }'

} done

