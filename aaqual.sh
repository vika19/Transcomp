#!/bin/bash
# aaqual.sh any*.aa.gz : count protein sizes with faCount (dna), with tr for XXX gap count
# output 3-col table: ID, aatotal, XXgaps
# add on >aalen qual columns; ** Only good for evigene cdna_bestorf.aa files
# >whitefly:vel2k25Loc100074t1 aalen=41,63%,complete; clen=199; strand=-; offs=144-19; 
# ** FIXME2: w/ new aaqual, set col2 == size - gap, not sizetotal, and let user add gaps in

## option: offs= column, strand= ??
dostat=0; if [ "X$stat" != "X" ]; then dostat=$stat; fi
doff=0; if [ "X$off" != "X" ]; then doff=$off; fi
export doff=$doff;

inaa=$*
for az in $inaa; do {
  TCAT=cat; 
  nogz=`echo $az | sed 's/.gz//;'`; if [ $az != $nogz ]; then TCAT="gunzip -c"; fi
  nam=`echo $az | sed 's/.gz//; s/\.qual//;'`; 

  if [ ! -f $nam.qual ]; then 
  $TCAT $az | perl -ne \
'if(/^>(\S+)/) { puta() if($d); $d=$1; ($al)=m/aalen=([^;\s]+)/; $al||="na";
($cl)=m/clen=(\d+)/; $cl||=0; $aat=$aag=0; 
if($doff){ ($ofs)=m/offs=([\d-]+)/; $ofs||=0; $cl.="\t$ofs"; } }
else { s/\*$//; $aat += tr/A-WYZa-wyz/A-WYZa-wyz/; $aag += tr/Xx\*/Xx\*/; }
END{ puta(); } BEGIN{ $doff=$ENV{doff}; }
sub puta { $al=($aat+$aag).",na" if($al eq "na"); print join("\t",$d,$aat,$aag,$al,$cl)."\n"; }' \
  > $nam.qual

  fi

  # add stat here if desired:
  if [ $dostat != 0  ]; then
  pnam=`basename $nam | sed 's/\.aa//; s/.allcd//;'`
echo "# aa-quality for $pnam : top 1k summary"
cat $nam.qual | egrep -v '^#|^total' | sort -k2,2nr | head -1000 | env nam=$pnam perl -ne \
'($aw,$nn)=(split)[1,2]; $n++; $sw+=$aw; $sn+=$nn; push @aw,$aw;
END{ $aw=int($sw/$n); $an=int(10*$sn/$n)/10; @aw=sort{$b <=> $a}@aw ; ($mx,$md,$mi)=@aw[0,int($n/2),-1];
print "#$ENV{nam}\t n=$n; average=$aw; median=$md; min,max=$mi,$mx; sum=$sw; gaps=$sn,$an\n"; }'
  fi

} done

