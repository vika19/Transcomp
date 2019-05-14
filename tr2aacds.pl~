#!/usr/bin/perl
# tr2aacds.pl

=item about

  EvidentialGene/tr2aacds.pl: 
  convert mRNA assembly sets of too many transcripts, to best open
  reading frame coding sequences, with filtering of identical coding sequences, 
  and alternate CDS transcripts found by high identity subset alignments. 

  inputs : one file of all transcripts.fasta[.gz] ; optional input transcripts.aa,.cds
  outputs: best orf aa and cds, classified into okay and drop sets.

  This omnibus script comes out of many tests and refinements for
  selecting best transcript assemblies from a superset of many, many
  assemblies of same RNAseq data.  The focus on CDS of mRNA transcripts
  allows for several of those refinements in balancing too few/too long mistakes,
  and too many subset/nearly same assemblies.  The goal is a biologically
  valid transcript assembly set.

  The input mRNA assembly set is presumed a large collection from many
  assemblers, data slices and options from the same transcript source.
  Typically 1+ million input transcripts are processed, and reduced to a 
  biologically realistic set of 20k - 50k main plus similar number of
  alternate transcripts.
  
  CDS sequence is primary quality of mRNA transcripts, after trying/discarding 
  AA sequence (for those missed silent codon paralog changes).  CDS and AA
  are guessed as ~longest ORF of transcript (with sizes options for complete/partial).
  Long UTRs are also scanned for ORFs, of joined genes (common w/ hi-express neighbors,
  some assemblers).   
  
  CDS are classified by identity/alignment to each other, with fastanrdb (identicals), 
  cd-hit-est (ident. fragments) and blastn (local hi-identity alternates).
  High identity redundant assemblies are discarded, alternates classified
  by hi-identity local alignments are classified and CDS-duplicative ones
  discarded.   A final okay set contains main transcripts with informative
  alternates, and unique CDS (no alternate or low identity).  Minimum CDS/AA sizes
  are used. Protein completeness/partial, and UTR-poor qualities also
  score and filter excess assemblies.
  
  Transcript as whole is ignored due to many UTR misassemblies and difficulty
  in assessing quality, but for CDS/UTR ratio as quality measure.  Test show
  selection of longest CDS-ORF has strong +correlation with highest 
  protein orthology score.

=item example

    $evigene/prot/tr2aacds.pl -debug -NCPU $ncpu -MAXMEM $maxmem -log -cdna $trset
    See below tr2aacds_qsub.sh for cluster script.

=item requirements

  See source, but currently:
  EvidentialGene source tree: http://arthropods.eugenes.org/EvidentialGene/evigene/
            or ftp://arthropods.eugenes.org/evigene/  (all: evigene.tar)            
  blastn,makeblastdb : NCBI C++ blast (tested ncbi2227)
  cd-hit-est, cd-hit  : http://cd-hit.org/
  fastanrdb : exonerate fastanrdb, create non-redundant fasta sequence database

  These need to be in PATH or ENV. 

=item more

  Use of high-identity CDS filtering (blastn identity >=98%; nrdb/cdhit-est 100%) 
  seems to balance well removal of true same-CDS and retention of paralogs and clonal
  different transcripts, however checking may be needed.
  
  > inputs : one file of all transcripts.fasta[.gz]
  In theory this can be run subset transcript assemblies, then
  re-run on okay outputs to combine.  This script currently uses
  -NCPU for parallelization on clusters.  The "slow" computes are
  for blastn-self of large seq set, and cd-hit-est of same.  Both use NCPU
  effectively on one compute node w/ many cores.  Running this on subset
  assemblies (as generated) may be effectively same, and more parallel/useful.
  
  TODO: Step 0. cdna_bestorf can be time-consuming, and can be parallelized by
  splitting input tr set.
  
  Runtime is ~3hr for ~1 Million transcripts on 32-core node.
    (0:bestorf 1hr, 3:blastn 1hr,  )
  
  
=item see also

  pipelining of various scripts from evigene/scripts/rnaseq/  
  
  okayset/*.tr may be submitted to NCBI TSA with
    evigene/scripts/rnaseq/asmrna2ncbitsa.pl
    : fixme, needs update for tr2aacds output formats
    
=item author
  
  don gilbert, gilbertd near indiana edu, 2013
  part of EvidentialGene, http://arthropods.eugenes.org/EvidentialGene/

=cut

use constant VERSION => '2013.07.27'; # '04.15'; #04.07; 03.14'; # '.03.11'; '.03.06';
## 11mar.FIXME: need blastn -ungapped otherwise miss perfect match of parts ; e.g. alt-tr half=perfect, other=imperf
## 14mar: blast_ncpu using fasplit(query.fa,ncpu) instead of blastn -num_tasks ..
## 07apr: add -ablastab homolog scores for classifier
## 15apr: lastz not ready to replace blastn; testing option..
## 27jul: lastz added, maybe ready/better than blastn, maybe not .. compare more w/ blastn

use FindBin;
use lib ("$FindBin::Bin"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...
# this script is now in evigene/scripts/prot/..

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname);
# use cdna_proteins; #?
# use cdna_evigenesub; #?

## FIND evigene path = self path
# my $EVIGENES="$FindBin::Bin/..";  
#patch for wacky cray bigdog2# do we also need 2nd use lib($EVIGENES) ??
my $EVIGENES=$ENV{EVIGENES} || "$FindBin::Bin/.."; #?? cray-bigdog2 failing at this ??? stupid machine

my $MINCDS = $ENV{MINCDS} || 90; # what? #maketraa3.sh: $AAMIN=40; $AMINPOO=100; $AMINBAD=200
my $CDSBLAST_IDENT= 98; # 98 + -ungapped ; was 95; # or 98; # or 95? ** drop this back to 95 default; at least for high okay count species
## .. may be getting high okay main/noclass gene count from mixed strain/heterozygotes? need some checking.
my $CDSBLAST_EVALUE= 1e-19; # is this ok?
my $NCPU= 1;
my $MAXMEM= 1000; # in Mb always 

# my $prefix = $ENV{prefix} || ""; # what?
my ($debug,$dryrun,$tidyup,$runsteps,$USE_LASTZ)= (0) x 9;  
my ($logfile,$aaseq,$aasize,$cdnaseq,$cdsseq,$aacdseq,$aaclstr,$aablast)= (undef) x 20; 

#?? rename?  outsetokay, outsetdrop, inputset, tmpfiles
my (@okayset,@dropset,@inputset,@tmpfiles,@erasefiles); # tidyup file sets
# my @input=();

use constant { LOG_NOTE => 0, LOG_WARN => 1, LOG_DIE => -1, };
my $logh= undef;
sub loggit{ my $dowarn=shift; 
  my $s= join(' ',@_); chomp($s); $s="FATAL $s" if($dowarn == LOG_DIE);
  if($logh){ print $logh "#t2ac: $s\n"; } elsif($dowarn>0||$debug){ warn "#t2ac: $s\n"; }
  if($dowarn == LOG_DIE) { die "#t2ac: $s\n" ; }
}

my $optok= GetOptions( 
  "cdnaseq|mrnaseq|trinput=s", \$cdnaseq, ## \@input,  # one only for this?
  "aaseq|aainput:s", \$aaseq,
  "cdsseq|cdsinput:s", \$cdsseq,
  "logfile:s", \$logfile,
  "ablastab=s", \$aablast,   # option traa-refaa.blastp.tall4 table for asmrna_dupfilter2.pl classifier
  # "runsteps=s", \$runsteps,  ## general options here?  -runsteps=lastz,xxx,
  # "format=s", \$format, 
  # "prefix=s", \$prefix, ## == prefix for new IDs from trclass
  "MINCDS=i", \$MINCDS,  
  "CDSBLAST_IDENT=i", \$CDSBLAST_IDENT, "CDSBLAST_EVALUE=i", \$CDSBLAST_EVALUE,  
  "NCPU=i", \$NCPU, "MAXMEM=i", \$MAXMEM,  
  "uselastz|lastz!", \$USE_LASTZ, # instead of blastn for cds-align
  "tidyup!", \$tidyup, 
  "dryrun|n!", \$dryrun,
  "debug!", \$debug,
);

#?? push @input, @ARGV; # FIXME: do something w/ remaining @ARGV .. warn?

die "EvidentialGene tr2aacds.pl VERSION ",VERSION,"
  convert large, redundant mRNA assembly set to best protein coding sequences, 
  filtering by quality of duplicates, fragments, and alternate transcripts.
  See http://eugenes.org/EvidentialGene/about/EvidentialGene_trassembly_pipe.html
Usage: tr2aacds.pl -mrnaseq transcripts.fasta[.gz] 
  opts: -MINCDS=$MINCDS -NCPU=$NCPU -MAXMEM=$MAXMEM.Mb -logfile -tidyup -dryrun -debug 
" unless($optok and $cdnaseq); 

if(not $logfile and defined $logfile) { # use output name
  $logfile= makename($cdnaseq,".tr2aacds.log");  
}
if($logfile) { open(LOG, ">>$logfile") or die $logfile; $logh= *LOG; }

$tidyup= 1 unless($dryrun||$debug); # default on unless debug|dryrun ?

loggit(1, "EvidentialGene tr2aacds.pl VERSION",VERSION);
loggit(1, "ERR: unused arguments:",@ARGV) if(@ARGV>0);

my $APPblastn=    findapp("blastn"); 
my $APPmakeblastdb= findapp("makeblastdb");
my $APPlastz="echo MISSING_lastz"; 
 	 $APPlastz= findapp("lastz") if ($USE_LASTZ || $debug); 
	 ## lastz 04.15 testing still; 2013.03.24 : replace blastn default for basic local align: better at finding perfect local aligns
my $APPfastanrdb= findapp("fastanrdb");
my $APPcdhitest=  findapp("cd-hit-est"); 
my $APPcdhit=     findapp("cd-hit");  
## .. these should call findevigeneapp() to warn/die if missing
my $APPcdnabest= findevigeneapp("cdna_bestorf.pl"); # allow ENV/path substitutions?
my $APPtraa2cds= findevigeneapp("prot/traa2cds.pl");
my $APPtrdupfilter= findevigeneapp("rnaseq/asmrna_dupfilter2.pl");  
my $APPaaqual=    findevigeneapp("prot/aaqual.sh");
## FAIL at this point if any apps missing?
#-------------------------------------

=item  tr2aacds pipeline algorithm

  0. make/collect input asm_name.{tr,aa,cds}, working mostly on .cds here

  1. perfect redundant removal:  fastanrdb  input.cds > input_nr.cds
    .. use aa.qual info for choosing top cds among identicals: pCDS only? ie need trsize along w/ cdssize

  2. perfect fragment removal: cd-hit-est -c 1.0 -l $MINCDS ..

  3. blastn, basic local align hi-ident subsequences for alternate tr.

  4. classify main/alternate cds, okay & drop subsets, using evigene/rnaseq/asmrna_dupfilter2.pl
     .. merges alignment table, protein-quality and identity, to score okay-main, ok-alt, and drop sets.

  5. make final output files from outclass: okay-main, okay-alts, drops 
      okayset is for public consumption, drops for data-overload enthusiasts (may contain valids).

  -- maybe option to run parts, w/ diff inputs, per offload blast,cdhit steps to cluster computers
     currently checks for intermediate files, so can be rerun to next step.

=cut

sub MAIN_start {}
MAIN: {
  # final output is classify table from asmrna_dupfilter2
  # and seq filesets classified: perfect dups, perfect fragments, alt partial dups : keep,drop class

#  my $check_outclass= makename($cdnaseq,".trclass");
#  if(-s $check_outclass and -d inputset) {
#    # skip/warn/..
#  }
  
  my($cdsseqnr,$cdsseqnrcd1,$cdsblast,$outclass,$outaln); # intermediate/output files
  my(%perfect_dups,%perfect_frags);
  
  loggit(0, "BEGIN with cdnaseq=",$cdnaseq,"date=",`date`);
  # fail unless -s $cdnaseq; .. is STDIN ok for this?

  # 0. make/collect input asm_name.{tr,aa,cds}, working mostly on .cds here
  # FIXMEd: parallelize this for NCPU by splitting input cdnaseq to NCPU parts; cat parts> whole.
  ($cdsseq,$aaseq) = get_bestorf($cdnaseq,$aaseq,$cdsseq);
  loggit(0, "bestorf_cds=",$cdsseq,"nrec=",facount($cdsseq));
  # FIXME: Check DUP IDs in aaseq, FAIL if found ..
  ($aasize)= make_aaqual($aaseq);
  if(my $ndup= fadupids($aaseq)) { loggit(LOG_DIE,"ERR: $ndup duplicate ids in $aaseq\n"); }
  
  # these subtasks will check existence of output file before run task
  # 1. nonredundant removal:  fastanrdb  input.cds 
  ($cdsseqnr)= nonredundant_cds($cdsseq);
  loggit(0, "nonredundant_cds=",$cdsseqnr,"nrec=",facount($cdsseqnr));

  ## 1.1. reassign redundant cds from aaqual = pCDS to reduce utrbad/utrpoor set.
  ## .. need only change header ID in $cdsseqnr ; fastanrdb puts all redundant ids on header.
  ## Note this rewrites $cdsseqnr file, to same name
  my($nbest,$ndups)= nonredundant_reassignbest($cdsseqnr,$aasize);
  loggit(0,"nonredundant_reassignbest=",$nbest,"of",$ndups); 
  
  # 2. perfect fragment removal  : clusterize : $NCPU, $MAXMEM
  ($cdsseqnrcd1)= nofragments_cds($cdsseqnr);
  loggit(0, "nofragments_cds=",$cdsseqnrcd1,"nrec=",facount($cdsseqnrcd1));
  # 2.1. FIXME need to check the cdsseqnrcd1.clstr clusters for useful alternate tr (?) 
  #   .. slightly shorter CDS may be valid alts

	if($USE_LASTZ) {
  # 3z. lastz replace blastn for alignments of hi-ident subsequences  
	($cdsblast)= lastz_cds($cdsseqnrcd1);  # asmrna_dupfilter2.pl can now parse lastz.general.format
 	loggit(0, "lastz_cds=",$cdsblast);
	} else {
  # 3. blastn alignments of hi-ident subsequences : clusterized  : $NCPU, $MAXMEM
  ($cdsblast)= blastn_cds($cdsseqnrcd1);
  loggit(0, "blastn_cds=",$cdsblast);
  }

  # 4.1 cdhit -c 0.9 -i $aaseq -o $aaseqcd .. for aaseqcd.clstr input to dupfilter ..
  #  .. filter uses aaclstr only to drop aa-hiident alts as uninformative; eg 1/2 of 13k althi.
  ($aacdseq,$aaclstr)= aacluster($aaseq); # but use only aaseq IDs from cdsseqnrcd1
   
  # 4. classify main/alternate cds, okay & drop subsets, using evigene/rnaseq/asmrna_dupfilter2.pl
  ($outclass,$outaln)= asmdupfilter_cds($cdsblast,$aasize,$aaclstr);
  loggit(0, "asmdupfilter_cds=",$outclass);
  # FIXME? here, or in asmrna_dupfilter2.pl, create ID main,alt table from outclass/trclass, 
  #    with new numeric IDs, old/cur ids, main/alt num
  ## Better leave to other script, evgrna2genbanktsa.pl, merging trclass, naming table, other annots ?
  
  my $OUTH= (defined $logh) ? $logh : *STDOUT;
  asmdupclass_sum($OUTH,$outclass,$aasize); 
  
  # 5. make final output files from outclass: okay-main, okay-alts, drops 
  # 5. fixme: add hdr classinfo from drops: %perfect_dups, %perfect_frags, 
  %perfect_dups = redundant_idset($cdsseqnr); # read headers after rebest, for dupclass_fileset
  %perfect_frags= fragment_idset("$cdsseqnrcd1.clstr"); # read headers after rebest, for dupclass_fileset

  my @outfiles= asmdupclass_fileset($outclass,$cdnaseq,$aaseq,$cdsseq,\%perfect_dups,\%perfect_frags); 
  loggit(0, "asmdupfilter_fileset=", @outfiles);

  ## turn off tidyup unless( -s $outclass);
  if($tidyup and -s $outclass) {
    loggit(0, "tidyup output folders: okayset dropset inputset tmpfiles");
    ## tidyup needs tod/basename($fn) ?? use perl:rename($fn,"$tod/$tf") ? log nfiles moved not each filename?
    sub tidyup{ my($tod,@td)= @_; mkdir($tod); foreach my $fn (@td) { my $tf=basename($fn); runcmd("mv $fn $tod/$tf") if(-f $fn); } }
    tidyup("okayset",@okayset);  
    tidyup("dropset",@dropset);  
    tidyup("inputset",@inputset);  
    tidyup("tmpfiles",@tmpfiles);  
    ## my $rmlist= join" ",grep{ -f $_ } @erasefiles;
    ## loggit(0,"tidyup erase:",$rmlist) if($rmlist); #too long for log .. chop
    my @rmlist;
    foreach my $fn (@erasefiles) { if(-f $fn) { unlink($fn); push @rmlist,$fn; } } # too verbose: runcmd("rm $fn") 
    if(@rmlist) { my $nrm=@rmlist; my $rml=join" ",@rmlist[0..4]; loggit(0,"tidyup erase: n=$nrm, $rml .."); } 
  }
  loggit(0, "DONE at date=",`date`);
  loggit(0, "======================================"); # log may append runs
}


#-------------------------------------

sub nonredundant_cds
{
  my($cdsseq)=@_;

  my $cdsnrseq = makename($cdsseq,"nr.cds");
  my $cmd="$APPfastanrdb $cdsseq > $cdsnrseq";
  unless(-s $cdsnrseq) {
  my $runerr= runcmd($cmd);
  }
  push @tmpfiles, $cdsnrseq;
  return($cdsnrseq);
}

  ## 1.1. reassign redundant cds from aaqual = pCDS to reduce utrbad/utrpoor set.
  ## .. need only change header ID in $cdsseqnr ; fastanrdb puts all redundant ids on header.

# %perfect_dups= redundant_idset($cdsseqnr); # read headers after rebest, for dupclass_fileset
sub redundant_idset
{
  my($cdsnrseq)=@_;
  my %pdup=();
  open(F,$cdsnrseq) or return %pdup; 
  while(<F>) { if(/^>/) { s/>//; my($d1,@dp)=split; map{ $pdup{$_}=$d1; }@dp; } } close(F);
  return %pdup;
}

sub nonredundant_reassignbest
{
  my($cdsnrseq, $aaqual)=@_;
  my (%aq, %better, ); 
  my ($nbetter,$nrec)=(0,0);
  my $flagfile= "$cdsnrseq.isbest";
  unless( -s $cdsnrseq and -s $aaqual ) {
    loggit(1,"ERR: nonredundant_reassignbest missing cdsnr:$cdsnrseq or aaqual:$aaqual"); 
    return($nbetter,$nrec);
  }
  return($nbetter,$nrec) if( -f $flagfile or $dryrun);
  
#   my $cdsnrhdr = makename($cdsnrseq,".hdr");
#   my $cmd="grep '^>' $cdsnrseq > $cdsnrhdr"; # dont really need .hdr file; combine below
  
  use constant PCDS_DIFF => 5; # what? should change any small diff?
  sub checkd{ my($d1,@d2)=@_; $nrec++; my($aw1,$aww1,$pc1,$aq1)=split",",$aq{$d1}; my $idbest=$d1; 
    foreach my $d (@d2) { my($aw,$aww,$pc,$aq)=split",",$aq{$d}; my $pdif=$pc - $pc1; 
    if($pdif >= PCDS_DIFF) { $pc1= $pc; $idbest=$d; } }
    if($idbest ne $d1) { $nbetter++; @d2= grep { $_ ne $idbest } @d2;
      $better{$d1}= join(" ",$idbest,$d1,@d2); }
  }
  
  runcmd("touch $flagfile"); push @tmpfiles, $flagfile; #??
  open(F,$aaqual);   while(<F>) { my($id,$aw,$gp,$aq,$tw)=split; $aq{$id}="$aw,$aq" if($aq); } close(F);
  open(F,$cdsnrseq); while(<F>) { if(/^>/) { s/>//; my @dp=split; checkd(@dp) if(@dp>1);} } close(F);
  if($nbetter>0) { # rewrite $cdsnrseq headers 
    my $ok= open(B,">$cdsnrseq.best");
    if($ok) { 
      open(F,$cdsnrseq); 
      while(<F>) { if(/^>(\S+)/) { if(my $hdr= $better{$1}) { s/>.*$/>$hdr/; } } print B $_; } 
      close(B); close(F); 
      my $cmd="mv $cdsnrseq $cdsnrseq.old; mv $cdsnrseq.best $cdsnrseq";
      system($cmd); #? runcmd($cmd);
      push @tmpfiles, "$cdsnrseq.old";
    }
  }
   
  return($nbetter,$nrec); # what?
}

=item 1.1 reassign redundant cds using aaqual  

  grep '^>' shrimt1trin1nr.cds > shrimt1trin1nr.cds.hdr
  grep '  ' shrimt1trin1nr.cds.hdr | cat shrimt1trin1.aa.qual - | perl -ne \
'if(/^>/) { s/>//; my @dp=split; checkd(@dp) if(@dp>1); } elsif(/^(\w+)\t/) { ($id,$aw,$gp,$aq,$tw)=split; 
$aq{$id}="$aw,$aq"; } END{ warn "done ncheck=$ncheck nbad=$nbad nbetter=$nbet\n"; } 
sub checkd{ my($d1,@d2)=@_; ($aw1,$aww1,$pc1,$aq1)=split",",$aq{$d1}; $ncheck++; $nbad++ if($pc1<1); 
foreach $d (@d2) { ($aw,$aww,$pc,$aq)=split",",$aq{$d}; $nbad++ if($pc<1); $pdif=$pc - $pc1; 
if($pdif>9) {  print "$d:$aw,$pc,$aq\tbetter\t$d1:$aw1,$pc1,$aq1\n"; $nbet++; } }}' 

done ncheck=15630 nbad= nbetter=8389  # all diff
done ncheck=15630 nbad= nbetter=2292  # pdif > 9

  must changes:  good vs utrbad
shrimt1trin1loc21564c0t5:1385,81%,complete      better  shrimt1trin1loc21564c0t2:1385,47%,complete-utrbad
shrimt1trin1loc9285c0t1:389,65%,complete        better  shrimt1trin1loc8886c0t1:389,33%,complete-utrbad
shrimt1trin1loc468087c0t1:42,50%,complete-utrpoor       better  shrimt1trin1loc229219c0t1:42,31%,complete-utrbad
shrimt1trin1loc43780c0t2:95,54%,complete-utrpoor        better  shrimt1trin1loc43780c0t1:95,37%,complete-utrbad

  maybe changes: same utr-qual but minor pCDS improvement
shrimt1trin1loc21426c0t4:555,83%,complete       better  shrimt1trin1loc21426c0t3:555,66%,complete
  = 15% pcds improve
shrimt1trin1loc65367c0t1:55,30%,complete-utrbad better  shrimt1trin1loc36046c0t1:55,7%,complete-utrbad
  = 23% pcds improve -- still bad but 30% >> 7% cds  
  
  ignore small changes?
shrimt1trin1loc58089c0t1:523,73%,complete       better  shrimt1trin1loc45529c0t1:523,70%,complete

  
=cut
    

  # 2. perfect fragment removal
sub nofragments_cds
{
  my($cdsseq)=@_;

  my $cdsnrseq = makename($cdsseq,"cd1.cds");
  my $cdlog = makename($cdsseq,"cd1.log");
  ## *** clusterize this .. -T NCPU works ok
  ## cd-hit-est -M 1500  -c 1.00 -d 0;  -l  length of throw_away_sequences, default 10
  ## .. redo opts so can use ENV replacements ..
  my $opts=" -c 1.00"; 
  $opts.=" -T $NCPU" if($NCPU>1);
  $opts.=" -M $MAXMEM" if($MAXMEM>1); # M in Megabytes
  $opts.=" -l ".($MINCDS-1) if($MINCDS>10); 
  if(my $cdopt=$ENV{CDHITESTOPT}) { $opts .= " ".$cdopt; }
  
## pogonus1all3nr.cds : cd-hit-est failing w/ bizzare mem request; maybe cds > 1.2G file too big?
## opt here to fasplit(2+), then cd-hit-est-2d -i pt1.cds -i2 pt2.cds, combine parts..

  my $cmd="$APPcdhitest $opts -d 0 -i $cdsseq -o $cdsnrseq 1> $cdlog 2>&1"; # log to ...
  unless(-s $cdsnrseq) {
  my $runerr= runcmd($cmd);
  unlink("$cdsnrseq.bak.clstr") if (-f "$cdsnrseq.bak.clstr");
  # push @erasefiles,"$cdsnrseq.bak.clstr";
  }
  # FIXME: grep '^>' $cdsnrseq > $cdsnrseq.hdr for later ref
  push @tmpfiles, $cdsnrseq, "$cdsnrseq.clstr", $cdlog ;
  return($cdsnrseq);
}


  # %perfect_frags= fragment_idset("$cdsseqnrcd1.clstr"); # read headers after rebest, for dupclass_fileset
sub fragment_idset
{
  my($cdhitclstr)=@_;
  my(%fragids,@fragids,$mainid);
  open(F, $cdhitclstr) or return %fragids;
  while(<F>) { # cd-hit clstr format
    if(/^>/) { if($mainid and @fragids) { map{ $fragids{$_}= $mainid } @fragids; } @fragids=(); $mainid=0; }
    elsif(/^(\d+)/) { my $i=$1;
      m/(\d+)(nt|aa), >(.+)\.\.\. (.+)$/;  
      my($tlen,$typ,$tid,$pinfo)=($1,$2,$3,$4); 
      my $ismain=($pinfo =~ /\*/)?1:0;  # FIXME: drop /$/; new merge fnum at end of line now;
      unless($tid) {
        # $nerr++;
      } elsif($ismain) {
        $mainid= $tid; ## $mainlen=$tlen; $pi=100; # push @cluster, "$tid,$pi";
      } else {
        push @fragids, $tid;
      }
    }
  } close(F);
  return %fragids;
}

sub aacluster
{
  my($aaseq)=@_;
  return () unless(-s $aaseq and $aaseq !~ /\.gz$/);
  
  my $aacdseq = makename($aaseq,"_cd90.aa");
  my $cdlog   = makename($aaseq,"_cd90.log");
  ## cd-hit  -M 1500  -c 0.9 -d 0;  -l  length of throw_away_sequences, default 10
  my $opts=""; 
  $opts.=" -T $NCPU" if($NCPU>1);
  $opts.=" -M $MAXMEM" if($MAXMEM>1); # M in Megabytes
  $opts.=" -l ".int($MINCDS/3) if($MINCDS>30); 

  my $cmd="$APPcdhit $opts -c 0.90 -d 0 -i $aaseq -o $aacdseq 1> $cdlog 2>&1"; # log to ...
  unless(-s $aacdseq) {
  my $runerr= runcmd($cmd);
  unlink("$aacdseq.bak.clstr") if (-f "$aacdseq.bak.clstr");
  # push @erasefiles,"$aacdseq.bak.clstr";
  }
  my $aaclstr= "$aacdseq.clstr";
  push @tmpfiles, $aacdseq, $aaclstr, $cdlog ;
  return($aacdseq,$aaclstr); ## want only "$aacdseq.clstr"
}


sub asmdupclass_sum
{
  my($OUTH,$outclass,$aasize)=@_;
  my($nt,%tab,%idclass);
  
  # open(OC,"cat $outclass | cut -f2,3 | sort | uniq -c |") or warn "ERR: asmdupclass_sum $outclass";
  unless( open(OC,$outclass) ) { loggit(1,"ERR: asmdupclass_sum missing: $outclass"); return -1; }

  print $OUTH  "# Class Table for $outclass \n";
  while(<OC>) {
    s/maybeok/okay/; s/altmidfrag/altmfrag/; # was amfrag
    my($id,$ac,$cl)=split; 
    for my $ct ($cl,"total") { $tab{$ct}{$ac}+=1; } $nt+=1;
    # my $fclass= ($ac =~ /drop/) ? "drop" : ($cl =~ /^alt/) ? "okalt" : "okay"; # is this bad default? problem was "amfrag" not altmidfrag!
    my $fclass= "none"; if($ac =~ /^drop/) { $fclass="drop"; } elsif($ac =~ /^okay/) { $fclass= ($cl =~ /^alt|^amfrag/)?"okalt":"okay"; }
    $idclass{$id}= $fclass; 
  } close(OC);
  
  my @ac=qw(okay drop); 
  printf $OUTH "%-9s\t","class"; print $OUTH join("\t",@ac,@ac)."\n";
  foreach my $cl (sort keys %tab) { 
    my @pv=(); my @nv=();
    foreach my $ac (@ac) { 
      my $n= $tab{$cl}{$ac}||0; push @nv,$n; 
      my $p= int(1000*$n/$nt)/10; push @pv,$p; 
    }
    if($cl eq "total"){ print $OUTH (("-") x 45); print $OUTH "\n"; } ;
    printf $OUTH "%-9s\t",$cl; print $OUTH join("\t",@pv,@nv)."\n"; 
  } 
  print $OUTH (("=") x 45); print $OUTH "\n";
  
  ##  aastat for top1k of okay idclass;  evigene/scripts/prot/aastat.sh
  if(-s $aasize) {
    my $top=1000; my($nok,$n,$sw,$sn,@aw,@nn);
    #?? add same stats for okalt, drop sets? : @$aw{$class} @$nn{class}, $nok{class}
    ##? count aaqual: complete+utrok+gapok vs partial/utrbad/gapbad 
    open(AA,"sort -k2,2nr $aasize |");
    while(<AA>) { next unless(/^\w/); my($id,$aw,$nn,$aqual)= split;
      if($idclass{$id} =~ /okay/) { $nok++; push @aw,$aw; push @nn,$nn; } 
    } close(AA);
    
    # stat top1k and allokay
    sub aastat { my ($name,$n)= @_; $n=@aw if($n > @aw);
      my($sw,$sn)=(0,0); for my $i (0..$n-1) { $sw+=$aw[$i]; $sn+=$nn[$i]; }
      my $aw=int($sw/$n); my $an=int(10*$sn/$n)/10;  my($mx,$md,$mi)= @aw[0,int($n/2),$n-1];
      print $OUTH  "$name\t n=$n; average=$aw; median=$md; min,max=$mi,$mx; sum=$sw; gaps=$sn,$an\n";     
    }
    print $OUTH  "# AA-quality for okay set of $aasize (no okalt): all and longest $top summary \n";
    aastat("okay.top",$top);
    aastat("okay.all",$nok) unless($nok<$top);
    # bug for cacao3all7.aa.qual: okay.all n=51551 when only n=31341 are -drop -alt in cacao3all7.trclass
    # problem was "^amfrag" vs ^altmidfrag above .. no dupids in aa.qual; cacao3all7.aa n=3397790;
  }
  
}


# @outfiles= asmdupclass_fileset($outclass,$cdnaseq,$aaseq,$cdsseq); 
# fixme: add hdr classinfo from drops: %perfect_dups, %perfect_frags, 
sub asmdupclass_fileset
{
  my($outclass,$cdnaseq,$aaseq,$cdsseq,$perfect_dups,$perfect_frags)=@_;
  my($nt,%idclass,%classinfo,@outfiles,%outfiles);
  
  # check/skip for existing fileset
  foreach my $inf ($cdnaseq,$aaseq,$cdsseq) {
    my($suf)= $inf =~ m,(\.[\.\w]+)$,;  $suf =~ s/\.gz//;
    foreach my $class (qw(okay okalt drop)) {  # 
      my $cname= makename($inf,".$class$suf");
      push @outfiles, $cname if($dryrun or -s $cname);  
    }
  }
  return @outfiles if(@outfiles > 3);

  if(ref $perfect_dups) { foreach my $pid (keys %$perfect_dups) { 
    my $mainid= $perfect_dups->{$pid}; $classinfo{$pid}= "perfectdup,drop,match:$mainid"; } }
  if(ref $perfect_frags) { foreach my $pid (keys %$perfect_frags) { 
    my $mainid= $perfect_frags->{$pid}; $classinfo{$pid}= "perfectfrag,drop,match:$mainid"; } }
      
  unless( open(OC,$outclass) ) { loggit(1,"ERR: asmdupclass_fileset missing: $outclass"); return -1; }
  while(<OC>) {
    next unless(/^\w/); 
    s/maybeok/okay/; # s/altmidfrag/amfrag/; 
    my($id,$ac,$cl,$bestid,$pia,$aaqual,$flags)=split; 
    my $fclass= ($ac =~ /drop/) ? "drop" : ($cl =~ /^alt/) ? "okalt" : "okay";
    $idclass{$id}= $fclass; $nt++;
    my $match= ($bestid eq $id)?"":",match:$bestid,pct:$pia";
    # fixme: match bestid == thisid ; ie no match for noclass...
    $classinfo{$id}= "$cl,$ac$match; aalen=$aaqual;";
    # .. above classinfo for perfectdups from 1.fastanrdb, 2. perfectfragments cdhitest
  } close(OC);
  if($nt < 2) { loggit(1,"ERR: asmdupclass_fileset classes=$nt"); return -1; }
  
  @outfiles= (); %outfiles= ();
  foreach my $inf ($cdnaseq,$aaseq,$cdsseq) {
    my($suf)= $inf =~ m,(\.[\.\w]+)$,;  $suf =~ s/\.gz//;
    foreach my $class (qw(okay okalt drop)) {
      my $cname= makename($inf,".$class$suf");
      if(-s $cname) { rename($cname,"$cname.old"); } # or what?
      my $chandle=undef;
      my $okw= open( $chandle, ">$cname");
      if($okw) {
        push @outfiles, $cname;
        $outfiles{"$inf.$class"}= $chandle;
        push @okayset, $cname unless($class eq "drop");
        push @dropset, $cname if($class eq "drop");
      } else {
        loggit(1,"ERR: writing $cname");
      }
    }
  }
  
  foreach my $inf ($cdnaseq,$aaseq,$cdsseq) {
    my $okin= 0; 
    if($inf =~ /\.gz$/) { $okin= open(IN,"gunzip -c $inf|"); }
    else { $okin= open(IN,$inf); }
    loggit(1,"ERR: asmdupclass_fileset reading: $inf") unless($okin); 
    my $chandle=undef;
    ## BUT insure no re-rev : UPDATE.maybe: revcomp(cdnaseq) if aaseq strand=- ??? or instead traa2cds -trout
    while(<IN>) {
      if(/^>(\S+)/) { my $id=$1; 
        my $class= $idclass{$id} || "drop"; # NOTE: 1,2 perfect dups are not in idclass; dropped already
        $chandle= $outfiles{"$inf.$class"}; # fail if missing?
        ## FIXME: want classinfo also for 1. perfdups, 2.perffrags, need ids from those files
        if(my $classinfo= $classinfo{$id}) { 
          $classinfo =~ s/ aalen=\S+// if(/aalen=/);
          s/evgclass=\S+//; s/$/ evgclass=$classinfo/; 
          }
        loggit(1,"ERR: missing output id=$id; $inf.$class") unless(ref $chandle);        
      }
      print $chandle $_ if(ref $chandle); ## bad chandle
    } close(IN);
  }

## add aa.qual per okayset, dropset .. pull from aaseq, not orig.aa.qual, add evgclass= info
## grep -v drop $pt.trclass | cut -f1,2 | sed 's/okay//' | ggrep -F -f - inputset/$pt.aa.qual > okayset/$pt.aa.qual

  foreach my $inf ($cdnaseq,$aaseq,$cdsseq) {
    foreach my $class (qw(okay okalt drop)) {
      my $ohand= $outfiles{"$inf.$class"};
      close($ohand) if(ref $ohand);
    }
  }
  
  return(@outfiles);
}



sub asmdupfilter_cds
{
  my($cdsblast,$aasize,$aacdhit)=@_;

  my $outclass= makename($aasize,".trclass");
  my $outaln  = makename($aasize,".alntab"); # change to .cds.alntab ?
  my $aflog   = makename($aasize,".adupfilt.log");
  ## my $aacdhit = makename($aasize,".aa.clstr");
  my $aaclstropt= ( $aacdhit and -s $aacdhit ) ? "-acdhit $aacdhit" : "";
  my $dbg=($debug)? "-debug" : "";

	## opt for cdsblast ==  lastz : $USE_LASTZ global, file=name-self97.lastz
	my $cdsblastop= ( $cdsblast =~ m/blast/ ) ? "-blastab $cdsblast"
			: ( $cdsblast =~ m/lastz/ or $USE_LASTZ ) ? "-lastz $cdsblast" : "-blastab $cdsblast";

	##  opt for -aablast=$aablast, table of trid refid bitscore identity align
	##  for homology scoring to preserve from drops..
	## FIXME: allow -ablastab to be .names table in asmrna_dupfilter2.pl, ie decide by file name
  my $aablastopt= ( $aablast and -s $aablast ) ? "-ablastab $aablast" : "";

  ##... simplest..
  # $evigene/scripts/rnaseq/asmrna_dupfilter2.pl -debug -CDSALIGN -tinyaln 35 
  #  -aasize $pt.aa.qual -blast slf95-$pt*blastn -outeq $pt.baln5 -outclass $pt.bclass5 > & log.dupfltb5.$pt
  ##... more inputs.. skip -dupids
  # $evigene/scripts/rnaseq/asmrna_dupfilter2.pl -debug -CDSALIGN -tinyaln 35 -aasize $pt.aa.qual 
  #  -acdhit $pt.aa.clstr -dupids $pt.dupids -blast outz/sd-slf95-$pt*blastn.gz 
  #  -outeq $pt.baln5c -outclass $pt.bclass5c
  
  my $cmd="$APPtrdupfilter $dbg -tinyaln 35 -aasize $aasize -CDSALIGN $cdsblastop $aablastopt $aaclstropt"
    ." -outeqtab $outaln -outclass $outclass >$aflog 2>&1";
    
  unless(-s $outclass) {
  my $runerr= runcmd($cmd);
  }
  push @tmpfiles, $outaln, $aflog; ## outclass is main output file?
  return($outclass,$outaln);
}
  
sub make_aaqual
{
  my($aaseq)=@_;
  my $aasize= makename($aaseq,".aa.qual"); 
  my $cmd="$APPaaqual $aaseq"; 
  unless(-s $aasize) {
  my $runerr= runcmd($cmd);
  }
  
  
  push @inputset, $aasize;
  return($aasize);
}


sub fadupids { 
  my $fa=shift; my($ok,$ndup)=(0,0); my %ids; 
  if($fa =~ /\.gz$/) { $ok= open(F,"gunzip -c $fa|"); }
  else { $ok= open(F,$fa); }
  while(<F>) { if(/^>(\S+)/) { my $id=$1; my $ni= ++$ids{$id}; 
    if($ni>1) { $ndup++; loggit(1,"ERR: dup id:$id in $fa"); } 
  } } close(F);
  # die "ERR: $ndup duplicate ids in $fa\n" if($ndup>0); # leave to caller ?
  return (wantarray) ? ($ndup,\%ids) : $ndup;
}

sub facount {
  my($infile)=@_;
  my $nrec=0;
  if(-s $infile) { 
    if($infile=~/\.gz$/) { $nrec=`gunzip -c $infile | grep -c '^>'`; }
    else { $nrec=`grep -c '^>' $infile`; }
    chomp($nrec); }
  return $nrec; # NOT: "nrec=".$nrec
}


sub fasize { 
  my $fa=shift; my $b=0; my $ok;
  #if($fa=~/\.gz$/) { $b=`gunzip -c $fa | grep -v '^>' | wc -c | sed 's/ .*//'`; } 
  #else { $b=`grep -v '^>' $fa | wc -c | sed 's/ .*//'`; }
  if($fa =~ /\.gz$/) { $ok= open(F,"gunzip -c $fa|"); }
  else { $ok= open(F,$fa); }
  while(<F>) { $b += length($_) unless(/^>/); } close(F);
  chomp($b); return $b; 
}

sub fasplit {
  my($fa, $spldir, $npart, $splsize)=@_;
  my @splist= (); my $ok= 0;
  my($atsize,$atfile)=(0,0);
  my $fabase= basename($fa);
  if($fa =~ /\.gz$/) { $ok= open(F,"gunzip -c $fa|"); }
  else { $ok= open(F,$fa); }
  unless($ok) { loggit(LOG_DIE,"ERR: fasplit $fa $spldir $npart $splsize"); return @splist; }
  mkdir($spldir) unless(-d $spldir);
  while(<F>) {
    if($atfile==0 or ($atsize > $splsize && /^>/)) {
      $atfile++; if($atfile>$npart) { } # finishup???
      close(SPL) if($atfile>1);
      # my $spname= "$spldir$fabase.split$atfile.fa";
      my $spname= $spldir . makename($fabase,".split$atfile.fa");
      $ok= open(SPL,'>',$spname);  $atsize= 0;
      unless($ok) { loggit(LOG_DIE,"ERR: fasplit $atfile $spname"); return @splist; }
      push @splist, $spname;
      }
    print SPL;
    $atsize+= length($_) unless(/^>/);
  } 
  close(SPL); close(F);
  return @splist;
}

  # FIXME: parallelize this for NCPU by splitting input cdnaseq to NCPU parts; cat parts> whole.
sub make_bestorf_ncpu
{
  my($npart,$cdnaseq,$cdsseq,$aaseq)=@_;
  return unless(-s $cdnaseq);
  
  my $MINAA=int($MINCDS/3);
  my $csize= fasize($cdnaseq);
  my $spldir= makename($cdnaseq,"_split/");
  my $splsize= 1 + int($csize/$npart);
  
  mkdir($spldir); # dryrun?
  my @splset= fasplit( $cdnaseq, $spldir, $npart, $splsize); 
  my (@cdsset,@aaset);
  my $icpu= 0; 
  foreach my $cdna1 (@splset) {
    my $aa1 = makename($cdna1,".aa");
    my $cds1= makename($cdna1,".cds");
    my $cmd1="$APPcdnabest -nostop -minaa=$MINAA -cdna $cdna1 -aaseq $aa1 -cdsseq $cds1";
    push @aaset, $aa1; push @cdsset, $cds1;
    my $pid= forkcmd($cmd1);    
    if(++$icpu > $npart) { while (wait() != -1) { }; $icpu= 0; }
  }
  ## THIS FAILED TO WAIT for pids ...
  ## not enough: wait(); # waitpid(@pids); ??
  while (wait() != -1) { };
  
  my $cmd;
  $cmd= "cat ".join(' ',@aaset)." > $aaseq";   runcmd($cmd);
  $cmd= "cat ".join(' ',@cdsset)." > $cdsseq"; runcmd($cmd);
  
  #?? are these useful parts to keep for later?
  if($debug||$dryrun) {
    push @erasefiles, @splset, @aaset, @cdsset; # which?
  } else {
    # runcmd("rm -r $spldir");
    foreach my $fn (@splset, @aaset, @cdsset) {  unlink($fn) if(-f $fn); } 
    rmdir($spldir); 
  }
  
  return($cdsseq,$aaseq);
}


sub make_bestorf
{
  my($cdnaseq)=@_;
  my $MINAA=int($MINCDS/3);
  my $aaseq = makename($cdnaseq,".aa");
  my $cdsseq= makename($cdnaseq,".cds");
  # FIXME: parallelize this for NCPU by splitting input cdnaseq to NCPU parts; cat parts> whole.
  my $cmd="$APPcdnabest -nostop -minaa=$MINAA -cdna $cdnaseq -aaseq $aaseq -cdsseq $cdsseq";
  unless(-s $cdsseq) {
    if($NCPU>1 and not $dryrun) { ($cdsseq,$aaseq)= make_bestorf_ncpu($NCPU,$cdnaseq,$cdsseq,$aaseq); }
    else { runcmd($cmd); }
  }
  return($cdsseq,$aaseq);
}


sub make_cdsFromTrAA
{
  my($cdnaseq,$aaseq)=@_;
  my $cdsseq= makename($cdnaseq,".cds");
  my $cmd="$APPtraa2cds -cdna $cdnaseq -aaseq $aaseq -cdsseq $cdsseq";
  unless(-s $cdsseq) {
  my $runerr= runcmd($cmd);
  }
  return($cdsseq);
}

sub get_bestorf
{
  my($cdnaseq,$aaseq,$cdsseq)=@_;
  $aaseq = makename($cdnaseq,".aa") unless($aaseq); # if(defined $aaseq and not $aaseq);
  $cdsseq= makename($cdnaseq,".cds") unless($cdsseq); # if(defined $cdsseq and not $cdsseq);
  if( $cdsseq and -s $cdsseq ) {
    # have already ..  
  } elsif( $aaseq and -s $aaseq ) {
    ($cdsseq) = make_cdsFromTrAA($cdnaseq,$aaseq);
  } else {
    ($cdsseq,$aaseq) = make_bestorf($cdnaseq);
  }
  push @inputset, $cdsseq,$aaseq ;
  return($cdsseq,$aaseq);
}


sub blastn_cds
{
  my($cdsseq)=@_;
  if($NCPU>1 and not $dryrun) { return blastn_cds_ncpu($NCPU,$cdsseq); }
 
  ## *** clusterize this .. NCPU works ok .. NOT!
  my $cdsdb= makename($cdsseq,"_db");
  my $cdsbltab= makename($cdsseq,"-self$CDSBLAST_IDENT.blastn");
  my $blog= makename($cdsdb,".log","xxxsuf");
  return($cdsbltab) unless(-s $cdsseq);

## 11mar.FIXME: need blastn -ungapped otherwise miss perfect match of parts ; e.g. alt-tr half=perfect, other=imperf
## 14mar : Very Sloooowwwww for daphmag, cluster not using NCPU effectively w/ -num_threads (e.g get 2x speedup for 32 cpu)
##   .. redo this w/ fasplit(query) and ncpu x 1-cpu parts.

  my $fmtcmd="$APPmakeblastdb -in $cdsseq -dbtype nucl -out $cdsdb -logfile $blog";
  my $opts=""; 
  $opts.=" -num_threads $NCPU" if($NCPU>1); # -num_threads isnt effective
  my $cmd="$APPblastn -task megablast $opts -ungapped -perc_identity $CDSBLAST_IDENT -evalue $CDSBLAST_EVALUE -dust no"
    ." -outfmt 7 -db $cdsdb -query $cdsseq -out $cdsbltab";
  unless(-s $cdsbltab) {
  runcmd($fmtcmd);
  my $runerr= runcmd($cmd);
  }
  push @erasefiles, "$cdsdb.nsq","$cdsdb.nin","$cdsdb.nhr";
  push @tmpfiles, $cdsbltab,$blog;
  return($cdsbltab);    
}

  # parallelize blastn by splitting query.fa to NCPU parts; blastn -num_threads not effective.
sub blastn_cds_ncpu
{
  my($npart,$cdsseq)=@_;
  
  my $cdsdb= makename($cdsseq,"_db");
  my $cdsbltab= makename($cdsseq,"-self$CDSBLAST_IDENT.blastn");
  my $blog= makename($cdsdb,".log","xxxsuf");
  return($cdsbltab) unless(-s $cdsseq);
  
  my $opts=""; 
  ## -dust no/yes ? does it make a diff?  -ungapped may need to be option..
  if(my $blopt=$ENV{BLASTNOPT}) { $opts .= " ".$blopt; }
  else { $opts.= " -ungapped -dust no "; } #  -evalue $CDSBLAST_EVALUE ??
  my $blcmd0="$APPblastn -task megablast -perc_identity $CDSBLAST_IDENT -evalue $CDSBLAST_EVALUE"
    ." $opts -outfmt 7 -db $cdsdb "; # add parts: -query $cdsseq -out $cdsbltab
  my $fmtcmd="$APPmakeblastdb -in $cdsseq -dbtype nucl -out $cdsdb -logfile $blog";

  unless(-s $cdsbltab) {
    runcmd($fmtcmd);
    push @erasefiles, "$cdsdb.nsq","$cdsdb.nin","$cdsdb.nhr";

    my $csize= fasize($cdsseq);
    my $spldir= makename($cdsseq,"_blsplit/");
    my $splsize= 1 + int($csize/$npart);

    mkdir($spldir); # dryrun?
    my @splset= fasplit( $cdsseq, $spldir, $npart, $splsize); 
    my (@bloset);
    my $icpu= 0; my $ipart=0;
    foreach my $cds1 (@splset) {
      $ipart++;
      my $cdsbltab1= makename($cds1,"-self$CDSBLAST_IDENT.blastn$ipart");
      my $cmd1= $blcmd0 . " -query $cds1 -out $cdsbltab1";  # add parts: -query $cdsseq -out $cdsbltab
      push @bloset, $cdsbltab1;
      my $pid= forkcmd($cmd1);    
      if(++$icpu > $npart) { while (wait() != -1) { }; $icpu= 0; }
    }
    while (wait() != -1) { };
    
    my $cmd= "cat ".join(' ',@bloset)." > $cdsbltab"; 
    my $runerr= runcmd($cmd); #??

    if($debug||$dryrun) {
      push @erasefiles, @splset, @bloset;  
    } else {
      foreach my $fn (@splset, @bloset) {  unlink($fn) if(-f $fn); } 
      rmdir($spldir); 
    }
    
  }

  push @tmpfiles, $cdsbltab,$blog;
  return($cdsbltab);    
}

sub lastz_cds
{
  my($cdsseq)=@_;
  if($NCPU>1 and not $dryrun) { return lastz_cds_ncpu($NCPU,$cdsseq); }
 
  my $cdsbltab= makename($cdsseq,"-self$CDSBLAST_IDENT.lastz");
  my $blog= makename($cdsbltab,".log","lastz");
  return($cdsbltab) unless(-s $cdsseq);

  my $opts=""; 
  if(my $blopt=$ENV{LASTZOPT}) { $opts = $blopt; }
	else { 
	$opts="--identity=$CDSBLAST_IDENT --coverage=20 --strand=plus --step=10 --seed=match12 --filter=nmatch:90 --gfextend --notransition --exact=20 --match=1,5 --nochain --nogapped --ambiguous=n";
	}
  ##my $cmd="$APPlastz \"$cdsseq\[multiple\]\" $cdsseq  $opts --format=general --output=$cdsbltab"; 
  my $cmd= $APPlastz .' "' . $cdsseq . '[multiple]" ' . $cdsseq . ' ' . $opts .' --format=general --output='.$cdsbltab; 

  unless(-s $cdsbltab) {
  my $runerr= runcmd($cmd);
  }
  push @tmpfiles, $cdsbltab, $blog;
  return($cdsbltab);    
}

sub lastz_cds_ncpu
{
  my($npart,$cdsseq)=@_;
  
  my $cdsbltab= makename($cdsseq,"-self$CDSBLAST_IDENT.lastz");
  my $blog= makename($cdsbltab,".log","lastz");
  return($cdsbltab) unless(-s $cdsseq);
  
 	## lastz opts tested; not sure how many are needed yet. PI,PA main ones.
 	## FIXME: ** before use, this needs more tests w/ other options; got wacky results e.g. self partial aligns, not full. 
	## PI=97 == CDSBLAST_IDENT; PA=20;
	## lseedopt="--step=10 --seed=match12"
	## lopt="--identity=$PI --coverage=$PA --filter=nmatch:90 --gfextend --notransition --exact=20 --match=1,5 --nochain --nogapped --ambiguous=n"
  ## $bbin/lastz "$trf[multiple]" $qfile $lseedopt $lopt --format=general --output=$onam.lzout  &
  ##
  ## ?? add --strand since we have stranded CDS-seq
  ## --strand=plus   search + strand only (matching strand of query spec)
  ## --coverage=$CDSLASTZ_COVER ??
  ## maybe not? --exact=20 --match=1,5 .. problems for snp,indel? 
  my $opts=""; 
  if(my $blopt=$ENV{LASTZOPT}) { $opts = $blopt; }
	else { 
	$opts="--identity=$CDSBLAST_IDENT --coverage=20 --strand=plus --step=10 --seed=match12 --filter=nmatch:90 --gfextend --notransition --exact=20 --match=1,5 --nochain --nogapped --ambiguous=n";
	}
	
	## lastz target[xxx] query  options : is recommended form, does option order mater?
  ##my $blcmd0="$APPlastz \"$cdsseq\[multiple\]\" QUERYHERE $opts --format=general "; # per part :  --output=$onam.lzout $qfile
  my $blcmd0= $APPlastz .' "' . $cdsseq . '[multiple]" QUERYHERE ' . $opts .' --format=general '; # per part :  --output=$onam.lzout $qfile

  unless(-s $cdsbltab) {
    my $csize= fasize($cdsseq);
    my $spldir= makename($cdsseq,"_blsplit/");
    my $splsize= 1 + int($csize/$npart);

    mkdir($spldir); # dryrun?
    my @splset= fasplit( $cdsseq, $spldir, $npart, $splsize); 
    my (@bloset);
    my $icpu= 0; my $ipart=0;
    foreach my $cds1 (@splset) {
      $ipart++;
      my $cdsbltab1= makename($cds1,"-self$CDSBLAST_IDENT.lastz$ipart");
      my $cmd1= $blcmd0 . " --output=$cdsbltab1";  $cmd1 =~ s/QUERYHERE/$cds1/;  
      push @bloset, $cdsbltab1;
      my $pid= forkcmd($cmd1);    
      if(++$icpu > $npart) { while (wait() != -1) { }; $icpu= 0; }
    }
    while (wait() != -1) { };
    
    my $cmd= "cat ".join(' ',@bloset)." > $cdsbltab"; 
    my $runerr= runcmd($cmd); #??

    if($debug||$dryrun) {
      push @erasefiles, @splset, @bloset;  
    } else {
      foreach my $fn (@splset, @bloset) {  unlink($fn) if(-f $fn); } 
      rmdir($spldir); 
    }
    
  }

  push @tmpfiles, $cdsbltab,$blog;
  return($cdsbltab);    
}


sub findapp
{
  my($aname)=@_;
  my $app="";
  $app=$ENV{uc($aname)} if(not $app and $ENV{uc($aname)});  
  $app=$ENV{$aname} if(not $app and $ENV{$aname});  
  $app=`which $aname` unless($app); 
  chomp($app);
  ## #tr2aacds: app=blastn, path=no blastn in 
  my $dol=0; if(not $app or $app =~ /^no $aname/) { 
    $app="echo MISSING_$aname"; $dol=($dryrun||$debug)? LOG_WARN : LOG_DIE; 
    }
  loggit( $dol, "app=$aname, path=$app");
  return($app);
}

sub findevigeneapp
{
  my($aname)=@_;
  my $app= $aname;
  # my $EVIGENES="$FindBin::Bin/.."; # ok?
  # my $APPcdnabest= findevigeneapp("$EVIGENES/cdna_bestorf.pl"); # allow ENV/path substitutions?
  $app="$EVIGENES/$aname" unless(-x $app);
  my $dol=0; 
  unless( -x $app) { 
    $app="echo MISSING_$aname"; 
    $dol=($dryrun||$debug)? LOG_WARN : LOG_DIE; 
    }
  loggit( $dol, "evigeneapp=$aname, path=$app");
  return($app);
}

sub runcmd
{
  my @cmd= @_;
  loggit( ($dryrun) ? 1 : 0,"CMD=",@cmd);  
  ## fail if $cmd[0] =~ /MISSING_/
  my $err= ($dryrun) ? 0 : system(@cmd);  
  if($err) { loggit(1,"ERR=$err ",$cmd[0]); } # ..
  return $err;
}

sub forkcmd
{
  my @cmd= @_;
  loggit( ($dryrun) ? 1 : 0,"forkCMD=",@cmd);  
  unless($dryrun) {
    my $pid= fork();
    if($pid) { # parent
      return $pid;
    } else { # child
      my $err= system(@cmd);
      exit($err);
    }
  }
  # if($err) { loggit(1,"ERR=$err ",$cmd[0]); } # ..
}

sub makename
{
  my($infile,$osuf,$insuf)=@_;
  $insuf ||= 'aa|blast|cds|cdna|tr|fasta|fa';  ## fixme need insuf: tr|fasta|fa
  my $outfile= $infile; $outfile =~ s/\.gz$//;  ## drop \.gz always 
  $outfile =~ s,\.($insuf)[^\/\s]*$,,; 
  $outfile.= $osuf if($osuf); 
  $outfile.= "_out" if($outfile eq $infile);
  return $outfile;
}


__END__

=item test1

  /bio/bio-grid/aabugs4/crusts/shrimpt/best1test
  log.tr2aacds1                     shrimt1trin1.drop.cds             shrimt1trin1nr.cds
  log.tr2aacds2                     shrimt1trin1.drop.tr              shrimt1trin1nrcd1.blastn
  log.tr2aacds3                     shrimt1trin1.okalt.aa             shrimt1trin1nrcd1.cds
  log.tr2aacds4                     shrimt1trin1.okalt.cds            shrimt1trin1nrcd1.cds.bak.clstr
  shrimt1trin1.aa                   shrimt1trin1.okalt.tr             shrimt1trin1nrcd1.cds.clstr
  shrimt1trin1.aa.qual              shrimt1trin1.okay.aa              shrimt1trin1nrcd1.log
  shrimt1trin1.adupfilt.log         shrimt1trin1.okay.cds             shrimt1trin1nrcd1_db.nhr
  shrimt1trin1.alntab               shrimt1trin1.okay.tr              shrimt1trin1nrcd1_db.nin
  shrimt1trin1.cds                  shrimt1trin1.tr.gz                shrimt1trin1nrcd1_db.nsq
  shrimt1trin1.drop.aa              shrimt1trin1.trclass

  >> cleanup
  dropset/               log.tr2aacds2          okayset/               tmpfiles/
  inputset/              log.tr2aacds3          shrimt1trin1.tr.gz
  log.tr2aacds1          log.tr2aacds4          shrimt1trin1.trclass
  
  OUTPUT seq:
    shrimt1trin1.{okay,okalt}.{tr,cds,aa} >> okayset/
    shrimt1trin1.drop.{tr,cds,aa} >> dropset/

  INPUT/made seq: >> inputset/
    shrimt1trin1.aa
    shrimt1trin1.aa.qual
    shrimt1trin1.cds
    shrimt1trin1nr.cds
    shrimt1trin1nrcd1.cds
      
  TEMP files to erase >> tmpfiles/ # (move to tmpfile subdir?)
    shrimt1trin1nrcd1_db.{nsq,nhr,nin} : drop
    shrimt1trin1nrcd1.blastn  : probably drop  
    shrimt1trin1nrcd1.cds.bak.clstr : always erase junk
    shrimt1trin1nrcd1.log : drop
    shrimt1trin1nrcd1.cds.clstr : probably drop; leave for inspection?
    shrimt1trin1.adupfilt.log : leave for inspection?
    shrimt1trin1.alntab : for inspect, reruns?
    
=cut

=item tr2aacds_qsub.sh

  #! /bin/bash
  ### env trset=allstrimpt1.tr datad=`pwd` qsub -q normal tr2aacds_qsub.sh
  #PBS -N tr2cds
  #PBS -l nodes=1:ppn=32,walltime=18:55:00
  #PBS -o tr2cds.$$.out
  #PBS -e tr2cds.$$.err
  #PBS -V

  ncpu=30
  maxmem=50000
  evigene=$HOME/bio/evigene/scripts
  #t2ac: app=cd-hit-est, path= echo MISSING_cd-hit-est
  #v45bad32kseq#export PATH=$HOME/bio/cdhit/bin:$PATH
  export PATH=$HOME/bio/cdhit461/bin:$PATH
  #t2ac: app=fastanrdb, path= echo MISSING_fastanrdb
  export fastanrdb=$HOME/bio/exonerate/bin/fastanrdb
  #t2ac: app=blastn, path= echo MISSING_blastn
  export PATH=$HOME/bio/ncbi2227/bin:$PATH

  if [ "X" = "X$trset" ]; then
    echo "missing env trset=xxxx.tr"; exit -1
  fi
  if [ "X" = "X$datad" ]; then
    echo "missing env datad=/path/to/data"; exit -1
  fi

  cd $datad/
  echo $evigene/prot/tr2aacds.pl -debug -NCPU $ncpu -MAXMEM $maxmem -log -cdna $trset
  $evigene/prot/tr2aacds.pl -debug -NCPU $ncpu -MAXMEM $maxmem -log -cdna $trset

=cut

=item blastn fasplit speedup

speedier, same trclass results  
  c: 142 min total, 110 min blastn, blastn-threads
  d:  46 min total,  22 min blastn, blastn-fasplit, 5+x faster
            
catfish1all4.c : using -num_threads 30
#t2ac: BEGIN with cdnaseq= catfish1all4.tr.gz date= Wed Mar 13 09:50:37 PDT 2013
# Class Table for catfish1all4.trclass 
class           okay    drop    okay    drop
althi           5.5     19.5    46419   163572
althi1          0.8     5.8     6947    49022
althia2         0       5.1     0       43011
altmfrag        0.2     0.5     2419    4317
altmfraga2      0       0       291     417
altmid          0.6     0.9     5188    7669
altmida2        0       0       799     809
main            4.3     9.8     36322   82723
maina2          0.2     0.2     1852    2450
noclass         2.1     25.9    18326   217631
noclassa2       0       0       37      191
parthi          0       13.4    0       112502
parthi1         0       2.1     0       17972
parthia2        0       1.9     0       16530
---------------------------------------------
total           14.1    85.8    118600  718816
=============================================
# AA-quality for okay set of catfish1all4.aa.qual (no okalt): all and longest 1000 summary 
okay.top	 n=1000; average=2205; median=1921; min,max=1476,13376; sum=2205718; gaps=5109,5.1
okay.all	 n=56537; average=273; median=125; min,max=40,13376; sum=15437828; gaps=130303,2.3
#t2ac: asmdupfilter_fileset= catfish1all4.okay.tr catfish1all4.okalt.tr catfish1all4.drop.tr catfish1all4.okay.aa catfish1all4.okalt.aa catfish1all4.drop.aa catfish1all4.okay.cds catfish1all4.okalt.cds catfish1all4.drop.cds
#t2ac: tidyup output folders: okayset dropset inputset tmpfiles
#t2ac: DONE at date= Wed Mar 13 12:12:14 PDT 2013
#t2ac: ======================================

catfish1all4.d : using blastn-fasplit(query,30)
#t2ac: BEGIN with cdnaseq= catfish1all4.tr.gz date= Thu Mar 14 12:00:19 PDT 2013
# Class Table for catfish1all4.trclass 
class           okay    drop    okay    drop
althi           5.5     19.5    46420   163576
althi1          0.8     5.8     6947    49019
althia2         0       5.1     0       43012
altmfrag        0.2     0.5     2419    4318
altmfraga2      0       0       291     417
altmid          0.6     0.9     5188    7669
altmida2        0       0       799     809
main            4.3     9.8     36321   82724
maina2          0.2     0.2     1852    2450
noclass         2.1     25.9    18326   217632
noclassa2       0       0       37      191
parthi          0       13.4    0       112498
parthi1         0       2.1     0       17971
parthia2        0       1.9     0       16530
---------------------------------------------
total           14.1    85.8    118600  718816
=============================================
# AA-quality for okay set of catfish1all4.aa.qual (no okalt): all and longest 1000 summary 
okay.top	 n=1000; average=2205; median=1921; min,max=1476,13376; sum=2205718; gaps=5109,5.1
okay.all	 n=56536; average=273; median=125; min,max=40,13376; sum=15437776; gaps=130303,2.3
#t2ac: asmdupfilter_fileset= catfish1all4.okay.tr catfish1all4.okalt.tr catfish1all4.drop.tr catfish1all4.okay.aa catfish1all4.okalt.aa catfish1all4.drop.aa catfish1all4.okay.cds catfish1all4.okalt.cds catfish1all4.drop.cds
#t2ac: tidyup output folders: okayset dropset inputset tmpfiles
#t2ac: DONE at date= Thu Mar 14 12:46:10 PDT 2013
#t2ac: ======================================


=cut

=item cd-hit-est MEM failure: V4.6.1 cures it

** dont know why this failed, other runs w/ more cds seq ok ..
>> CD-HIT version bug; update to V4.6.1 cures it; cause is CDS-size > 32k (some 56k real cds = titin)

Program: CD-HIT, V4.5.7 (+OpenMP), Jan 20 2013, 09:22:23
Command: /home/ux455375/bio/cdhit/bin/cd-hit-est -T 30 -M
         50000 -l 89 -c 1.00 -d 0 -i pogonus1all3nr.cds -o
         pogonus1all3nrcd1.cds

Started: Mon Mar  4 22:22:03 2013
================================================================
total seq: 1129597
longest and shortest : 56202 and 90
Total letters: 1250540301
Sequences have been sorted

Fatal Error:
not enough memory, please set -M option greater than 18446744060060

Program halted !!

Approximated minimal memory consumption:
Sequence        : 1400M
Buffer          : 30 X 18446744073201M = 18446744058472M
Table           : 2 X 34M = 69M
Miscellaneous   : 17M
Total           : 18446744059960M

=cut
