#!/usr/bin/perl
# tr2aacds.pl

use constant VERSION => '2019.02.01';
use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...

use constant { kAAQUAL_MAX => 3, kAAQUAL_MIN => -3, kAAQUAL_NONE => 0, };
use constant { LOG_NOTE => 0, LOG_WARN => 1, LOG_DIE => -1, };
our $logh= undef;

use Getopt::Long;
use File::Basename qw(basename dirname);
#use cdna_evigenesub; # the use of this was not complete, I'm making a stand-alone for tr2aacds

# cdna_evigenesub globals:
our $EVIGENES=$ENV{EVIGENES} || "$FindBin::Bin/..";  
our $EGAPP='tr2aacds_ssbuild';  
our $EGLOG='t2ac';
our $dryrun=0;
our $DEBUG= $ENV{debug}|| 0;
our $EVGLOGH;

my $MINCDS = $ENV{MINCDS} || 90; # what? #maketraa3.sh: $AAMIN=40; $AMINPOO=100; $AMINBAD=200
my $CDSBLAST_IDENT= 98;
my $CDSBLAST_EVALUE= 1e-19;
my $NCPU= 1;
my $MAXMEM= 1000; # in Mb

our $okayclass  ='main|noclass';
our $okaltclass ='alt|part';
our $dropclass  ='drop|cull';

my ($debug,$tidyup,$runsteps,$USE_LASTZ)= (0) x 9;
my ($logfile,$aaseq,$aasize,$cdnaseq,$cdsseq,$aacdseq,$aaclstr,$aablast,$table)= (undef) x 20; 
my (@okayset,@dropset,@inputset,@tmpfiles,@erasefiles); # tidyup file sets

my $optok= GetOptions( 
  "cdnaseq|mrnaseq|trinput=s", \$cdnaseq, ## \@input,  # one only for this?
  "aaseq|aainput:s", \$aaseq,
  "cdsseq|cdsinput:s", \$cdsseq,
  "logfile:s", \$logfile,
  "ablastab=s", \$aablast,   # option traa-refaa.blastp.tall4 table for asmrna_dupfilter2.pl classifier
  "MINCDS=i", \$MINCDS,  
  "CDSBLAST_IDENT=i", \$CDSBLAST_IDENT, "CDSBLAST_EVALUE=i", \$CDSBLAST_EVALUE,  
  "NCPU=i", \$NCPU, "MAXMEM=i", \$MAXMEM,  
  "uselastz|lastz!", \$USE_LASTZ, # instead of blastn for cds-align
  "tidyup!", \$tidyup, 
  "dryrun|n!", \$dryrun,
  "debug!", \$DEBUG,
  "table:s", \$table,
);

$tidyup= 1 unless($dryrun||$debug);

openloggit( $logfile, $cdnaseq); 
loggit(1, "EvidentialGene tr2aacds.pl VERSION",VERSION);
loggit(1, "input cdnaseq=",$cdnaseq);
loggit(1, "ERR: unused arguments:",@ARGV) if(@ARGV>0); # do something w/ remaining @ARGV .. warn

my $APPblastn=    findapp("blastn"); 
my $APPmakeblastdb= findapp("makeblastdb");
#my $APPlastz="echo MISSING_lastz"; 
my $APPlastz= findapp("lastz") if ($USE_LASTZ || $debug); 
my $APPfastanrdb= findapp("fastanrdb");
my $APPcdhitest=  findapp("cd-hit-est"); 
my $APPcdhit=     findapp("cd-hit");
  
my $APPtranseq=   findevigeneapp("transeq");
my $APPselectORF= findevigeneapp("selectORF.ba");
my $APPreformataaheader= findevigeneapp("reformat_aa_header.ba");

## .. these should call findevigeneapp() to warn/die if missing
my $APPcdnabest= findevigeneapp("cdna_bestorf.pl");
my $APPtraa2cds= findevigeneapp("prot/traa2cds.pl");
my $APPtrdupfilter= findevigeneapp("rnaseq/asmrna_dupfilter2.pl");  
my $APPaaqual=    findevigeneapp("prot/aaqual.sh"); # replace w/ cdna_evigenesub:makeAaQual()

=item  tr2aacds pipeline algorithm

  I added a frame option - this is for use with ORFfinder and transeq.  This is the first step
  in allowing for flexibility in reading frame.

  0. make/collect input asm_name.{tr,aa,cds}, working mostly on .cds here

  1. perfect redundant removal:  fastanrdb  input.cds > input_nr.cds
    .. use aa.qual info for choosing top cds among identicals: pCDS only? ie need trsize along w/ cdssize

  2. perfect fragment removal: cd-hit-est -c 1.0 -l $MINCDS ..

  3. blastn, basic local align hi-ident subsequences for alternate tr.

  4. classify main/alternate cds, okay & drop subsets, using evigene/rnaseq/asmrna_dupfilter2.pl
     .. merges alignment table, protein-quality and identity, to score okay-main, ok-alt, and drop sets.

  5. make final output files from outclass: okay-main, okay-alts, drops 
      okayset is for public consumption, drops for data-overload enthusiasts (may contain valids).

=cut


#----------------------------------------------------
#MAIN
#----------------------------------------------------

sub MAIN_start {}
MAIN: {
  unless($ENV{TMPDIR}) {  ## FIXME sort TMPDIR, cant use /tmp
  require Cwd; if(my $cdir=Cwd::getcwd()) { $ENV{TMPDIR}=$cdir; } 
  if(my $cdir=FindBin::cwd2()) { $ENV{TMPDIR}=$cdir; } 
  }
 
  #check table used is valid
  my @table = (1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,31);
  if ( grep( /$table/, @table) ) {
     print "using ncbi amino acid table #$table\n";
  } else {
     print "please use a valid amino acid table:\n";
     print "1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,31\n";
     exit;
  }

  #make intermediate/output files
  my($cdsseqnr,$cdsseqnrcd1,$cdsblast,$outclass,$outaln);
  my(%perfect_dups,%perfect_frags);
  
  loggit(0, "BEGIN with cdnaseq=",$cdnaseq,"date=",`date`);


  # 0. make/collect input assembly_name.{tr,aa,cds}, we will work mostly with cds though

  #Required to have cdnaseq - may or may not have aaseq, cdsseq
  ($cdsseq,$aaseq) = get_bestorf($cdnaseq,$aaseq,$cdsseq);
  loggit(0, "bestorf_cds=",$cdsseq,"nrec=",facount($cdsseq));
  ($aasize)= make_aaqual($aaseq);

  # 1. nonredundant removal:  fastanrdb  input.cds 
  ($cdsseqnr)= nonredundant_cds($cdsseq);
  loggit(0, "nonredundant_cds=",$cdsseqnr,"nrec=",facount($cdsseqnr));

  ## 1.1. reassign redundant cds from aaqual = pCDS to reduce utrbad/utrpoor set.
  ## .. need only change header ID in $cdsseqnr ; fastanrdb puts all redundant ids on header.
  ## Note this rewrites $cdsseqnr file, to same name
  my($nbest,$ndups)= nonredundant_reassignbest($cdsseqnr,$aasize);
  loggit(0,"nonredundant_reassignbest=",$nbest,"of",$ndups); 

  }

#----------------------------------------------------
#Subprocesses - Set up
#----------------------------------------------------

sub findapp
{
  my($aname, $nodie)=@_;
  my $app=""; $nodie ||= 0;
  $app=$ENV{$aname} if(not $app and $ENV{$aname});  
  $app=$ENV{uc($aname)} if(not $app and $ENV{uc($aname)});  
  $app=`which $aname` unless($app); 
  chomp($app);

  my $dol=0; if(not $app or $app =~ /^no $aname/) { 
    $app="echo MISSING_$aname"; 
    $dol=($dryrun||$DEBUG||$nodie)? LOG_WARN : LOG_DIE; 
    }

  loggit( $dol, "app=$aname, path=$app");
  return($app);
}

sub findevigeneapp
{
  my($aname, $nodie)=@_;
  my $app= $aname; $nodie ||= 0;
  # my $EVIGENES="$FindBin::Bin/.."; # ok?
  # my $APPcdnabest= findevigeneapp("$EVIGENES/cdna_bestorf.pl"); # allow ENV/path substitutions?
  $app="$EVIGENES/$aname" unless(-x $app);
  my $dol=0; 
  unless( -x $app) { 
    $app="echo MISSING_$aname"; 
    $dol=($dryrun||$DEBUG||$nodie)? LOG_WARN : LOG_DIE; 
    }
  loggit( $dol, "evigeneapp=$aname, path=$app");
  return($app);
}

sub makename
{
  my($infile,$osuf,$insuf)=@_;
  $insuf ||= 'aa|blast|cdna|mrna|cds|tr|trclass|tbl|fasta|faa|fsa|fa';
  
  unless($infile) { warn "makename MISSING infile"; $infile="Noname$$"; }
  my $outfile= $infile; $outfile =~ s/\.gz$//; # bad for infile empty/undef .. do what?
  
  $outfile =~ s,\.($insuf)[^\/\s]*$,,; 
  $outfile.= $osuf if($osuf); 
  $outfile.= "_out" if($outfile eq $infile);
  
  return $outfile;
}

sub openloggit {
  my($logfile,$trname)= @_;
  
  if(not $logfile and defined $logfile) { # use output name
    $logfile= $trname || $EGLOG;
    $logfile= makename($logfile,".$EGAPP.log");  # need program suffix??
  }
  
  if($logfile) { open($logh, '>>', $logfile) or die $logfile; } 
}

sub loggit{ 
  my($dowarn,@msg)= @_; return unless($dowarn or @msg);
  my $s= join(' ',@msg); chomp($s); $s="FATAL $s" if($dowarn == LOG_DIE);
  if($logh){ print $logh "#$EGLOG: $s\n"; } elsif($dowarn>0||$DEBUG){ warn "#$EGLOG: $s\n"; }
  if($dowarn == LOG_DIE) { die "#$EGLOG: $s\n" ; }
}

sub openRead {
  my($fna)= @_; 
  my($ok,$hin)= (0,undef);
  $ok= ($fna =~ /\.gz$/) ? open($hin,"gunzip -c $fna|") : open($hin,$fna);  
  loggit(1,"ERR: openRead $fna") unless($ok);
  return ($ok,$hin);
}

sub runcmd
{
  my @cmd= @_;
  loggit( ($dryrun) ? 1 : 0,"CMD=",@cmd);  
  ## fail if $cmd[0] =~ /MISSING_/
  my $err= ($dryrun) ? 0 : system(@cmd);  ## ?? add option run: @result=`$cmd`;
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
}

#----------------------------------------------------
#Subprocesses - Fasta handling
#----------------------------------------------------

sub fasize
{ 
   #Reports the number of sequence characters in a file
   my $fa=shift; my $b=0; my $ok;

   if($fa =~ /\.gz$/) { $ok= open(F,"gunzip -c $fa|"); }
   else { $ok= open(F,$fa); }

   while(<F>) { $b += length($_) unless(/^>/); } close(F);
   chomp($b); return $b; 
}

sub facount 
{
  #Reports the number of sequence entries in a file
  my $fa=shift; my $n=0;
  my($ok,$hin)= openRead($fa);

  if($ok) { while(<$hin>) { $n++ if(/^>/); } close($hin); }
  return $n;  
}

sub fasplit
{
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

#----------------------------------------------------
#Subprocesses - Introduced in Step 0
#----------------------------------------------------

sub get_bestorf
{
  my($cdnaseq,$aaseq,$cdsseq)=@_;

  #create names for files if file not provided in main options
  $aaseq = makename($cdnaseq,".aa") unless($aaseq);
  $cdsseq= makename($cdnaseq,".cds") unless($cdsseq);

  # if cds exists - grab cds, make aa if not provided; else if only aa provided, make cds from AA, 
  # else if only cdna is provided, make best orf/aa from cdna
  
  if( $cdsseq and -s $cdsseq ) {
    if ( $aaseq and -s $aaseq ) {
       #do nothing, we have all three!
    } else {
       ($aaseq) = make_aaseq_from_transeq($cdsseq, $aaseq);
    }
  } elsif( $aaseq and -s $aaseq ) {
       #Alternate not needed here, as the cds is pulled from the cdna and offset in aa headers
       #However, need to reformat aa input from ORFfinder...
       ($cdsseq) = make_cdsFromTrAA($cdnaseq,$aaseq);
  } else {
    if ($table eq 1) {
       ($cdsseq,$aaseq) = make_bestorf($cdnaseq);
    } else {    
       ($cdsseq,$aaseq) = make_cdsFromORFfinder($cdnaseq, $cdsseq, $aaseq);
    }
  }

  push @inputset, $cdsseq,$aaseq ;
  return($cdsseq,$aaseq);
}

sub make_cdsFromTrAA
{
  my($cdnaseq,$aaseq)=@_;
  my $cdsseq= makename($aaseq,".cds");

  #Run prot/traa2cds script
  my $cmd="$APPtraa2cds -cdna $cdnaseq -aaseq $aaseq -cdsseq $cdsseq";
  unless(-s $cdsseq) {
  my $runerr= runcmd($cmd);
  }

  #returns cds file
  return($cdsseq);
}

sub make_cdsFromORFfinder
{
  my($cdnaseq, $cdsseq, $aaseq)=@_;
  
  my $cmd="ORFfinder -in $cdnaseq -g $table -n true -out $cdsseq -outfmt 1";
  runcmd($cmd);

  #bash to fix header for cdsseq is needed for aa offset
  my $cmd2="./$APPselectORF $cdnaseq $cdsseq $table";
  runcmd($cmd2);

  ($aaseq) = make_aaseq_from_transeq($cdsseq, $aaseq);

  return($cdsseq, $aaseq);
}

sub make_bestorf
{
  my($cdnaseq)=@_;

  #MINCDS is user defined or default 90BP, MINAA scales to 30AA
  my $MINAA=int($MINCDS/3);

  my $aaseq = makename($cdnaseq,".aa");
  my $cdsseq= makename($cdnaseq,".cds");
 
  unless(-s $cdsseq) {
    if($NCPU>1 and not $dryrun) {
         #Run in paralell if more than 1 CPU given 
         ($cdsseq,$aaseq)= make_bestorf_ncpu($NCPU,$cdnaseq,$cdsseq,$aaseq);
    }
    else { 
          #Run cdna_bestorf.pl
          my $cmd="$APPcdnabest -nostop -minaa=$MINAA -cdna $cdnaseq -aaseq $aaseq -cdsseq $cdsseq";
          runcmd($cmd); 
    }
  }
  return($cdsseq,$aaseq);
}

sub make_bestorf_ncpu
{
  #Parallel version of above script
  my($npart,$cdnaseq,$cdsseq,$aaseq)=@_;
  return unless(-s $cdnaseq);
  
  my $MINAA=int($MINCDS/3);
  my $csize= fasize($cdnaseq);
  my $spldir= makename($cdnaseq,"_split/");
  my $splsize= 1 + int($csize/$npart);
  
  mkdir($spldir);
  my @splset= fasplit( $cdnaseq, $spldir, $npart, $splsize); 
  my (@cdsset,@aaset);
  my $icpu= 0; 
  foreach my $cdna1 (@splset) {
    my $aa1  = makename($cdna1,".aa$icpu"); 
    my $cds1 = makename($cdna1,".cds$icpu"); 
    my $cmd1="$APPcdnabest -nostop -minaa=$MINAA -cdna $cdna1 -aaseq $aa1 -cdsseq $cds1";
    push @aaset, $aa1; push @cdsset, $cds1;
    my $pid= forkcmd($cmd1);    
    if(++$icpu > $npart) { while (wait() != -1) { }; $icpu= 0; }
  }

  while (wait() != -1) { };
  
  my $cmd;
  $cmd= "cat ".join(' ',@aaset)." > $aaseq"; runcmd($cmd);
  $cmd= "cat ".join(' ',@cdsset)." > $cdsseq"; runcmd($cmd);
  
  if($debug||$dryrun) {
    push @erasefiles, @splset, @aaset, @cdsset; # which?
  } else {
    foreach my $fn (@splset, @aaset, @cdsset) {  unlink($fn) if(-f $fn); } 
    rmdir($spldir); 
  }
  
  return($cdsseq,$aaseq);
}

sub make_aaseq_from_transeq
{
  my($cdsseq, $aaseq)=@_;

  my ($cmd, $cmd2);

  #assume frame 1, as the cds is being provided, allow for table definition
  #can't - partials!
  #$cmd = "$APPtranseq -sequence $cdsseq -outseq $aaseq -table $table -frame 1";
  $cmd = "./best_aa.ba $cdsseq $aaseq $table";
  my $runerr= runcmd($cmd);

  #just renames the type - could do this in perl, but leaving in case I need to change header
  $cmd2 = "cat $aaseq | ./$APPreformataaheader > tmp; mv tmp $aaseq";
  runcmd($cmd2);

  push @inputset, $aaseq;
  return($aaseq);
}

sub make_aaqual
{
  my($aaseq)=@_;
  
  my $aasize= makename($aaseq,".aa.qual"); 

  #Run prot/aaqual.sh
  my $cmd="$APPaaqual $aaseq";

  unless(-s $aasize) { my $runerr= runcmd($cmd); }
  push @inputset, $aasize;
  return($aasize);
}

#----------------------------------------------------
#Subprocesses - Introduced in Step 1
#----------------------------------------------------

sub nonredundant_cds
{
  #Runs Fastanrdb on sequences -any that have a match gets a second name in the header
  my($cdsseq)=@_;

  my $cdsnrseq = makename($cdsseq,"nr.cds");
  my $cmd="$APPfastanrdb $cdsseq > $cdsnrseq";
  unless(-s $cdsnrseq) {
  my $runerr= runcmd($cmd);
  }
  push @tmpfiles, $cdsnrseq;
  return($cdsnrseq);
}

#old version from tr2aacds.pl #doesn't work either.  Redesign this whole bit?
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

=item 1.1 new version of redundant clean up - not working and I don't know enough perl to fix...
sub nonredundant_reassignbest
{
  #determines which of the redundant copies is better
  my($cdsnrseq, $aaqual)=@_;
  my ( %better, );
  my ($nbetter,$nrec)=(0,0);
  my $flagfile= "$cdsnrseq.isbest";

  unless( -s $cdsnrseq and -s $aaqual ) {
    loggit(1,"ERR: nonredundant_reassignbest missing cdsnr:$cdsnrseq or aaqual:$aaqual"); 
    return($nbetter,$nrec);
  }

  return($nbetter,$nrec) if( -f $flagfile or $dryrun);

  #just let's you know it's getting this far    
  runcmd("touch $flagfile"); push @tmpfiles, $flagfile; 
  
  our $aaqualh= make_aaqual($aaqual);

  #Read in output from fastanrdb
  my($ok,$hin)= openRead($cdsnrseq); 

  #split headers, if more than one (meaning duplicate), must determine "best"

  #for each header, split into an array @dp
  #if more than one entry (duplicates)
  #  increase number of records, determin number better (nrcheckaaqual returns 1 for better, 0 for not)
  #  mv cdsnrseq to old, and best to cdsnrseq (main)
  
  #if number better is more than 0 (not best)
  #open new file for bests
  #for each header, print first entry in better hash, replace header in cdsnrseq
  #clean up files

  while(<$hin>) { if(/^>/){ s/>//; my @dp=split; 
     if(@dp>1) { $nrec++; $nbetter += nrcheckaaqual('nrdup', $aaqualh, \%better, @dp);} } 
  } close($hin);
  
  if($nbetter>0) { # rewrite $cdsnrseq headers 
    my $ok= open(B,">$cdsnrseq.best");
    if($ok) { 
      ($ok,$hin)= openRead($cdsnrseq); 
      while(<$hin>) { if(/^>(\S+)/) { print $better{$1}; if(my $hdr= $better{$1}) { s/>.*$/>$hdr/; } } print B $_; } 
      close(B); close($hin); 
      my $cmd="mv $cdsnrseq $cdsnrseq.old; mv $cdsnrseq.best $cdsnrseq";
      system($cmd); #? runcmd($cmd);
      push @tmpfiles, "$cdsnrseq.old";
    }
  }
   
  return($nbetter,$nrec);
}

use constant fPCDS_DIFF => 1; # what? should change any small diff?
use constant dPCDS_DIFF => 5; # what? should change any small diff?

sub nrcheckaaqual 
{

  my($nrtest, $aaqualh, $better, $dmain, @d2)=@_;  
  my $dofrag= ($nrtest =~ /frag/)?1:0;
  my $isbetter= 0;  
  print "$aaqualh\n"; #this is the file name, not the hash...
  my $aqmain= "";
  #my $aqmain= $aaqualh->{$dmain} or return 0;
  my($awm,$apm,$acm,$aqm)=split",",$aqmain;  
  my $okmain=($acm < kAAQUAL_MAX  or $dmain =~ /utrorf/)?0:1;
  my $idbest= $dmain; 
  my $awmMIN= 0.90 * $awm;
  return 0 if($dofrag and $okmain);
  
  ## sort @d2 by best aaqual ? using aaqualh acm scores ?
  foreach my $d2 (@d2) { 
    my $aqfrag= $aaqualh->{$d2} or next;
    my($awf,$apf,$acf,$aqf)=split",", $aqfrag;
    my $pdif= $apf - $apm; 
     
    if($dofrag) {
      my $okfrag= ($acf < kAAQUAL_MAX or $d2 =~ /utrorf/)?0:1;
      if($okfrag and $awf >= $awmMIN and $pdif >= fPCDS_DIFF) { 
        # delete $better->{$idbest}; $better->{$d2}= $dmain;
        $idbest= $d2; $apm= $apf;  $isbetter=1; 
        #not# return $isbetter; # ONLy return here if @d2 sorted best
      }
    } else { 
      ## cds perfect-duplicate, so aaquals are same, only utrs differ
      if($pdif >= dPCDS_DIFF) { $apm= $apf; $idbest= $d2; $isbetter=1; }  
    }
  }
  
  if($idbest ne $dmain) { 
    $isbetter=1; 
    if($dofrag) { delete $better->{$dmain}; $better->{$idbest}= $dmain; }
    else { @d2= grep { $_ ne $idbest } @d2;  $better->{$dmain}= join(" ",$idbest,$dmain,@d2);  }
  }

  return $isbetter;
}
=cut
#----------------------------------------------------
#Subprocesses - Introduced in Step 2
#----------------------------------------------------

#----------------------------------------------------
#Subprocesses - Introduced in Step 3
#----------------------------------------------------

#----------------------------------------------------
#Subprocesses - Introduced in Step 4
#----------------------------------------------------

#----------------------------------------------------
#Subprocesses - Introduced in Step 5
#----------------------------------------------------


