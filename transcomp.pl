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
our $EGAPP='transcomp.pl';  
our $EGLOG='trc';
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
  
my $APPtranseq=   findapp("transeq");
#my $APPselectORF= findapp("selectORF.ba");
#my $APPselectORF = findapp("TransDecoder.LongOrfs");
my $APPreformataaheader= findapp("reformat_aa_header.ba");

## .. these should call findevigeneapp() to warn/die if missing
my $APPcdnabest= findevigeneapp("cdna_bestorf.pl");
my $APPtraa2cds= findevigeneapp("prot/traa2cds.pl");
my $APPtrdupfilter= findevigeneapp("rnaseq/asmrna_dupfilter2.pl");  
#my $APPaaqual=    findevigeneapp("prot/aaqual.sh"); # replace w/ cdna_evigenesub:makeAaQual()
my $APPaaqual=    findapp("aaqualwTransd.sh"); 


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
  my $stub = $cdnaseq;
  $stub =~ s/\.f.*//g;
  #print "STUB = $stub\n";

  ($cdsseq,$aaseq) = get_bestorf($cdnaseq,$aaseq,$cdsseq);
  loggit(0, "bestorf_cds=",$cdsseq,"nrec=",facount($cdsseq));
  ($aasize)= make_aaqual($aaseq);
  my(%qualhash)= create_qualhash($aasize);
  
  print "1. nonredundant removal:  fastanrdb  input.cds\n"; 
  ($cdsseqnr)= nonredundant_cds($cdsseq);
  loggit(0, "nonredundant_cds=",$cdsseqnr,"nrec=",facount($cdsseqnr));

  print "1.1. reassign redundant cds from aaqual = pCDS to reduce utrbad/utrpoor set.\n";
  ## .. need only change header ID in $cdsseqnr ; fastanrdb puts all redundant ids on header.
  ## Note this rewrites $cdsseqnr file, to same name
  my(%perfduphash) = redundancy_resolve($cdsseqnr, "perfectdup", \%qualhash);
  #dump_hash(\%perfduphash);
  #NEED TO DO:loggit(0,"nonredundant_reassignbest=",$nbest,"of",$ndups); 

  print "2. perfect fragment removal  : clusterize : \$NCPU, \$MAXMEM\n";
  my $cdsnrids= {}; # hash ref now, or id file?
  ($cdsseqnrcd1)= nofragments_cds($cdsseqnr, "cd1");
  my (%fragduphash) = classify_cd($cdsseqnrcd1, \%qualhash, "perfectfrag");
  #loggit(0, "nofragments_cds=",$cdsseqnrcd1,"nrec=",facount($cdsseqnrcd1));

  #(%perfduphash) = reconcile_duplicates(\%perfduphash, \%fragduphash, $cdsseq);
  print_file_from_hash(\%perfduphash, 3, "dupclean", $cdsseq, $aaseq);

  print "3. BLAST\n";
  # 3. blastn alignments of hi-ident subsequences : clusterized  : $NCPU, $MAXMEM
  ($cdsblast)= blastn_cds($cdsseq);
  my(%althash)=blast2alt($cdsblast, \%perfduphash, \%qualhash);
  #dump_hash(\%althash);   

  print "4. cdhit -c 0.9 -i \$aaseq -o \$aaseqcd .. for aaseqcd.clstr input to dupfilter\n";
  #  .. filter uses aaclstr only to drop aa-hiident alts as uninformative; eg 1/2 of 13k althi.
  my($cdsseqnrcd90)= nofragments_cds($cdsseqnr, "cd90");
  my(%cdalthash) = classify_cd($cdsseqnrcd90, \%qualhash, "cd90");
  #dump_hash(\%cdalthash);   

  print "5. Classify - okay, okalt, drop\n";
  # 5.1 - merge alts from blast and cdhit
  %althash = combinealts(\%althash, \%cdalthash);
  #dump_hash(\%althash);

  # 5.2 - make okay, okalt files
  classify(\%althash,\%perfduphash);
  label_output_files("okalt",$cdsseq,$cdnaseq,$aaseq,\%althash, \%qualhash);
  #label_output_files("okay",$cdsseq,$cdnaseq,$aaseq,\%althash, \%qualhash);

  # 6. Tag files
  #make drop
  #make drop files
  #mark in files
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
    #if ($table eq 1) {
    #   ($cdsseq,$aaseq) = make_bestorf($cdnaseq);
    #} else {    
       ($cdsseq,$aaseq) = make_cdsFromORFfinder($cdnaseq, $cdsseq, $aaseq);
    #}
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
  
  #my $cmd="ORFfinder -in $cdnaseq -g $table -n true -out $cdsseq -outfmt 1";
  #runcmd($cmd);

  #bash to fix header for cdsseq is needed for aa offset
  #my $cmd2="$APPselectORF $cdnaseq $cdsseq $table";
  my $code;
  if ($table == "1") {$code = "universal";
  } elsif ($table == "10") {$code = "Euplotes";
  } elsif ($table == "6") {$code = "Tetrahymena";
  } elsif ($table == "2") {$code = "Mitochondrial-Vertebrates";
  } elsif ($table == "5") {$code = "Mitochondrial-Arthropods";
  } elsif ($table == "9") {$code = "Mitochondrial-Enchinoderms";
  } elsif ($table == "13") {$code = "Mitochondrial-Ascidians";
  } elsif ($table == "14") {$code = "Mitochondrial-Platyhelminths";
  } elsif ($table == "3") {$code = "Mitochondrial-Yeasts";
  } elsif ($table == "11") {$code = "Mitochondrial-Protozoans";
  }

  my $cmd2="TransDecoder.LongOrfs -t $cdnaseq -G $code";
  runcmd($cmd2);
  $cmd2 = "grep \"\>\" *transdecoder_dir/longest_orfs.cds \| sed 's/::/ #/g' \| sort -uk2,2 \> cds.list\;";
  runcmd($cmd2);
  #$cmd2 = "sed -i 's/\ #/::/g' cds.list\; subset_fasta.pl -i cds.list \< $cdnaseq.transdecoder_dir/longest_orfs.cds \> $cdsseq";
  #runcmd($cmd2);
  #$cmd2 = "sed -i -e 's/Gene.[0-9]*:://g' -e 's/::/ /g' $cdsseq";
  my $cmd = "sed -i 's/\ #/::/g' cds.list\; subset_fasta.pl -i cds.list \< $cdnaseq.transdecoder_dir/longest_orfs.cds \> $cdsseq\; sed -i -e 's/Gene.[0-9]*:://g' -e 's/::/ /g' $cdsseq\;";
  runcmd($cmd);
  $cmd = "subset_fasta.pl -i cds.list \< $cdnaseq.transdecoder_dir/longest_orfs.pep \> $aaseq\; sed -i -e 's/Gene.[0-9]*:://g' -e 's/::/ /g' $aaseq\;";
  runcmd($cmd);
 
  #($aaseq) = make_aaseq_from_transeq($cdsseq, $aaseq);

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
  $cmd = "best_aa.ba $cdsseq $aaseq $table";
  my $runerr= runcmd($cmd);

  #just renames the type - could do this in perl but leaving in case I need to change header
  $cmd2 = "cat $aaseq | $APPreformataaheader > tmp; mv tmp $aaseq";
  runcmd($cmd2);

  push @inputset, $aaseq;
  return($aaseq);
}

sub make_aaqual
{
  my($aaseq)=@_;
  
  my $aasize= makename($aaseq,".aa.qual"); 

  #Run prot/aaqual.sh
  #my $cmd="$APPaaqual $aaseq";
  my $cmd="$APPaaqual < $cdsseq > $aaseq.qual";

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
  print "#####################$cmd#################\n";
  unless(-s $cdsnrseq) {
  my $runerr= runcmd($cmd);
  }
  push @tmpfiles, $cdsnrseq;
  return($cdsnrseq);
}

sub create_qualhash
{
 my ($aaqual) = @_;

 my $data = "$aaqual";
 my %qualhash;

 open(my $fh, '<', $data) or die "Could not open '$data' $!\n";

 while (my $line = <$fh>) {
   chomp $line;
   my @fields = split "\t" , $line;

   my $key = shift @fields;
   my @rest = @fields;
   my $rests = join(",", @rest);
   $qualhash{$key} = $rests ;
 }

 return %qualhash;
}

sub redundancy_resolve
{
 my $file = $_[0];
 my $mode = $_[1];
 my %qualhash = %{$_[2]};

 my(%perfduphash, %duphash);

 my $tmp = makename($file,".tmp"); 

 open(O,'>',$tmp);
 open(F,$file); while(my $line = <F>) { 
	if($line =~ m/^>/) { 
		chomp $line;
		$line=~ s/>//; 
		my @duplist = split(/\  */,$line);

		foreach my $transcript_name (@duplist) {
		    	$duphash{$transcript_name} = $qualhash{$transcript_name};
		}

		my($best, $others);
	 	($best, $others) = sort_duphash(\%duphash, $mode);	

		$perfduphash{$best} = $others;

		%{duphash}= ();		

		#if ($mode eq "perfectdup") {
			print O ">$best\n";
  		#}
	} else {
		print O "$line";
	}
  }
  close(F);
  close(O); 

  my $cmd= "mv $file $file.$mode; mv $tmp $file";
  runcmd($cmd);

  return(%perfduphash)
}

sub sort_duphash 
{
 my %inhash = %{$_[0]};
 my $mode = $_[1];

 my %outhash;
 my $best;
 my $others;

 my @keys = keys %inhash;
 my $l = scalar @keys; 

 if ($l == "1") {
	if (%inhash) { #true if not empty
		$best = $keys[0];
		$best =~ s/\s+$//;
		$keys[0] = "no $mode";}
 } else {
	#@qual = split %duphash, ",";
	#need a way to sort...
 	#sort %inhash by aalen, aalen_wogap, type?, clen
   	#(see https://perlmaven.com/how-to-sort-a-hash-in-perl)

	$best = $keys[0]; #move to bottom for cd sorting ease - "longest is best!"!
	splice(@keys,0,1);

	for my $i (0 .. $#keys) {
		$keys[$i] =~ s/\s+$//;
		$keys[$i] ="$keys[$i]/$mode";
	}
 }

 $others = join("," , @keys);

 return($best, $others)
}


#----------------------------------------------------
#Subprocesses - Introduced in Step 2 - remove perfect frags
#----------------------------------------------------

sub nofragments_cds
{
  my($seqnr, $mode)=@_;
  
  my $cdlog = makename($seqnr,"$mode.log");
  my($opts, $app, $seqnrcd);

  ## *** clusterize this .. -T NCPU works ok
  ## cd1: cd-hit-est -M 1500  -c 1.00 -d 0;  -l  length of throw_away_sequences, default 10
  ## cd9: cdhit -c 0.9 -i $aaseq -o $aaseqcd .. for aaseqcd.clstr input to dupfilter ..
  if ($mode eq "cd1") {
	$app= $APPcdhitest;
	$opts=" -c 1.00"; 
  #} elsif ($mode eq "cd90") {
  } elsif (index($mode, "cd90") !=-1) {
 	$app= $APPcdhit;
 	$opts=" -c 0.9";
 	$seqnrcd = makename($seqnr,"$mode.aa");
  }

  $opts.=" -T $NCPU" if($NCPU>1);
  $opts.=" -M $MAXMEM" if($MAXMEM>1); # M in Megabytes
  $opts.=" -l ".($MINCDS-1) if($MINCDS>10); 

  if(my $cdopt=$ENV{CDHITESTOPT}) { $opts .= " ".$cdopt; }
  
  my $cmd="$app $opts -d 0 -i $seqnr -o $seqnrcd 1> $cdlog 2>&1"; # log to ...
  #print "\nours:\n$cmd\n";
  #print "/N/soft/rhel7/cd-hit/4.6.8/cd-hit-est  -c 1.00 -M 1000 -l 89 -d 0 -i test2/cdna2nr.cds -o test2/cdna2nrcd1.cds 1> test2/cdna2nrcd1.log 2>&1\n";
  #runcmd($cmd);
	
  unless(-s $seqnrcd) {
    my $runerr= runcmd($cmd); # $idops="update";
    unlink("$seqnrcd.bak.clstr") if (-f "$seqnrcd.bak.clstr");
  }
  
  return($seqnrcd);
}

sub classify_cd
{
 my $file = $_[0];
 my %qualhash = %{$_[1]};
 my $mode = $_[2];

 my(%duphash, %outhash, @cluster);
 my($best,$others,$tag, $out);

 $file = $file . ".clstr";
 open(F, $file); while(<F>) {
 	if (/Cluster/) {
		if(@cluster) {
			my @list;
			foreach my $nametag (@cluster) {
				my @set = split(" ", $nametag);

				if ($set[1] eq "best") { #assumes longest == best
					$best = $set[0];
					$best =~ s/\s+$//;
				} else {
					push @list, "$set[0]/$mode$set[1]";  
				}
			}

			if (@list) {
				$others = join("," , @list);
			} else {
				$others = "no $mode";
			}

			$outhash{$best} = $others;
			%duphash = (); @cluster = ();	#reset
		}
	} else {
                s/>//;
		s/\.\.\.//g;
		my @line = split;
		my $name = $line[2];
		chomp $name;

		if ( $line[3] eq "*") {
			$tag = "best";
		} else {
		        $tag = $line[4];
			$tag =~ s/%//;
			if ($mode eq "cd90") {
				if ($tag >= 98) { $tag = "-a2"; } else { $tag = "-default";}
			} else {
				$tag = "";
			}
		}
		$name = "$name $tag";				
		push @cluster, "$name\n";
	}
 }
 close(F);

 #my @keys= keys %outhash; foreach my $k (@keys) { print "k: $k e: $outhash{$k}\n";}

 return (%outhash);

}

sub reconcile_duplicates
{
 my %perfectdups = %{$_[0]};
 my %fragdups = %{$_[1]};
 my $cdsseq = $_[2];
 my $mode ="cleanupdups";

 my @keys = keys %fragdups;
 my @pvalues = values %perfectdups;
 my @remove;

 foreach my $key (@keys) {
	if ($fragdups{$key} eq "no perfectfrag") {
		#skip it
		#print "$fragdups{$key} skipped\n";
	} elsif (exists($perfectdups{$key})) { 
		#it is a frag of a best
		#print "\nfragdups key $key is a hit to a best: $key";

		if ($perfectdups{$key} eq "no perfectdup") { 
			#print " but no perfect dups are found\n";
			$perfectdups{$key} = $fragdups{$key}; 
		} else { 
			#print " so it is added to the perfect dups\n";
			$perfectdups{$key} = $perfectdups{$key} . "," . $fragdups{$key};
                }

		my @frags = split /,/, $fragdups{$key};
		foreach my $f (@frags) { 
			$f =~ s/\/.*//g; 
			#print "fragment $f removed from perdups as it is a fragment.\n"; 
			push @remove, $f 
		};
				
	} elsif ( "$key/perfectdup" ~~ @pvalues ) {
		#print "this is a fragment of a duplicate";
		#it is a fragment of a duplicate
	
		my ($pkey) = grep{ $perfectdups{$_} eq '$key/perfectdup'} keys %perfectdups;
		#print " and the fragments matching $fragdups{$key} have been added to $pkey.\n";

		$perfectdups{$pkey} = $perfectdups{$pkey} . "," . $fragdups{$key}; 
	}
 }

 #remove fragdup key from perfectdup keys, as it now merged to better hit
 foreach my $d (@remove) { delete $perfectdups{$d}; }

 return(%perfectdups); #MAIN -> all keys are either okay or okalt
}

sub dump_hash
{
  my %hash = %{$_[0]};

  my @keys = keys %hash;
  foreach my $k (@keys) { print "k: $k e: $hash{$k}\n";}
}

sub print_file_from_hash
{
  my %hash = %{$_[0]};
  my $which = $_[1];
  my $mode = $_[2];
  my $cdsseq = $_[3];
  my $aaseq = $_[4];
 
  my @keys = keys %hash;
  my $list = makename($cdsseq,".list");
  my $cdsout = makename($cdsseq, ".cds.$mode");
  my $aasout = makename($cdsseq, ".aa.$mode");

  my $cmd;

  open(O,'>',$list);

  foreach my $k (@keys) {
	print O ">$k\n";
  }
  close(O);

  if ($which == 1) { 
	$cmd = "subset_fasta.pl -i $list < $cdsseq > $cdsout; mv $cdsseq $cdsseq.pre$mode; mv $cdsout $cdsseq";
	runcmd($cmd);
  } elsif ($which == 2) {
	$cmd = "subset_fasta.pl -i $list < $aaseq > $aasout; mv $aaseq $aaseq.pre$mode; mv $aasout $aaseq";
	runcmd($cmd);
  } elsif ($which == 3) {
	$cmd = "subset_fasta.pl -i $list < $cdsseq > $cdsout; mv $cdsseq $cdsseq.pre$mode; mv $cdsout $cdsseq";
	runcmd($cmd);


	$cmd = "subset_fasta.pl -i $list < $aaseq > $aasout; mv $aaseq $aaseq.pre$mode; mv $aasout $aaseq";
	runcmd($cmd);
  }

  #no return, just prints files  
}

#----------------------------------------------------
#Subprocesses - Introduced in Step 3 - BLAST Alternates
#----------------------------------------------------

sub blastn_cds
{
  my($cdsseq)=@_;
  if($NCPU>1 and not $dryrun) { return blastn_cds_ncpu($NCPU,$cdsseq); }
 
  ## *** clusterize this .. NCPU works ok .. NOT!
  my $cdsdb= makename($cdsseq,"_db");
  my $cdsbltab= makename($cdsseq,"-self$CDSBLAST_IDENT.blastn");
  my $blog= makename($cdsdb,".log","xxxsuf");
  return($cdsbltab) unless(-s $cdsseq);

  my $fmtcmd="$APPmakeblastdb -in $cdsseq -dbtype nucl -out $cdsdb -logfile $blog";
  my $opts="";
  $opts.=" -num_threads $NCPU" if($NCPU>1); # -num_threads isnt effective
  my $cmd="$APPblastn -task megablast $opts -ungapped -perc_identity $CDSBLAST_IDENT -evalue $CDSBLAST_EVALUE -dust no" . ' -outfmt "6 std qlen slen" ' . "-db $cdsdb -query $cdsseq -out $cdsbltab";
  unless(-s $cdsbltab) {
	  runcmd($fmtcmd);
	  my $runerr= runcmd($cmd);
  }
  push @erasefiles, "$cdsdb.nsq","$cdsdb.nin","$cdsdb.nhr";
  push @tmpfiles, $cdsbltab,$blog;
  return($cdsbltab);
}

sub blast2alt
{
  my $cdsbltab = $_[0];
  my %perfduphash = %{$_[1]};
  my %qualhash = %{$_[2]};
  my(%althash, %duphash, %blasthash);

  #load file into hash
  open(O, $cdsbltab); while(my $line = <O>) {
      chomp $line;
      my @blastentry = split(/\t/,$line);
      my $query = shift @blastentry;

      #keep only non-self hits, id>=90, bitscore >= 0.7*qlen, alignment length >=200bp (~exon size)
      my $bitcutoff = 0.7 * $blastentry[11];
      my $alencutoff = 200;

      if ($blastentry[2] >= 90 && $blastentry[10] >= $bitcutoff && $blastentry[2] >= $alencutoff) {
	  if ($query eq $blastentry[0]) { #skip as it is self to self;
	  } else {
	      if ($blastentry[2]>=98) { $blastentry[13] = "/blast-hi1";
	      } elsif ($blastentry[2]>=93) { $blastentry[13] = "/blast-mid";
	      } else { $blastentry[13] = "/blast-default";}

	      if (exists($blasthash{$query}) && $blasthash{$query}) {
		   $blasthash{$query} = $blasthash{$query} . "," . $blastentry[0] . $blastentry[13];
  	      } else {
                   $blasthash{$query} = $blastentry[0] . $blastentry[13];
              }
          }
      }
  }
  close(O);

  #collapse %blasthash to $best: @others
  my ($best, $others);
  my @keys = keys %blasthash;
  foreach my $query (@keys) { 
        my (@list, @duplist);
	my $hits=$blasthash{$query};
	my @hit_list=split(",",$hits);	
	push @hit_list, $query;
	
	foreach my $name (@hit_list) {
		$name =~ s/\/blast-.*//;
		push @duplist, $name;
	}
	

        foreach my $transcript (@duplist) {
               $duphash{$transcript} = $qualhash{$transcript};
	}
	($best, $others) = sort_duphash(\%duphash, "blast");
	%duphash = ();

	if (exists($althash{$best})) { #skip;
	} else {
		$althash{$best} = $blasthash{$best};	
	}
  }

  return(%althash);
}

#----------------------------------------------------
#Subprocesses - Introduced in Step 5 - Alts Combined
#----------------------------------------------------

sub combinealts 
{
  my %althash = %{$_[0]};
  my %cdalthash = %{$_[1]};
  
  #update althash with details from cdalthash
  my @keysa = keys %althash;
  my @keysc = keys %cdalthash;

  foreach my $key (@keysa) {
      if (exists($cdalthash{$key}) && $cdalthash{$key} ne "no cd90") {
	  my @cdhits = split(",",$cdalthash{$key});
	  my @althits = split(",",$althash{$key});
	  
	  my @cdhit;
	 
	  foreach my $cdhitname (@cdhits) {
	     @cdhit = split("/", $cdhitname);

	     if (index($althash{$key}, $cdhit[0]) != -1) {
		my $find = "$cdhit[0]/";
		my $replace = $cdhit[0] . "/" . $cdhit[1] . "/";
	        $althash{$key} =~ s/$find/$replace/;
	     } else {
               $althash{$key} = $althash{$key} . "," . $cdalthash{$key};
	     }
	  }
      }
  }

  foreach my $key (@keysc) {
      if ($cdalthash{$key} ne "no cd90") {
          if(exists($althash{$key})) { #skip;

	  } else {
	      $althash{$key} = $cdalthash{$key};
          }

      }
  }

  return (%althash)
}

sub classify
{
  my %althash = %{$_[0]};  
  my %perfduphash = %{$_[1]};  
  my (@okay, @okalt, @main);

  #main set - this is the keys from the perfduphash, as perfduphash is main:dup
  @main = keys %perfduphash;  

  #okalt set - this is the entries from the althash, as althash is okay:okalt
  @okalt = values %althash;

  #okay set - all main that are not in alt
  for (@okalt) {
      s/\/.*//;
  }

  my %tmp = map { $_ => 1 } @okalt;


  foreach my $transcript (@main) {
      if (exists($tmp{$transcript})) {#skip;
      } else {
          push @okay, $transcript;
      }
  }
 
  #print scalar(@main);
  #print "\n";
  #print scalar(@okay);
  #print "\n";
  #print scalar(@okalt);
  #print "\n";

#make files

   my $list = makename($cdsseq,".list");
   my $cdsout = makename($cdsseq, ".cds.okay");
   my $aaout = makename($cdsseq, ".aa.okay");
   my $cdnaout = makename($cdsseq, ".fa.okay");

   my $lista = makename($cdsseq,".alt.list");
   my $cdsouta = makename($cdsseq, ".cds.okalt");
   my $aaouta = makename($cdsseq, ".aa.okalt");
   my $cdnaouta = makename($cdsseq, ".fa.okalt");


   open(O, ">", $list);
   foreach my $k (@okay) {
       print O ">$k\n";
   }
   close(O);

   open(O, ">", $lista);
   foreach my $k (@okalt) {
       print O ">$k\n";
   }
   close(O);
 
   my $cmd;

   $cmd = "subset_fasta.pl -i $list < $cdsseq > $cdsout"; 
   #print "$cmd\n";
   runcmd($cmd);

   $cmd = "subset_fasta.pl -i $list < $cdnaseq > $cdnaout"; 
   #print "$cmd\n";
   runcmd($cmd);

   $cmd = "subset_fasta.pl -i $list < $aaseq  > $aaout"; 
   #print "$cmd\n";
   runcmd($cmd);

   $cmd = "subset_fasta.pl -i $lista < $cdsseq > $cdsouta"; 
   #print "$cmd\n";
   runcmd($cmd);

   $cmd = "subset_fasta.pl -i $lista < $cdnaseq > $cdnaouta"; 
   #print "$cmd\n";
   runcmd($cmd);

   $cmd = "subset_fasta.pl -i $lista < $aaseq > $aaouta"; 
   #print "$cmd\n";
   runcmd($cmd);
}

sub label_output_files
{
  my $mode=$_[0]; #okay
  my $cdsout=$_[1]; #. "." . $mode; #x.cds.okay
  my $cdnaout=$_[2]; #. "." . $mode;  #x.cdna.okay
  my $aaout=$_[3]; #. "." . $mode;    #x.aa.okay
  my %althash = %{$_[4]};
  my %qualhash = %{$_[5]};

  my ($name,$class,$match,@aametric,$aalen,@okay,$offs,$clen,$strand);


  my $stub = $cdsout;
  $stub =~ s/.cds//; 
  my $blast = $stub . "-self98.blastn";
  print "$blast\n";

  if ($mode eq "okalt") {
      my @keys = keys %althash;
      foreach my $key (@keys) {
          #print "key: $key\n";
          #print "alts: $althash{$key}\n";

          my @alts = split(",", $althash{$key});
          foreach my $entry (@alts) {
              my @labels = split ("/", $entry);
              #print "name: $labels[0]  tag: $labels[1]$labels[2]\n";
	      $name = $labels[0];
	      $class = "$labels[1]$labels[2],okalt";
	      $match = "$key";

	      #grab aa metrics
	      @aametric = split(",",$qualhash{$name});
              $aalen = "$aametric[0],$aametric[3],$aametric[4]";
	      $clen = $aametric[5];

              $offs = "";
	      $strand = "+";

               #change header
	       #my $cmdcdna = ">$name evgclass=$class,match:$match,pct:\\\; aalen=$aalen\\\;";
	       #my $cmdcds = ">$name type=cds\\\; aalen=$aalen\\\; clen=$clen\\\; strand=$strand; offs=$offs\\\; evgclass=$class,match:$match,pct:\\\;";
	       #my $cmdaa = ">$name aalen=$aalen\\\; clen=$clen\\\; strand=$strand; offs=$offs; evgclass=$class,match:$match,pct:\\\;";

	       my $cmdmatch = ">$name evgclass=$class,match:$match,pct:";

	       my $cmd = "sed \-i \'s|>$name|$cmdmatch|g\' $stub.fa.$mode";
	       runcmd($cmd);

	       $cmd = "sed \-i \'s|>$name|$cmdmatch|g\' $stub.cds.$mode";
               runcmd($cmd);

	       $cmd = "sed \-i \'s|>$name|$cmdmatch|g\' $stub.aa.$mode";
               runcmd($cmd);

	       #make command to pull blast value - do after header made
	       $cmd = "P=`grep '$key.*$name' $blast | awk '\{printf \"%.2f%/%s\",\$4/\$13*100,\$3\}'`; sed -Ei \"s|(>$name.*pct:)(.*)|\\1\$P\\2|g\" $stub.fa.$mode";
               runcmd($cmd);
	       $cmd = "P=`grep '$key.*$name' $blast | awk '\{printf \"%.2f%/%s\",\$4/\$13*100,\$3\}'`; sed -Ei \"s|(>$name.*pct:)(.*)|\\1\$P\\2|g\" $stub.aa.$mode";
               runcmd($cmd);
	       $cmd = "P=`grep '$key.*$name' $blast | awk '\{printf \"%.2f%/%s\",\$4/\$13*100,\$3\}'`; sed -Ei \"s|(>$name.*pct:)(.*)|\\1\$P\\2|g\" $stub.cds.$mode";
               runcmd($cmd);
          }

       }
  } elsif ($mode eq "okay") {
      my $list = $stub . ".list";
      open(O, $list); while (my $line = <O>) {
           chomp $line; 
	   $line =~ s/>//;
           push @okay, $line;
      }

      foreach my $key (@okay) {
          $name = $key;
	  $class = "okay/main";

	  #grab aa metrics
	  @aametric = split(",",$qualhash{$name});
          $aalen = "$aametric[0],$aametric[3],$aametric[4]";

          #change header
	  my $cmdpart = ">$name evgclass=$class\\\; aalen=$aalen\\\;";
	  my $cmd = "sed -i 's|>$name|$cmdpart|g' $cdsout.$mode";
	  #print "$cmd\n";
          runcmd($cmd);
       }
  }
}

#make drop file
#add information
	#main - everything in @okay, @okalt
	#alt link - from %althash
	#duplicates - from %perduphash
	#no reading frame - ?
	#low quality - ?


