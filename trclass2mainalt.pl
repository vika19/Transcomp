#!/usr/bin/perl
# trclass2mainalt.pl

=item trclass2mainalt

  $evigene/scripts/prot/trclass2mainalt.pl -idpre Anofunz4iEVm  -trclass  evg2anofunz4h.tgclass3 
  cut from evigene/scripts/evgmrna2tsa2.pl  
  
=cut

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/.."); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname);
use cdna_evigenesub; #added # 

our $EVIGENES="$FindBin::Bin";  
our $dryrun=0; ## $DRYRUN ?
our $DEBUG= $ENV{debug}|| 0;

my $IDPREFIX= $ENV{idprefix} || 'EVGm'; 
my $SHOWDROPS=0; 
my ($trclass,$output,$logfile);  
my ($trpath,$trname, $sradatah, %settings);
my ($pubid_format,$altid_format,$GBPROID);
my $pubidnum_start=0;
my $CULLXEQ=0;

my $optok= GetOptions(
  "class|trclass=s", \$trclass,
  "dropshow!", \$SHOWDROPS,  
  "idprefix=s", \$IDPREFIX,  # FIXME: idpre option  overwritten by spppref
  "dryrun|n!", \$dryrun, 
  "CULLXEQ!", \$CULLXEQ,  
  "debug!", \$DEBUG, 
  );

die "EvidentialGene trclass2mainalt -trclass evigenes.trclass [-idprefix $IDPREFIX ] 
  makes tables of main-alt linkage and public ids, from tr2aacds output trclass table"
  unless($optok and $trclass);  

MAIN();

sub MAIN
{
  my($upstatus,$upfiles,$uptemp,$upokids)= (0) x 9; 

  # loggit(0, "BEGIN with input=",$trclass||$cdnaseq,"date=",`date`);
	# do_settings("restore",$trclass);  

	# ($cdnaseq,$trpath,$trname)= get_evgtrset($trclass,$cdnaseq,"publicset"); # subs.pm; cdnaseq => mrnaseq here
	#	loggit(0, "get_evgtrset=",$cdnaseq,$trpath,$trname); ## facount($cdnaseq)
	#	loggit(LOG_DIE, "Missing -mrna",$cdnaseq) unless($cdnaseq and -s $cdnaseq);

	($pubid_format,$altid_format)= make_IDPREFIX(); # default abbrev of $organism now

	my($maintab,$pubids,$nmaintr,$nalltr)= trclass2maintab($trclass,"publicset",$upokids);
	# loggit(0, "trclass2maintab primary n=",$nmaintr,"allntr=",$nalltr,$pubids); 

  # skip this?
  # evigene/scripts/rnaseq/asmrna_altreclass.pl -trclass $trclass -altrenum -out -debug
  # loggit(0, "DONE at date=",`date`);
}

# add for updated flags from genomap data, althi1.exonequal can be culled
# here maybe cull short, 1xon alts also
sub cullExonEq {
  my($md,$ads,$alt,$notes)=@_;
  my %culled; my $ncul=0;
  for my $td (@$ads) {
    my $cl= $alt->{$md}{$td};
    my $note= $notes->{$td};
    if($cl=~/^althi/ and not($note=~/refbest/) and $note=~/feq:([^;:\s]+)/ ) { 
      my $feq=$1;
      my @xe= grep /altmapxe/, split",",$feq;
      for my $xe (@xe) {
        my($xd,$xev)= $xe=~m,(\w+)/altmapxe(.+),;
        my $xcl= $alt->{$md}{$xd};
        if($xev>50 and not($xcl =~ /^cull/)) { #  and $xcl=~/^althi/
          $alt->{$md}{$td}="cull".$cl;
          $culled{$td}=$xe; $ncul++; last;
        }
      }
    }
  }
  return \%culled;
}

sub trclass2maintab
{
  my($trclass,$pubdir, $okids)=@_;
  my $ntr=0;  my $nerr=0;
  my $mainindex= $pubidnum_start;
  ## okids = \%validids after merge filesets
  my $hasokids= ($okids and ref($okids) and scalar(%$okids))?1:0;
  
  my $maintab = makename($trclass,".mainalt.tab","trclass");  # > $pt.mainalt.tab
  my $pubidtab= makename($trclass,".pubids","trclass");   
  if(not -f $pubidtab and $pubdir and -d $pubdir) {
  	my($pubd,$ft);
  	($pubd,$ft)= getFileset($pubdir,'pubids',$pubd);  $pubidtab=$ft if($ft);  
  	($pubd,$ft)= getFileset($pubdir,'mainalt.tab',$pubd);  $maintab=$ft if($ft);  
  }
  return($maintab,$pubidtab,$mainindex,$ntr) if( -s $maintab and -s $pubidtab);# or dryrun ..

  my($ok,%main,%mainsize,%alt,%altsize,%maindrops,%altdrops,%balt,%drop,$outh,$outpubidh,$inh);
  my(%aaqual,%piad,%notes);
  
  ($ok,$inh)= openRead($trclass);
  $ok= open($outh,'>',$maintab) if($ok);
  $ok= open($outpubidh,'>',$pubidtab) if($ok);
  unless($ok) { loggit(1,"ERR: parse $trclass TO $maintab"); return; }

  ## FIXME: only althi are reliably locus alternates; altmid .. are more likely paralogs
  while(<$inh>) {
    next unless(/^\w/); chomp;
    my($td,$ok,$cl,$md,$piad,$aq,$fl)=split "\t";
    
    ## new data bug, md == 0, piad == 0, from asmrna_dupfilter2b; keep? call these 'noclass' for now?
    # Anofunz4hEVm000033t23	okay	main	0	0	137,100%,partial	aaref:161,AGAP002273-PA,pflag:0
    # Anofunz4hEVm000033t26	okay	althi1	0	0	121,99%,partial	aaref:145,AGAP002273-PA,pflag:0,feq:Anofunz4hEVm000033t23/altpar48.0.0,Anofunz4hEVm000033t24/altpar35.0.0,Anofunz4hEVm000033t27/altpar50.0.0,Anofunz4hEVm000033t28/altpar50.0.0,Anofunz4hEVm000033t29/altpar100.0.0,Anofunz4hEVm000033t30/altpar50.0.0
    if($DEBUG and not $md) { $md="$td.miss"; $nerr++; }
    
    unless($cl and $md and $aq) { $nerr++; loggit(1,"ERR: trclass inline:$_"); next; } # what?? report?

		my $dropit=0;
		## maybeok needs handling, drop? or keep/mark/cull ??
		## should revise asm dupfliter to avoid these maybes, includes 'refbest', others some good/bad
		## now all are from exoneq flag, all 'althi1|part|frag' classes; 5702 refbest(dups?) of 13059
		## * maybe keep all maybes ; found some refbest/good being dropped here..
		## * but maybe is large subset, most althi1: 236768 drop, 195819 maybeok, 90660 okay in evg4anoalb.trclass 
    ## .. use hoscore if available, else keep til have hoscores?
		if($ok eq 'maybeok') { 
		   $ok="okay"; $cl.="maybe"; # if($fl=~/refbest/) { $ok="okay"; $cl.="maybe"; } 
		}
		
    if($ok ne 'okay') { $drop{$td}=1; $cl=$ok.$cl; $dropit=1;  } # next unless($SHOWDROPS);  OPTION: include drops?
    elsif($hasokids and not $okids->{$td}) { $drop{$td}=1;  $dropit=1; $cl='dropid'.$cl; } # next unless($SHOWDROPS);  unless $td =~ /utrorf/ ??

    ## 98/100/-sense/PitaEaR000975t24 << buggers.
    ## all piad entries should have pd and sense/asense fields, placeholder where needed '.' ?
    ## NEW SYNTAX tgclass for piad: 99/96/./altmap95/Anofunz4hEVm003223t1
    # need more sensible way to stick in genomap flags into trclass struct
    my $tgalt="";
    my($pi,$pa,$asense,$pd,@px)=split"/",$piad; # asense before pd ** NEED To revise asm dupfilter to clarify
    if($pd=~/^(alt|main|par|nocla)/) { $tgalt=$pd; $pd=(@px<1)?"":shift @px; } 
    if(@px) {
       # @px is what else?
    }
     
    if($asense =~ /sense/) { $pd="" unless($pd =~ /^\w/); }
    elsif($asense =~ /^\w/) { $pd=$asense; $asense=""; }
    $md=$pd if($pd =~ /^\w/);
   	$cl=~s/a2$//;  #? dont need 

    ## DROP alt intermediate problems: creates new nomain due to dropd mid link, need showdrops link hash ??
 		#OFF# if($dropit and not $SHOWDROPS) {
 			#OFF# next unless($cl =~ /main|noclass/); # this is enough, skip dropalts : FIXME, DONT skip here..
 		#OFF# }

 		$aaqual{$td}= $aq; # save for pubtable
 		$piad{$td}= $piad; # save for pubtable
 		$notes{$td}= $fl;
 		
 		my $aasize= ($aq =~ m/(\d+).(\d*)/) ? "$1.$2" : 0; # size.pCDS for best sort
    if($cl =~ /^main|^noclass/) { 
    	$main{$td}=$cl; $balt{$td}=$td;
    	$mainsize{$td}= $aasize; # ($aq =~ m/(\d+).(\d*)/) ? "$1.$2" : 0;
    } else { 
    	$alt{$md}{$td}= $cl; $balt{$td}=$md; 
    	$altsize{$td}= $aasize; 
    	# NOMAIN fix: always add dropmain here 
    }  
  }

  
  # if($SHOWDROPS) {
  #   # drops from dropset/*.aa headers for perfect_dups, perfect_frags info not in trclass ..
  #   my $ndr=0;
  #   my($dset,$droptr)= getFileset("$trpath/dropset",'drop.tr|drop.aa');  
  #   my($ok,$hin)= ($droptr) ? openRead($droptr) : (0,0);
  #   if($ok) { 
  #     while(<$hin>) { if(/^>(\S+)/) {  my $td=$1;
  #       my($cl,$ok1,$md)= m/evgclass=(\w+),(\w+),match:([^\s;,]+)/;
  #       $ok1="drop"; # ensure no bad cases
  #       if($cl and $md and not $balt{$td}) { $drop{$td}=1; $altdrops{$md}{$td}= $ok1.$cl; $ndr++; } # not: $balt{$td}=$md; 
  #     } 
  #   } close($hin); }
  #   loggit(0,"trclass2maintab: dropset adds $ndr"); 
  # }
	
  
  my %hasmain;
  my @amain= grep { not $main{$_} } sort keys %alt; # dropmain here now only for SHOWDROPS !
  foreach my $am (@amain) { 
    my $md= $balt{$am} || $am; ## $md=$am if($drop{$md});
  	if(!$main{$md} and $drop{$md}) { my $md1= $balt{$md}||""; if($md1 and $main{$md1}) { $md=$md1; } }
  	
    if($main{$md}) { my @at=keys %{$alt{$am}}; map{ $alt{$md}{$_}=$alt{$am}{$_}; $balt{$_}=$md; } @at; } 
    elsif($md) { 
     $main{$am}="NOMAIN";  # FIXME: get rid of these by finding alt main
     } 
  }
  
  foreach my $td (keys %balt) {
    my $md= $balt{$td} || $td; 
    $main{$md}="NOMAIN" unless($main{$md});# FIXME: get rid of these by finding alt main
  }


     
  ## add headers to these:
  #originalID     MainClass  Alternates
  #Public_mRNA_ID         originalID      PublicGeneID    AltNum
  ##FIXME: use extended realt format now?
  # #Public_mRNA_ID originalID      PublicGeneID    AltNum  Class   AAqual  pIdAln  Notes

  print $outh '#'.join("\t",qw(originalID MainClass Alternates))."\n";
  #oprint $outpubidh '#'.join("\t",qw(Public_mRNA_ID originalID PublicGeneID AltNum))."\n"
  print $outpubidh '#'.join("\t",qw(Public_mRNA_ID originalID PublicGeneID AltNum Class AAqual pIdAln Notes))."\n"
    if($outpubidh);

	my %doneid=();
  my @mainlist= sort{ $mainsize{$b} <=> $mainsize{$a} or $a cmp $b } keys %main;
  
  foreach my $md (@mainlist) { 

    my @ad= sort{$alt{$md}{$a} cmp $alt{$md}{$b}
      or $altsize{$b} <=> $altsize{$a} or $a cmp $b } keys %{$alt{$md}}; 

    my $culls= ($CULLXEQ) ? cullExonEq($md,\@ad,\%alt,\%notes) : {}; #?? here

    my $ad=join",",map{ "$_/".$alt{$md}{$_} } @ad; 
    my $mc=$main{$md}; 
    if($SHOWDROPS) {
    	my @add= sort{$altdrops{$md}{$a} cmp $altdrops{$md}{$b}} keys %{$altdrops{$md}};  
    	my $add=join",",map{ "$_/".$altdrops{$md}{$_} } @add; 
    	$ad .= ",$add" if ($add);
    }
    
    print $outh join("\t",$md,$mc,$ad)."\n";  # mainalt.tab
    
    if($outpubidh) { # should be required ??
      my $ialt= 0; my $needmain=0;
      my($cla,$aaq,$pida,$nots);
      $cla= $mc;  # cla=$main{$td}=$cl;  
      $aaq= $aaqual{$md}||"noaa";
      $pida=$piad{$md}||0;  
      $nots=$notes{$md}||"nonote";  
      if($mc eq "NOMAIN") { $cla=(@ad>0)?"main":"noclass"; } ## needs to change, to main? to noclass?
      elsif($mc =~ /^alt/) { } # is this were nomain show up? or @ad?

      if($drop{$md} or $doneid{$md}) { $needmain=1; }
      else {
      	$mainindex++; $needmain=0; # BUG: move below drop{}
      	my ($pubmrnaid,$pubgeneid)= make_pubid($md, $mainindex, ++$ialt);
      	#o# print $outpubidh join("\t",$pubmrnaid,$md,$pubgeneid,$ialt)."\n"; 
      	print $outpubidh join("\t",$pubmrnaid,$md,$pubgeneid,$ialt,$cla,$aaq,$pida,$nots)."\n";  #n
      	$ntr++; $doneid{$md}++;
      	}
      
      my @sad= sort{ $altsize{$b} <=> $altsize{$a} or $a cmp $b } @ad;      
      foreach my $ad (@sad) {
        unless($drop{$ad} or $doneid{$ad}) {
        $cla=  $alt{$md}{$ad}||"nocl"; 
        $aaq=  $aaqual{$ad}||"noaa";
        $pida= $piad{$ad}||0; # fixme == piad above
        $nots= $notes{$ad}||"nonote"; # fixme
        # $cull= $culls->{$ad}||""; # $cla == "cullalt..";
      	if($needmain) { 
      	  $mainindex++; $needmain=0; 
      	  if($cla=~/^alt/){ $cla=(@sad>1)?"main":"noclass"; } 
      	}  
        my ($altmrnaid,$altgeneid)= make_pubid($ad, $mainindex, ++$ialt);
        #o# print $outpubidh join("\t",$altmrnaid,$ad,$altgeneid,$ialt)."\n"; 
        print $outpubidh join("\t",$altmrnaid,$ad,$altgeneid,$ialt,$cla,$aaq,$pida,$nots)."\n"; 
        $ntr++; $doneid{$ad}++;
        }
      }
    }
  }
  close($inh); close($outh);

	# push @publicset, $maintab, $pubidtab; #?
  return($maintab,$pubidtab,$mainindex,$ntr);  # return main,alt id hashes ....
}


sub make_pubid
{
  my($oid, $mainindex, $altnum)= @_;
  my $pubidnum= $mainindex;
  $pubidnum_start= $pubidnum; #?
  my $pubgene = sprintf( $pubid_format, $pubidnum); 
  my $pubid   = $pubgene . sprintf( $altid_format, $altnum);
  return($pubid,$pubgene);
}

sub make_IDPREFIX
{
  my($existingID)= @_;
  my $digits=6; # or 7?
  my $ChangeDefaultIdPrefix= 0; # ($IDPREFIX eq $DEFAULT_SETTINGS{'IDPREFIX'}) ? 1:0;
  
  if($existingID and $existingID !~ m/^$IDPREFIX/) {
    $existingID=~ s/t\d+$//;
    my($prefix,$nums)= $existingID =~ m/^(\w+)(\d\d\d+)$/; ## prefix may have numbers, how many? 
    if($prefix) { $IDPREFIX= $prefix; $digits= length($nums) || 6; }
    $ChangeDefaultIdPrefix=0;
  }

  my $nd= ( $IDPREFIX =~ s/(0\d+)$// ) ? length($1) : $digits;
  my $pubid_format = $IDPREFIX.'%0'.$nd.'d';  
  my $altid_format = 't%d'; 
  # my $GBPROID= $IDPREFIX."_".$DATE;  
  return($pubid_format,$altid_format); # ,$GBPROID);
}
