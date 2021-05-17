#!/usr/bin/perl
#
#  AMIGOS
#  Carlos M Duarte, Copyright 1998, 1999, 2003
#  A program to measure standard and non-standard torsion angles for 
#  nucleotides in standard PDB formatted files
#  Last modified 08/30/03;	version 0.3
#
$PDBSFX = ".pdb";
$pi = atan2(1,1) * 4;
$PDBSPRD = "rna_sprd.txt";
$PDBAREA = "rna_area.txt";
open (ALLOUT, ">>all_sprd.txt")|| die "coupldn't open all_sprd.txt\n";
printf ALLOUT"ResNo\tResID\tType\tETA\tTHETA\tCHI\tALPHA\tBETA\tGAMMA\tDELTA\tEPSILON\tZETA\tNUTWO\tTYPEA\n";
open (ALLAREA, ">>all_area.txt")|| die "couldn't open all_area.txt\n";
printf ALLAREA"ResNo\tResID\tType\tETA\tTHETA\tCHI\tALPHA\tBETA\tGAMMA\tDELTA\tEPSILON\tZETA\tNUTWO\tTYPEA\n";
#
#  The section below will open each file in a directory and then 
#  open output files in whatever directory you are in. There are 4 types
#  of output files
#  1) a large single file of all nts from all files
#  2) a large single file of all nts which fall within user spec'd range
#  3&4) for each RNA PDB file files there are corresponding output files
#	of the above type
#
print STDOUT "What directory are your structures file(s) in?\n";
$dir = <STDIN>;
chop ($dir);
print $dir;
opendir (THISDIR, $dir);
@allfiles = readdir THISDIR;
foreach $FILE (@allfiles) { 
        $_ = $FILE;
        if (/pdb$/ || /ent$/) {
                print "$_\n";
		$PDB = $dir . "/" . $FILE;
		print SPRDOUT "$FILE\n";
		$FILEAREA = $FILE."_area.txt";
		open (AREAOUT, ">>$FILEAREA")|| die "coupldn't open $FILEAREA \n";
		printf AREAOUT"ResID\tETA\tTHETA\n";
		print AREAOUT "$FILE\n";
		print ALLAREA "$FILE\n";
#
#  The subroutine call "&MEASURE" calls the subroutine that processes the file
#  It is followed by measurement file output sequences
#
		&MEASURE($PDB);
		print AREAOUT"This area has $AREACOUNT out of $RESCOUNT residues\n";
		$AREATOTAL += $AREACOUNT;
		$RESTOTAL += $RESCOUNT;
		$FILEOUT = $FILE."_sprd.txt";
		open (OUTFILE, ">>$FILEOUT")|| die "couldn't open $FILEOUT \n";
		printf OUTFILE"ResNo\tResID\tType\tETA\tTHETA\tCHI\tALPHA\tBETA\tGAMMA\tDELTA\tEPSILON\tZETA\tNUTWO\tTYPEA\n";
		for ($i = 1; $i <= $NUMRES; $i++) {
			if ($CPCP[$i] != 0 && $PCPC[$i] != 0) {
				printf OUTFILE"%-5d\t%-7s\t%-4s\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f%\t%6.1f\t%6.1f\t%3d\n", $i, $RESNAM[$i], $RESTYP[$i], $CPCP[$i], $PCPC[$i], $CHI[$i], $ALPHA[$i], $BETA[$i], $GAMMA[$i], $DELTA[$i], $EPSILON[$i], $ZETA[$i], $NUTWO[$i], $AHLX[$i];
				printf ALLOUT"%-5d\t%-7s\t%-4s\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f%\t%6.1f\t%3d\n", $i, $RESNAM[$i], $RESTYP[$i], $CPCP[$i], $PCPC[$i], $CHI[$i], $ALPHA[$i], $BETA[$i], $GAMMA[$i], $DELTA[$i], $EPSILON[$i], $ZETA[$i], $NUTWO[$i], $AHLX[$i];
				}
        		}
                }
	}
print ALLAREA"This area has $AREATOTAL out of $RESTOTAL residues\n";
#
#  Start of the processing subroutine
#
#  File is parsed according to format specified by the PDB (see PDB website)
#
sub MEASURE {
open (PDBFILE, $PDB)|| die "Couldn't open $PDB\n";
@ATOM = 0;
@RNA = 0;
$ATNUM = 0;
$NUMRES = 0;
$RESCOUNT = 0;
$NUMCHN = 0;
$AREACOUNT = 0;
while (<PDBFILE>) {
	last if (/^END/);
	if (/^ATOM/) {
		chop $_;
		$FOO = substr($_,12,4);
		$BAR = substr($_,21,6);
                if ($FOO eq $ATNAM && $BAR eq $ATRES[$ATNUM]) {
		  $INFORM = substr($_,0,27);
		  print "$INFORM skipped - duplicate atom model\n";
		  next;
		}
		$ATNUM++;
		$ATOM[$ATNUM] = $_;
		$ATNAM = substr($_, 12, 4);
		$ATRES[$ATNUM] = substr($_, 21, 5);
		$CHTM[$ATNUM] = substr($_, 21, 1);
		$ATRTP[$ATNUM] = substr($_, 17, 3);
		if ($ATRES[$ATNUM] ne $ATRES[$ATNUM-1]) {
			$NUMRES++;
			$RESCOUNT++;
			$RESNAM[$NUMRES] = $ATRES[$ATNUM];
			$RESSTRT[$NUMRES] = $ATNUM;
			$RESTYP[$NUMRES] = $ATRTP[$ATNUM];
			}
#
#  The following section assigns each atom type to arrays that will later 
#  be used to run measurements
#
#  examples of acceptable ribose atom names: C4'  -or-  C4*
#
		if ($ATNAM eq " P  ") {
			$PIDX[$NUMRES] = $ATNUM;
			}
		if ($ATNAM eq " C4'" || $ATNAM eq " C4*") {
			$C4IDX[$NUMRES] = $ATNUM;
			}
		if ($ATNAM eq " O4'" || $ATNAM eq " O4*") {
			$O4IDX[$NUMRES] = $ATNUM;
			}
		if ($ATNAM eq " C1'" || $ATNAM eq " C1*") {
			$C1IDX[$NUMRES] = $ATNUM;
			}
		if ($ATNAM eq " N1 " || $ATNAM eq " N1 ") {
			$N1IDX[$NUMRES] = $ATNUM;
			}
		if ($ATNAM eq " N9 " || $ATNAM eq " N9 ") {
			$N9IDX[$NUMRES] = $ATNUM;
			}
		if ($ATNAM eq " C2 " || $ATNAM eq " C2 ") {
			$C2IDX[$NUMRES] = $ATNUM;
			}
		if ($ATNAM eq " C4 " || $ATNAM eq " C4 ") {
			$C42IDX[$NUMRES] = $ATNUM;
			}
		if ($ATNAM eq " O3'" || $ATNAM eq " O3*") {
			$O3IDX[$NUMRES] = $ATNUM;
			}
		if ($ATNAM eq " O5'" || $ATNAM eq " O5*") {
			$O5IDX[$NUMRES] = $ATNUM;
			}
		if ($ATNAM eq " C5'" || $ATNAM eq " C5*") {
			$C5IDX[$NUMRES] = $ATNUM;
			}
		if ($ATNAM eq " C3'" || $ATNAM eq " C3*") {
			$C3IDX[$NUMRES] = $ATNUM;
			}
		if ($ATNAM eq " C2'" || $ATNAM eq " C2*") {
			$C21IDX[$NUMRES] = $ATNUM;
			}
		if ($ATNAM eq " O2'" || $ATNAM eq " O2*") {
			$RNA[$NUMRES] = 1;
			}
#
# The chain designation is based on the PDB file. Each residue has a 
#  corresponding residue number. If the first character is a letter
#  the start and end of a new chain are considered to be where the letter
#  changes (e.g. A15 followed by B1). The first and last nucleotide of a
#  chain are not measured
#
		if ($NUMCHN == 0) {
			$NUMCHN++;
			$STRTATM[$NUMCHN] = $ATNUM;
			$CHAIN[$NUMCHN] = $CHTM;
			$STRTRES[$NUMCHN] = $NUMRES;
			}
		elsif ($CHTM[$ATNUM] ne $CHTM[$ATNUM-1]) {
			$NUMCHN++;
			$STRTATM[$NUMCHN] = $ATNUM;
			$ENDATM[$NUMCHN-1] = $ATNUM-1;
			$STRTRES[$NUMCHN] = $NUMRES;
			$ENDRES[$NUMCHN-1] = $NUMRES-1;
			$CHAIN[$NUMCHN] = $CHTM;
			}
		}
	}
$ENDATM[$NUMCHN] = $ATNUM;
$ENDRES[$NUMCHN] = $NUMRES;
#
# Following section does the actual parsing of atoms for each individual 
# nucleotide measurement 
#
for ($i = 1; $i <= $NUMCHN; $i++) {
	print "chain number is $i\n";
	$RESCOUNT -= 2;
	$TRES = $ENDRES[$i] - $STRTRES[$i] + 1;
	for ($j = 1; $j <= $TRES; $j++) {
		$NRES = $j + $STRTRES[$i] - 1;
		if ($j == 1 || $j == $TRES || $RNA[$NRES] == 0 || $RNA[$NRES-1] == 0 || $RNA[$NRES+1] == 0) {
			$ALPHA[$NRES] = $BETA[$NRES] = $GAMMA[$NRES] = $DELTA[$NRES] = $EPSILON[$NRES] = $ZETA[$NRES] = $PCPC[$NRES] = $CPCP[$NRES] = $PCPPLATE[$NRES] = $CPCPLATE[$NRES] = $NUTWO[$NRES] = 0;
			}
		else {
			$P1NAM = $PIDX[$NRES];
			$C1NAM = $C4IDX[$NRES];
			$O4NAM = $O4IDX[$NRES];
			$CC1NAM = $C1IDX[$NRES];
			$O5NAM = $O5IDX[$NRES];
			$C5NAM = $C5IDX[$NRES];
			$C3NAM = $C3IDX[$NRES];
			$O3NAM = $O3IDX[$NRES];
			$C21NAM = $C21IDX[$NRES];
			@P1 = (substr($ATOM[$P1NAM],30,8),substr($ATOM[$P1NAM],38,8),substr($ATOM[$P1NAM],46,8));
			@C1 = (substr($ATOM[$C1NAM],30,8),substr($ATOM[$C1NAM],38,8),substr($ATOM[$C1NAM],46,8));
			@O4 = (substr($ATOM[$O4NAM],30,8),substr($ATOM[$O4NAM],38,8),substr($ATOM[$O4NAM],46,8));
			@CC1 = (substr($ATOM[$CC1NAM],30,8),substr($ATOM[$CC1NAM],38,8),substr($ATOM[$CC1NAM],46,8));
			@O5 = (substr($ATOM[$O5NAM],30,8),substr($ATOM[$O5NAM],38,8),substr($ATOM[$O5NAM],46,8));
			@C5 = (substr($ATOM[$C5NAM],30,8),substr($ATOM[$C5NAM],38,8),substr($ATOM[$C5NAM],46,8));
			@C3 = (substr($ATOM[$C3NAM],30,8),substr($ATOM[$C3NAM],38,8),substr($ATOM[$C3NAM],46,8));
			@O3 = (substr($ATOM[$O3NAM],30,8),substr($ATOM[$O3NAM],38,8),substr($ATOM[$O3NAM],46,8));
			@C21 = (substr($ATOM[$C21NAM],30,8),substr($ATOM[$C21NAM],38,8),substr($ATOM[$C21NAM],46,8));
			&MEAS(@P1,@C1);
			$DISPC[$NRES] = $DIST;
			&DIHED(@P1, @O5, @C5, @C1, 1);
			$BETA[$NRES] = $TOR;
			&DIHED(@O5, @C5, @C1, @C3, 2);
			$GAMMA[$NRES] = $TOR;
			&DIHED(@C5, @C1, @C3, @O3, 3);
			$DELTA[$NRES] = $TOR;
			&DIHED(@CC1, @C21, @C3, @C1, 3);
			$NUTWO[$NRES] = $TOR;
			if ($RESTYP[$NRES] eq "  G" || $RESTYP[$NRES] eq "  A"){
				$N9NAM = $N9IDX[$NRES];
				$CC4NAM = $C42IDX[$NRES];
				@N9 = (substr($ATOM[$N9NAM],30,8),substr($ATOM[$N9NAM],38,8),substr($ATOM[$N9NAM],46,8));
				@CC4 = (substr($ATOM[$CC4NAM],30,8),substr($ATOM[$CC4NAM],38,8),substr($ATOM[$CC4NAM],46,8));
				&DIHED(@O4, @CC1, @N9, @CC4, 4);
				$CHI[$NRES] = $TOR;
				}
			elsif ($RESTYP[$NRES] eq "  U" || $RESTYP[$NRES] eq "  C" || $RESTYP[$NRES] eq "  T"){
				$N1NAM = $N1IDX[$NRES];
				$CC2NAM = $C2IDX[$NRES];
	                        @N1 = (substr($ATOM[$N1NAM],30,8),substr($ATOM[$N1NAM],38,8),substr($ATOM[$N1NAM],46,8));
				@CC2 = (substr($ATOM[$CC2NAM],30,8),substr($ATOM[$CC2NAM],38,8),substr($ATOM[$CC2NAM],46,8));
				&DIHED(@O4, @CC1, @N1, @CC2, 5);
				$CHI[$NRES] = $TOR;
				}
			else {
				print "you have a freak ribonucleotide, position $NRES, type $RESTYP[$NRES]\n";
				}
			$P2NAM = $PIDX[$NRES+1];
			$C2NAM = $C4IDX[$NRES+1];
			$O52NAM = $O5IDX[$NRES+1];
			$P0NAM = $PIDX[$NRES - 1];
			@P2 = (substr($ATOM[$P2NAM],30,8),substr($ATOM[$P2NAM],38,8),substr($ATOM[$P2NAM],46,8));
			@C2 = (substr($ATOM[$C2NAM],30,8),substr($ATOM[$C2NAM],38,8),substr($ATOM[$C2NAM],46,8));
			@O52 = (substr($ATOM[$O52NAM],30,8),substr($ATOM[$O52NAM],38,8),substr($ATOM[$O52NAM],46,8));
			@P0 = (substr($ATOM[$P0NAM],30,8),substr($ATOM[$P0NAM],38,8),substr($ATOM[$P0NAM],46,8));
			&MEAS(@C1,@P2);
			$DISCP[$NRES] = $DIST;
			&ANG(@P1,@C1,@P2);
			$CANG[$NRES] = $THETA;
			&DIHED(@C3,@O3,@P2,@O52, 6);
			$ZETA[$NRES] = $TOR;
			&DIHED(@P1,@C1,@P2,@C2, 7);
			$PCPC[$NRES] = $TOR;
			&DIHED(@C1,@C3,@O3,@P2, 8);
			$EPSILON[$NRES] = $TOR;
			$C0NAM = $C4IDX[$NRES - 1];
			$O30NAM = $O3IDX[$NRES - 1];
			@C0 = (substr($ATOM[$C0NAM],30,8),substr($ATOM[$C0NAM],38,8),substr($ATOM[$C0NAM],46,8));
			@O30 = (substr($ATOM[$O30NAM],30,8),substr($ATOM[$O30NAM],38,8),substr($ATOM[$O30NAM],46,8));
			&ANG(@C0,@P1,@C1);
			$PANG[$NRES] = $THETA;
			&DIHED(@C0,@P1,@C1,@P2, 9);
			$CPCP[$NRES] = $TOR;
			&DIHED(@O30,@P1,@O5,@C5, 10);
			$ALPHA[$NRES] = $TOR;
			}
		}
	}
#
#  At this point the measurements for the file have been processed
#  The following sections check to see if nucleotides within this file
#  fall within user specified ranges
#
#  Call routine TYPEA 
# 
for ($i = 1; $i <= $NUMRES; $i++) {
	&TYPEA($ALPHA[$i], $BETA[$i], $GAMMA[$i], $DELTA[$i], $EPSILON[$i], $ZETA[$i]);
	if ($TYPA == 1) {
		$AHLX[$i] = 1;
		$GOOD++;
		}
	else {
		$AHLX[$i] = 0;
		$BAD++;
		}
#
#  Call routine AREA 
#
	&AREA($CPCP[$i], $PCPC[$i]);
        if ($i == $ENDRES[$#ENDRES]) {
                $end = $ATNUM;
                }
        else {
                $end = $RESSTRT[$i+1] - 1;
                }
	}
	}
#
#  A subroutine to measure distances between atoms
#
sub MEAS {
	local (@HERE) = @_;
	$X2 = ($HERE[0] - $HERE[3])**2;
	$Y2 = ($HERE[1] - $HERE[4])**2;
	$Z2 = ($HERE[2] - $HERE[5])**2;
	$DIST = ($X2 + $Y2 + $Z2)**(1/2);
	}
#
#  A subroutine to measure the angle between three atoms
#
sub ANG {
	local (@NOW) = @_;
	local ($SIN, $COS);
	foreach $foo (@NOW)   {
                if (/!\d/) {
                        $THETA = 400;
                        print "$foo in sub ANG, res $NRES of $FILE is not a number\n";
                        return;
                        }
                }
	$AX = ($NOW[0] - $NOW[3]);
	$BX = ($NOW[6] - $NOW[3]);
	$AY = ($NOW[1] - $NOW[4]);
	$BY = ($NOW[7] - $NOW[4]);
	$AZ = ($NOW[2] - $NOW[5]);
	$BZ = ($NOW[8] - $NOW[5]);
	$COS = (($AX * $BX) + ($AY * $BY) + ($AZ * $BZ))/ ((($AX**2 + $AY**2 + $AZ**2)*($BX**2 + $BY**2 +$BZ**2))**(1/2));
	$SIN = (1 - ($COS**2))**(1/2);
	$THETA = atan2($SIN, $COS)*(180/$pi);
	}
#
#  A subroutine to measure the torsion defined by four atoms
#
sub DIHED {
	local (@DO) = @_;
	local ($SIN, $COS, $NUCOS);
	foreach $foo (@DO) {
		if (/!\d/ || $foo == 0) {
			$TOR = 400;
			}
		}
	$PX = ($DO[0] - $DO[3]);
	$PY = ($DO[1] - $DO[4]);
	$PZ = ($DO[2] - $DO[5]);
	$QX = ($DO[6] - $DO[3]);
	$QY = ($DO[7] - $DO[4]);
	$QZ = ($DO[8] - $DO[5]);
	$RX = ($DO[9] - $DO[3]);
	$RY = ($DO[10] - $DO[4]);
	$RZ = ($DO[11] - $DO[5]);
        $LNX = ($QY * $RZ) - ($QZ * $RY);
        $LNY = ($QZ * $RX) - ($QX * $RZ);
        $LNZ = ($QX * $RY) - ($QY * $RX);
        $MNX = ($QY * $PZ) - ($QZ * $PY);
        $MNY = ($QZ * $PX) - ($QX * $PZ);
        $MNZ = ($QX * $PY) - ($QY * $PX);
        $COS = (($LNX*$MNX)+($LNY*$MNY)+($LNZ*$MNZ))/((($LNX**2+$LNY**2+$LNZ**2)*($MNX**2+$MNY**2+$MNZ**2))**(1/2));
	$SIN = (1 - ($COS**2))**(1/2);
	$NUCOS = (($MNX*$RX)+($MNY*$RY)+($MNZ*$RZ))/((($MNX**2+$MNY**2+$MNZ**2)*($RX**2+$RY**2+$RZ**2))**(1/2));
	$TOR = (atan2($SIN, $COS)*(180/$pi))*($NUCOS/(($NUCOS**2)**(1/2)));
	if ($TOR < 0){
		$TOR = 360 + $TOR
		}
	}
#
#  Subroutine to determine if nucleotides fall within user-specified ranges
#  of eta and theta. Nucleotides which do are added to the files 
#  all_area.txt AND xxx.pdb_area.txt
#
sub AREA {
	local (@MSMT) = @_;
	@UPLIM = ( 
		190, 	# change to high limit of eta range
		260 	# change to high limit of theta range
		);
	@LOWLIM = ( 
		150,	# change to lower limit of eta range
		190	# change to lower limit of theta range
		);
	for ($z=0; $z <= $#MSMT; $z++) {
		if (($MSMT[$z] == 0) && ($MSMT[$z+1] == 0)) {
			return;
			}
#
#  The code now finds nucleotides outside of the specified region
#  To change it to find nucleotides inside the specified region comment 
#  the next executable line and remove the comment (#) from the next line
#             if (($MSMT[$z] < $LOWLIM[$z}) || ($MSMT[$z] > $UPLIM[$z])) {
#
              if (($MSMT[$z] > $LOWLIM[$z]) && ($MSMT[$z] < $UPLIM[$z])) {
			return;
			}
		}
	printf AREAOUT "%-7s\t%6.1f\t%6.1f\n", $RESNAM[$i], $MSMT[0], $MSMT[1];
	printf ALLAREA"%-5d\t%-7s\t%-4s\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f%\t%6.1f\t%3d\n", $i, $RESNAM[$i], $RESTYP[$i], $CPCP[$i], $PCPC[$i], $CHI[$i], $ALPHA[$i], $BETA[$i], $GAMMA[$i], $DELTA[$i], $EPSILON[$i], $ZETA[$i], $NUTWO[$i], $AHLX[$i];
	$AREACOUNT++;
	}
#
#  Subroutine to determine if nucleotides fall within user-specified ranges
#  of the standard torsions. Nucleotides which do are noted with a "1"
#  in the TYPEA column of output files all_sprd.txt and xxx.pdb_sprd.txt
#
sub TYPEA {
	$TYPA = 0;
	local (@MSMT) = @_;
	@UPLIM = ( 

                360,	# change to high limit for alpha 
                210, 	# change to high limit for beta
                360, 	# change to high limit for gamma
                100, 	# change to high limit for delta
                225, 	# change to high limit for epsilon
                340	# change to high limit for zeta
                );
        @LOWLIM = ( 
                155,	# change to low limit for alpha
                149, 	# change to low limit for beta
                0,	# change to low limit for gamma
                65,	# change to low limit for delta
                155,	# change to low limit for epsilon
                250	# change to low limit for zeta
                );
	$TYPA = 1;
	}
