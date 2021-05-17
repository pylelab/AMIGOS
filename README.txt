A newer version of this program, AMIGOS II, is now available at
https://pylelab.org/software

-------------------------------------------------------------------------------

Algorithmic Method for Identifying Groupings of Overall Structure
(AMIGOS)
By Carlos M Duarte, Copyright 1998, 1999, 2003
Modified 05/17/2021 by Chengxin Zhang
Modified 08/30/2003 by Leven Wadley & Carlos Duarte
	 (Modification listed @ end of file)
Modified 03/08/1999 (Modifications listed @ end of file)
Modified 09/08/1998

AMIGOS is a perl script which outputs tables of torsion angles from 
nucleic acid PDB files. Since it is a perl script it requires no compilation
as long as perl is properly installed on the machine you are working on. 
To run it, make the file executable, and execute it.

By default AMIGOS measures the standard backbone torsion angles (alpha,
beta, gamma, delta, epsilon & zeta), the sugar pucker torsion (nu 2), chi,
and the pseudotorsions eta & theta. It will also output tables of nucleotides
whose measurements fall within user-specified ranges.

REQUIREMENTS
-  perl4 (or later) compiler

-  A directory with some nucleic acid PDB files

(YOU CAN RUN THE PROGRAM WITH THE INFORMATION ABOVE. IF YOU WANT TO KNOW
MORE, READ ON) 

EXECUTING THE PROGRAM

-  place AMIGOS.pl somewhere in your path.

-  make sure you are in a directory in which you have write access.

-  type "AMIGOS.pl".

-  if that fails, try "perl AMIGOS.pl" in the directory where AMIGOS.pl
   resides.

-  enter the name of the directory containing the nucleic acid PDB file(s).
	When specifying the directory containing the pdb files, use the full 
	path of the directory (Home directory shortcuts, e.g. ~bozo, 
	don't work). 

-  output will scroll across the screen telling you which file is being 
	processed and whether the files have residues with unusual names.

-  The directory you are working in will start to fill with a variety of files 
	which are explained below.

KNOWN BUGS OR ISSUES

PDB files that do not conform to the PDB Format Description (available
at http://www.pdb.org) may not be read correctly. This situation
arises most frequently when third-party software is used to output
files in PDB format. Complain to them if you have this problem.

HETATM entries in PDB files are ignored. This can mean that modified
bases in some PDB files (e.g. some tRNAs) will not be
considered. Likewise, bases that are adjacent to bases that contain
HETATM entries will not have torsions calculated. This choice in
program design was intentional.

OUTPUT:
	In its current form the program creates four different types of files.
1) An output file of the measurements of all nucleotides form all files
	(all_sprd.txt)

2) An output file with all nucleotides which fall within user specified ranges
   of eta & theta (all_area.txt)

For each PDB file there are two files created
3) An output file of the measurements of all nucleotides (filename_sprd.txt)

4) An output file with all nucleotides which fall within user specified ranges
   of eta & theta (filename_area.txt)

	The total number of files created is 2n+2, where n is the number of
PDB files in the directory you are processing.

	The output does not contain measurements for nucleotides at the start
or end of a nucleotide chain. The start and end of chains are determined 
by the PDB file CHAIN ID (column 22 of ATOM entry) designation for each 
residue. If there are intervening nucleotides which are in the PDB file as
HET (heteroatoms) they will not be considered in the measurement and this
can introduce error (This case is seen in some tRNA files) 

	While the program is running the on screen report (standard output) 
tells you what file is being processed and how many chains the script thinks 
it has. The program also reports residues with unusual names, and chi is 
not calculated for these.

INPUT (PDB FIles):
	The program assumes that any file ending in "ent" or "pdb" is a 
properly formatted PDB nucleic acid file. Any residue which does not
contain an "O2'" or "O2*" is ignored for geometric calculations.  This
effectively screens out all amino acids and DNA.  It assumes that properly 
named residues are "A", "C", "G", "U", or "T", and reports the name of 
all other residues which contain an O2'.  Backbone atom names can either 
end in "'" or "*".

PROGRAM:
	The program has routines to measure distances, and angles and can
also be expanded to measure any four atoms the user desires, if one cares
to pick through the code. The variable names are not straightforward when
it comes to the arrays for each atom used to make the measurements. This
was not purposefully done, but is historical, as the program grew to do
much more than originally intended.

ADVANCED MODIFICATIONS:
	By default the area files contain the measurements of all nucleotides
which fall outside the Helical region (Duarte & Pyle, 1998). This can be
modified in the script, but you've got to edit the file. I've included 
instructions on how to do this in the comments of the script. It is not 
difficult to modify this, and if the comments don't explain it well enough
e-mail or phone me and I'll be happy to walk you through it.

        The program has routines to measure distances, and angles and can
also be expanded to measure any four atoms the user desires, if one cares
to pick through the code. The variable names are not straightforward when
it comes to the arrays for each atom used to make the measurements. This
was not purposefully done, but is historical, as the program grew to do
much more than originally intended.

	Good luck, and thanks for using the code.

	Please e-mail or me with any problems or questions

Carlos Duarte	(cmd63@columbia.edu) 

	This program was used in the following paper 

	Duarte, C. & Pyle, A.M., (1998) Stepping through an RNA structure:
a novel approach to conformational analysis, Journal of Molecular Biology,
284(5):1465-78. 

*** Modifications of 5/17/21 ***

- Program modified to handle division of 0.

*** Modifications of 8/30/03 ***

- Program modified to permit multiple conformations of the same
  atom. Currently, only the first conformation is used and remaining
  conformations are ignored.


*** Modifications of 3/8/99 ***

- Program altered to only calculate measurements for residues with O2'
or O2* atoms.  This effectively screens out all amino acids and DNA
nucleotides.

- Program changed to stop reading PDB file when "END" is encountered. 

