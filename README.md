

DEPENDENCIES
------------

depends on lapack, and  assumes the lapack is installed in /usr/lib/ 
(sudo apt-get install liblapack-dev)


INSTALLING
----------

after unpacking in $STRUCT_HOME directory
```bash 
> cd $STRUCT_HOME/10_objects/
> make
> make clean
```

the executable called struct should now be in $STRUCT_HOME



USAGE
-----
```bash
$STRUCT_HOME/struct -in/-from <pdb/db tgt file> [-c1 <tgt chain>] \
		    [ -to <pdb/db qry file>] [ -c2 <qry chain>] \
		    [-max_out <number of almts to output>] [ -v] \
		    [-no_bb] [ -p <parameter file>]


  -from/to  input files in db or pdb format (to do the alignment 
          on backbone level, both files should be pdb
  -c1/2   chains for the pdb inputs 1 and 2
  -max_out number of alignments (pdb and allignemnt files) to output
  -no_bb  do not do the alignment on the backbone level
  -v      verbose output
  -p      parameter file
``````


TEST CASES
----------

## 1-on-1 comparison

```bash
> cd $STRUCT_HOME/11_tests/01_2d8bA_1d0nA/
> $STRUCT_HOME/struct -from 2d8bA.pdb  -to 1d0nA.pdb
```

view the results:
```bash 
> cat   2d8bA_1d0nA.struct_out
> pymol 1d0nA.pdb 2d8bA.to_1d0nA.*.pdb
```
If you do not use pymol, you can view these files in your
favorite molecular viewer by finding them under open->file menu.


## Making a db file 
A db file is a file containing directions and center points for
   each tentative, heuristically determined, element of secondary
   structure; it is used for fast searching though a large set 
   of structures, such as PDB itself.
```bash
> cd $STRUCT_HOME/11_tests/01_2d8bA_1d0nA/
> $STRUCT_HOME/struct -in  1d0nA.pdb
```

Note: it is possible to use the full pdb file here,
provided you have it available in your directory:
```bash
> $STRUCT_HOME/struct -in 1d0n.pdb -c1 A 
```

view the output:
```bash
> cat 1d0nA.db
```

## Database search

```bash
> cd $STRUCT_HOME/11_tests/02_db_search
> $STRUCT_HOME/struct -from small_test.db -to ../01_2d8bA_1d0nA/2d8bA.db
```
View the output:
```bash
> cat digest.struct_out
```

View the output sorted by the direction score
```bash
> awk '$1 != "%" && $1 != "done"' digest.struct_out | sort -grk 5
```
Note: the "database" we are searching here is a  
   concatenation of the db files, like the one produced in  
   the testcase (b); 
   direction is the only info we have in this case.



