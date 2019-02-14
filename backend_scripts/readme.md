## Program connections file
The csv files stores the connections between all the programs that are present. 
The way to read it is:
  - rows indicate output, as in the output of STAR is indicated in the STAR row
  - columns indicate inputs, as in the output of this row as input of the program in this column
 
 There are different possible values:
  - S indicates the same program for input and output, this will throw an error
  - -1 indicates no connection between the two programs, this happens with programs of 
     the same type or if an advanced step is fed into a previous one
  - 0 indicates direct connection between the two programs
  - an INT indicates connection between the two programs through the program represented by the INT
  - an INT + C indicates connection through multiple steps
