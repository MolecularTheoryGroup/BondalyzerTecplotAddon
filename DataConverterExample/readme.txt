------------------------------------------------
*** Note for WINDOWS users *** : This example has
changed. To build it, create a new project of
type "win32 dll" in VC++ and add the following files:

cnvss.c
cnvss.h
------------------------------------------------


cnvss -- Addon Example Dataset Converter

This example reads a simple space-delimited spreadsheet
file. The file must be in the following form:

  Variable 1  Variable 2   ...   Variable N

  1.1 1.2 1.3
  2.1 2.2 2.3
  .
  .
  .
  n.1 n.2 n.3

  [Blank lines are ignored]

  Constraints:

  1) all of the variables names must be on the same line
  2) subsequent lines must have *exactly* that number of variables

Illustrates:
-- Reading ascii field data into Tecplot
