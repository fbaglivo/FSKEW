SKEW NORMAL (Matlab version)

To install Skew normal procedures do the following steps:
1) choose, or create, a folder (example: c:\skewn\library) 
2) copy snmatlab.zip in this folder and extract all its files
3) open Matlab and add the folder containing the procedures in
	the existing path, with the command:
	addpath "skew normal path" (example: addpath c:\skewn\library)
   Note: if you want to add this path permanently, choose 'Set path...'
	 from 'File' menu; then choose 'Add to path...', insert the
	 "skew normal path" (example: c:\skewn\library) and use
	 the 'Save settings' option.

Now, you are able to use all procedures of skew normal library.
See help for details.

How to use help

In this version of Matlab skew normal library, three types of help are 
available:

1) online help: this can be seen in the command window typing  

>> help "procedure-name"

example: the command

>> help sn_cumulants

shows help for sn_cumulants procedure.

2) Help window: from menu 'Help' of Matlab choose 'Help Window' and then
'Skew normal distribution library'.

3) Help HTML: there is also a HTML version of the help that can be seen
with a Web Browser. You can load this help with the command:

>> web "skew normal path"\sn.html

example:

>> web c:\skewn\library\sn.html
