NOTE: this documentation file and the associated software are
      never sent unsollicited; if you get this material without
      any request from your side, please disregard it, and 
      inform the authur at the address below if the problem 
      persists, including the header(s) of the e-mail message(s).


The `SN' library for the SKEW NORMAL distribution (Matlab version)

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


Requirements:
- The procedures have been developed and tested on Matlab 5.
- The statistical toolbox stats and the optimisation toolbox optim 
  are required. 

Home page:
the home page of the SN distribution, where you can find related material, is 
http://www.stat.unipd.it/dip/homes/azzalini/SN

Licence: GNU General Public Linense, see file COPYING.

