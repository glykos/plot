# plot
Lightweight plotting program for GNU/Linux and MacOS that can read standard input and can be used in unix pipelines
____________________________________________________________________

0. WINDOZE

   If you work with a windoze machine, forget it. It won't happen.
   Delete this file, its folder, everything. End of story. Really.

1. GNU/LINUX

   If you work with GNU/Linux you are in good shape. If it doesn't
   work for you, I'll fix it. But before getting me to fix it, 
   confirm that it is indeed broken :

      * Open a shell   
      * Go to `plot/bin/linux/`
      * Type :   `./plot -cc < ../../example/matrix.dat`
      * Move the cursor in the graphics window, press 'c'.
	After a while press 'n'.
      * Press 'q' to quit.

   If you get a segmentation fault any time you try to run plot,
   your stack size is limited. Type in the terminal 
   `ulimit -s unlimited` if you use bash, or,
   `limit stacksize unlimited` if you work with tcsh.

   If you get any other error message, or you don't see the graphs,
   it is indeed broken and you might as well send a bug report 
   with a problem description and details of your machine to
   glykos@mbg.duth.gr . If it is not broken, copy the program
   somewhere in the users' path (eg. `/usr/local/bin/`).

   Famous last words : the GNU/Linux executable that I provide has
   been tested with various distributions ranging from RedHat 7.3
   (released May 2002) to Ubuntu 18.04 (released April 2018). If
   you can find a 32bit GNU/Linux machine that this executable 
   fails to run, I'd be grateful if you could let me know.
   For 64bit machines you either need the 32bit libraries, or
   you can try the 64bit executable that I provide.

 

2. MacOSX 
   
   See (1) above. The only difference is in the names of the 
   directories containing the executables. Please note that I do
   not have access to a MacOSX machine and I despise proprietary
   operating systems (the MaxOSX executables were cross-compiled
   on a GNU/Linux machine). I'd love to help if I can, but I do
   have better things to do than dealing with Apple's binary
   compatibility issues [see (3) below for other solutions].

3. ADVENTUROUS UNIX USERS

   If you have gcc and g77 (or gfortran) and you can build Ygl, 
   you might as well try to compile 'plot' from the source code
   located in src/. The file `2comp_gnu` contains the somewhat
   minimalistic and out-of-date instructions for building an
   executable. Of course, it will not work. But I'll be happy to
   assist if I can ...

4. ALL USERS WITH A WORKING EXECUTABLE

   This program has documentation !
____________________________________________________________________

