
## Unix

for Mac users you need to edit the first part of the Makefile.def, uncomment
the first line and comment or delete the second line!


there are a number of example classes and example scripts to test these
classes. In addition for elements and materials there are example c and 
fortran routines.

typing make in the DEVELOPER directory will build them all. the resulting
shared object libraries (dynamic link libraries) will be found in the same
directory as the source code.

to test them you need to ensure that the LD_LIBRARY_PATH is set to include ./

to test them run 'OpenSees example1.tcl', note OpenSees must be on your path.

## Windows

1. Open ViualStudio
2. File -> New Project
3. Give it name of New class and location of current class files (ElasticPPcpp and C:\???\OpenSeesDeveloper\material (VisualStudio2010 want VisualStudio C++ Win32 - Win32 Project)
4. Select OK (visual Studio20
5. New windows pops up, Select Application Settings
6. Select DLL AND select Empty Project
7. Select Finish. A new project pops up in the workspace with class name (ElasticPPcpp)

8.  Right click on source files, select add Existing Item.
9.  browse to the class file (should be up 1 directory, select it (`elasticPPcpp.cpp`)
10. Right click on heared files, select add existing item.
11. browse to header file, select it (elasticPPcpp.h)

12. right click on Project (ElasticPPcpp), select build.

It fails in compilation, as cannot find file `elementAPI.h`
13. right click on Project (ElasticPPCPP), select Properties.
14. select C/C++ folder icon
15. In Configuration pull down menu, select all configurations
16. Click 3 dots to right of Additional Include Directories
17. In window that pops up, select folder
18. add to line `..\..\core` (this directory contains the `elementAPI.h` file). (keep adding ..\ till it works!)
19. select ok.
20. right click on project (ElasticPPcpp), select build

It fails in Linking, a lot of unresolved external symbols.
21. Right click on source files, select add Existing Item.
22. browse to the core directory (should be 1 directory up)
23. select all the `.cpp` files here

24. right click on Project (ElasticPPCPP), select build.


Now test it, 

25. using windows finder, copy the `.dll` (C:\??\OpenSeesDeveloper\material\ElasticPPCPP\Debug) into the class directory (C:\???\material\)
26. run OpenSees (C:\??\OpenSeesDeveloper\bin)
27. `cd` into directory containing example (`cd ..\material`)
28. source an example script containing your new command, i.e. type: `source example1.tcl`

If it fails and you are using the x64 bit version you need to change the code generated
30. Under build, choose the configuration manager, you woill see the project marked as win32
31. in the pull down platform menu, choose new, x64 should pop up .. select that
32. build again


