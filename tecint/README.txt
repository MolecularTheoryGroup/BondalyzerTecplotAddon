These directories contain the Unix and Windows distribution of Tecplot 
Integration Library. /doc contains subdirectories PDF and HTML, each of which 
contain the documentation for the library. /include contains the C include file 
to be included in any source files that use the library API's. /lib contains 
subdirectories that each contain the library, libtecint.a (tecint.lib for Windows)
for the platform indicated by the directory name (these are the same platform 
names used by Tecplot). Directory win32 contains the release version of the Windows
library, created with multi-threaded DLL settings.

To use the library with your Tecplot add-on, include file tecint.h in your source
and modify the link command to include the library for your particular platform.
Refer to the documentation for use of the Integrate and related routines.
