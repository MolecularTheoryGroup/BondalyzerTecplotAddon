This directory contains a set of sample code used to create a general purpose add-on.  It include a sidebar (entitled "My Sidebar") and a modeless dialog (entitled "General Purpose dialog").

Refer to the ADK User's Manual for details on how to create a general purpose add-on.


--------------------------------------------------------------------------------------


While debugging your DLL, it may be easier for you
to run Tecplot directly from you Visual C++ DLL project.
To do this:

1. Select Project/Properties/Configuration Properties/Debugging.
2. Set "Command" to the Tecplot executable (include the full path if necessary).
3. Set "Command Arguments" to "-loadaddon GENERALPURPOSESAMPLE"
4. Set "Working Directory" to "Debug".
