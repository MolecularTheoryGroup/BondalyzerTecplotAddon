/*
****************************************************************
****************** BEGIN DEVELOPMENT NOTES *********************
****************************************************************

D ****************************************************************
D *                 Build 1.0 9-04-98                            *
D ****************************************************************





****************************************************************
****************** END DEVELOPMENT NOTES ***********************
****************************************************************
****************************************************************
*  D in column 1 marks date information.                       *
*  C in column 1 marks notes on new changes.                   *
*  B in column 1 marks notes on bug fixes.                     *
****************************************************************

*/

#define ADDON_NAME "Simple Spreadsheet Converter"
#include "TecplotVersion.h"
#define ADDON_VERSION TecVersionId
#ifdef __DATE__
# define ADDON_DATE __DATE__
#else
# define ADDON_DATE ""
#endif /* __DATE__ */
#define MinTecplotVersionAllowed 740503

