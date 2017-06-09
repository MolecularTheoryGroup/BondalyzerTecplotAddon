#include <cstdlib>
#include <cstdio>
#include "parsearguments.h"
using namespace std;

void ParseArguments(int argc,
		    char**argv,
		    char*& outputdir)
{
  // parse arguments
  if(argc>1)
    {
      for(int arg_idx=1; arg_idx<argc; arg_idx++)
        {
	  if(argv[arg_idx][0]=='-')
            {
	      if(argv[arg_idx][1]=='o' && ++arg_idx < argc)
		{
		  outputdir=argv[arg_idx];		  
		}
	      else
		{
		  printf("usage: %s [-o outputdir] \n -o outputdir : use outputdir as output directory [default: mesh_output]\n\n",
			 argv[0]);
		  exit(1);
		}
            }
	  else
            {
	      printf("usage: %s [-o outputdir] \n -o outputdir : use outputdir as output directory [default: mesh_output]\n\n",
		     argv[0]);
	      exit(1);
            }
        }
    }
  return;
}
