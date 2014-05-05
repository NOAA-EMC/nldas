#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "preproc.h"

/* slimmed down version of the VIC code */
FILE *open_file(char string[], char type[])
{
  /* open a file */

  FILE *stream;

  stream = fopen(string,type);

  if(stream == NULL){
    fprintf(stderr, "\nUnable to open file \"%s\"\n\n", string);
    exit(EXIT_FAILURE);
  }

#if VERBOSE
  fprintf(stderr,"\n \"%s\" has been",string);
  if(strcmp(type,"r") == 0)
    fprintf(stderr," opened for reading.\n");
  else if(strcmp(type,"w") == 0)
    fprintf(stderr," truncated or created for writing.\n");
#endif

  return stream;

}
