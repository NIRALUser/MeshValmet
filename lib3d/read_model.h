
#ifndef READ_MODEL_FILE
#define READ_MODEL_FILE

#include <model_in.h>

#include <string.h>

/* Reads a model from file 'fname' and returns the model read. If an error
 * occurs a message is printed and the program exists. 
 */
static struct model *read_model_file(const char *fname)
{

  int rcode;
  struct model *m;
  const char *errstr;
  
  int character = '.';
  const char ext1[] = ".raw";
  const char ext2[] = ".byu";
  const char ext3[] = ".iv";
  
  const char* pointer = strrchr(fname,character);
  
  if(strcmp(pointer,ext1) == 0 || strcmp(pointer,ext3) == 0)
    rcode = read_fmodel(&m,fname,MESH_FF_AUTO,1);
  else if(strcmp(pointer,ext2) == 0)
    rcode = read_fmodel(&m,fname,MESH_FF_BYU,1);
  
  if (rcode <= 0) {
    switch (rcode) {
    case 0:
      errstr = "no models in file";
      break;
    case MESH_NO_MEM:
      errstr = "no memory";
      break;
    case MESH_CORRUPTED:
      errstr = "corrupted file or I/O error";
      break;
    case MESH_MODEL_ERR:
      errstr = "model error";
      break;
    case MESH_NOT_TRIAG:
      errstr = "not a triangular mesh model";
      break;
    case MESH_BAD_FF:
      errstr = "unrecognized file format";
      break;
    case MESH_BAD_FNAME:
      errstr = strerror(errno);
      break;
    default:
      errstr = "unknown error";
    }
    fprintf(stderr,"ERROR: %s: %s\n",fname,errstr);
    exit(1);
  } else if (m->num_faces == 0) {
    fprintf(stderr,"ERROR: %s: empty model (no faces)\n",fname);
    exit(1);
  }
  return m;
}
#endif
