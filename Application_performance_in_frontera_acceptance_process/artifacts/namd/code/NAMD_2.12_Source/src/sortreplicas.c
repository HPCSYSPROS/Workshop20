
#include "largefiles.h"  /* must be first! */

#include "vmdplugin.h"

extern int molfile_dcdplugin_init(void);
extern int molfile_dcdplugin_register(void *, vmdplugin_register_cb);
extern int molfile_dcdplugin_fini(void);
/*
extern int molfile_jsplugin_init(void);
extern int molfile_jsplugin_register(void *, vmdplugin_register_cb);
extern int molfile_jsplugin_fini(void);
*/

#include "molfile_plugin.h"

static int register_cb(void *v, vmdplugin_t *p) {
  *((vmdplugin_t **)v) = p;
  return VMDPLUGIN_SUCCESS;
}

#include <sys/types.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

int main(int argc, char **argv) {

  molfile_timestep_t frame;
  molfile_plugin_t *plugin;
  char *output_root;
  char *filename;
  int num_replicas;
  int runs_per_frame;
  long long int final_step = -1;
  long long int checkstep = -1;
  int colvars;
  FILE **hist_in;
  FILE **hist_out;
  FILE **colv_in;
  FILE **colv_out;
  void **traj_in;
  void **traj_out;
  int natoms=MOLFILE_NUMATOMS_UNKNOWN;
  int i, i_run;

  molfile_dcdplugin_init();
  molfile_dcdplugin_register(&plugin, register_cb);

  if ( argc < 4 || argc > 5 ) {
    fprintf(stderr, "args: <job_output_root> <num_replicas> <runs_per_frame> [final_step]\n");
    exit(-1);
  }
  output_root = argv[1];
  num_replicas = atoi(argv[2]);
  runs_per_frame = atoi(argv[3]);
  if ( argc > 4 ) {
    sscanf(argv[4], "%lld", &final_step);
  }

  filename = (char*) malloc(strlen(output_root)+100);
  hist_in = (FILE**) malloc(num_replicas*sizeof(FILE*));
  hist_out = (FILE**) malloc(num_replicas*sizeof(FILE*));
  colv_in = (FILE**) malloc(num_replicas*sizeof(FILE*));
  colv_out = (FILE**) malloc(num_replicas*sizeof(FILE*));
  traj_in = (void**) malloc(num_replicas*sizeof(FILE*));
  traj_out = (void**) malloc(num_replicas*sizeof(FILE*));

  for ( i=0; i<num_replicas; ++i ) {
    char *root_end;
    if ( strstr(output_root,"%s") ) {
      char istr[10];
      sprintf(istr,"%d",i);
      sprintf(filename,output_root,istr);
    } else {
      sprintf(filename,output_root,i);
    }
    root_end = filename + strlen(filename);

    sprintf(root_end,".%d.history",i);
    hist_in[i] = fopen(filename,"r");
    if ( ! hist_in[i] ) {
      fprintf(stderr, "error opening input file %s: %s\n",
					filename, strerror(errno));
      exit(-1);
    }
    sprintf(root_end,".%d.dcd",i);
    traj_in[i] = plugin->open_file_read(filename,"dcd",&natoms);
    if ( ! traj_in[i] ) {
      fprintf(stderr, "error opening input file %s: %s\n",
					filename, strerror(errno));
      exit(-1);
    }
    sprintf(root_end,".%d.colvars.traj",i);
    colv_in[i] = fopen(filename,"r");
    if ( colv_in[i] ) {
      if ( i == 0 ) {
        printf("Found first input colvars trajectory file %s.\n", filename);
        colvars = 1;
      } else if ( ! colvars ) {
        fprintf(stderr, "missing input colvars trajectory files before %s\n", filename);
        exit(-1);
      }
    } else {
      if ( i == 0 ) {
        colvars = 0;
      } else if ( colvars ) {
        fprintf(stderr, "error opening input colvars trajectory file %s: %s\n",
					filename, strerror(errno));
        exit(-1);
      }
    }
  }

  for ( i=0; i<num_replicas; ++i ) {
    char *root_end;
    if ( strstr(output_root,"%s") ) {
      char istr[10];
      sprintf(istr,"%d",i);
      sprintf(filename,output_root,istr);
    } else {
      sprintf(filename,output_root,i);
    }
    root_end = filename + strlen(filename);

    sprintf(root_end,".%d.sort.history",i);
    hist_out[i] = fopen(filename,"w");
    if ( ! hist_out[i] ) {
      fprintf(stderr, "error opening output file %s: %s\n",
					filename, strerror(errno));
      exit(-1);
    }
    sprintf(root_end,".%d.sort.dcd",i);
    traj_out[i] = plugin->open_file_write(filename,"dcd",natoms);
    if ( ! traj_out[i] ) {
      fprintf(stderr, "error opening output file %s: %s\n",
					filename, strerror(errno));
      exit(-1);
    }
    if ( colvars ) {
      sprintf(root_end,".%d.sort.colvars.traj",i);
      colv_out[i] = fopen(filename,"w");
      if ( ! colv_out[i] ) {
        fprintf(stderr, "error opening output file %s: %s\n",
					filename, strerror(errno));
        exit(-1);
      }
    }
  }

  frame.coords = (float*) malloc(3*natoms*sizeof(float));
  frame.velocities = (float*) NULL;

#define LINE_MAX 10000

  i_run = 0;
  for ( ; 1; ++i_run ) { /* loop until read fails */
    char line[LINE_MAX];
    for ( i=0; i<num_replicas; ++i ) {
      char *r;
      char sav;
      int f1,f2;
      int rc;
      int rep_id = -1;
      long long int step;
      r = fgets(line, LINE_MAX, hist_in[i]);
      if ( ! r ) { break; }
      rc = sscanf(line, "%lld %n%d%n", &step, &f1, &rep_id, &f2);
      if ( rc != 2 ) {
        fprintf(stderr,"Format error for replica %d at line %d: %s",
							i, i_run, line);
        exit(-1);
      }
      if ( i == 0 ) {
        if ( step <= checkstep ) {
          fprintf(stderr,"Step out of order for replica %d at line %d: %s",
							i, i_run, line);
          exit(-1);
        }
        checkstep = step;
        if ( final_step >= 0 && checkstep > final_step ) {
          printf("Stopping after final step %lld.\n", final_step);
          break;
        }
      } else if ( step != checkstep ) {
        fprintf(stderr,"Step mismatch for replica %d at line %d: %s",
							i, i_run, line);
        exit(-1);
      }
      if ( rep_id < 0 || rep_id >= num_replicas ) {
        fprintf(stderr,"Invalid replica ID for replica %d at line %d: %s",
							i, i_run, line);
        exit(-1);
      }
      sav = line[f1];
      line[f1] = 0;
      fprintf(hist_out[rep_id],"%s%d%s",line,i,line+f2);
      line[f1] = sav;
      if ( colvars ) {
       long long int oldcstep = -1;
       while ( 1 ) {
        long long int cstep;
        char cline[LINE_MAX];
#ifdef WIN32
        __int64 oldpos = _ftelli64(colv_in[i]);
#else
        off_t oldpos = ftello(colv_in[i]);
#endif
        r = fgets(cline, LINE_MAX, colv_in[i]);
        if ( ! r ) { break; }
        if ( cline[0] == '#' ) {
          fprintf(colv_out[rep_id],"%s",cline);
          continue;
        } 
        rc = sscanf(cline, "%lld", &cstep);
        if ( rc != 1 ) {
          fprintf(stderr,"Format error in colvar trajectory for replica %d: %s",
							i, cline);
          exit(-1);
        }
        if ( cstep == oldcstep ) continue;  /* filter out repeats */
        if ( cstep < oldcstep ) {
          fprintf(stderr,"Step out of order in colvar trajectory for replica %d: %s",
							i, cline);
          exit(-1);
        }
        if ( cstep > step ) {
#ifdef WIN32
          _fseeki64(colv_in[i], oldpos, SEEK_SET);
#else
          fseeko(colv_in[i], oldpos, SEEK_SET);
#endif
          break;
        }
        if ( i_run != 0 || oldcstep != -1 ) {  /* skip first entry */
          fprintf(colv_out[rep_id],"%s",cline);
        }
        oldcstep = cstep;
       }
      }
      if ( (i_run+1) % runs_per_frame ) continue;
      rc = plugin->read_next_timestep(traj_in[i],natoms,&frame);
      if ( rc == MOLFILE_SUCCESS ) {
        plugin->write_timestep(traj_out[rep_id],&frame);
      } else {
        fprintf(stderr,"Unable to read frame for replica %d at line %d: %s",
							i, i_run, line);
        break;
      }
    }
    if ( i < num_replicas ) {
      printf("Processed %d runs.\n",i_run);
      if ( i ) fprintf(stderr,"Uneven input lengths for replica %d at line %d: %s",
							i, i_run, line);
      break;
    }
  }

  free(frame.coords);

  for ( i=0; i<num_replicas; ++i ) {
    if ( fclose(hist_in[i]) ) {
      fprintf(stderr, "error closing history input file %d: %s\n", i, strerror(errno));
    }
    plugin->close_file_read(traj_in[i]);
    if ( fclose(hist_out[i]) ) {
      fprintf(stderr, "error closing history output file %d: %s\n", i, strerror(errno));
    }
    plugin->close_file_write(traj_out[i]);
    if ( colvars ) {
      if ( fclose(colv_in[i]) ) {
        fprintf(stderr, "error closing colvars input file %d: %s\n", i, strerror(errno));
      }
      if ( fclose(colv_out[i]) ) {
        fprintf(stderr, "error closing colvars output file %d: %s\n", i, strerror(errno));
      }
    }
  }

  molfile_dcdplugin_fini();
  exit(0);
}

