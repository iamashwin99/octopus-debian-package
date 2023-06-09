/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <libgen.h>

#include <fortran_types.h>

#ifdef HAVE_DIRENT_H
#include <dirent.h>
#endif

#include "string_f.h" /* Fortran <-> c string compatibility issues */

/* *********************** interface functions ********************** */

void FC_FUNC_(oct_mkdir, OCT_MKDIR)
		 (STR_F_TYPE name STR_ARG1)
{
  struct stat buf;
  char *name_c;

  TO_C_STR1(name, name_c);
  if(!*name_c) return;

  if(stat(name_c, &buf) == 0){
    free(name_c);
    return;
  }

#ifndef _WIN32
  mkdir(name_c, 0775);
#else
  mkdir(name_c);
#endif

  free(name_c);
}

void FC_FUNC_(oct_stat, OCT_STAT)
     (fint *ierr, STR_F_TYPE name, STR_F_TYPE mod_time STR_ARG2)
{
  char *name_c, *mod_time_c;
  struct stat statbuf;
  time_t mtime;
  struct tm * timeinfo;

  TO_C_STR1(name, name_c);
  *ierr = stat(name_c, &statbuf);
  free(name_c);

  if(*ierr == 0)
  {
    mtime = statbuf.st_mtime; /* last modification time */
    timeinfo = localtime(&mtime);
    mod_time_c = asctime(timeinfo);
  }
  else
  {
    perror(name_c); /* what is the problem? */
    mod_time_c = malloc(sizeof(char));
    mod_time_c[0] = '\0';
  }

  TO_F_STR2(mod_time_c, mod_time);
  if(*ierr != 0)
  {
    printf("ierr = %i\n", *ierr);
    free(mod_time_c);
  }
  /* otherwise, do not do this since 'mod_time_c' points at static data of asctime */
}

int FC_FUNC_(oct_dir_exists, OCT_DIR_EXISTS)
     (STR_F_TYPE name STR_ARG1)
{
  int ierr;
  char *name_c;
  struct stat statbuf;

  TO_C_STR1(name, name_c);
  ierr = stat(name_c, &statbuf);
  free(name_c);

  if(ierr == 0)
  {
    return S_ISDIR(statbuf.st_mode);
  }
  else
  {
    return 0;
  }
}


void FC_FUNC_(oct_rm, OCT_RM)
     (STR_F_TYPE name STR_ARG1)
{
  char *name_c;

  TO_C_STR1(name, name_c);
  unlink(name_c);
  free(name_c);
}

void FC_FUNC_(oct_getcwd, OCT_GETCWD)
  (STR_F_TYPE name STR_ARG1)
{
  char s[256];
  getcwd(s, 256);
  TO_F_STR1(s, name);
}

void FC_FUNC_(oct_realpath, OCT_REALPATH)
     (STR_F_TYPE fnam, STR_F_TYPE rnam STR_ARG2)
{
  char *fn=NULL, *rn=NULL;
  TO_C_STR1(fnam, fn);
  rn = realpath(fn, NULL);
  free(fn);
  if(rn!=NULL){
    TO_F_STR2(rn, rnam);
  }else{
    TO_F_STR2("", rnam);
  }
  free(rn);
  return;
}

void FC_FUNC_(oct_dirname, OCT_DIRNAME)
     (STR_F_TYPE fnam, STR_F_TYPE dnam STR_ARG2)
{
  char *fn=NULL, *dn=NULL;
  TO_C_STR1(fnam, fn);
  dn = dirname(fn);
  if(dn!=NULL){
    TO_F_STR2(dn, dnam);
  }else{
    TO_F_STR2("", dnam);
  }
  free(fn);
  return;
}

void FC_FUNC_(oct_basename, OCT_BASENAME)
     (STR_F_TYPE fnam, STR_F_TYPE bnam STR_ARG2)
{
  char *fn=NULL, *bn=NULL;
  TO_C_STR1(fnam, fn);
  bn = basename(fn);
  free(fn);
  if(bn!=NULL){
    TO_F_STR2(bn, bnam);
  }else{
    TO_F_STR2("", bnam);
  }
  return;
}

void FC_FUNC_(oct_getenv, OCT_GETENV)
  (STR_F_TYPE var, STR_F_TYPE value STR_ARG2)
{
  char *name_c, *var_c;

  TO_C_STR1(var, name_c);
  var_c = getenv(name_c);
  free(name_c);

  if(var_c != NULL){
    TO_F_STR2(var_c, value);
  }else{
    TO_F_STR2("", value);
  }
}

/* this function gets a string of the form '1-12, 34' and fills
	 array l with the 1 if the number is in the list, or 0 otherwise */
void FC_FUNC_(oct_wfs_list, OCT_WFS_LIST)
		 (STR_F_TYPE str, fint l[16384] STR_ARG1)
{
  int i, i1, i2;
  char c[20], *c1, *str_c, *s;

  TO_C_STR1(str, str_c);
  s = str_c;
  
  /* clear list */
  for(i=0; i<16384; i++)
    l[i] = 0;
  
  while(*s){
    /* get integer */
    for(c1 = c; isdigit(*s) || isspace(*s); s++)
      if(isdigit(*s)) *c1++ = *s;
    *c1 = '\0';
    i1 = atoi(c) - 1;
    
    if(*s == '-'){ /* range */
      s++;
      for(c1 = c; isdigit(*s) || isspace(*s); s++)
	if(isdigit(*s)) *c1++ = *s;
      *c1 = '\0';
      i2 = atoi(c) - 1;
    } else /* single value */
      i2 = i1;
    
    for(i=i1; i<=i2; i++)
      if(i>=0 && i<16384)
				l[i] = 1;

    if(*s) s++;
  }

  free(str_c);
}

/* ------------------------------ from varia.c ------------------------------- */
#include "varia.h"

void FC_FUNC_(oct_progress_bar, OCT_PROGRESS_BAR)
  (fint *a, fint *max)
{
  if (*max == 0) return; /* Skip the bar if the length is 0 */
  progress_bar(*a, *max);
}

/* ------------------------------ some stuff  -------------------------------- */
void FC_FUNC_(oct_gettimeofday, OCT_GETTIMEOFDAY)
  (fint *sec, fint *usec)
{
#ifdef HAVE_GETTIMEOFDAY
  struct timeval tv;
  struct timezone tz;

  gettimeofday(&tv, &tz);

  /* The typecast below should use long. However, this causes incompatibilities
     with Fortran integers. 
     Using int will cause wrong results when tv.tv_sec exceeds INT_MAX=2147483647 */
  *sec  = (int) tv.tv_sec;
  *usec = (int) tv.tv_usec;
/*
  char str[sizeof("HH:MM:SS")];
  time_t local;
  local = tv.tv_sec;
  strftime(str, sizeof(str), "%T", localtime(&local));
  printf("%s.%06ld \n", str, (long) tv.tv_usec);
  printf("%ld.%06ld \n", (long) tv.tv_sec, (long) tv.tv_usec);
*/
#else
  *sec  = 0;
  *usec = 0; 
#endif
}

double FC_FUNC_(oct_clock, OCT_CLOCK)
  ()
{
#ifdef HAVE_GETTIMEOFDAY
  int sec, usec;
  FC_FUNC_(oct_gettimeofday, OCT_GETTIMEOFDAY) (&sec, &usec);
  return sec + 1.0e-6*usec;
#else
  return (double) clock()/CLOCKS_PER_SEC;
#endif
}

void FC_FUNC_(oct_nanosleep, OCT_NANOSLEEP)
	(fint *sec, fint *nsec)
{
#ifdef HAVE_NANOSLEEP
  /* Datatypes should be long instead of int (see comment in gettimeofday) */
  struct timespec req;
  req.tv_sec  = (time_t) *sec;
  req.tv_nsec = (long)   *nsec;
  nanosleep(&req, NULL);
#endif
}


void FC_FUNC_(oct_sysname, OCT_SYSNAME)
		 (STR_F_TYPE name STR_ARG1)
{
  char *name_c;
  
  sysname(&name_c);
  TO_F_STR1(name_c, name);
  free(name_c);
}


int FC_FUNC_(oct_number_of_lines, OCT_NUMBER_OF_LINES)
  (STR_F_TYPE name STR_ARG1)
{

  FILE *pf;
  int c, i;
  char *name_c;

  TO_C_STR1(name, name_c);
  pf = fopen(name_c, "r");
  free(name_c);
  
  if (pf != NULL) {
    i = 0;
    while ((c = getc(pf)) != EOF) {
      if (c == '\n') i++;
    }
    fclose(pf); 
    return i;
  }else{
    return -1;
  }
}


/* Given a string in C, it breaks it line by line and returns each 
   as a Fortran string. Returns 0 if string does not have more lines. 
*/
void FC_FUNC_(oct_break_c_string, OCT_BREAK_C_STRING)
  (char **str, char **s, STR_F_TYPE line_f STR_ARG1)
{
  char *c, line[256]; /* hopefully no line is longer than 256 characters ;) */

  if(*s == NULL) *s = *str;

  if(*s == NULL || **s == '\0'){
    *s = (char *)(0);
    return;
  }

  for(c=line; **s!='\0' && **s!='\n'; (*s)++, c++)
    *c = **s;
  *c = '\0';
  if(**s=='\n') (*s)++;

  TO_F_STR1(line, line_f);
}

/* 

This function searches in directory given by dirname all files that have the following name:

*_<real_number>_<integer>*

It returns the value of <real_number> found that is closest to freq (or abs(freq))
and for which the value of <integer> matches with the tag argument.

The value found is returned in the freq argument.

ierr results:
0 : value found
1 : no matching file found
2 : cannot open the directory or function not available
 
*/

void FC_FUNC_(oct_search_file_lr, OCT_SEARCH_FILE_LR)
     (double * freq, const fint * tag, fint * ierr, STR_F_TYPE dirname STR_ARG1)
{
#if defined(HAVE_DIRENT_H) && defined(HAVE_CLOSEDIR) && defined(HAVE_READDIR) && defined(HAVE_STRCHR) && defined(HAVE_STRTOD)

  DIR * dir;
  struct dirent * ent;
  char * name_c;
  char * num_start, * num_end;
  double read_value, min;
  int found_something, read_tag;

  TO_C_STR1(dirname, name_c);
  dir = opendir(name_c);

  if(dir == NULL){
    *ierr = 2;
    return;
  }
  free(name_c);

  ent = NULL;
  found_something = 0;

  while(1) {
    ent = readdir(dir);
    if( ent == NULL ) break;

    num_start = strchr(ent -> d_name, '_');

    if(num_start != NULL) {
      num_start++; /* now this points to the beginning of the number */

      /* take the numerical value from the string */
      read_value = strtod(num_start, &num_end);
      
      if ( num_end == num_start ) continue; /* no number found */

      /* check that we have the correct tag */
      if(num_end[0] == '_') {

	num_start = num_end + 1;
	read_tag = (int) strtol(num_start, &num_end, 10);
	if ( num_end == num_start ) continue; /* no tag found */
	if ( read_tag != *tag ) continue; /* tag does not match */
 
      } else continue;


      /* if this is the first number we found */
      if ( !found_something ) {
	min = read_value;
	found_something = 1;
      } else if (fabs(fabs(min)-fabs(*freq)) > fabs(fabs(read_value)-fabs(*freq)) ) {
      /* if the value is closer than previous */
	min = read_value;
      }
    }

  }

  closedir(dir);

  if(found_something) {
    *ierr = 0;
    *freq = min ;
  } else {
    *ierr = 1;
  }

#else
#warning directory search not compiled
  fprintf(stderr, "Warning: Directory search not available since certain C functions are not available.\n");
  *ierr = 2;
#endif

}

void * FC_FUNC_(oct_get_memory_usage, OCT_GET_MEMORY_USAGE)()
{
#if defined(HAVE_SYSCONF) && defined(HAVE_GETPID)
  static size_t pagesize = 0;
  FILE *f;
  int pid;
  unsigned long mem;
  char s[256];
  
  if(pagesize == 0)
    pagesize = sysconf(_SC_PAGESIZE);
  
  pid = getpid();
  sprintf(s, "%s%d%s", "/proc/", pid, "/statm");
  if((f = fopen(s, "r")) == (FILE *)NULL) return (void *)(-1);
  fscanf(f, "%lu", &mem);
  fclose(f);
  
  return (void *) (mem*pagesize);
#elif HAVE_SBRK
  return sbrk(0);
#else
  return 0;
#endif
}

void FC_FUNC_(oct_exit_failure, OCT_EXIT_FAILURE)(){
  exit(EXIT_FAILURE);
}
