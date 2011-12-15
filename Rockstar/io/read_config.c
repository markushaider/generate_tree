#undef _POSIX_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include "read_config.h"

void add_to_string_array(char ***array, char *string, size_t length, int num_entries) {
  char **elem;
  if (!(num_entries % 10)) {
    *array = (char **) realloc(*array, sizeof(char *)*(num_entries+10));
    if (!(*array)) {
      fprintf(stderr, "[Error] Couldn't allocate memory for keys/values!\n");
      exit(1);
    }
  }

  elem = (*array)+num_entries;
  *elem = (char *)malloc(sizeof(char)*(length+1));
  if (!(*elem)) {
    fprintf(stderr, "[Error] Couldn't allocate memory for keys/values!\n");
    exit(1);
  }

  memcpy(*elem, string, length*sizeof(char));
  (*elem)[length] = 0;
}

inline static char *_config_to_string(struct configfile *c, char *key) {
  int i;
  for (i=0; i<c->num_entries; i++)
    if (strcmp(key,c->keys[i])==0) return c->values[i];
  return 0;
}

char *config_to_string(struct configfile *c, char *key, char *def_val) {
  char *val = _config_to_string(c, key);
  if (!val) {
    add_to_string_array(&(c->keys), key, strlen(key), c->num_entries);
    add_to_string_array(&(c->values), def_val, strlen(def_val), c->num_entries);
    c->num_entries++;
    return strdup(def_val);
  }
  return strdup(val);
}

double config_to_real(struct configfile *c, char *key, double def_val) {
  char buffer[100];
  char *val = _config_to_string(c, key);
  if (!val) {
    sprintf(buffer, "%.10e", def_val);
    add_to_string_array(&(c->keys), key, strlen(key), c->num_entries);
    add_to_string_array(&(c->values), buffer, strlen(buffer), c->num_entries);
    c->num_entries++;
    return def_val;
  }
  return strtod(val, (char **)NULL);
}

void config_to_real3(struct configfile *c, char *key, double *res, char *def_val) {
  char *val = _config_to_string(c, key), *cur, *endp;
  int64_t i, length;
  if (!val) {
    add_to_string_array(&(c->keys), key, strlen(key), c->num_entries);
    add_to_string_array(&(c->values), def_val, strlen(def_val), c->num_entries);
    c->num_entries++;
    val = def_val;
  }
  length = strlen(val);
  cur = val;
  for (i=0; i<3; i++) {
    while ((cur < val+length) && !isdigit(*cur)) cur++;
    if (cur >= (val+length)) { res[i] = 0; continue; }
    res[i] = strtod(cur, &endp);
    cur = endp;
  }
}


inline static void trim(char **a, char **b) {
  if ((*a)>(*b)) return;
  while ((*a) <= (*b)) {
    if (!(((**a)==' ') || ((**a)=='\t'))) {
      if (((**a)=='\'') || ((**a)=='\"')) 
	(*a)+=1;
      break;
    }   
    (*a)+=1;
  }
  while ((*a) <= (*b)) {
    if (!(((**b)==' ') || ((**b)=='\t') || ((**b)=='\n') || ((**b)=='\0'))) {
      if (((**b)=='\'') || ((**b)=='\"'))
	(*b)-=1;
      break;
    }
    (*b)-=1;
  }
}

void load_config(struct configfile *c, char *filename) {
  char buffer[1024];
  FILE *input;
  char *key_start, *key_end;
  char *val_start, *val_end;

  c->keys = c->values = 0;
  input = fopen(filename, "r");
  if (!input) {
    c->num_entries = 0;
    return;
  }
  
  while (fgets(buffer, 1024, input)) {
    val_start = buffer;
    strsep(&val_start, "#"); //Ignore comments;
    key_start = val_start = buffer;
    val_end = buffer + strlen(buffer);
    if (val_end == key_start) continue; //Blank line
    strsep(&val_start, "=");
    if (!val_start) continue; //No = sign found.
    key_end = val_start - 1;
    trim(&key_start, &key_end);
    trim(&val_start, &val_end);
    if (key_end < key_start) continue; //No key found.
    if (val_end < val_start) continue; //No value found.

    add_to_string_array(&(c->keys), key_start, key_end-key_start+1, c->num_entries);
    add_to_string_array(&(c->values), val_start, val_end-val_start+1, c->num_entries);
    c->num_entries++;
  }
  fclose(input);
}

void write_config(struct configfile c, char *filename) {
  FILE *output;
  int i;
  output = fopen(filename, "w");
  if (!output) {
    fprintf(stderr, "[Error] Couldn't open file %s for writing!\n", filename);
    exit(1);
  }
  
  for (i=0; i<c.num_entries; i++)
    fprintf(output, "\"%s\" = \"%s\"\n", c.keys[i], c.values[i]);
  fclose(output);
}

void free_config(struct configfile c) {
  int i;
  for (i=0; i<c.num_entries; i++) {
    free(c.keys[i]);
    free(c.values[i]);
  }
  free(c.keys);
  free(c.values);
  c.num_entries = 0;
  c.keys = c.values = 0;
}

