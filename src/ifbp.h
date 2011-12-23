/*

 Copyright (C) 2011 Sun Ruping <ruping@molgen.mpg.de>

 This file is part of Breakpointer.

 Breakpointer is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License.

*/

#include <cstdio>
#include <getopt.h>
#include <cstdlib>
#include <cstring>

struct parameters {
  char* mapping_f;
  unsigned int windowsize;
  unsigned int readlen;
  unsigned int unique;
  unsigned int indiprint;
  char* tag_uniq;
  unsigned int val_uniq;
};

struct parameters* interface(struct parameters* param,int argc, char *argv[]);
void delete_param(struct parameters* param);
void usage(void);

const char* program_name;

struct parameters* interface(struct parameters* param, int argc, char *argv[]){

  program_name = argv[0];
  int c;     // the next argument
  int help = 0;

  if (argc < 2){
    usage();
    exit(0);
  }

  param = new struct parameters;
  param->mapping_f = new char;
  param->tag_uniq = new char;
  //param->val_uniq = new char;

  const struct option long_options[] ={
    {"unique",0,0,'u'},
    {"readlen",1,0,'l'},
    {"mapping",1,0,'m'},
    {"windowsize",1,0,'w'},
    {"indiprint",0,0,'i'},
    {"tag_uniq",1,0,'t'},
    {"val_uniq",1,0,'v'},
    {"help",0,0,'h'},
    {0,0,0,0}
  };

  while (1) {

    int option_index = 0;
    c = getopt_long_only (argc, argv,"hium:w:l:v:t:",long_options, &option_index);

    if (c == -1){
      break;
    }

    switch(c) {
    case 0:
      break;
    case 'u':
      param->unique = 1;
      break;
    case 'm':
      param->mapping_f = optarg;
      break;
    case 'w':
      param->windowsize = atoi(optarg);
      break;
    case 'l':
      param->readlen = atoi(optarg);
      break;
    case 't':
      param->tag_uniq = optarg;
      break;
    case 'v':
      param->val_uniq = atoi(optarg);
      break;
    case 'i':
      param->indiprint = 1;
      break;
    case 'h':
      help = 1;
      break;
    case '?':
      help = 1;
      break;
    default:
      help = 1;
      break;
    }
  }

  if(help){
    usage();
    delete_param(param);
    exit(0);
  }

  return param;
}

void usage()
{
  fprintf(stdout, "\nbreakpointer (search for depth skewed regions) @ BreakPointer v0.1 2011 Sun Ruping <ruping@molgen.mpg.de>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "Usage: %s options >output\n\n", program_name);
  fprintf(stdout, "-m --mapping     <string> A BAM alignment file or a file containing the filenames of multiple BAM files (one file per line). MUST be according to the chromosome and start position.\n                          In case of multiple BAM files, make sure these BAM files are sorted samely (require header tag: \"@HD\tVN:1.0\tSO:coordinate\").\n");
  fprintf(stdout, "-w --windowsize  <int>    the size in bp of the sliding window (default: 10 for reads < 50bp, 20 for longer and variable read length).\n");
  fprintf(stdout, "-l --readlen     <int>    the size in bp of the read length (default: allowing variable read length).\n");
  fprintf(stdout, "-u --unique               take only uniquelly mapped reads (default: take all mapped reads).\n                          since different mappers generate different tags for uniqueness, if -q is set, user shoule provide unique tag info (see tag/val_uniq). \n                          we recommand not to set this option if the mapping file only contain a few multiple location reads, in case users are not sure about the unique tags\n");
  fprintf(stdout, "-t --tag_uniq    <string> the tag in the bam file denotating whether a read is uniquely mapped (default \"XT\" is taken as from BWA).\n");
  fprintf(stdout, "-v --val_uniq    <int>    the value for the above tag of uniquely mapped reads (default value is taken as from the output from BWA).\n");
  fprintf(stdout, "-h --help                 print the help message.\n");
  fprintf(stdout, "\n");
}

void delete_param(struct parameters* param)
{
  delete(param->mapping_f);
  delete(param->tag_uniq);
  delete(param);
}
