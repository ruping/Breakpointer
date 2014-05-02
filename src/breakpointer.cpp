/*****************************************************************************

  breakpointer.cpp @ Breakpointer
  using a sliding window to get region with skewed depth toward ending reads.

  (c) 2011 - Sun Ruping
  Dept. Vingron (Computational Mol. Bio.)
  Max-Planck-Institute for Molecular Genetics
  Ihnestr. 73, D-14195, Berlin, Germany

  current: Department of Systems Biology, New York, USA
  EMAIL: rs3412@columbia.edu

  Breakpointer is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License.

******************************************************************************/


#include <api/BamReader.h>
#include <api/BamMultiReader.h>
#include "mathstats.h"
#include "ifbp.h"
using namespace BamTools;

#include <cstring>
#include <vector>
#include <deque>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
using namespace std; 

//merged print
unsigned int indipr    = 0;
string last_chr        = "SRP";
unsigned int last_end  = 0;
unsigned int ol_start  = 0;
unsigned int ol_end    = 0;
unsigned int ol_dis    = 0;
float ol_number = 0;
float ol_depth  = 0;
float ol_ratio1 = 0;
float ol_ratio2 = 0;
float ol_score  = 0;

//for window storage
struct bucket {
  unsigned int bdepth;  //the depth of this window
  unsigned int bdeps;   //the start depth
  unsigned int bdepe;   //the end depth 
};

struct window {
  map <unsigned int, struct bucket> buckets;
  unsigned int end;    //the end of the window = start  + winsize - 1 we don't need it here
  unsigned int depth;  //the depth of this window
  unsigned int deps;   //the start depth
  unsigned int depe;   //the end depth 
};

//for read storage
struct read {
  unsigned int start;  // start of the read
  unsigned int end;    // end of the read
};

inline void ParseCigar(const vector<CigarOp> &cigar, vector<int> &blockStarts, vector<int> &blockEnds, unsigned int &alignmentEnd); //deprecated
inline void print_endepth(const string &chr, unsigned int winstart, const float &winsize, struct window &window, const float &prob);
inline void splitstring(const string &str, vector<string> &elements, const string &delimiter);
inline string int2str(unsigned int &);

int main (int argc, char **argv){

  struct parameters *param = 0;
  param = interface(param, argc, argv);
   
  // check the arguments

  unsigned int read_length = 0;
  if ( param->readlen ) read_length = param->readlen;  //argument readlength
  else cerr << "no readlength argument is given, using variable read length setting" << endl;

  unsigned int windowsize;
  if ( param->windowsize ) windowsize = param->windowsize;   //argument windowsize
  else {
    if (read_length <= 50 && read_length != 0 ) windowsize = 10;
    if (read_length > 50)                       windowsize = 20;
    if (read_length == 0)                       windowsize = 20;
  }
  cerr << "windowsize is: " << windowsize << endl;

  indipr = param->indiprint;           // argument print indi

  float winsize = windowsize;
  float readlen = read_length;
  float prob = 0.;
  if (read_length != 0) {
    prob  = (2 * winsize) / (winsize + readlen);
    cerr << "binomial prob: " << prob << endl;
  }
  else {
    cerr << "binomial probs will be decided for each length group" << endl;
  }

  //string tag_uniq = param->tag_uniq;
  //unsigned int val_uniq;
  //if (tag_uniq == "") tag_uniq = "XT";  //default for BWA alignment
  //if (param->val_uniq) val_uniq = param->val_uniq;
  //else val_uniq = 85;                   //default for BWA alignment
  //cerr << "unique tag is: " << tag_uniq << "\t" << val_uniq << endl; 
  
  string oldchr;                        //for checking the chromosome
  unsigned int oldstart = 0;            //compare start (piling up)

  map <unsigned int, struct window> windows; //MAP container of windows
  deque <struct read> reads;            //DEQUE container of reads
  set <string> pileup;                  //SET container of piling-up reads

//-------------------------------------------------------------------------------------------------------+
// BAM input (file or filenames?)                                                                        |
//-------------------------------------------------------------------------------------------------------+
  char *fof = param->mapping_f;
  FILE *IN=NULL;
  char line[5000];
  int filecount=0;
  vector <string> fnames;
 
  if (strchr(fof,' ')!=NULL) {
    char *ptr;
    ptr=strtok(fof," ");
    while (ptr!=NULL) {
      fnames.push_back(ptr);
      filecount++;
      ptr=strtok(NULL," ");
    }
  } else {
    IN=fopen(fof,"rt");
    if (IN!=NULL) {
      long linecount=0;
      while (fgets(line,5000-1,IN)!=NULL) {
        linecount++;
        if (line[0]!='#' && line[0]!='\n') {
          char *ptr=strchr(line,'\n');
          if (ptr!=NULL && ptr[0]=='\n') {
            ptr[0]='\0';
          }
          FILE *dummy=NULL;
          dummy=fopen(line,"rt");
          if (dummy!=NULL) {     // seems to be a file of filenames...
            fclose(dummy);
            fnames.push_back(line);
            filecount++;
          } else if (filecount==0 || linecount>=1000-1) {  // seems to be a single file
            fnames.push_back(fof);
            filecount++;
            break;
          }
        }
      }
      fclose(IN);
    }
  }  //file or file name decided and stored in vector "fnames"

  cerr << "the input mapping files are:" << endl;
  vector <string>::iterator fit = fnames.begin();
  for(; fit != fnames.end(); fit++) {
    cerr << *fit << endl;
  }

//-------------------------------------------------------------------------------------------------------+
// end of file or filenames                                                                              |
//-------------------------------------------------------------------------------------------------------+

  // open the BAM file(s)  
  BamMultiReader reader;
  reader.Open(fnames);
  
  // get header & reference information
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();

  if ( ! reader.LocateIndexes() )     // opens any existing index files that match our BAM files
     reader.CreateIndexes();         // creates index files for BAM files that still lack one

  BamAlignment bam;

  while (reader.GetNextAlignment(bam)) {  //getting each alignment

    if (bam.IsMapped() == false) continue; //skip unaligned reads

    string read_qual =  bam.Qualities;
    unsigned int real_length = read_qual.size();

    if (read_length != 0) {     // the length is preset
      if (real_length != read_length) //skip the read with different length
         continue;
    }

    if (refs.at(bam.RefID).RefName != oldchr && !windows.empty()) {  //a new chr, the windows should be printed out and then clean up

      map <unsigned int,window>::iterator iter = windows.begin();
      for (; iter != windows.end() ; iter++) { //print the last windows of the old chr
        print_endepth(oldchr, (*iter).first, winsize, (*iter).second, prob);
      }

      windows.clear(); //clear windows
      reads.clear();   //clear reads

      oldstart = 0;
    }
   

    unsigned int unique = 0;
    if ( bam.HasTag("NH") ) {
      bam.GetTag("NH", unique);                   // uniqueness
    } else if (bam.HasTag("XT")) {
      string xt;
      bam.GetTag("XT", xt);                       // bwa aligner
      xt = xt.substr(0,1);
      if (xt != "R") {
        unique = 1;
      }
    } else {
      if (bam.MapQuality > 10) {                   // bowtie2
        unique = 1;
      }
    }

    if (param->unique == 1) {
      if (unique != 1) {                         // skipe uniquelly mapped reads 
        continue;
      }
    }

    string strand = "+";
    if ( bam.IsReverseStrand() ) strand = "-";

    unsigned int alignmentStart, alignmentEnd;
    //CIGAR string and figure out the alignment blocks
    //vector<int> blockLengths;
    //vector<int> blockStarts;
    //blockStarts.push_back(0);
    //extract the block starts and lengths from the CIGAR string
    //ParseCigar(bam.CigarData, blockStarts, blockLengths, alignmentEnd);
    alignmentStart = bam.Position+1;
    alignmentEnd   = bam.GetEndPosition();

    //if (alignmentEnd_new != alignmentEnd) cerr << bam.Name << "\t" << alignmentStart << "\t" << alignmentEnd << "\t" << alignmentEnd_new  << endl;

    //skip piling up reads (pileup based no starts+ends+strand)
    string alignSum = int2str(alignmentStart) + int2str(alignmentEnd) + strand;
    if (alignmentStart != oldstart) {
      pileup.clear();          //clear pileup set
      pileup.insert(alignSum); //insert the new read
    }
    else if (alignmentStart == oldstart){
      if   (pileup.count(alignSum)) continue; //if find pileup read
      else pileup.insert(alignSum); //insert
    }   

    //insert two ends into the windows map, default depths are 0
    bucket tmpb = {0, 0, 0};  
    map <unsigned int, struct bucket> tmpm;
    tmpm.insert( pair <unsigned int, struct bucket> (real_length, tmpb));
    window tmpw1 = {tmpm, alignmentStart+windowsize-1, 0, 0, 0};    
    window tmpw2 = {tmpm, alignmentEnd+windowsize-1,   0, 0, 0};
    windows.insert( pair < unsigned int, struct window > (alignmentStart, tmpw1) );
    windows.insert( pair < unsigned int, struct window > (alignmentEnd, tmpw2) );

    //add the current read into the reads deque
    struct read tmp3 = {alignmentStart, alignmentEnd}; 
    reads.push_back(tmp3);

    map <unsigned int, struct window>::iterator iter = windows.begin();

    while (iter != windows.end()) { //iterate the windows

      if ((iter->second).depth != 0) {  //old windows, only compare with the new read

        if ((iter->second).end < alignmentStart){ //the window is beyond the new read start, print the window and delete it
          print_endepth(oldchr, (*iter).first, winsize, (*iter).second, prob);
          windows.erase(iter++);
          continue;
        }

        if ((iter->second).end >= alignmentStart && iter->first <= alignmentEnd){
          (iter->second).depth++;    //depth++
          if ((iter->second).buckets.count(real_length) > 0)
            ((iter->second).buckets)[real_length].bdepth++;
          else {
            ((iter->second).buckets).insert(pair <unsigned int, struct bucket> (real_length, tmpb));
            ((iter->second).buckets)[real_length].bdepth++;
          }
          if (iter->first <= alignmentStart){
             (iter->second).deps++;  //depth_start++
             ((iter->second).buckets)[real_length].bdeps++;
          }
          if ((iter->second).end >= alignmentEnd){
             (iter->second).depe++;  //depth_end++
             ((iter->second).buckets)[real_length].bdepe++;
          }
        }

      } //old window

      else if ((iter->second).depth == 0){ //new window, check all the current reads
         
        deque <struct read>::iterator iter2 = reads.begin();

        for(;iter2 != reads.end();){ // loop over the deque of reads 

          if (iter2->end < alignmentStart){ //read beyond the current alignment start, the read should be deleted
            iter2 = reads.erase(iter2);
            continue;
          }

          if (iter2->end >= iter->first && iter2->start <= (iter->second).end){
            (iter->second).depth++;    //depth++
            if ((iter->second).buckets.count(real_length) > 0)
              ((iter->second).buckets)[real_length].bdepth++;
            else {
              ((iter->second).buckets).insert(pair <unsigned int, struct bucket> (real_length, tmpb));
              ((iter->second).buckets)[real_length].bdepth++;
            }
            if (iter2->end <= (iter->second).end){
               (iter->second).depe++;  //depth_end++
               ((iter->second).buckets)[real_length].bdepe++;
            }
            if (iter2->start >= iter->first){
               (iter->second).deps++;  //depth_start++
               ((iter->second).buckets)[real_length].bdeps++;
            } 
          } //overlap

          iter2++;

        } //loop the reads

      } //new window

      iter++; //increment of windows pointer

    } //iterate the windows

    oldchr    = refs.at(bam.RefID).RefName;
    oldstart  = alignmentStart;

  } //getting each alignment

  reader.Close();
  

  //print out the windows in the pool
  if (!windows.empty()) {
    map <unsigned int,window>::iterator iter = windows.begin();
    for (; iter != windows.end() ; iter++){ //print the last windows of the old chr
      print_endepth(oldchr, (*iter).first, winsize, (*iter).second, prob);
    }
    windows.clear(); //clear windows
    reads.clear();   //clear reads
  }

  //print the last guy
  if (last_chr != "SRP") {
    ol_dis = ol_end - ol_start + 1;
    float av_depth  = ol_depth/ol_number;
    float av_ratio1 = ol_ratio1/ol_number;
    float av_ratio2 = ol_ratio2/ol_number;
    float av_score  = ol_score/ol_number;
    printf("%s\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\n", last_chr.c_str(), ol_start, ol_end,
           ol_dis, av_depth, av_ratio1, av_ratio2, av_score);
    last_chr = "SRP";  // ending
  }

  cerr << "step1 of @Breakpointer done." << endl;

}

inline string int2str(unsigned int &i){
  string s;
  stringstream ss(s);
  ss << i;
  return ss.str();
}

inline void splitstring(const string &str, vector<string> &elements, const string &delimiter) {
  string::size_type lastPos = str.find_first_not_of(delimiter, 0);
  string::size_type pos     = str.find_first_of(delimiter, lastPos);

  while (string::npos != pos || string::npos != lastPos) {
    elements.push_back(str.substr(lastPos, pos - lastPos));
    lastPos = str.find_first_not_of(delimiter, pos);
    pos = str.find_first_of(delimiter, lastPos);
  }
}

inline void ParseCigar(const vector<CigarOp> &cigar, vector<int> &blockStarts, vector<int> &blockLengths, unsigned int &alignmentEnd) {

  int currPosition = 0;
  int blockLength  = 0;

  //  Rip through the CIGAR ops and figure out if there is more
  //  than one block for this alignment
  vector<CigarOp>::const_iterator cigItr = cigar.begin();
  vector<CigarOp>::const_iterator cigEnd = cigar.end();
  for (; cigItr != cigEnd; ++cigItr) {
    switch (cigItr->Type) {
    case ('M') :
      blockLength  += cigItr->Length;
      currPosition += cigItr->Length;
    case ('I') : break;
    case ('S') : break;
    case ('D') :
      blockLength  += cigItr->Length;
      currPosition += cigItr->Length;
      break;
    case ('P') : break;
    case ('N') :
      blockStarts.push_back(currPosition + cigItr->Length);
      blockLengths.push_back(blockLength);
      currPosition += cigItr->Length;
      blockLength = 0;
      break;
    case ('H') : break;                             // for 'H' - do nothing, move to next op
    default    :
      printf("ERROR: Invalid Cigar op type\n");   // shouldn't get here
      exit(1);
    }
  }
  // add the kast block and set the
  // alignment end (i.e., relative to the start)
  blockLengths.push_back(blockLength);
  alignmentEnd = currPosition;
}

inline void print_endepth(const string &chr, unsigned int winstart, const float &winsize, struct window &window, const float &prob){

  float depth  = window.depth;
  float starts = window.deps;
  float ends   = window.depe;
  float ratio1 = (starts+ends)/depth;
  float ratio2 = starts/(starts+ends);

  if (depth > 1 && ratio1 > prob) {

    double score;

    if (prob != 0.) {  // single
      score = pbinom((starts+ends), depth, prob, 0);
      score = -log10(score);
    }
    else {             // multiple
      float bcs = 0.;
      map <unsigned int, struct bucket>::iterator bit = window.buckets.begin();
      for (; bit != window.buckets.end(); bit++) {
        if ((bit->second).bdepth < 2) continue;
        else {
          float bdp = (bit->second).bdepth;
          float bds = (bit->second).bdeps;
          float bde = (bit->second).bdepe;
          float bprob = (2 * winsize) / (winsize + (bit->first));
          float bscore = (bdp/depth)*pbinom((bds+bde), bdp, bprob, 0);
          bcs += bscore;
        }
      }
      if (bcs > 0.) score = -log10(bcs);
      else score = 0.;
      //cerr << winstart << "\t" << window.buckets.size() << "\t" << score << endl;
    }

    if (score > 1.) { // merge and print      
 
      if (indipr == 1) {
        fprintf(stderr, "%s\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.5f\n", chr.c_str(), winstart, window.end,
                window.depth, window.deps, window.depe, ratio1, ratio2, score);
      }

     
      if (chr != last_chr) {
        
        if (last_chr != "SRP") {
          ol_dis = ol_end - ol_start + 1;
          float av_depth  = ol_depth/ol_number;
          float av_ratio1 = ol_ratio1/ol_number;
          float av_ratio2 = ol_ratio2/ol_number;
          float av_score  = ol_score/ol_number;
          printf("%s\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\n", last_chr.c_str(), ol_start, ol_end,
                 ol_dis, av_depth, av_ratio1, av_ratio2, av_score);
        }
 
        //reset everything
        last_end  = 0;
        ol_start  = 0;
        ol_end    = 0;
        ol_dis    = 0;
        ol_number = 0;
        ol_depth  = 0;
        ol_ratio1 = 0;
        ol_ratio2 = 0;
        ol_score  = 0;
      }

      if (winstart <= last_end){ //overlapping : put this window into vector
        ol_end     = window.end;
        ol_depth  += window.depth;
        ol_ratio1 += ratio1;
        ol_ratio2 += ratio2;
        ol_score  += score;
        ol_number += 1;
      }
     
      if (winstart > last_end){  //non-overlapping

        if (last_end != 0){
          ol_dis = ol_end - ol_start + 1;
          float av_depth  = ol_depth/ol_number;
          float av_ratio1 = ol_ratio1/ol_number;
          float av_ratio2 = ol_ratio2/ol_number;
          float av_score  = ol_score/ol_number;
          printf("%s\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\n", chr.c_str(), ol_start, ol_end,
                 ol_dis, av_depth, av_ratio1, av_ratio2, av_score);
        }

        // reset ol
        ol_start  = winstart;
        ol_end    = window.end;
        ol_dis    = 0;
        ol_number = 1;
        ol_depth  = window.depth;
        ol_ratio1 = ratio1;
        ol_ratio2 = ratio2;
        ol_score  = score;

      }
     
      last_chr = chr;
      last_end = window.end;

    } //merging
  }
}
