/*****************************************************************************

  breakmis.cpp @ Breakpointer
  Mismatch screening for each merged region obtained using breakpointer.

  (c) 2011 - Sun Ruping
  Dept. Vingron (Computational Mol. Bio.)
  Max-Planck-Institute for Molecular Genetics
  Ihnestr. 73, D-14195, Berlin, Germany
  EMAIL: ruping@molgen.mpg.de

  Breakpointer is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License.

******************************************************************************/

#include <api/BamReader.h>
#include <api/BamMultiReader.h>
using namespace BamTools;

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <deque>
#include <map>
#include <set>
#include <string>
#include <cstring>
#include <sstream>
#include <iomanip>
#include "ifbm.h"
#include "mathstats.h"

using namespace std;

struct SVread {
  string name;
  unsigned int start;
  unsigned int end;
  string strand;
  string seq;
  set <unsigned int> me;
};

struct region {
  string chr;
  unsigned int start;
  unsigned int end;
  unsigned int dis;
  float depth;
  float ratio1;
  float ratio2;
  float score;
  unsigned int coverage;
  unsigned int mismatch;
  map <unsigned int, unsigned int> mispos;
  set <unsigned int> forbid;
  map <unsigned int, unsigned int> posendcov;
  vector <struct SVread> SVreads;
};


inline void ParseCigar(const vector<CigarOp> &cigar, vector<int> &blockStarts, vector<int> &blockEnds, unsigned int &alignmentEnd);
inline void splitstring(const string &str, vector<string> &elements, const string &delimiter);
inline void splitgfftag(const string &str, map <string, string> &elements, const string &delimiter, vector<string> &tag_want);
inline void eatline(const string &str, deque <struct region> &regions_ref);
inline string int2str(unsigned int &i);
inline string flo2str(float &f);
inline void print_mismatch(struct region &region);
inline void finished(const unsigned int &where);

int main (int argc, char *argv[]) {
 
  struct parameters *param = 0;
  param = interface(param, argc, argv);

  unsigned int readlen = 0;
  if ( param->readlen ) readlen = param->readlen;  //argument readlength
  else cerr << "no readlength argument is given, using variable read length setting" << endl;

  unsigned int endlen;
  if (readlen == 0) endlen = 10;
  else if (readlen < 50) endlen = 10;
  else if (readlen >=50 && readlen <= 100) endlen = 15;
  else endlen = 20;

  string qual_clip = param->qual_clip;
  if (qual_clip == "no")  cerr << "quality clipping is turned off." << endl;
  else if (qual_clip == "phred33" || qual_clip == "phred64" || qual_clip == "solexa64")  cerr << "read quality type for clipping is set to: " << qual_clip << ".\n"; 
  else{
    qual_clip = "phred33";
    cerr << "quality clipping is on, quality taken default: " << qual_clip << ".\n";
  }

  unsigned int unique = param->unique;  // argument unique
  string tag_uniq = param->tag_uniq;
  unsigned int val_uniq;
  if (tag_uniq == "") tag_uniq = "XT"; //default for BWA alignment
  if (param->val_uniq) val_uniq = param->val_uniq;
  else val_uniq = 85;   //default for BWA alignment
  cerr << "unique tag is: " << tag_uniq << "\t" << val_uniq << endl;

  string mistag = param->mistag;  // tag for mismatch
  if (mistag == "") mistag = "MD";
  cerr << "mismatch tag in the bam file is: " << mistag << endl;

  string old_chr = "SRP";
  unsigned int oldstart = 0;  
  map <string, unsigned int> pileup;             //SET container of piling-up reads
 
//-------------------------------------------------------------------------------------------------------+
// end of file or filenames                                                                              |
//-------------------------------------------------------------------------------------------------------+
  char *fof = param->mapping_f;
  FILE *IN=NULL;
  char fline[5000];
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
      long flinecount=0;
      while (fgets(fline,5000-1,IN)!=NULL) {
        flinecount++;
        if (fline[0]!='#' && fline[0]!='\n') {
          char *ptr=strchr(fline,'\n');
          if (ptr!=NULL && ptr[0]=='\n') {
            ptr[0]='\0';
          }
          FILE *dummy=NULL;
          dummy=fopen(fline,"rt");
          if (dummy!=NULL) {     // seems to be a file of filenames...
            fclose(dummy);
            fnames.push_back(fline);
            filecount++;
          } else if (filecount==0 || flinecount>=1000-1) {  // seems to be a single file
            fnames.push_back(fof);
            filecount++;
            break;
          }
        }
      }
      fclose(IN);
    }
  } //file or file name decided and stored in vector "fnames"

  cerr << "the input mapping files are:" << endl;
  vector <string>::iterator fit = fnames.begin();
  for(; fit != fnames.end(); fit++){
    cerr << *fit << endl;
  }

//-------------------------------------------------------------------------------------------------------+
// end of file or filenames                                                                              |
//-------------------------------------------------------------------------------------------------------+

  BamMultiReader reader;
  reader.Open(fnames);
  BamAlignment bam;

  // get header & reference information
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();

  if ( !reader.LocateIndexes() )     // opens any existing index files that match our BAM files
     reader.CreateIndexes();         // creates index files for BAM files that still lack one

  //regions for the input of region file
  ifstream region_f;
  string line;
  region_f.open(param->region_f, ios_base::in);  // the file is opened

  deque <struct region> regions;
  getline(region_f, line); //get the first line
  eatline(line,regions);

  deque <struct region>::iterator it = regions.begin();

  while ( it->chr != old_chr ) {     // a new chr come from the region file

    old_chr = it->chr;  // set the current chr as old chr

    int chr_id  = reader.GetReferenceID(it->chr);

    if (chr_id == -1) {  //reference not found

      for (; it != regions.end() && it->chr == old_chr; ) {
        print_mismatch(*it);       // print the old region info
        it = regions.erase(it);             // erase the current region
      }

      while ( regions.empty() ) {
        getline(region_f, line);
        if ( region_f.eof() ) finished(1);
        eatline(line, regions);
        it = regions.begin();
        if (it->chr == old_chr){
          print_mismatch(*it);
          regions.clear();
          continue;
        }
      }
      continue;
    }

    // set to new chr
    int chr_len = refs.at(chr_id).RefLength;
    if ( !reader.SetRegion(chr_id, 1, chr_id, chr_len) ) {
        cerr << "bamtools count ERROR: Jump region failed " << it->chr << endl;
        reader.Close();
        exit(1);
    }

    while (reader.GetNextAlignment(bam)) {    //reading each alignment

      if (bam.IsMapped() == false) continue;  //skip unaligned reads

      string read_qual =  bam.Qualities;
      unsigned int real_length = read_qual.size();

      if (readlen != 0) {     // the length is preset
        if (real_length != readlen) //skip the read with different length
          continue;
      }
 
      string queryB   = bam.QueryBases;
    
      //skip multiple location reads
      if (unique == 1){
        unsigned int XT;
        if (bam.GetTag(tag_uniq, XT)) {            // only for bwa mapping reporting bam file
          if (XT != val_uniq) continue;              // skip non-uniquely mapped reads 
        }
      }

      string strand = "+";
      if (bam.IsReverseStrand()) strand = "-";

      unsigned int alignmentStart, alignmentEnd;
      alignmentStart = bam.Position+1;
      alignmentEnd   = bam.GetEndPosition();
      unsigned int readends1 = alignmentStart+endlen;
      unsigned int readends2 = alignmentEnd-endlen;

      // get clipping infomation
      unsigned int clipleft = 0;
      unsigned int clipright = 0;
      bool clipStatus = false;
      if (qual_clip != "no"){
        clipleft = 10000;
        clipright = 10000;
        unsigned int qpos = 0;
        unsigned int qsize = 0;
        for (; qpos < real_length; qpos++){
          int qual_now;
          if (qual_clip == "phred33")  qual_now = int(read_qual[qpos]) - 33;
          else if (qual_clip == "phred64")  qual_now = int(read_qual[qpos]) - 64;
          else if (qual_clip == "solexa64") qual_now = int(read_qual[qpos]) - 64;
          else qual_now = int(read_qual[qpos]) - 33;
          // cerr << qual_now << " ";        

          if (qual_now >= 5) {
            qsize++;
            if (qsize >= 5) {
              if (clipleft == 10000) {clipleft = (qpos + 1) - qsize;}
              clipright = real_length - (qpos + 1);
            } 
          } //qual_clip above threshold
          else
            qsize = 0;
        }
        if (clipleft == 10000)  clipleft  = real_length;
        if (clipright == 10000) clipright = real_length;      
        if (clipleft != 0 || clipright != 0) clipStatus = true;
      }
      // end: get clipping infomation

      // decode the MD tag to get the mismatches
      string MD;
      vector <string> tagMD;
      vector <unsigned int> mismatch;
      vector <unsigned int> fbpos;
      bool MisStatus = false;
      if (bam.GetTag(mistag, MD)) {
        splitstring(MD, tagMD, "ACGTN^");
        if (tagMD.size() > 1) {
          tagMD.pop_back();
          unsigned int pos = 0;
          vector <string>::iterator mditer = tagMD.begin();
          for (; mditer != tagMD.end(); mditer++) {

            pos += (atoi((*mditer).c_str()) + 1);              //get the position in the alignement

            if (clipStatus == false) {  //clipStatus == false
              if ( pos < endlen || pos > (real_length - endlen + 1) ) {           //mismatch in the ends
                unsigned int mis = pos + (alignmentStart - 1);                // to genomic coordinates
                mismatch.push_back(mis);
                if (MisStatus == false) MisStatus = true;
              }
              else {                                                        //forbid this pos as a real mis pos
                unsigned int forbidpos = pos + (alignmentStart - 1);
                fbpos.push_back(forbidpos);
                if (MisStatus == false) MisStatus = true;
              }
            }  //noclip

            else {                     //clipStatus == true
              if (pos > clipleft && pos < (real_length - clipright + 1)) {       // not in clipped region
                if (pos < endlen || pos > (real_length - endlen + 1)) {
                  unsigned int mis = pos + (alignmentStart - 1);             // to genomic coordinates
                  mismatch.push_back(mis);
                  if (MisStatus == false) MisStatus = true;
                }
                else {                                                       //forbid this pos as a real mis pos
                  unsigned int forbidpos = pos + (alignmentStart - 1);
                  fbpos.push_back(forbidpos);
                  if (MisStatus == false) MisStatus = true;
                }
              } // not in clipped region
            } // yes clip
            
          } //iterator of tagMD
        }   //tagMD > 1
      }     //get MD
      // end: decode the MD tag to get the mismatches

      // skip piling up reads (taking into account: mismatches & clipping information)
      string alignSum = int2str(alignmentStart) + int2str(alignmentEnd) + strand;

      if (alignmentStart != oldstart){
        pileup.clear();                      //clear pileup set
        if (clipStatus == true)
          pileup.insert(map<string, unsigned int>::value_type(alignSum, 0));                                   // 0: a clipped read
        else{
          if (tagMD.size() > 1) pileup.insert(map<string, unsigned int>::value_type(alignSum, 2));             // 2: a quality read with mismatches
          else pileup.insert(map<string, unsigned int>::value_type(alignSum, 1));                              // 1: a quality read with no mismatches
        }
      }
      else if (alignmentStart == oldstart) {
        if (pileup.count(alignSum)) {                            // looks like a pileup
          if ( clipStatus == false ){                            // a perfect read
            if ( tagMD.size() > 1 ) {                            // a perfect read with mismatch
              if (pileup[alignSum] == 2) continue;               // already got, skip
              if (pileup[alignSum] == 1) pileup[alignSum] += 1;  // let this read in
              if (pileup[alignSum] == 0) pileup[alignSum] += 2;  // let this read in
            }
            else {                                               // a perfect read with out mismatch
              if (pileup[alignSum] == 2) continue;
              if (pileup[alignSum] == 1) continue;
              if (pileup[alignSum] == 0) pileup[alignSum] += 1;
            }
          }
          else                                                   // a clipped read
            continue;
        }                                     // found pileup key
        else{                                 // not found key
          if ( clipStatus == true )           // a clipped read
            pileup.insert(map<string, unsigned int>::value_type(alignSum, 0));
          else{                               // a perfect read
            if ( tagMD.size() > 1 ) pileup.insert(map<string, unsigned int>::value_type(alignSum, 2));
            else pileup.insert(map<string, unsigned int>::value_type(alignSum, 1));
          }
        }
      }
      // end: skip piling up reads (taking into account: mismatches & clipping information)


      oldstart = alignmentStart;
      
      // loop for all the regions
      struct SVread tmpread;            // storing current read
      tmpread.name  = bam.Name;
      tmpread.start = alignmentStart;
      tmpread.end   = alignmentEnd;
      tmpread.strand= strand;
      tmpread.seq   = queryB;

      deque <struct region>::iterator iter = regions.begin();

      if ( iter->start > alignmentEnd ) continue;    // skip reads not overlapping with the first region

      while ( iter->chr == old_chr && iter->start <= alignmentEnd && iter != regions.end() ) {

        bool endinside = false;
        bool misinside = false;

        if (iter->end < alignmentStart) {            // the region end is beyond the alignmentStart
          print_mismatch(*iter);            // print out 
          iter = regions.erase(iter);                // this region should be removed 
          if ( regions.empty() ){                    // regions is empty
            getline(region_f, line);                 // get a line of region file
            if ( ! region_f.eof() ){
              eatline(line, regions);                // eat a line and put it into the duque
              iter = regions.begin();
            }
            else finished(2);
          }
          continue;
        }

        if (iter->end >= alignmentStart && iter->start <= alignmentEnd) {  //overlapping, should add some coverage

          iter->coverage++;          
          tmpread.me.clear();

          //add end coverage to each pos, two pairs: alignmentStart-readends1 & readends2-alignmentEnd
          unsigned int region_pos;
          if (alignmentStart <= iter->start) {
            if (readends1 < iter->start && readends2 > iter->start) region_pos = readends2;
            else region_pos = iter->start;
          }
          else  region_pos = alignmentStart;
          for (; region_pos <= iter->end; region_pos++){
            if (region_pos > alignmentEnd) 
              break;
            if ((region_pos >= alignmentStart && region_pos <= readends1) || (region_pos >= readends2 && region_pos <= alignmentEnd)) {
              if (! (iter->posendcov).count(region_pos) ) 
                iter->posendcov.insert( map <unsigned int, unsigned int>::value_type(region_pos, 1) );
              else
                (iter->posendcov)[region_pos] ++;
            }
          } //add end coverage

          if ((alignmentStart >= iter->start && alignmentStart <= iter->end) || (alignmentEnd >= iter->start && alignmentEnd <= iter->end)) {
            endinside   = true;     // current read is ending inside the region
          }          

          if (MisStatus == true) {  // do some mismatch counting

            vector <unsigned int>::iterator misiter = mismatch.begin();
            for(; misiter != mismatch.end(); misiter++) {
              if (*misiter >= iter->start && *misiter <= iter->end) {  //mismatches found

                iter->mismatch++;

                tmpread.me.insert(*misiter);
                misinside = true;        // current read has ME

                if (! (iter->mispos).count(*misiter) )
                  iter->mispos.insert( map <unsigned int, unsigned int>::value_type(*misiter, 1) );
                else
                  (iter->mispos)[*misiter] ++;
              }
            }  //iterator of mismatch(es)

            vector <unsigned int>::iterator fbiter = fbpos.begin();     // forbidden positions
            for(; fbiter != fbpos.end(); fbiter++){
              if (*fbiter >= iter->start && *fbiter <= iter->end) {
                if (! (iter->forbid).count(*fbiter) )
                  iter->forbid.insert(*fbiter);     // insert the forbidden pos 
              }
            }  // iterator of fbpos
          }    // the read has mismatch(es)
        }      // overlapping

        if (endinside == true && misinside == true) {
          (iter->SVreads).push_back(tmpread);         // store current read if end inside and mis inside
        }

        if ( (iter+1) != regions.end() ) 
          iter++;                                     // if this region is not the last element in the deque
        else {                                        // the last element
          getline(region_f, line);                    // get a line of region file
          if ( ! region_f.eof() ){
            eatline(line, regions);                   // eat a line and put it into the duque
            iter = regions.end();
            iter--;
          }
          else  finished(3);
        }
      
      } //while the read is overlapping a region

    } //read a new bam

    
    // bam alignments in this region have been read, need to loop back
    it = regions.begin();                    //reset to beginning
    for (; it != regions.end() && it->chr == old_chr; ) {  // there are some regions left
      print_mismatch(*it);          // print the old region info
      it = regions.erase(it);                // erase the current region
    }

    while ( regions.empty() ) {

      getline(region_f, line);
      if ( region_f.eof() )  finished(4);
      eatline(line, regions);
      it = regions.begin();
      if (it->chr == old_chr){
        print_mismatch(*it);
        regions.clear();
        continue;
      }
    } // region is empty

  } //new chromosome from region file


  //close everything
  regions.clear();
  reader.Close();
  region_f.close();

  return 0;

} //main

inline string int2str(unsigned int &i){
  string s;
  stringstream ss(s);
  ss << i;
  return ss.str();
}

inline string flo2str(float &f){
  string s;
  stringstream ss(s);
  ss << f;
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

inline void splitgfftag(const string &str, map <string, string> &elements, const string &delimiter, vector<string> &tag_want) {

  string::size_type lastPos = str.find_first_not_of(delimiter, 0);
  string::size_type pos     = str.find_first_of(delimiter, lastPos);

  while (string::npos != pos || string::npos != lastPos) {
    string tmp_tag = str.substr(lastPos, pos - lastPos);
    vector <string> tag_now;
    splitstring(tmp_tag, tag_now, "=");

    if( tag_now.size() == 2 ){  //skip strange tags
      vector <string>::iterator tagiter = tag_now.begin();
      string tagtitle = *tagiter;
      string tagthing = *++tagiter;
      vector <string>::iterator wantiter = tag_want.begin();
      while (wantiter != tag_want.end()){
        if (tagtitle == *wantiter) {    
          elements.insert(map <string, string>::value_type(tagtitle,tagthing));
        }
        wantiter++;
      }
    }

    lastPos = str.find_first_not_of(delimiter, pos);
    pos = str.find_first_of(delimiter, lastPos);
  }
}

inline void eatline(const string &str, deque <struct region> &regions_ref) {
  
  vector <string> line_content;
  //split line and then put it into a deque

  splitstring(str, line_content, "\t");
  vector <string>::iterator iter = line_content.begin();
  unsigned int i;
  struct region tmp;

  for(i = 1; iter != line_content.end(); iter++, i++){
    switch (i) {
    case 1: // chr
      tmp.chr = *iter;
      continue;
    case 2: // start
      tmp.start = atoi((*iter).c_str());
      continue;
    case 3: // end
      tmp.end = atoi((*iter).c_str());
      continue;
    case 4: // dis
      tmp.dis = atoi((*iter).c_str());
      continue;
    case 5: // depth
      tmp.depth = atof((*iter).c_str());
      continue;
    case 6: // ratio1
      tmp.ratio1 = atof((*iter).c_str());
      continue;
    case 7: // ratio2
      tmp.ratio2 = atof((*iter).c_str());
      continue;
    case 8: // score
      tmp.score = atof((*iter).c_str());
      continue;
    default:
      break;
    }
  }
  tmp.coverage = 0;
  tmp.mismatch = 0;
  regions_ref.push_back(tmp);

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
    case ('D') : break;
      blockLength  += cigItr->Length;
      currPosition += cigItr->Length;
    case ('P') : break;
    case ('N') :
      blockStarts.push_back(currPosition + cigItr->Length);
      blockLengths.push_back(blockLength);
      currPosition += cigItr->Length;
      blockLength = 0;
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

inline void print_mismatch(struct region &region){  // do some mismatch screening thresholding to reach high accuracy

  unsigned int realmis      = 0;
  unsigned int totalmispos  = 0;
  float        totalendbase = 0;  

  map <unsigned int, unsigned int>::iterator posenditer = region.posendcov.begin();
  for (; posenditer != region.posendcov.end(); posenditer++){
    totalendbase += posenditer->second;
  }
  float leasterr = 0.01;
  float n_emis = region.mismatch;
  float localerr = n_emis/totalendbase;
  float baserr = max(leasterr, localerr);
  double mismatch_score = 0.;

  map <unsigned int, unsigned int>::iterator iter = region.mispos.begin();  
  for(;iter != region.mispos.end(); iter++){
    if (! region.forbid.count(iter->first) ) { // not found in forbidden pos
      totalmispos++; 
      mismatch_score += pbinom(iter->second, (region.posendcov)[iter->first], baserr, 0);
      if ( iter->second >= 2 ) {
        realmis++;
      } // frequency > 2
    }   // forbidden pos
  }     // for each mis pos

  unsigned int max_nme = 0;
  map <unsigned int, vector <struct SVread> > topreads;
  vector <struct SVread>::iterator riter = region.SVreads.begin();
  for (; riter != region.SVreads.end(); riter++){
    set <unsigned int>::iterator miter = riter->me.begin();
    unsigned int nme = 0;
    while ( miter != riter->me.end() ) {
      if ( ! region.forbid.count(*miter) ) {
        nme ++;
        miter++;
      }
      else {
        riter->me.erase(miter++);  // clean non-end mis
      }
    }
    topreads[nme].push_back(*riter);    
    if ( nme > max_nme ) max_nme = nme;
  } // current SVread

  string seedseq = "RME";    // decide the seed sequence 
  if (max_nme > 0) {
    if (topreads[max_nme].size() == 1){
      vector <struct SVread>::iterator topiter = topreads[max_nme].begin();
      set <unsigned int>::iterator meiter = (topiter->me).begin();
      if ( (*meiter - topiter->start) < (topiter->end - *meiter) ) {
        if (topiter->seq == "NA") seedseq = (topiter->name)+"["+(topiter->strand)+"p]";
        else seedseq = (topiter->seq).substr(0,25);
      }
      else {
        if (topiter->seq == "NA") seedseq = (topiter->name)+"["+(topiter->strand)+"s]";
        else seedseq = (topiter->seq).substr((topiter->seq).length()-25,25);
      }
    }  
    else {
      vector <struct SVread>::iterator topiter  = topreads[max_nme].begin();
      vector <struct SVread>::iterator topiter2 = topreads[max_nme].end();
      topiter2--;
      if (topiter->start >= region.start && topiter->end > region.end) { // the first top read is starting in the region
         if (topiter->seq == "NA") seedseq = (topiter->name)+"["+(topiter->strand)+"p]";
         else seedseq = (topiter->seq).substr(0,25);       
      }
      else if (topiter2->start < region.start && topiter2->end <= region.end) {
         if (topiter->seq == "NA") seedseq = (topiter->name)+"["+(topiter->strand)+"s]";
         else seedseq = (topiter2->seq).substr((topiter2->seq).length()-25,25);
      }
      else {
        // where is the changing point? if ratio2 > 0.5 take the start, else take the end
        vector <struct SVread>::iterator toprem = topiter;
        for (; topiter != topreads[max_nme].end(); topiter++) {

          set <unsigned int>::iterator meiter  = (topiter->me).begin(); // the first mis
          set <unsigned int>::iterator meiter2 = (topiter->me).end();   // the last  mis
          meiter2--;

          if ( topiter->end <= region.end && (*meiter - topiter->start) > (topiter->end - *meiter) ) {
            toprem = topiter;
          }
          if ( topiter->start >= region.start && (*meiter2 - topiter->start) < (topiter->end - *meiter2)) {  // now its turning
              if ( region.ratio2 > 0.5 ) {
                if (topiter->seq == "NA") seedseq = (topiter->name)+"["+(topiter->strand)+"p]";
                else seedseq = (topiter->seq).substr(0,25);
              }
              else {
                if (topiter->seq == "NA") seedseq = (topiter->name)+"["+(topiter->strand)+"s]";
                else seedseq = (toprem->seq).substr((toprem->seq).length()-25,25);
              }
              break;
          } // turning
        } // for top iter

        if (seedseq == "RME"){
          if (toprem->seq == "NA") seedseq = (toprem->name)+"["+(toprem->strand)+"s]";
          else seedseq = (toprem->seq).substr((toprem->seq).length()-25,25);
        }

      } //else
    } // > 1
  } //max_nme > 0

  float fdepth = region.coverage;
  float edepth = fdepth*(region.ratio1);
  float mirate = (region.mismatch)/edepth;

  string chrom        = region.chr.c_str();
  if (chrom.substr(0,3) != "chr") chrom = "chr"+chrom;
  string source       = "Breakpointer";
  string type         = "Depth-Skewed";
  unsigned int start  = region.start;
  unsigned int end    = region.end;
  float confi         = mismatch_score;
  string strand       = "+";
  string phase        = ".";
  string tag          = "ID="+chrom+":"+int2str(start)+";SIZE="+int2str(region.dis)+";DEPTH="+int2str(region.coverage)+";EndsRatio="+flo2str(region.ratio1)+";StartsRatio="+flo2str(region.ratio2)+";BinomialScore="+flo2str(region.score)+";MIS="+int2str(region.mismatch)+";realMIS="+int2str(realmis)+";MISRATE="+flo2str(mirate)+";seedseq="+seedseq;

  if ( totalmispos > 1 )  //screen for number of positions where there are mismatches
    if ( region.coverage >= 5 )  //screen for coverage
     //if ( !(region.coverage > 100 && region.score < 1.1) ) // filter out hard to say stuff just for XLMR project!!!!!
      //if ( realmis > 0 || (mirate > 1 && region.coverage >= 10 && region.mismatch < 50) ) //screen for realmis and mirate
        cout << chrom <<"\t"<< source <<"\t"<< type <<"\t"<< start <<"\t"<< end <<"\t"<< setprecision(3) << confi <<"\t"<< strand <<"\t"<< phase <<"\t"<< tag << endl;
}


inline void finished(const unsigned int &where){
  cerr << "Finished: end of region file, Zone: " << where << endl;
  exit(0);
}
