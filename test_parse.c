#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include "mol_general.h"
#include "GeneralHashFunctions.h"

bool opt_bitstring;
bool opt_none; 
bool opt_verbose;
bool opt_text;
bool opt_code;
bool opt_bin;
bool opt_stdin;
bool opt_exact;
bool opt_debug;
bool opt_molout;
bool opt_molstat;
bool opt_molstat_X;
bool opt_molstat_v;  // new in v0.4d
bool opt_xmdlout;
bool opt_strict;  // new in v0.2f
bool opt_nochirality; // new in v0.4e (thanks to Robert Bruccoleri for this suggestion)
bool opt_metalrings;  // new in v0.3
bool opt_geom;  // new in v0.3d
bool opt_chiral;  // new in v0.3f
bool opt_iso;  // new in v0.3p
bool opt_chg;  // new in v0.3p 
bool opt_rad;  // new in v0.3p
int  opt_rs;  // new in v0.3i
bool opt_fp;  // new in v0.3m
int  fpformat;  // new in v0.3m
bool opt_hfp;  // new in v0.4
bool opt_matchnum;  // new in v0.4a
bool opt_matchnum1;  // new in v0.4a
bool opt_pos;  // v0.5
bool opt_list;

int  hfpformat;
char* molfilename;
char* ndl_molfilename;
char* molname;
char* molcomment;
int  n_atoms;
int  n_bonds;
int  n_rings;    // the number of rings we determined ourselves
int  n_countablerings; // v0.3n; number of rings which are not envelope rings
int  n_cmrings;  // the number of rings we read from a (CheckMol-tweaked) MDL molfile
int  n_charges;  // number of charges    // unnecessary!!!!
int  n_heavyatoms;
int  n_trueheavyatoms;  // v0.4b
int  n_heavybonds;
int  ndl_n_atoms;
int  ndl_n_bonds;
int  ndl_n_rings;
int  ndl_n_heavyatoms;
int  ndl_n_trueheavyatoms;  // v0.4b
int  ndl_n_heavybonds;
  //cm_mdlmolfile  : boolean;
bool  found_arominfo;
bool  found_querymol;
bool  ndl_querymol;  // v0.3p
int  tmp_n_atoms;      // v0.3m
int  tmp_n_bonds;      // v0.3m
int  tmp_n_rings;      // v0.3m
int  tmp_n_heavyatoms; // v0.3m
int  tmp_n_trueheavyatoms;  // v0.4b
int  tmp_n_heavybonds; // v0.3m

typedef struct ATOM_REC {
   char element[3];
   char atype[4];
   float x;
   float y;
   float z;
   int formal_charge;
   float real_charge;
   short Hexp;
   short Htot;
   int neighbor_count;
   int ring_count;
   bool arom;
   bool q_arom;
   bool stereo_care;
   bool heavy;
   bool metal;
   int nvalences;
   bool tag;
   int nucleon_number;
   int radical_type;
}ATOM_REC;

typedef struct BOND_REC {
  int a1;
  int a2;
  char btype;
  int ring_count;
  bool arom;
  bool q_arom;
  short topo;
  short stereo;
  short mdl_stereo;
}BOND_REC;

typedef struct MOLSTAT_REC{
              int n_QA, n_QB, n_chg;                  // number of query atoms, query bonds, charges
              int n_C1, n_C2, n_C;                    // number of sp, sp2 hybridized, and total no. of carbons
              int n_CHB1p, n_CHB2p, n_CHB3p, n_CHB4;  // number of C atoms with at least 1, 2, 3 hetero bonds
              int n_O2, n_O3;                         // number of sp2 and sp3 oxygens
              int n_N1, n_N2, n_N3;                   // number of sp, sp2, and sp3 nitrogens
              int n_S, n_SeTe;                        // number of sulfur atoms and selenium or tellurium atoms
              int n_F, n_Cl, n_Br, n_I;               // number of fluorine, chlorine, bromine, iodine atoms
              int n_P, n_B;                           // number of phosphorus and boron atoms
              int n_Met, n_X;                         // number of metal and "other" atoms (not listed elsewhere); v0.3l
              int n_b1, n_b2, n_b3, n_bar;            // number single, double, triple, and aromatic bonds
              int n_C1O, n_C2O, n_CN, n_XY;           // number of C-O single bonds, C=O double bonds, CN bonds (any type), hetero/hetero bonds
              int n_r3, n_r4, n_r5, n_r6, n_r7, n_r8; // number of 3-, 4-, 5-, 6-, 7-, and 8-membered rings
              int n_r9, n_r10, n_r11, n_r12, n_r13p;  // number of 9-, 10-, 11-, 12-, and 13plus-membered rings
              int n_rN, n_rN1, n_rN2, n_rN3p;         // number of rings containing N (any number), 1 N, 2 N, and 3 N or more
              int n_rO, n_rO1, n_rO2p;                // number of rings containing O (any number), 1 O, and 2 O or more
              int n_rS, n_rX, n_rar, n_rBz;          // number of rings containing S (any number), any heteroatom (any number), 
                                                                // number of aromatic rings, number of benzene rings
              int n_br2p;                             // number of bonds belonging to more than one ring (v0.3n)
//                  {$IFDEF extended_molstat}
              int n_psg01, n_psg02, n_psg13, n_psg14; // number of atoms belonging to elements of group 1, 2, etc. 
              int n_psg15, n_psg16, n_psg17, n_psg18; // number of atoms belonging to elements of group 15, 16, etc. 
              int n_pstm, n_psla;                     // number of transition metals, lanthanides/actinides
//                  {$ENDIF}
}MOLSTAT_REC;

typedef struct RINGPROP_REC { // new in v0.3
                int    size     ;
                bool    arom    ;
                bool    envelope;
} RINGPROP_REC;

typedef struct CONVAL_REC { // new in v0.3j
int def;         // better as longint for large molecules?
int tmp;
}CONVAL_REC;

  // new in v0.3q: periodic table lookup
typedef struct PT_REC{ 
char el[3];
float am;    // average isotopic mass; v0.4d  
int count;   // v0.4d
}PT_REC;

typedef struct FRAG_REC{
int size;
int a_atomicnum[max_ringfragpath_length];
char b_code[max_ringfragpath_length];
bool ring;
}FRAG_REC;

typedef int NEIGHBOR_REC[max_neighbors];
typedef int RINGPATH_TYPE[max_ringsize];
typedef char FRAGSTR[max_fragstr_length];
//typedef RINGPATH_TYPE RINGLIST[max_rings]; //i

ATOM_REC *atom[max_atoms];  
BOND_REC *bond[max_bonds];
RINGPROP_REC *ringprop[max_rings];
CONVAL_REC *cv[max_atoms];  // new in v0.3j
int *ring[max_rings];
int *fgloc[used_fg+1]; // fgloc[0] is empty.
bool fg[max_fg];
bool *gmm[max_ndl_gmmsize];

MOLSTAT_REC molstat;
PT_REC pt[max_atomicnum+1];

char molbuf[max_atoms+max_bonds+8192][BUFSIZ];
int molbufindex;
int li;  // the line number of a molecule.
bool mol_in_queue;
int mol_count;
int ringsearch_mode;

bool  mol_OK;  // new in v0.2i
int n_ar     ;  // new in v0.3
int  prev_n_ar ;  // new in v0.3
bool  ez_search ;  // new in v0.3d
bool  rs_search ;  // new in v0.3f
bool  ez_flag ;  // new in v0.3f
bool  chir_flag ;  // new in v0.3f
bool  rs_strict ;  // new in v0.3j

int  n_Ctot;
int n_Otot;
int n_Ntot;
int ndl_n_Ctot;
int ndl_n_Otot;
int ndl_n_Ntot;  // new in v0.3g
int tmp_n_Ctot;
int tmp_n_Otot;
int tmp_n_Ntot;  // new in v0.3m

bool  ether_generic;   // v0.3j
bool  amine_generic;   // v0.3j  
bool  hydroxy_generic; // v0.3j
long  fpdecimal     ;  // v0.3m
long  fpincrement   ;
int  fpindex        ;
bool  fp_exacthit;
bool  tmfmismatch   ;  // v0.3m; rely on tweaked MDL molfiles only if same level
bool  auto_ssr      ;  // v0.3n; indicates that SAR -> SSR fallback has happened
long  recursion_level; // v0.3p
long  recursion_depth; // v0.3p
bool  keep_DT        ; // v0.3p  molfile output: Deuterium/Tritium as D/T or 2H/3H
bool  rfile_is_open ;
bool  fp_exactblock ;
int  tmfcode ;  // v0.3m; version number for tweaked MDL molfiles (tweaklevel)
  // new in v0.4: periodic table lookup, fragments for hashed fingerprints
int max_vringsize;
bool  hfp[hfpsize];
bool fg[max_fg];
int  hfpformat;  // end of new v0.4 variables
  // new in v0.4a: global match matrix (+ cumulative copy), various auxiliary variables
bool  use_gmm ;
bool  valid_gmm;
/*  overall_match  : boolean;
  tmp_tag        : array[1..max_atoms] of boolean;
  found_untagged : boolean;  // end of new v0.4a variables
  n_matches      : integer;  // v0.4b
  match_string   : ansistring;  // v0.4b
  sep_label      : string;  // v0.4b
  mol_formula    : string;  // v0.4d
  mol_weight     : single;  // v0.4d */
int  fglang;  // v0.5

int ATOI(char *str) { //return -1 if the atoi failed, else return the int value;
static int val;
str=strcat(str,"0");
val=atoi(str);
if (val>0) { return val/10; } else {return -1;} 
} 

void init_globals() {
static int i;
opt_verbose     = false;
opt_debug    = false;
opt_exact       = false;
opt_stdin       = false;
opt_text        = false;
opt_code        = false;
opt_bin         = false;
opt_bitstring   = false;
opt_molout      = false;
opt_molstat     = false;
opt_molstat_X   = false;
opt_molstat_v   = false;  // new in v0.4d
opt_xmdlout     = false;
opt_strict      = false;  // new in v0.2f
opt_metalrings  = false;  // new in v0.3
opt_geom        = false;  // new in v0.3d
opt_chiral      = false;  // new in v0.3f
opt_fp          = false;  // new in v0.3m
opt_iso         = false;  // new in v0.3p
opt_chg         = false;  // new in v0.3p
opt_rad         = false;  // new in v0.3p
  //cm_mdlmolfile   := false;
found_arominfo  = false;
found_querymol  = false;
ndl_querymol    = false;
opt_rs          = rs_sar;  // v0.3i
ringsearch_mode = rs_sar;
rfile_is_open   = false;  // new in v0.2g
ez_search       = false;  // new in v0.3d
rs_search       = false;  // new in v0.3f
ez_flag         = false;  // new in v0.3f
chir_flag       = false;  // new in v0.3f
rs_strict       = false;  // new in v0.3j
n_Ctot = 0;
n_Otot = 0;
n_Ntot = 0;             // new in v0.3g
ndl_n_Ctot = 0;
ndl_n_Otot = 0;
ndl_n_Ntot = 0; // new in v0.3g
  for(i=1;i<=max_fg;i++) {fg[i] = false;}
/*  try
    molbuf[0]=(MOLSTAT_REC *)malloc(sizeof(MOLSTAT_REC));
    getmem(molbuf,sizeof(molbuftype));
  except
    on e:Eoutofmemory do
      begin
        writeln('Not enough memory');
        halt(4);
      end;
  end; */
ether_generic   = false;       // v0.3j
amine_generic   = false;       // v0.3j
hydroxy_generic = false;       // v0.3j
fpformat        = fpf_decimal; // v0.3m
fpindex         = 0;           // v0.3m
fp_exacthit     = false;       // v0.3m
fp_exactblock   = false;       // v0.3m
tmfcode         = 0;           // v0.3m
tmfmismatch     = false;       // v0.3m
auto_ssr        = false;       // v0.3n
keep_DT         = true;        // v0.3p
opt_hfp         = false;       // v0.4
hfpformat       = fpf_decimal; // v0.4
opt_matchnum    = false;       // v0.4a
opt_matchnum1   = false;       // v0.4a
opt_pos         = false;       // v0.5
  use_gmm         = false;        // v0.4a
  valid_gmm       = false;        // v0.4a
  fglang          = -1;           // v0.5
}

void readbuffer(char *rfile) {
int i,pos_0,pos_1; 
static char line_buffer[BUFSIZ];
static int line_number;
line_number=0;
pos_0=0; 
for(i=0;i<strlen(rfile);i++) {
  if(rfile[i]=='\n'){ 
    line_number++;
    pos_1=i;
//    printf("%c, %d\n",rfile[i+1], i);
    memset(line_buffer,0,sizeof(line_buffer));
    strncpy(line_buffer,rfile+pos_0+1,pos_1-pos_0);
    line_buffer[strlen(line_buffer)-1]='\0';
    strcpy(molbuf[line_number],line_buffer);
    pos_0=pos_1;
  }
}
//last line:
    memset(line_buffer,0,sizeof(line_buffer));
    strncpy(line_buffer,rfile+pos_0,strlen(rfile)-pos_0);
    line_buffer[strlen(line_buffer)-1]='\0';
    strcpy(molbuf[line_number+1],line_buffer);
molbufindex=line_number+1;
//    printf("%s: %d\n",molbuf[line_number+1],molbufindex);
}

void readfile() { 
static char *eom;
int i;
static FILE *infile; 
static char line_buffer[150]; /* BUFSIZ is defined if you include stdio.h */ 
static int line_number; 
if (!opt_stdin) {
  infile = fopen(molfilename, "r"); 
  if (!infile) { 
    printf("Couldn't open file %s for reading.\n", molfilename); 
    exit(0);
  } 
  line_number = 0; 
  while (fgets(line_buffer, sizeof(line_buffer), infile)) { 
  if (line_number < (max_atoms+max_bonds+64)) {        
    ++line_number; /* note that the newline is in the buffer */ 
    for (i=strlen(line_buffer)-1;i>=0;i--){
      if(line_buffer[i]=='\n'||line_buffer[i]=='\r'){line_buffer[i]='\0'; } else {break;} }
    strcpy(molbuf[line_number],line_buffer);
//   printf("%s: %d\n",molbuf[line_number],line_number);
  } else {
    printf("Not enough memories for molfile!");
    exit(0);
  }
//printf("%4d: %s", line_number, molbuf[line_number]); 
  eom=strstr(line_buffer,molend);
  if(eom){molbufindex=line_number;}
  } fclose(infile);
//printf("\nTotal number of lines = %d\n", molbufindex); 
} else {  line_number = 0; 
  while (fgets(line_buffer,sizeof(line_buffer), stdin)!=NULL) {
  if (line_number < (max_atoms+max_bonds+64)) {        
    ++line_number; /* note that the newline is in the buffer */ 
  eom=strstr(line_buffer,molend);
  if(eom){molbufindex=line_number;}
    for (i=strlen(line_buffer)-1;i>=0;i--){
      if(line_buffer[i]=='\n'||line_buffer[i]=='\r'){line_buffer[i]='\0'; } else {break;} }
//    line_buffer[strlen(line_buffer)]='\0';
    strcpy(molbuf[line_number],line_buffer);
//    printf("%s: %d\n",molbuf[line_number],line_number);
  } else {
    printf("Not enough memories for molfile!");
    exit(0);
  }
//printf("%4d: %s", line_number, molbuf[line_number]); 
  } 
}
} 

void zap_molecule() {
    if (atom != NULL) {
        free(*atom);
        *atom=NULL;
    }
    if (bond != NULL) {
        free(*bond);
        *bond=NULL;
    }
    if (ring != NULL) {
        free(*ring);
        *ring=NULL;
    }
    if (ringprop != NULL) {
        free(*ringprop);
        *ringprop=NULL;
    }
    if (fgloc!= NULL) {
        free(*fgloc);
        *fgloc=NULL;
    }
  n_atoms = 0;
  n_bonds = 0;
  n_rings = 0;
}

bool is_heavyatom(int id) {
  bool r;
  char el[2];
  r=true;
  strcpy(el,atom[id]->element);
  if (!strcmp(el,"DU")||!strcmp(el,"LP")) {r = false;}
  if (!strcmp(el,"H ")) { 
      if (opt_iso == false) {r = false; } else
        {
          if (atom[id]->nucleon_number < 2) { r= false;}
        }
  }
  return r;  
}
  // note: deuterium is regarded as a heavy atom if -i option is used;
  // this applies only to matchmol, in checkmol D/T is always non-heavy
bool is_trueheavyatom(int id) {   // v0.4b
  bool  r;
  char el[2];
  r  = true;
  strcpy(el,atom[id]->element);
  if (!strcmp(el,"DU") || !strcmp(el,"LP")) { r = false; }
  if (!strcmp(el,"H ") || !strcmp(el,"D ") || !strcmp(el,"T ")) {r = false;}
  return r;  
}

bool is_metal(int id){
bool  r;
char  el[3];
  r = false;
  strcpy(el,atom[id]->element);
  if (!strcmp(el,"LI") || !strcmp(el,"NA") || !strcmp(el,"K ") || !strcmp(el,"RB") ||
     !strcmp(el,"CS") || !strcmp(el,"BE") || !strcmp(el,"MG") || !strcmp(el,"CA") || 
     !strcmp(el,"SR") || !strcmp(el,"BA") || !strcmp(el,"TI") || !strcmp(el,"ZR") || 
     !strcmp(el,"CR") || !strcmp(el,"MO") || !strcmp(el,"MN") || !strcmp(el,"FE") ||
     !strcmp(el,"CO") || !strcmp(el,"NI") || !strcmp(el,"PD") || !strcmp(el,"PT") ||
     !strcmp(el,"SN") || !strcmp(el,"CU") || !strcmp(el,"AG") || !strcmp(el,"AU") ||
     !strcmp(el,"ZN") || !strcmp(el,"CD") || !strcmp(el,"HG") || !strcmp(el,"AL") ||
     !strcmp(el,"SN") || !strcmp(el,"PB") || !strcmp(el,"SB") || !strcmp(el,"BI")) 
  { r = true; }                                // etc. etc.
return r;
}

int get_nvalences(char a_el[3]) {  // changed name and position in v0.3m
// preliminary version; should be extended to element/atomtype
int  res;
  res = 1;
  if (!strcmp(a_el,"H ")) {res = 1;}
  if (!strcmp(a_el,"D ")) {res = 1;} // v0.3n
  if (!strcmp(a_el,"C ")) { res = 4;}
  if (!strcmp(a_el,"N ")) { res = 3;}
  if (!strcmp(a_el,"O ")) { res = 2;}
  if (!strcmp(a_el,"S ")) { res = 2;}
  if (!strcmp(a_el,"SE")) { res = 2;}
  if (!strcmp(a_el,"TE")) { res = 2;}
  if (!strcmp(a_el,"P ")) { res = 3;}
  if (!strcmp(a_el,"F ")) { res = 1;}
  if (!strcmp(a_el,"CL")) { res = 1;}
  if (!strcmp(a_el,"BR")) { res = 1;}
  if (!strcmp(a_el,"I ")) { res = 1;}
  if (!strcmp(a_el,"B ")) { res = 3;}
  if (!strcmp(a_el,"LI")) { res = 1;}
  if (!strcmp(a_el,"NA")) { res = 1;}
  if (!strcmp(a_el,"K ")) { res = 1;}
  if (!strcmp(a_el,"CA")) { res = 2;}
  if (!strcmp(a_el,"SR")) { res = 2;}
  if (!strcmp(a_el,"MG")) { res = 2;}
  if (!strcmp(a_el,"FE")) { res = 3;}
  if (!strcmp(a_el,"MN")) { res = 2;}
  if (!strcmp(a_el,"HG")) { res = 2;}
  if (!strcmp(a_el,"SI")) { res = 4;}
  if (!strcmp(a_el,"SN")) { res = 4;}
  if (!strcmp(a_el,"ZN")) { res = 2;}
  if (!strcmp(a_el,"CU")) { res = 2;}
  if (!strcmp(a_el,"A ")) { res = 4;}
  if (!strcmp(a_el,"Q ")) { res = 4;}
  return res;
}

void count_neighbors() {// counts heavy-atom neighbors and explicit hydrogens
int i;
  if ((n_atoms < 1) || (n_bonds < 1)) {return;}
  for (i = 0; i< n_bonds; i++) {
      if (atom[bond[i]->a1]->heavy) {atom[(bond[i]->a2)]->neighbor_count++;}
      if (atom[bond[i]->a2]->heavy) {atom[(bond[i]->a1)]->neighbor_count++;}
      if (!strcmp(atom[(bond[i]->a1)]->element,"H ")) {atom[(bond[i]->a2)]->Hexp++;}
      if (!strcmp(atom[(bond[i]->a2)]->element,"H ")) {atom[(bond[i]->a1)]->Hexp++;}
      // plausibility check (new in v02.i)
      if ((atom[(bond[i]->a1)]->neighbor_count > max_neighbors) || 
         (atom[(bond[i]->a2)]->neighbor_count > max_neighbors)) {
           mol_OK = false;
      }
  }
}

bool odd(unsigned num)
{
 return num%2==1;
}

bool rc_identical(int rc_int) {
  if (rc_int == 0) { return true;} else { return false;}
}


bool rc_1in2(int rc_int) {
  if (odd(rc_int)) {return true;} else {return false;}
}

bool rc_2in1(int rc_int) {
  rc_int = rc_int/2;
  if (odd(rc_int)) {return true;} else {return false;}
}

bool rc_different(int rc_int){
  rc_int = rc_int/4;
  if (odd(rc_int)){return true;} else {return false;}
}

bool rc_independent(int rc_int) {
  rc_int = rc_int/8;
  if (odd(rc_int)) {return true;} else {return false;}
}

int path_pos(int id, RINGPATH_TYPE a_path) {  // new version in v0.3l 
int  i, pp;
  pp = 0;
  for (i = 0;i<max_ringsize;i++) {
      if (a_path[i] == id) {
          pp = i;
          break;
      }
  } 
          return pp; //get the i where a_path[i]==id
}

int ringcompare(RINGPATH_TYPE rp1, RINGPATH_TYPE rp2){
int i, j, rc, rs1, rs2;
int  n_common, max_cra;
  rc = 0;
  n_common = 0;
  rs1 = path_pos(0,rp1);
  rs2 = path_pos(0,rp2);
  if (rs1 < rs2) {max_cra = rs1;} else {max_cra = rs2;}
  for (i=0;i<rs1;i++) {
    for (j =0;j<rs2;j++){
      if (rp1[i] == rp2[j]){n_common++;}
    }
  }
  if ((rs1 == rs2) && (n_common == max_cra)) {rc = 0;} else {
    if (n_common == 0) {rc=rc+8;}
    if (n_common < max_cra) {rc=rc+4;} else {
          if (rs1 < rs2) {rc++;} else {rc=rc+2;}
    }
  }
  return rc;
}

bool is_newring(RINGPATH_TYPE n_path) {
int  i, j;
bool  nr, same_ring;
RINGPATH_TYPE  tmp_path;
int  rc_result;
bool  found_ring;
int  pl;  // new in v0.3
  nr = true;
  pl = path_pos(0,n_path);  // new in v0.3
  if (n_rings > 0) {
    if(ringsearch_mode==rs_sar) {
      found_ring = false;
      i=0;
      while ((i < n_rings) && (!found_ring)) {
        if (pl ==ringprop[i]->size){  // compare only rings of same size
          same_ring = true;   
          for (j = 0;j<pl;j++){
            if (*(*(ring+i)+j)!=n_path[j]) {same_ring = false;}
          }
          if (same_ring) {
            nr= false;
            found_ring = true;
            return nr;
          }
        }
        i++;
      }  // while
    } else {
      for (i=0;i<n_rings;i++) {
        for (j=0;j<max_ringsize;j++){tmp_path[j] = *(*(ring+i)+j);}
        rc_result = ringcompare(n_path,tmp_path);
                        if (rc_identical(rc_result)) {nr = false;}
                        if (rc_1in2(rc_result)) {
                            // exchange existing ring by smaller one
                            for (j=0; j<max_ringsize;j++) {*(*(ring+i)+j) = n_path[j];}
                            // update ring property record (new in v0.3)
                            ringprop[i]->size = pl;
                            nr = false;
/*                            {$IFDEF debug}
                            debugoutput('replacing ring '+inttostr(i)+' by smaller one (ringsize: '+inttostr(path_length(n_path))+')');
                            {$ENDIF}*/
                        }
                        if (rc_2in1(rc_result)){ // new ring contains existing one, but is larger ==> discard
                            nr= false;
                        }
        }
    }
  }
  return nr;
}

void add_ring(RINGPATH_TYPE n_path) {
// new in v0.3: store rings in an ordered way (with ascending ring size)
int i, j, k, s, pl ;
  pl = path_pos(0,n_path);
  if (pl < 1) {return;}
  if (n_rings < max_rings) {
      n_rings++;
      j = 0;
      if (n_rings > 1) {
          for (i=1; i<n_rings;i++) {
              s= ringprop[i-1]->size;
              if (pl >= s) {j = i;} //The last added n_path.
          }
      }
                // the next position is ours
      if (j < n_rings) {
// push up the remaining rings by one position
          for (k=(n_rings-1); k>j; k--) {
              ringprop[k]->size = ringprop[k-1]->size;
              //ringprop^[k].arom := ringprop^[(k-1)].arom;
              for (i = 0; i<max_ringsize; i++) {
                  *(*(ring+k)+i) = *(*(ring+k-1)+i);
              }
           }
       }
      ringprop[j]->size = pl;  
      for (i=0;i<path_pos(0,n_path);i++) {
          *(*(ring+j)+i) = n_path[i];
          //inc(atom^[(n_path[i])].ring_count);
       }
  } 
}

bool is_ringpath(RINGPATH_TYPE s_path){
int i, j, nb_count,a_ref, a_tmp, a_left,a_right;
NEIGHBOR_REC  nb;
bool  rp, new_atom;
int  a_last, pl;
RINGPATH_TYPE l_path;
  rp = false;
  new_atom = false;
  pl = path_pos(0,s_path);
  memset(l_path,0,max_ringsize*sizeof(int));
  if (pl>0) {
    for(i=0;i<pl;i++) {l_path[i]=s_path[i];}
  // check if the last atom is a metal and stop if opt_metalrings is not set (v0.3)
    if (opt_metalrings == false) {
      if (atom[l_path[pl-1]]->metal) {rp=false;
      return false; 
      }
    }
  // check if ring is already closed
    if ((pl > 2) && (l_path[pl-1] == l_path[0])) { 
  //order ring path
     if (pl < 3) {return rp;}
     a_ref = n_atoms;  // start with highest possible value for an atom number
     for (i =0;i< pl;i++) {
        if (l_path[i] < a_ref) {a_ref = l_path[i];}  // find the minimum value ==> reference atom
      }
      if (path_pos(a_ref,l_path) < (pl-1)) {a_right = l_path[(path_pos(a_ref,l_path)+1)];} else {a_right = l_path[0];}//1,9,3->left 3, right 9
      if (path_pos(a_ref,l_path) > 0) {a_left = l_path[(path_pos(a_ref,l_path))-1];} else {a_left = l_path[pl-2];} 

      if (a_right == a_left) {return false;}  // should never happen
      if (a_right < a_left) {// correct ring numbering direction, only shift of the reference atom to the left end required
        while (path_pos(a_ref,l_path) > 0) {
            a_tmp = l_path[1];
            for (i =0;i< (pl -1);i++) {l_path[i] = l_path[(i+1)];}
            l_path[pl-1] = a_tmp;
        }
      } 
      if (a_right > a_left) {  // wrong ring numbering direction, two steps required
        while (l_path[pl-1]!=a_ref) {  // step one: create "mirrored" ring path with reference atom at right end
            a_tmp = l_path[pl-2];
            for (i=(pl-1); i>0; i--) {l_path[i] = l_path[(i-1)];}
            l_path[0] = a_tmp;
        }  
        for (i=0;i<(pl/2);i++) {// one more mirroring
            a_tmp = l_path[i];
            l_path[i] = l_path[(pl-1)-i];
            l_path[(pl-1)-i] = a_tmp;
        }//9,1,3
      }  
      l_path[pl-1]=0; // remove last entry (redundant!)}
      if (is_newring(l_path)) { //printf("test search_path fourth2: %d, %d, %d,%d, %d, %d, %d, %d\n",a_ref, l_path[0],l_path[1],l_path[2], l_path[3], l_path[4], l_path[5], l_path[6]);
            if (n_rings < max_rings) {add_ring(l_path);
            } else{ rp=false;
            }
        rp=true;
      } 
      return rp;
   }
  // any other case: ring is not (yet) closed; 
  // if the ring path exists, increase one atom to the search path. 
    a_last = l_path[pl-1];
    memset(nb,0,sizeof(NEIGHBOR_REC));
    nb_count=0;
    for (i= 0;i< n_bonds;i++ ){
     if ((bond[i]->a1 == a_last) && (nb_count < max_neighbors) && (atom[bond[i]->a2]->heavy)) {
        nb[nb_count]=bond[i]->a2;
        nb_count++;
      }
      if ((bond[i]->a2== a_last)&& (nb_count < max_neighbors)&& (atom[bond[i]->a1]->heavy)) {
        nb[nb_count]=bond[i]->a1;
        nb_count++;
      }
    }
    if (atom[a_last]->neighbor_count > 1) {
      if ((rp == false) && (n_rings < max_rings)) { // added in v0.2: check if max_rings is reached
 // if ring is not closed, continue searching
          for (i= 0; i<atom[a_last]->neighbor_count;i++) {
              new_atom = true;//printf("%d,%d,%d,%d,%d,%d\n",a_last,atom[a_last]->neighbor_count,nb[0],nb[1],nb[2],n_rings);
              for (j =1; j<pl; j++) {
                if (nb[i]==l_path[j]) {      // v0.3k
                  new_atom = false; // v0.3k
                  break;
                }
              }
              // added in v0.1a: check if max_rings not yet reached
              // added in v0.2:  limit ring size to max_vringsize instead of max_ringsize
              if ((new_atom) && (pl < max_vringsize) && (n_rings < max_rings)) {
                  l_path[pl] = nb[i];
                  if (pl < max_ringsize-1) {l_path[pl+1] = 0; } // just to be sure
                  recursion_level++;                             // v0.3p (begin)
                  if (recursion_level > max_recursion_depth) {
                      n_rings = max_rings;
                      return false;
                  }       
                                          // v0.3p (end)
                  if (is_ringpath(l_path)) {rp=true;}
              }
        }
      }
    }
  }
  return rp;
}

void chk_ringbonds(){ //check rings 
int i, j,n, a1rc, a2rc;
int ra1, ra2;
int nb_count;
NEIGHBOR_REC nb;
RINGPATH_TYPE search_path; 
  memset(search_path,0,max_ringsize*sizeof(int));
if (n_bonds>0) {
  for(i=0; i<n_bonds; i++){
    a1rc=atom[bond[i]->a1]->ring_count;
    a2rc=atom[bond[i]->a2]->ring_count;
    if ((n_rings==0) || (a1rc<n_rings) ||(a2rc<n_rings)) {
      recursion_level=0;
      ra1=bond[i]->a1;
      ra2=bond[i]->a2;
      nb_count = 0;
      for (j= 0;j< n_bonds;j++ ){ //get ra2's neighbors!
        if ((bond[j]->a1 == ra2) && (nb_count < max_neighbors) && (atom[bond[j]->a2]->heavy)) {
          nb[nb_count]=bond[j]->a2;
          nb_count++;
        }
        if ((bond[j]->a2== ra2)&& (nb_count < max_neighbors)&& (atom[bond[j]->a1]->heavy)) {
          nb[nb_count]=bond[j]->a1;
          nb_count++;
        }
      }
      if((atom[ra2]->neighbor_count > 1)&&(atom[ra1]->neighbor_count>1)) {
        search_path[0]=ra1;
        search_path[1]=ra2;
//        printf("gethere:%d,%d,%d,%d\n",n_rings, ra1, ra2, nb[0]);
        for(n=0; n< atom[ra2]->neighbor_count; n++){
          if ((nb[n]!=ra1)&&(atom[nb[n]]->heavy)) {
            search_path[2]=nb[n];
            is_ringpath(search_path);
          }
        }
      }
    }
  }
}
}

void clear_rings() {
int i,j;
  n_rings = 0;
for(i=0;i<max_rings;i++) {
  for(j=0;j<max_ringsize;j++) *(*(ring+i)+j)=0;
}
  for (i=0;i<max_rings;i++) {  // new in v0.3
      ringprop[i]->size     = 0;
      ringprop[i]->arom     = false;
      ringprop[i]->envelope = false;
  }
  if (n_atoms > 0) {
      for (i=1;i<=n_atoms;i++){atom[i]->ring_count= 0;}
  }
  if (n_bonds > 0) {
      for (i=0;i<n_bonds;i++){bond[i]->ring_count = 0;}
  }
}

int ring_lastpos(RINGPATH_TYPE s){
int i, rc, rlp;
  rlp = 0;
  if (n_rings > 0) {
      for (i =0; i<n_rings;i++) {
          rc = ringcompare(s, ring[i]);
          if (rc_identical(rc)) { rlp = i;}
       }
   }
return rlp;
}

void remove_redundant_rings() {
int i, j, k, rlp;
int *tmp_path; // : ringpath_type;
  if (n_rings>1) {
    for (i=0;i<(n_rings-1);i++) {
      tmp_path= *(ring+i);
      rlp = ring_lastpos(tmp_path);
      while (rlp > i) {
          for (j=rlp; j<(n_rings-1); j++) {
              ring[j] = ring[(j+1)];
              ringprop[j]->size = ringprop[(j+1)]->size;  // new in v0.3
              ringprop[j]->arom = ringprop[(j+1)]->arom;
              ringprop[j]->envelope = ringprop[(j+1)]->envelope;
           }
          for (k =0;k<max_ringsize;k++) {ring[n_rings-1][k] = 0;}
          n_rings--;
          rlp = ring_lastpos(tmp_path);
      }
    }
  }
}

int get_bond(int ba1, int ba2){
int i, b_id;
  b_id = 0;
  if (n_bonds > 0) {
    for (i=0; i<n_bonds; i++) {
        if (((bond[i]->a1 == ba1) && (bond[i]->a2 == ba2)) ||
           ((bond[i]->a1 == ba2) && (bond[i]->a2 == ba1)))  b_id = i;
    }
  }
return b_id;
}

void chk_envelopes() {  // new in v0.3d
// checks if a ring completely contains one or more other rings
int  a,i,j,k,l,pl,pli ;
bool  found_atom, found_all_atoms, found_ring;
if (n_rings > 1) {
  for (i=1;i<n_rings;i++) {
      found_ring = false;
      j = 0;
      pli = ringprop[i]->size;  // path_length(ring^[i]);
      while ((j < i) && (found_ring == false)) {
          found_all_atoms = true;
          pl = ringprop[j]->size;  // path_length(ring^[j]);
          for (k=0;k<pl;k++) {
              found_atom = false;
              a= *(*(ring+j)+k);
              for (l =0;l<pli;l++) {
                  if (*(*(ring+i)+l) == a) {found_atom = true;}
              }
              if (found_atom == false) {found_all_atoms = false;}
          }
          if(found_all_atoms){found_ring=true;}
          j++;
        }
      if (found_ring){ringprop[i]->envelope = true;}
  }
}
}

void update_ringcount() {
int i, j, a1, a2, b, pl;
  if (n_rings > 0) {
      chk_envelopes();
      for (i=0;i<n_rings;i++) {
        if (ringprop[i]->envelope == false) {
              pl = ringprop[i]->size;  // path_length(ring^[i]);  // v0.3d
              a2 = *(*(ring+i)+pl-1);
              for (j=0; j<pl; j++) {
                  a1= *(*(ring+i)+j);
                  atom[a1]->ring_count++;
                  b= get_bond(a1,a2);
                  bond[b]->ring_count++;
                  a2= a1;
              }
        }
     }
  }
}

bool is_oxo_C(int id) {
int  i, nb_count, k;
bool  r  ;
int  nb[max_neighbors];
r=false;
  if ((id < 1) ||(id > n_atoms)) {return r;}
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == id) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2 == id) && (nb_count < max_neighbors) && (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  if ((!strcmp(atom[id]->element,"C ")) && (atom[id]->neighbor_count > 0)) {
      for (i=0; i<atom[id]->neighbor_count; i++) {
          if ((bond[get_bond(id,nb[i])]->btype == 'D') &&
             (!strcmp(atom[(nb[i])]->element,"O ") )) {   // no N, amidines are different...
             r = true;
          }
      }
  }
  return r;
}

bool is_thioxo_C(int id) {
int  i, nb_count, k;
bool  r  ;
int  nb[max_neighbors];
r=false;
  if ((id < 1) || (id > n_atoms)) {return r;}
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == id) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== id)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  if ((!strcmp(atom[id]->element,"C ")) && (atom[id]->neighbor_count > 0)) {
      for (i=0; i<atom[id]->neighbor_count; i++) {
          if ((bond[get_bond(id,nb[i])]->btype == 'D') && (!strcmp(atom[(nb[i])]->element,"S ") ||
             !strcmp(atom[(nb[i])]->element,"SE"))) {   // no N, amidines are different...
             r = true;
          }
      }
  }
  return r;
}

bool is_true_imino_C(int id) {
int i, b, a_n, nb_count, k;
bool  r  ;
int  nb[max_neighbors];  // v0.3j
  r = false;
  a_n = 0;
  if ((id < 1) || (id > n_atoms)) return r;
  memset(nb, 0, max_neighbors*sizeof(int));
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == id) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== id)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  if (!strcmp(atom[id]->element,"C ") && (atom[id]->neighbor_count > 0)) {
      for (i=0;i<atom[id]->neighbor_count;i++) {
          b = get_bond(id,nb[i]);  // v0.3j
          if ((bond[b]->btype == 'D') && (bond[b]->arom == false) &&  // v0.3j
             !strcmp(atom[(nb[i])]->element,"N ")) a_n = nb[i];
      }
      if (a_n > 0) {
          r = true;
  memset(nb, 0, max_neighbors*sizeof(int));
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_n) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_n)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
          for (i=0; i<atom[a_n]->neighbor_count; i++) {
              if (strcmp(atom[(nb[i])]->element,"C ") && strcmp(atom[(nb[i])]->element,"H ") &&
              strcmp(atom[(nb[i])]->element,"D ")) r = false; // v0.3n: D
          }
      }
  }
  return r;
}

bool is_diazonium(int a_view, int a_ref) {
bool r ;
int  nb_next[max_neighbors];
int  bond_count, nb_next_count;
int  i,chg_count;
int  n1, n2;
  r = false;
  bond_count = 0;
  chg_count = 0;
  n1 = 0; n2 = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->element,"N ") && (atom[a_ref]->neighbor_count == 2)) {
          n1 = a_ref;
          chg_count = atom[n1]->formal_charge;
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == n1)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == n1)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          if (!strcmp(atom[(nb_next[0])]->element,"N ")) {
              n2 = nb_next[0];
              chg_count = chg_count + atom[n2]->formal_charge;                      
              if (bond[get_bond(n1,n2)]->btype == 'S') bond_count++;
              if (bond[get_bond(n1,n2)]->btype == 'D') bond_count=bond_count+2;
              if (bond[get_bond(n1,n2)]->btype == 'T') bond_count=bond_count+3;
          }
          if ((n1 > 0) && (n2 > 0) && (atom[n2]->neighbor_count == 1) &&
             (bond_count >= 2) && (chg_count > 0)) r = true;
        }
    }
  return r;
}

bool is_exocyclic_imino_C(int id,int r_id){
int  i,j,k,nb_count;
bool  r;
int  nb[max_neighbors];
RINGPATH_TYPE  testring;
int  ring_size;
  r = false;
  if ((id < 1) || (id > n_atoms)) return r;
  memset(nb, 0, max_neighbors * sizeof(int));
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == id) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
     nb[nb_count]=bond[k]->a2;
     nb_count++;
    }
    if ((bond[k]->a2== id)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  } //NEIGHBOR copy the neighbor part here!!
  memset(testring,0,sizeof(RINGPATH_TYPE));
  ring_size = ringprop[r_id]->size;  // v0.3j
  for (j=0;j<ring_size;j++) testring[j]= *(*(ring+r_id)+j);  // v0.3j
  if (!strcmp(atom[id]->element,"C ") && (atom[id]->neighbor_count > 0)) {
      for (i=0;i<atom[id]->neighbor_count;i++) {
          if ((bond[get_bond(id,nb[i])]->btype == 'D') && !strcmp(atom[(nb[i])]->element,"N ")) {
                 r = true;
               for (j=0;j<ring_size;j++){
                   if (nb[i] == *(*(ring+r_id)+j)) r = false;
               }
          }
      }
  }
  return r;
}

int find_exocyclic_methylene_C(int id, int r_id) { // renamed and rewritten in v0.3j
int  i,j,k,nb_count,r;
int  nb[max_neighbors];
RINGPATH_TYPE  testring;
int  ring_size;
  r = 0;
  if ((id < 1) || (id > n_atoms)) {return 0;}
      memset(nb, 0, max_neighbors * sizeof(int));
      nb_count=0;
      for (k= 0;k< n_bonds;k++ ){
        if ((bond[k]->a1 == id) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
          nb[nb_count]=bond[k]->a2;
          nb_count++;
        }
        if ((bond[k]->a2== id)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
          nb[nb_count]=bond[k]->a1;
          nb_count++;
        }
      } //NEIGHBOR copy the neighbor part here!!
  memset(testring,0,sizeof(RINGPATH_TYPE));
  ring_size = ringprop[r_id]->size;  // v0.3j
  for (j=0;j<ring_size;j++) testring[j] = *(*(ring+r_id)+j);  // v0.3j
  if (!strcmp(atom[id]->element,"C ") && (atom[id]->neighbor_count > 0)) {
      for (i=0;i<atom[id]->neighbor_count;i++) {
          if ((bond[get_bond(id,nb[i])]->btype == 'D') && !strcmp(atom[(nb[i])]->element,"C ")) {
                 r = nb[i];
                 for (j =0; j<ring_size; j++) {
                   if (nb[i] == *(*(ring+r_id)+j)) r = 0;
                 }
          } 
      }
  }
  return r;
}

void update_atypes() {
int i, j, k, b_id;
int nb[max_neighbors];
int single_count, double_count, triple_count, arom_count, acyl_count;
int C_count, O_count, nb_count;
int total_bonds;
int NdO_count;
int NdC_count;
int Htotal;
  if(n_atoms < 1){return;}
  for (i=1;i<=n_atoms;i++) {
      single_count = 0;
      double_count = 0;
      triple_count = 0;
      arom_count   = 0;
      total_bonds  = 0;
      acyl_count   = 0;
      C_count      = 0;
      O_count      = 0;
      NdO_count = 0;
      NdC_count = 0;
      Htotal    = 0;
      memset(nb, 0, max_neighbors * sizeof(int));
      nb_count=0;
      for (k= 0;k< n_bonds;k++ ){
        if ((bond[k]->a1 == i) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
          nb[nb_count]=bond[k]->a2;
          nb_count++;
        }
        if ((bond[k]->a2== i)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
          nb[nb_count]=bond[k]->a1;
          nb_count++;
        }
      } //NEIGHBOR copy the neighbor part here!!

      if (atom[i]->neighbor_count > 0){ // count single, double, triple, and aromatic bonds to all neighbor atoms
          for (j= 0;j<atom[i]->neighbor_count;j++){
              if ((is_oxo_C(nb[j])) || (is_thioxo_C(nb[j]))) {acyl_count++;}
              if (!strcmp(atom[(nb[j])]->element,"C ")){C_count++;}
              if (!strcmp(atom[(nb[j])]->element,"O ")) {O_count++;}
              b_id = get_bond(i,nb[j]);
              if (b_id >= 0){
                  if (bond[b_id]->btype == 'S') {single_count++;}
                  if (bond[b_id]->btype == 'D') {double_count++;}
                  if (bond[b_id]->btype == 'T') {triple_count++;}
                  if (bond[b_id]->btype == 'A') {    // v0.3n: special treatment for acyclic bonds
                                               // flagged as "aromatic" (in query structures)
                  if (bond[b_id]->ring_count > 0) {arom_count++;} else {double_count++;}
                  }
                  if ((!strcmp(atom[i]->element,"N ") && !strcmp(atom[(nb[j])]->element,"O ")) ||
                     (!strcmp(atom[i]->element,"O ") && !strcmp(atom[(nb[j])]->element,"N "))) 
                  {           // check if it is an N-oxide drawn with a double bond ==> should be N3
                       if (bond[b_id]->btype == 'D') {NdO_count++;}
                  }
                  if ((!strcmp(atom[i]->element,"N ") && !strcmp(atom[(nb[j])]->element,"C ")) ||
                     (!strcmp(atom[i]->element,"C ") && !strcmp(atom[(nb[j])]->element,"N "))) {
                       if (bond[b_id]->btype == 'D'){NdC_count++;}
                  }
              }
          }
          total_bonds = single_count + 2*double_count + 3*triple_count + (int)(1.5*arom_count);  
          // calculate number of total hydrogens per atom
          //Htotal := nvalences(atom^[i].element) - total_bonds + atom^[i].formal_charge;
          Htotal = atom[i]->nvalences - total_bonds + atom[i]->formal_charge;
          if (Htotal < 0) {Htotal = 0;}  // e.g., N in nitro group
          atom[i]->Htot = Htotal;
          // refine atom types, based on bond types
          if (!strcmp(atom[i]->element,"C ")) {
              if (arom_count > 1) {strcpy(atom[i]->atype,"CAR");}
              if ((triple_count == 1) || (double_count == 2)) {strcpy(atom[i]->atype,"C1 ");}
              if (double_count == 1) {strcpy(atom[i]->atype,"C2 ");}
              if ((triple_count == 0) && (double_count == 0) && (arom_count < 2)) {
                 strcpy(atom[i]->atype,"C3 ");
              }
//              printf("TEST ATYPE: %d, %s, %s\n",i,atom[i]->element, atom[i]->atype);
          }  
          if (!strcmp(atom[i]->element,"O ")) {
              if (double_count == 1) {strcpy(atom[i]->atype,"O2 ");}
              if (double_count == 0) {strcpy(atom[i]->atype,"O3 ");}
          }
          if (!strcmp(atom[i]->element,"N ")) {
              if (total_bonds > 3){
                  if (O_count == 0) {
                      if ((single_count > 3) ||((single_count == 2) && (double_count == 1) && (C_count >=2))) { atom[i]->formal_charge = 1;} 
                  } else { // could be an N-oxide -> should be found elsewhere 
                      if ((O_count == 1) && (atom[i]->formal_charge == 0)) { 
                        strcpy(atom[i]->atype,"N3 ");}  // v0.3m
                      if ((O_count == 2) && (atom[i]->formal_charge == 0)) { 
                          if (atom[i]->neighbor_count > 2) {strcpy(atom[i]->atype,"N2 ");}  // nitro v0.3o
                          if (atom[i]->neighbor_count == 2) {strcpy(atom[i]->atype,"N1 ");}  // NO2   v0.3o
                       }
                      // the rest is left empty, so far....
                  }
              }
              if ((triple_count == 1) || ((double_count == 2) && (atom[i]->neighbor_count == 2))) {
                   strcpy(atom[i]->atype,"N1 ");} // v0.3n
              if (double_count == 1) {
                  //if NdC_count > 0 then atom^[i].atype := 'N2 ';
                  if ((NdC_count == 0) && (NdO_count > 0) && (C_count >= 2)){
                    strcpy(atom[i]->atype,"N3 ");}  // N-oxide is N3 except in hetarene etc.
                  else {strcpy(atom[i]->atype,"N2 ");}                   // fallback, added in v0.3g 
               }
              if ((arom_count > 1) || (atom[i]->arom == true)){strcpy(atom[i]->atype,"NAR");}  // v0.3n
              if ((triple_count == 0) && (double_count == 0)) { 
                  if ((atom[i]->formal_charge == 0)) {
                      if (acyl_count ==0){strcpy(atom[i]->atype,"N3 ");}
                      if (acyl_count > 0) {strcpy(atom[i]->atype,"NAM");}
                  }
                  if (atom[i]->formal_charge == 1) {strcpy(atom[i]->atype,"N3+");}
               }
          } 
          if (!strcmp(atom[i]->element,"P ")) {
              if (single_count > 4) {strcpy(atom[i]->atype,"P4 ");}
              if ((single_count <= 4) && (double_count == 0)) {strcpy(atom[i]->atype,"P3 ");}
              if (double_count == 2) {strcpy(atom[i]->atype,"P3D");}
           }
          if (!strcmp(atom[i]->element,"S ")) {
              if ((double_count == 1) && (single_count == 0)) { strcpy(atom[i]->atype,"S2 ");}
              if (double_count == 0) {strcpy(atom[i]->atype,"S3 ");}
              if ((double_count == 1) && (single_count > 0)) {strcpy(atom[i]->atype,"SO ");}
              if ((double_count == 2) && (single_count > 0)) {strcpy(atom[i]->atype,"SO2");}
          }
          // further atom types should go here
       }
  }
}

void update_Htotal() {
int  i, i2, j, k, b_id;
int  nb[max_neighbors], nb_next[max_neighbors];
int  single_count, double_count, triple_count, arom_count;
int  nb_count, nb_next_count, total_bonds;
int  Htotal ;
int  nval;   // new in v0.3
bool  diazon;       // new in v0.3j
int  a1, a2, a3;   // new in v0.3j
  if (n_atoms < 1) {return;}
  diazon = false;
  for (i = 1; i<=n_atoms; i++) {
      single_count = 0;
      double_count = 0;
      triple_count = 0;
      arom_count   = 0;
      total_bonds  = 0;
      Htotal    = 0;
      memset(nb, 0, max_neighbors * sizeof(int));
      nb_count=0;
      for (k= 0;k< n_bonds;k++ ){
        if ((bond[k]->a1 == i) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
          nb[nb_count]=bond[k]->a2;
          nb_count++;
        }
        if ((bond[k]->a2== i)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
          nb[nb_count]=bond[k]->a1;
          nb_count++;
        }
      } //NEIGHBOR copy the neighbor part here!!

      if (atom[i]->neighbor_count > 0) {  // count single, double, triple, and aromatic bonds to all neighbor atoms
          for (j=0; j<atom[i]->neighbor_count; j++) {
              b_id = get_bond(i,nb[j]);
              if (b_id >= 0) {
                  if (bond[b_id]->btype == 'S') single_count++;
                  if (bond[b_id]->btype == 'D') double_count++;
                  if (bond[b_id]->btype == 'T') triple_count++;
                  if (bond[b_id]->btype == 'A') arom_count++;
              }
          }
          //check for diazonium salts
          a1 = i; //id
          a2 = nb[0]; //prev_id
          if (!strcmp(atom[a1]->element,"N ") && !strcmp(atom[a2]->element,"N ")) {
              if (atom[a2]->neighbor_count == 2) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i2=0; i2<n_bonds;i2++) {
        if ((bond[i2]->a1 == a2)&& (bond[i2]->a2 != a1) && (nb_next_count < max_neighbors) 
        && (atom[bond[i2]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i2]->a2;
          nb_next_count++;
        }
        if ((bond[i2]->a2 == a2)&& (bond[i2]->a1 != a1) && (nb_next_count < max_neighbors)  
        && (atom[bond[i2]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i2]->a1;
          nb_next_count++;
        }
      }
                  a3= nb_next[0];
                  if (!strcmp(atom[a3]->element,"C ") && is_diazonium(a3,a2)) diazon = true;
              }
          }
      }
      total_bonds = single_count + 2*double_count + 3*triple_count + (int)(1.5*arom_count);  
      // calculate number of total hydrogens per atom
      //nval := nvalences(atom^[i].element);    // new in v0.3
      nval = atom[i]->nvalences;    // new in v0.3m
      if (!strcmp(atom[i]->element,"P ")) { 
          if ((total_bonds-atom[i]->formal_charge) > 3) nval = 5;  // refined in v0.3n
      }  
      if (!strcmp(atom[i]->element,"S ")) {      // v0.3h
          if ((total_bonds > 2) && (atom[i]->formal_charge < 1)) nval = 4;  // updated in v0.3j
          if (total_bonds > 4) nval=6;  // this will need some refinement...
      }
      Htotal = nval - total_bonds + atom[i]->formal_charge;
      if ((atom[i]->radical_type == 1) || (atom[i]->radical_type == 3)) Htotal = Htotal - 2; // v0.3p
      if (atom[i]->radical_type == 2) Htotal = Htotal - 1; // v0.3p
      if (diazon) Htotal = 0;      // v0.3j
      if (Htotal < 0) Htotal = 0;  // e.g., N in nitro group
      atom[i]->Htot = Htotal;
      if (atom[i]->Hexp > atom[i]->Htot) atom[i]->Htot = atom[i]->Hexp;  // v0.3n; just to be sure...
    }
}

void chk_arom() {
int  i, j, pi_count, ring_size;
int  b, a1, a2;  // v0.3n
RINGPATH_TYPE  testring;
int  a_ref, a_prev, a_next;
int  b_bk, b_fw, b_exo;
char bt_bk, bt_fw;
bool  ar_bk, ar_fw, ar_exo ;  // new in v0.3
bool  conj_intr, ko, aromatic;
bool  aromatic_bt ;  // v0.3n
int  n_db, n_sb, n_ar ;
bool  cumul ;
int  exo_mC ;        // v0.3j
int  arom_pi_diff ;  // v0.3j
  if (n_rings < 1) return;     
  // first, do a very quick check for benzene, pyridine, etc.
  for (i =0; i<n_rings; i++) {
      ring_size = ringprop[i]->size;
      if (ring_size == 6) {
          memset(testring,0,sizeof(RINGPATH_TYPE));
          for(j=0;j<ring_size;j++) {testring[j] = *(*(ring+i)+j);}
          cumul = false;
          n_sb = 0;
          n_db = 0;
          n_ar = 0;
          a_prev = testring[ring_size-1];
          for (j = 0;j<ring_size;j++) {
              a_ref = testring[j];
              if (j < ring_size-1) {a_next = testring[j+1];} else {a_next = testring[0];}
              b_bk  = get_bond(a_prev,a_ref);
              b_fw  = get_bond(a_ref,a_next);
              bt_bk = bond[b_bk]->btype;
              bt_fw = bond[b_fw]->btype;
              if (bt_fw == 'S') n_sb++;
              if (bt_fw == 'D') n_db++;
              if (bt_fw == 'A') n_ar++;
              if ((bt_fw != 'A') && (bt_bk == bt_fw)) cumul = true;
              a_prev = a_ref;
          }
          if ((n_ar == 6) || ((n_sb == 3) && (n_db == 3) && (cumul == false))) { // this ring is aromatic
              a_prev = testring[ring_size-1];
              for (j =0; j<ring_size;j++) {
                  a_ref = testring[j];
                  b_bk  = get_bond(a_prev,a_ref);
                  bond[b_bk]->arom = true;
                  a_prev = a_ref;
              }
              ringprop[i]->arom = true;            
          }
      }
  } 
  for (i = 0; i<n_rings; i++) {
      if (ringprop[i]->arom == false) {// do the hard work only for those rings which are not yet flagged aromatic
          memset(testring,0,sizeof(RINGPATH_TYPE));
          ring_size = ringprop[i]->size;  // v0.3j
          for (j= 0;j<ring_size;j++) testring[j] = *(*(ring+i)+j);  // v0.3j
          pi_count = 0;
          arom_pi_diff = 0;  // v0.3j
          conj_intr = false;
          ko        = false;
          a_prev    = testring[ring_size-1];
          for (j = 0;j<ring_size;j++){
              a_ref = testring[j];
              if (j < ring_size-1) {a_next= testring[j+1];} else {a_next= testring[0];}
              b_bk  = get_bond(a_prev,a_ref);
              b_fw  = get_bond(a_ref,a_next);
              bt_bk = bond[b_bk]->btype;
              bt_fw = bond[b_fw]->btype;
              ar_bk = bond[b_bk]->arom;
              ar_fw = bond[b_fw]->arom;
              if ((bt_bk == 'S') && (bt_fw == 'S') && (ar_bk == false) && (ar_fw == false)) {
                  // first, assume the worst case (interrupted conjugation)
                  conj_intr= true;  
                  // conjugation can be restored by hetero atoms
                  if (!strcmp(atom[a_ref]->atype,"O3 ") || !strcmp(atom[a_ref]->atype,"S3 ") ||
                     !strcmp(atom[a_ref]->element,"N ") || !strcmp(atom[a_ref]->element,"SE") ||
                     !strcmp(atom[a_ref]->element,"Q ")) {    // v0.4: query atom of type Q is OK
                       conj_intr = false;
                       pi_count=pi_count+2;  // lone pair adds for 2 pi electrons
                     }
                  // conjugation can be restored by a formal charge at a methylene group
                  if (!strcmp(atom[a_ref]->element,"C ") && (atom[a_ref]->formal_charge != 0)) {
                      conj_intr = false;
                      pi_count  = pi_count - atom[a_ref]->formal_charge;  // neg. charge increases pi_count!
                  }
                  // conjugation can be restored by carbonyl groups etc.
                  if ((is_oxo_C(a_ref)) || (is_thioxo_C(a_ref)) || (is_exocyclic_imino_C(a_ref,i)))
                  {conj_intr= false;}
                  // conjugation can be restored by exocyclic C=C double bond,
                  // adds 2 pi electrons to 5-membered rings, not to 7-membered rings (CAUTION!)
                  // apply only to non-aromatic exocyclic C=C bonds
                  exo_mC = find_exocyclic_methylene_C(a_ref,i);  // v0.3j
                  if ((exo_mC > 0) && odd(ring_size)) {
                      b_exo = get_bond(a_ref,exo_mC);           // v0.3j 
                      ar_exo = bond[b_exo]->arom;
                      if ((ring_size - 1)%4 == 0) {  // 5-membered rings and related
                          conj_intr = false;
                          pi_count=pi_count+2;
                      } else {                          // 7-membered rings and related
                          if (!ar_exo) conj_intr= false;
                      }
                  }
                  // if conjugation is still interrupted ==> knock-out
                  if (conj_intr) ko = true;
                } else {
                  if ((bt_bk == 'S') && (bt_fw == 'S') && (ar_bk == true) && (ar_fw == true)) {
                      if (!strcmp(atom[a_ref]->atype,"O3 ") || !strcmp(atom[a_ref]->atype,"S3 ") ||
                         !strcmp(atom[a_ref]->element,"N ") || !strcmp(atom[a_ref]->element,"SE")) {
                           pi_count=pi_count+2;  // lone pair adds for 2 pi electrons
                      }
                      if (!strcmp(atom[a_ref]->element,"C ") && (atom[a_ref]->formal_charge != 0)) 
                      {pi_count = pi_count - atom[a_ref]->formal_charge;}  // neg. charge increases pi_count!
                      exo_mC = find_exocyclic_methylene_C(a_ref,i);  // v0.3j
                      if ((exo_mC > 0) && odd(ring_size)) {       // v0.3j
                          b_exo = get_bond(a_ref,exo_mC);            // v0.3j
                          ar_exo = bond[b_exo]->arom;
                          if ((ring_size - 1)%4 == 0) pi_count=pi_count+2;  // 5-membered rings and related
                      }
                   } else {   // any other case: increase pi count by one electron
                      pi_count++;  // v0.3j; adjustment for bridgehead N: see below
                      if (((bt_bk == 'S') && (bt_fw == 'S')) &&
                          (((ar_bk == true) && (ar_fw == false)) ||
                           ((ar_bk == false) && (ar_fw == true)))) { 
                          // v0.3j; if a bridgehead N were not aromatic, it could 
                          // contribute 2 pi electrons --> try also this variant
                          // (example: CAS 32278-54-9)
                          if (!strcmp(atom[a_ref]->element,"N ")) arom_pi_diff++;
                      }
                   }
                }
              // last command:
              a_prev = a_ref;
            }  // for j := 1 to ring_size
          // now we can draw our conclusion
          //if not ((ko) or (odd(pi_count))) then
          if (!ko) {   // v0.3j; odd pi_count might be compensated by arom_pi_diff
          // apply Hueckel's rule
              if ((abs(ring_size-pi_count) < 2) && (((pi_count - 2)%4 == 0) ||
                  (((pi_count + arom_pi_diff) - 2)%4 == 0))) {
                  // this ring is aromatic
                  ringprop[i]->arom= true;
                  // now mark _all_ bonds in the ring as aromatic
                  a_prev= testring[ring_size-1];
                  for(j=0; j<ring_size; j++) {
                      a_ref = testring[j];
                      bond[get_bond(a_prev,a_ref)]->arom= true;
                      a_prev= a_ref;
                  }
              }
          }
      }
  }  // (for i := 1 to n_rings)
  // finally, mark all involved atoms as aromatic
  for (i=1; i<=n_atoms; i++) {atom[i]->arom=false;} //initialize atom->arom
  for (i=0; i<n_bonds; i++) {
      if (bond[i]->arom) {
          a1 = bond[i]->a1;  // v0.3n
          a2 = bond[i]->a2;  // v0.3n
          atom[a1]->arom = true;
          atom[a2]->arom = true;
          // v0.3n: update atom types if applicable (C and N)
          if (!strcmp(atom[a1]->element,"C ")) strcpy(atom[a1]->atype,"CAR");
          if (!strcmp(atom[a2]->element,"C ")) strcpy(atom[a2]->atype,"CAR");
          if (!strcmp(atom[a1]->element,"N ")) strcpy(atom[a1]->atype,"NAR");
          if (!strcmp(atom[a2]->element,"N ")) strcpy(atom[a2]->atype,"NAR");
      }
  }
  // update aromaticity information in ringprop
  // new in v0.3n: accept rings as aromatic if all bonds are of type 'A'
  for (i=0;i<n_rings;i++){
      ring_size = ringprop[i]->size;
      memset(testring,0,sizeof(RINGPATH_TYPE));
      for(j=0;j<ring_size;j++) {testring[j] = *(*(ring+i)+j);}
      //ring_size := path_length(testring);
      ring_size = ringprop[i]->size;  // v0.3j
      aromatic = true;
      aromatic_bt = true;  // v0.3n
      a_prev = testring[ring_size-1];
      for (j=0; j<ring_size; j++) {
          a_ref = testring[j];
          b = get_bond(a_prev,a_ref);  // v0.3n
          if (!bond[b]->arom) {aromatic = false;}        
          if (bond[b]->btype != 'A') aromatic_bt = false;  // v0.3n
          a_prev = a_ref;
      }

      if ((aromatic_bt) && (!aromatic)) {  // v0.3n: update aromaticity flag
          a_prev = testring[ring_size-1];
          for (j =0;j<ring_size;j++) {
              a_ref = testring[j];
              b = get_bond(a_prev,a_ref);
              bond[b]->arom = true;
              if (!strcmp(atom[a_ref]->element,"C ")) strcpy(atom[a_ref]->atype,"CAR");
              if (!strcmp(atom[a_ref]->element,"N ")) strcpy(atom[a_ref]->atype,"NAR");
              a_prev = a_ref;
          }
          aromatic = true;
      }                                  // end v0.3n block  
      if (aromatic) {ringprop[i]->arom = true;   } else {ringprop[i]->arom = false;}
   }  
}

int count_aromatic_rings() {
int i, n;
  n = 0;
    if (n_rings > 0) {
        for (i=0;i<n_rings;i++){
          if (ringprop[i]->arom) n++;
        }
    }
  return n;
}

void write_mol() {
int  i, j;
RINGPATH_TYPE  testring;
int  ring_size;
//  if progmode = pmCheckMol then
  printf("Molecule name: %s\n",molname);
//  else
//    writeln('Molecule name (haystack): ',molname);
  printf("atoms: %d bonds: %d rings: %d\n",n_atoms,n_bonds,n_rings);
  if (n_atoms < 1) return;
  if (n_bonds < 1) return;
  for (i = 1; i<=n_atoms; i++) {
    if (i <   10) printf(" ");
    if (i <  100) printf(" ");
    if (i < 1000) printf(" ");
    printf("%d %s %s %9.4f %9.4f %9.4f",i,atom[i]->element,atom[i]->atype,atom[i]->x,atom[i]->y,atom[i]->z);
    printf("  (%d heavy-atom neighbors, Hexp: %d, Htot: %d)",atom[i]->neighbor_count,atom[i]->Hexp,atom[i]->Htot);
    if (atom[i]->formal_charge != 0) printf("  charge: %d",atom[i]->formal_charge);
    if (atom[i]->arom) printf(" aromatic");
    printf("\n");
  }
  for (i=0; i<n_bonds;i++) {
      if (i <   10) printf(" ");
      if (i <  100) printf(" ");
      if (i < 1000) printf(" ");
      printf("%d %d %d %c",i,bond[i]->a1,bond[i]->a2,bond[i]->btype);
      if (bond[i]->ring_count > 0) printf(", contained in %d rings",bond[i]->ring_count);
      if (bond[i]->arom) printf(" (aromatic) ");
      printf("\n");
  }
  if (n_rings > 0) {
      for (i =0; i<n_rings; i++) {
          printf("ring %d: ",i);
          //aromatic := true;
          memset(testring,0,sizeof(RINGPATH_TYPE));
          ring_size=ringprop[i]->size;  // v0.3j
          //for j := 1 to max_ringsize do if ring^[i,j] > 0 then testring[j] := ring^[i,j];
          for (j=0; j<ring_size; j++) {testring[j]= *(*(ring+i)+j);}  // v0.3j
          //ring_size := path_length(testring);
          //a_prev := testring[ring_size];
          for (j=0; j<ring_size;j++) {
              printf("%d ",testring[j]);
              //a_ref := testring[j];
              //if (not bond^[get_bond(a_prev,a_ref)].arom) then aromatic := false;
              //a_prev := a_ref;
          }
          //if aromatic then write(' (aromatic)');
          if (ringprop[i]->arom) printf(" (aromatic)");
          if (ringprop[i]->envelope) printf(" (env)");
          printf("\n");
       }
  }
}

void add2fgloc(int row, int a) {  // v0.5
int  i, n ;
bool  a_found;
  if ((row < 0) || (row > used_fg)) return;
  n = *(*(fgloc+row)+0);
  a_found = false;
  if (n > 0) {
      for (i = 0; i< n;i++) { 
          if (*(*(fgloc+row)+i) == a) a_found = true;
      }
  }
  if ((a_found = false) && (n < max_fgpos)) {
      n++;
      *(*(fgloc+row)+0)= n;
      *(*(fgloc+row)+n)= a;
  }
}

bool is_electroneg(char a_el[3]){
// new in v0.3j
bool  res = false;
  if (!strcmp(a_el,"N ")) res = true;
  if (!strcmp(a_el,"P ")) res = true;
  if (!strcmp(a_el,"O ")) res = true;
  if (!strcmp(a_el,"S ")) res = true;
  if (!strcmp(a_el,"SE")) res = true;
  if (!strcmp(a_el,"TE")) res = true;
  if (!strcmp(a_el,"F ")) res = true;
  if (!strcmp(a_el,"CL")) res = true;
  if (!strcmp(a_el,"BR")) res = true;
  if (!strcmp(a_el,"I ")) res = true;
  return res;
}

int raw_hetbond_count(int a){   // new in v0.2j, ignores bond order
int  i,nb_count,k ;
int  nb[max_neighbors];
int  hbc;
  hbc = 0;
  if ((a > 0) && (a <= n_atoms)) {
      memset(nb, 0, max_neighbors * sizeof(int));
      nb_count=0;
      for (k= 0;k< n_bonds;k++ ){
        if ((bond[k]->a1 == a) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
          nb[nb_count]=bond[k]->a2;
          nb_count++;
        }
        if ((bond[k]->a2== a)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
          nb[nb_count]=bond[k]->a1;
          nb_count++;
        }
      } //NEIGHBOR copy the neighbor part here!!
      if (atom[a]->neighbor_count > 0) {
          for (i= 0; i<atom[a]->neighbor_count; i++) {
  // added 'D ' in v0.3n
              if (strcmp(atom[(nb[i])]->element,"C ") && strcmp(atom[(nb[i])]->element,"A ") &&
              strcmp(atom[(nb[i])]->element,"H ") && strcmp(atom[(nb[i])]->element,"D ") &&
              strcmp(atom[(nb[i])]->element,"LP") && strcmp(atom[(nb[i])]->element,"DU")) hbc++;
          }
      }
  }
  return hbc;
}

int hetbond_count(int a) {
int  i, nb_count,k;
int  nb[max_neighbors];
char  bt;   // v0.4
float  hbc;
int  hbc_int;  // v0.4
  hbc = 0.0;
  if ((a > 0) && (a <= n_atoms)) {
      memset(nb, 0, max_neighbors * sizeof(int));
      nb_count=0;
      for (k= 0;k< n_bonds;k++ ){
        if ((bond[k]->a1 == a) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
          nb[nb_count]=bond[k]->a2;
          nb_count++;
        }
        if ((bond[k]->a2== a)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
          nb[nb_count]=bond[k]->a1;
          nb_count++;
        }
      } //NEIGHBOR copy the neighbor part here!!
      if (atom[a]->neighbor_count > 0) {
          for (i=0; i<atom[a]->neighbor_count; i++) {
              if (strcmp(atom[(nb[i])]->element,"C ") && strcmp(atom[(nb[i])]->element,"A ") &&
              strcmp(atom[(nb[i])]->element,"H ") && strcmp(atom[(nb[i])]->element,"D ") &&
              strcmp(atom[(nb[i])]->element,"LP") && strcmp(atom[(nb[i])]->element,"DU")) {  // added 'D ' in v0.3n
                  bt = bond[(get_bond(a,nb[i]))]->btype;  // v0.4
                  if (bt == 'S') hbc = hbc + 1;
                  if (bt == 'A') hbc = hbc + 1.5;
                  if (bt == 'D') hbc = hbc + 2;
                  if (bt == 'T') hbc = hbc + 3;
              }
          }
      }
  }
  hbc_int = (int)(hbc+0.5);  // v0.4
  return hbc_int;  // v0.4
}

bool is_arene(int r_id) {
int  i,j  ;
bool  r;
RINGPATH_TYPE  testring ;
int  ring_size;
int  a_prev, a_ref;
  r = false;
  if ((r_id < 0) || (r_id > n_rings-1)) return false;
  memset(testring,0,sizeof(RINGPATH_TYPE));
  ring_size= ringprop[r_id]->size;  // v0.3j
  //for j := 1 to max_ringsize do if ring^[r_id,j] > 0 then testring[j] := ring^[r_id,j];
  for (j=0; j<ring_size;j++) testring[j] = *(*(ring+r_id)+j);  // v0.3j
  //ring_size := path_length(testring);
  if (ring_size > 2) {
      r = true;
      a_prev = testring[ring_size-1];
      for (i =0;i<ring_size;i++){
          a_ref = testring[i];
          if (bond[get_bond(a_prev,a_ref)]->arom == false) {r= false;}
          a_prev = a_ref;
      }
  }
  return r;
}

bool is_heterocycle(int r_id) {
int  i,j;
bool  r;
RINGPATH_TYPE  testring;
int  ring_size;
int  a_ref;
  r= false;
  if ((r_id < 0) || (r_id >= n_rings)) return r; 
  memset(testring,0,sizeof(RINGPATH_TYPE));
  ring_size= ringprop[r_id]->size;  // v0.3j
  //for j := 1 to max_ringsize do if ring^[r_id,j] > 0 then testring[j] := ring^[r_id,j];
  for (j =0;j<ring_size;j++) testring[j] =  *(*(ring+r_id)+j);  // v0.3j
  //ring_size := path_length(testring);
  if (ring_size > 2) {
      for (i=0;i<ring_size;i++) {
          a_ref = testring[i];
          if (strcmp(atom[a_ref]->element,"C ")) r=true;
      }
  }
  return r;
}

bool is_hydroxy(int a_view, int a_ref) {
bool  r=false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->atype,"O3 ") && (atom[a_ref]->neighbor_count == 1)) r = true;
  }
  return r;
}

bool is_amino(int a_view, int a_ref) {
bool  r = false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if ((!strcmp(atom[a_ref]->atype,"N3 ") || !strcmp(atom[a_ref]->atype,"N3+")) &&
        (atom[a_ref]->neighbor_count == 1)) r= true;
  }
  return r;
}

bool is_alkyl(int a_view, int a_ref) {
int  i;
bool  r;
int  nb_next[max_neighbors];
int  het_count, nb_next_count;
  r = false;
  het_count = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S') &&
     !strcmp(atom[a_ref]->atype,"C3 ") && (atom[a_ref]->arom == false)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
      for (i=0; i<(atom[a_ref]->neighbor_count-1);i++) {
          if (strcmp(atom[(nb_next[i])]->element,"C ") && strcmp(atom[(nb_next[i])]->element,"H ") &&
            strcmp(atom[(nb_next[i])]->element,"D ") && strcmp(atom[(nb_next[i])]->element,"DU") && 
            strcmp(atom[(nb_next[i])]->element,"LP")) het_count++;  // added 'D ' in v0.3n
      }
      if (het_count <= 1) r = true;  // we consider (e.g.) alkoxyalkyl groups as alkyl
  }
  return r;
}

bool is_alkoxy(int a_view,int a_ref) {
bool  r  ;
int  nb_next_count, i ;
int nb_next[max_neighbors];
  r = false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->atype,"O3 ") && (atom[a_ref]->neighbor_count == 2)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          if (is_alkyl(a_ref,nb_next[0])) r = true;
      }
  }
  return r;
}

bool is_alkynyl(int a_view, int a_ref) {  // new in v0.3j
int  i, nb_next_count  ;
bool  r  ;
int  nb_next[max_neighbors];
int  c1_count;
  r = false;
  c1_count  = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S') &&
     !strcmp(atom[a_ref]->atype,"C1 ") && (atom[a_ref]->arom == false)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
      for (i=0;i<(atom[a_ref]->neighbor_count - 1);i++) {
          if (!strcmp(atom[(nb_next[i])]->atype,"C1 ")) c1_count++;  
      }
      if (c1_count == 1) r = true;
  }
  return r;
}

bool is_aryl(int a_view,int a_ref) {
bool  r= false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S') &&
     !strcmp(atom[a_ref]->element,"C ") && (atom[a_ref]->arom == true)) r = true;
  return r;
}

bool is_aryloxy(int a_view, int a_ref) {
bool  r ;
int  nb_next[max_neighbors];
int nb_next_count, i;
  r = false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->atype,"O3 ") && (atom[a_ref]->neighbor_count == 2)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          if (is_aryl(a_ref,nb_next[0])) r = true;
      }
  }
  return r;
}

bool is_alkylamino(int a_view, int a_ref) {
bool  r ;
int  nb_next[max_neighbors];
int nb_next_count, i;
int  alkyl_count ;
  r = false;
  alkyl_count = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->element,"N ") && (atom[a_ref]->neighbor_count == 2)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          if (is_alkyl(a_ref,nb_next[0])) alkyl_count++;
          if (alkyl_count == 1) r = true ; 
      }
  }
  return r;
}

bool is_dialkylamino(int a_view, int a_ref) {
bool  r ;
int  nb_next[max_neighbors];
int nb_next_count, i;
int  alkyl_count;
  r = false;
  alkyl_count = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->element,"N ") && (atom[a_ref]->neighbor_count == 3)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          for (i = 0; i< 2; i++) {
              if (is_alkyl(a_ref,nb_next[i])) alkyl_count++;
          }
          if (alkyl_count == 2) r = true;  
      }
  }
  return r;
}

bool is_arylamino(int a_view, int a_ref){
bool  r ;
int  nb_next[max_neighbors];
int nb_next_count, i;
int  aryl_count;
  r = false;
  aryl_count = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) { 
      if (!strcmp(atom[a_ref]->element,"N ") && (atom[a_ref]->neighbor_count == 2)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          if (is_aryl(a_ref,nb_next[0])) aryl_count++;
          if (aryl_count == 1) r = true;  
      }
  }
  return r;
}

bool is_diarylamino(int a_view, int a_ref) {
bool  r ;
int  nb_next[max_neighbors];
int nb_next_count, i;
int  aryl_count;
  r = false;
  aryl_count = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) { 
      if (!strcmp(atom[a_ref]->element,"N ") && (atom[a_ref]->neighbor_count == 3)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          for (i =0;i< 2; i++) {
              if (is_aryl(a_ref,nb_next[i])) aryl_count++;
          }
          if (aryl_count == 2) r = true;
      }
  }
  return r;
}

bool is_alkylarylamino(int a_view, int a_ref) {
bool  r ;
int  nb_next[max_neighbors];
int nb_next_count, i;
int  alkyl_count, aryl_count;
  r = false;
  alkyl_count = 0;  aryl_count  = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) { 
      if (!strcmp(atom[a_ref]->element,"N ") && (atom[a_ref]->neighbor_count == 3)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          for (i=0;i<2;i++) {
              if (is_alkyl(a_ref,nb_next[i])) alkyl_count++;
              if (is_aryl(a_ref,nb_next[i])) aryl_count++;
          }
          if ((alkyl_count == 1) && (aryl_count == 1))  r = true;  
      }
  }
return r;
}

bool is_thiocarbamoyl(int a_view, int a_ref) {
bool  r ;
int  nb_next[max_neighbors];
int nb_next_count, i;
  r= false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) { 
      if ((is_thioxo_C(a_ref)) && (atom[a_ref]->neighbor_count == 3)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          for (i=0;i<2;i++) {
           if (!strcmp(atom[(nb_next[i])]->atype,"N3 ") || !strcmp(atom[(nb_next[i])]->atype,"NAM")) r= true;
          }
      }
  }
return r;
}

bool is_true_alkyl(int a_view, int a_ref) {
bool  r ;
int  nb_next[max_neighbors];
int nb_next_count, i;
int  het_count;
  r= false;
  het_count = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S') &&
     !strcmp(atom[a_ref]->atype,"C3 ") && (atom[a_ref]->arom == false)) {
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
      for (i = 0;i<(atom[a_ref]->neighbor_count - 1);i++) {
        if (strcmp(atom[(nb_next[i])]->element,"C ") &&
            strcmp(atom[(nb_next[i])]->element,"H ") &&
            strcmp(atom[(nb_next[i])]->element,"D ") &&
            strcmp(atom[(nb_next[i])]->element,"DU")) het_count++;  // added 'D ' in v0.3n
      }
      if (het_count == 0) r= true;  //
  }
  return r;
}

bool is_true_alkoxy(int a_view, int a_ref) {
bool  r ;
int  nb_next[max_neighbors];
int nb_next_count, i;
  r = false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->atype,"O3 ") && (atom[a_ref]->neighbor_count == 2)) {
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          if (is_true_alkyl(a_ref,nb_next[0])) r= true;
      }
  }
  return r;
}

bool is_true_alkylsulfanyl(int a_view, int a_ref) {
bool  r ;
int  nb_next[max_neighbors];
int nb_next_count, i;
  r = false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->atype,"S3 ") && (atom[a_ref]->neighbor_count == 2)) {
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          if (is_true_alkyl(a_ref,nb_next[0])) r = true;
      }
  }
  return r;
}

bool is_true_alkylamino(int a_view, int a_ref) {
bool  r ;
int  nb_next[max_neighbors];
int nb_next_count, i, alkyl_count;
  r = false;
  alkyl_count = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if ((!strcmp(atom[a_ref]->atype,"N3 ") || !strcmp(atom[a_ref]->atype,"N3+")) &&
         (atom[a_ref]->neighbor_count == 2)) {
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          if (is_true_alkyl(a_ref,nb_next[0])) alkyl_count++;
          if (alkyl_count == 1) r = true;
      }
  }
  return r;
}

bool is_true_dialkylamino(int a_view, int a_ref) {
bool  r ;
int  nb_next[max_neighbors];
int nb_next_count, i, alkyl_count;
  r = false;
  alkyl_count = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if ((!strcmp(atom[a_ref]->atype,"N3 ") || !strcmp(atom[a_ref]->atype,"N3+")) &&
         (atom[a_ref]->neighbor_count == 3)) {
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          for (i = 0;i<2;i++) {
              if (is_true_alkyl(a_ref,nb_next[i])) alkyl_count++;
          }
          if (alkyl_count == 2)  r = true; 
      }
  }
  return r;
}

bool is_alkenyl(int a_view, int a_ref) {  // new in v0.3j; rewritten in v0.4b
int  i, nb_next_count;
bool  r;
int  nb_bond;       // v0.4b
int nb_next[max_neighbors];          // v0.4b
int  cc_dbl_count;  // v0.4b
  r = false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S') &&
     !strcmp(atom[a_ref]->atype,"C2 ") && (atom[a_ref]->arom == false)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
      cc_dbl_count = 0;
      for (i = 0;i<(atom[a_ref]->neighbor_count - 1); i++) {
          nb_bond = get_bond(a_ref,nb_next[i]);
          if ((bond[nb_bond]->btype == 'D') && !strcmp(atom[(nb_next[i])]->element,"C ")) {
            cc_dbl_count++;
          }
      }
      if (cc_dbl_count > 0) r = true;       // we consider (e.g.) alkoxyalkenyl groups as alkenyl
  }                                        // v0.3k: changed c2_count = 1 into c2_count >= 1
  return r;
}

bool is_alkenyloxy(int a_view, int a_ref) {  // v0.3j
bool  r;
int nb_next_count, i;
int nb_next[max_neighbors];
  r= false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->atype,"O3 ") && (atom[a_ref]->neighbor_count == 2)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          if (is_alkenyl(a_ref,nb_next[0])) r= true;
      }
  }
  return r;
}

bool is_alkynyloxy(int a_view, int a_ref) {  // v0.3j
bool  r  ;
int nb_next_count, i;
int nb_next[max_neighbors];
  r = false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->atype,"O3 ") && (atom[a_ref]->neighbor_count == 2)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          if (is_alkynyl(a_ref,nb_next[0])) r= true;
      }
  }
  return r;
}

bool is_alkylsulfanyl(int a_view, int a_ref){
bool  r ;
int nb_next_count, i;
int nb_next[max_neighbors];
  r = false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->atype,"S3 ") && (atom[a_ref]->neighbor_count == 2)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          if (is_alkyl(a_ref,nb_next[0])) r= true;
      }
  }
  return r;
}

bool is_sulfanyl(int a_view, int a_ref){
bool  r ;
  r = false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->atype,"S3 ") && (atom[a_ref]->neighbor_count == 1)) r = true;
  }
  return r;
}

bool is_arylsulfanyl(int a_view, int a_ref) {
bool  r ;
int  nb_next_count, i;
int  nb_next[max_neighbors];
  r = false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->atype,"S3 ") && (atom[a_ref]->neighbor_count == 2)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          if (is_aryl(a_ref,nb_next[0])) r= true;
      }
  }
  return r;
}

bool is_alkenylsulfanyl(int a_view, int a_ref) {  // v0.3j
bool  r ;
int  nb_next_count, i;
int  nb_next[max_neighbors];
  r = false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')){
      if (!strcmp(atom[a_ref]->atype,"S3 ") && (atom[a_ref]->neighbor_count == 2)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          if (is_alkenyl(a_ref,nb_next[0])) r = true;
      }
   }
  return r;
}

bool is_alkynylsulfanyl(int a_view,int a_ref) {  // v0.3j
bool  r ;
int  nb_next_count, i;
int  nb_next[max_neighbors];
  r = false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')){
      if (!strcmp(atom[a_ref]->atype,"S3 ") && (atom[a_ref]->neighbor_count == 2)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          if (is_alkynyl(a_ref,nb_next[0])) r= true;
      }
  }
  return r;
}

bool is_siloxy(int a_view, int a_ref) {
bool  r ;
int  nb_next_count, i;
int  nb_next[max_neighbors];
  r = false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->atype,"O3 ") && (atom[a_ref]->neighbor_count == 2)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          if (!strcmp(atom[(nb_next[0])]->element,"SI")) r = true;
      }
  }
  return r;
}

bool is_alkanoyl(int a_view, int a_ref) {
bool  r ;
int  nb_next_count, i;
int  nb_next[max_neighbors];
  r = false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if ((is_oxo_C(a_ref)) && (atom[a_ref]->neighbor_count == 3)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          for (i=0;i<2;i++) {
              if (is_alkyl(a_ref,nb_next[i])) r = true;
          }
      }    
  }
  return r;  
}

bool is_aroyl(int a_view, int a_ref) {
bool  r ;
int  nb_next_count, i;
int  nb_next[max_neighbors];
  r = false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if ((is_oxo_C(a_ref)) && (atom[a_ref]->neighbor_count == 3)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          for (i=0;i<2;i++){
              if (is_aryl(a_ref,nb_next[i])) r = true;
          }
      }    
   }
  return r;  
}

bool is_acyl(int a_view, int a_ref) {
bool r= false;
  if ((is_alkanoyl(a_view,a_ref)) || (is_aroyl(a_view,a_ref)))  r= true;
  return r;
}

bool is_acyl_gen(int a_view, int a_ref) {  // new in v0.3j
bool  r = false;
  if (is_oxo_C(a_ref)) r = true;
  return r;
}

bool is_acylamino(int a_view, int a_ref) {
bool  r ;
int  nb_next_count, i;
int  nb_next[max_neighbors];
int  acyl_count ;
  r = false;
  acyl_count = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->element,"N ") && (atom[a_ref]->neighbor_count == 2)) {
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          if (is_acyl(a_ref,nb_next[0])) acyl_count++;
          if (acyl_count == 1) r = true;  
      }
  }
return r;
}

bool is_carbamoyl(int a_view, int a_ref) {
bool  r ;
int  nb_next_count, i;
int  nb_next[max_neighbors];
  r = false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {

      if (is_oxo_C(a_ref) && (atom[a_ref]->neighbor_count == 3)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          for (i=0;i<2;i++){
              if (!strcmp(atom[(nb_next[i])]->atype,"N3 ") || 
                 !strcmp(atom[(nb_next[i])]->atype,"NAM")) r = true;
          }
      }    
  }
  return r;  
}

bool is_subst_acylamino(int a_view, int a_ref) {
// may be substituted _or_ unsubstituted acylamino group!
bool  r ;
int  nb_next_count, i;
int  nb_next[max_neighbors];
int  acyl_count ;
  r = false;
  acyl_count = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->element,"N ") && (atom[a_ref]->neighbor_count >= 2)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          for (i=0;i<(atom[a_ref]->neighbor_count - 1);i++) {
              if (is_acyl_gen(a_ref,nb_next[i])) acyl_count++;  // v0.3j
          }
          if (acyl_count > 0) r= true;  
      }
  }
return r;
}

bool is_subst_amino(int a_view, int a_ref) {
bool r= false;
  if ((is_amino(a_view,a_ref)) || (is_alkylamino(a_view,a_ref)) ||
     (is_arylamino(a_view,a_ref)) || (is_dialkylamino(a_view,a_ref)) ||
     (is_alkylarylamino(a_view,a_ref)) || (is_diarylamino(a_view,a_ref))) r =true;
return r;
}

void fgloc_set_hydrazino(int a_ref){   // v0.5
int  i, nb_count, k, a_n, b_id;
int  nb[max_neighbors];  // v0.3j
  if (opt_pos==false) return;
  a_n = 0;  // v0.5
  if (atom[a_ref]->heavy) {
      if (!strcmp(atom[a_ref]->element,"N ") && (atom[a_ref]->neighbor_count >= 2)) {
  memset(nb, 0, max_neighbors * sizeof(int));
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_ref) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
     nb[nb_count]=bond[k]->a2;
     nb_count++;
    }
    if ((bond[k]->a2== a_ref)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  } //NEIGHBOR copy the neighbor part here!! 
  for (i=0; i<atom[a_ref]->neighbor_count; i++) {
      if ((is_amino(a_ref,nb[i])) ||
                 (is_subst_amino(a_ref,nb[i]))) a_n = nb[i];
      }
      if (a_n > 0) {
          b_id = get_bond(a_ref,a_n);
          add2fgloc(fg_hydrazine,b_id);
      }
    }
  }
}

bool is_hydroximino_C(int id) {
int  i,k,nb_count;
bool  r;
int  nb[max_neighbors];
int  a_het ;
  r = false;
  a_het = 0;
  if ((id < 1) || (id > n_atoms)) return r;
  memset(nb, 0, max_neighbors * sizeof(int));
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == id) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
     nb[nb_count]=bond[k]->a2;
     nb_count++;
    }
    if ((bond[k]->a2== id)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  } //NEIGHBOR copy the neighbor part here!!
  if (!strcmp(atom[id]->element,"C ") && (atom[id]->neighbor_count > 0)) {
      for (i=0;i<atom[id]->neighbor_count; i++) {
          if ((bond[get_bond(id,nb[i])]->btype == 'D') &&
             !strcmp(atom[(nb[i])]->element,"N ") && (hetbond_count(nb[i]) == 3)) a_het = nb[i];
      }
      if (a_het > 0) {
  memset(nb, 0, max_neighbors * sizeof(int));
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_het) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
     nb[nb_count]=bond[k]->a2;
     nb_count++;
    }
    if ((bond[k]->a2== a_het)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  } //NEIGHBOR copy the neighbor part here!!
          if (!strcmp(atom[a_het]->element,"N ") && (atom[a_het]->neighbor_count > 0)) {
              for (i = 0;i<atom[a_het]->neighbor_count;i++) {
                  if (is_hydroxy(a_het,nb[i])) r = true;
              }
          }
      }
  }
return r;
}

bool is_hydrazino(int a_view, int a_ref) {
bool  r ;
int  nb_next_count, i, NR_count;
int  nb_next[max_neighbors];
  r = false;
  NR_count = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->element,"N ") && (atom[a_ref]->neighbor_count >= 2)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          for (i=0;i<(atom[a_ref]->neighbor_count - 1);i++) { // fixed in v0.3c
              if ((is_amino(a_ref,nb_next[i])) ||
                 (is_subst_amino(a_ref,nb_next[i]))) NR_count++;
          }
          if (NR_count == 1) r= true ; 
      }
  }
return r;
}

bool is_hydroxylamino(int a_view, int a_ref) {
bool  r ;
int  nb_next_count, i;
int  nb_next[max_neighbors];
int  OH_count, het_count;
  r = false;
  OH_count  = 0;
  het_count = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->element,"N ") && (atom[a_ref]->neighbor_count >= 2)) { // v0.3c
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
        for (i=0;i<(atom[a_ref]->neighbor_count - 1);i++) { // v0.3c
         if (is_hydroxy(a_ref,nb_next[i])) OH_count++;
         if (strcmp(atom[nb_next[i]]->element,"C ") && strcmp(atom[nb_next[i]]->element,"H ") &&   // v0.3k
           strcmp(atom[nb_next[i]]->element,"D ") && strcmp(atom[nb_next[i]]->element,"DU") &&
           strcmp(atom[nb_next[i]]->element,"LP")) het_count++; // v0.3n: D
         }
        if ((OH_count == 1) && (het_count == 1)) r= true;  
      }
  }
return r;
}

bool is_alkoxycarbonyl(int a_view, int a_ref) {
bool  r ;
int  nb_next_count, i;
int  nb_next[max_neighbors];
  r = false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if ((is_oxo_C(a_ref)) && (atom[a_ref]->neighbor_count == 3)) {
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          for (i=0;i<2;i++) {
              if (is_alkoxy(a_ref,nb_next[i])) r= true;
          }
      }    
  }
  return r;  
}

bool is_aryloxycarbonyl(int a_view, int a_ref) {
bool  r=false ;
int  nb_next_count, i;
int  nb_next[max_neighbors];
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if ((is_oxo_C(a_ref)) && (atom[a_ref]->neighbor_count == 3)) {
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          for (i=0;i<2;i++) {
              if (is_aryloxy(a_ref,nb_next[i])) r = true;
          }
      }    
  }
  return r;  
}

bool is_alkoxythiocarbonyl(int a_view, int a_ref) {
bool  r ;
int  nb_next_count, i;
int  nb_next[max_neighbors];
  r= false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if ((is_thioxo_C(a_ref)) && (atom[a_ref]->neighbor_count == 3)) {
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          for (i=0;i<2;i++) {
              if (is_alkoxy(a_ref,nb_next[i])) r = true;
          }
      }    
  }
  return r;  
}

bool is_aryloxythiocarbonyl(int a_view, int a_ref) {
bool  r ;
int  nb_next_count, i;
int  nb_next[max_neighbors];
  r = false;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if ((is_thioxo_C(a_ref)) && (atom[a_ref]->neighbor_count == 3)) {
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          for (i=0;i<2;i++) {
              if (is_aryloxy(a_ref,nb_next[i]))  r= true;
          }
      }    
  }
  return r;  
}

bool is_imino_C(int id){
bool r;
int  i, k, nb_count;
int  nb[max_neighbors];
  r = false;
  if ((id < 1) || (id > n_atoms)) return r;
    nb_count=0;
    for (k= 0;k< n_bonds;k++ ){
      if ((bond[k]->a1 == id) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
        nb[nb_count]=bond[k]->a2;
        nb_count++;
      }
      if ((bond[k]->a2== id)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
        nb[nb_count]=bond[k]->a1;
        nb_count++;
      }
    }
  if (!strcmp(atom[id]->element,"C ") && (atom[id]->neighbor_count > 0)) {
      for (i=0;i<atom[id]->neighbor_count;i++) {
          if ((bond[get_bond(id,nb[i])]->btype == 'D') && (!strcmp(atom[(nb[i])]->element,"N "))) {
             r= true;
          }
      }
  }
  return r;
}

bool is_hydrazono_C(int id) {
int  i,k,nb_count;
bool  r;
int  nb[max_neighbors];
int  a_het;
  r = false;
  a_het = 0;
  if ((id < 1) || (id > n_atoms)) return r;
  memset(nb, 0, max_neighbors * sizeof(int));
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == id) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
     nb[nb_count]=bond[k]->a2;
     nb_count++;
    }
    if ((bond[k]->a2== id)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  } //NEIGHBOR copy the neighbor part here!!
  if (!strcmp(atom[id]->element,"C ") && (atom[id]->neighbor_count > 0)) {
      for (i= 0;i<atom[id]->neighbor_count;i++) { 
          if ((bond[get_bond(id,nb[i])]->btype == 'D') &&
           !strcmp(atom[(nb[i])]->element,"N ")) a_het = nb[i];
      }
      if (a_het > 0) {
  memset(nb, 0, max_neighbors * sizeof(int));
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_het) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
     nb[nb_count]=bond[k]->a2;
     nb_count++;
    }
    if ((bond[k]->a2== a_het)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  } //NEIGHBOR copy the neighbor part here!!
          if (!strcmp(atom[a_het]->element,"N ") && (atom[a_het]->neighbor_count > 0)) {
              for (i=0;i<atom[a_het]->neighbor_count;i++) {
                  if (is_amino(a_het,nb[i]) || is_alkylamino(a_het,nb[i]) ||
                     is_alkylarylamino(a_het,nb[i]) || is_arylamino(a_het,nb[i]) ||
                     is_dialkylamino(a_het,nb[i]) || is_diarylamino(a_het,nb[i])) r = true;
              }
          }
      }
  }
  return r;
}

bool is_C_monosubst_amino(int a_view, int a_ref) {  // new in v0.3j, a_ref = N
bool  r ;
int  nb_next_count, i;
int  nb_next[max_neighbors];
int  C_count ;
  r = false;
  C_count = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if ((!strcmp(atom[a_ref]->atype,"N3 ") || !strcmp(atom[a_ref]->atype,"NAM")) && 
        (atom[a_ref]->neighbor_count == 2)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          if (!strcmp(atom[(nb_next[0])]->element,"C ")) C_count++;
          if (C_count == 1) r = true; 
      }
  }
  return r;
}

bool is_C_disubst_amino(int a_view, int a_ref) {  // new in v0.3j
bool  r ;
int  nb_next_count, i;
int  nb_next[max_neighbors];
int  b, C_count ;
  r = false;
  C_count = 0;
  b = get_bond(a_view,a_ref);
  if ((atom[a_view]->heavy) && (bond[b]->btype == 'S') && (bond[b]->arom == false)) {
      if ((!strcmp(atom[a_ref]->atype,"N3 ") || !strcmp(atom[a_ref]->atype,"NAM")) && 
        (atom[a_ref]->neighbor_count == 3)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          for (i =0;i< 2;i++) {
              if (!strcmp(atom[(nb_next[i])]->element,"C ")) C_count++;
          }
          if (C_count == 2) r = true;  
      }
  }
return r;
}

bool is_nitro(int a_view, int a_ref) {
bool  r ;
int  nb_next_count, i, bond_count;
int  nb_next[max_neighbors];
int  O_count ;
  r = false;
  O_count = 0;
  bond_count = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->element,"N ") && (atom[a_ref]->neighbor_count == 3)) {
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          for (i=0;i<2;i++) {
              if (!strcmp(atom[(nb_next[i])]->element,"O ")) O_count++;
              if (bond[get_bond(a_ref,nb_next[i])]->btype == 'S') bond_count++;
              if (bond[get_bond(a_ref,nb_next[i])]->btype == 'D') bond_count=bond_count+2;
          }
          if ((O_count == 2) && (bond_count >= 3))  r= true;  
      }
  }
  return r;
}

bool is_azido(int a_view, int a_ref) {
bool  r ;
int  nb_next_count, i, bond_count, n1, n2, n3;
int  nb_next[max_neighbors];
  r = false;
  bond_count = 0;
  n1 = 0; n2 = 0; n3 = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->element,"N ") && (atom[a_ref]->neighbor_count == 2)) { // v0.3c
          n1 = a_ref;
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 ==n1)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == n1)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          if (!strcmp(atom[(nb_next[0])]->element,"N ")) {
              n2 = nb_next[0];
              if (bond[get_bond(n1,n2)]->btype == 'S') bond_count++;
              if (bond[get_bond(n1,n2)]->btype == 'D') bond_count=bond_count+2;
              if (bond[get_bond(n1,n2)]->btype == 'T') bond_count=bond_count+3;
          }
          if ((n2 > 0) && (atom[n2]->neighbor_count == 2)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 ==n2)&& (bond[i]->a2 != n1) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == n2)&& (bond[i]->a1 != n1) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
              if (!strcmp(atom[(nb_next[0])]->element,"N ")) {
                  n3 = nb_next[0];
                  if (bond[get_bond(n2,n3)]->btype == 'S') bond_count++;
                  if (bond[get_bond(n2,n3)]->btype == 'D') bond_count=bond_count+2;
                  if (bond[get_bond(n2,n3)]->btype == 'T') bond_count=bond_count+3;
              }
          } 
          if ((n1 > 0) && (n2 > 0) && (n3 > 0) && (atom[n3]->neighbor_count == 1) &&
             (bond_count > 3)) r= true;  
      }
  }
return r;
}

bool is_nitroso(int a_view, int a_ref) { // new in v0.3j
bool  r ;
int  nb_next_count, i, O_count, a2;
int  nb_next[max_neighbors];
  r = false;
  O_count = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->element,"N ") && (atom[a_ref]->neighbor_count == 2)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 ==a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          a2 = nb_next[0];
          if (!strcmp(atom[a2]->element,"O ") && (bond[get_bond(a_ref,a2)]->btype == 'D')) O_count++;
          if (O_count == 1) r= true;  
      }
  }
  return r;
}

bool is_subst_hydrazino(int a_view, int a_ref) {
int  i, nb_next_count;
bool  r ;
int  nb_next[max_neighbors];
int  NR_count, a2;
  r = false;
  NR_count = 0;
  if ((atom[a_view]->heavy) && (bond[get_bond(a_view,a_ref)]->btype == 'S')) {
      if (!strcmp(atom[a_ref]->element,"N ") && (atom[a_ref]->neighbor_count >= 2)) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 ==a_ref)&& (bond[i]->a2 != a_view) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_ref)&& (bond[i]->a1 != a_view) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
          for (i=0;i<(atom[a_ref]->neighbor_count-1);i++) {
            a2 = nb_next[i];
            if (!strcmp(atom[a2]->element,"N ") && !(is_nitroso(a_ref,a2))) NR_count++;  // v0.3j
          }
          if (NR_count == 1) r= true;  
      }
  }
  return r;
}

void chk_ccx(int a_view, int a_ref) {
int  i, nb_count, k, b_id;
int  nb[max_neighbors];
int  OH_count, OR_count, N_count ;  // v0.5
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_ref) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_ref)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  OH_count = 0; OR_count = 0; N_count = 0;
  for (i=0;i<atom[a_ref]->neighbor_count;i++) {
      if (bond[(get_bond(a_ref,nb[i]))]->btype == 'S') {
          if (is_hydroxy(a_ref,nb[i])) OH_count++;
          if ((is_alkoxy(a_ref,nb[i])) || (is_aryloxy(a_ref,nb[i])) ||
             (is_siloxy(a_ref,nb[i]))) OR_count++;
          if (!strcmp(atom[(nb[i])]->atype,"N3 ") || !strcmp(atom[(nb[i])]->atype,"NAM")) N_count++;
      }
  }
  if (OH_count == 1) {
      fg[fg_enol]     = true;
      if (opt_pos) add2fgloc(fg_enol,a_ref);  // v0.5
  }
  if (OR_count == 1) {
      fg[fg_enolether] = true;  
      if (opt_pos) add2fgloc(fg_enolether,a_ref);  // v0.5
  }
  if (N_count == 1)  {
      fg[fg_enamine]  = true;
      if (opt_pos) add2fgloc(fg_enamine,a_ref);  // v0.5
  }
  // new in v0.2f   (regard anything else as an alkene)
  if ((OH_count + OR_count + N_count) == 0) {
      fg[fg_alkene] = true;
        if (opt_pos) {
            b_id = get_bond(a_view,a_ref);
            add2fgloc(fg_alkene,b_id);
        }
  }
}

void chk_xccx(int a_view,int a_ref) {
int  i, nb_count, k, b_id;
int  nb[max_neighbors];
int  OH_count, OR_count, N_count;  // v0.5
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_view) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_view)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  OH_count = 0; OR_count = 0; N_count = 0;
  for (i=0; i<atom[a_view]->neighbor_count;i++) {
      if (bond[(get_bond(a_view,nb[i]))]->btype == 'S') {
          if (is_hydroxy(a_view,nb[i])) OH_count++;
          if ((is_alkoxy(a_view,nb[i])) || (is_aryloxy(a_view,nb[i])) ||
             (is_siloxy(a_view,nb[i]))) OR_count++;
          if (!strcmp(atom[(nb[i])]->atype,"N3 ") ||
             !strcmp(atom[(nb[i])]->atype,"NAM")) N_count++;
      }
  }
  memset(nb, 0, max_neighbors* sizeof(int));
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_ref) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_ref)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  for (i=0;i<atom[a_ref]->neighbor_count;i++) {
      if (bond[(get_bond(a_ref,nb[i]))]->btype == 'S') {
          if (is_hydroxy(a_ref,nb[i])) OH_count++;
          if ((is_alkoxy(a_ref,nb[i])) ||
             (is_aryloxy(a_ref,nb[i])) ||
             (is_siloxy(a_ref,nb[i]))) OR_count++;
          if (!strcmp(atom[(nb[i])]->atype,"N3 ") ||
             !strcmp(atom[(nb[i])]->atype,"NAM")) N_count++;
      }
  }
  if (OH_count == 2) { 
      fg[fg_enediol]     = true;
      if (opt_pos) {  // v0.5
          b_id = get_bond(a_view,a_ref);
          add2fgloc(fg_enediol,b_id);
      }
  }
  // new in v0.2f   (regard anything else as an alkene)
  if ((OH_count + OR_count + N_count) == 0) { 
      fg[fg_alkene] = true;
      if (opt_pos) {
          b_id = get_bond(a_view,a_ref);
          add2fgloc(fg_alkene,b_id);
      }
  }
}

void chk_n_o_dbl(int a1, int a2) {
int  i, nb_count, k;
int  nb[max_neighbors];
int OR_count, N_count, C_count;
int b, het_count;          // v0.3j
char  bt ;            // v0.3k
float  bo_sum;      // v0.3k
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a1) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a1)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  OR_count = 0; N_count = 0; C_count = 0; bo_sum = 0.0;    // v0.3k
  het_count = 0; // v0.3k
  for (i=0; i<atom[a1]->neighbor_count;i++) {
      if (nb[i] != a2) { 
          b = get_bond(a1,nb[i]);         // v0.3j
          bt = bond[b]->btype;                           // v0.3k
          if (strcmp(atom[(nb[i])]->element,"C ") && strcmp(atom[(nb[i])]->element,"H ") &&
           strcmp(atom[(nb[i])]->element,"D ") && strcmp(atom[(nb[i])]->element,"DU") &&
           strcmp(atom[(nb[i])]->element,"LP") && (bond[b]->arom == false)) het_count++; 
// added 'D ' in v0.3n; // v0.3k: ignore hetero atoms in aromatic rings like isoxazole 
          if (bt == 'S') bo_sum = bo_sum + 1;           
          if (bt == 'D') bo_sum = bo_sum + 2;
          if (bt == 'A') bo_sum = bo_sum + 1.5;
          if (!strcmp(atom[(nb[i])]->element,"O ")) OR_count++;
          if (!strcmp(atom[(nb[i])]->element,"N ")) N_count++;
          if (!strcmp(atom[(nb[i])]->element,"C ") && (bond[b]->btype == 'S')) C_count++;  // v0.3k
          // if (is_alkyl(a1,nb[i])) or (is_aryl(a1,nb[i])) then inc(c_count);
      }
  }
  if (((OR_count+N_count+C_count)== 1) && (atom[a1]->neighbor_count == 2)) {  // excludes nitro etc.
      if (OR_count == 1) { 
          fg[fg_nitrite]       = true;
          if (opt_pos) add2fgloc(fg_nitrite,a1);   // v0.5
      }
      if (C_count == 1) { 
          fg[fg_nitroso_compound] = true;
          if (opt_pos) add2fgloc(fg_nitroso_compound,a1);   // v0.5
      }
      if (N_count == 1) {
          fg[fg_nitroso_compound]  = true; // instead of nitrosamine  v0.3j
          if (opt_pos) add2fgloc(fg_nitroso_compound,a1);   // v0.5
      }
      //if (n_count = 1) then fg[fg_nitrosamine]   := true;  // still missing
  }
  //if ((c_count > 1) and (or_count = 0) and (n_count = 0)) then
  //  begin
  //    fg[fg_n_oxide] := true;
  //  end;
  // new approach in v0.3k
  if ((het_count == 0) && (bo_sum > 2)) {  // =O does not count!
      fg[fg_n_oxide] = true;
      if (opt_pos) add2fgloc(fg_n_oxide,a1);   // v0.5
  }
}

void chk_c_hal(int a1,int a2) {
  fg[fg_halogen_deriv] = true;
  if (opt_pos) add2fgloc(fg_halogen_deriv,a2);  // v0.5
  if (atom[a1]->arom) {
    fg[fg_aryl_halide] = true;
    if (!strcmp(atom[a2]->element,"F ")) {
        fg[fg_aryl_fluoride] = true;
        if (opt_pos) add2fgloc(fg_aryl_fluoride,a2);  // v0.5
    }
    if (!strcmp(atom[a2]->element,"CL")) {
        fg[fg_aryl_chloride] = true;
        if (opt_pos) add2fgloc(fg_aryl_chloride,a2);  // v0.5
    }
    if (!strcmp(atom[a2]->element,"BR")) {
        fg[fg_aryl_bromide] = true;
        if (opt_pos) add2fgloc(fg_aryl_bromide,a2);  // v0.5
    }
    if (!strcmp(atom[a2]->element,"I ")) {
        fg[fg_aryl_iodide]   = true;
        if (opt_pos) add2fgloc(fg_aryl_iodide,a2);  // v0.5
    }
  } else {
    if (!strcmp(atom[a1]->atype,"C3 ") && (hetbond_count(a1) <= 2)) {// alkyl halides
        fg[fg_alkyl_halide] = true;
        if (!strcmp(atom[a2]->element,"F ")) {
            fg[fg_alkyl_fluoride] = true;
            if (opt_pos) add2fgloc(fg_alkyl_fluoride,a2);  // v0.5
        }
        if (!strcmp(atom[a2]->element,"CL")) {
            fg[fg_alkyl_chloride] = true;
            if (opt_pos) add2fgloc(fg_alkyl_chloride,a2);  // v0.5
        }
        if (!strcmp(atom[a2]->element,"BR")) {
            fg[fg_alkyl_bromide]  = true;
            if (opt_pos) add2fgloc(fg_alkyl_bromide,a2);  // v0.5
        }
        if (!strcmp(atom[a2]->element,"I ")) {
            fg[fg_alkyl_iodide]   = true;                      
            if (opt_pos) add2fgloc(fg_alkyl_iodide,a2);  // v0.5
        }
    }
    if (!strcmp(atom[a1]->atype,"C2 ") && (hetbond_count(a1) == 3)) { // acyl halides and related compounds
        if (is_oxo_C(a1)) {
            fg[fg_acyl_halide] = true;
            if (opt_pos) add2fgloc(fg_acyl_halide,a1);  // v0.5
            if (!strcmp(atom[a2]->element,"F ")) {
                fg[fg_acyl_fluoride] = true;
                if (opt_pos) add2fgloc(fg_acyl_fluoride,a1);  // v0.5
            }
            if (!strcmp(atom[a2]->element,"CL")) { 
                fg[fg_acyl_chloride] = true;
                if (opt_pos) add2fgloc(fg_acyl_chloride,a1);  // v0.5
            }
            if (!strcmp(atom[a2]->element,"BR")) {
                fg[fg_acyl_bromide]  = true;
                if (opt_pos) add2fgloc(fg_acyl_bromide,a1);  // v0.5
            }
            if (!strcmp(atom[a2]->element,"I ")) {
                fg[fg_acyl_iodide]   = true;                      
                if (opt_pos) add2fgloc(fg_acyl_iodide,a1);  // v0.5
            }
        }
        if (is_thioxo_C(a1)) {
            fg[fg_thiocarboxylic_acid_deriv] = true;
            if (opt_pos) add2fgloc(fg_thiocarboxylic_acid_deriv,a1);  // v0.5
        }
        if (is_imino_C(a1)) {
            fg[fg_imidoyl_halide] = true;
            if (opt_pos) add2fgloc(fg_imidoyl_halide,a1);  // v0.5
        }
    }
    if (!strcmp(atom[a1]->atype,"C2 ") && (hetbond_count(a1) == 4)) {// chloroformates etc.
        fg[fg_co2_deriv] = true;
        if (opt_pos) add2fgloc(fg_co2_deriv,a1);  // v0.5
        if (is_oxo_C(a1)) {
            fg[fg_carbonic_acid_deriv] = true;
            if (opt_pos) add2fgloc(fg_carbonic_acid_deriv,a1);  // v0.5
            if ((is_alkoxycarbonyl(a2,a1)) || (is_aryloxycarbonyl(a2,a1))) {
                fg[fg_carbonic_acid_ester_halide] = true;
                if (opt_pos) add2fgloc(fg_carbonic_acid_ester_halide,a1);  // v0.5
            }
            if (is_carbamoyl(a2,a1)) {
                fg[fg_carbamic_acid_deriv] = true;
                if (opt_pos) add2fgloc(fg_carbamic_acid_deriv,a1);  // v0.5
                fg[fg_carbamic_acid_halide] = true;
                if (opt_pos) add2fgloc(fg_carbamic_acid_halide,a1);  // v0.5
             }
        }
        if (is_thioxo_C(a1)) {
            fg[fg_thiocarbonic_acid_deriv] = true;
            if (opt_pos) add2fgloc(fg_thiocarbonic_acid_deriv,a1);  // v0.5
            if ((is_alkoxythiocarbonyl(a2,a1)) || (is_aryloxythiocarbonyl(a2,a1))) {
                fg[fg_thiocarbonic_acid_ester_halide] = true;
                if (opt_pos) add2fgloc(fg_thiocarbonic_acid_ester_halide,a1);  // v0.5
            }
            if (is_thiocarbamoyl(a2,a1)) {
                fg[fg_thiocarbamic_acid_deriv] = true;
                if (opt_pos) add2fgloc(fg_thiocarbamic_acid_deriv,a1);  // v0.5
                fg[fg_thiocarbamic_acid_halide] = true;
                if (opt_pos) add2fgloc(fg_thiocarbamic_acid_halide,a1);  // v0.5
            }
         }
     }
    // still missing: polyhalogen compounds (-CX2H, -CX3)
  }    // end of non-aromatic halogen compounds
}


void chk_c_o(int a1, int a2) {// a1 = C, a2 = O
  // ignore heteroaromatic rings (like furan, thiophene, etc.)
  if (bond[get_bond(a1,a2)]->arom == true) return;
  if (is_true_alkyl(a2,a1) && is_hydroxy(a1,a2)) {
      fg[fg_hydroxy] = true;
      fg[fg_alcohol] = true;
      if (opt_pos) {
          add2fgloc(fg_hydroxy,a2);
          add2fgloc(fg_alcohol,a2);
      }
      if (atom[a1]->neighbor_count <= 2) {
          fg[fg_prim_alcohol] = true;
          if (opt_pos) add2fgloc(fg_prim_alcohol,a2);   // v0.5
      }
      if (atom[a1]->neighbor_count == 3) {
          fg[fg_sec_alcohol]  = true;
          if (opt_pos) add2fgloc(fg_sec_alcohol,a2);   // v0.5
      }
      if (atom[a1]->neighbor_count == 4) {
          fg[fg_tert_alcohol] = true;
          if (opt_pos) add2fgloc(fg_tert_alcohol,a2);   // v0.5
      }
  }
  if (is_aryl(a2,a1) && is_hydroxy(a1,a2)) {
      fg[fg_hydroxy] = true;
      fg[fg_phenol]  = true;
      if (opt_pos) {
          add2fgloc(fg_hydroxy,a2);
          add2fgloc(fg_phenol,a2);
      }
  }
  if (is_true_alkyl(a2,a1) && is_true_alkoxy(a1,a2)) {
      fg[fg_ether]        = true;
      fg[fg_dialkylether] = true;
      if (opt_pos) {
          add2fgloc(fg_ether,a2);
          add2fgloc(fg_dialkylether,a2);
      }
  }
  if ((is_true_alkyl(a2,a1) && is_aryloxy(a1,a2)) || 
     (is_aryl(a2,a1) && is_true_alkoxy(a1,a2))) {
      fg[fg_ether]       = true;
      fg[fg_alkylarylether]= true;
      if (opt_pos) {
          add2fgloc(fg_ether,a2);
          add2fgloc(fg_alkylarylether,a2);
      }
  }
  if (is_aryl(a2,a1) && is_aryloxy(a1,a2)) {
      fg[fg_ether]        = true;
      fg[fg_diarylether] = true;
      if (opt_pos) {
          add2fgloc(fg_ether,a2);
          add2fgloc(fg_diarylether,a2);
      }
  }
  if ((is_true_alkyl(a2,a1) || is_aryl(a2,a1)) && is_alkynyloxy(a1,a2)) {
      fg[fg_ether]        = true;
      ether_generic       = true;
      if (opt_pos) add2fgloc(fg_ether,a2);   // v0.5
  }
  if (is_alkynyl(a2,a1) && is_hydroxy(a1,a2)) {
      fg[fg_hydroxy]  = true;
      hydroxy_generic = true;
      if (opt_pos) add2fgloc(fg_hydroxy,a2);   // v0.5
  }
}

void chk_c_s(int a1, int a2) {// a1 = C, a2 = S
int  i, nb_count, k;
int  nb[max_neighbors];
int  O_count, OH_count, OR_count, N_count, C_count, hal_count;
  // ignore heteroaromatic rings (like furan, thiophene, etc.)
  if (bond[get_bond(a1,a2)]->arom == true) return;
  if (is_alkyl(a2,a1) && is_sulfanyl(a1,a2)) {
      fg[fg_thiol] = true;
      fg[fg_alkylthiol] = true;
      if (opt_pos) {
          add2fgloc(fg_thiol,a2);
          add2fgloc(fg_alkylthiol,a2);
      }
  }
  if (is_aryl(a2,a1) && is_sulfanyl(a1,a2)) {
      fg[fg_thiol]       = true;
      fg[fg_arylthiol]   = true;
      if (opt_pos) {
          add2fgloc(fg_thiol,a2);
          add2fgloc(fg_arylthiol,a2);
      }
  }
  if (is_true_alkyl(a2,a1) && is_true_alkylsulfanyl(a1,a2)) {
      fg[fg_thioether]  = true;
      if (opt_pos) add2fgloc(fg_thioether,a2);   // v0.5
  }
  if ((is_true_alkyl(a2,a1) && is_arylsulfanyl(a1,a2)) || 
     (is_aryl(a2,a1) && is_true_alkylsulfanyl(a1,a2))) {
      fg[fg_thioether] = true;
      if (opt_pos) add2fgloc(fg_thioether,a2);   // v0.5
  }
  if (is_aryl(a2,a1) && is_arylsulfanyl(a1,a2)) {
      fg[fg_thioether]    = true;
      if (opt_pos) add2fgloc(fg_thioether,a2);   // v0.5
  }
  // check for sulfinic/sulfenic acid derivatives
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a2) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a2)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  O_count = 0; OH_count = 0; OR_count = 0; N_count = 0; C_count = 0; hal_count= 0;
  for (i =0;i<atom[a2]->neighbor_count;i++) {
      if (is_alkyl(a2,nb[i]) || is_aryl(a2,nb[i])) C_count++;
      if (is_hydroxy(a2,nb[i])) OH_count++;        
      if (is_alkoxy(a2,nb[i]) || is_aryloxy(a2,nb[i])) OR_count++;
      if (is_amino(a2,nb[i]) || is_subst_amino(a2,nb[i])) N_count++;
      if (!strcmp(atom[(nb[i])]->element,"F ") || !strcmp(atom[(nb[i])]->element,"CL") ||
        !strcmp(atom[(nb[i])]->element,"BR") || !strcmp(atom[(nb[i])]->element,"I ")) hal_count++;
      if (!strcmp(atom[(nb[i])]->element,"O ")) O_count++;
  }
  if (C_count == 1) {
      if ((atom[a2]->neighbor_count == 3) && ((O_count-(OH_count + OR_count)) == 1)) { // sulfinic acid & derivs
          fg[fg_sulfinic_acid_deriv]  = true;
          if (opt_pos) add2fgloc(fg_sulfinic_acid_deriv,a2);   // v0.5
          if (OH_count == 1)  {
              fg[fg_sulfinic_acid]       = true;
              if (opt_pos) add2fgloc(fg_sulfinic_acid,a2);   // v0.5
          }
          if (OR_count == 1)  {
              fg[fg_sulfinic_acid_ester]  = true;
              if (opt_pos) add2fgloc(fg_sulfinic_acid_ester,a2);   // v0.5
          }
          if (hal_count == 1) {
              fg[fg_sulfinic_acid_halide]= true;
              if (opt_pos) add2fgloc(fg_sulfinic_acid_halide,a2);   // v0.5
          }
          if (N_count == 1) {
              fg[fg_sulfinic_acid_amide]  = true;
              if (opt_pos) add2fgloc(fg_sulfinic_acid_amide,a2);   // v0.5
          }
      }
      if ((atom[a2]->neighbor_count == 2) && ((O_count-(OH_count+OR_count)) == 0)) { // sulfenic acid & derivs
          fg[fg_sulfenic_acid_deriv] = true;
          if (opt_pos) add2fgloc(fg_sulfenic_acid_deriv,a2);   // v0.5
          if (OH_count == 1) {
              fg[fg_sulfenic_acid]       = true;
              if (opt_pos) add2fgloc(fg_sulfenic_acid,a2);   // v0.5
          }
          if (OR_count == 1) {
              fg[fg_sulfenic_acid_ester]  = true;
              if (opt_pos) add2fgloc(fg_sulfenic_acid_ester,a2);   // v0.5
          }
          if (hal_count == 1) {
              fg[fg_sulfenic_acid_halide] = true;
              if (opt_pos) add2fgloc(fg_sulfenic_acid_halide,a2);   // v0.5
          }
          if (N_count == 1)  {
              fg[fg_sulfenic_acid_amide]  = true;
              if (opt_pos) add2fgloc(fg_sulfenic_acid_amide,a2);   // v0.5
          }
      }
  }
}

void chk_c_n(int a1, int a2) {// a1 = C, a2 = N
  // ignore heteroaromatic rings (like furan, thiophene, pyrrol, etc.)
  if (atom[a2]->arom == true) return;
  if (is_true_alkyl(a2,a1) && is_amino(a1,a2)) {
      fg[fg_amine]            = true;
      fg[fg_prim_amine]       = true;
      fg[fg_prim_aliph_amine] = true;
      if (opt_pos) {
          add2fgloc(fg_amine,a2);
          add2fgloc(fg_prim_amine,a2);
          add2fgloc(fg_prim_aliph_amine,a2);
      }
  }
  if (is_aryl(a2,a1) && is_amino(a1,a2)) {
      fg[fg_amine]            = true;
      fg[fg_prim_amine]       = true;
      fg[fg_prim_arom_amine]  = true;
      if (opt_pos) {
          add2fgloc(fg_amine,a2);
          add2fgloc(fg_prim_amine,a2);
          add2fgloc(fg_prim_arom_amine,a2); 
      }
  }
  if (is_true_alkyl(a2,a1) && is_true_alkylamino(a1,a2)) {
      fg[fg_amine]            = true;
      fg[fg_sec_amine]        = true;
      fg[fg_sec_aliph_amine]  = true;
      if (opt_pos) {
          add2fgloc(fg_amine,a2);
          add2fgloc(fg_sec_amine,a2);
          add2fgloc(fg_sec_aliph_amine,a2);
      }
  }
  if (is_aryl(a2,a1) && is_true_alkylamino(a1,a2)) {
      fg[fg_amine]            = true;
      fg[fg_sec_amine]        = true;
      fg[fg_sec_mixed_amine]  = true;
      if (opt_pos) {
          add2fgloc(fg_amine,a2);
          add2fgloc(fg_sec_amine,a2);
          add2fgloc(fg_sec_mixed_amine,a2);
      }
  }
  if (is_aryl(a2,a1) && is_arylamino(a1,a2)) {
      fg[fg_amine]            = true;
      fg[fg_sec_amine]        = true;
      fg[fg_sec_arom_amine]   = true;
      if (opt_pos) {
          add2fgloc(fg_amine,a2);
          add2fgloc(fg_sec_amine,a2);
          add2fgloc(fg_sec_arom_amine,a2);
      }
  }
  if (is_true_alkyl(a2,a1) && is_true_dialkylamino(a1,a2)) {
      fg[fg_amine]            = true;
      fg[fg_tert_amine]       = true;
      fg[fg_tert_aliph_amine] = true;
      if (opt_pos) {
          add2fgloc(fg_amine,a2);
          add2fgloc(fg_tert_amine,a2);
          add2fgloc(fg_tert_aliph_amine,a2);
      }
  }
  if ((is_true_alkyl(a2,a1) && is_diarylamino(a1,a2)) ||
     (is_aryl(a2,a1) && is_true_dialkylamino(a1,a2))){
      fg[fg_amine]            = true;
      fg[fg_tert_amine]       = true;
      fg[fg_tert_mixed_amine] = true;
      if (opt_pos) {
          add2fgloc(fg_amine,a2);
          add2fgloc(fg_tert_amine,a2);
          add2fgloc(fg_tert_mixed_amine,a2);
      }
  }
  if (is_aryl(a2,a1) && is_diarylamino(a1,a2)) {
      fg[fg_amine]            = true;
      fg[fg_tert_amine]       = true;
      fg[fg_tert_arom_amine]  = true;
      if (opt_pos) {
          add2fgloc(fg_amine,a2);
          add2fgloc(fg_tert_amine,a2);
          add2fgloc(fg_tert_arom_amine,a2);
      }
  }
  if ((is_alkyl(a2,a1) || is_aryl(a2,a1) ||    // v0.3k
      is_alkenyl(a2,a1) || is_alkynyl(a2,a1)) && (is_hydroxylamino(a1,a2) &&
      (is_acyl_gen(a2,a1)==false))) {
      fg[fg_hydroxylamine]      = true;        // v0.3k 
      if (opt_pos) add2fgloc(fg_hydroxylamine,a2);  // v0.5
  }
//printf("1:%d,%d,%d\n",is_nitro(a1,a2),is_alkyl(a2,a1),is_hydrazino(a1,a2));
//printf("4:%d,%d,%d,%d,%d,%d",a1,a2,is_nitro(a1,a2),is_alkyl(a2,a1),is_aryl(a2,a1),is_alkenyl(a2,a1));
  if ((is_alkyl(a2,a1) || is_aryl(a2,a1) || is_acyl(a2,a1) ||
      is_alkenyl(a2,a1) || is_alkynyl(a2,a1)) && (is_hydrazino(a1,a2))) {
      fg[fg_hydrazine]          = true;
      if (opt_pos) fgloc_set_hydrazino(a2);
  }
  if ((is_alkyl(a2,a1) || is_aryl(a2,a1) ||    // v0.3k
      is_alkenyl(a2,a1) || is_alkynyl(a2,a1)) && (is_azido(a1,a2))) {
      fg[fg_azide]           = true;
      if (opt_pos) add2fgloc(fg_azide,a2);   // v0.5
  }
  if ((is_alkyl(a2,a1) || is_aryl(a2,a1) ||    // v0.3k
      is_alkenyl(a2,a1) || is_alkynyl(a2,a1)) && (is_diazonium(a1,a2))) {
      fg[fg_diazonium_salt]           = true;
      if (opt_pos) add2fgloc(fg_diazonium_salt,a2);   // v0.5
  }
  if (((is_alkyl(a2,a1) || is_aryl(a2,a1) ||    // v0.3k
      is_alkenyl(a2,a1) || is_alkynyl(a2,a1))) && (is_nitro(a1,a2))) {
      fg[fg_nitro_compound]     = true;
      if (opt_pos) add2fgloc(fg_nitro_compound,a2);   // v0.5
  }
  if (is_alkynyl(a2,a1) && (is_amino(a1,a2) || is_C_monosubst_amino(a1,a2) || 
    is_C_disubst_amino(a1,a2)) && (!is_acylamino(a1,a2))) {  // v0.4c: fixed parentheses
      fg[fg_amine]            = true;
      amine_generic           = true;    
      if (opt_pos) add2fgloc(fg_amine,a2);   // v0.5
  }
}

void chk_c_c(int a1, int a2) {
int  i, nb_count, k, b_id;
int  nb[max_neighbors];
int  OH_count, NHR_count ;
bool  a1oh ;   // v0.5
  // ignore aromatic rings
  if (atom[a2]->arom == true) return;
  //check for 1,2-diols and 1,2-aminoalcoholes
  if (!strcmp(atom[a1]->atype,"C3 ") && !strcmp(atom[a2]->atype,"C3 ")) {
      if ((hetbond_count(a1) == 1) && (hetbond_count(a2) ==1)) {
          OH_count = 0; NHR_count = 0; a1oh = false;
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a1) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a1)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
          for (i=0;i<atom[a1]->neighbor_count;i++) {
              if (nb[i] != a2){
                  if (is_hydroxy(a1,nb[i])) {
                      OH_count++;        
                      a1oh = true;   // v0.5
                  }
                  if ((is_amino(a1,nb[i])) || (is_alkylamino(a1,nb[i])) 
                      || (is_arylamino(a1,nb[i])))  NHR_count++;
              }
          }  
  memset(nb,0,max_neighbors*sizeof(int));
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a2) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a2)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
          for (i=0;i<atom[a2]->neighbor_count;i++) { 
              if (nb[i] != a1) {
                  if (is_hydroxy(a2,nb[i])) OH_count++;        
                  if (is_amino(a2,nb[i]) || is_alkylamino(a2,nb[i]) ||
                    is_arylamino(a2,nb[i])) NHR_count++;
              }
          }
          if (OH_count == 2) { 
              fg[fg_1_2_diol] = true;
              if (opt_pos) {   // v0.5
                  b_id = get_bond(a1,a2);
                  add2fgloc(fg_1_2_diol,b_id);
              }
          }
          if ((OH_count == 1) && (NHR_count == 1)) {
              fg[fg_1_2_aminoalcohol] = true;
              if (opt_pos) {
                  if (a1oh) { add2fgloc(fg_1_2_aminoalcohol,a1); 
                  } else { add2fgloc(fg_1_2_aminoalcohol,a2);}
              }
          }
      }
  }
  // check for alpha-aminoacids and alpha-hydroxyacids
  if (!strcmp(atom[a1]->atype,"C3 ") && !strcmp(atom[a2]->atype,"C2 ")) {
      if ((hetbond_count(a1) == 1) && (hetbond_count(a2) == 3)) {
          OH_count = 0; NHR_count = 0;
  memset(nb,0,max_neighbors*sizeof(int));
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a1) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a1)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
          for (i=0;i<atom[a1]->neighbor_count;i++) {
              if (nb[i] != a2) {
                  if (is_hydroxy(a1,nb[i])) OH_count++;        
                  if (is_amino(a1,nb[i]) || is_alkylamino(a1,nb[i]) 
                      || is_arylamino(a1,nb[i])) NHR_count++;
              }
          }  
  memset(nb,0,max_neighbors*sizeof(int));
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a2) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a2)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
          for (i =0;i<atom[a2]->neighbor_count;i++) {
              if (nb[i] != a1) {
                  if (is_hydroxy(a2,nb[i])) OH_count++;        
              }
          }  
          if ((OH_count == 2) && (is_oxo_C(a2))) {
              fg[fg_alpha_hydroxyacid] = true;
              if (opt_pos) add2fgloc(fg_alpha_hydroxyacid,a2);   // v0.5
          }
          if ((OH_count == 1) && (NHR_count == 1) && (is_oxo_C(a2))) {
              fg[fg_alpha_aminoacid] = true;
              if (opt_pos) add2fgloc(fg_alpha_aminoacid,a2);   // v0.5
          }
      }
  }    
}

void chk_sulfoxide(int a1, int a2) {
int  i, nb_count, k;
int  nb[max_neighbors];
int  O_count, C_count;
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a1) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a1)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  O_count = 0; C_count = 0;
  for (i= 0;i<atom[a1]->neighbor_count;i++) {
      if (!strcmp(atom[(nb[i])]->element,"O ")) O_count++;
      if (is_alkyl(a1,nb[i]) || is_aryl(a1,nb[i]) || 
         is_alkenyl(a1,nb[i]) || is_alkynyl(a1,nb[i])) C_count++;   // v0.5
  }
  if ((O_count == 1) && (C_count == 2)) { 
      fg[fg_sulfoxide] = true;
      if (opt_pos) add2fgloc(fg_sulfoxide,a1);   // v0.5
  }
}

void chk_x_y_single(int a_view, int a_ref) {
int  b_id;  // v0.5
  if (!strcmp(atom[a_view]->atype,"O3 ") && !strcmp(atom[a_ref]->atype,"O3 ")) {
      if (is_hydroxy(a_ref,a_view) || is_hydroxy(a_view,a_ref)) {
          fg[fg_hydroperoxide] = true;
          if (opt_pos) {
              if (is_hydroxy(a_view,a_ref)) {add2fgloc(fg_hydroperoxide,a_view);
              } else {
                add2fgloc(fg_hydroperoxide,a_ref);}
          }
      }
      if ((is_alkoxy(a_ref,a_view) || is_aryloxy(a_ref,a_view) || is_siloxy(a_ref,a_view)) &&
         (is_alkoxy(a_view,a_ref) || is_aryloxy(a_view,a_ref) || is_siloxy(a_view,a_ref))) {
          fg[fg_peroxide] = true;
          if (opt_pos) {
              b_id = get_bond(a_view,a_ref);
              add2fgloc(fg_peroxide,b_id);
          }
      }
  }  // still missing: peracid
  if (!strcmp(atom[a_view]->atype,"S3 ") && !strcmp(atom[a_ref]->atype,"S3 ")) {
      if ((atom[a_view]->neighbor_count == 2) && (atom[a_ref]->neighbor_count== 2)) {
          fg[fg_disulfide] = true;
          if (opt_pos) {
              b_id = get_bond(a_view,a_ref);
              add2fgloc(fg_disulfide,b_id);
          }
      }
  }
  if (!strcmp(atom[a_view]->element,"N ") && !strcmp(atom[a_ref]->element,"N ") &&
     (hetbond_count(a_view) == 1) && (hetbond_count(a_ref) == 1)) {
      //if ((is_amino(a_ref,a_view)) or 
      //    (is_subst_amino(a_ref,a_view)) or
      //    (is_acylamino(a_ref,a_view))) and
      //   ((is_amino(a_view,a_ref)) or 
      //    (is_subst_amino(a_view,a_ref)) or
      //    (is_acylamino(a_ref,a_view))) then 
      if (bond[get_bond(a_view,a_ref)]->arom == false) {
          fg[fg_hydrazine] = true;
          if (opt_pos) {
              b_id = get_bond(a_view,a_ref);
              add2fgloc(fg_hydrazine,b_id);
          }
      }
  } 
  if (!strcmp(atom[a_view]->element,"N ") && !strcmp(atom[a_ref]->atype,"O3 ")) {
  // bond is in "opposite" direction
      if ((is_alkoxy(a_view,a_ref) || is_aryloxy(a_view,a_ref)) &&
         is_nitro(a_ref,a_view)) {
          fg[fg_nitrate] = true;
          if (opt_pos) add2fgloc(fg_nitrate,a_view);   // v0.5
      }
      if (((is_nitro(a_ref,a_view)==false) && (atom[a_view]->arom==false)) &&
         (is_amino(a_ref,a_view) || is_subst_amino(a_ref,a_view)) &&
         (is_acylamino(a_ref,a_view)==false)) {
          fg[fg_hydroxylamine] = true;    // new in v0.3c
          if (opt_pos) add2fgloc(fg_hydroxylamine,a_view);   // v0.5
      }
  }
  if (!strcmp(atom[a_view]->element,"S ") && !strcmp(atom[a_ref]->element,"O ")) {  
    chk_sulfoxide(a_view,a_ref);
  }
}

void chk_carboxyl_deriv(int a_view, int a_ref) {
int  i,k, nb_count;
int   O_count, N_count, S_count;
int  a_o=0, a_n=0, a_s=0;
int  nb[max_neighbors];
  O_count = 0; N_count = 0; S_count = 0;
  nb_count=0;
  memset(nb,0,max_neighbors*sizeof(int));
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_view) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_view)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  for (i=0; i<atom[a_view]->neighbor_count; i++) { 
      if (bond[(get_bond(a_view,nb[i]))]->btype == 'S') {
          if (strcmp(atom[(nb[i])]->element,"C ")) { 
              if (!strcmp(atom[(nb[i])]->element,"O ")) {
                O_count++; a_o = nb[i]; }
              if (!strcmp(atom[(nb[i])]->element,"N ")) {
                N_count++; a_n = nb[i]; }
              if (!strcmp(atom[(nb[i])]->element,"S ")) {
                S_count++; a_s = nb[i]; }
          }
      }
  }    
  if (is_oxo_C(a_view)) {
      if (O_count == 1) {  // anhydride is checked somewhere else
          if (bond[get_bond(a_view,a_o)]->arom == false) { 
              fg[fg_carboxylic_acid_deriv]  = true;
              if (opt_pos) add2fgloc(fg_carboxylic_acid_deriv,a_view);   // v0.5
          }
          if (is_hydroxy(a_view,a_o)) { 
              if (atom[a_o]->formal_charge == 0) { 
                  fg[fg_carboxylic_acid] = true;
                  if (opt_pos) add2fgloc(fg_carboxylic_acid,a_view);   // v0.5
              }
              if (atom[a_o]->formal_charge == -1) { 
                  fg[fg_carboxylic_acid_salt] = true;            
                  if (opt_pos) add2fgloc(fg_carboxylic_acid_salt,a_view);   // v0.5
              }
          }
          if (is_alkoxy(a_view,a_o) || is_aryloxy(a_view,a_o) ||
            is_alkenyloxy(a_view,a_o) || is_alkynyloxy(a_view,a_o)) { 
              if (bond[get_bond(a_view,a_o)]->arom == false) { 
                  fg[fg_carboxylic_acid_ester] = true;
                  if (opt_pos) add2fgloc(fg_carboxylic_acid_ester,a_view);   // v0.5
              }
              if (bond[get_bond(a_view,a_o)]->ring_count > 0) {
                  if (bond[get_bond(a_view,a_o)]->arom == true) {
                      fg[fg_oxohetarene] = true;
                      if (opt_pos) add2fgloc(fg_oxohetarene,a_view);   // v0.5
                  } else {
                      fg[fg_lactone] = true;
                      if (opt_pos) add2fgloc(fg_lactone,a_view);   // v0.5
                  }
              }
          }
      }
      if (N_count == 1) { 
          if (bond[get_bond(a_view,a_n)]->arom == false) {
              fg[fg_carboxylic_acid_deriv] = true;
              if (opt_pos) add2fgloc(fg_carboxylic_acid_deriv,a_view);   // v0.5
          } else { //fg[fg_lactam_heteroarom] := true;  // catches also pyridazines, 1,2,3-triazines, etc.
              fg[fg_oxohetarene] = true;
              if (opt_pos) add2fgloc(fg_oxohetarene,a_view);   // v0.5
          }
          if ((is_amino(a_view,a_n)) || 
            (!strcmp(atom[a_n]->atype,"NAM") && (atom[a_n]->neighbor_count == 1))) {
              fg[fg_carboxylic_acid_amide] = true;
              if (opt_pos) add2fgloc(fg_carboxylic_acid_amide,a_view);   // v0.5
              fg[fg_carboxylic_acid_prim_amide] = true;
              if (opt_pos) add2fgloc(fg_carboxylic_acid_prim_amide,a_view);   // v0.5
          }
          //if (is_alkylamino(a_view,a_n)) or (is_arylamino(a_view,a_n)) then 
          if ((is_C_monosubst_amino(a_view,a_n)) && (!(is_subst_acylamino(a_view,a_n)))) { // v0.3j
              if (bond[get_bond(a_view,a_n)]->arom == false) {
                  fg[fg_carboxylic_acid_amide] = true;
                  if (opt_pos) add2fgloc(fg_carboxylic_acid_amide,a_view);   // v0.5
              }
              if (bond[get_bond(a_view,a_n)]->arom == false) { 
                  fg[fg_carboxylic_acid_sec_amide] = true;
                  if (opt_pos) add2fgloc(fg_carboxylic_acid_sec_amide,a_view);   // v0.5
              }
              if (bond[get_bond(a_view,a_n)]->ring_count > 0) {
                  if (bond[get_bond(a_view,a_n)]->arom == true) {
                     fg[fg_oxohetarene] = true;
                     if (opt_pos) add2fgloc(fg_oxohetarene,a_view);   // v0.5
                  } else {
                        fg[fg_lactam] = true;
                        if (opt_pos) add2fgloc(fg_lactam,a_view);   // v0.5
                  }
              }
          }          
          //if (is_dialkylamino(a_view,a_n)) or (is_alkylarylamino(a_view,a_n)) or
          //   (is_diarylamino(a_view,a_n)) then 
          if ((is_C_disubst_amino(a_view,a_n)) && !(is_subst_acylamino(a_view,a_n))) { // v0.3j
              if (bond[get_bond(a_view,a_n)]->arom == false) {
                  fg[fg_carboxylic_acid_amide] = true;
                  if (opt_pos) add2fgloc(fg_carboxylic_acid_amide,a_view);   // v0.5
              }
              if (bond[get_bond(a_view,a_n)]->arom == false) { 
                  fg[fg_carboxylic_acid_tert_amide] = true;
                  if (opt_pos) add2fgloc(fg_carboxylic_acid_tert_amide,a_view);   // v0.5
              }
              if (bond[get_bond(a_view,a_n)]->ring_count > 0) {
                  if (bond[get_bond(a_view,a_n)]->arom == true) {
                    //fg[fg_lactam_heteroarom]    := true else 
                      fg[fg_oxohetarene]   = true;
                      if (opt_pos) add2fgloc(fg_oxohetarene,a_view);   // v0.5
                  } else {
                        fg[fg_lactam]  = true;
                        if (opt_pos) add2fgloc(fg_lactam,a_view);   // v0.5
                  }
              }
          }
          if (is_hydroxylamino(a_view,a_n)) { 
              fg[fg_hydroxamic_acid]  = true;
              if (opt_pos) add2fgloc(fg_hydroxamic_acid,a_view);   // v0.5
          }
          if (is_hydrazino(a_view,a_n)) { 
              fg[fg_carboxylic_acid_hydrazide] = true;
              if (opt_pos) add2fgloc(fg_carboxylic_acid_hydrazide,a_view);   // v0.5
          }
          if (is_azido(a_view,a_n)) { 
              fg[fg_carboxylic_acid_azide] = true;
              if (opt_pos) add2fgloc(fg_carboxylic_acid_azide,a_view);   // v0.5
          }
      }
      if (S_count == 1) {  // anhydride is checked somewhere else
          if (bond[get_bond(a_view,a_s)]->arom == false) {
              fg[fg_thiocarboxylic_acid_deriv]  = true;
              if (opt_pos) add2fgloc(fg_thiocarboxylic_acid_deriv,a_view);   // v0.5
          }
          if (is_sulfanyl(a_view,a_s)) { 
              fg[fg_thiocarboxylic_acid] = true;
              if (opt_pos) add2fgloc(fg_thiocarboxylic_acid,a_view);   // v0.5
          }
          if ((is_alkylsulfanyl(a_view,a_s)) || (is_arylsulfanyl(a_view,a_s))) { 
              if (bond[get_bond(a_view,a_s)]->arom == false) { 
                  fg[fg_thiocarboxylic_acid_ester] = true;
                  if (opt_pos) add2fgloc(fg_thiocarboxylic_acid_ester,a_view);   // v0.5
              }
              if (bond[get_bond(a_view,a_s)]->ring_count > 0) {
                  if (bond[get_bond(a_view,a_s)]->arom == true) {            
                      fg[fg_oxohetarene] = true;
                      if (opt_pos) add2fgloc(fg_oxohetarene,a_view);   // v0.5
                  } else {
                      fg[fg_thiolactone] = true;
                      if (opt_pos) add2fgloc(fg_thiolactone,a_view);   // v0.5
                  }
              }
          }            
      }
  }  // end oxo-C
  if (is_thioxo_C(a_view)) {
      if (O_count == 1) {  // anhydride is checked somewhere else
          if (bond[get_bond(a_view,a_o)]->arom == false) { 
              fg[fg_thiocarboxylic_acid_deriv] = true;
              if (opt_pos) add2fgloc(fg_thiocarboxylic_acid_deriv,a_view);   // v0.5
          }
          if (is_hydroxy(a_view,a_o)) { 
              fg[fg_thiocarboxylic_acid] = true;   // fixed in v0.3c
              if (opt_pos) add2fgloc(fg_thiocarboxylic_acid,a_view);   // v0.5
          }
          if ((is_alkoxy(a_view,a_o)) || (is_aryloxy(a_view,a_o))) { 
              if ((a_s>0) && (bond[get_bond(a_view,a_s)]->arom == false)) { 
                  fg[fg_thiocarboxylic_acid_ester] = true;
                  if (opt_pos) add2fgloc(fg_thiocarboxylic_acid_ester,a_view);   // v0.5
              }
              if (bond[get_bond(a_view,a_o)]->ring_count > 0) {
                  if (bond[get_bond(a_view,a_o)]->arom == true) {
                      fg[fg_thioxohetarene] = true;
                      if (opt_pos) add2fgloc(fg_thioxohetarene,a_view);   // v0.5
                  } else {
                      fg[fg_thiolactone]= true;
                      if (opt_pos) add2fgloc(fg_thiolactone,a_view);   // v0.5
                  }
              }
          }            
      }
      if (N_count == 1) { 
          if (bond[get_bond(a_view,a_n)]->arom == false) {
              fg[fg_thiocarboxylic_acid_deriv] = true;
              if (opt_pos) add2fgloc(fg_thiocarboxylic_acid_deriv,a_view);   // v0.5
          } else {
                fg[fg_thioxohetarene] = true;  // catches also pyridazines, 1,2,3-triazines, etc.
                if (opt_pos) add2fgloc(fg_thioxohetarene,a_view);   // v0.5
          }
          if ((is_amino(a_view,a_n)) || (atom[a_n]->neighbor_count == 1)) {  // v0.4a
              fg[fg_thiocarboxylic_acid_amide] = true;
              if (opt_pos) add2fgloc(fg_thiocarboxylic_acid_amide,a_view);   // v0.5
              // fg[fg_thiocarboxylic_acid_prim_amide] := true;
          }
          //if (is_alkylamino(a_view,a_n)) or (is_arylamino(a_view,a_n)) then 
          if ((is_C_monosubst_amino(a_view,a_n)) && !(is_subst_acylamino(a_view,a_n))) {// v0.3j
              if (bond[get_bond(a_view,a_n)]->arom == false) {
                  fg[fg_thiocarboxylic_acid_amide] = true;
                  if (opt_pos) add2fgloc(fg_thiocarboxylic_acid_amide,a_view);   // v0.5
              }
              //fg[fg_thiocarboxylic_acid_sec_amide]  := true;
              if (bond[get_bond(a_view,a_n)]->ring_count > 0) {
                  if (bond[get_bond(a_view,a_n)]->arom == true) {
                      fg[fg_thioxohetarene] = true;
                      if (opt_pos) add2fgloc(fg_thioxohetarene,a_view);   // v0.5
                  } else {
                      fg[fg_thiolactam] = true;
                      if (opt_pos) add2fgloc(fg_thiolactam,a_view);   // v0.5
                  }
              }
          }          
          //if (is_dialkylamino(a_view,a_n)) or (is_alkylarylamino(a_view,a_n)) or
          //   (is_diarylamino(a_view,a_n)) then 
          if ((is_C_disubst_amino(a_view,a_n)) && !(is_subst_acylamino(a_view,a_n))) {   // v0.3j
              if (bond[get_bond(a_view,a_n)]->arom == false) {
                  fg[fg_thiocarboxylic_acid_amide] = true;
                  if (opt_pos) add2fgloc(fg_thiocarboxylic_acid_amide,a_view);   // v0.5
              }
              //fg[fg_thiocarboxylic_acid_tert_amide] := true;
              if (bond[get_bond(a_view,a_n)]->ring_count > 0) {
                  if (bond[get_bond(a_view,a_n)]->arom == true) {
                    //fg[fg_thiolactam_heteroarom] := true else fg[fg_thiolactam] := true;
                      fg[fg_thioxohetarene] = true;
                      if (opt_pos) add2fgloc(fg_thioxohetarene,a_view);   // v0.5
                  } else {
                      fg[fg_thiolactam] = true;
                      if (opt_pos) add2fgloc(fg_thiolactam,a_view);   // v0.5
                  }
              }
          }          
      }
      if (S_count == 1) { // anhydride is checked somewhere else
          if ((a_s>0) && (bond[get_bond(a_view,a_s)]->arom == false)) { 
              fg[fg_thiocarboxylic_acid_deriv] = true;
              if (opt_pos) add2fgloc(fg_thiocarboxylic_acid_deriv,a_view);   // v0.5
          }
          if (is_sulfanyl(a_view,a_s)) {
              fg[fg_thiocarboxylic_acid] = true;
              if (opt_pos) add2fgloc(fg_thiocarboxylic_acid,a_view);   // v0.5
          }
          if ((is_alkylsulfanyl(a_view,a_s)) || (is_arylsulfanyl(a_view,a_s))) { 
              if (bond[get_bond(a_view,a_s)]->arom == false) {
                  fg[fg_thiocarboxylic_acid_ester] = true;
                  if (opt_pos) add2fgloc(fg_thiocarboxylic_acid_ester,a_view);   // v0.5
              }
              if (bond[get_bond(a_view,a_s)]->ring_count > 0) {
                  if (bond[get_bond(a_view,a_s)]->arom == true) {
                      fg[fg_thioxohetarene] = true;
                      if (opt_pos) add2fgloc(fg_thioxohetarene,a_view);   // v0.5
                  } else {
                      fg[fg_thiolactone] = true;
                      if (opt_pos) add2fgloc(fg_thiolactone,a_view);   // v0.5
                  }
              }
          }            
      }
  } // end thioxo-C
  if (is_true_imino_C(a_view)) {
      if (O_count == 1) { 
          if (bond[get_bond(a_view,a_o)]->arom == false) { 
              fg[fg_carboxylic_acid_deriv] = true;
              if (opt_pos) add2fgloc(fg_carboxylic_acid_deriv,a_view);   // v0.5
          }
          if ((is_alkoxy(a_view,a_o)) || (is_aryloxy(a_view,a_o))) {
              if (bond[get_bond(a_view,a_o)]->arom == false) {
                  fg[fg_imido_ester] = true;
                  if (opt_pos) add2fgloc(fg_imido_ester,a_view);   // v0.5
              }
          }            
      }
      if ((N_count == 1) && (bond[get_bond(a_view,a_n)]->arom == false)) { 
          if (bond[get_bond(a_view,a_n)]->arom == false) { 
              fg[fg_carboxylic_acid_deriv] = true;
              if (opt_pos) add2fgloc(fg_carboxylic_acid_deriv,a_view);   // v0.5
          }
          if ((is_amino(a_view,a_n)) || (is_subst_amino(a_view,a_n))) {
              if (bond[get_bond(a_view,a_n)]->arom == false) { 
                  fg[fg_carboxylic_acid_deriv] = true;
                  if (opt_pos) add2fgloc(fg_carboxylic_acid_deriv,a_view);   // v0.5
                  fg[fg_carboxylic_acid_amidine]  = true;
                  if (opt_pos) add2fgloc(fg_carboxylic_acid_amidine,a_view);   // v0.5
              }
          }          
          if (is_hydrazino(a_view,a_n)) { 
              if (bond[get_bond(a_view,a_n)]->arom == false) { 
                  fg[fg_carboxylic_acid_amidrazone] = true;
                  if (opt_pos) {
                      add2fgloc(fg_carboxylic_acid_amidrazone,a_view);   // v0.5
                      fgloc_set_hydrazino(a_n);   // v0.5
                  }
              }
          }
      }
      if ((N_count == 1) && (bond[get_bond(a_view,a_n)]->arom == true)) { // catches also pyridazines, 1,2,3-triazines, etc
          fg[fg_iminohetarene] = true;
          if (opt_pos) add2fgloc(fg_iminohetarene,a_view);   // v0.5
      }
      if (S_count == 1) { 
          if (bond[get_bond(a_view,a_s)]->arom == false)  { 
              fg[fg_carboxylic_acid_deriv]  = true;
              if (opt_pos) add2fgloc(fg_carboxylic_acid_deriv,a_view);   // v0.5
          }
          if ((is_alkylsulfanyl(a_view,a_s)) || (is_arylsulfanyl(a_view,a_s))) { 
              if (bond[get_bond(a_view,a_s)]->arom == false) {
                  fg[fg_imido_thioester] = true;
                  if (opt_pos) add2fgloc(fg_imido_thioester,a_view);   // v0.5
              }
          }            
      }
  }   // is true_imino_c
  if (is_hydroximino_C(a_view)) {
      if ((bond[get_bond(a_view,a_n)]->arom == false) && (a_n>0)) {
          fg[fg_carboxylic_acid_deriv]  = true;
          if (opt_pos) add2fgloc(fg_carboxylic_acid_deriv,a_view);   // v0.5
      } 
      if (O_count == 1) { 
          if (is_hydroxy(a_view,a_o)) { 
              fg[fg_hydroxamic_acid] = true;
              if (opt_pos) add2fgloc(fg_hydroxamic_acid,a_view);   // v0.5
          }          
      }
  }   // is hydroximino_c
  if (is_hydrazono_C(a_view)) {
      if ((bond[get_bond(a_view,a_n)]->arom == false) && (a_n>0)) { 
          fg[fg_carboxylic_acid_deriv]  = true;
          if (opt_pos) add2fgloc(fg_carboxylic_acid_deriv,a_view);   // v0.5
      } 
      if (N_count == 1) { 
          if ((is_amino(a_view,a_n)) || (is_subst_amino(a_view,a_n))) { 
              fg[fg_carboxylic_acid_amidrazone] = true;
              if (opt_pos) add2fgloc(fg_carboxylic_acid_amidrazone,a_view);   // v0.5
          }          
      }
  }   // is hydrazono_c
}

void chk_co2_sp2(int a_view, int a_ref) {
int  i, k, nb_count;
int  nb[max_neighbors];
int  O_count, OR_count, N_count, NN_count, NNX_count, S_count, SR_count;
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_view) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_view)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  O_count = 0; OR_count = 0; N_count = 0;
  NN_count = 0; NNX_count = 0; S_count = 0; SR_count = 0;
  for (i=0;i<atom[a_view]->neighbor_count;i++) {
      if (bond[(get_bond(a_view,nb[i]))]->btype == 'S') {
          if (strcmp(atom[(nb[i])]->element,"C ")) {
              if (!strcmp(atom[(nb[i])]->element,"O ")) {
                  O_count++;
                  if ((is_alkoxy(a_view,nb[i])) || (is_alkenyloxy(a_view,nb[i])) ||  // v0.3j
                     (is_aryloxy(a_view,nb[i]))) OR_count++;
              }
              if (!strcmp(atom[(nb[i])]->element,"N ")) { 
                  N_count++; 
                  if (is_hydrazino(a_view,nb[i])) NN_count++;
                  if (is_subst_hydrazino(a_view,nb[i])) NNX_count++;  // more general...
              }
              if (!strcmp(atom[(nb[i])]->element,"S ")) {
                  S_count++; 
                  if ((is_alkylsulfanyl(a_view,nb[i])) ||
                     (is_arylsulfanyl(a_view,nb[i]))) SR_count++;
              }
          }
      }
  }     
  if (is_oxo_C(a_view)) {
      if (O_count == 2) { 
          fg[fg_carbonic_acid_deriv] = true;
          if (opt_pos) add2fgloc(fg_carbonic_acid_deriv,a_view);   // v0.5
          if (OR_count == 1) { 
              fg[fg_carbonic_acid_monoester] = true;
              if (opt_pos) add2fgloc(fg_carbonic_acid_monoester,a_view);   // v0.5
          }
          if (OR_count == 2) {
              fg[fg_carbonic_acid_diester] = true;
              if (opt_pos) add2fgloc(fg_carbonic_acid_diester,a_view);   // v0.5
          }
      }
      if ((O_count == 1) && (S_count == 1)) {   
          fg[fg_thiocarbonic_acid_deriv] = true;
          if (opt_pos) add2fgloc(fg_thiocarbonic_acid_deriv,a_view);   // v0.5
          if (OR_count +SR_count == 1) {
              fg[fg_thiocarbonic_acid_monoester] = true;
              if (opt_pos) add2fgloc(fg_thiocarbonic_acid_monoester,a_view);   // v0.5
          }
          if (OR_count + SR_count == 2) { 
              fg[fg_thiocarbonic_acid_diester] = true;
              if (opt_pos) add2fgloc(fg_thiocarbonic_acid_diester,a_view);   // v0.5
          }
      }
      if (S_count == 2) { 
          fg[fg_thiocarbonic_acid_deriv] = true;
          if (opt_pos) add2fgloc(fg_thiocarbonic_acid_deriv,a_view);   // v0.5
          if (SR_count == 1) {
              fg[fg_thiocarbonic_acid_monoester] = true;
              if (opt_pos) add2fgloc(fg_thiocarbonic_acid_monoester,a_view);   // v0.5
          }
          if (SR_count == 2) { 
              fg[fg_thiocarbonic_acid_diester] = true;
              if (opt_pos) add2fgloc(fg_thiocarbonic_acid_diester,a_view);   // v0.5
          }
      }
      if ((O_count == 1) && (N_count == 1)) {  
          fg[fg_carbamic_acid_deriv] = true;
          if (opt_pos) add2fgloc(fg_carbamic_acid_deriv,a_view);   // v0.5
          if (OR_count == 0) { 
              fg[fg_carbamic_acid] = true;
              if (opt_pos) add2fgloc(fg_carbamic_acid,a_view);   // v0.5
          }
          if (OR_count == 1) { 
              fg[fg_carbamic_acid_ester] = true;
              if (opt_pos) add2fgloc(fg_carbamic_acid_ester,a_view);   // v0.5
          }
      }
      if ((S_count == 1) && (N_count == 1)) { 
          fg[fg_thiocarbamic_acid_deriv] = true;
          if (opt_pos) add2fgloc(fg_thiocarbamic_acid_deriv,a_view);   // v0.5
          if (SR_count == 0) {
              fg[fg_thiocarbamic_acid] = true;
              if (opt_pos) add2fgloc(fg_thiocarbamic_acid,a_view);   // v0.5
          }
          if (SR_count == 1) { 
              fg[fg_thiocarbamic_acid_ester] = true;
              if (opt_pos) add2fgloc(fg_thiocarbamic_acid_ester,a_view);   // v0.5
          }
      }
      if (N_count == 2) {
          if (NN_count == 1) {
              fg[fg_semicarbazide]  = true;
              if (opt_pos) add2fgloc(fg_semicarbazide,a_view);   // v0.5
          } else {
              if (NNX_count == 0) { 
                  fg[fg_urea] = true;  // excludes semicarbazones
                  if (opt_pos) add2fgloc(fg_urea,a_view);   // v0.5
              }
          }                     
      }  
  }  // end oxo-C
  if (is_thioxo_C(a_view)) {
      if (O_count == 2) {   
          fg[fg_thiocarbonic_acid_deriv] = true;
          if (opt_pos) add2fgloc(fg_thiocarbonic_acid_deriv,a_view);   // v0.5
          if (OR_count == 1) { 
              fg[fg_thiocarbonic_acid_monoester] = true;
              if (opt_pos) add2fgloc(fg_thiocarbonic_acid_monoester,a_view);   // v0.5
          }
          if (OR_count == 2) { 
              fg[fg_thiocarbonic_acid_diester] = true;
              if (opt_pos) add2fgloc(fg_thiocarbonic_acid_diester,a_view);   // v0.5
          }
      }
      if ((O_count == 1) && (S_count == 1)) {
          fg[fg_thiocarbonic_acid_deriv]  = true;
          if (opt_pos) add2fgloc(fg_thiocarbonic_acid_deriv,a_view);   // v0.5
          if (OR_count + SR_count == 1) { 
              fg[fg_thiocarbonic_acid_monoester]  = true;
              if (opt_pos) add2fgloc(fg_thiocarbonic_acid_monoester,a_view);   // v0.5
          }
          if (OR_count + SR_count == 2) { 
              fg[fg_thiocarbonic_acid_diester] = true;
              if (opt_pos) add2fgloc(fg_thiocarbonic_acid_diester,a_view);   // v0.5
          }
      }
      if (S_count == 2) { 
          fg[fg_thiocarbonic_acid_deriv]  = true;
          if (opt_pos) add2fgloc(fg_thiocarbonic_acid_deriv,a_view);   // v0.5
          if (SR_count == 1) {
              fg[fg_thiocarbonic_acid_monoester] = true;
              if (opt_pos) add2fgloc(fg_thiocarbonic_acid_monoester,a_view);   // v0.5
          }
          if (SR_count == 2) { 

              fg[fg_thiocarbonic_acid_diester] = true;
              if (opt_pos) add2fgloc(fg_thiocarbonic_acid_diester,a_view);   // v0.5
          }
      }
      if ((O_count == 1) && (N_count == 1)) {   
          fg[fg_thiocarbamic_acid_deriv] = true;
          if (opt_pos) add2fgloc(fg_thiocarbamic_acid_deriv,a_view);   // v0.5
          if (OR_count == 0) {
              fg[fg_thiocarbamic_acid] = true;
              if (opt_pos) add2fgloc(fg_thiocarbamic_acid,a_view);   // v0.5
          }
          if (OR_count == 1) { 
              fg[fg_thiocarbamic_acid_ester] = true;
              if (opt_pos) add2fgloc(fg_thiocarbamic_acid_ester,a_view);   // v0.5
          }
      }
      if ((S_count == 1) && (N_count == 1)) { 
          fg[fg_thiocarbamic_acid_deriv]  = true;
          if (opt_pos) add2fgloc(fg_thiocarbamic_acid_deriv,a_view);   // v0.5
          if (SR_count == 0) {
              fg[fg_thiocarbamic_acid]  = true;
              if (opt_pos) add2fgloc(fg_thiocarbamic_acid,a_view);   // v0.5
          }
          if (SR_count == 1) {
              fg[fg_thiocarbamic_acid_ester] = true;
              if (opt_pos) add2fgloc(fg_thiocarbamic_acid_ester,a_view);   // v0.5
          }
      }
      if (N_count == 2) {
          if (NN_count == 1) {
              fg[fg_thiosemicarbazide] = true;
              if (opt_pos) add2fgloc(fg_thiosemicarbazide,a_view);   // v0.5
          } else {
              if (NNX_count == 0) {
                  fg[fg_thiourea] = true;  // excludes thiosemicarbazones
                  if (opt_pos) add2fgloc(fg_thiourea,a_view);   // v0.5
              }
          }                    
      } 
  }  // end thioxo-C
  if ((is_true_imino_C(a_view)) && (bond[get_bond(a_view,a_ref)]->arom == false)) {
      if ((O_count == 1) && (N_count == 1)) {
          fg[fg_isourea]  = true;
          if (opt_pos) add2fgloc(fg_isourea,a_view);   // v0.5
      }
      if ((S_count == 1) && (N_count == 1)) {
          fg[fg_isothiourea]  = true;
          if (opt_pos) add2fgloc(fg_isothiourea,a_view);   // v0.5
      }
      if (N_count == 2) {
          fg[fg_guanidine] = true;
          if (opt_pos) add2fgloc(fg_guanidine,a_view);   // v0.5
      }  
  }  // end Imino-C
}

void chk_co2_sp(int a_view, int a_ref) {
int  i, k, nb_count;
int  nb[max_neighbors];
int  O_count, N_count, S_count ;
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_view) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_view)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  O_count = 0; N_count = 0; S_count = 0;
  for (i =0; i<atom[a_view]->neighbor_count;i++) {
      if (bond[(get_bond(a_view,nb[i]))]->btype == 'D') {
          if (strcmp(atom[(nb[i])]->element,"C ")) { 
              if (!strcmp(atom[(nb[i])]->element,"O ")) O_count++; 
              if (!strcmp(atom[(nb[i])]->element,"N ")) N_count++; 
              if (!strcmp(atom[(nb[i])]->element,"S ")) S_count++; 
          }
      }
  }     
  if (O_count + S_count == 2) {
      fg[fg_co2_deriv] = true;  // new in v0.3b
      if (opt_pos) add2fgloc(fg_co2_deriv,a_view);   // v0.5
  }
  if ((O_count == 1) && (N_count == 1)) { 
      fg[fg_isocyanate] = true;
      if (opt_pos) add2fgloc(fg_isocyanate,a_view);   // v0.5
  }
  if ((S_count == 1) && (N_count == 1)) {
      fg[fg_isothiocyanate] = true;
      if (opt_pos) add2fgloc(fg_isothiocyanate,a_view);   // v0.5
  }
  if (N_count == 2) {
      fg[fg_carbodiimide] = true;
      if (opt_pos) add2fgloc(fg_carbodiimide,a_view);   // v0.5
  }
}

void chk_ion(int a_ref) {
int  i, charge, k, nb_count;
int  nb[max_neighbors];
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_ref) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_ref)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  charge = atom[a_ref]->formal_charge;
  if (charge != 0) {  // check if charge is neutralized by an adjacent opposite charge
      for (i =0; i<atom[a_ref]->neighbor_count;i++) {
          charge = charge + atom[(nb[i])]->formal_charge;
      }
      if (charge > 0) {
          fg[fg_cation] = true;
          if (opt_pos) add2fgloc(fg_cation,a_ref);  // v0.5
      }
      if (charge < 0) { 
          fg[fg_anion]  = true;
          if (opt_pos) add2fgloc(fg_anion,a_ref);  // v0.5
      }
  }
}

bool is_cyano(int a_view, int a_ref) {// a_view = C, a_ref = N
bool  r= false;
  if (!strcmp(atom[a_view]->atype,"C1 ") && (bond[get_bond(a_view,a_ref)]->btype == 'T') &&
     !strcmp(atom[a_ref]->atype,"N1 ") && (atom[a_ref]->neighbor_count == 1)) r = true;
  return r;
}

bool is_cyano_c(int a_ref){
bool r;
int  i, k, nb_count;
int  nb[max_neighbors];
  r = false;
  if (!strcmp(atom[a_ref]->atype,"C1 ") && (atom[a_ref]->neighbor_count > 0)) {
    nb_count=0;
    for (k= 0;k< n_bonds;k++ ){
      if ((bond[k]->a1 == a_ref) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
        nb[nb_count]=bond[k]->a2;
        nb_count++;
      }
      if ((bond[k]->a2== a_ref)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
        nb[nb_count]=bond[k]->a1;
        nb_count++;
      }
    }
    for (i=0;i<atom[a_ref]->neighbor_count;i++) {
        if (is_cyano(a_ref,nb[i])) r = true;
    }
  }
return r;
}

void chk_imine(int a_ref, int a_view) { // a_ref = C, a_view = N
int  i, k, nb_count;
int  nb[max_neighbors];
int  a_het, a_c, het_count, C_count, O_count;  // v0.3k
  het_count = 0; C_count = 0; O_count = 0;  // v0.3k
  if (atom[a_view]->neighbor_count == 1) {
      if (atom[a_ref]->arom == false) { 
          fg[fg_imine] = true;
          if (opt_pos) add2fgloc(fg_imine,a_ref);  // v0.5
      }
  } else {
    nb_count=0;
    for (k= 0;k< n_bonds;k++ ) {
      if ((bond[k]->a1 == a_view) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
        nb[nb_count]=bond[k]->a2;
        nb_count++;
      }
      if ((bond[k]->a2== a_view)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
        nb[nb_count]=bond[k]->a1;
        nb_count++;
      }
    }
      if (atom[a_view]->neighbor_count > 1) {
          for (i =0; i<atom[a_view]->neighbor_count;i++) {
              if ((nb[i] != a_ref) && (bond[(get_bond(a_view,nb[i]))]->btype == 'S')) {
                if (!strcmp(atom[(nb[i])]->element,"C ")) { 
                      a_c = nb[i];
                      C_count++;
                  }
                if (!strcmp(atom[(nb[i])]->element,"O ") || !strcmp(atom[(nb[i])]->element,"N ")){ 
                      a_het = nb[i];
                      het_count++;
                }
                if (!strcmp(atom[(nb[i])]->element,"O ") && (atom[(nb[i])]->neighbor_count == 1) &&
                 (bond[(get_bond(a_view,nb[i]))]->arom == false)) O_count++;
              } // v0.3k
              if ((nb[i] != a_ref) && (bond[(get_bond(a_view,nb[i]))]->btype == 'D' )) {
   // v0.3k; make sure we do not count nitro groups in "azi" form etc.
                if (!strcmp(atom[(nb[i])]->element,"O ") || !strcmp(atom[(nb[i])]->element,"N ") ||
                !strcmp(atom[(nb[i])]->element,"S ")) { 
                      a_het = nb[i];  // v0.3m
                      het_count++;
                }
                if (!strcmp(atom[(nb[i])]->element,"O ") && (atom[(nb[i])]->neighbor_count == 1) &&
                 (bond[(get_bond(a_view,nb[i]))]->arom == false)) O_count++;// v0.3k
              }
          }
          if (C_count == 1) { 
              if (((is_alkyl(a_view,a_c)) || (is_aryl(a_view,a_c)) ||
                  (is_alkenyl(a_view,a_c)) || (is_alkynyl(a_view,a_c))) &&
                 (atom[a_ref]->arom == false) && (het_count == 0)) { 
                     fg[fg_imine] = true;   // v0.3k
                     if (opt_pos) add2fgloc(fg_imine,a_ref);  // v0.5
              }
          }
          if (het_count == 1) {
              if (!strcmp(atom[a_het]->element,"O ")) {
                  if (is_hydroxy(a_view,a_het)) {
                      fg[fg_oxime] = true;
                      if (opt_pos) add2fgloc(fg_oxime,a_ref);  // v0.5
                  }
                  if ((is_alkoxy(a_view,a_het)) || (is_aryloxy(a_view,a_het)) ||
                     (is_alkenyloxy(a_view,a_het)) || (is_alkynyloxy(a_view,a_het))) {
                      fg[fg_oxime_ether] = true;                  
                      if (opt_pos) add2fgloc(fg_oxime_ether,a_ref);  // v0.5
                  }
              }
              if (!strcmp(atom[a_het]->element,"N ")) {
                  if ((is_amino(a_view,a_het)) || (is_alkylamino(a_view,a_het)) ||
                     (is_dialkylamino(a_view,a_het)) || (is_alkylarylamino(a_view,a_het)) ||
                     (is_arylamino(a_view,a_het)) || (is_diarylamino(a_view,a_het))) { 
                      fg[fg_hydrazone] = true;
                      if (opt_pos) add2fgloc(fg_hydrazone,a_ref);  // v0.5
                  } else { // check for semicarbazone or thiosemicarbazone
                  memset(nb, 0, max_neighbors*sizeof(int));
                  nb_count=0;
                  for (k= 0;k< n_bonds;k++ ) {
                  if ((bond[k]->a1 == a_het) && (nb_count < max_neighbors) &&
                    (atom[bond[k]->a2]->heavy)) {
                    nb[nb_count]=bond[k]->a2;
                    nb_count++;
                  }
                  if ((bond[k]->a2== a_het)&& (nb_count < max_neighbors)&&
                  (atom[bond[k]->a1]->heavy)) {
                    nb[nb_count]=bond[k]->a1;
                    nb_count++;
                  }
                  }
                    if (atom[a_het]->neighbor_count > 1) {
                      for (i=0;i<atom[a_het]->neighbor_count;i++){ 
                        if (nb[i] != a_view) {
                          if (is_carbamoyl(a_het,nb[i])) { 
                              fg[fg_semicarbazone] = true;
                               if (opt_pos) add2fgloc(fg_semicarbazone,a_ref);  // v0.5
                          }
                         if (is_thiocarbamoyl(a_het,nb[i])) {
                              fg[fg_thiosemicarbazone] = true;                                
                             if (opt_pos) add2fgloc(fg_thiosemicarbazone,a_ref);  // v0.5
                         }
                       }
                     }
                   } 
                 }
              }
          }     // v0.3k: nitro groups in "azi" form
          if ((het_count == 2) && (O_count == 2)) { 
              fg[fg_nitro_compound] = true;
              if (opt_pos) add2fgloc(fg_nitro_compound,a_view);  // v0.5
          }
      }
  }
}

bool is_nitrile(int a_view, int a_ref) { // a_view = C, a_ref = N
bool r;
int nb_next[max_neighbors];
int nb_next_count, i;
  r = false;
  if (is_cyano(a_view,a_ref)) {
    if ((atom[a_view]->neighbor_count == 1) && (atom[a_view]->formal_charge == 0)) { 
      r = true; 
    } else {  // HCN is also a nitrile!
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_view)&& (bond[i]->a2 != a_ref) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_view)&& (bond[i]->a1 != a_ref) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
      if (!strcmp(atom[(nb_next[0])]->element,"C ") || !strcmp(atom[(nb_next[0])]->element,"H ") ||
        !strcmp(atom[(nb_next[0])]->element,"D ")) r = true;  // v0.3n: D
    }
  }
  return r;
}

void chk_carbonyl_deriv(int a_view, int a_ref) { // a_view = C
int  i, k;
int  nb[max_neighbors];
int  C_count, CN_count;
char  bt;     // new in v0.3b
int  n_db, nb_count;  // new in v0.3b
    nb_count=0;
    for (k= 0;k< n_bonds;k++ ){
      if ((bond[k]->a1 == a_view) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
        nb[nb_count]=bond[k]->a2;
        nb_count++;
      }
      if ((bond[k]->a2== a_view)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
        nb[nb_count]=bond[k]->a1;
        nb_count++;
      }
    }
  C_count = 0; CN_count = 0;
  n_db = 0;  // new in v0.3b
  for (i = 0;i<atom[a_view]->neighbor_count; i++) {
      bt = bond[(get_bond(a_view,nb[i]))]->btype;
      if (bt == 'S') {
          if (!strcmp(atom[(nb[i])]->element,"C ")) { 
              if (is_cyano_c(nb[i])) {
                CN_count++;
              } else {
                C_count++;
              } 
          }
      } else {  // new in v0.3b
          if (bt == 'D') n_db++;
      }
  }
  if (is_oxo_C(a_view)) {
      fg[fg_carbonyl]  = true;
      if (opt_pos) add2fgloc(fg_carbonyl,a_view);  // v0.5
      if ((C_count + CN_count) < 2) {  // new in v0.3b (detection of ketenes)
          if (n_db <= 1) { 
              fg[fg_aldehyde] = true;
              if (opt_pos) add2fgloc(fg_aldehyde,a_view);  // v0.5
          } else {
              fg[fg_ketene] = true;
              if (opt_pos) add2fgloc(fg_ketene,a_view);  // v0.5
          }
      }
      if (C_count == 2) { 
          if (atom[a_view]->arom) {
              fg[fg_oxohetarene] = true;
              if (opt_pos) add2fgloc(fg_oxohetarene,a_view);  // v0.5
          } else {
              fg[fg_ketone]   = true;
              if (opt_pos) add2fgloc(fg_ketone,a_view);  // v0.5
          }
      }
      if (CN_count > 0) { 
          fg[fg_acyl_cyanide] = true;
          if (opt_pos) add2fgloc(fg_acyl_cyanide,a_view);  // v0.5
      }
  }
  if (is_thioxo_C(a_view)) {
      fg[fg_thiocarbonyl]  = true;
      if (opt_pos) add2fgloc(fg_thiocarbonyl,a_view);  // v0.5
      if (C_count < 2) {
          fg[fg_thioaldehyde]  = true;
          if (opt_pos) add2fgloc(fg_thioaldehyde,a_view);  // v0.5
      }
      if (C_count == 2) {
          if (atom[a_view]->arom) {
              fg[fg_thioxohetarene] = true;
              if (opt_pos) add2fgloc(fg_thioxohetarene,a_view);  // v0.5
          } else {
              fg[fg_thioketone]   = true;
              if (opt_pos) add2fgloc(fg_thioketone,a_view);  // v0.5
          }
      }
  }
  if (is_imino_C(a_view)) {
      chk_imine(a_view,a_ref);
  }
}

bool is_isonitrile(int a_view,int a_ref) {   // only recognized with CN triple bond!
// a_view = C, a_ref = N
bool  r= false;
  if (!strcmp(atom[a_view]->atype,"C1 ") && (bond[get_bond(a_view,a_ref)]->btype == 'T') &&
     !strcmp(atom[a_ref]->atype,"N1 ") && (atom[a_ref]->neighbor_count == 2) &&
     (atom[a_view]->neighbor_count == 1))  r = true;
  return r;
}

bool is_cyanate(int a_view, int a_ref) { // a_view = C, a_ref = N
bool  r = false;
int nb_next[max_neighbors];
int nb_next_count, i;
  if (is_cyano(a_view,a_ref)) {
    if (atom[a_view]->neighbor_count == 2) {
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_view)&& (bond[i]->a2 != a_ref) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_view)&& (bond[i]->a1 != a_ref) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
        if ((is_alkoxy(a_view,nb_next[0])) ||
           (is_aryloxy(a_view,nb_next[0]))) r = true;
    }
  }
  return r;
}

bool is_thiocyanate(int a_view, int a_ref) {
bool  r = false;
int nb_next[max_neighbors];
int nb_next_count, i;
  if (is_cyano(a_view,a_ref)) {
    if (atom[a_view]->neighbor_count == 2) { 
      memset(nb_next, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == a_view)&& (bond[i]->a2 != a_ref) && (nb_next_count < max_neighbors) 
        && (atom[bond[i]->a2]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == a_view)&& (bond[i]->a1 != a_ref) && (nb_next_count < max_neighbors)  
        && (atom[bond[i]->a1]->heavy)) {
          nb_next[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
        if ((is_alkylsulfanyl(a_view,nb_next[0])) || 
           (is_arylsulfanyl(a_view,nb_next[0]))) r = true;
    }
  }
  return r;
}

bool is_true_exocyclic_imino_C(int id, int r_id) {  // v0.3j
int  r    ;
int  i, j, nb_count,k;
int  nb[max_neighbors];
RINGPATH_TYPE  testring ;
int  b, ring_size ;
  r = false;
  if ((id < 1) || (id > n_atoms)) return r;
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == id) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== id)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  memset(testring,0,sizeof(RINGPATH_TYPE));
  ring_size = ringprop[r_id]->size;  // v0.3j
  for (j=0;j<ring_size;j++) testring[j] = *(*(ring+r_id)+j);  // v0.3j
  if (!strcmp(atom[id]->element,"C ") && (atom[id]->neighbor_count > 0)) {
      for (i=0;i<atom[id]->neighbor_count;i++) {
          b= get_bond(id,nb[i]);
          if ((bond[b]->btype == 'D') && (bond[b]->arom == false) &&
             !strcmp(atom[(nb[i])]->element,"N ")) {
                 r = true;
                 for (j =0;j<ring_size;j++) {
                   if (nb[i] == *(*(ring+r_id)+j)) r= false; 
                 }
          }
      }
  }
  return r;
}

void chk_triple(int a1, int a2) {
int  b_id;  // v0.5
  if (!strcmp(atom[a1]->element,"C ") && !strcmp(atom[a2]->element,"C ") &&
    (bond[get_bond(a1,a2)]->arom == false)) {
      fg[fg_alkyne] = true;
      if (opt_pos) {
          b_id = get_bond(a1,a2);
          add2fgloc(fg_alkyne,b_id);
      }
  }
  if (is_nitrile(a1,a2)) { 
      fg[fg_nitrile]     = true;
      if (opt_pos) add2fgloc(fg_nitrile,a1);
  }
  if (is_isonitrile(a1,a2)) { 
      fg[fg_isonitrile]  = true;
      if (opt_pos) add2fgloc(fg_isonitrile,a2);  // a2 is the N atom!
  }
  if (is_cyanate(a1,a2)) { 
      fg[fg_cyanate]     = true;
      if (opt_pos) add2fgloc(fg_cyanate,a1);
  }
  if (is_thiocyanate(a1,a2)) { 
      fg[fg_thiocyanate] = true;
      if (opt_pos) add2fgloc(fg_thiocyanate,a1);
  }
}

void chk_double(int a1, int a2) {
int  b_id;  // v0.5
  if (!strcmp(atom[a1]->element,"C ") && strcmp(atom[a2]->element,"C ") &&
    (bond[get_bond(a1,a2)]->arom == false)) {
      if (hetbond_count(a1) == 2) chk_carbonyl_deriv(a1,a2);
      if (hetbond_count(a1) == 3) chk_carboxyl_deriv(a1,a2);
      if (hetbond_count(a1) == 4) {
          if (!strcmp(atom[a1]->atype,"C2 ")) chk_co2_sp2(a1,a2);
          if (!strcmp(atom[a1]->atype,"C1 ")) chk_co2_sp(a1,a2);
      }
  }  // end C=X
  if (!strcmp(atom[a1]->atype,"C2 ") && !strcmp(atom[a2]->atype,"C2 ")
    && (bond[get_bond(a1,a2)]->arom == false)) {
      if ((hetbond_count(a1) == 0) && (hetbond_count(a2) == 2)) {
          fg[fg_ketene_acetal_deriv] = true;
          if (opt_pos) add2fgloc(fg_ketene_acetal_deriv,a2);  // v0.5
      }
      if ((hetbond_count(a1) == 0) && (hetbond_count(a2) == 1)) chk_ccx(a1,a2);
      if ((hetbond_count(a1) == 1) && (hetbond_count(a2) == 1)) chk_xccx(a1,a2);
      if (((hetbond_count(a1) == 0) && (hetbond_count(a2) == 0)) &&
         (atom[a1]->arom == false) && (atom[a2]->arom == false)) { 
          fg[fg_alkene] = true;
          if (opt_pos) {   // v0.5
              b_id = get_bond(a1,a2);
              add2fgloc(fg_alkene,b_id);
          }
      }
  }
  if (!strcmp(atom[a1]->element,"N ") && !strcmp(atom[a2]->element,"N ") && 
     (hetbond_count(a1) == 2) && (hetbond_count(a2) == 2) &&
     (bond[get_bond(a1,a2)]->arom == false) &&
     (atom[a1]->neighbor_count == 2) && (atom[a2]->neighbor_count == 2)) {
           fg[fg_azo_compound] = true;
           if (opt_pos) {  // v0.5
               b_id = get_bond(a1,a2);
               add2fgloc(fg_azo_compound,b_id);
           }
      }
  if (!strcmp(atom[a1]->element,"N ") && !strcmp(atom[a2]->element,"O ")) chk_n_o_dbl(a1,a2);
  if (!strcmp(atom[a1]->element,"S ") && !strcmp(atom[a2]->element,"O ")) chk_sulfoxide(a1,a2);
}

void chk_single(int a1, int a2) {
  if (!strcmp(atom[a1]->element,"C ") && 
    (!strcmp(atom[a2]->element,"F ") || !strcmp(atom[a2]->element,"CL") || 
    !strcmp(atom[a2]->element,"BR") || !strcmp(atom[a2]->element,"I "))) chk_c_hal(a1,a2);
  if (!strcmp(atom[a1]->element,"C ") && !strcmp(atom[a2]->element,"O ")) chk_c_o(a1,a2);   
  if (!strcmp(atom[a1]->element,"C ") && !strcmp(atom[a2]->element,"S ")) chk_c_s(a1,a2);   
  if (!strcmp(atom[a1]->element,"C ") && !strcmp(atom[a2]->element,"N ")) chk_c_n(a1,a2);       
  if (!strcmp(atom[a1]->element,"C ") && (atom[a2]->metal && (is_cyano_c(a1) == false))) {
      fg[fg_organometallic] = true;
      if (opt_pos) add2fgloc(fg_organometallic,a2);   // v0.5
      if (!strcmp(atom[a2]->element,"LI")) { 
          fg[fg_organolithium] = true;
          if (opt_pos) add2fgloc(fg_organolithium,a2);   // v0.5
      }
      if (!strcmp(atom[a2]->element,"MG")) { 
           fg[fg_organomagnesium] = true;
           if (opt_pos) add2fgloc(fg_organomagnesium,a2);   // v0.5
      }
  } 
  if (!strcmp(atom[a1]->element,"C ") && !strcmp(atom[a2]->element,"C ")) chk_c_c(a1,a2);       
  if (strcmp(atom[a1]->element,"C ") && strcmp(atom[a2]->element,"C ")) chk_x_y_single(a1,a2);
}

void chk_carbonyl_deriv_sp3(int a_ref) {
int  i, nb_count,k;
int  nb[max_neighbors];
int  OH_count, OR_count, N_count, SH_count, SR_count;
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_ref) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_ref)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  OH_count = 0; OR_count = 0; N_count = 0; SH_count = 0; SR_count = 0;
  for (i=0;i<atom[a_ref]->neighbor_count; i++) {
      if (is_hydroxy(a_ref,nb[i])) OH_count++;
      if ((is_alkoxy(a_ref,nb[i])) || (is_aryloxy(a_ref,nb[i])) ||
         (is_alkenyloxy(a_ref,nb[i])) || (is_alkynyloxy(a_ref,nb[i]))) OR_count++;
      if (is_sulfanyl(a_ref,nb[i])) SH_count++;
      if ((is_alkylsulfanyl(a_ref,nb[i])) || (is_arylsulfanyl(a_ref,nb[i])) ||
         (is_alkenylsulfanyl(a_ref,nb[i])) || (is_alkynylsulfanyl(a_ref,nb[i]))) SR_count++;
      if (!strcmp(atom[(nb[i])]->atype,"N3 ") || !strcmp(atom[(nb[i])]->atype,"NAM")) N_count++;
  }
  if (OH_count == 2) {
      fg[fg_carbonyl_hydrate] = true;
      if (opt_pos) add2fgloc(fg_carbonyl_hydrate,a_ref);  // v0.5
  }
  if ((OH_count == 1) && (OR_count == 1)) { 
      fg[fg_hemiacetal] = true;
      if (opt_pos) add2fgloc(fg_hemiacetal,a_ref);  // v0.5
  }
  if (OR_count == 2) { 
      fg[fg_acetal] = true;
      if (opt_pos) add2fgloc(fg_acetal,a_ref);  // v0.5
  }
  if (((OH_count == 1) || (OR_count == 1)) && (N_count == 1)) { 
      fg[fg_hemiaminal] = true;  
      if (opt_pos) add2fgloc(fg_hemiaminal,a_ref);  // v0.5
  }
  if (N_count == 2) {
      fg[fg_aminal] = true;  
      if (opt_pos) add2fgloc(fg_aminal,a_ref);  // v0.5
  }
  if (((SH_count == 1) || (SR_count == 1)) && (N_count == 1)) { 
      fg[fg_thiohemiaminal] = true;  
      if (opt_pos) add2fgloc(fg_thiohemiaminal,a_ref);  // v0.5
  }
  if ((SR_count == 2) || ((OR_count == 1) && (SR_count == 1))) { 
      fg[fg_thioacetal] = true;  
      if (opt_pos) add2fgloc(fg_thioacetal,a_ref);  // v0.5
  }
}

void chk_carboxyl_deriv_sp3(int a_ref) {
int  i, nb_count, k ;
int  nb[max_neighbors];
int  OR_count, OH_count, N_count;  // oh_count new in v0.3c
int  electroneg_count;             // new in v0.3j
int  hal_count       ;
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_ref) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_ref)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  OR_count = 0; OH_count = 0; N_count = 0;
  electroneg_count = 0;
  hal_count        = 0;
  for (i =0; i<atom[a_ref]->neighbor_count; i++) {// v0.3j
      if (is_electroneg(atom[(nb[i])]->element)) electroneg_count++; 
      if (!strcmp(atom[(nb[i])]->element, "F ") || !strcmp(atom[(nb[i])]->element,"CL") ||
        !strcmp(atom[(nb[i])]->element,"BR") || !strcmp(atom[(nb[i])]->element,"I ")) hal_count++;
      if ((is_alkoxy(a_ref,nb[i])) || (is_aryloxy(a_ref,nb[i])) || (is_siloxy(a_ref,nb[i]))) {
        OR_count++;}
      if (is_hydroxy(a_ref,nb[i])) OH_count++;  // new in v0.3c   
      if (!strcmp(atom[(nb[i])]->atype,"N3 ") || !strcmp(atom[(nb[i])]->atype,"NAM")) N_count++;
  }
  //if (or_count + n_count > 1) then fg[fg_orthocarboxylic_acid_deriv] := true;  // until v0.3i
  if ((electroneg_count == 3) && (hal_count < 3)) {
      fg[fg_orthocarboxylic_acid_deriv] = true;  // v0.3j
      if (opt_pos) add2fgloc(fg_orthocarboxylic_acid_deriv,a_ref);   // v0.5
  }
  if (OR_count == 3) { 
      fg[fg_carboxylic_acid_orthoester] = true;
      if (opt_pos) add2fgloc(fg_carboxylic_acid_orthoester,a_ref);   // v0.5
  }
  if ((OR_count == 2) && (N_count == 1)) { 
      fg[fg_carboxylic_acid_amide_acetal] = true;
      if (opt_pos) add2fgloc(fg_carboxylic_acid_amide_acetal,a_ref);   // v0.5
  }
  if ((OH_count > 0) && (OH_count + OR_count + N_count == 3)) { 
      fg[fg_orthocarboxylic_acid_deriv] = true;  // new in v0.3c
      if (opt_pos) add2fgloc(fg_orthocarboxylic_acid_deriv,a_ref);   // v0.5
  }
}

void chk_anhydride(int a_ref) {
int  i, nb_count, k;
int  nb[max_neighbors];
int  acyl_count;
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_ref) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_ref)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  acyl_count = 0;
  for (i=0;i<atom[a_ref]->neighbor_count;i++){
      if (is_acyl_gen(a_ref,nb[i])) acyl_count++; // v0.4b replaced is_acyl() by is_acyl_gen()
  }
  if ((acyl_count == 2) && !strcmp(atom[a_ref]->atype,"O3 ")) {
      fg[fg_carboxylic_acid_deriv]     = true;
      fg[fg_carboxylic_acid_anhydride] = true;
      if (opt_pos) {   // v0.5
          add2fgloc(fg_carboxylic_acid_anhydride,a_ref);
          for (i=0;i<atom[a_ref]->neighbor_count;i++) {
              if (is_acyl(a_ref,nb[i])) add2fgloc(fg_carboxylic_acid_deriv,nb[i]);  // use is_acyl, not is_acyl_gen here!
          }
       }
   }
}

void chk_imide(int a_ref) {
int  i, nb_count,k;
int  nb[max_neighbors];
int  acyl_count;
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_ref) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_ref)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  acyl_count = 0;
  for (i=0; i<atom[a_ref]->neighbor_count;i++) {
      if ((is_acyl_gen(a_ref,nb[i])) || (is_carbamoyl(a_ref,nb[i]))) acyl_count++;  // v0.3j
  }
  if ((acyl_count >= 2) && !strcmp(atom[a_ref]->element,"N ")) {  // v0.3j: accept also N-acyl-imides    begin
      fg[fg_carboxylic_acid_deriv] = true;
      fg[fg_carboxylic_acid_imide] = true;
      if (atom[a_ref]->neighbor_count == 2) fg[fg_carboxylic_acid_unsubst_imide] = true;
      if (atom[a_ref]->neighbor_count == 3) fg[fg_carboxylic_acid_subst_imide] = true;
      if (opt_pos) {   // v0.5
          add2fgloc(fg_carboxylic_acid_imide,a_ref);
          if (atom[a_ref]->neighbor_count == 2) add2fgloc(fg_carboxylic_acid_unsubst_imide,a_ref);
          if (atom[a_ref]->neighbor_count == 3) add2fgloc(fg_carboxylic_acid_subst_imide,a_ref);
          for (i=0;i<atom[a_ref]->neighbor_count;i++) {
              if (is_acyl(a_ref,nb[i])) add2fgloc(fg_carboxylic_acid_deriv,nb[i]);  // use is_acyl, not is_acyl_gen here!
          }
      }
  }
}

void get_ringstat(int r_id){
int  i,j;
RINGPATH_TYPE  testring;
int  ring_size;
int  a_ref;
int  nN, nO, nS, nX;
  nN = 0; nO = 0; nS = 0; nX = 0;
  if ((r_id < 0) || (r_id > n_rings-1)) return;
  memset(testring,0,sizeof(RINGPATH_TYPE));
  ring_size= ringprop[r_id]->size;  // v0.3j
  for (j =0;j<ring_size;j++) testring[j] =  *(*(ring+r_id)+j);  // v0.3j
  if (ring_size > 2) {
      for (i=0;i<ring_size;i++) {
          a_ref = testring[i];
          if (strcmp(atom[a_ref]->element,"C ") && strcmp(atom[a_ref]->element,"A ")){ 
              nX++;   // general heteroatom count
              if (!strcmp(atom[a_ref]->element,"N ")) nN++;
              if (!strcmp(atom[a_ref]->element,"O ")) nO++;
              if (!strcmp(atom[a_ref]->element,"S ")) nS++;
          }
      }
      if (nN > 0){
          molstat.n_rN++;
          if (nN == 1) molstat.n_rN1++;
          if (nN == 2) molstat.n_rN2++;
          if (nN > 2) molstat.n_rN3p++;
      }
      if (nO > 0){ 
          molstat.n_rO++;
          if (nO == 1) molstat.n_rO1++;
          if (nO == 2) molstat.n_rO2p++;
      }
      if (nS > 0) molstat.n_rS++;
      if (nX > 0) molstat.n_rX++;
      // general ringsize descriptors; v0.3m
    if (ring_size>=13) {molstat.n_r13p++;} else {
      switch (ring_size)
      {
       case 3 : molstat.n_r3++;break;
       case 4 : molstat.n_r4++;break;
       case 5 : molstat.n_r5++;break;
       case 6 : molstat.n_r6++;break;
       case 7 : molstat.n_r7++;break;
       case 8 : molstat.n_r8++;break;
       case 9 : molstat.n_r9++;break;
       case 10: molstat.n_r10++;break;
       case 11: molstat.n_r11++;break;
       case 12: molstat.n_r12++;break;
       default: break;
      }
    }
  }
}

void get_molstat() {
int  i;
int  a1, a2 ;
char  btype;
int  hbc;
int  n_b2formal;  // new in v0.2e
  if (n_atoms == 0) return;
//atoms:
  for (i=1;i<=n_atoms;i++) {
     if (atom[i]->heavy) {
          if (!strcmp(atom[i]->atype,"C1 ")) molstat.n_C1++;
          if (!strcmp(atom[i]->atype,"C2 ") || !strcmp(atom[i]->atype,"CAR")) molstat.n_C2++;
          if (!strcmp(atom[i]->element,"C ")) molstat.n_C++;
          if (!strcmp(atom[i]->atype,"O2 ")) molstat.n_O2++;
          if (!strcmp(atom[i]->atype,"O3 ")) molstat.n_O3++;
          if (!strcmp(atom[i]->atype,"N1 ")) molstat.n_N1++;
          if ((!strcmp(atom[i]->atype,"N2 ")) || (!strcmp(atom[i]->atype,"NAR")) || 
            (!strcmp(atom[i]->atype,"NAM") && (atom[i]->arom == true))) molstat.n_N2++;  // v0.3n
          if (!strcmp(atom[i]->atype,"N3 ") || !strcmp(atom[i]->atype,"NPL") || !strcmp(atom[i]->atype,"N3+") || 
            (!strcmp(atom[i]->atype,"NAM") && (atom[i]->arom == false))) molstat.n_N3++;  // v0.3n
          if (!strcmp(atom[i]->element,"A ")) molstat.n_QA++;  // query atom
          if (!strcmp(atom[i]->element,"Q ")) molstat.n_QA++;  // query atom
          if (!strcmp(atom[i]->element,"X ")) molstat.n_QA++;  // query atom (halogen)  // v0.3p
          if (!strcmp(atom[i]->element,"S ")) molstat.n_S++;
          if (!strcmp(atom[i]->element,"SE")) molstat.n_SeTe++;
          if (!strcmp(atom[i]->element,"TE")) molstat.n_SeTe++;
          if (!strcmp(atom[i]->element,"F ")) molstat.n_F++;
          if (!strcmp(atom[i]->element,"CL")) molstat.n_Cl++;
          if (!strcmp(atom[i]->element,"BR")) molstat.n_Br++;
          if (!strcmp(atom[i]->element,"I ")) molstat.n_I++;
          if (!strcmp(atom[i]->element,"P ")) molstat.n_P++;
          if (!strcmp(atom[i]->element,"B ")) molstat.n_B++;
    // check for known metals
          if (atom[i]->metal) molstat.n_Met++;  // v0.3l
          // still missing: unknown elements
          // check number of heteroatom bonds per C atom
          if (!strcmp(atom[i]->element,"C ")) {
             hbc = raw_hetbond_count(i);   // new in v0.2j (replaces hetbond_count)
             if (hbc >= 1) molstat.n_CHB1p++;
             if (hbc >= 2) molstat.n_CHB2p++;
             if (hbc >= 3) molstat.n_CHB3p++;
             if (hbc == 4) molstat.n_CHB4++;
          }
          if (atom[i]->formal_charge != 0) { 
            molstat.n_chg++;
            n_charges++;
          }
   // check for "other" elements;  v0.3l
   //(elem = 'F ') or (elem = 'CL') or (elem = 'BR') or (elem = 'I ') or  // leave halogens as type X, v0.3m
         if (!(atom[i]->metal || !strcmp(atom[i]->element,"C ") || 
         !strcmp(atom[i]->element,"N ") || 
         !strcmp(atom[i]->element,"O ") || !strcmp(atom[i]->element,"S ") ||
         !strcmp(atom[i]->element,"SE") || !strcmp(atom[i]->element,"TE") || 
         !strcmp(atom[i]->element,"P ") || !strcmp(atom[i]->element,"B ") ||
         !strcmp(atom[i]->element,"A ") || !strcmp(atom[i]->element,"Q "))) molstat.n_X++;
/*         {$IFDEF extended_molstat}
         if (!strcmp(elem = 'LI') or (!strcmp(elem = 'NA') or (!strcmp(elem = 'K ') or (!strcmp(elem = 'RB') or
            (!strcmp(elem = 'CS') or (!strcmp(elem = 'FR') then inc(n_psg01); 
         if (!strcmp(elem = 'BE') or (!strcmp(elem = 'MG') or (!strcmp(elem = 'CA') or (!strcmp(elem = 'SR') or
            (elem = 'BA') or (elem = 'RA') then inc(n_psg02); 
         if (elem = 'B ') or (elem = 'AL') or (elem = 'GA') or (elem = 'IN') or
            (elem = 'TL') then inc(n_psg13); 
         if (elem = 'C ') or (elem = 'SI') or (elem = 'GE') or (elem = 'SN') or
            (elem = 'PB') then inc(n_psg14); 
         if (elem = 'N ') or (elem = 'P ') or (elem = 'AS') or (elem = 'SB') or
            (elem = 'BI') then inc(n_psg15); 
         if (elem = 'O ') or (elem = 'S ') or (elem = 'SE') or (elem = 'TE') or
            (elem = 'PO') then inc(n_psg16); 
         if (elem = 'F ') or (elem = 'CL') or (elem = 'BR') or (elem = 'I ') or
            (elem = 'AT') then inc(n_psg17); 
         if (elem = 'HE') or (elem = 'NE') or (elem = 'AR') or (elem = 'KR') or
            (elem = 'XE') or (elem = 'RN') then inc(n_psg18); 
         if (elem = 'SC') or (elem = 'Y ') or (elem = 'LU') or (elem = 'LR') or
            (elem = 'TI') or (elem = 'ZR') or (elem = 'HF') or (elem = 'RF') or
            (elem = 'V ') or (elem = 'NB') or (elem = 'TA') or (elem = 'DB') or
            (elem = 'CR') or (elem = 'MO') or (elem = 'W ') or (elem = 'SG') or
            (elem = 'MN') or (elem = 'TC') or (elem = 'RE') or (elem = 'BH') or
            (elem = 'FE') or (elem = 'RU') or (elem = 'OS') or (elem = 'HS') or
            (elem = 'CO') or (elem = 'RH') or (elem = 'IR') or (elem = 'MT') or
            (elem = 'NI') or (elem = 'PD') or (elem = 'PT') or (elem = 'DS') or
            (elem = 'CU') or (elem = 'AG') or (elem = 'AU') or (elem = 'RG') or
            (elem = 'ZN') or (elem = 'CD') or (elem = 'HG') then inc(n_pstm); 
         if (elem = 'LA') or (elem = 'CE') or (elem = 'PR') or (elem = 'ND') or
            (elem = 'PM') or (elem = 'SM') or (elem = 'EU') or (elem = 'GD') or
            (elem = 'TB') or (elem = 'DY') or (elem = 'HO') or (elem = 'ER') or
            (elem = 'TM') or (elem = 'YB') or
            (elem = 'AC') or (elem = 'TH') or (elem = 'PA') or (elem = 'U ') or
            (elem = 'NP') or (elem = 'PU') or (elem = 'AM') or (elem = 'CM') or
            (elem = 'BK') or (elem = 'CF') or (elem = 'ES') or (elem = 'FM') or
            (elem = 'MD') or (elem = 'NO') then inc(n_psla); 
         {$ENDIF}*/
    } // is heavy
 }  // atoms
 if (n_bonds > 0) {
    for (i =0; i<n_bonds;i++) {
        a1 = bond[i]->a1; a2 = bond[i]->a2;
        btype = bond[i]->btype;
        if (bond[i]->arom) {molstat.n_bar++;} else { // v0.3n: ignore bonds to (explicit) hydrogens
            if ((btype == 'S') && (atom[a1]->heavy && atom[a2]->heavy)) molstat.n_b1++;
            if (btype == 'D') molstat.n_b2++;
            if (btype == 'T') molstat.n_b3++;
        }
        if ((!strcmp(atom[a1]->element,"C ") && !strcmp(atom[a2]->element,"O ")) || 
        (!strcmp(atom[a1]->element,"O ") && !strcmp(atom[a2]->element,"C "))) {
            if (btype == 'S') molstat.n_C1O++;
            if (btype == 'D') molstat.n_C2O++;
        }
        if ((!strcmp(atom[a1]->element,"C ") && !strcmp(atom[a2]->element,"N ")) || 
        (!strcmp(atom[a1]->element,"N ") && !strcmp(atom[a2]->element,"C "))) molstat.n_CN++;
        if (strcmp(atom[a1]->element,"C ") && (atom[a1]->heavy) 
        && strcmp(atom[a2]->element,"C ") && (atom[a2]->heavy)) molstat.n_XY++;
              // new in v0.3n: number of bonds belonging to more than one ring
        if (bond[i]->ring_count > 1) molstat.n_br2p++;
     }
  } // bonds
  if (n_rings > 0) {
    n_b2formal = 0;        // v0.3n
    n_countablerings = 0;  // v0.3n
    for (i =0;i<n_rings;i++) {
        if (ringprop[i]->envelope == false) n_countablerings++;  // v0.3n
        if (is_arene(i) && (ringprop[i]->envelope == false)) {   // v0.3n: ignore envelope rings
            molstat.n_rar++;
              if ((ringprop[i]->size == 6) && (is_heterocycle(i) == false)) molstat.n_rBz++;  // v0.3l
        }
        get_ringstat(i);
        if ((ringprop[i]->arom == true) && (ringprop[i]->envelope == false)) n_b2formal++;  // new in v0.3n; replaces assignment below
    //n_b2formal := n_rar;  // new in v0.2e; adds 1 formal double bond for each aromatic ring
    // in order to allow an isolated double bond in the needle
    // to be matched as a ring fragment of an aromatic ring
    }
      if (n_b2formal > (molstat.n_bar/2)) n_b2formal = molstat.n_bar/2;
      molstat.n_b2 = molstat.n_b2 + n_b2formal;
  } // rings
}

void fix_ssr_ringcounts() { // new in v0.3n
  // if SAR -> SSR fallback happens, set some molstat values
  // to a maximum (ring counts for various ring sizes);
  // this should be necessary only for ring sizes which
  // are a) too large for the SSR (depending on ssr_vringsize)
  // and b) which are likely to contain "envelope rings"
  // (size 6 and above)
//  if (molstat.n_r3 = 0) then molstat.n_r3 := max_rings;
//  if (molstat.n_r4 = 0) then molstat.n_r4 := max_rings;
//  if (molstat.n_r5 = 0) then molstat.n_r5 := max_rings;
  if (molstat.n_r6 == 0)   molstat.n_r6 = max_rings;
  if (molstat.n_r7 == 0)   molstat.n_r7 = max_rings;
  if (molstat.n_r8 == 0)   molstat.n_r8 = max_rings;
  if (molstat.n_r9 == 0)   molstat.n_r9 = max_rings;
  if (molstat.n_r10 == 0)  molstat.n_r10 = max_rings;
  if (molstat.n_r11 == 0)  molstat.n_r11 = max_rings;
  if (molstat.n_r12 == 0)  molstat.n_r12 = max_rings;
  if (molstat.n_r13p == 0) molstat.n_r13p = max_rings;
}

void write_molstat() {
  if (auto_ssr) fix_ssr_ringcounts();  // v0.3n
    printf("n_atoms:%d;",n_heavyatoms);  // count only non-H atoms (some molfiles contain explicit H's)
    if (opt_molstat_v == false) {  // v0.4d
      if (n_bonds > 0) printf("n_bonds:%d;",n_heavybonds);  // count only bonds between non-H atoms
      if (n_rings > 0) printf("n_rings:%d;",n_rings);  // changed to non-envelope rings in v0.3n
      //      if n_QA    > 0 then write('n_QA:',n_QA,';');
      //      if n_QB    > 0 then write('n_QB:',n_QB,';');
      if (opt_chg && (molstat.n_chg > 0)) printf("n_chg:%d;",molstat.n_chg);   // v0.3p
      if (molstat.n_C1 > 0) printf("n_C1:%d;",molstat.n_C1);
      if (molstat.n_C2    > 0) printf("n_C2:%d;",molstat.n_C2);
      // requirement of a given number of sp3 carbons might be too restrictive,
      // so we use the total number of carbons instead  (initially used variable n_C3 is now n_C)
      if (molstat.n_C     > 0) printf("n_C:%d;",molstat.n_C);
      if (molstat.n_CHB1p > 0) printf("n_CHB1p:%d;",molstat.n_CHB1p);
      if (molstat.n_CHB2p > 0) printf("n_CHB2p:%d;",molstat.n_CHB2p);
      if (molstat.n_CHB3p > 0) printf("n_CHB3p:%d;",molstat.n_CHB3p);
      if (molstat.n_CHB4  > 0) printf("n_CHB4:%d;",molstat.n_CHB4);
      if (molstat.n_O2    > 0) printf("n_O2:%d;",molstat.n_O2);
      if (molstat.n_O3    > 0) printf("n_O3:%d;",molstat.n_O3);
      if (molstat.n_N1    > 0) printf("n_N1:%d;",molstat.n_N1);
      if (molstat.n_N2    > 0) printf("n_N2:%d;",molstat.n_N2);
      if (molstat.n_N3    > 0) printf("n_N3:%d;",molstat.n_N3);
      if (molstat.n_S     > 0) printf("n_S:%d;",molstat.n_S);
      if (molstat.n_SeTe  > 0) printf("n_SeTe:%d;",molstat.n_SeTe);
      if (molstat.n_F     > 0) printf("n_F:%d;",molstat.n_F);
      if (molstat.n_Cl    > 0) printf("n_Cl:%d;",molstat.n_Cl);
      if (molstat.n_Br    > 0) printf("n_Br:%d;",molstat.n_Br);
      if (molstat.n_I     > 0) printf("n_I:%d;",molstat.n_I);
      if (molstat.n_P     > 0) printf("n_P:%d;",molstat.n_P);
      if (molstat.n_B     > 0) printf("n_B:%d;",molstat.n_B);
      if (molstat.n_Met   > 0) printf("n_Met:%d;",molstat.n_Met);
      if (molstat.n_X     > 0) printf("n_X:%d;",molstat.n_X);
      if (molstat.n_b1    > 0) printf("n_b1:%d;",molstat.n_b1);
      if (molstat.n_b2    > 0) printf("n_b2:%d;",molstat.n_b2);
      if (molstat.n_b3    > 0) printf("n_b3:%d;",molstat.n_b3);
      if (molstat.n_bar   > 0) printf("n_bar:%d;",molstat.n_bar);
      if (molstat.n_C1O   > 0) printf("n_C1O:%d;",molstat.n_C1O);
      if (molstat.n_C2O   > 0) printf("n_C2O:%d;",molstat.n_C2O);
      if (molstat.n_CN    > 0) printf("n_CN:%d;",molstat.n_CN);
      if (molstat.n_XY    > 0) printf("n_XY:%d;",molstat.n_XY);
      if (molstat.n_r3    > 0) printf("n_r3:%d;",molstat.n_r3);
      if (molstat.n_r4    > 0) printf("n_r4:%d;",molstat.n_r4);
      if (molstat.n_r5    > 0) printf("n_r5:%d;",molstat.n_r5);
      if (molstat.n_r6    > 0) printf("n_r6:%d;",molstat.n_r6);
      if (molstat.n_r7    > 0) printf("n_r7:%d;",molstat.n_r7);
      if (molstat.n_r8    > 0) printf("n_r8:%d;",molstat.n_r8);
      if (molstat.n_r9    > 0) printf("n_r9:%d;",molstat.n_r9);
      if (molstat.n_r10   > 0) printf("n_r10:%d;",molstat.n_r10);
      if (molstat.n_r11   > 0) printf("n_r11:%d;",molstat.n_r11);
      if (molstat.n_r12   > 0) printf("n_r12:%d;",molstat.n_r12);
      if (molstat.n_r13p  > 0) printf("n_r13p:%d;",molstat.n_r13p);
      if (molstat.n_rN    > 0) printf("n_rN:%d;",molstat.n_rN);
      if (molstat.n_rN1   > 0) printf("n_rN1:%d;",molstat.n_rN1);
      if (molstat.n_rN2   > 0) printf("n_rN2:%d;",molstat.n_rN2);
      if (molstat.n_rN3p  > 0) printf("n_rN3p:%d;",molstat.n_rN3p);
      if (molstat.n_rO    > 0) printf("n_rO:%d;",molstat.n_rO);
      if (molstat.n_rO1   > 0) printf("n_rO1:%d;",molstat.n_rO1);
      if (molstat.n_rO2p  > 0) printf("n_rO2p:%d;",molstat.n_rO2p);
      if (molstat.n_rS    > 0) printf("n_rS:%d;",molstat.n_rS);
      if (molstat.n_rX    > 0) printf("n_rX:%d;",molstat.n_rX);
      if (molstat.n_rar   > 0) printf("n_rar:%d;",molstat.n_rar);
/*      {$IFDEF extended_molstat}
          if n_rBz   > 0 then write('n_rbz:',n_rBz,';');
          if n_br2p  > 0 then write('n_br2p:',n_br2p,';');
          if n_psg01 > 0 then write('n_psg01:',n_psg01,';');
          if n_psg02 > 0 then write('n_psg02:',n_psg02,';');
          if n_psg13 > 0 then write('n_psg13:',n_psg13,';');
          if n_psg14 > 0 then write('n_psg14:',n_psg14,';');
          if n_psg15 > 0 then write('n_psg15:',n_psg15,';');
          if n_psg16 > 0 then write('n_psg16:',n_psg16,';');
          if n_psg17 > 0 then write('n_psg17:',n_psg17,';');
          if n_psg18 > 0 then write('n_psg18:',n_psg18,';');
          if n_pstm > 0 then write('n_pstm:',n_pstm,';');
          if n_psla > 0 then write('n_psla:',n_psla,';');
          {$ENDIF}*/
    } else {
   // full output, even if field is 0;  v0.4d
      printf("n_bonds:%d;",n_heavybonds);   // count only bonds between non-H atoms
      printf("n_rings:%d;",n_rings);  // changed to non-envelope rings in v0.3n
      printf("n_QA:%d;",molstat.n_QA);
      printf("n_QB:%d;",molstat.n_QB);
      if (opt_chg) printf("n_chg:%d;",molstat.n_chg);    // v0.3p
      printf("n_C1:%d;",molstat.n_C1);
      printf("n_C2:%d;",molstat.n_C2);
      // requirement of a given number of sp3 carbons might be too restrictive,
      // so we use the total number of carbons instead  (initially used variable n_C3 is now n_C)
      printf("n_C:%d;",molstat.n_C);
      printf("n_CHB1p:%d;",molstat.n_CHB1p);
      printf("n_CHB2p:%d;",molstat.n_CHB2p);
      printf("n_CHB3p:%d;",molstat.n_CHB3p);
      printf("n_CHB4:%d;",molstat.n_CHB4);
      printf("n_O2:%d;",molstat.n_O2);
      printf("n_O3:%d;",molstat.n_O3);
      printf("n_N1:%d;",molstat.n_N1);
      printf("n_N2:%d;",molstat.n_N2);
      printf("n_N3:%d;",molstat.n_N3);
      printf("n_S:%d;",molstat.n_S);
      printf("n_SeTe:%d;",molstat.n_SeTe);
      printf("n_F:%d;",molstat.n_F);
      printf("n_Cl:%d;",molstat.n_Cl);
      printf("n_Br:%d;",molstat.n_Br);
      printf("n_I:%d;",molstat.n_I);
      printf("n_P:%d;",molstat.n_P);
      printf("n_B:%d;",molstat.n_B);
      printf("n_Met:%d;",molstat.n_Met);
      printf("n_X:%d;",molstat.n_X);
      printf("n_b1:%d;",molstat.n_b1);
      printf("n_b2:%d;",molstat.n_b2);
      printf("n_b3:%d;",molstat.n_b3);
      printf("n_bar:%d;",molstat.n_bar);
      printf("n_C1O:%d;",molstat.n_C1O);
      printf("n_C2O:%d;",molstat.n_C2O);
      printf("n_CN:%d;",molstat.n_CN);
      printf("n_XY:%d;",molstat.n_XY);
      printf("n_r3:%d;",molstat.n_r3);
      printf("n_r4:%d;",molstat.n_r4);
      printf("n_r5:%d;",molstat.n_r5);
      printf("n_r6:%d;",molstat.n_r6);
      printf("n_r7:%d;",molstat.n_r7);
      printf("n_r8:%d;",molstat.n_r8);
      printf("n_r9:%d;",molstat.n_r9);
      printf("n_r10:%d;",molstat.n_r10);
      printf("n_r11:%d;",molstat.n_r11);
      printf("n_r12:%d;",molstat.n_r12);
      printf("n_r13p:%d;",molstat.n_r13p);
      printf("n_rN:%d;",molstat.n_rN);
      printf("n_rN1:%d;",molstat.n_rN1);
      printf("n_rN2:%d;",molstat.n_rN2);
      printf("n_rN3p:%d;",molstat.n_rN3p);
      printf("n_rO:%d;",molstat.n_rO);
      printf("n_rO1:%d;",molstat.n_rO1);
      printf("n_rO2p:%d;",molstat.n_rO2p);
      printf("n_rS:%d;",molstat.n_rS);
      printf("n_rX:%d;",molstat.n_rX);
      printf("n_rar:%d;",molstat.n_rar);
/*      {$IFDEF extended_molstat}
          write('n_rbz:',n_rBz,';');
          write('n_br2p:',n_br2p,';');
          write('n_psg01:',n_psg01,';');
          write('n_psg02:',n_psg02,';');
          write('n_psg13:',n_psg13,';');
          write('n_psg14:',n_psg14,';');
          write('n_psg15:',n_psg15,';');
          write('n_psg16:',n_psg16,';');
          write('n_psg17:',n_psg17,';');
          write('n_psg18:',n_psg18,';');
          write('n_pstm:',n_pstm,';');
          write('n_psla:',n_psla,';');
          {$ENDIF}*/
    }
    printf("\n");
}

void write_molstat_X() {
  if (auto_ssr) fix_ssr_ringcounts() ;  // v0.3n
      printf("%d,%d,%d,",n_heavyatoms,n_heavybonds,n_rings); 
      printf("%d,%d,",molstat.n_QA,molstat.n_QB);
      if (opt_chg) {printf("%d,",molstat.n_chg);} else {printf("0,");}  // v0.3p
      printf("%d,%d,%d,%d,",molstat.n_C1,molstat.n_C2,molstat.n_C,molstat.n_CHB1p);
      printf("%d,%d,%d,%d,",molstat.n_CHB2p,molstat.n_CHB3p,molstat.n_CHB4,molstat.n_O2);
      printf("%d,%d,%d,%d,",molstat.n_O3,molstat.n_N1,molstat.n_N2,molstat.n_N3);
      printf("%d,%d,%d,%d,",molstat.n_S,molstat.n_SeTe,molstat.n_F,molstat.n_Cl); 
      printf("%d,%d,%d,%d,",molstat.n_Br,molstat.n_I,molstat.n_P,molstat.n_B);
      printf("%d,%d,%d,%d,",molstat.n_Met,molstat.n_X,molstat.n_b1,molstat.n_b2); 
      printf("%d,%d,%d,%d,",molstat.n_b3,molstat.n_bar,molstat.n_C1O,molstat.n_C2O); 
      printf("%d,%d,%d,%d,",molstat.n_CN,molstat.n_XY,molstat.n_r3,molstat.n_r4); 
      printf("%d,%d,%d,%d,",molstat.n_r5,molstat.n_r6,molstat.n_r7,molstat.n_r8); 
      printf("%d,%d,%d,%d,",molstat.n_r9,molstat.n_r10,molstat.n_r11,molstat.n_r12);
      printf("%d,%d,%d,%d,",molstat.n_r13p,molstat.n_rN,molstat.n_rN1,molstat.n_rN2);
      printf("%d,%d,%d,%d,",molstat.n_rN3p,molstat.n_rO,molstat.n_rO1,molstat.n_rO2p);
      printf("%d,%d,%d\n",molstat.n_rS,molstat.n_rX,molstat.n_rar);
}

bool normalize_ionic_bonds() {  // v0.3k // changed from a procedure into a function in v0.3m
int  i ;
int  a1, a2 ;
int  fc1, fc2;
char  bt;
bool  res; // v0.3m
  res = false; // v0.3m
  if (n_bonds == 0) return res;
  for (i=0;i<n_bonds;i++) {
      a1= bond[i]->a1;
      a2= bond[i]->a2;
      bt= bond[i]->btype;
      fc1= atom[a1]->formal_charge;
      fc2= atom[a2]->formal_charge;
      if ((fc1 * fc2 == -1) && ((bt == 'S') || (bt == 'D'))) {
          atom[a1]->formal_charge = 0;
          atom[a2]->formal_charge = 0;
          if (!strcmp(atom[a1]->atype,"N3+")) strcpy(atom[a1]->atype,"N3 ");  // v0.3m
          if (!strcmp(atom[a2]->atype,"N3+")) strcpy(atom[a2]->atype,"N3 ");  // v0.3m
          if (bt == 'D') bond[i]->btype = 'T';
          if (bt == 'S') bond[i]->btype = 'D';
          res = true;  // v0.3m
       }
  }
  return res;  // v0.3m (return true if any change was made
}

void chk_so2_deriv(int a_ref){
int  i, nb_count, k;
int  nb[max_neighbors];
int  het_count, O_count, OR_count, hal_count, N_count, C_count;
  if (strcmp(atom[a_ref]->atype,"SO2")) return;
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_ref) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_ref)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  het_count = 0; O_count = 0; OR_count = 0; hal_count = 0; N_count = 0; C_count = 0;
  for (i=0;i<atom[a_ref]->neighbor_count;i++) {
      if (bond[(get_bond(a_ref,nb[i]))]->btype == 'S') {
          if (strcmp(atom[(nb[i])]->element,"C ") && strcmp(atom[(nb[i])]->element,"H ") &&
        strcmp(atom[(nb[i])]->element,"D ") && strcmp(atom[(nb[i])]->element,"DU") &&
 // added 'D ' in v0.3n
        strcmp(atom[(nb[i])]->element,"LP")) het_count++;
          if (!strcmp(atom[(nb[i])]->element,"O ")) {
              O_count++;
              if (is_alkoxy(a_ref,nb[i]) || is_aryloxy(a_ref,nb[i])) OR_count++;
          }
          if (!strcmp(atom[(nb[i])]->element,"N ")) N_count++;
          if (!strcmp(atom[(nb[i])]->element,"C ")) C_count++;
          if (!strcmp(atom[(nb[i])]->element,"F ") || !strcmp(atom[(nb[i])]->element,"CL") ||
        !strcmp(atom[(nb[i])]->element,"BR") || !strcmp(atom[(nb[i])]->element,"I ")) hal_count++;
      }
  }
  if (het_count == 2) {  // sulfuric acid derivative
      fg[fg_sulfuric_acid_deriv] = true;
      if (opt_pos) add2fgloc(fg_sulfuric_acid_deriv,a_ref);   // v0.5
      if (O_count == 2) {
          if (OR_count == 0) {
              fg[fg_sulfuric_acid] = true;
              if (opt_pos) add2fgloc(fg_sulfuric_acid,a_ref);   // v0.5
          }
          if (OR_count == 1) { 
              fg[fg_sulfuric_acid_monoester] = true;
              if (opt_pos) add2fgloc(fg_sulfuric_acid_monoester,a_ref);   // v0.5
          }
          if (OR_count == 2) {
              fg[fg_sulfuric_acid_diester] = true;
              if (opt_pos) add2fgloc(fg_sulfuric_acid_diester,a_ref);   // v0.5
          }
      }
      if (O_count == 1) {
          if ((OR_count == 1) && (N_count == 1)) { 
              fg[fg_sulfuric_acid_amide_ester] = true;
              if (opt_pos) add2fgloc(fg_sulfuric_acid_amide_ester,a_ref);   // v0.5
          }
          if ((OR_count == 0) && (N_count == 1)) { 
              fg[fg_sulfuric_acid_amide] = true;
              if (opt_pos) add2fgloc(fg_sulfuric_acid_amide,a_ref);   // v0.5
          }
      }
      if (N_count == 2) { 
          fg[fg_sulfuric_acid_diamide] = true;
          if (opt_pos) add2fgloc(fg_sulfuric_acid_diamide,a_ref);   // v0.5
      }
      if (hal_count > 0) { 
          fg[fg_sulfuryl_halide] = true;
          if (opt_pos) add2fgloc(fg_sulfuryl_halide,a_ref);   // v0.5
      }
  }
  if ((het_count == 1) && (C_count == 1)) {   // sulfonic acid derivative
      fg[fg_sulfonic_acid_deriv] = true;
      if (opt_pos) add2fgloc(fg_sulfonic_acid_deriv,a_ref);   // v0.5
      if ((O_count == 1) && (OR_count == 0)) { 
          fg[fg_sulfonic_acid] = true;
          if (opt_pos) add2fgloc(fg_sulfonic_acid,a_ref);   // v0.5
      }
      if ((O_count == 1) && (OR_count == 1)) { 
          fg[fg_sulfonic_acid_ester] = true;
          if (opt_pos) add2fgloc(fg_sulfonic_acid_ester,a_ref);   // v0.5
      }
      if (N_count == 1) {
          fg[fg_sulfonamide] = true;
          if (opt_pos) add2fgloc(fg_sulfonamide,a_ref);   // v0.5
      }
      if (hal_count == 1) { 
          fg[fg_sulfonyl_halide] = true;
          if (opt_pos) add2fgloc(fg_sulfonyl_halide,a_ref);   // v0.5
      }
  }
  if ((het_count == 0) && (C_count == 2)) {   // sulfone
      fg[fg_sulfone] = true;
      if (opt_pos) add2fgloc(fg_sulfone,a_ref);   // v0.5
  }
}

void chk_p_deriv(int a_ref) {
int  i, nb_count, k;
int  nb[max_neighbors];
int  het_count, OH_count, OR_count, hal_count, N_count, C_count;
char dbl_het[3];
  if (strcmp(atom[a_ref]->element,"P ")) return;
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_ref) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_ref)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  het_count = 0; OH_count = 0; OR_count= 0; hal_count= 0; N_count= 0; C_count= 0;
  for (i =0; i<atom[a_ref]->neighbor_count;i++) {
      if (bond[(get_bond(a_ref,nb[i]))]->btype == 'D') strcpy(dbl_het,atom[(nb[i])]->element);
      if (bond[(get_bond(a_ref,nb[i]))]->btype == 'S') {
          if (!strcmp(atom[(nb[i])]->element,"C ")) C_count++;
          if (is_hydroxy(a_ref,nb[i])) OH_count++;
          if ((is_alkoxy(a_ref,nb[i])) || (is_aryloxy(a_ref,nb[i]))) OR_count++;
          if (!strcmp(atom[(nb[i])]->element,"N ")) N_count++;
          if (!strcmp(atom[(nb[i])]->element,"F ") || !strcmp(atom[(nb[i])]->element,"CL") ||
        !strcmp(atom[(nb[i])]->element,"BR") || !strcmp(atom[(nb[i])]->element,"I ")) hal_count++;
      }
  }
  het_count = OH_count + OR_count + hal_count + N_count;
  if (!strcmp(atom[a_ref]->atype,"P3D") || !strcmp(atom[a_ref]->atype,"P4 ")) {
      if (!strcmp(dbl_het,"O ")) {
          if (C_count == 0) {
              fg[fg_phosphoric_acid_deriv] = true;
              if (opt_pos) add2fgloc(fg_phosphoric_acid_deriv,a_ref);  // v0.5
              if (OH_count == 3) {
                  fg[fg_phosphoric_acid]        = true;
                  if (opt_pos) add2fgloc(fg_phosphoric_acid,a_ref);  // v0.5
              }
              if (OR_count > 0) { 
                  fg[fg_phosphoric_acid_ester] = true;
                  if (opt_pos) add2fgloc(fg_phosphoric_acid_ester,a_ref);  // v0.5
              }
              if (hal_count > 0) { 
                  fg[fg_phosphoric_acid_halide] = true;            
                  if (opt_pos) add2fgloc(fg_phosphoric_acid_halide,a_ref);  // v0.5
              }
              if (N_count > 0) { 
                  fg[fg_phosphoric_acid_amide] = true;
                  if (opt_pos) add2fgloc(fg_phosphoric_acid_amide,a_ref);  // v0.5
              }
          }
          if (C_count == 1) {
              fg[fg_phosphonic_acid_deriv] = true;
              if (opt_pos) add2fgloc(fg_phosphonic_acid_deriv,a_ref);  // v0.5
              if (OH_count == 2) {
                  fg[fg_phosphonic_acid]        = true;
                  if (opt_pos) add2fgloc(fg_phosphonic_acid,a_ref); } // v0.5
              if (OR_count > 0) { 
                  fg[fg_phosphonic_acid_ester]  = true;
                  if (opt_pos) add2fgloc(fg_phosphonic_acid_ester,a_ref);  // v0.5
              }
              //if (hal_count > 0)  then fg[fg_phosphonic_acid_halide] := true;            
              //if (n_count > 0)    then fg[fg_phosphonic_acid_amide]  := true;
          }
          if (C_count == 3) { 
              fg[fg_phosphinoxide] = true;  
              if (opt_pos) add2fgloc(fg_phosphinoxide,a_ref);  // v0.5
          }
      }
      if (!strcmp(dbl_het,"S ")) {
          if (C_count == 0) {
              fg[fg_thiophosphoric_acid_deriv] = true;
              if (opt_pos) add2fgloc(fg_thiophosphoric_acid_deriv,a_ref);  // v0.5
              if (OH_count == 3) { 
                  fg[fg_thiophosphoric_acid]        = true;
                  if (opt_pos) add2fgloc(fg_thiophosphoric_acid,a_ref);  // v0.5
              }
              if (OR_count > 0) { 
                  fg[fg_thiophosphoric_acid_ester]  = true;
                  if (opt_pos) add2fgloc(fg_thiophosphoric_acid_ester,a_ref);  // v0.5
              }
              if (hal_count > 0) {
                  fg[fg_thiophosphoric_acid_halide] = true;            
                  if (opt_pos) add2fgloc(fg_thiophosphoric_acid_halide,a_ref);  // v0.5
              }
              if (N_count > 0) { 
                  fg[fg_thiophosphoric_acid_amide]  = true;
                  if (opt_pos) add2fgloc(fg_thiophosphoric_acid_amide,a_ref);  // v0.5
              }
          }
      }
  }
//  if (atom^[a_ref].atype = 'P4 ') then fg[fg_phosphoric_acid_deriv] := true;
  if (!strcmp(atom[a_ref]->atype,"P3 ")) {    // changed P3D into P3 in v0.3b
      if ((C_count == 3) && (het_count == 0)) { 
          fg[fg_phosphine] = true;
          if (opt_pos) add2fgloc(fg_phosphine,a_ref);  // v0.5
      }
      if ((C_count == 3) && (OH_count == 1)) { 
          fg[fg_phosphinoxide] = true;
          if (opt_pos) add2fgloc(fg_phosphinoxide,a_ref);  // v0.5
      }
   }
}

void chk_b_deriv(int a_ref) {
int  i, nb_count, k;
int  nb[max_neighbors];
int  het_count, OH_count, OR_count, hal_count, N_count, C_count;
  if (strcmp(atom[a_ref]->element,"B ")) return;
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_ref) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_ref)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  het_count = 0; OH_count = 0; OR_count = 0; hal_count = 0; N_count = 0; C_count = 0;
  for (i=0;i<atom[a_ref]->neighbor_count;i++) {
      if (bond[(get_bond(a_ref,nb[i]))]->btype == 'S') {
          if (!strcmp(atom[(nb[i])]->element,"C ")) {
            C_count++;
          } else if ((strcmp(atom[(nb[i])]->element,"H ") && strcmp(atom[(nb[i])]->element,"D ")
          && strcmp(atom[(nb[i])]->element,"LP"))) {
            het_count++; }// v0.3n: D     
          if (is_hydroxy(a_ref,nb[i])) OH_count++;
          if ((is_alkoxy(a_ref,nb[i])) || (is_aryloxy(a_ref,nb[i]))) OR_count++;  // fixed in v0.3b
          if (!strcmp(atom[(nb[i])]->element,"N ")) N_count++;
          if (!strcmp(atom[(nb[i])]->element,"F ") || !strcmp(atom[(nb[i])]->element,"CL") ||
          !strcmp(atom[(nb[i])]->element,"BR") || !strcmp(atom[(nb[i])]->element,"I ")) hal_count++;
      }
  }
  het_count = OH_count + OR_count + hal_count + N_count;  // fixed in v0.3b
  if ((C_count == 1) && (het_count == 2)) {
      fg[fg_boronic_acid_deriv] = true;
      if (opt_pos) add2fgloc(fg_boronic_acid_deriv,a_ref);  // v0.5
      if (OH_count == 2) {
          fg[fg_boronic_acid]        = true;
          if (opt_pos) add2fgloc(fg_boronic_acid,a_ref);  // v0.5
      }
      if (OR_count > 0) {
          fg[fg_boronic_acid_ester]  = true;
          if (opt_pos) add2fgloc(fg_boronic_acid_ester,a_ref);  // v0.5
      }
  }
}


void chk_ammon(int a_ref) {
int  i, nb_count, k;
int  nb[max_neighbors];
int  het_count, O_count, OR_count, r_count;
char  bt;  // v0.3k
float  bo_sum;
bool  ha;
  if (strcmp(atom[a_ref]->atype,"N3+") && (atom[a_ref]->formal_charge == 0)) return;
  if (strcmp(atom[a_ref]->element,"N ")) return;   // just to be sure;  v0.3i
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_ref) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_ref)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  het_count = 0; O_count = 0; OR_count = 0; r_count = 0;
  bo_sum = 0;
  for (i=0;i<atom[a_ref]->neighbor_count;i++) {
      bt = bond[(get_bond(a_ref,nb[i]))]->btype;  // v0.3k
      ha = atom[nb[i]]->heavy;                   // v0.3k
      if (bt == 'S') {
          if (ha) bo_sum = bo_sum + 1;
          if (strcmp(atom[(nb[i])]->element,"C ") && strcmp(atom[(nb[i])]->element,"H ") &&
             strcmp(atom[(nb[i])]->element,"D ") && strcmp(atom[(nb[i])]->element,"DU")) {   // added 'D ' in v0.3n
               het_count++;
               if (!strcmp(atom[(nb[i])]->element,"O ")) {
                 O_count++;
                 if (atom[(nb[i])]->neighbor_count > 1) OR_count++;
               }
          }
          if ((is_alkyl(a_ref,nb[i])) || (is_aryl(a_ref,nb[i]))   // v0.3k
              || (is_alkenyl(a_ref,nb[i])) || (is_alkynyl(a_ref,nb[i]))) r_count++;
      }
      if (bt == 'D') {
          if (ha) bo_sum = bo_sum + 2;
          if (strcmp(atom[(nb[i])]->element,"C ")) {
              het_count=het_count+2;
              if (!strcmp(atom[(nb[i])]->element,"O ")) { 
                  O_count=O_count+2;
              }
          }
          if (!strcmp(atom[(nb[i])]->element,"C ")) r_count++;
      }
      if ((bt == 'A') && ha) bo_sum = bo_sum + 1.5;
  }   // v0.3k: corrected end of "for ..." loop
  if ((het_count == 0) && (r_count == 4)) { 
      fg[fg_quart_ammonium] = true;
      if (opt_pos) add2fgloc(fg_quart_ammonium,a_ref);   // v0.5
  }
  if ((het_count == 1) && (atom[a_ref]->neighbor_count >= 3)) {
      if ((O_count == 1) && (OR_count == 0) && (bo_sum > 3)) { 
          fg[fg_n_oxide] = true;  // finds only aliphatic N-oxides!
          if (opt_pos) add2fgloc(fg_n_oxide,a_ref);   // v0.5
      }
      if ((((O_count == 1) && (OR_count == 1)) || (O_count == 0))
       && (atom[a_ref]->arom == true)) { 
           fg[fg_quart_ammonium] = true;
           if (opt_pos) add2fgloc(fg_quart_ammonium,a_ref);   // v0.5
      }
  }
}

void chk_12diphenol(int a_view, int a_ref) {
// a_view and a_ref are adjacent ring atoms
int  i, nb_count, k;
int  nb[max_neighbors];
int  OH_count, b_id ;  // v0.5
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_view) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_view)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  OH_count = 0;
  for (i= 0;i<atom[a_view]->neighbor_count;i++) {
      if (bond[(get_bond(a_view,nb[i]))]->btype == 'S') {
          if (is_hydroxy(a_view,nb[i])) OH_count++;
      }
  }
  memset(nb,0,max_neighbors*sizeof(int));
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_ref) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_ref)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  for (i= 0;i< atom[a_ref]->neighbor_count;i++) {
      if (bond[(get_bond(a_ref,nb[i]))]->btype == 'S') {
          if (is_hydroxy(a_ref,nb[i])) OH_count++;
      }
  }
  if (OH_count == 2) {
      fg[fg_1_2_diphenol] = true;
      if (opt_pos) {
          b_id = get_bond(a_view,a_ref);
          add2fgloc(fg_1_2_diphenol,b_id);
      }
  }
}

void chk_arom_fg(int a1, int a2) {
  if ((hetbond_count(a1) == 1) && (hetbond_count(a2) == 1)) chk_12diphenol(a1,a2);
}

void chk_oxo_thioxo_imino_hetarene(int r_id){
int  i,j  ;
RINGPATH_TYPE  testring ;
int  ring_size, a_ref;
  if ((r_id < 0) || (r_id >= n_rings)) return;
  memset(testring,0,sizeof(RINGPATH_TYPE));
  ring_size= ringprop[r_id]->size;  // v0.3j
  //for j := 1 to max_ringsize do if ring^[r_id,j] > 0 then testring[j] := ring^[r_id,j];
  for (j= 0;j<ring_size;j++) testring[j] = *(*(ring+r_id)+j);  // v0.3j
  //ring_size := path_length(testring);
  //if (is_arene(r_id)) and (odd(ring_size) = false) then
  if (is_arene(r_id)) {
      for (i=0;i<ring_size;i++) {
          a_ref = testring[i];
          if (is_oxo_C(a_ref)) {
              fg[fg_oxohetarene]  = true;
              if (opt_pos) add2fgloc(fg_oxohetarene,a_ref);   // v0.5
          }
          if (is_thioxo_C(a_ref)) {
              fg[fg_thioxohetarene] = true;
              if (opt_pos) add2fgloc(fg_thioxohetarene,a_ref);   // v0.5
          }
          if (is_true_exocyclic_imino_C(a_ref,r_id)) {
              fg[fg_iminohetarene] = true;
              if (opt_pos) add2fgloc(fg_iminohetarene,a_ref);   // v0.5
          }
      }
  }
}

char* mkfglabel(int fgnum) {
//char* fp;
static char f[100];
if ((fgnum < 1) || (fgnum > used_fg)) return f;
switch (fgnum) {
  case 1: strcpy(f,"cation");break;
  case 2: strcpy(f,"anion");break;
  case 3: strcpy(f,	"carbonyl compound");break;
  case 4: strcpy(f,	"aldehyde");break;
  case 5: strcpy(f,	"ketone");break;
  case 6: strcpy(f,	"thiocarbonyl compound");break;
  case 7: strcpy(f,	"thioaldehyde");break;
  case 8: strcpy(f,	"thioketone");break;
  case 9: strcpy(f,	"imine");break;
  case 10: strcpy(f,	"hydrazone");break;
  case 11: strcpy(f,	"semicarbazone");break;
  case 12: strcpy(f,	"thiosemicarbazone");break;
  case 13: strcpy(f,	"oxime");break;
  case 14: strcpy(f,	"oxime ether");break;
  case 15: strcpy(f,	"ketene");break;
  case 16: strcpy(f,	"ketene acetal or derivative");break;
  case 17: strcpy(f,	"carbonyl hydrate");break;
  case 18: strcpy(f,	"hemiacetal");break;
  case 19: strcpy(f,	"acetal");break;
  case 20: strcpy(f,	"hemiaminal");break;
  case 21: strcpy(f,	"aminal");break;
  case 22: strcpy(f,	"hemithioaminal");break;
  case 23: strcpy(f,	"thioacetal");break;
  case 24: strcpy(f,	"enamine");break;
  case 25: strcpy(f,	"enol");break;
  case 26: strcpy(f,	"enol ether");break;
  case 27: strcpy(f,	"hydroxy compound");break;
  case 28: strcpy(f,	"alcohol");break;
  case 29: strcpy(f,	"primary alcohol");break;
  case 30: strcpy(f,	"secondary alcohol");break;
  case 31: strcpy(f,	"tertiary alcohol");break;
  case 32: strcpy(f,	"1,2-diol");break;
  case 33: strcpy(f,	"1,2-aminoalcohol");break;
  case 34: strcpy(f,	"phenol or hydroxyhetarene");break;
  case 35: strcpy(f,	"1,2-diphenol");break;
  case 36: strcpy(f,	"enediol");break;
  case 37: strcpy(f,	"ether");break;
  case 38: strcpy(f,	"dialkyl ether");break;
  case 39: strcpy(f,	"alkyl aryl ether");break;
  case 40: strcpy(f,	"diaryl ether");break;
  case 41: strcpy(f,	"thioether");break;
  case 42: strcpy(f,	"disulfide");break;
  case 43: strcpy(f,	"peroxide");break;
  case 44: strcpy(f,	"hydroperoxide");break;
  case 45: strcpy(f,	"hydrazine derivative");break;
  case 46: strcpy(f,	"hydroxylamine");break;
  case 47: strcpy(f,	"amine");break;
  case 48: strcpy(f,	"primary amine");break;
  case 49: strcpy(f,	"primary aliphatic amine (alkylamine)");break;
  case 50: strcpy(f,	"primary aromatic amine");break;
  case 51: strcpy(f,	"secondary amine");break;
  case 52: strcpy(f,	"secondary aliphatic amine (dialkylamine)");break;
  case 53: strcpy(f,	"secondary aliphatic/aromatic amine (alkylarylamine)");break;
  case 54: strcpy(f,	"secondary aromatic amine (diarylamine)");break;
  case 55: strcpy(f,	"tertiary amine");break;
  case 56: strcpy(f,	"tertiary aliphatic amine (trialkylamine)");break;
  case 57: strcpy(f,	"tertiary aliphatic/aromatic amine (alkylarylamine)");break;
  case 58: strcpy(f,	"tertiary aromatic amine (triarylamine)");break;
  case 59: strcpy(f,	"quaternary ammonium salt");break;
  case 60: strcpy(f,	"N-oxide");break;
  case 61: strcpy(f,	"halogen derivative");break;
  case 62: strcpy(f,	"alkyl halide");break;
  case 63: strcpy(f,	"alkyl fluoride");break;
  case 64: strcpy(f,	"alkyl chloride"); break;
  case 65: strcpy(f,	"alkyl bromide");break;
  case 66: strcpy(f,	"alkyl iodide");break;
  case 67: strcpy(f,	"aryl halide"); break;
  case 68: strcpy(f,	"aryl fluoride");break;
  case 69: strcpy(f,	"aryl chloride");break;
  case 70: strcpy(f,	"aryl bromide");break;
  case 71: strcpy(f,	"aryl iodide"); break;
  case 72: strcpy(f,	"organometallic compound");break;
  case 73: strcpy(f,	"organolithium compound");break;
  case 74: strcpy(f,	"organomagnesium compound");break;
  case 75: strcpy(f,	"carboxylic acid derivative");break;
  case 76: strcpy(f,	"carboxylic acid");break;
  case 77: strcpy(f,	"carboxylic acid salt");break;
  case 78: strcpy(f,	"carboxylic acid ester");break;
  case 79: strcpy(f,	"lactone");     break;
  case 80: strcpy(f,	"carboxylic acid amide");break;
  case 81: strcpy(f,	"primary carboxylic acid amide");break;
  case 82: strcpy(f,	"secondary carboxylic acid amide");break;
  case 83: strcpy(f,	"tertiary carboxylic acid amide");break;
  case 84: strcpy(f,	"lactam");      break;
  case 85: strcpy(f,	"carboxylic acid hydrazide");break;
  case 86: strcpy(f,	"carboxylic acid azide");break;
  case 87: strcpy(f,	"hydroxamic acid");break;
  case 88: strcpy(f,	"carboxylic acid amidine");break;
  case 89: strcpy(f,	"carboxylic acid amidrazone");break;
  case 90: strcpy(f,	"carbonitrile");break;
  case 91: strcpy(f,	"acyl halide"); break;
  case 92: strcpy(f,	"acyl fluoride");break;
  case 93: strcpy(f,	"acyl chloride");break;
  case 94: strcpy(f,	"acyl bromide");break;
  case 95: strcpy(f,	"acyl iodide"); break;
  case 96: strcpy(f,	"acyl cyanide");break;
  case 97: strcpy(f,	"imido ester"); break;
  case 98: strcpy(f,	"imidoyl halide");break;
  case 99: strcpy(f,	"thiocarboxylic acid derivative");break;
  case 100: strcpy(f,	"thiocarboxylic acid");break;
  case 101: strcpy(f,	"thiocarboxylic acid ester");break;
  case 102: strcpy(f,	"thiolactone");break;
  case 103: strcpy(f,	"thiocarboxylic acid amide");break;
  case 104: strcpy(f,	"thiolactam"); break;
  case 105: strcpy(f,	"imidothioester");break;
  case 106: strcpy(f,	"oxo(het)arene");break;
  case 107: strcpy(f,	"thioxo(het)arene");break;
  case 108: strcpy(f,	"imino(het)arene");break;
  case 109: strcpy(f,	"orthocarboxylic acid derivative");break;
  case 110: strcpy(f,	"orthoester"); break;
  case 111: strcpy(f,	"amide acetal");break;
  case 112: strcpy(f,	"carboxylic acid anhydride");break;
  case 113: strcpy(f,	"carboxylic acid imide");break;
  case 114: strcpy(f,	"carboxylic acid imide, N-unsubstituted");break;
  case 115: strcpy(f,	"carboxylic acid imide, N-substituted");break;
  case 116: strcpy(f,	"CO2 derivative (general)");break;
  case 117: strcpy(f,	"carbonic acid derivative");break;
  case 118: strcpy(f,	"carbonic acid monoester");break;
  case 119: strcpy(f,	"carbonic acid diester");break;
  case 120: strcpy(f,	"carbonic acid ester halide (alkyl/aryl haloformate)"); break;
  case 121: strcpy(f,	"thiocarbonic acid derivative");break;
  case 122: strcpy(f,	"thiocarbonic acid monoester");break;
  case 123: strcpy(f,	"thiocarbonic acid diester");break;
  case 124: strcpy(f,	"thiocarbonic acid ester halide (alkyl/aryl halothioformate");break;
  case 125: strcpy(f,	"carbamic acid derivative");break;
  case 126: strcpy(f,	"carbamic acid");break;
  case 127: strcpy(f,	"carbamic acid ester (urethane)");break;
  case 128: strcpy(f,	"carbamic acid halide (haloformic acid amide)");break;
  case 129: strcpy(f,	"thiocarbamic acid derivative");break;
  case 130: strcpy(f,	"thiocarbamic acid");break;
  case 131: strcpy(f,	"thiocarbamic acid ester");break;
  case 132: strcpy(f,	"thiocarbamic acid halide (halothioformic acid amide)");break;
  case 133: strcpy(f,	"urea");       break;
  case 134: strcpy(f,	"isourea");    break;
  case 135: strcpy(f,	"thiourea");   break;
  case 136: strcpy(f,	"isothiourea");break;
  case 137: strcpy(f,	"guanidine");  break;
  case 138: strcpy(f,	"semicarbazide");break;
  case 139: strcpy(f,	"thiosemicarbazide");break;
  case 140: strcpy(f,	"azide");      break;
  case 141: strcpy(f,	"azo compound");break;
  case 142: strcpy(f,	"diazonium salt");break;
  case 143: strcpy(f,	"isonitrile"); break;
  case 144: strcpy(f,	"cyanate");    break;
  case 145: strcpy(f,	"isocyanate"); break;
  case 146: strcpy(f,	"thiocyanate");break;
  case 147: strcpy(f,	"isothiocyanate");break;
  case 148: strcpy(f,	"carbodiimide");break;
  case 149: strcpy(f,	"nitroso compound");break;
  case 150: strcpy(f,	"nitro compound");break;
  case 151: strcpy(f,	"nitrite");    break;
  case 152: strcpy(f,	"nitrate");    break;
  case 153: strcpy(f,	"sulfuric acid derivative");break;
  case 154: strcpy(f,	"sulfuric acid");break;
  case 155: strcpy(f,	"sulfuric acid monoester");break;
  case 156: strcpy(f,	"sulfuric acid diester");break;
  case 157: strcpy(f,	"sulfuric acid amide ester");break;
  case 158: strcpy(f,	"sulfuric acid amide");break;
  case 159: strcpy(f,	"sulfuric acid diamide");break;
  case 160: strcpy(f,	"sulfuryl halide");break;
  case 161: strcpy(f,	"sulfonic acid derivative");break;
  case 162: strcpy(f,	"sulfonic acid");break;
  case 163: strcpy(f,	"sulfonic acid ester");break;
  case 164: strcpy(f,	"sulfonamide");break;
  case 165: strcpy(f,	"sulfonyl halide");break;
  case 166: strcpy(f,	"sulfone");    break;
  case 167: strcpy(f,	"sulfoxide");  break;
  case 168: strcpy(f,	"sulfinic acid derivative");break;
  case 169: strcpy(f,	"sulfinic acid");break;
  case 170: strcpy(f,	"sulfinic acid ester");break;
  case 171: strcpy(f,	"sulfinic acid halide");break;
  case 172: strcpy(f,	"sulfinic acid amide");break;
  case 173: strcpy(f,	"sulfenic acid derivative");break;
  case 174: strcpy(f,	"sulfenic acid");break;
  case 175: strcpy(f,	"sulfenic acid ester");break;
  case 176: strcpy(f,	"sulfenic acid halide");break;
  case 177: strcpy(f,	"sulfenic acid amide");break;
  case 178: strcpy(f,	"thiol (sulfanyl compound)");break;
  case 179: strcpy(f,	"alkylthiol"); break;
  case 180: strcpy(f,	"arylthiol (thiophenol)");break;
  case 181: strcpy(f,	"phosphoric acid derivative");break;
  case 182: strcpy(f,	"phosphoric acid");break;
  case 183: strcpy(f,	"phosphoric acid ester");break;
  case 184: strcpy(f,	"phosphoric acid halide");break;
  case 185: strcpy(f,	"phosphoric acid amide");break;
  case 186: strcpy(f,	"thiophosphoric acid derivative");break;
  case 187: strcpy(f,	"thiophosphoric acid");break;
  case 188: strcpy(f,	"thiophosphoric acid ester");break;
  case 189: strcpy(f,	"thiophosphoric acid halide");break;
  case 190: strcpy(f,	"thiophosphoric acid amide");break;
  case 191: strcpy(f,	"phosphonic acid derivative");break;
  case 192: strcpy(f,	"phosphonic acid");break;
  case 193: strcpy(f,	"phosphonic acid ester");break;
  case 194: strcpy(f,	"phosphine");  break;
  case 195: strcpy(f,	"phosphine oxide");break;
  case 196: strcpy(f,	"boronic acid derivative");break;
  case 197: strcpy(f,	"boronic acid");break;
  case 198: strcpy(f,	"boronic acid ester");break;
  case 199: strcpy(f,	"alkene");     break;
  case 200: strcpy(f,	"alkyne");     break;
  case 201: strcpy(f,	"aromatic compound");break;
  case 202: strcpy(f,	"heterocyclic compound");break;
  case 203: strcpy(f,	"alpha-aminoacid");break;
  case 204: strcpy(f,	"alpha-hydroxyacid");break;
  default: break;
}//printf("HERE: %s\n", f);
//fp=f;
return f;
}

void write_fg_text() {
int  i;
  for (i= 1;i<=used_fg;i++) { // fg[0] = 0; leave an empty functional group
  // first define some exceptions
      if (fg[i] == true) { 
          if (
               (i == fg_carbonyl) || 
               (i == fg_thiocarbonyl) || 
               (i == fg_alcohol) ||
               (i == fg_hydroxy) || 
               (i == fg_ether) || 
               (i == fg_amine) || 
               (i == fg_halogen_deriv) ||
               (i == fg_alkyl_halide) || 
               (i == fg_aryl_halide) || 
               (i == fg_carboxylic_acid_deriv) ||
               (i == fg_carboxylic_acid_amide) || 
               (i == fg_acyl_halide) || 
               (i == fg_thiocarboxylic_acid_deriv) ||
               (i == fg_carboxylic_acid_imide) || 
               (i == fg_carbonic_acid_deriv) || 
               (i == fg_carbamic_acid_deriv) ||
               (i == fg_thiocarbamic_acid_deriv) || 
               (i == fg_sulfuric_acid_deriv) || 
               (i == fg_sulfonic_acid_deriv) ||
               (i == fg_sulfinic_acid_deriv) || 
               (i == fg_sulfenic_acid_deriv) || 
               (i == fg_phosphoric_acid_deriv) ||
               (i == fg_thiophosphoric_acid_deriv)
             ) {
              if ((i == fg_hydroxy) && hydroxy_generic)    printf("%s\n",mkfglabel(fg_hydroxy));
              if ((i == fg_ether) && ether_generic)        printf("%s\n",mkfglabel(fg_ether));
              if ((i == fg_amine) && amine_generic)        printf("%s\n",mkfglabel(fg_amine));
              if (i == fg_halogen_deriv) {
                  if ((!fg[fg_alkyl_halide]) && (!fg[fg_aryl_halide]) && (!fg[fg_acyl_halide])) {
                  printf("%s\n",mkfglabel(fg_halogen_deriv));}
              }
              if ((i == fg_carbonic_acid_deriv) && !(fg[fg_carbonic_acid_monoester] || 
                 fg[fg_carbonic_acid_diester] || fg[fg_carbonic_acid_ester_halide]))   {
                    printf("%s\n",mkfglabel(fg_carbonic_acid_deriv)); }
              if ((i == fg_carbamic_acid_deriv) && !(fg[fg_carbamic_acid] || 
                fg[fg_carbamic_acid_ester] || fg[fg_carbamic_acid_halide])) {
                    printf("%s\n",mkfglabel(fg_carbamic_acid_deriv));}
              if ((i == fg_thiocarbamic_acid_deriv) && !(fg[fg_thiocarbamic_acid] ||
                  fg[fg_thiocarbamic_acid_ester] || fg[fg_thiocarbamic_acid_halide])) {
                    printf("%s\n",mkfglabel(fg_thiocarbamic_acid_deriv));}
          } else { printf("%s\n",mkfglabel(i)); } // now treat the rest normally
      } // if fg[i]...
  }  // for i
}

void write_fg_binary () {  // v0.4
const int bsize=32;
int   i, j, n1;
long int  fgincrement ;
long int  fgdecimal ;
  n1 = 0;
  for (i = 0;i<(max_fg/bsize);i++) {
      fgdecimal = 0;
      for (j =1; j<=bsize; j++) {
          fgincrement = 1;
          if (fg[((bsize*i)+j)]) {//printf("'%d'",(bsize*i)+j); 
              n1++;
              fgincrement = fgincrement<<(j-1);
              fgdecimal = fgdecimal + fgincrement;
          }
      }
      if (i > 0) printf(",");
      printf("%ld",fgdecimal);
  }
  printf(";%d\n",n1);
}

void chk_functionalgroups() {
int  i,n ;
int  a1,a2,a_tmp;
char  bt;
int  pos_chg, neg_chg;
  if ((n_atoms < 1) || (n_bonds < 1)) return;
  if (opt_pos) { // v0.5
 // allocate memory for the location table
    if (fgloc!=NULL) { 
      for(n=0;n<used_fg;n++) fgloc[n]=(int *)malloc(max_fgpos*sizeof(int));
      if (fgloc==NULL) {printf("Error allocating memory for fgs!");exit(1);}
    }
  }
  pos_chg = 0; neg_chg = 0;
  for (i =1; i<=n_atoms; i++) {  // a few groups are best discovered in the atom list
      if (!strcmp(atom[i]->atype,"SO2")) chk_so2_deriv(i);
      //if (atom^[i].atype = 'SO ') then fg[fg_sulfoxide] := true;  // do another check in the bond list!!
      if (!strcmp(atom[i]->element,"P ")) chk_p_deriv(i);
      if (!strcmp(atom[i]->element,"B ")) chk_b_deriv(i);
      if (!strcmp(atom[i]->atype,"N3+") || (atom[i]->formal_charge > 0)) chk_ammon(i);
      if (!strcmp(atom[i]->atype,"C3 ") && (hetbond_count(i) == 2)) chk_carbonyl_deriv_sp3(i);
      if (!strcmp(atom[i]->atype,"C3 ") && (hetbond_count(i) == 3)) chk_carboxyl_deriv_sp3(i);
      if (!strcmp(atom[i]->atype,"O3 ") && (atom[i]->neighbor_count == 2)) chk_anhydride(i);
      if ((!strcmp(atom[i]->atype,"N3 ") || !strcmp(atom[i]->atype,"NAM")) &&
         (atom[i]->neighbor_count >= 2)) chk_imide(i);
      if (atom[i]->formal_charge > 0) pos_chg = pos_chg + atom[i]->formal_charge;
      if (atom[i]->formal_charge < 0) neg_chg = neg_chg + atom[i]->formal_charge;
      chk_ion(i);
  }
  for (i=0;i<n_bonds;i++) {  // most groups are best discovered in the bond list
      a1 = bond[i]->a1;
      a2 = bond[i]->a2;
      bt = bond[i]->btype;
      if ((atom[a1]->heavy) && (atom[a2]->heavy)) {// orient bond first!
  if (!strcmp(atom[a1]->element,"H ") || !strcmp(atom[a2]->element,"H ") ||  
    !strcmp(atom[a1]->element,"D ") || !strcmp(atom[a2]->element,"D ")) return; // v0.3n: D
  if (!strcmp(atom[a2]->element, "C ") && strcmp(atom[a1]->element,"C ")) {
    a_tmp = a1;
    a1 = a2;
    a2 = a_tmp;}
  if (!strcmp(atom[a2]->element, atom[a1]->element)) {
      if (hetbond_count(a1) > hetbond_count(a2)) {
        a_tmp = a1;
        a1 = a2;
        a2 = a_tmp;}
  }
  if (strcmp(atom[a2]->element,"C ") && strcmp(atom[a1]->element,"C ") &&
    strcmp(atom[a1]->element, atom[a2]->element)) {
      if (!strcmp(atom[a1]->element, "O ") || !strcmp(atom[a2]->element,"O ")) {
          if (!strcmp(atom[a1]->element,"O ")) {
            a_tmp = a1;
            a1 = a2;
            a2 = a_tmp;}
      }
  }
  if (strcmp(atom[a2]->element,"C ") && strcmp(atom[a1]->element,"C ") &&
    !strcmp(atom[a1]->element,atom[a2]->element)) {
      if ((atom[a2]->neighbor_count - hetbond_count(a2)) > 
          (atom[a1]->neighbor_count - hetbond_count(a1))) {
        a_tmp = a1;
        a1 = a2;
        a2 = a_tmp;}
  }// orient bond
          if (bt == 'T') chk_triple(a1,a2);
          if (bt == 'D') chk_double(a1,a2);
          if (bt == 'S') chk_single(a1,a2);
          if (bond[i]->arom) chk_arom_fg(a1,a2);
      }
  }
  if (n_rings > 0) {
      for (i=0;i<n_rings;i++) {
          chk_oxo_thioxo_imino_hetarene(i);
          if (is_arene(i)) { 
              fg[fg_aromatic] = true;
              if (opt_pos && (ringprop[i]->envelope == false)) add2fgloc(fg_aromatic,i);    // v0.5
          }
          if (is_heterocycle(i)) {
              fg[fg_heterocycle] = true;
              if (opt_pos && (ringprop[i]->envelope == false)) add2fgloc(fg_heterocycle,i);    // v0.5
          }
      }
  }
  if ((pos_chg + neg_chg) > 0) fg[fg_cation] = true;
  if ((pos_chg + neg_chg) < 0) fg[fg_anion] = true;
}

void init_pt() {  // v0.4, v0.4d, pt[0] is empty
  // fills the periodic table with element symbols; index = atomic number (Z)
  strcpy(pt[1].el,"H "); strcpy(pt[2].el,"HE"); strcpy(pt[3].el,"LI"); strcpy(pt[4].el,"BE");  strcpy(pt[5].el,"B "); 
  pt[1].am = 1.0079; pt[2].am = 4.0026;  pt[3].am = 6.941;  pt[4].am = 9.01218;  pt[5].am = 10.81; 

  strcpy(pt[6].el,"C "); strcpy(pt[7].el,"N "); strcpy(pt[8].el,"O "); strcpy(pt[9].el,"F ");  strcpy(pt[10].el,"NE"); 
  pt[6].am = 12.011; pt[7].am = 14.00674;  pt[8].am = 15.9994;  pt[9].am = 18.999840;  pt[10].am = 20.179; 

  strcpy(pt[11].el,"NA"); strcpy(pt[12].el,"MG"); strcpy(pt[13].el,"AL"); strcpy(pt[14].el,"SI");  strcpy(pt[15].el,"P "); 
  pt[11].am = 22.98977; pt[12].am = 24.305;  pt[13].am = 26.98154;  pt[14].am = 28.086;  pt[15].am = 30.97376; 

  strcpy(pt[16].el,"S "); strcpy(pt[17].el,"CL"); strcpy(pt[18].el,"AR"); strcpy(pt[19].el,"K ");  strcpy(pt[20].el,"CA"); 
  pt[16].am = 32.06; pt[17].am = 35.453;  pt[18].am = 39.948;  pt[19].am = 39.098;  pt[20].am = 40.08; 

  strcpy(pt[21].el,"SC"); strcpy(pt[22].el,"TI"); strcpy(pt[23].el,"V "); strcpy(pt[24].el,"CR");  strcpy(pt[25].el,"MN"); 
  pt[21].am = 44.9559; pt[22].am = 47.90;  pt[23].am = 50.9414;  pt[24].am = 51.996;  pt[25].am = 54.9380; 

  strcpy(pt[26].el,"FE"); strcpy(pt[27].el,"CO"); strcpy(pt[28].el,"NI"); strcpy(pt[29].el,"CU");  strcpy(pt[30].el,"ZN"); 
  pt[26].am = 55.847; pt[27].am = 58.9332;  pt[28].am = 58.70;  pt[29].am = 63.546;  pt[30].am = 65.38; 

  strcpy(pt[31].el,"GA"); strcpy(pt[32].el,"GE"); strcpy(pt[33].el,"AS"); strcpy(pt[34].el,"SE");  strcpy(pt[35].el,"BR"); 
  pt[31].am = 	69.72; pt[32].am = 72.59;  pt[33].am = 74.9216;  pt[34].am = 78.96;  pt[35].am = 79.904; 

  strcpy(pt[36].el,"KR"); strcpy(pt[37].el,"RB");  strcpy(pt[38].el,"SR");  strcpy(pt[39].el,"Y ");  strcpy(pt[40].el,"Zr"); 
  pt[36].am = 83.80; pt[37].am = 85.4678;  pt[38].am = 87.62;  pt[39].am = 88.9059;  pt[40].am = 91.22; 

  strcpy(pt[41].el,"NB"); strcpy(pt[42].el,"MO");  strcpy(pt[43].el,"TC");  strcpy(pt[44].el,"RU");  strcpy(pt[45].el,"RH"); 
  pt[41].am = 92.9064; pt[42].am = 95.94;  pt[43].am = 97;  pt[44].am = 101.07;  pt[45].am = 102.9055; 

  strcpy(pt[46].el,"PD"); strcpy(pt[47].el,"AG");  strcpy(pt[48].el,"CD");  strcpy(pt[49].el,"IN");  strcpy(pt[50].el,"SN"); 
  pt[46].am = 106.4; pt[47].am = 106.4;  pt[48].am = 112.40;  pt[49].am = 114.82;  pt[50].am = 118.69; 

  strcpy(pt[51].el,"SB"); strcpy(pt[52].el,"TE");  strcpy(pt[53].el,"I ");  strcpy(pt[54].el,"XE");  strcpy(pt[55].el,"CS"); 
  pt[51].am = 121.75; pt[52].am = 127.60;  pt[53].am = 126.9045;  pt[54].am = 131.30;  pt[55].am = 132.9054; 

  strcpy(pt[56].el,"BA"); strcpy(pt[57].el,"LA");  strcpy(pt[58].el,"CE");  strcpy(pt[59].el,"PR");  strcpy(pt[60].el,"ND"); 
  pt[56].am = 137.34; pt[57].am = 138.9055;  pt[58].am = 140.12;  pt[59].am = 140.9077;  pt[60].am = 144.24; 

  strcpy(pt[61].el,"PM"); strcpy(pt[62].el,"SM"); strcpy(pt[63].el,"EU"); strcpy(pt[64].el,"GD"); strcpy(pt[65].el,"TB"); 
  pt[61].am = 145; pt[62].am = 150.4;  pt[63].am = 151.96;  pt[64].am = 157.25;  pt[65].am = 158.9254; 

  strcpy(pt[66].el,"DY"); strcpy(pt[57].el,"HO");  strcpy(pt[68].el,"ER"); strcpy(pt[69].el,"TM"); strcpy(pt[70].el,"YB"); 
  pt[66].am = 162.50; pt[57].am = 164.9304;  pt[68].am = 167.26;  pt[69].am = 168.9342;  pt[70].am = 173.04; 

  strcpy(pt[71].el,"LU"); strcpy(pt[72].el,"HF");  strcpy(pt[73].el,"TA");  strcpy(pt[74].el,"W ");  strcpy(pt[75].el,"RE"); 
  pt[71].am = 174.96; pt[72].am = 178.49;  pt[73].am = 180.9479;  pt[74].am = 183.5;  pt[75].am = 186.207; 

  strcpy(pt[76].el,"OS"); strcpy(pt[77].el,"IR");  strcpy(pt[78].el,"PT");  strcpy(pt[79].el,"AU"); strcpy(pt[80].el,"HG"); 
  pt[76].am = 190.2; pt[77].am = 192.22;  pt[78].am = 195.09;  pt[79].am = 196.9665;  pt[80].am = 200.59; 

  strcpy(pt[81].el,"TL"); strcpy(pt[82].el,"PB");  strcpy(pt[83].el,"BI");  strcpy(pt[84].el,"PO");  strcpy(pt[85].el,"AT"); 
  pt[81].am = 204.37; pt[82].am = 207.2;  pt[83].am = 208.9804;  pt[84].am = 209;  pt[85].am = 210; 

  strcpy(pt[86].el,"RN"); strcpy(pt[87].el,"FR");  strcpy(pt[88].el,"RA");  strcpy(pt[89].el,"AC");  strcpy(pt[90].el,"TH"); 
  pt[86].am = 222; pt[87].am = 223;  pt[88].am = 226.0254;  pt[89].am = 227;  pt[90].am = 232.0381; 

  strcpy(pt[91].el,"PA"); strcpy(pt[92].el,"U ");  strcpy(pt[93].el,"NP");  strcpy(pt[94].el,"PU");  strcpy(pt[95].el,"AM"); 
  pt[91].am = 231.0359; pt[92].am = 238.029;  pt[93].am = 237.0482;  pt[94].am = 244;  pt[95].am = 243; 

  strcpy(pt[96].el,"CM"); strcpy(pt[97].el,"BK");  strcpy(pt[98].el,"CF");  strcpy(pt[99].el,"ES");  strcpy(pt[100].el,"FM"); 
  pt[96].am = 247; pt[97].am = 247;  pt[98].am = 251;  pt[99].am = 254;  pt[100].am = 257; 

 strcpy(pt[101].el,"MD"); strcpy(pt[102].el,"NO"); strcpy(pt[103].el,"LR"); strcpy(pt[104].el,"RF");  strcpy(pt[105].el,"DB"); 
  pt[101].am = 258; pt[102].am = 259;  pt[103].am = 262;  pt[104].am = 261.11;  pt[105].am = 	268; 

  strcpy(pt[106].el,"SG"); strcpy(pt[107].el,"BH");strcpy(pt[108].el,"HS");strcpy(pt[109].el,"MT");  strcpy(pt[110].el,"DS"); 
  pt[106].am = 271; pt[107].am = 270;  pt[108].am = 269;  pt[109].am = 278;  pt[110].am = 281; 

  strcpy(pt[111].el,"RG"); strcpy(pt[112].el,"CN"); strcpy(pt[113].el,"UT");strcpy(pt[114].el,"UQ");  strcpy(pt[115].el,"UP"); 
  pt[111].am = 281; pt[112].am = 285;  pt[113].am = 286 ;  pt[114].am = 289;  pt[115].am = 289; 

  strcpy(pt[116].el,"UH"); strcpy(pt[117].el,"US"); strcpy(pt[118].el,"UO");  // from 112 there should be 3 letters Uub ...
  pt[116].am = 293; pt[117].am = 294;  pt[118].am = 294; // from 112 there should be 3 letters: Uub ...
}

int atomicnumber(char q[3]) {  // v0.4
int  i, res ;
  res = 0;
  if (strlen(q)>0) {
      if (strlen(q) < 2) strcat(q," ");
//      q[2] = lowercase(q[2]);  // v0.4e
      for (i= 1;i<=max_atomicnum;i++) {
          if (!strcmp(q,pt[i].el)) res = i;
      }
  }
  if (!strcmp(q,"Q ") || !strcmp(q,"A ") || !strcmp(q,"X ")) res = 999;
return res;
}

char* mk_fragstr(FRAG_REC frag){
int  i ;
static char  res[max_fragstr_length], tmpstr[max_fragstr_length];
memset(res, 0, max_fragstr_length);
  for (i = 0; i<frag.size;i++) {
      sprintf(tmpstr,"%d",frag.a_atomicnum[i]);
      strcat(res,tmpstr);
      memset(tmpstr, 0, max_fragstr_length);
      if ((i < frag.size-1) || (frag.ring == true)) {
        tmpstr[0]=frag.b_code[i];
        tmpstr[1]='\0';
        strcat(res,tmpstr);
      }
  } //printf("gethere,%s\n",res);
  return res;  
}

char* assemble_frag(RINGPATH_TYPE fp, bool is_ring) {
int  i, a1, a2, b, ref, a_tmp;
short int  pl;
char  el1[3];
char fstr[max_fragstr_length];
char  bt, b_tmp, nbt[2];
FRAG_REC  lfrag;
bool  valid;
  valid = true;
//  memset(lfrag,0,sizeof(FRAG_REC));
  pl = path_pos(0,fp);
  a1 = fp[0];
  memset(fstr, 0, max_fragstr_length);
  strcpy(el1,atom[a1]->element);
  lfrag.size = 1;
  lfrag.a_atomicnum[0] = atomicnumber(el1);
  if (el1[1] == ' ') el1[1]='\0';
  strcat(fstr,el1);
  if (pl > 0) {
      lfrag.size = pl;
      lfrag.ring = is_ring;
      for (i =1;i<pl;i++) {
          a2 = a1;
          a1 = fp[i];
          b = get_bond(a1,a2);
          strcpy(el1,atom[a1]->element);
          bt = bond[b]->btype;
          if (bt == 'S') strcpy(nbt,"-");
          if (bt == 'D') strcpy(nbt,"=");
          if ((bt == 'A') || (bond[b]->arom == true)) strcpy(nbt,"~");  // v0.4b  switched lines with the next one
          if (bt == 'T') strcpy(nbt,"#");  // v0.4b  'T' remains 'T' even in aromatic rings
          if ((bt=='l') ||(bt=='s') || (bt=='d') || (bt=='a')) valid = false;  // reject fragments with wildcard bonds
          lfrag.a_atomicnum[i] = atomicnumber(el1);
          if (lfrag.a_atomicnum[i] > 900) valid = false;  // reject fragments with wildcard atoms
          lfrag.b_code[(i-1)] = nbt[0];
          if (el1[1] == ' ') el1[1]='\0';
          strcat(fstr,nbt);
          strcat(fstr,el1);
      }
      if (is_ring) {   // add the last bond (between first and last atom)
          a1 = fp[pl-1];
          a2 = fp[0];
          b = get_bond(a1,a2);
          bt = bond[b]->btype;
          if (bt == 'S') strcpy(nbt,"-");
          if (bt == 'D') strcpy(nbt,"=");
          if (bt == 'T') strcpy(nbt,"#");
          if ((bt == 'A') || (bond[b]->arom == true)) strcpy(nbt,"~");
          if ((bt=='l') ||(bt=='s') || (bt=='d') || (bt=='a')) valid = false;  // reject fragments with wildcard bonds
          lfrag.b_code[pl-1] = nbt[0];
          strcat(fstr,nbt);
      }//printf("gethere,%s\n",fstr);
  }
  if (valid) {
      if (is_ring == false) {//order frag
//order_frag(lfrag);
          ref = 0;
          for (i = 0;i<(lfrag.size/2);i++) {
            ref = ref + lfrag.a_atomicnum[i];
            ref = ref - lfrag.a_atomicnum[lfrag.size-1-i];
            if (ref != 0) break;
          }
          if (ref == 0) {
            for (i=0;i<(lfrag.size/2);i++) {
              ref = ref + (int)(lfrag.b_code[i]);
              ref = ref - (int)(lfrag.b_code[lfrag.size-2-i]);
              if (ref != 0) break;
            }
          }
          if (ref < 0) {//revert frag
            if (lfrag.ring == false) {
              for (i = 0;i<(lfrag.size/2);i++) {
                a_tmp = lfrag.a_atomicnum[i];
                b_tmp =  lfrag.b_code[i];
                lfrag.a_atomicnum[i] = lfrag.a_atomicnum[lfrag.size-1-i];
                lfrag.a_atomicnum[lfrag.size-1-i] = a_tmp;
                lfrag.b_code[i] = lfrag.b_code[lfrag.size-2-i];
                lfrag.b_code[lfrag.size-2-i] = b_tmp;
              }
            }
          }
      }
      return mk_fragstr(lfrag);
   } else { return "\0"; }
}

void mk_hfp(char fstr[max_fragstr_length]){
unsigned int  h1, h2, h3, h4;
int len;
  len=strlen(fstr);
  h1 = (FNVHash(fstr,len) % hfpsize) + 1;
  h2 = (DJBHash(fstr,len) % hfpsize) + 1;
  if (n_hfpbits > 2) h3 = (BKDRHash(fstr,len) % hfpsize) + 1;
  if (n_hfpbits > 3) h4 = (DEKHash(fstr,len) %hfpsize) + 1;
  hfp[h1] = true;
  hfp[h2] = true;
  if (n_hfpbits > 2) hfp[h3] = true;
  if (n_hfpbits > 3) hfp[h4] = true;
//  printf("gethere: %d, %d, %s\n",h1,h2,fstr);
}

void collect_frags(RINGPATH_TYPE fp) {  // recursive procedure
RINGPATH_TYPE  lfp ;
int  nb[max_neighbors];
short int  pl ;
short int  n_nb ;
int  i, j, a_last, nb_count,nb_next_count, k ;
char  thisfrag[max_fragstr_length] ;
  pl = path_pos(0,fp);
  memset(lfp,0, sizeof(RINGPATH_TYPE));
  memset(nb,0,max_neighbors*sizeof(int));
  memset(thisfrag,0,max_fragstr_length);
    for(i=0;i<pl;i++) {lfp[i]=fp[i];}
  a_last = lfp[pl-1];
  if (pl > max_fragpath_length) return;
  if (atom[a_last]->heavy == false) return;
  // check if the last atom is a query atom and stop if yes
  if (!strcmp(atom[a_last]->element,"Q ") || !strcmp(atom[a_last]->element,"A ") ||
     !strcmp(atom[a_last]->element,"X ")) {
      lfp[pl-1] = 0;  // remove last atom (this could be the only one...)
      return;
  }
  // check if fragment forms a ring; if yes, discard it (rings are treated separately
  if (pl > 2) {
      for (i=0;i<pl-1;i++) { 
          if (lfp[pl-1] == lfp[i]) return;
      }
  }
  if ((pl >= min_fragpath_length) && (pl <= max_fragpath_length)) {
      //assemble the fragment and put it into the list
      strcpy(thisfrag,assemble_frag(lfp,false));
      if (strlen(thisfrag)!=0) mk_hfp(thisfrag);
  }
  if (pl == 1) {
  memset(nb,0,max_neighbors*sizeof(int));
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == a_last) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== a_last)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  } else {
      memset(nb, 0, max_neighbors * sizeof(int));
      nb_next_count= 0;
      for (k=0; k<n_bonds;k++) {
        if ((bond[k]->a1 == a_last) && (bond[k]->a2 != lfp[(pl-2)]) &&
           (nb_next_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
          nb[nb_next_count] = bond[k]->a2;
          nb_next_count++;
        }
        if ((bond[k]->a2 == a_last) && (bond[k]->a1 != lfp[(pl-2)]) &&
           (nb_next_count < max_neighbors) && (atom[bond[k]->a1]->heavy)) {
          nb[nb_next_count] = bond[k]->a1;
          nb_next_count++;
        }
      } 
  }
  // any other case: fragment path is not yet complete
  if (pl == 1) {n_nb = atom[a_last]->neighbor_count;} else {n_nb = atom[a_last]->neighbor_count-1;}
  if (n_nb > 0) { 
      for (j=0;j<n_nb;j++) {         
          lfp[pl] = nb[j];
          collect_frags(lfp);
      }
  }
//  free(lfp); free(nb);
}  // collect_frags

void make_linear_frags() {
int  a;
RINGPATH_TYPE  lfp;
  memset(lfp,0,sizeof(RINGPATH_TYPE));
  if ((n_atoms < 2) || (n_bonds < 1)) return;
  for (a = 1;a<=n_atoms;a++) { 
      lfp[0] = a;
      collect_frags(lfp);
   }
}

void make_ring_frags() {
int  i, j,k ,a ;
RINGPATH_TYPE  lfp ;
long int  h1, h2, h3, h4 ;
char  fstr[max_fragstr_length];
int  rs, len ;
  if (n_rings < 1) return;
  for (i=0;i<n_rings;i++) {
      rs = ringprop[i]->size;
      if ((rs <= max_ringfragpath_length) && (ringprop[i]->envelope == false)) {
          memset(lfp,0,sizeof(RINGPATH_TYPE));
          for (j=0;j<rs;j++) lfp[j] = *(*(ring+i)+j); 
          strcpy(fstr,assemble_frag(lfp,true));
          len=strlen(fstr);
          if (len==0) return;
          h1 = ((FNVHash(fstr,len)%hfpsize)) + 1;
          h2 = ((DJBHash(fstr,len)%hfpsize)) + 1;
          if (n_hfpbits > 2) h3 = ((BKDRHash(fstr,len)%hfpsize)) + 1;
          if (n_hfpbits > 3) h4 = ((DEKHash(fstr,len)%hfpsize)) + 1;
          //rotate the ring into every possible position
          for (j = 0;j<(rs-1);j++) {// other rs-1 positions. 
              a = lfp[rs-1]; //rotate the ring
              for (k = rs-1;k>0;k--) lfp[k] = lfp[(k-1)];
              lfp[0] = a;
              strcpy(fstr,assemble_frag(lfp,true));
              len=strlen(fstr);
              if (len==0) return;
              h1 = h1 + ((FNVHash(fstr,len)%hfpsize) + 1);
              h2 = h2 + ((DJBHash(fstr,len)%hfpsize) + 1);
              if (n_hfpbits > 2) h3 = h3 + ((BKDRHash(fstr,len)%hfpsize) + 1);
              if (n_hfpbits > 3) h4 = h4 + ((DEKHash(fstr,len)%hfpsize) + 1);
          }

          //now reverse the ring and repeat the whole procedure
          memset(lfp,0,sizeof(RINGPATH_TYPE));
          for (j=0;j<rs;j++) lfp[j] = *(*(ring+i)+j); 
          for (k=0;k<(rs/2);k++) {
            a = lfp[k];
            lfp[k]=lfp[(rs-1-k)];
            lfp[(rs-1-k)] = a;
          }//reverse_ring(lfp);
          strcpy(fstr,assemble_frag(lfp,true));
          len=strlen(fstr);
          if (len==0) return;
          h1 = h1 + ((FNVHash(fstr,len)%hfpsize) + 1);
          h2 = h2 + ((DJBHash(fstr,len)%hfpsize) + 1);
          if (n_hfpbits > 2) h3= h3 + ((BKDRHash(fstr,len)%hfpsize) + 1);
          if (n_hfpbits > 3) h4= h4 + ((DEKHash(fstr,len)%hfpsize) + 1);
          for (j =0;j<(rs-1);j++) {
              a = lfp[rs-1]; //rotate the ring
              for (k = rs-1;k>0;k--) lfp[k] = lfp[(k-1)];
              lfp[0] = a;
              strcpy(fstr,assemble_frag(lfp,true));
              len=strlen(fstr);
              if (len==0) return;
              h1 = h1 + ((FNVHash(fstr,len)%hfpsize) + 1);
              h2 = h2 + ((DJBHash(fstr,len)%hfpsize) + 1);
              if (n_hfpbits > 2) h3 = h3 + ((BKDRHash(fstr,len)%hfpsize) + 1);
              if (n_hfpbits > 3) h4 = h4 + ((DEKHash(fstr,len)%hfpsize) + 1);
          }

          //accumulated hash values have to be broken down again
          //writeln(' sum: ',h1,'  ',h2);
          h1 = (h1%hfpsize) + 1;
          h2 = (h2%hfpsize) + 1;
          if (n_hfpbits > 2) h3 = (h3%hfpsize) + 1;
          if (n_hfpbits > 3) h4 = (h4%hfpsize) + 1;
          
//          printf("==cumulative ring hash for ring %d ========> %d, %d\n",i, h1, h2);

          hfp[h1] = true;
          hfp[h2] = true;
          if (n_hfpbits > 2) hfp[h3] = true;
          if (n_hfpbits > 3) hfp[h4] = true;
      }
  } 
}

/*void reorder_frag(FRAG_REC frag) {
// order branched fragments by atomic number + bond type
int  i, j, k, newpos;
long  score1, score2 ;
FRAG_REC  lfrag ;
  memset(lfrag,0,sizeof(FRAG_REC));
  lfrag.size = 0;
  lfrag.ring = frag.ring;
  for (i=0;i<frag.size;i++) do
    score1 = frag.a_atomicnum[i] * 10;
    switch(frag.b_code[i]){
      case '-' : score1++;        break;
      case '=' : score1=score1+2; break;
      case '#' : score1=score1+3; break;
      case '~' : score1=score1+4; break;
      default  : break;
    }
    lfrag.size++;
    if (lfrag.size == 1) {
          lfrag.a_atomicnum[0] = frag.a_atomicnum[0];
          lfrag.b_code[0] = frag.b_code[0];          
    } else {
          // sort any additional entry
      newpos = 0;
      for (j = 0; j<(lfrag.size -1);j++) {
        score2 = lfrag.a_atomicnum[j] * 10;
        switch(lfrag.b_code[j]){
          case '-' : score2++;        break;
          case '=' : score2=score2+2; break;
          case '#' : score2=score2+3; break;
          case '~' : score2=score2+4; break;
          default  : break;
        }
        if (score2 > score1) newpos++;
      }
      if (newpos < lfrag.size) {
        for (k = lfrag.size-1; k>newpos; k--) {
          lfrag.a_atomicnum[k] = lfrag.a_atomicnum[k-1];
          lfrag.b_code[k] = lfrag.b_code[k-1];
        }
      }
      lfrag.a_atomicnum[newpos] = frag.a_atomicnum[i];
      lfrag.b_code[newpos] = frag.b_code[i];
    }
  } // for i ...
  // now copy the ordered elements back to frag
  for (i = 0;i<frag.size;i++) {
      frag.a_atomicnum[i] = lfrag.a_atomicnum[i];
      frag.b_code[i] = lfrag.b_code[i];
  }
}*/

char* assemble_brfragstr(int center_atomicnum, FRAG_REC frag) {
int   i;
static char  fstr[max_fragstr_length], tmpstr[max_fragstr_length];
//  memset(fstr, 0, max_fragstr_length);
  sprintf(fstr,"%d",center_atomicnum);
  strcat(fstr,"(");
  for (i =0;i<frag.size;i++) {
      if (i > 0) strcat(fstr,",");
      memset(tmpstr, 0, max_fragstr_length);
      tmpstr[0]=frag.b_code[i];
      tmpstr[1]='\0';
      strcat(fstr,tmpstr);
      memset(tmpstr, 0, max_fragstr_length);
      sprintf(tmpstr,"%d",frag.a_atomicnum[i]);
      strcat(fstr,tmpstr);
  }
  strcat(fstr,")");
  return fstr;
}

void make_branched_frags() {
int   i, j, k, o1,o2,o3,o4, a1, a2, b, can;
NEIGHBOR_REC  nb;
int  nb_count,newpos;
long  score1, score2 ;
bool  valid;
FRAG_REC  lfrag, lfrag2, lfrag3;
char  bt, nbt;
char  brfragstr[max_fragstr_length];
  if (n_atoms > 3) {
      for (i =1;i<=n_atoms;i++) {
          a1 = i;
          if ((!strcmp(atom[i]->element,"C ") || !strcmp(atom[i]->element,"N ")) &&
            ((atom[i]->neighbor_count == 4) || (atom[i]->neighbor_count == 3))) {
              valid = true;
              memset(nb,0,sizeof(NEIGHBOR_REC));
              lfrag.size = atom[i]->neighbor_count;
              lfrag.a_atomicnum[i]=0;
//              lfrag.b_code[i]=' ';
              lfrag.ring = false;
  nb_count=0;
  for (k= 0;k< n_bonds;k++ ){
    if ((bond[k]->a1 == i) && (nb_count < max_neighbors) && (atom[bond[k]->a2]->heavy)) {
      nb[nb_count]=bond[k]->a2;
      nb_count++;
    }
    if ((bond[k]->a2== i)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }// nb = get_neighbors(i);
              can = atomicnumber(atom[i]->element);
              for (j =0;j<atom[i]->neighbor_count;j++) {
                  a2 = nb[j];
                  b = get_bond(a1,a2);
                  bt = ' ';
                  nbt = ' ';
                  if (b >= 0) {bt = bond[b]->btype;} else {valid = false;}
                  if (bt == 'S') nbt = '-';
                  if (bt == 'D') nbt = '=';
                  if ((bt == 'A') || (bond[b]->arom == true)) nbt = '~';
                  if (bt == 'T') nbt = '#';
                  if ((bt=='l') || (bt=='s') || (bt=='d') || (bt=='a')) valid = false;  // reject fragments with wildcard bonds
                  lfrag.a_atomicnum[j] = atomicnumber(atom[a2]->element);
                  if (lfrag.a_atomicnum[j] > 900) valid = false;  // reject fragments with wildcard atoms
                  if (!strcmp(atom[a2]->element,"A ") || !strcmp(atom[a2]->element,"Q ") ||
                   !strcmp(atom[a2]->element,"X ")) valid = false;  // reject fragments with wildcard atoms
                  lfrag.b_code[j] = nbt;
              }  // for j...
              if (valid) {
                  //brfragstr := assemble_brfragstr(can,lfrag);
//printf("  branched atom:%d  %s\n",i,brfragstr);
  lfrag2.size = 0;
  lfrag2.ring = lfrag.ring;
  for (o1=0;o1<lfrag.size;o1++) {
    score1 = lfrag.a_atomicnum[o1] * 10;
    switch(lfrag.b_code[o1]){
      case '-' : score1++;        break;
      case '=' : score1=score1+2; break;
      case '#' : score1=score1+3; break;
      case '~' : score1=score1+4; break;
      default  : break;
    }
    lfrag2.size++;
    if (lfrag2.size == 1) {
          lfrag2.a_atomicnum[0] = lfrag.a_atomicnum[0];
          lfrag2.b_code[0] = lfrag.b_code[0];          
    } else {
          // sort any additional entry
      newpos = 0;
      for (o2 = 0; o2<(lfrag2.size -1);o2++) {
        score2 = lfrag2.a_atomicnum[o2] * 10;
        switch(lfrag2.b_code[o2]){
          case '-' : score2++;        break;
          case '=' : score2=score2+2; break;
          case '#' : score2=score2+3; break;
          case '~' : score2=score2+4; break;
          default  : break;
        }
        if (score2 > score1) newpos++;
      }
      if (newpos < lfrag2.size) {
        for (o3 = lfrag2.size-1; o3>newpos; o3--) {
          lfrag2.a_atomicnum[o3] = lfrag2.a_atomicnum[o3-1];
          lfrag2.b_code[o3] = lfrag2.b_code[o3-1];
        }
      }
      lfrag2.a_atomicnum[newpos] = lfrag.a_atomicnum[o1];
      lfrag2.b_code[newpos] = lfrag.b_code[o1];
    }
  } // for i ...
  // now copy the ordered elements back to frag
  for (o4 = 0;o4<lfrag.size;o4++) {
      lfrag.a_atomicnum[o4] = lfrag2.a_atomicnum[o4];
      lfrag.b_code[o4] = lfrag2.b_code[o4];
  }                  //reorder_frag(lfrag);
                  strcpy(brfragstr,assemble_brfragstr(can,lfrag));
                  //writeln('        primary fragment after reordering: ',brfragstr);
                  mk_hfp(brfragstr);
                  // if this is a quaternary C or N, make all possible tertiary subgraphs
                  if (atom[i]->neighbor_count == 4) {
                      for (j=0;j<4;j++) { // outer loop determines which substituent NOT to include
//                          memset(lfrag3,0,sizeof(FRAG_REC));
                          lfrag3.size = 0;
                          lfrag3.ring = false;
                          for (k = 0;k<4;k++) {  // inner loop copies the remaining substituents from lfrag to lfrag3
                              if (j != k) {
                                  lfrag3.size++;
                                  lfrag3.a_atomicnum[lfrag3.size-1] = lfrag.a_atomicnum[k];
                                  lfrag3.b_code[lfrag3.size-1] = lfrag.b_code[k];
                              }
                          }
                          // and now the hash
  lfrag2.size = 0;
  lfrag2.ring = lfrag3.ring;
  for (o1=0;o1<lfrag3.size;o1++) {
    score1 = lfrag3.a_atomicnum[o1] * 10;
    switch(lfrag3.b_code[o1]){
      case '-' : score1++;        break;
      case '=' : score1=score1+2; break;
      case '#' : score1=score1+3; break;
      case '~' : score1=score1+4; break;
      default  : break;
    }
    lfrag2.size++;
    if (lfrag2.size == 1) {
          lfrag2.a_atomicnum[0] = lfrag3.a_atomicnum[0];
          lfrag2.b_code[0] = lfrag3.b_code[0];          
    } else {
          // sort any additional entry
      newpos = 0;
      for (o2 = 0; o2<(lfrag2.size -1);o2++) {
        score2 = lfrag2.a_atomicnum[o2] * 10;
        switch(lfrag2.b_code[o2]){
          case '-' : score2++;        break;
          case '=' : score2=score2+2; break;
          case '#' : score2=score2+3; break;
          case '~' : score2=score2+4; break;
          default  : break;
        }
        if (score2 > score1) newpos++;
      }
      if (newpos < lfrag2.size) {
        for (o3 = lfrag2.size-1; o3>newpos; o3--) {
          lfrag2.a_atomicnum[o3] = lfrag2.a_atomicnum[o3-1];
          lfrag2.b_code[o3] = lfrag2.b_code[o3-1];
        }
      }
      lfrag2.a_atomicnum[newpos] = lfrag3.a_atomicnum[o1];
      lfrag2.b_code[newpos] = lfrag3.b_code[o1];
    }
  } // for i ...
  // now copy the ordered elements back to frag
  for (o4 = 0;o4<lfrag3.size;o4++) {
      lfrag3.a_atomicnum[o4] = lfrag2.a_atomicnum[o4];
      lfrag3.b_code[o4] = lfrag2.b_code[o4];
  }//reorder_frag(lfrag3);
                          strcpy(brfragstr,assemble_brfragstr(can,lfrag3));
                          //writeln('        sub-branch after reordering: ',brfragstr);
                          mk_hfp(brfragstr);
                      }
                  }
              } // if valid
          }  // elem = C ...
       }  // for i...
   }
}

void write_hfp_dec() {
const int bsize = 32;
int  i, j, n1 ;
long int hfpincrement ;
long int hfpdecimal;
  n1 = 0;
  for (i = 0 ; i< ((hfpsize/bsize));i++) {
      hfpdecimal = 0;
      for (j =1;j<=bsize;j++){
          hfpincrement = 1;
          if (hfp[((bsize*i)+j)]) { 
              n1++;
              hfpincrement = hfpincrement<<(j-1);
              hfpdecimal = hfpdecimal + hfpincrement;
//          printf("test: %d, %d\n",j, hfpdecimal);
          }
      }
      if (i > 0) printf(",");
      printf("%ld",hfpdecimal);
  }
  printf(";%d",n1);
  printf("\n");
}

void make_hashed_fp() {
  init_pt();
  memset(hfp,0,hfpsize*sizeof(bool));
  // get the fragments
  make_linear_frags();
  make_ring_frags();
//  {$IFDEF extended_hfp}
  make_branched_frags();
//  {$ENDIF}
  // now write the output
//  if (hfpformat == fpf_boolean) write_hfp();
  if (hfpformat == fpf_decimal) write_hfp_dec(); 
}

void read_molfile() { //statistical analysis: step 1
int  n, v, tmp_n_atoms, tmp_n_bonds;  // v0.3l
char  rline[BUFSIZ];
char tmpstr[BUFSIZ];
char xstr[9];
char ystr[9];
char zstr[9];
char chgstr[4];
float  xval, yval, zval;
int chgval, chg_id, chg, chg_num;
char  a1str[4];
char a2str[4];
char elemstr[3], atypestr[4];
int  a1val, a2val;
int  ri,rc,bt,bs;
int  sepcount;
int  i;              // v0.3j
bool  clearcharges;   // v0.3j
  found_querymol = false;  // v0.3p
  clearcharges = true;     // v0.3j
//  if n_atoms > 0 then zap_molecule; //memeory clean to be tested later
//  rline='';
  ri = 1;
  molname = molbuf[ri];            // line 1
//  molname[strlen(molname) - 1] = '\0';
  memset(rline,0,BUFSIZ); 
  if (ri < molbufindex) {ri++;}  // line 2
  strcpy(rline,molbuf[ri]);
  if (strstr(rline,"CheckMol")||strstr(rline,"Tweaked")) //to test if it is a valid tweaked molfile 
    {
      found_arominfo = true;
      tmfcode = 1;  // v0.3m (begin)
      if ((strlen(rline) >= 39) && (strstr(rline,"TMF")))  // v0.3m; encoding of tweaklevel
        { 
          strncpy(tmpstr, rline+37, 2); //tweak levels
          tmfcode=atoi(tmpstr); //tmfcode = int(tmpstr)
          memset(tmpstr,0,BUFSIZ); 
        }
      if (tmfcode<0) {tmfmismatch = true; } else {tmfmismatch = false;}
      if (((strstr(rline,":r0")) && (ringsearch_mode != rs_sar)) ||
         ((strstr(rline,":r1")) && (ringsearch_mode != rs_ssr))) { tmfmismatch = true;}
      if (((strstr(rline,":m0")) && (opt_metalrings == true)) ||
         ((strstr(rline,":m1")) && (opt_metalrings == false))) { tmfmismatch = true;}
/*      {$IFDEF debug}
      if tmfmismatch then 
        debugoutput('"tweaked" molfile: version mismatch') else
        debugoutput('"tweaked" molfile: version OK');
      {$ENDIF}
      // v0.3m (end)*/
    }
  memset(rline,0,BUFSIZ); 
  if (ri < molbufindex) {ri++;}  // line 3
  strcpy(rline,molbuf[ri]);
  molcomment = molbuf[ri];
//  molcomment[strlen(molcomment) - 1] = '\0';
  memset(rline,0,BUFSIZ); 
  if (ri < molbufindex) {ri++;}  // line 4
  strcpy(rline,molbuf[ri]);
  memset(tmpstr,0,BUFSIZ);
  strncpy(tmpstr,rline,3); // take the first 3 chars. 
  n_atoms=atoi(tmpstr);
  memset(tmpstr,0,BUFSIZ); // must clear the pointer to tmpstr.
  strncpy(tmpstr,rline+3,3);
  n_bonds=atoi(tmpstr);
  memset(tmpstr,0,BUFSIZ); 
  strncpy(tmpstr,rline+9,3);   // if it is a CheckMol-tweaked molfile, this is the number of rings
  n_cmrings=atoi(tmpstr);
  memset(tmpstr,0,BUFSIZ); 
  if (n_cmrings<0) { n_cmrings = 0; } // In case of bad transfer of string to int.
  // do some range checking for n_atoms, n_bonds; new in v0.3l
  tmp_n_atoms = n_atoms;
  if (n_atoms > max_atoms) {n_atoms = max_atoms;}
  if (n_atoms < 0) {n_atoms = 0;}
  tmp_n_bonds = n_bonds;
  if (n_bonds > max_bonds) {n_bonds = max_bonds;}
  if (n_bonds < 0) {n_bonds = 0;}
  if ((n_atoms == 0) && opt_verbose) {  // v0.3l
      printf("WARNING: Possibly no struct read!");
      exit(0);
    }
  for(n=0;n<=n_atoms;n++) {
    atom[n]=(ATOM_REC *)malloc(sizeof(ATOM_REC));} //initialize atom records; atom[1-n], avoiding atom[0];
  if ( atom==NULL ) {printf("Error allocating memory for ATOM!");exit(1);}
  for(n=0;n<n_bonds;n++) {
    bond[n]=(BOND_REC *)malloc(sizeof(BOND_REC)); }     // this would be only one calloc() in C;  v0.3l
  if ( bond==NULL ) {printf("Error allocating memory for BOND!");exit(1);}
//    ring[n]=(RINGLIST *)malloc(sizeof(RINGLIST));
  for(n=0;n<max_rings;n++) {
    ringprop[n]=(RINGPROP_REC *)malloc(sizeof(RINGPROP_REC));}
  if ( ringprop==NULL ) {printf("Error allocating memory for RINGPROP!");exit(1);}
  for(n=0;n<max_rings;n++) ring[n]=(int *)malloc(max_ringsize*sizeof(int));
  if ( ring==NULL) {printf("Error allocating memory for ring!");exit(1);}
  // check for the chirality flag
  if ((strlen(rline) > 15) && (rline[14]=='1')) { chir_flag = true;}  // new in v0.3f
  memset(rline,0,BUFSIZ); 
//start reading the stereo structure
  n_heavyatoms = 0;
  n_trueheavyatoms = 0;  // v0.4b
  n_heavybonds = 0;
  n_Ctot       = 0;  // v0.3g
  n_Otot       = 0;  // v0.3g
  n_Ntot       = 0;  // v0.3g

  if (n_atoms > 0) {  // v0.3l
      for (n=1; n<=tmp_n_atoms; n++) {
          if (n <= max_atoms) {v = n;} else {v = max_atoms; } 
              atom[v]->x= 0; atom[v]->y= 0; atom[v]->z = 0;  // v0.3g
              atom[v]->formal_charge  = 0;
              atom[v]->real_charge    = 0;
              atom[v]->Hexp           = 0;
              atom[v]->Htot           = 0;
              atom[v]->neighbor_count = 0;
              atom[v]->ring_count     = 0;
              atom[v]->arom           = false;
              atom[v]->q_arom         = false;
              atom[v]->stereo_care    = false;
              atom[v]->metal          = false;
              atom[v]->heavy          = false;
              atom[v]->tag            = false;
              atom[v]->nucleon_number = 0;
              atom[v]->radical_type   = 0;
          ri++;
          strcpy(rline,molbuf[ri]);
          memset(elemstr,0,3);
          strncpy(elemstr,rline+31,2); 
          for(i=0;i<sizeof(elemstr);i++) elemstr[i] = toupper(elemstr[i]); 
          strncpy(atypestr,rline+31,3);
//          elemstr  := get_MDLelement(atomtype); str3 to str2!
          if (!strcmp(elemstr,"C ")) {n_Ctot++;}
          if (!strcmp(elemstr,"O ")) {n_Otot++;}
          if (!strcmp(elemstr,"N ")) {n_Ntot++;}
//          newatomtype = convert_MDLtype(atomtype); //??
          strncpy(xstr,rline,10);   // fixed in v0.3k (was: 2,9 etc.)
          strncpy(ystr,rline+10,10);
          strncpy(zstr,rline+20,10);
          memset(chgstr,0,3);
          strncpy(chgstr,rline+36,3);   // new in v0.3j
          chgstr[4]='\0';
          chgval=atoi(chgstr);
          if (chgval > 0) {
              if ((chgval >= 1) && (chgval <= 7)) {chgval=4-chgval;} else {chgval = 0;}
          }                      // end (v0.3j) charge: +3~-3
          xval=atof(xstr);
          yval=atof(ystr);
          zval=atof(zstr);  // v0.3k: removed superfluous val(chgstr,chgval,code)
//          mstat[0]->n_C1=10; 
              strcpy(atom[v]->element,elemstr);
//convert old atype to new atype
  if (!strcmp(atypestr,"O  ")) strcpy(atypestr,"O3 ");
  if (!strcmp(atypestr,"N  ")) strcpy(atypestr,"N3 ");
  if (!strcmp(atypestr,"Cl ")) strcpy(atypestr,"CL ");
  if (!strcmp(atypestr,"Br ")) strcpy(atypestr,"BR ");
  if (!strcmp(atypestr,"Al ")) strcpy(atypestr,"AL ");
  if (!strcmp(atypestr,"ANY")) strcpy(atypestr,"A ");
  if (!strcmp(atypestr,"Ca ")) strcpy(atypestr,"CA ");
  if (!strcmp(atypestr,"Du ")) strcpy(atypestr,"DU ");
  if (!strcmp(atypestr,"Li ")) strcpy(atypestr,"LI ");
  if (!strcmp(atypestr,"Na ")) strcpy(atypestr,"NA ");
  if (!strcmp(atypestr,"S  ")) strcpy(atypestr,"S3 ");
  if (!strcmp(atypestr,"Si ")) strcpy(atypestr,"SI ");
  if (!strcmp(atypestr,"P  ")) strcpy(atypestr,"P4 ");
              strcpy(atom[v]->atype, atypestr);
              if ((!strcmp(elemstr,"A ")) || (!strcmp(elemstr,"Q ")) || (!strcmp(elemstr,"X "))) {
                found_querymol = true;  // 'X ' added in v0.3n ??
                strcpy(atom[v]->atype,elemstr);
                strcat(atom[v]->atype," ");
              }
              if (!strcmp(elemstr,"D ")) {   // v0.3p
                  keep_DT = true;
                  strcpy(atom[v]->element,"H ");
                  atom[v]->nucleon_number = 2;
                  strcpy(atom[v]->atype,"DU ");
                }
              if (!strcmp(elemstr,"T ")) {   // v0.3p
                  keep_DT = true;
                  strcpy(atom[v]->element,"H ");
                  atom[v]->nucleon_number = 3; 
                  strcpy(atom[v]->atype,"DU ");
                }
              atom[v]->x = xval; atom[v]->y = yval; atom[v]->z = zval; 
              atom[v]->formal_charge = (int)(chgval+0.5); atom[v]->real_charge = 0; // v0.3j
              // read aromaticity flag from CheckMol-tweaked MDL molfile
              if ((strlen(rline) > 37) && (rline[37]=='0')) {
                  atom[v]->arom = true;
                  found_arominfo = true;
                }
              // new in v0.3d: read stereo care flag
              if ((strlen(rline) > 47) && (rline[47]=='1')) { atom[v]->stereo_care = true;}
              if (is_heavyatom(n)) { 
                  n_heavyatoms++;
                  atom[v]->heavy = true;
                  if (is_metal(n)) {atom[v]->metal = true;} else {atom[v]->metal=false;}
                  if (is_trueheavyatom(n)) { n_trueheavyatoms++; } // v0.4b
              atom[v]->nvalences = get_nvalences(elemstr);  // v0.3m                
            } 
          memset(rline,0,BUFSIZ); 
        }
    }  // if (n_atoms > 0)...
   if (n_bonds > 0) {  // v0.3l
      for (n=0;n<tmp_n_bonds;n++) {
          if (n <= max_bonds) {v = n;} else {v = max_bonds;}  // just for safety; v0.3l
          ri++;
          memset(a1str,0,4);
          memset(a2str,0,4);
          strcpy(rline,molbuf[ri]);
          strncpy(a1str,rline,3);
          strncpy(a2str,rline+3,3);
          a1val=atoi(a1str);
          if (a1val < 0) {a1val = 1; }  // v0.3l
          a2val=atoi(a2str);
          if (a2val < 0) {a2val = 1;}  // v0.3l
          bond[v]->a1 = a1val; bond[v]->a2 = a2val;
              if (rline[8] == '1') {bond[v]->btype = 'S';}  // single bond
              if (rline[8] == '2') {bond[v]->btype = 'D';}  // double bond
              if (rline[8] == '3') {bond[v]->btype = 'T';}  // triple bond
              if (rline[8] == '4') {bond[v]->btype = 'A';}  // aromatic bond
              if (rline[8] == '5') {bond[v]->btype = 'l';}  // single or double
              if (rline[8] == '6') {bond[v]->btype = 's';}  // single or aromatic
              if (rline[8] == '7') {bond[v]->btype = 'd';}  // double or aromatic
              if (rline[8] == '8') {bond[v]->btype = 'a';}  // any
              if ((bond[v]->btype =='l') || (bond[v]->btype=='s') || (bond[v]->btype=='d')
 || (bond[v]->btype=='a')) { found_querymol = true;}
                bond[v]->arom = false;
                bond[v]->q_arom = false;  // v0.3p
              // read aromaticity flag from CheckMol-tweaked MDL molfile
              if ((bond[v]->btype == 'A') || (rline[7] == '0')) {
                  bond[v]->arom = true;
                  if (rline[7] == '0') { found_arominfo = true;}
              }
              memset(tmpstr,0,BUFSIZ);
              strncpy(tmpstr,rline+12,3);  // new in v0.3d: read ring_count from tweaked molfile
              rc=atoi(tmpstr);
              if ((rc<0) || tmfmismatch) { 
                bond[v]->ring_count = 0;} else {bond[v]->ring_count = rc;}  // v0.3n: added tmfmismatch check
              memset(tmpstr,0,BUFSIZ);
              strncpy(tmpstr,rline+15,3);  // new in v0.3d: read bond topology;
              bt=atoi(tmpstr);         // extended features are encoded by leading zero
              if ((bt < 0) || (bt > 2)) {bond[v]->topo = btopo_any;} else   // v0.3n changed >5 into >2
                {
                  if (tmpstr[1] == '0') {bond[v]->topo = bt + 3;} else {bond[v]->topo = bt;}
                }
              // new in v0.3d: add stereo property from MDL "stereo care" flag in atom block
              bond[v]->stereo = bstereo_any;
              if (bond[v]->btype =='D') {
                if (atom[a1val]->stereo_care && atom[a2val]->stereo_care) { // this is the MDL-conformant encoding,
                      bond[v]->stereo = bstereo_xyz;   // for an alternative see below
                      ez_flag = true;         // v0.3f
                  } else {
                    memset(tmpstr,0,BUFSIZ);
                    strncpy(tmpstr,rline+9,3);  // new in v0.3d: read bond stereo specification;
                    bs=atoi(tmpstr);         // this extended feature is encoded by a leading zero
                    bond[v]->mdl_stereo = bs;            // v0.3n
                    if ((bs <= 0) || (bs > 2)) {
                      bond[v]->stereo = bstereo_any;
                    } else {
                      bond[v]->stereo = bstereo_xyz;
                    }
                    if (tmpstr[1] == '0') {bond[v]->stereo = bstereo_xyz;}
                  }
                }
              //if stereo <> bstereo_any then ez_search := true;
              if (bond[v]->stereo != bstereo_any) {ez_flag = true;}  // changed in v0.3f
              if ((bond[v]->btype =='S') && (strlen(rline)>11) && (rline[11]=='1')) 
              {bond[v]->stereo = bstereo_up;}
              if ((bond[v]->btype =='S') && (strlen(rline)>11) && (rline[11]=='6')) 
              {bond[v]->stereo = bstereo_down;}
              memset(tmpstr,0,BUFSIZ);
              strncpy(tmpstr,rline+9,3);  // new in v0.3n: save original bond stereo specification;
              bs=atoi(tmpstr);         // v0.3n
              bond[v]->mdl_stereo = bs;            // v0.3n
          }          //if atom^[a1val].heavy and atom^[a2val].heavy then inc(n_heavybonds);  // moved down; v0.4b
       }  // if (n_bonds > 0)...
  sepcount = 0;

  while ((ri < molbufindex) && (sepcount < 1)) {
      ri++;
      strcpy(rline,molbuf[ri]);
      if (strstr(rline,"M  CHG")) { 
                           // new in v0.3j
        if (clearcharges) {  // "M  CHG" supersedes all "old-style" charge values
              for (i = 0; i<=n_atoms; i++) {atom[i]->formal_charge = 0;}
          }
              memset(tmpstr,0,BUFSIZ);
              strncpy(tmpstr,rline+6,3);
              chg_num=atoi(tmpstr);  
        for(i=0;i<chg_num;i++) {
              memset(tmpstr,0,BUFSIZ);
              strncpy(tmpstr,rline+5+8*i+4,4);  
              chg_id=atoi(tmpstr); // atom array start by 0. 
              memset(tmpstr,0,BUFSIZ);
              strncpy(tmpstr,rline+5+8*i+8,4);  
              chg=atoi(tmpstr); // chg<0
              if((chg_id!=0) && (chg!=0)) {
              atom[chg_id]->formal_charge = chg;}
//         printf("FORMAL_CHARGE: %d, %d\n",chg_id, atom[chg_id]->formal_charge); 
         }
          clearcharges = false;  // subsequent "M  CHG" lines must not clear previous values
      }
      if (strstr(rline,"M  ISO")) {
              memset(tmpstr,0,BUFSIZ);
              strncpy(tmpstr,rline+6,3);  
              chg_num=atoi(tmpstr);
        for(i=0;i<chg_num;i++) {
              memset(tmpstr,0,BUFSIZ);
              strncpy(tmpstr,rline+5+8*i+4,4);  
              chg_id=atoi(tmpstr);
              memset(tmpstr,0,BUFSIZ);
              strncpy(tmpstr,rline+5+8*i+8,4);  
              chg=atoi(tmpstr);
          if((chg_id!=0) && (chg>0)) { 
              atom[chg_id]->nucleon_number = chg;
              if(!strcmp(atom[chg_id]->element,"H ") && (chg > 1)) {
                  keep_DT = false;
                  if(opt_iso) {
                      atom[chg_id]->heavy = true;
                      n_heavyatoms++;
                   }
                  strcpy(atom[chg_id]->atype,"DU ");
              }
           }
         }
       }
      if (strstr(rline,"M  RAD")) {
              memset(tmpstr,0,BUFSIZ);
              strncpy(tmpstr,rline+6,3);  
              chg_num=atoi(tmpstr);
        for(i=0;i<chg_num;i++) {
              memset(tmpstr,0,BUFSIZ);
              strncpy(tmpstr,rline+5+8*i+4,4);  
              chg_id=atoi(tmpstr);
              memset(tmpstr,0,BUFSIZ);
              strncpy(tmpstr,rline+5+8*i+8,4);  
              chg=atoi(tmpstr);
              if((chg_id!=0) && (chg!=0)) {
              atom[chg_id]->radical_type = chg;}
         }  // v0.3p
      }
      if (strstr(rline,"$$$$")) {
        sepcount++;
//        sep_label = get_sep_label(rline);   // v0.4c
        if (molbufindex > (ri + 2)) { mol_in_queue = true;}  // we assume this is an SDF file
        }
    }
  if (n_bonds > 0) {  // v0.4b  (must be done after "M  ISO" check)
      for (n = 0; n< n_bonds; n++) {
        a1val=bond[n]->a1;
        a2val=bond[n]->a2;
        if (atom[a1val]->heavy && atom[a2val]->heavy) {n_heavybonds++;}
        }
    }
//  fillchar(ring^,sizeof(ringlist),0); //initialize ring
  for (n = 0; n<n_cmrings; n++) {  // new in v0.3
      ringprop[n]->size     = 0;
      ringprop[n]->arom     = false;
      ringprop[n]->envelope = false;
  }
  li = ri + 1;

}

void write_MDLmolfile() {
int i;
int  a_chg, a_iso, a_rad;
  printf("%80s\n",molname);
  printf("  Tweaked                         TMF0%d",tweaklevel);  // v0.3m
  if (ringsearch_mode == rs_sar) printf(":r0");         // v0.3m
  if (ringsearch_mode == rs_ssr) printf(":r1");         // v0.3m
  if (opt_metalrings) {printf(":m1");} else {printf(":m0"); }// v0.3m
  printf("\n");
  printf("%s\n",molcomment);
  printf("%3d%3d%3d%3d%3d               999 V2000\n",n_atoms,n_bonds,0,n_rings,chir_flag);
  for (i= 1;i<=n_atoms;i++) {
    printf("%10.4f%10.4f%10.4f ",atom[i]->x,atom[i]->y,atom[i]->z);
    if (!strcmp(atom[i]->element,"H ") && (atom[i]->nucleon_number == 2) && keep_DT) { 
      printf("D "); } else 
    if (!strcmp(atom[i]->element,"H ") && (atom[i]->nucleon_number == 3) && keep_DT) {// v0.3p
      printf("T "); } else {
      printf("%s",atom[i]->element);
    }
    printf("  0"); //mass difference (isotopes) not used yet!
    if (atom[i]->arom) {printf(" 00"); } else { printf("  0");}
    printf("  0  0  0  0  0  0  0  0  0  0\n");
  }
  for (i=0;i<n_bonds;i++) {
    printf("%3d%3d",bond[i]->a1,bond[i]->a2);
    if (bond[i]->arom) {printf(" 0");} else {printf("  ");}
    if (bond[i]->btype == 'S') printf("1");
    if (bond[i]->btype == 'D') printf("2");
    if (bond[i]->btype == 'T') printf("3");
    if (bond[i]->btype == 'A') printf("4");
    if (bond[i]->btype == 'l') printf("5");
    if (bond[i]->btype == 's') printf("6");
    if (bond[i]->btype == 'd') printf("7");
    if (bond[i]->btype == 'a') printf("8");
    printf("  %d%3d  0  0\n",bond[i]->mdl_stereo,bond[i]->ring_count);
  }
  for (i=1;i<=n_atoms;i++) {
      a_chg= atom[i]->formal_charge;
      if (a_chg != 0) {
          printf("M  CHG  1 %3d %3d\n",i,a_chg);
      }
  }
  for (i= 1;i<=n_atoms;i++) {
      a_iso= atom[i]->nucleon_number;
      if (a_iso != 0) {
          if (strcmp(atom[i]->element,"H ") || (keep_DT == false)) {
              printf("M  ISO  1 %3d %3d\n",i,a_iso);
          }
      }
  }
  for (i= 1;i<=n_atoms;i++) {
      a_rad = atom[i]->radical_type;
      if (a_rad != 0) {
          printf("M  RAD  1 %3d %3d\n",i,a_rad);
      }
  }
  printf("M  END\n");
}

void list_molstat_codes() {
  printf("n_atoms:     number of heavy atoms\n");
  printf("n_bonds:     number of bonds between non-H atoms\n");
  printf("n_rings:     number of rings\n");
  printf("n_QA:        number of query atoms\n");
  printf("n_QB:        number of query bonds\n");
  printf("n_chg:       number of charges\n");
  printf("n_C1:        number of sp-hybridized carbon atoms\n");
  printf("n_C2:        number of sp2-hybridized carbon atoms\n");
  printf("n_C:         total number of carbon atoms\n");
  printf("n_CHB1p:     number of carbon atoms with at least 1 bond to a hetero atom\n");
  printf("n_CHB2p:     number of carbon atoms with at least 2 bonds to a hetero atom\n");
  printf("n_CHB3p:     number of carbon atoms with at least 3 bonds to a hetero atom\n");
  printf("n_CHB4:      number of carbon atoms with 4 bonds to a hetero atom\n");
  printf("n_O2:        number of sp2-hybridized oxygen atoms\n");
  printf("n_O3:        number of sp3-hybridized oxygen atoms\n");
  printf("n_N1:        number of sp-hybridized nitrogen atoms\n");
  printf("n_N2:        number of sp2-hybridized nitrogen atoms\n");
  printf("n_N3:        number of sp3-hybridized nitrogen atoms\n");
  printf("n_S:         number of sulfur atoms\n");
  printf("n_SeTe:      total number of selenium and tellurium atoms\n");
  printf("n_F:         number of fluorine atoms\n");
  printf("n_Cl:        number of chlorine atoms\n");
  printf("n_Br:        number of bromine atoms\n");
  printf("n_I:         number of iodine atoms\n");
  printf("n_P:         number of phosphorus atoms\n");
  printf("n_B:         number of boron atoms\n");
  printf("n_Met:       total number of metal atoms\n");
  printf("n_X:         total number of 'other' atoms (not listed above) and halogens\n");
  printf("n_b1:        number of single bonds\n");
  printf("n_b2:        number of double bonds\n");
  printf("n_b3:        number of triple bonds\n");
  printf("n_bar:       number of aromatic bonds\n");
  printf("n_C1O:       number of C-O single bonds\n");
  printf("n_C2O:       number of C=O double bonds\n");
  printf("n_CN:        number of C/N bonds (any type)\n");
  printf("n_XY:        number of heteroatom/heteroatom bonds (any type)\n");
  printf("n_r3:        number of 3-membered rings\n");
  printf("n_r4:        number of 4-membered rings\n");
  printf("n_r5:        number of 5-membered rings\n");
  printf("n_r6:        number of 6-membered rings\n");
  printf("n_r7:        number of 7-membered rings\n");
  printf("n_r8:        number of 8-membered rings\n");
  printf("n_r9:        number of 9-membered rings\n");
  printf("n_r10:       number of 10-membered rings\n");
  printf("n_r11:       number of 11-membered rings\n");
  printf("n_r12:       number of 12-membered rings\n");
  printf("n_r13p:      number of 13-membered or larger rings\n");
  printf("n_rN:        number of rings containing nitrogen (any number)\n");
  printf("n_rN1:       number of rings containing 1 nitrogen atom\n");
  printf("n_rN2:       number of rings containing 2 nitrogen atoms\n");
  printf("n_rN3p:      number of rings containing 3 or more nitrogen atoms\n");
  printf("n_rO:        number of rings containing oxygen (any number)\n");
  printf("n_rO1:       number of rings containing 1 oxygen atom\n");
  printf("n_rO2p:      number of rings containing 2 or more oxygen atoms\n");
  printf("n_rS:        number of rings containing sulfur (any number)\n");
  printf("n_rX:        number of heterocycles (any type)\n");
  printf("n_rar:       number of aromatic rings (any type)\n");
/*  {$IFDEF extended_molstat}
  printf("n_rbz:       number of benzene rings");
  printf("n_br2p:      number of bonds belonging to two or more rings");
  printf("n_psg01:     number of atoms belonging to group 1 of the periodic system");
  printf("n_psg02:     number of atoms belonging to group 2 of the periodic system");
  printf("n_psg13:     number of atoms belonging to group 13 of the periodic system");
  printf("n_psg14:     number of atoms belonging to group 14 of the periodic system");
  printf("n_psg15:     number of atoms belonging to group 15 of the periodic system");
  printf("n_psg16:     number of atoms belonging to group 16 of the periodic system");
  printf("n_psg17:     number of atoms belonging to group 17 of the periodic system");
  printf("n_psg18:     number of atoms belonging to group 18 of the periodic system");
  printf("n_pstm:      number of atoms belonging to the transition metals");
  printf("n_psla:      number of atoms belonging to the lanthanides or actinides");
  {$ENDIF} */
}

void show_usage() {
      printf("molstat 2014--for MOLBASE\n");
      printf("Usage: molstat [options] <filename>\n");
      printf("  where options can be:\n");
      printf("    -a  count charges in fingerprint\n");
      printf("    -b  bitstring (in decimal format) representing the presence of each group\n");
      printf("    -l  print a list of fingerprint codes + explanation and exit\n");
      printf("    -m  write MDL molfile (tweaked for aromatic atoms/bonds)\n");
      printf("    -r  force SSR (set of small rings) ring search mode\n");
      printf("    -t  list functional groups\n");
      printf("    -x  print molecular statistics (number of various atom types, bond types,\n");
      printf("    -v  extended molecular statistics\n");
      printf("    -s  standard in mode\n");
      printf("        ring sizes, etc.), 0s are omitted\n");
      printf("    -X  same as above, listing all records (even if 0) as comma-separated list\n");
      printf("    -H  hashed fingerprint mode with decimal output\n");
      printf("    -M  accept metal atoms as ring members\n");
    // the "debug" option (-D) remains undocumented
}

//=======================================================================================

int main (int argc, char** argv)
{
  int c;
      int digit_optind = 0;
  char *rfile;
init_globals();
  while (1)
  {
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    static struct option long_options[] =
    {
      {"add", 1, 0, 0},
      {"append", 0, 0, 0},
      {"delete", 1, 0, 0},
      {"verbose", 0, 0, 0},
      {"create", 1, 0, 'c'},
      {"file", 1, 0, 0},{0, 0, 0, 0},
      {0, 0, 0, 0},{0, 0, 0, 0},{0, 0, 0, 0},{0, 0, 0, 0},{0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "abhmrtlxvsMXH", long_options, &option_index);
    if (c == -1)
        break;

          switch (c)
            {
            case 0:
              printf ("option %s", long_options[option_index].name);
              if (optarg)
                printf (" with arg %s", optarg);
              printf ("\n");
              break;
            case 'a':
              opt_chg=true;
              break;
            case 'b':
              opt_bin=true;
              break;
            case 'h':
              opt_hfp= true;
              hfpformat= fpf_boolean;
              break;
            case 'm':
              opt_text=false;
              opt_bin= false;
              opt_bitstring= false;
              opt_code= false;
              opt_pos= false;  // v0.5
              opt_molstat= false;
              opt_hfp= false;
              opt_xmdlout= true;
              break;
            case 'r':
              opt_rs=rs_ssr;
              break;
            case 't':
              opt_verbose=true;
              opt_text=true;
//              fglang=1;
              break;
            case 'l':
              opt_list=true;
              break;
            case 'x':
              opt_molstat=true;
              break;
            case 's':
              opt_stdin=true;
              break;
            case 'v':
              opt_molstat=true;
              opt_molstat_v=true;
              break;
            case 'M':
              opt_metalrings=true;
              break;
            case 'X':
              opt_molstat =true;
              opt_molstat_X = true;
              break;
            case 'H':
              opt_hfp=true;
              hfpformat=fpf_decimal;
              break;
            case '?':
              break;
            default:
              printf ("?? getopt returned character code 0%o ??\n", c);
            }
        }
if(argc<2) {show_usage(); return 0; }
if(opt_list) {list_molstat_codes(); return 0;};
if(!opt_stdin) {
  if (strlen(argv[2])>30) { 
    rfile=argv[2]; 
    readbuffer(rfile);//read from "$mol" 
  } else {
    molfilename=argv[2];
    readfile();  
  }
} else {
  readfile(); 
}
if (ringsearch_mode == rs_sar) {max_vringsize = max_ringsize;} else {
                                   max_vringsize = ssr_vringsize;}
mol_OK=true;
read_molfile();
count_neighbors();
if ((!mol_OK) || (n_atoms < 1)){  // v0.3g; check if this is a valid query structure
  printf("Invalid molecule!\n");
  exit(EXIT_SUCCESS);
}
  chk_ringbonds();
  if (ringsearch_mode == rs_ssr) remove_redundant_rings();
      if (n_rings == max_rings) {
          ringsearch_mode = rs_ssr;
          auto_ssr = true;  // v0.3n
          clear_rings();
          max_vringsize = ssr_vringsize;   // v0.3n (was: 10)
          chk_ringbonds();
          remove_redundant_rings();
      }
  update_ringcount();
      // new in v0.3k: if output is a molfile, leave the original
      // representation of N-oxides, S-oxides, nitro groups, etc.
      // unchanged (ionic or non-ionic), in any other case make covalent bonds
  if (!opt_xmdlout) {normalize_ionic_bonds();}  // v0.3k
  update_atypes();
  update_Htotal();    // added in v0.3
  chk_arom();
  if (ringsearch_mode == rs_ssr) {  // new in v0.3
          while((prev_n_ar - n_ar) != 0){
            prev_n_ar = count_aromatic_rings();
            chk_arom();
            n_ar = count_aromatic_rings();
          }
  } 
/*} else { //matchmol routine
      if (opt_xmdlout==false) {normalize_ionic_bonds(); } // v0.3k
      update_atypes();  // added in v0.2f
      update_Htotal();  // end v0.2b snippet
}*/
if (opt_verbose) {write_mol();}
get_molstat();
if (opt_molstat || opt_hfp) {
          if (opt_molstat) { 
              if (opt_molstat_X) {printf("\r\r");write_molstat_X();} else {write_molstat();}
          }
} else {
          if (found_querymol) {
              printf("input structure contains query atom or query bond!\n");
              exit(0);
          }
}
      if (opt_none) { opt_text = true;}
      if (opt_text || opt_code || opt_bin || opt_bitstring || opt_pos) {// v0.5
          chk_functionalgroups();
          if (opt_text) write_fg_text();
//          if opt_text_de   then write_fg_text(lang_de);
//          if opt_code      then write_fg_code;
          if (opt_bin) write_fg_binary();
//          if (opt_bitstring) write_fg_bitstring;
//          if (opt_pos) write_fg_pos;  // v0.5
      }
      if (opt_xmdlout) {
          init_pt();
//          calc_mf_mw();  // v0.4d
          write_MDLmolfile();
      }
      if (opt_hfp) make_hashed_fp();
      zap_molecule();
return 0;
}
