#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
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

int fin;
FILE *infile;
int  hfpformat;
char* molfilename;
char* ndl_molfilename;
char* hay_molfilename;
char* molname;
char* ndl_molname;
char* tmp_molname;
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
bool  tmp_querymol;
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

typedef struct P_3D{
double x;
double y;
double z;
}P_3D;

typedef int NEIGHBOR_REC[max_neighbors];
typedef int RINGPATH_TYPE[max_ringsize];
typedef char FRAGSTR[max_fragstr_length];
typedef int MATCHPATH_TYPE[max_matchpath_length];
typedef bool MATCHMATRIX[max_ndl_gmmsize][max_matchpath_length];
typedef bool GLOBAL_MATCHMATRIX[max_ndl_gmmsize][max_hay_gmmsize];
typedef int CHIRPATH_TYPE[4];
//typedef RINGPATH_TYPE RINGLIST[max_rings]; //i

ATOM_REC *atom[max_atoms]; 
ATOM_REC *ndl_atom[max_atoms]; 
ATOM_REC *tmp_atom[max_atoms]; 
BOND_REC *bond[max_bonds];
BOND_REC *ndl_bond[max_bonds];
BOND_REC *tmp_bond[max_bonds];
RINGPROP_REC *ringprop[max_rings];
RINGPROP_REC *ndl_ringprop[max_rings];
RINGPROP_REC *tmp_ringprop[max_rings];
CONVAL_REC *cv[max_atoms];  // new in v0.3j
int *ring[max_rings];
int *ndl_ring[max_rings];
int *tmp_ring[max_rings];
int *fgloc[used_fg+1]; // fgloc[0] is empty.
bool fg[max_fg];
bool *gmm[max_ndl_gmmsize];
bool *gmm_total[max_ndl_gmmsize];
int ndl_ref_atom;
MOLSTAT_REC molstat;
MOLSTAT_REC ndl_molstat;
MOLSTAT_REC tmp_molstat;
PT_REC pt[max_atomicnum+1];
bool matchresult;
bool matchsummary;
MATCHPATH_TYPE ndl_matchpath;
MATCHPATH_TYPE hay_matchpath;

char molbuf[max_atoms+max_bonds+8192][BUFSIZ];
int molbufindex;
int li;  // the start line number of a molecule.
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
bool  found_untagged;
int  tmp_tag[max_atoms+1];
int n_matches;
bool  overall_match;
/*  match_string   : ansistring;  // v0.4b
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
  //ringsearch_mode := rs_sar;
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

void zap_needle() {
    if (ndl_atom != NULL) {
        free(*ndl_atom);
        *ndl_atom=NULL;
    }
    if (ndl_bond != NULL) {
        free(*ndl_bond);
        *ndl_bond=NULL;
    }
    if (ndl_ring != NULL) {
        free(*ndl_ring);
        *ndl_ring=NULL;
    }
    if (ndl_ringprop != NULL) {
        free(*ndl_ringprop);
        *ndl_ringprop=NULL;
    }
  ndl_n_atoms = 0;
  ndl_n_bonds = 0;
  ndl_n_rings = 0;
}

void zap_tmp() {
    if (tmp_atom != NULL) {
        free(*tmp_atom);
        *tmp_atom=NULL;
    }
    if (tmp_bond != NULL) {
        free(*tmp_bond);
        *tmp_bond=NULL;
    }
    if (tmp_ring != NULL) {
        free(*tmp_ring);
        *tmp_ring=NULL;
    }
    if (tmp_ringprop != NULL) {
        free(*tmp_ringprop);
        *tmp_ringprop=NULL;
    }
  tmp_n_atoms = 0;
  tmp_n_bonds = 0;
  tmp_n_rings = 0;
}

void readfile(char *filename) { 
char *inname = filename; 
int i;
static char *eom;
static char line_buffer[BUFSIZ]; /* BUFSIZ is defined if you include stdio.h */ 
static char numb[BUFSIZ];
static int line_number; 
if (!opt_stdin) {
  infile = fopen(inname, "r"); 
  if (!infile) { 
    printf("Couldn't open file %s for reading.\n", inname); 
    exit(0);
  } 
  line_number = 0; 
  mol_in_queue=false;
//skip fin lines of the haystack file!
  for (i=0;i<fin;i++) fgets(numb, sizeof(numb), infile);
  while (fgets(line_buffer, sizeof(line_buffer), infile)) { 
    if (line_number < (max_atoms+max_bonds+64)) {        
      ++line_number; /* note that the newline is in the buffer */ 
      for (i=strlen(line_buffer)-1;i>=0;i--){
        if(line_buffer[i]=='\n'||line_buffer[i]=='\r') {
          line_buffer[strlen(line_buffer)-1]='\0';} else {break;} 
      }//remove last '\n'
      strcpy(molbuf[line_number],line_buffer);
//    printf("%s: %d\n",molbuf[line_number],line_number);
    } else {
      printf("Not enough memories for molfile!");
      exit(0);
    }
//printf("%4d: %s", line_number, molbuf[line_number]); 
    eom=strstr(line_buffer,molend);
    if(eom){molbufindex=line_number;}
    if(strstr(line_buffer,"$$$$")) {mol_in_queue=true; fin=fin+line_number; break;}
  }
//printf("\nTotal number of lines = %d %d\n", fin, molbufindex); 
} else {
  line_number = 0; 
  mol_in_queue=false;
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
    eom=strstr(line_buffer,molend);
    if(eom){molbufindex=line_number;}
    if(strstr(line_buffer,"$$$$")) {mol_in_queue=true; fin=fin+line_number; break;}
  }
} 
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
  memset(rline,0,BUFSIZ); 
  if (ri < molbufindex) {ri++;}  // line 2
  strcpy(rline,molbuf[ri]);
  if (strstr(rline,"CheckMol") || strstr(rline,"Tweaked")) //to test if it is a valid tweaked molfile 
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
  molcomment = molbuf[ri];
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
                bond[v]->ring_count = 0;
              } else {
              bond[v]->ring_count = rc;}  // v0.3n: added tmfmismatch check
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
              if (bond[v]->stereo != bstereo_any) {ez_flag = true;   }  // changed in v0.3f
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

void count_neighbors() {// counts heavy-atom neighbors and explicit hydrogens
int i;
  if ((n_atoms < 1) || (n_bonds < 1)) {return;}
//initialize
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

void clear_atom_tags() {
int i;
  if (n_atoms > 0) {
    for (i= 1;i<=n_atoms;i++) atom[i]->tag = false;
  }
}

void set_atom_tags() {
int i;
  if (n_atoms > 0){
    for (i= 1;i<=n_atoms;i++) atom[i]->tag = true;
  }
}

void clear_ndl_atom_tags() {
int i;
  if (ndl_n_atoms > 0) {
    for (i= 1;i<=ndl_n_atoms;i++) ndl_atom[i]->tag = false;
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
  for (i=1;i<rs1;i++) {
    for (j =1;j<rs2;j++){
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
      if (is_newring(l_path)) { 
            if (n_rings < max_rings) {
              add_ring(l_path); 
            } else{
              return false;
            }
      } 
      return true;
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
              new_atom = true;
              for (j =1; j<pl; j++) {
                if (nb[i]==l_path[j]) {      // v0.3k
                  new_atom = false; // v0.3k
                  break;
                }
              }
              if ((new_atom) && (pl < max_vringsize) && (n_rings < max_rings)) {
                  l_path[pl] = nb[i];
                  if (pl < max_ringsize-1) {l_path[pl+1] = 0; } // just to be sure
                  recursion_level++;                             // v0.3p (begin)
                  if (recursion_level > max_recursion_depth) {
                      n_rings = max_rings;
                      return rp;
                  } 
                  if (is_ringpath(l_path)) {rp=true; }      
                                          // v0.3p (end)
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
int i;
  n_rings = 0;
  memset(ring, 0, max_rings*max_ringsize * sizeof(int));
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
  b_id = -1;
  if (n_bonds > 0) {
    for (i=0; i<n_bonds; i++) {
        if (((bond[i]->a1 == ba1) && (bond[i]->a2 == ba2)) ||
           ((bond[i]->a1 == ba2) && (bond[i]->a2 == ba1)))  b_id = i;
    }
  }
return b_id;
}

int get_ndl_bond(int ba1, int ba2){
int i, b_id;
  b_id = -1;
  if (ndl_n_bonds > 0) {
    for (i=0; i<ndl_n_bonds; i++) {
        if (((ndl_bond[i]->a1 == ba1) && (ndl_bond[i]->a2 == ba2)) ||
           ((ndl_bond[i]->a1 == ba2) && (ndl_bond[i]->a2 == ba1)))  b_id = i;
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
    if ((bond[k]->a2== id)&& (nb_count < max_neighbors)&& (atom[bond[k]->a1]->heavy)) {
      nb[nb_count]=bond[k]->a1;
      nb_count++;
    }
  }
  if ((!strcmp(atom[id]->element,"C ")) && (atom[id]->neighbor_count > 0)) {
      for (i=0; i<atom[id]->neighbor_count; i++) {
          if ((bond[get_bond(id,nb[i])]->btype == 'D') &&
             (!strcmp(atom[(nb[i])]->element,"O "))) {   // no N, amidines are different...
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

void update_atypes_quick() {  // v0.4b
int i;
  if (n_atoms < 1) return;
  for (i = 1;i<=n_atoms;i++) {
      if (!strcmp(atom[i]->element,"C ") && (atom[i]->arom)) strcpy(atom[i]->atype,"CAR");
      if (!strcmp(atom[i]->atype,"N2 ") && (atom[i]->arom)) strcpy(atom[i]->atype,"NAR");
      if (atom[i]->Hexp > atom[i]->Htot) atom[i]->Htot = atom[i]->Hexp;  
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
          if (!ko) {    // v0.3j; odd pi_count might be compensated by arom_pi_diff
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

void copy_mol_to_needle() {
int i, j, n;
  if (n_atoms == 0) return;
  for(n=0;n<=n_atoms;n++) {
    ndl_atom[n]=(ATOM_REC *)malloc(sizeof(ATOM_REC));} //initialize atom records; atom[1-n], avoiding atom[0];
  if ( ndl_atom==NULL ) {printf("Error allocating memory for ndl ATOM!");exit(1);}
  for(n=0;n<n_bonds;n++) {
    ndl_bond[n]=(BOND_REC *)malloc(sizeof(BOND_REC)); }     // this would be only one calloc() in C;  v0.3l
  if ( ndl_bond==NULL ) {printf("Error allocating memory for ndl BOND!");exit(1);}
//    ring[n]=(RINGLIST *)malloc(sizeof(RINGLIST));
  for(n=0;n<max_rings;n++) {
    ndl_ringprop[n]=(RINGPROP_REC *)malloc(sizeof(RINGPROP_REC));}
  if ( ndl_ringprop==NULL ) {printf("Error allocating memory for ndl RINGPROP!");exit(1);}
  for(n=0;n<max_rings;n++) ndl_ring[n]=(int *)malloc(max_ringsize*sizeof(int));
  if ( ndl_ring==NULL) {printf("Error allocating memory for ndl ring!");exit(1);}

  ndl_n_atoms = n_atoms;
  ndl_n_bonds = n_bonds;
  ndl_n_rings = n_rings;
  ndl_n_heavyatoms = n_heavyatoms;
  ndl_n_trueheavyatoms = n_trueheavyatoms;  // v0.4b
  ndl_n_heavybonds = n_heavybonds;
  ndl_molname = molname;
  ndl_n_Ctot = n_Ctot;
  ndl_n_Otot = n_Otot;
  ndl_n_Ntot = n_Ntot;
  for (i = 1;i<=n_atoms;i++) {
      strcpy(ndl_atom[i]->element,atom[i]->element);
      strcpy(ndl_atom[i]->atype,atom[i]->atype);
      ndl_atom[i]->x              = atom[i]->x;
      ndl_atom[i]->y              = atom[i]->y;
      ndl_atom[i]->z              = atom[i]->z;
      ndl_atom[i]->formal_charge  = atom[i]->formal_charge;
      ndl_atom[i]->real_charge    = atom[i]->real_charge;
      ndl_atom[i]->Hexp           = atom[i]->Hexp;
      ndl_atom[i]->Htot           = atom[i]->Htot;
      ndl_atom[i]->neighbor_count = atom[i]->neighbor_count;
      ndl_atom[i]->ring_count     = atom[i]->ring_count;
      ndl_atom[i]->arom           = atom[i]->arom;
      ndl_atom[i]->stereo_care    = atom[i]->stereo_care;
      ndl_atom[i]->heavy          = atom[i]->heavy;  // v0.3l
      ndl_atom[i]->metal          = atom[i]->metal;  // v0.3l
      ndl_atom[i]->tag            = atom[i]->tag;  // v0.3o
      ndl_atom[i]->nucleon_number = atom[i]->nucleon_number;  // v0.3p
      ndl_atom[i]->radical_type   = atom[i]->radical_type;  // v0.3p
  }
  if (n_bonds > 0) {
      for (i =0;i<n_bonds;i++){
          ndl_bond[i]->a1         = bond[i]->a1;
          ndl_bond[i]->a2         = bond[i]->a2;
          ndl_bond[i]->btype      = bond[i]->btype;
          ndl_bond[i]->arom       = bond[i]->arom;
          ndl_bond[i]->ring_count = bond[i]->ring_count;  // new in v0.3d
          ndl_bond[i]->topo       = bond[i]->topo;        // new in v0.3d
          ndl_bond[i]->stereo     = bond[i]->stereo;      // new in v0.3d
          ndl_bond[i]->q_arom     = bond[i]->q_arom; 
      }
  }
  if (n_rings > 0) {
      for (i=0;i<n_rings;i++){
          for (j=0;j<max_ringsize;j++) *(*(ndl_ring+i)+j) = *(*(ring+i)+j);
      }
      for (i=0;i<max_rings;i++){   // new in v0.3
          ndl_ringprop[i]->size     = ringprop[i]->size;
          ndl_ringprop[i]->arom     = ringprop[i]->arom;
          ndl_ringprop[i]->envelope = ringprop[i]->envelope;
      }
  }
  ndl_molstat.n_QA = molstat.n_QA;       ndl_molstat.n_QB = molstat.n_QB; 
  ndl_molstat.n_chg = molstat.n_chg;     ndl_molstat.n_C1 = molstat.n_C1; 
  ndl_molstat.n_C2 = molstat.n_C2;       ndl_molstat.n_C = molstat.n_C;
  ndl_molstat.n_CHB1p = molstat.n_CHB1p; ndl_molstat.n_CHB2p = molstat.n_CHB2p;
  ndl_molstat.n_CHB3p = molstat.n_CHB3p; ndl_molstat.n_CHB4 = molstat.n_CHB4;
  ndl_molstat.n_O2 = molstat.n_O2;       ndl_molstat.n_O3  = molstat.n_O3;
  ndl_molstat.n_N1 = molstat.n_N1;       ndl_molstat.n_N2 = molstat.n_N2; 
  ndl_molstat.n_N3 = molstat.n_N3;       ndl_molstat.n_S = molstat.n_S; 
  ndl_molstat.n_SeTe = molstat.n_SeTe;   ndl_molstat.n_F = molstat.n_F; 
  ndl_molstat.n_Cl = molstat.n_Cl;       ndl_molstat.n_Br = molstat.n_Br;
  ndl_molstat.n_I = molstat.n_I;         ndl_molstat.n_P = molstat.n_P; 
  ndl_molstat.n_B = molstat.n_B;         ndl_molstat.n_Met = molstat.n_Met; 
  ndl_molstat.n_X = molstat.n_X;         ndl_molstat.n_b1 = molstat.n_b1; 
  ndl_molstat.n_b2 = molstat.n_b2;       ndl_molstat.n_b3 = molstat.n_b3; 
  ndl_molstat.n_bar = molstat.n_bar;     ndl_molstat.n_C1O = molstat.n_C1O; 
  ndl_molstat.n_C2O = molstat.n_C2O;     ndl_molstat.n_CN = molstat.n_CN; 
  ndl_molstat.n_XY = molstat.n_XY;       ndl_molstat.n_r3 = molstat.n_r3; 
  ndl_molstat.n_r4 = molstat.n_r4;       ndl_molstat.n_r5 = molstat.n_r5; 
  ndl_molstat.n_r6 = molstat.n_r6;       ndl_molstat.n_r7 = molstat.n_r7; 
  ndl_molstat.n_r8 = molstat.n_r8;       ndl_molstat.n_r9 = molstat.n_r9; 
  ndl_molstat.n_r10 = molstat.n_r10;     ndl_molstat.n_r11 = molstat.n_r11; 
  ndl_molstat.n_r12 = molstat.n_r12;     ndl_molstat.n_r13p = molstat.n_r13p;
  ndl_molstat.n_rN = molstat.n_rN;       ndl_molstat.n_rN1 = molstat.n_rN1; 
  ndl_molstat.n_rN2 = molstat.n_rN2;     ndl_molstat.n_rN3p = molstat.n_rN3p;
  ndl_molstat.n_rO = molstat.n_rO;       ndl_molstat.n_rO1 = molstat.n_rO1; 
  ndl_molstat.n_rO2p = molstat.n_rO2p;   ndl_molstat.n_rS = molstat.n_rS; 
  ndl_molstat.n_rX = molstat.n_rX;       ndl_molstat.n_rar = molstat.n_rar; 
  ndl_molstat.n_rBz = molstat.n_rBz;     ndl_molstat.n_br2p = molstat.n_br2p;  // v0.3n
  // make sure some modes can be switched on only by the query file
  // and not by subsequent haystack file(s)
  if (ez_flag) ez_search = true;    // new in v0.3f
  if (chir_flag) rs_search = true;    // new in v0.3f
  ndl_querymol = found_querymol;  // v0.3p
}

void copy_mol_to_tmp() {
int i, j, n;
  if (n_atoms == 0) return;
  for(n=0;n<=n_atoms;n++) {
    tmp_atom[n]=(ATOM_REC *)malloc(sizeof(ATOM_REC));} //initialize atom records; atom[1-n], avoiding atom[0];
  if ( tmp_atom==NULL ) {printf("Error allocating memory for tmp ATOM!");exit(1);}
  for(n=0;n<n_bonds;n++) {
    tmp_bond[n]=(BOND_REC *)malloc(sizeof(BOND_REC)); }     // this would be only one calloc() in C;  v0.3l
  if ( tmp_bond==NULL ) {printf("Error allocating memory for tmp BOND!");exit(1);}
//    ring[n]=(RINGLIST *)malloc(sizeof(RINGLIST));
  for(n=0;n<max_rings;n++) {
    tmp_ringprop[n]=(RINGPROP_REC *)malloc(sizeof(RINGPROP_REC));}
  if ( tmp_ringprop==NULL ) {printf("Error allocating memory for tmp RINGPROP!");exit(1);}
  for(n=0;n<max_rings;n++) tmp_ring[n]=(int *)malloc(max_ringsize*sizeof(int));
  if ( tmp_ring==NULL) {printf("Error allocating memory for tmp ring!");exit(1);}

  tmp_n_atoms = n_atoms;
  tmp_n_bonds = n_bonds;
  tmp_n_rings = n_rings;
  tmp_n_heavyatoms = n_heavyatoms;
  tmp_n_trueheavyatoms = n_trueheavyatoms;  // v0.4b
  tmp_n_heavybonds = n_heavybonds;
  tmp_molname = molname;
  tmp_n_Ctot = n_Ctot;
  tmp_n_Otot = n_Otot;
  tmp_n_Ntot = n_Ntot;
  for (i = 1;i<=n_atoms;i++) {
      strcpy(tmp_atom[i]->element,atom[i]->element);
      strcpy(tmp_atom[i]->atype,atom[i]->atype);
      tmp_atom[i]->x              = atom[i]->x;
      tmp_atom[i]->y              = atom[i]->y;
      tmp_atom[i]->z              = atom[i]->z;
      tmp_atom[i]->formal_charge  = atom[i]->formal_charge;
      tmp_atom[i]->real_charge    = atom[i]->real_charge;
      tmp_atom[i]->Hexp           = atom[i]->Hexp;
      tmp_atom[i]->Htot           = atom[i]->Htot;
      tmp_atom[i]->neighbor_count = atom[i]->neighbor_count;
      tmp_atom[i]->ring_count     = atom[i]->ring_count;
      tmp_atom[i]->arom           = atom[i]->arom;
      tmp_atom[i]->stereo_care    = atom[i]->stereo_care;
      tmp_atom[i]->heavy          = atom[i]->heavy;  // v0.3l
      tmp_atom[i]->metal          = atom[i]->metal;  // v0.3l
      tmp_atom[i]->tag            = atom[i]->tag;  // v0.3o
      tmp_atom[i]->nucleon_number = atom[i]->nucleon_number;  // v0.3p
      tmp_atom[i]->radical_type   = atom[i]->radical_type;  // v0.3p
  }
  if (n_bonds > 0) {
      for (i =0;i<n_bonds;i++){
          tmp_bond[i]->a1         = bond[i]->a1;
          tmp_bond[i]->a2         = bond[i]->a2;
          tmp_bond[i]->btype      = bond[i]->btype;
          tmp_bond[i]->arom       = bond[i]->arom;
          tmp_bond[i]->ring_count = bond[i]->ring_count;  // new in v0.3d
          tmp_bond[i]->topo       = bond[i]->topo;        // new in v0.3d
          tmp_bond[i]->stereo     = bond[i]->stereo;      // new in v0.3d
          tmp_bond[i]->q_arom     = bond[i]->q_arom; 
      }
  }
  if (n_rings > 0) {
      for (i=0;i<n_rings;i++){
          for (j=0;j<max_ringsize;j++) *(*(tmp_ring+i)+j) = *(*(ring+i)+j);
      }
      for (i=0;i<max_rings;i++){   // new in v0.3
          tmp_ringprop[i]->size     = ringprop[i]->size;
          tmp_ringprop[i]->arom     = ringprop[i]->arom;
          tmp_ringprop[i]->envelope = ringprop[i]->envelope;
      }
  }
  tmp_molstat.n_QA = molstat.n_QA;       tmp_molstat.n_QB = molstat.n_QB; 
  tmp_molstat.n_chg = molstat.n_chg;     tmp_molstat.n_C1 = molstat.n_C1; 
  tmp_molstat.n_C2 = molstat.n_C2;       tmp_molstat.n_C = molstat.n_C;
  tmp_molstat.n_CHB1p = molstat.n_CHB1p; tmp_molstat.n_CHB2p = molstat.n_CHB2p;
  tmp_molstat.n_CHB3p = molstat.n_CHB3p; tmp_molstat.n_CHB4 = molstat.n_CHB4;
  tmp_molstat.n_O2 = molstat.n_O2;       tmp_molstat.n_O3  = molstat.n_O3;
  tmp_molstat.n_N1 = molstat.n_N1;       tmp_molstat.n_N2 = molstat.n_N2; 
  tmp_molstat.n_N3 = molstat.n_N3;       tmp_molstat.n_S = molstat.n_S; 
  tmp_molstat.n_SeTe = molstat.n_SeTe;   tmp_molstat.n_F = molstat.n_F; 
  tmp_molstat.n_Cl = molstat.n_Cl;       tmp_molstat.n_Br = molstat.n_Br;
  tmp_molstat.n_I = molstat.n_I;         tmp_molstat.n_P = molstat.n_P; 
  tmp_molstat.n_B = molstat.n_B;         tmp_molstat.n_Met = molstat.n_Met; 
  tmp_molstat.n_X = molstat.n_X;         tmp_molstat.n_b1 = molstat.n_b1; 
  tmp_molstat.n_b2 = molstat.n_b2;       tmp_molstat.n_b3 = molstat.n_b3; 
  tmp_molstat.n_bar = molstat.n_bar;     tmp_molstat.n_C1O = molstat.n_C1O; 
  tmp_molstat.n_C2O = molstat.n_C2O;     tmp_molstat.n_CN = molstat.n_CN; 
  tmp_molstat.n_XY = molstat.n_XY;       tmp_molstat.n_r3 = molstat.n_r3; 
  tmp_molstat.n_r4 = molstat.n_r4;       tmp_molstat.n_r5 = molstat.n_r5; 
  tmp_molstat.n_r6 = molstat.n_r6;       tmp_molstat.n_r7 = molstat.n_r7; 
  tmp_molstat.n_r8 = molstat.n_r8;       tmp_molstat.n_r9 = molstat.n_r9; 
  tmp_molstat.n_r10 = molstat.n_r10;     tmp_molstat.n_r11 = molstat.n_r11; 
  tmp_molstat.n_r12 = molstat.n_r12;     tmp_molstat.n_r13p = molstat.n_r13p;
  tmp_molstat.n_rN = molstat.n_rN;       tmp_molstat.n_rN1 = molstat.n_rN1; 
  tmp_molstat.n_rN2 = molstat.n_rN2;     tmp_molstat.n_rN3p = molstat.n_rN3p;
  tmp_molstat.n_rO = molstat.n_rO;       tmp_molstat.n_rO1 = molstat.n_rO1; 
  tmp_molstat.n_rO2p = molstat.n_rO2p;   tmp_molstat.n_rS = molstat.n_rS; 
  tmp_molstat.n_rX = molstat.n_rX;       tmp_molstat.n_rar = molstat.n_rar; 
  tmp_molstat.n_rBz = molstat.n_rBz;     tmp_molstat.n_br2p = molstat.n_br2p; 
// v0.3n
  // make sure some modes can be switched on only by the query file
  // and not by subsequent haystack file(s)
  if (ez_flag) ez_search = true;    // new in v0.3f
  if (chir_flag) rs_search = true;    // new in v0.3f
}

void copy_tmp_to_mol() {
int i, j, n;
  if (tmp_n_atoms == 0) return;
  for(n=0;n<=tmp_n_atoms;n++) {
    atom[n]=(ATOM_REC *)malloc(sizeof(ATOM_REC));} //initialize atom records; atom[1-n], avoiding atom[0];
  if ( atom==NULL ) {printf("Error allocating memory for tmp ATOM!");exit(1);}
  for(n=0;n<tmp_n_bonds;n++) {
    bond[n]=(BOND_REC *)malloc(sizeof(BOND_REC)); }     // this would be only one calloc() in C;  v0.3l
  if ( bond==NULL ) {printf("Error allocating memory for tmp BOND!");exit(1);}
//    ring[n]=(RINGLIST *)malloc(sizeof(RINGLIST));
  for(n=0;n<max_rings;n++) {
    ringprop[n]=(RINGPROP_REC *)malloc(sizeof(RINGPROP_REC));}
  if ( ringprop==NULL ) {printf("Error allocating memory for tmp RINGPROP!");exit(1);}
  for(n=0;n<max_rings;n++) ring[n]=(int *)malloc(max_ringsize*sizeof(int));
  if ( ring==NULL) {printf("Error allocating memory for tmp ring!");exit(1);}

  n_atoms = tmp_n_atoms;
  n_bonds = tmp_n_bonds;
  n_rings = tmp_n_rings;
  n_heavyatoms = tmp_n_heavyatoms;
  n_trueheavyatoms = tmp_n_trueheavyatoms;  // v0.4b
  n_heavybonds = tmp_n_heavybonds;
  molname = tmp_molname;
  n_Ctot = tmp_n_Ctot;
  n_Otot = tmp_n_Otot;
  n_Ntot = tmp_n_Ntot;
  for (i = 1;i<=tmp_n_atoms;i++) {
      strcpy(atom[i]->element,tmp_atom[i]->element);
      strcpy(atom[i]->atype,tmp_atom[i]->atype);
      atom[i]->x              = tmp_atom[i]->x;
      atom[i]->y              = tmp_atom[i]->y;
      atom[i]->z              = tmp_atom[i]->z;
      atom[i]->formal_charge  = tmp_atom[i]->formal_charge;
      atom[i]->real_charge    = tmp_atom[i]->real_charge;
      atom[i]->Hexp           = tmp_atom[i]->Hexp;
      atom[i]->Htot           = tmp_atom[i]->Htot;
      atom[i]->neighbor_count = tmp_atom[i]->neighbor_count;
      atom[i]->ring_count     = tmp_atom[i]->ring_count;
      atom[i]->arom           = tmp_atom[i]->arom;
      atom[i]->stereo_care    = tmp_atom[i]->stereo_care;
      atom[i]->heavy          = tmp_atom[i]->heavy;  // v0.3l
      atom[i]->metal          = tmp_atom[i]->metal;  // v0.3l
      atom[i]->tag            = tmp_atom[i]->tag;  // v0.3o
      atom[i]->nucleon_number = tmp_atom[i]->nucleon_number;  // v0.3p
      atom[i]->radical_type   = tmp_atom[i]->radical_type;  // v0.3p
  }
  if (tmp_n_bonds > 0) {
      for (i =0;i<tmp_n_bonds;i++){
          bond[i]->a1         = tmp_bond[i]->a1;
          bond[i]->a2         = tmp_bond[i]->a2;
          bond[i]->btype      = tmp_bond[i]->btype;
          bond[i]->arom       = tmp_bond[i]->arom;
          bond[i]->ring_count = tmp_bond[i]->ring_count;  // new in v0.3d
          bond[i]->topo       = tmp_bond[i]->topo;        // new in v0.3d
          bond[i]->stereo     = tmp_bond[i]->stereo;      // new in v0.3d
          bond[i]->q_arom     = tmp_bond[i]->q_arom;
      }
  }
  if (tmp_n_rings > 0) {
      for (i=0;i<tmp_n_rings;i++){
          for (j=0;j<max_ringsize;j++) *(*(ring+i)+j)=*(*(tmp_ring+i)+j) ;
      }
      for (i=0;i<max_rings;i++){   // new in v0.3
          ringprop[i]->size     = tmp_ringprop[i]->size;
          ringprop[i]->arom     = tmp_ringprop[i]->arom;
          ringprop[i]->envelope = tmp_ringprop[i]->envelope;
      }
  }
  molstat.n_QA = tmp_molstat.n_QA;       molstat.n_QB = tmp_molstat.n_QB; 
  molstat.n_chg = tmp_molstat.n_chg;     molstat.n_C1 = tmp_molstat.n_C1; 
  molstat.n_C2 = tmp_molstat.n_C2;       molstat.n_C = tmp_molstat.n_C;
  molstat.n_CHB1p = tmp_molstat.n_CHB1p; molstat.n_CHB2p = tmp_molstat.n_CHB2p;
  molstat.n_CHB3p = tmp_molstat.n_CHB3p; molstat.n_CHB4 = tmp_molstat.n_CHB4;
  molstat.n_O2 = tmp_molstat.n_O2;       molstat.n_O3  = tmp_molstat.n_O3;
  molstat.n_N1 = tmp_molstat.n_N1;       molstat.n_N2 = tmp_molstat.n_N2; 
  molstat.n_N3 = tmp_molstat.n_N3;       molstat.n_S = tmp_molstat.n_S; 
  molstat.n_SeTe = tmp_molstat.n_SeTe;   molstat.n_F = tmp_molstat.n_F; 
  molstat.n_Cl = tmp_molstat.n_Cl;       molstat.n_Br = tmp_molstat.n_Br;
  molstat.n_I = tmp_molstat.n_I;         molstat.n_P = tmp_molstat.n_P; 
  molstat.n_B = tmp_molstat.n_B;         molstat.n_Met = tmp_molstat.n_Met; 
  molstat.n_X = tmp_molstat.n_X;         molstat.n_b1 = tmp_molstat.n_b1; 
  molstat.n_b2 = tmp_molstat.n_b2;       molstat.n_b3 = tmp_molstat.n_b3; 
  molstat.n_bar = tmp_molstat.n_bar;     molstat.n_C1O = tmp_molstat.n_C1O; 
  molstat.n_C2O = tmp_molstat.n_C2O;     molstat.n_CN = tmp_molstat.n_CN; 
  molstat.n_XY = tmp_molstat.n_XY;       molstat.n_r3 = tmp_molstat.n_r3; 
  molstat.n_r4 = tmp_molstat.n_r4;       molstat.n_r5 = tmp_molstat.n_r5; 
  molstat.n_r6 = tmp_molstat.n_r6;       molstat.n_r7 = tmp_molstat.n_r7; 
  molstat.n_r8 = tmp_molstat.n_r8;       molstat.n_r9 = tmp_molstat.n_r9; 
  molstat.n_r10 = tmp_molstat.n_r10;     molstat.n_r11 = tmp_molstat.n_r11; 
  molstat.n_r12 = tmp_molstat.n_r12;     molstat.n_r13p = tmp_molstat.n_r13p;
  molstat.n_rN = tmp_molstat.n_rN;       molstat.n_rN1 = tmp_molstat.n_rN1; 
  molstat.n_rN2 = tmp_molstat.n_rN2;     molstat.n_rN3p = tmp_molstat.n_rN3p;
  molstat.n_rO = tmp_molstat.n_rO;       molstat.n_rO1 = tmp_molstat.n_rO1; 
  molstat.n_rO2p = tmp_molstat.n_rO2p;   molstat.n_rS = tmp_molstat.n_rS; 
  molstat.n_rX = tmp_molstat.n_rX;       molstat.n_rar = tmp_molstat.n_rar; 
  molstat.n_rBz = tmp_molstat.n_rBz;     molstat.n_br2p = tmp_molstat.n_br2p; 
// v0.3n
  // make sure some modes can be switched on only by the query file
  // and not by subsequent haystack file(s)
  if (ez_flag) ez_search = true;    // new in v0.3f
  if (chir_flag) rs_search = true;    // new in v0.3f
}

void chk_wildcard_rings() {  // new in v0.3p
// checks if there are any wildcard atom types or bond types
// in a ring of the needle; if yes ==> set the q_arom flag in the
// atom and bond record of all ring members in order to perform the 
// match a bit more generously
int  i, j, rs;
int  a1, a2, b ;
bool  wcr;
  if (ndl_querymol == false) return;
  if (ndl_n_rings == 0) return;
  // now look for any not-yet-aromatic rings which contain a wildcard
  for (i=0;i<ndl_n_rings;i++) {
      wcr = false;
      if (ndl_ringprop[i]->arom == false) {
          rs = ndl_ringprop[i]->size;
          a2 = *(*(ndl_ring+i)+rs-1);
          for (j =0;j<rs;j++) {
            a1 = *(*(ndl_ring+i)+j);
            b = get_ndl_bond(a1,a2);
            if (!strcmp(ndl_atom[a1]->atype,"A  ") || !strcmp(ndl_atom[a1]->atype,"Q  ")) wcr= true;
            if ((ndl_bond[b]->btype == 'l') || (ndl_bond[b]->btype == 's') ||
              (ndl_bond[b]->btype == 'd') || (ndl_bond[b]->btype == 'a')) wcr = true;
              a2 = a1;   
          }
          if (wcr) { // if yes, flag all atoms and bonds in this ring as "potentially" aromatic
              a2 = *(*(ndl_ring+i)+rs-1);
              for (j=0;j<rs;j++) {
                  a1 = *(*(ndl_ring+i)+j);
                  b= get_ndl_bond(a1,a2);
                  ndl_atom[a1]->q_arom = true;
                  ndl_bond[b]->q_arom = true;
                  a2 = a1;   
              }
          }
      }
  }
  // and now undo this flagging for all rings which contain no wildcard
  for (i=0;i<ndl_n_rings;i++) {
    wcr = false;
    rs = ndl_ringprop[i]->size;
    a2 = *(*(ndl_ring+i)+rs-1);
    for (j=0;j<rs;j++){
      a1= *(*(ndl_ring+i)+j);
      b = get_ndl_bond(a1,a2);
      if (!strcmp(ndl_atom[a1]->atype,"A  ") || !strcmp(ndl_atom[a1]->atype,"Q  ")) wcr = true;
      if ((ndl_bond[b]->btype == 'l') || (ndl_bond[b]->btype == 's') ||
        (ndl_bond[b]->btype == 'd') || (ndl_bond[b]->btype == 'a'))  wcr= true;
      a2 = a1;   
    }
    if (wcr == false) {  // if yes, unflag all atoms and bonds in this ring
        a2 = *(*(ndl_ring+i)+rs-1);
        for (j=0;j<rs;j++) {
            a1 = *(*(ndl_ring+i)+j);
            b = get_ndl_bond(a1,a2);
            ndl_atom[a1]->q_arom = false;
            ndl_bond[b]->q_arom = false;
            a2 = a1;   
        }
    }
  }
  // some further refinement would be necessary here in order to unflag everything
  // which contains a wildcard but which definitely cannot be aromatic
}

void set_ndl_atom_tags() {
int  i;
  if (ndl_n_atoms > 0) {
    for (i=1;i<=ndl_n_atoms;i++) ndl_atom[i]->tag = true;
  }
}

void cv_init() {  // new in v0.3j
int i;
  if (cv == NULL) return;
  for (i= 1;i<=ndl_n_atoms;i++) {
      cv[i]->def = ndl_atom[i]->neighbor_count;
  }
}

int cv_count() { // new in v0.3j, modified in v0.3m
int  i, j ;
int  cvlist[max_atoms];
int  cvdef;
bool  isnew;
int  entries;
  if (cv == NULL) {
      return 0;
  }
  memset(cvlist,0,max_atoms*sizeof(int));
  entries = 0;
  for (i = 1;i<= ndl_n_atoms;i++) {
      if (ndl_atom[i]->heavy == true) {
          cvdef = cv[i]->def;
          isnew = true;
          if (entries > 0) {
              for (j = 1;j<=entries;j++) {
                if (cvlist[j] == cvdef) isnew = false;
              }
          }
          if (isnew) {
              entries++;
              cvlist[entries] = cvdef;
          }
          // now we have a list of unique connection values
      }
  }
  return entries;
}

int cv_iterate(int n_cv_prev) {  // new in v0.3j, modified in v0.3m
int  i, j, nb_count;
NEIGHBOR_REC nb;
int  nnb, nsum, n_cv;
  if ((cv == NULL) || (ndl_n_atoms == 0)) return 0;
  // update the connection values (Morgan algorithm)
  for (i= 1;i<=ndl_n_atoms;i++){
    if (ndl_atom[i]->heavy == true) {
      nb_count = 0;
      for (j= 0;j< ndl_n_bonds;j++ ){ //get ra2's neighbors!
        if ((ndl_bond[j]->a1 == i) && (nb_count < max_neighbors) && (ndl_atom[ndl_bond[j]->a2]->heavy)) {
          nb[nb_count]=ndl_bond[j]->a2;
          nb_count++;
        }
        if ((ndl_bond[j]->a2== i)&& (nb_count < max_neighbors)&& (ndl_atom[ndl_bond[j]->a1]->heavy)) {
          nb[nb_count]=ndl_bond[j]->a1;
          nb_count++;
        }
      }
      nnb  = ndl_atom[i]->neighbor_count;
      nsum = 0;
      if (nnb > 0) {
          for (j=0;j<nnb;j++) {
            if (ndl_atom[(nb[j])]->heavy) nsum = nsum + cv[(nb[j])]->def;
          }  
      }
       cv[i]->tmp = nsum;
    }
  }
  n_cv = cv_count();
  if (n_cv > n_cv_prev) {
      for (i = 1;i<=ndl_n_atoms;i++) {
          cv[i]->def = cv[i]->tmp;
      }
  }    
  return n_cv;
}

bool ndl_alkene_C(int ba) {   // new in v0.3f
bool res;
int  i, ba2;
  res = false;
  if ((ndl_n_atoms > 0) && (ndl_n_bonds > 0)) {
      for (i=0;i<ndl_n_bonds;i++) {
          if ((ndl_bond[i]->a1 == ba) || (ndl_bond[i]->a2 == ba)) {
            if (ndl_bond[i]->a1 == ba) {ba2 = ndl_bond[i]->a2;} else {
                ba2 = ndl_bond[i]->a1;}
            if (!strcmp(ndl_atom[ba]->atype,"C2 ") && !strcmp(ndl_atom[ba2]->atype,"C2 ") &&
               (ndl_bond[i]->btype =='D') && (ndl_bond[i]->arom == false)) res = true;
          }
      }
  }
  return res;  
}

int ndl_hetatom_count(int a) {
int  i,j, nb_count;
NEIGHBOR_REC nb;
int  hac;
hac = 0;
if ((a > 0) && (a <= ndl_n_atoms)) {
  nb_count = 0;
  for (j= 0;j< ndl_n_bonds;j++ ){ //get ra2's neighbors!
    if ((ndl_bond[j]->a1 == a) && (nb_count < max_neighbors) && (ndl_atom[ndl_bond[j]->a2]->heavy)) {
      nb[nb_count]=ndl_bond[j]->a2;
      nb_count++;
    }
    if ((ndl_bond[j]->a2== a) && (nb_count < max_neighbors) && (ndl_atom[ndl_bond[j]->a1]->heavy)) {
      nb[nb_count]=ndl_bond[j]->a1;
      nb_count++;
    }
  }
  if (ndl_atom[a]->neighbor_count > 0) {
    for (i= 0;i<ndl_atom[a]->neighbor_count;i++) {
 // note: query atoms like 'A' should be present only in the needle
    if ((strcmp(ndl_atom[(nb[i])]->element,"C ") && strcmp(ndl_atom[(nb[i])]->element,"A ") &&
      strcmp(ndl_atom[(nb[i])]->element,"H ") && strcmp(ndl_atom[(nb[i])]->element,"D ") &&
      strcmp(ndl_atom[(nb[i])]->element,"LP") && strcmp(ndl_atom[(nb[i])]->element,"DU"))) {
      hac++;}
      }
    }
  }
  return hac;
}

int find_ndl_ref_atom() {
int  i, index;
int  score, score_max ; // v0.4b
int  n_nb, n_hc;
  // finds a characteristic atom in the needle molecule,
  // i.e., one with as many substituents as possible and
  // with as many heteroatom substitutents as possible;
  // added in v0.2d: make sure that reference atom is a heavy atom
  // and not (accidentally) an explicit hydrogen;
  // new in v0.3d: special treatment in case of E/Z geometry search
  // to ensure that the entire A-B=C-D fragment is enclosed in one
  // matchpath, regardless where the recursive search starts;
  // refined in v0.3f: exclude only alkene-C as reference atoms
  // added in v0.3o: needle atom must be "tagged" in order to be
  // selected (prevents unconnected fragments from being overlooked)
  if (ndl_n_atoms == 0) return 0;
  score = -1;
  score_max = -1;  // v0.4b
  index = 0;
  if ((ez_search) && (ndl_n_heavyatoms > 2)) {
      for (i=1;i<=ndl_n_atoms;i++) { // ignore sp2-carbons if not aromatic
          //if ((ndl_atom^[i].atype <> 'C2 ') or (ndl_atom^[i].arom = true)) then
          if ((ndl_alkene_C(i) == false) && (ndl_atom[i]->tag)) {
              n_nb = ndl_atom[i]->neighbor_count;
              n_hc = ndl_hetatom_count(i);
              if ((11*n_nb + 7*n_hc > score) && (ndl_atom[i]->heavy)) {
                  index = i;
                  score = 11*n_nb + 7*n_hc;  // changed in v0.3j
              }
          }
          if (score > score_max) score_max = score;  // v0.4b
      }
  }
  // it is possible that no suitable reference atom has been found here
  // (e.g., with "pure" polyenes), so we need a fallback option anyway
  if (index == 0) {
      ez_search = false;  // just in case it was true
      opt_geom = false;   // just in case it was true
      for (i = 1;i<=ndl_n_atoms;i++) {
          n_nb = ndl_atom[i]->neighbor_count;
          n_hc = ndl_hetatom_count(i);
          if ((11*n_nb + 7*n_hc > score) && (ndl_atom[i]->heavy) &&  // v0.3j
             (ndl_atom[i]->tag)) {
              index = i;
              score = 11*n_nb + 7*n_hc;  // changed in v0.3j
          }
          if (score > score_max) score_max = score;  // v0.4b
      }
  }
  // now index must be > 0 in any case (except for H2, or all tags have been cleared)
  if (score_max < 33) {// v0.4b   fallback for unbranched fragments (use terminal atom!)
      for (i = 1;i<=ndl_n_atoms;i++) {
          n_nb = ndl_atom[i]->neighbor_count;
          if ((n_nb == 1) && (ndl_atom[i]->heavy) && (ndl_atom[i]->tag)) index= i;   
      }                     // v0.4c (added heavy atom check); v0.4d (added tag check)
  }
  if (index == 0) index++;  // just to be sure...
  return index;  
}

int find_ndl_ref_atom_cv() {   // new in v0.3j, modified in v0.3m
int  i, n, res, it;
int  n_cv, n_cv_prev;
bool  finished;
int  cvmax;
  if (ndl_n_atoms == 0) {
      return 0;
  }
  res = 1;
  for(n=0;n<max_atoms;n++) {
    cv[n]=(CONVAL_REC *)malloc(sizeof(CONVAL_REC));}
  if ( cv==NULL ) {printf("Error allocating memory for cv!");return find_ndl_ref_atom();}
  cv_init();
  n_cv_prev = 0;
  it = 0;
  finished = false;
  while ((finished==false) && (it < 10000)) {
    it++;  // iteration counter (a safeguard against infinite loops)
    n_cv = cv_iterate(n_cv_prev);
    if (n_cv <= n_cv_prev) finished = true;
    n_cv_prev = n_cv;
  }
  // now that we have canonical connection values (Morgan algorithm),
  // pick the atom with the highest value
  // added in v0.3o: atom must be "tagged"
  cvmax = 0;
  for (i =1;i<=ndl_n_atoms;i++) {
//      printf("cv for atom %d: %d\n",i,cv[i]->def);
      if ((cv[i]->def > cvmax) && ((ndl_alkene_C(i) == false) || (ez_search == false))
          && (ndl_atom[i]->tag)) {
          cvmax = cv[i]->def;
          res  = i;
      }
  }
    if (cv != NULL) {
        free(*cv);
        *cv = NULL;
    }  
  return res;
}

bool gmm_collision() {   // v0.4b
// scans global match matrix for collisions: i.e., matches of a single haystack
//atom with (at least) two different needle atoms which do not have an alternative
// match (e.g., 4-methylheptane would match 1,4-dimethylcyclohexane
int  i, j, k, l, n_hits, hit_pos, n_alt, n_coll, dup_pos;
bool  res;
  res = false;
  if (!(use_gmm && valid_gmm)) return res;
  n_coll = 0;
  for (i= 1;i<=ndl_n_atoms;i++) {
      n_hits  = 0;
      hit_pos = 0;
      for (j = 1;j<=n_atoms;j++) {
          if (gmm[i][j]) {
              n_hits++;
              hit_pos = j;
          }
      }
      if (n_hits == 1) {
          dup_pos = 0;
          for (k = 1;k<=ndl_n_atoms;k++) {
              if ((k!=i) && gmm[k][hit_pos]) {
                  dup_pos = k;
                  n_alt = 0;  // number of alternative matches
                  for (l= 1;l<=n_atoms;l++) {
                      if ((l != hit_pos) && gmm[dup_pos][l]) n_alt++;
                  }
                  if (n_alt < 1) n_coll++;
              }
          }  
      }
  }
  if (n_coll > 0) res = true;
  return res;
}

int matchpath_pos(int id, MATCHPATH_TYPE a_path) {
int  i, pp;
  pp = 0;
  for (i = 0;i<max_matchpath_length;i++) {
      if (a_path[i] == id) {
          pp = i;
          break;
      }
  } 
  return pp;
}

//============================= geometry functions ==========================

double dist3d(P_3D p1, P_3D p2) {
double res;
  res    = sqrt(pow((p1.x-p2.x),2) + pow((p1.y-p2.y),2) + pow((p1.z-p2.z),2));
  return res;
}

P_3D subtract_3d(P_3D p1, P_3D p2) {
P_3D  p ;
  p.x = p1.x - p2.x;
  p.y = p1.y - p2.y;
  p.z = p1.z - p2.z;
  return p;
}

P_3D add_3d(P_3D p1, P_3D p2) {
P_3D  p;
  p.x = p1.x + p2.x;
  p.y = p1.y + p2.y;
  p.z = p1.z + p2.z;
  return p;
}

/*
void vec2origin(P_3D p1, P_3D p2) {
P_3D  p;
  p = subtract_3d(p2,p1);
  p2 = p;
  p1.x = 0; p1.y = 0; p1.z = 0;
}
*/

double scalar_prod(P_3D p1, P_3D p2, P_3D p3) {
P_3D  p;
double  res ;
  p = subtract_3d(p2,p1);
  p2 = p;
  p = subtract_3d(p3,p1);
  p3 = p;
  p1.x = 0; p1.y = 0; p1.z = 0;
  res = p2.x*p3.x + p2.y*p3.y + p2.z*p3.z;
  return res;
}

P_3D cross_prod(P_3D p1, P_3D p2, P_3D p3) {
P_3D  p, orig_p1;
  orig_p1 = p1;
  p = subtract_3d(p2,p1);
  p2 = p;
  p = subtract_3d(p3,p1);
  p3 = p;
  p.x = p2.y*p3.z - p2.z*p3.y;
  p.y = p2.z*p3.x - p2.x*p3.z;
  p.z = p2.x*p3.y - p2.y*p3.x;
  return add_3d(orig_p1,p);
}

double angle_3d(P_3D p1, P_3D p2, P_3D p3) {
P_3D  lp1,lp2,lp3,p;
double  res,  magn_1, magn_2,  cos_phi;
  res = 0;
  lp1 = p1; lp2 = p2; lp3 = p3;
  p = subtract_3d(lp2,lp1);
  lp2 = p;
  p = subtract_3d(lp3,lp1);
  lp3 = p;
  lp1.x = 0; lp1.y = 0; lp1.z = 0;
  magn_1 = dist3d(lp1,lp2);
  magn_2 = dist3d(lp1,lp3);
  if (magn_1 * magn_2 == 0) {   // emergency exit
      return pi ;
  }
  cos_phi = scalar_prod(lp1,lp2,lp3) / (magn_1 * magn_2);
  if (cos_phi < -1) cos_phi = -1;
  if (cos_phi > 1)  cos_phi = 1;
  res = acos(cos_phi);
  return res;
}

double torsion(P_3D p1, P_3D p2, P_3D p3, P_3D p4) {
P_3D  lp1,lp2,lp3,lp4,d1,c1,c2,c1xc2,c2xc1;
double  res, dist1,dist2,sign;
  // copy everything into local variables
  lp1 = p1; lp2 = p2; lp3 = p3; lp4 = p4;
  // get the vector between the two central atoms
  d1 = subtract_3d(p3,p2);
  // shift the first atom parallel to be attached to p3 instead of p2
  lp1 = add_3d(p1,d1);
  // now get the cross product vectors
  c1 = cross_prod(lp3,lp2,lp1);
  c2 = cross_prod(lp3,lp2,lp4);
  res = angle_3d(p3,c1,c2);
  //now check if it is clockwise or anticlockwise:
  //first, make the cross products of the two cross products c1 and c2 (both ways)
  c1xc2 = cross_prod(lp3,c1,c2);
  c2xc1 = cross_prod(lp3,c2,c1);
  //next, get the distances from these points to our refernce point lp2
  dist1 = dist3d(lp2,c1xc2);
  dist2 = dist3d(lp2,c2xc1);
  if (dist1 <= dist2) {sign = 1;} else {sign = -1;}
  return sign*res;
}


double ctorsion(P_3D p1, P_3D p2, P_3D p3, P_3D p4) {
// calculates "pseudo-torsion" defined by atoms 3 and 4, being both
// attached to atom 2, with respect to axis of atoms 1 and 2
P_3D  lp1,lp2,lp3,lp4,c1,c2;
double  res, dist1,dist2,sign;
P_3D  c1xc2, c2xc1;
  // copy everything into local variables
  lp1 = p1; lp2 = p2; lp3 = p3; lp4 = p4;
  // get the cross product vectors
  c1 = cross_prod(lp2,lp1,lp3);
  c2 = cross_prod(lp2,lp1,lp4);
  res = angle_3d(p2,c1,c2);
  //now check if it is clockwise or anticlockwise:
  //first, make the cross products of the two cross products c1 and c2 (both ways)
  c1xc2 = cross_prod(lp2,c1,c2);
  c2xc1 = cross_prod(lp2,c2,c1);
  //next, get the distances from these points to our refernce point lp1
  dist1 = dist3d(lp1,c1xc2);
  dist2 = dist3d(lp1,c2xc1);
  if (dist1 <= dist2) {sign = 1;} else {sign = -1;}
  return sign*res;
}


bool is_cis(P_3D p1, P_3D p2, P_3D p3, P_3D p4) {  // new in v0.3h, uses the dihedral angle
double  phi;
bool  res;
  res = false;
  phi = torsion(p1,p2,p3,p4);
  if (abs(phi) < pi/2) res = true;
  return res;
}

//====================== end of geometry functions ==========================

bool atomtypes_OK_strict(int ndl_a, int hay_a) { // new in v0.2f
int  ndl_nbc , hay_nbc;
short  ndl_Hexp , hay_Htot;
bool  res;
  res = false;   
  ndl_nbc = ndl_atom[ndl_a]->neighbor_count;
  ndl_Hexp = ndl_atom[ndl_a]->Hexp;
  hay_nbc = atom[hay_a]->neighbor_count;
  hay_Htot = atom[hay_a]->Htot;
  // v0.3o: formal charges must be the same
  if (ndl_atom[ndl_a]->formal_charge != atom[hay_a]->formal_charge) {
      return false;
  }
  // v0.3p: isotope nucleon numbers must be the same
  if (ndl_atom[ndl_a]->nucleon_number != atom[hay_a]->nucleon_number) {
      return false;
  }
  // v0.3p: radicals numbers must be the same
  if (ndl_atom[ndl_a]->radical_type != atom[hay_a]->radical_type) {
      return false;
  }
  if (!strcmp(ndl_atom[ndl_a]->atype,atom[hay_a]->atype)) {res= true; } else {
      if (!strcmp(ndl_atom[ndl_a]->element, atom[hay_a]->element)) {
          if ((ndl_atom[ndl_a]->arom) && (atom[hay_a]->arom)) res = true;
          if (ndl_querymol && (ndl_atom[ndl_a]->q_arom) &&
            (atom[hay_a]->arom)) res = true;   // v0.3p
      }
  }
  if (!strcmp(ndl_atom[ndl_a]->element, "A ") && (atom[hay_a]->heavy)) res = true;
  if (!strcmp(ndl_atom[ndl_a]->element, "Q ")) {
      if ((atom[hay_a]->heavy) && strcmp(atom[hay_a]->element,"C ")) res = true;
  }
  if (!strcmp(ndl_atom[ndl_a]->element,"X ")) {
      if (!strcmp(atom[hay_a]->element,"F ") || !strcmp(atom[hay_a]->element,"CL") ||
        !strcmp(atom[hay_a]->element,"BR") || !strcmp(atom[hay_a]->element,"I ")) res = true;
  }
  // if needle atom has more substituents than haystack atom ==> no match
  if (ndl_nbc > hay_nbc) res = false;
  // check for explicit hydrogens
  if (ndl_Hexp > hay_Htot) res = false;
  // new in v0.3m: in "fingerprint mode", also query atom symbols must match
  if (opt_fp) {
      if (strcmp(ndl_atom[ndl_a]->element, atom[hay_a]->element)) res = false;
  }
  return res;
}

bool atomtypes_OK(int ndl_a, int hay_a) {
int ndl_nbc, hay_nbc;
short  ndl_Hexp, hay_Htot;
bool  res;
  if ((ndl_a < 1) || (ndl_a > ndl_n_atoms) || (hay_a < 1) || (hay_a > n_atoms)) {
      return false;
  }
  // check for opposite charges;  v0.3l, refined in v0.3o, v0.3p
  // except in strict mode, matching pairs of charged+uncharged atoms 
  // are tolerated (this is a feature, not a bug)
  //if ((ndl_atom^[ndl_a].formal_charge <> 0) and (atom^[hst_a].formal_charge <> 0) 
  //and (ndl_atom^[ndl_a].formal_charge <> atom^[hst_a].formal_charge)) then    
  if (opt_chg) {
      if (ndl_atom[ndl_a]->formal_charge != atom[hay_a]->formal_charge) {
          return false;
      }
  }
  // v0.3p: isotopes must be the same
  if (opt_iso) {
      if (ndl_atom[ndl_a]->nucleon_number != atom[hay_a]->nucleon_number) {
          return false;
      }
  }
  // v0.3p: radicals must be the same
  if (opt_rad) {
      if (ndl_atom[ndl_a]->radical_type != atom[hay_a]->radical_type) {
          return false;
      }
  }
  // in exact mode, check if (disconnected) fragment is already tagged; v0.3o
  if ((opt_exact) && (atom[hay_a]->tag == true)) {
      return false;
  }
  if (opt_strict) {
      return atomtypes_OK_strict(ndl_a, hay_a);
  }  
  res = false;   
  ndl_nbc = ndl_atom[ndl_a]->neighbor_count;
  ndl_Hexp = ndl_atom[ndl_a]->Hexp;
  hay_nbc = atom[hay_a]->neighbor_count;
  hay_Htot = atom[hay_a]->Htot;
  if (!strcmp(ndl_atom[ndl_a]->element, atom[hay_a]->element)) res = true;  // very simplified...
  if (!strcmp(ndl_atom[ndl_a]->element,"A ") && (atom[hay_a]->heavy)) res = true;
  if (!strcmp(ndl_atom[ndl_a]->element,"Q ")) {
      if ((atom[hay_a]->heavy) && strcmp(atom[hay_a]->element,"C ")) res = true;
  }
  if (!strcmp(ndl_atom[ndl_a]->element,"X ")) {
      if (!strcmp(atom[hay_a]->element,"F ") || !strcmp(atom[hay_a]->element,"CL") || 
        !strcmp(atom[hay_a]->element,"BR") || !strcmp(atom[hay_a]->element,"I ")) res = true;
  }
  // v0.3o: in exact mode, check for identical neighbor_count
  if (opt_exact) {
      if (ndl_nbc != hay_nbc) {
          res = false;//'exact match failed: different number of neighbor atoms'
      }
  }
  // if needle atom has more substituents than haystack atom ==> no match
  if (ndl_nbc > hay_nbc) res = false;
  // check for explicit hydrogens
  if (ndl_Hexp > hay_Htot) res = false;
  return res;
}

bool bondtypes_OK_strict(int ndl_b, int hay_b) {
bool  ndl_arom, hay_arom;
char  ndl_btype, hay_btype;
int  ndl_rc, hay_rc;   // new in v0.3d
short  ndl_btopo ;  // new in v0.3d
bool  res ;
  res       = false;
  ndl_arom  = ndl_bond[ndl_b]->arom;
  ndl_btype = ndl_bond[ndl_b]->btype;
  ndl_rc    = ndl_bond[ndl_b]->ring_count;
  ndl_btopo = ndl_bond[ndl_b]->topo;
  hay_arom  = bond[hay_b]->arom;
  hay_btype = bond[hay_b]->btype;
  hay_rc    = bond[hay_b]->ring_count;
  if ((ndl_arom == true) && (hay_arom == true)) {
      res = true; 
      if ((ndl_btype == 'T') && (hay_btype != 'T')) res = false;  // v0.4b
      if ((hay_btype == 'T') && (ndl_btype != 'T') && (ndl_btype != 'a')) res= false; 
  }
  if ((ndl_arom == false) && (hay_arom == false)) {
      if (ndl_btype == hay_btype) res = true;
      if ((ndl_btype == 'l') && ((hay_btype == 'S') || (hay_btype == 'D'))) res = true;
      if ((ndl_btype == 's') && (hay_btype = 'S')) res = true;
      if ((ndl_btype == 'd') && (hay_btype = 'D')) res = true;
  }
  // a little exception:
  if ((ndl_arom == false) && (hay_arom == true)) {
      if (ndl_btype == 'A') res = true;
      if ((ndl_btype == 's') || (ndl_btype == 'd')) res = true;
      if (ndl_bond[ndl_b]->q_arom == true) res = true;  // v0.3p
  }
  if (ndl_btype == 'a') res = true;
    // new in v0.3d: strict comparison of topology (and even ring_count!)
  if ((ndl_btopo < btopo_always_any) || (ndl_btopo == btopo_exact_rc)) {
      if (ndl_rc != hay_rc) {
          res= false;     // this excludes further ring annulations as well as
      }
  } else {
      if ((ndl_btopo == btopo_excess_rc) && (hay_rc <= ndl_rc)) res = false;
  }
  return res;
}

bool bondtypes_OK(int ndl_b, int hay_b) {
bool  ndl_arom, hay_arom;
char  ndl_btype, hay_btype;
int  ndl_rc, hay_rc;   // new in v0.3d
short  ndl_btopo;  // new in v0.3d
bool  res;
int  a1, a2;
  if ((ndl_b < 0) || (ndl_b > ndl_n_bonds) ||
     (hay_b < 0) || (hay_b > n_bonds)) return false;
  if (opt_strict) {     
      return bondtypes_OK_strict(ndl_b, hay_b);
  }
  res = false;
  ndl_arom = ndl_bond[ndl_b]->arom;
  ndl_btype = ndl_bond[ndl_b]->btype;
  hay_arom = bond[hay_b]->arom;
  hay_btype = bond[hay_b]->btype;
  ndl_rc    = ndl_bond[ndl_b]->ring_count;
  hay_rc    = bond[hay_b]->ring_count;
  ndl_btopo = ndl_bond[ndl_b]->topo;
  if ((ndl_arom == true) && (hay_arom == true)) {
      res = true;
      if ((ndl_btype == 'T') && (hay_btype != 'T')) res = false;  // v0.4b
      if ((hay_btype == 'T') && (ndl_btype != 'T') && (ndl_btype != 'a')) res = false; 
  }
  if ((ndl_arom == false) && (hay_arom == false)) {
      if (ndl_btype == hay_btype) res = true;
      if ((ndl_btype == 'l') && ((hay_btype == 'S') || (hay_btype == 'D'))) res = true;
      if ((ndl_btype == 's') && (hay_btype == 'S')) res = true;
      if ((ndl_btype == 'd') && (hay_btype == 'D')) res = true;
  }
  // a little exception:
  if ((ndl_arom == false) && (hay_arom == true)) {
    if (ndl_btype == 'A') res = true;
    if ((ndl_btype == 's') || (ndl_btype == 'd')) res = true;
    if (ndl_btype == 'D') {   // added in 0.2d: do not accept C=O etc. as C-O/arom
      a1 = ndl_bond[ndl_b]->a1;
      a2 = ndl_bond[ndl_b]->a2;
      if (!(!strcmp(ndl_atom[a1]->element,"O ") || !strcmp(ndl_atom[a2]->element,"O ") ||
      !strcmp(ndl_atom[a1]->element,"S ") || !strcmp(ndl_atom[a2]->element,"S ") ||
      !strcmp(ndl_atom[a1]->element,"SE") || !strcmp(ndl_atom[a2]->element,"SE") ||
      !strcmp(ndl_atom[a1]->element,"TE") || !strcmp(ndl_atom[a1]->element,"TE"))) res= true;
    }
    if (ndl_bond[ndl_b]->q_arom == true) res= true;  // v0.3p
  }
  if (ndl_btype == 'a') res = true;//printf("%d,%d,%d,%d,%c,%d\n",ndl_b, ndl_arom, hay_arom,ndl_btopo,ndl_btype,res);
  // new in v0.3d: obey topology requirements in query structure
  if ((ndl_btopo != btopo_any) && (ndl_btopo != btopo_always_any)) {
      if ((ndl_btopo == btopo_ring)  && (hay_rc == 0)) res = false;
      if ((ndl_btopo == btopo_chain) && (hay_rc > 0)) res = false;
      if ((ndl_btopo == btopo_excess_rc) && (hay_rc <= ndl_rc)) res = false;
      if ((ndl_btopo == btopo_exact_rc)  && (hay_rc != ndl_rc)) res = false;
  }
  return res;
}

bool matrix_OK(MATCHMATRIX m, int ndl_dim, int hay_dim) {
  // new, recursive version in v0.2i: can handle up to max_neighbors substituents
bool  mr ;
MATCHMATRIX  lm;
int  i, ii, j, lhay_dim;
  if ((ndl_dim < 1) || (ndl_dim > max_neighbors) ||
     (hay_dim < 1) || (hay_dim > max_neighbors) ||
     (ndl_dim > hay_dim)) return false;
  mr = false;
  if (ndl_dim == 1) {
      for (i=0;i<hay_dim;i++) {if (m[0][i]) mr = true;}
  } else {
      for (i=0;i<hay_dim;i++) {
          if (m[0][i]) {  // write remaining fields into a new matchmatrix which is smaller by 1x1
              memset(lm,0,sizeof(MATCHMATRIX));
              for (j = 1;j<ndl_dim;j++) {
                  lhay_dim = 0;
                  for (ii=0;ii<hay_dim;ii++) {
                      if (ii != i) {
                          lm[(j-1)][lhay_dim] = m[j][ii];
                          lhay_dim++;
                      }
                  }
              }
              if (matrix_OK(lm,ndl_dim-1,lhay_dim)) {   // recursive call to matrix_OK
                  return true;  // stop any further work immediately
              } 
          }
      }
  }
  return mr;
}

bool is_flat(double angle_deg) {  // new in v0.3j
 if ((abs(angle_deg) > 5) && (abs(angle_deg) < 175)) {
   return false;} else {return true;}
}

bool chirality_OK(CHIRPATH_TYPE ndl_cp, CHIRPATH_TYPE hay_cp) {
bool  res;
double  ndl_ct, hay_ct, ndl_ct_deg, hay_ct_deg;
P_3D  np1,np2,np3,np4,hp1,hp2,hp3,hp4;
int  level, i;
bool  up, down, updown;
int  ta1,ta2,ta3,ta4,ba1,ba2;
  if (opt_nochirality) {
    return true;
  }
  res = true;
  // fill temporary atom variables
  ta1 = ndl_cp[0];  // this is the central atom
  ta2 = ndl_cp[1];
  ta3 = ndl_cp[2];
  ta4 = ndl_cp[3];
  // first, get the central atom of the needle
  np2.x = ndl_atom[ta1]->x;
  np2.y = ndl_atom[ta1]->y;
  np2.z = ndl_atom[ta1]->z;
  // next, do the same for all 3 substituent atoms
  np1.x = ndl_atom[ta2]->x;
  np1.y = ndl_atom[ta2]->y;
  np1.z = ndl_atom[ta2]->z;
  np3.x = ndl_atom[ta3]->x;
  np3.y = ndl_atom[ta3]->y;
  np3.z = ndl_atom[ta3]->z;
  np4.x = ndl_atom[ta4]->x;
  np4.y = ndl_atom[ta4]->y;
  np4.z = ndl_atom[ta4]->z;
  // now check all needle bonds if we should care about up/down bonds
  level  = 0;
  updown = false;
  up     = false;
  down   = false;
  if (ndl_n_bonds > 0) {
    for (i =0;i<ndl_n_bonds;i++) {
      if ((ndl_bond[i]->stereo == bstereo_up) || (ndl_bond[i]->stereo == bstereo_down)) {
          ba1 = ndl_bond[i]->a1; ba2 = ndl_bond[i]->a2;
          if ((ba1 == ta1) && (ndl_bond[i]->stereo == bstereo_up)) {
              up = true;
              if ((ba2 == ta2) || (ba2 == ta3) || (ba2 == ta4)) {
                   updown = true;
                   if (ba2 == ta2) np1.z = np1.z + 0.8;
                   if (ba2 == ta3) np3.z = np3.z + 0.8;
                   if (ba2 == ta4) np4.z = np4.z + 0.8;
              } else {level = level + 1;}
          }
          if ((ba1 == ta1) && (ndl_bond[i]->stereo == bstereo_down)) {
              down = true;
              if ((ba2 == ta2) || (ba2 == ta3) || (ba2 == ta4)) {
                  updown = true;
                  if (ba2 == ta2) np1.z = np1.z - 0.8;
                  if (ba2 == ta3) np3.z = np3.z - 0.8;
                  if (ba2 == ta4) np4.z = np4.z - 0.8;
              } else {level = level - 1;}
          }
          if ((ba2 == ta1) && (ndl_bond[i]->stereo == bstereo_up)) {
              down = true;
              if ((ba1 == ta2) || (ba1 == ta3) || (ba1 == ta4)) {
                  updown = true;
                  if (ba1 == ta2) np1.z = np1.z - 0.8;
                  if (ba1 == ta3) np3.z = np3.z - 0.8;
                  if (ba1 == ta4) np4.z = np4.z - 0.8;
              } else {level = level - 1; }
         }
         if ((ba2 == ta1) && (ndl_bond[i]->stereo == bstereo_down)) {
             up = true;
             if ((ba1 == ta2) || (ba1 == ta3) || (ba1 == ta4)) {
                  updown = true;
                  if (ba1 == ta2) np1.z = np1.z + 0.8;
                  if (ba1 == ta3) np3.z = np3.z + 0.8;
                  if (ba1 == ta4) np4.z = np4.z + 0.8;
              } else {level = level + 1;}
          }
       }  
    }  // for i ...
      if ((updown == false) && (level != 0)) {
          if (level > 0) np2.z = np2.z + 0.3;
          if (level < 0) np2.z = np2.z - 0.3;
      } else {
          if (up)   np2.z = np2.z + 0.1;
          if (down) np2.z = np2.z - 0.1;
      }
  }  
  // fill temporary atom variables again
  ta1 = hay_cp[0];
  ta2 = hay_cp[1];
  ta3 = hay_cp[2];
  ta4 = hay_cp[3];
  // then, get the central atom of the haystack
  hp2.x = atom[ta1]->x;
  hp2.y = atom[ta1]->y;
  hp2.z = atom[ta1]->z;
  // next, do the same for all 3 substituent atoms
  hp1.x = atom[ta2]->x;
  hp1.y = atom[ta2]->y;
  hp1.z = atom[ta2]->z;
  hp3.x = atom[ta3]->x;
  hp3.y = atom[ta3]->y;
  hp3.z = atom[ta3]->z;
  hp4.x = atom[ta4]->x;
  hp4.y = atom[ta4]->y;
  hp4.z = atom[ta4]->z;
  // now check all haystack bonds if we should care about up/down bonds
  level  = 0;
  updown = false;
  up     = false;
  down   = false;
  if (n_bonds > 0) {
      for (i=0;i<n_bonds;i++) {
          if ((bond[i]->stereo == bstereo_up) || (bond[i]->stereo == bstereo_down)) {
              ba1 = bond[i]->a1; ba2 = bond[i]->a2;
              if ((ba1 == ta1) && (bond[i]->stereo == bstereo_up)) {
                  up = true;
                  if ((ba2 == ta2) || (ba2 == ta3) || (ba2 == ta4)) {
                      updown = true;
                      if (ba2 == ta2) hp1.z = hp1.z + 0.8;
                      if (ba2 == ta3) hp3.z = hp3.z + 0.8;
                      if (ba2 == ta4) hp4.z = hp4.z + 0.8;
                   } else {level = level + 1;}
              }
              if ((ba1 == ta1) && (bond[i]->stereo == bstereo_down)) {
                  down = true;
                  if ((ba2 == ta2) || (ba2 == ta3) || (ba2 = ta4)) {
                      updown = true;
                      if (ba2 == ta2) hp1.z = hp1.z - 0.8;
                      if (ba2 == ta3) hp3.z = hp3.z - 0.8;
                      if (ba2 == ta4) hp4.z = hp4.z - 0.8;
                  } else { level = level - 1;}
              }
              if ((ba2 == ta1) && (bond[i]->stereo == bstereo_up)) {
                  down = true;
                  if ((ba1 == ta2) || (ba1 == ta3) || (ba1 == ta4)) {
                      updown = true;
                      if (ba1 == ta2) hp1.z = hp1.z - 0.8;
                      if (ba1 == ta3) hp3.z = hp3.z - 0.8;
                      if (ba1 == ta4) hp4.z = hp4.z - 0.8;
                  } else { level = level - 1;}
              }
              if ((ba2 == ta1) && (bond[i]->stereo == bstereo_down)) {
                  up = true;
                  if ((ba1 == ta2) || (ba1 == ta3) || (ba1 = ta4)) {
                      updown = true;
                      if (ba1 == ta2) hp1.z = hp1.z + 0.8;
                      if (ba1 == ta3) hp3.z = hp3.z + 0.8;
                      if (ba1 == ta4) hp4.z = hp4.z + 0.8;
                  } else {level = level + 1;}
              }
          }  
      }  // for i ...
      if ((updown == false) && (level != 0)) {
          if (level > 0) hp2.z = hp2.z + 0.3;
          if (level < 0) hp2.z = hp2.z - 0.3;
      } else {
          if (up )  hp2.z = hp2.z + 0.1;
          if (down) hp2.z = hp2.z - 0.1;
      }
  }  
  // get the pseudo-torsion angles
  ndl_ct = ctorsion(np1,np2,np3,np4);
  hay_ct = ctorsion(hp1,hp2,hp3,hp4);
  ndl_ct_deg = ndl_ct*180/pi;
  hay_ct_deg = hay_ct*180/pi;
  // now do a plausibility check and finally check the sense
  // (clockwise or counterclockwise)
  if ((! is_flat(ndl_ct_deg)) &&
     (! is_flat(hay_ct_deg)) &&
     (ndl_ct_deg * hay_ct_deg < 0)) res= false;
  if (rs_strict) {
      if ((is_flat(ndl_ct_deg) && (!is_flat(hay_ct_deg))) ||
         (is_flat(hay_ct_deg) && (!is_flat(ndl_ct_deg))) ||
         (ndl_ct_deg * hay_ct_deg < 0)) res = false;
  }
  return res;
}

bool ndl_maybe_chiral(int na) { // new in v0.3h
bool  res;
int  n_nb;
  res  = false;
  n_nb = ndl_atom[na]->neighbor_count;
  if (!strcmp(ndl_atom[na]->atype,"C3 ") && (n_nb > 2)) res = true;
  if (!strcmp(ndl_atom[na]->element,"N ")) {
      if (!strcmp(ndl_atom[na]->atype,"N3+") && (n_nb == 4)) res = true;
  }
  if (!strcmp(ndl_atom[na]->element,"S ")) {
      if ((n_nb == 3) && (ndl_hetatom_count(na) == 1)) res = true;
  }
  if (!strcmp(ndl_atom[na]->element,"P ") || !strcmp(ndl_atom[na]->element,"AS")) {  
 // "As" added in v0.3j
      if (n_nb > 3) res = true;  // are we missing something here?
      if (ndl_hetatom_count(na) >= 2) res = false;  // v0.3m; ignore phosphates etc.
  }
  return res;
}

bool is_matching(MATCHPATH_TYPE ndl_xmp, MATCHPATH_TYPE hay_xmp) {
int  i, j, k, l, m;
int  ndl_n_nb, n_nb, ndl_a, hay_a, ndl_b, hay_b, nb_count, nb_next_count;
int  prev_ndl_a, prev_hay_a, next_ndl_a, next_hay_a;
NEIGHBOR_REC  ndl_nb, hay_nb;
bool res;
MATCHMATRIX  mm;  
int  ndl_mp_len, hay_mp_len;
MATCHPATH_TYPE  ndl_mp;
MATCHPATH_TYPE  hay_mp;
bool  emptyline ;
bool  ndl_cis, hay_cis;
int  na1,na2,na3,na4;  // v0.3d
int  ha1,ha2,ha3,ha4;  // atom variables for E/Z check
int  prev_ndl_b, prev_hay_b;       //
P_3D  p1,p2,p3,p4;
  //hst_torsion, ndl_torsion : double;
CHIRPATH_TYPE  ncp,hcp;
int  n_hits;
int  n_singlehits;
  // initialize local matchpath variables
  memset(ndl_mp,0,sizeof(MATCHPATH_TYPE));
  memset(hay_mp,0,sizeof(MATCHPATH_TYPE));
  // copy content of external variables into local ones
  for (i=0;i<max_matchpath_length;i++) {
      ndl_mp[i] = ndl_xmp[i];
      hay_mp[i] = hay_xmp[i];
  }
  ndl_mp_len = matchpath_pos(0,ndl_mp);
  hay_mp_len = matchpath_pos(0,hay_mp);
  if (ndl_mp_len != hay_mp_len) {
      // this should never happen....
      res = false;
      return res;
  } 
  ndl_a = ndl_mp[ndl_mp_len-1];
  hay_a = hay_mp[hay_mp_len-1];
  ndl_atom[ndl_a]->tag = false; // new in v0.3o: mark the last needle atom as "visited"
  ndl_b = 0;
  hay_b = 0;
  prev_ndl_a = 0;
  prev_hay_a = 0;
  if (ndl_mp_len > 1) {
      prev_ndl_a = ndl_mp[ndl_mp_len-2];
      prev_hay_a = hay_mp[hay_mp_len-2];
  }
  // if geometry checking is on, check it here
  if ((ez_search == true) && (ndl_mp_len > 3)) {
      na1 = ndl_mp[ndl_mp_len-1];
      na2 = ndl_mp[ndl_mp_len-2];
      na3 = ndl_mp[ndl_mp_len-3];
      na4 = ndl_mp[ndl_mp_len-4];
      ha1 = hay_mp[hay_mp_len-1];
      ha2 = hay_mp[hay_mp_len-2];
      ha3 = hay_mp[hay_mp_len-3];
      ha4 = hay_mp[hay_mp_len-4];
      prev_ndl_b = get_ndl_bond(na2,na3);
      prev_hay_b = get_bond(ha2,ha3);
      if ((ndl_bond[prev_ndl_b]->btype == 'D') && (bond[prev_hay_b]->arom == false) &&
         (!strcmp(atom[ha2]->element,"C ") || !strcmp(atom[ha2]->element,"N ")) &&
         (!strcmp(atom[ha3]->element,"C ") || !strcmp(atom[ha3]->element,"N "))) {
  // v0.3g; check C=C, C=N, N=N bonds
          p1.x = atom[ha1]->x; p1.y = atom[ha1]->y; p1.z = atom[ha1]->z;
          p2.x = atom[ha2]->x; p2.y = atom[ha2]->y; p2.z = atom[ha2]->z;
          p3.x = atom[ha3]->x; p3.y = atom[ha3]->y; p3.z = atom[ha3]->z;
          p4.x = atom[ha4]->x; p4.y = atom[ha4]->y; p4.z = atom[ha4]->z;
          hay_cis = is_cis(p1,p2,p3,p4);
          //hst_torsion := torsion(p1,p2,p3,p4);
          p1.x = ndl_atom[na1]->x; p1.y = ndl_atom[na1]->y; p1.z = ndl_atom[na1]->z;
          p2.x = ndl_atom[na2]->x; p2.y = ndl_atom[na2]->y; p2.z = ndl_atom[na2]->z;
          p3.x = ndl_atom[na3]->x; p3.y = ndl_atom[na3]->y; p3.z = ndl_atom[na3]->z;
          p4.x = ndl_atom[na4]->x; p4.y = ndl_atom[na4]->y; p4.z = ndl_atom[na4]->z;
          //ndl_torsion := torsion(p1,p2,p3,p4);
          ndl_cis = is_cis(p1,p2,p3,p4);
          if (ndl_cis != hay_cis) return false;
      }
  }   // end of E/Z geometry check
  // check whatever can be checked as early as now:
  // e.g. different elements or more substituents on needle atom than on haystack
  if (!atomtypes_OK(ndl_a, hay_a)) return false;
  // positive scenarios, e.g. one-atom fragments  (v0.3o)
  if ((atom[hay_a]->neighbor_count == 0) && (ndl_atom[ndl_a]->neighbor_count == 0)) {
      if (atomtypes_OK(ndl_a, hay_a)) {
          res = true;
          atom[hay_a]->tag = true;
          if (use_gmm && valid_gmm) *(*(gmm+ndl_a)+hay_a) = true;  // v0.4a
      } else { res= false;}
      return res;
  }
  // and other possibilities:
  ndl_b = get_ndl_bond(prev_ndl_a,ndl_a);
  hay_b = get_bond(prev_hay_a,hay_a);  
  if ((ndl_b >= 0) && (hay_b >= 0)) {
      // do a quick check if bond types match
      if (!bondtypes_OK(ndl_b, hay_b)) return false;
  }//printf("matchresult: %d,%d,%d,%d\n",ndl_a,hay_a,ndl_mp_len,res);
  // v0.4b: reversed checks a) and a.1), see below
  // a.1) haystack fragment forms a ring, but needle does not;  v0.3m
  if ((matchpath_pos(ndl_a,ndl_mp) == matchpath_pos(0,ndl_mp)-1) &&
     (matchpath_pos(hay_a,hay_mp) < matchpath_pos(0,hay_mp)-1)) {
       res = false;
       return res;
  }
  // a) we reached the end of our needle fragment (and atom/bond types match)
  if ((ndl_atom[ndl_a]->neighbor_count == 1) && (atomtypes_OK(ndl_a, hay_a)) &&
     (bondtypes_OK(ndl_b, hay_b))) {  
       res= true;
       if (use_gmm && valid_gmm) *(*(gmm+ndl_a)+hay_a) = true;  // v0.4a
       return res;
  }   
//
  // b) a ring is formed (ndl_a is already in the path) and atom/bond types match
  if ((matchpath_pos(ndl_a,ndl_mp) >= 0) && 
    (matchpath_pos(ndl_a,ndl_mp) < matchpath_pos(0,ndl_mp)-1)) {
      if ((matchpath_pos(ndl_a,ndl_mp) == matchpath_pos(hay_a, hay_mp)) &&
         (atomtypes_OK(ndl_a, hay_a)) && (bondtypes_OK(ndl_b, hay_b))) {
           // 1st chirality check
        if ((matchpath_pos(ndl_a,ndl_mp) > 1) && (rs_search || ndl_atom[ndl_a]->stereo_care)
           && ndl_maybe_chiral(ndl_a)) {
               na1 = ndl_a;  // the (potential) chiral center (v0.3f)
               na2 = ndl_mp[(matchpath_pos(ndl_a,ndl_mp)-1)];
               na3 = ndl_mp[(matchpath_pos(ndl_a,ndl_mp)+1)];
               na4 = ndl_mp[(matchpath_pos(0,ndl_mp)-2)];
               ha1 = hay_a;  
               ha2 = hay_mp[(matchpath_pos(hay_a,hay_mp)-1)];
               ha3 = hay_mp[(matchpath_pos(hay_a,hay_mp)+1)];
               ha4 = hay_mp[(matchpath_pos(0,hay_mp)-2)];
               memset(ncp,0,sizeof(CHIRPATH_TYPE));
               memset(hcp,0,sizeof(CHIRPATH_TYPE));
               ncp[1] = na1; ncp[2] = na2; ncp[3] = na3; ncp[4] = na4;
               hcp[1] = ha1; hcp[2] = ha2; hcp[3] = ha3; hcp[4] = ha4;
               if (!chirality_OK(ncp,hcp)) { 
                   res=false;
                   return res;
               } 
           }  // end of 1st chirality check, chirality check succeeded at ring junction
           res = true;
           if (use_gmm && valid_gmm) *(*(gmm+ndl_a)+hay_a) = true;  // v0.4a
           return res;
       } else {
           res = false;
           return res;
       }
  }
  // in all other cases, do the hard work:
  // first, get all heavy-atom neighbors of needle and haystack;
  // at the beginning of the search, this means all neighbors, then it means
  // all but the previous atom (where we came from)
  memset(ndl_nb, 0, sizeof(NEIGHBOR_REC));
  memset(hay_nb, 0, sizeof(NEIGHBOR_REC));
  if (matchpath_pos(0,ndl_mp) == 1) {
      ndl_n_nb = ndl_atom[ndl_a]->neighbor_count;
      n_nb = atom[hay_a]->neighbor_count;
      nb_count = 0;
      for (j= 0;j< ndl_n_bonds;j++ ){ //get ndl_a's neighbors!
        if ((ndl_bond[j]->a1 == ndl_a) && (nb_count < max_neighbors) && (ndl_atom[ndl_bond[j]->a2]->heavy)) {
          ndl_nb[nb_count]=ndl_bond[j]->a2;
          nb_count++;
        }
        if ((ndl_bond[j]->a2== ndl_a) && (nb_count < max_neighbors) && (ndl_atom[bond[j]->a1]->heavy)) {
          ndl_nb[nb_count]=ndl_bond[j]->a1;
          nb_count++;
        }
      }
//      ndl_nb = get_ndl_neighbors(ndl_a);
      nb_count = 0;
      for (j= 0;j< n_bonds;j++ ){ //get hay_a's neighbors!
        if ((bond[j]->a1 == hay_a) && (nb_count < max_neighbors) && (atom[bond[j]->a2]->heavy)) {
          hay_nb[nb_count]=bond[j]->a2;
          nb_count++;
        }
        if ((bond[j]->a2== hay_a)&& (nb_count < max_neighbors)&& (atom[bond[j]->a1]->heavy)) {
          hay_nb[nb_count]=bond[j]->a1;
          nb_count++;
        }
      }
//      hay_nb = get_neighbors(hay_a);
  } else {
      ndl_n_nb = ndl_atom[ndl_a]->neighbor_count -1;
      n_nb = atom[hay_a]->neighbor_count -1;
      nb_next_count= 0;
      for (i=0; i<ndl_n_bonds;i++) {
        if ((ndl_bond[i]->a1 == ndl_a)&& (ndl_bond[i]->a2 != prev_ndl_a) && 
        (nb_next_count < max_neighbors) && (ndl_atom[ndl_bond[i]->a2]->heavy)) {
          ndl_nb[nb_next_count] = ndl_bond[i]->a2;
          nb_next_count++;
        }
        if ((ndl_bond[i]->a2 == ndl_a)&& (ndl_bond[i]->a1 != prev_ndl_a) &&
        (nb_next_count < max_neighbors) && (ndl_atom[ndl_bond[i]->a1]->heavy)) {
          ndl_nb[nb_next_count] = ndl_bond[i]->a1;
          nb_next_count++;
        }
      }
//      ndl_nb = get_ndl_nextneighbors(ndl_a,prev_ndl_a);
      nb_next_count= 0;
      for (i=0; i<n_bonds;i++) {
        if ((bond[i]->a1 == hay_a)&& (bond[i]->a2 !=prev_hay_a) && 
        (nb_next_count < max_neighbors) && (atom[bond[i]->a2]->heavy)) {
          hay_nb[nb_next_count] = bond[i]->a2;
          nb_next_count++;
        }
        if ((bond[i]->a2 == hay_a)&& (bond[i]->a1 != prev_hay_a) &&
        (nb_next_count < max_neighbors) && (atom[bond[i]->a1]->heavy)) {
          hay_nb[nb_next_count] = bond[i]->a1;
          nb_next_count++;
        }
      }
//      hay_nb = get_nextneighbors(hay_a,prev_hay_a);
  }          
  // v0.3o: mark all neighbor atoms as "visited"
  for (i = 0;i< ndl_n_nb;i++) ndl_atom[(ndl_nb[i])]->tag = false;
  // now that the neighbor-arrays are filled, get all
  // combinations of matches recursively;
  // first, initialize the match matrix
  memset(mm,0,sizeof(MATCHMATRIX));    // new in v0.2i
  // make sure there are not too many neighbors (max. max_neighbors)  
  if ((ndl_n_nb > max_neighbors) || (n_nb > max_neighbors)) return false;
  // check if matchpath is not already filled up
  if (matchpath_pos(0,ndl_mp) == max_matchpath_length) return false;
  // next, check which chain of the needle matches which chain of the haystack 
  for (i=0;i<ndl_n_nb;i++) {
      emptyline = true;
      next_ndl_a = ndl_nb[i];
      for (j =0;j<n_nb;j++) {
          next_hay_a = hay_nb[j];
          ndl_mp[ndl_mp_len] = next_ndl_a;
          hay_mp[hay_mp_len] = next_hay_a;
//printf("gethere,%d,%d,%d,%d,%d\n", ndl_a, hay_a, ndl_mp[1], ndl_mp[2],res);
          if (is_matching(ndl_mp,hay_mp)) {             // recursive function call
              recursion_depth++;     // new in v0.3p: limit the recursion depth
              if ((max_match_recursion_depth > 0) && 
              (recursion_depth > max_match_recursion_depth)) {
                  if (opt_verbose) {
                    printf("Warning: max. number of match recursions (%d) reached, reverting to non-exhaustive match", max_match_recursion_depth);}
                  res= true;
                  if (use_gmm && valid_gmm) *(*(gmm+ndl_a)+hay_a) = true;  // v0.4a
                  return res;
              }                    // end of v0.3p recursion depth checking 
              emptyline = false;
              mm[i][j] = true;
          }
      }
      // if a needle substituent does not match any of the haystack substituents,
      // stop any further work immediately
      if (emptyline) { return false;}
  } 

  // finally, check the content of the matrix
  res = matrix_OK(mm,ndl_n_nb,n_nb);
  // optional: chirality check
  if (res && (rs_search || (use_gmm && valid_gmm) || ndl_atom[ndl_a]->stereo_care) 
    && ndl_maybe_chiral(ndl_a)) {
      // first, we have to clean up the match matrix in order to remove
      // "impossible" multiple matches (new in v0.3h)
    for (i=0;i<3;i++) {
      for (j=0;j<max_neighbors;j++) {  // haystack dimension
        n_hits= 0;
        l    = 0;
        for (k =0;k<max_neighbors;k++) {  // needle dimension
            if (mm[k][j]) {
              n_hits++;
              l = k;
            }
        }
        if (n_hits == 1) {  // a unique match ==> kick out any other match at this pos.
          for (m =0;m< max_neighbors;m++) {
            if (m != j) {
              if (mm[l][m]) {
                if (use_gmm && valid_gmm) *(*(gmm+ndl_nb[l])+hay_nb[m])=false;              
              }
              mm[l][m] = false;
            }
          }
        }
      }
    }
      // end of match matrix clean-up
      if (prev_ndl_a > 0) {
          n_singlehits = 1;
          ncp[1] = prev_ndl_a;
          hcp[1] = prev_hay_a;
      } else {
          n_singlehits = 0;
      }
      ncp[0] = ndl_a;
      hcp[0] = hay_a;  
      i = 0; l = 0;
      while ((n_singlehits < 3) && (i < 3)) {
          n_hits= 0;
          for (k =0;k<n_nb;k++) {
              if (mm[i][k]) {
                  n_hits++;
                  l = k;
              }
          }
          if (n_hits == 1) {
              n_singlehits++;
              ncp[n_singlehits] = ndl_nb[i];  
              hcp[n_singlehits] = hay_nb[l];
          } 
          i++;
      }
      if (n_singlehits == 3) {
          if (!chirality_OK(ncp,hcp)) {
              res = false; 
          } 
      }
  }
  if (res && use_gmm && valid_gmm) gmm[ndl_a][hay_a]=true;  // v0.4a
  return res;
}

bool quick_match() {  // added in v0.2c
bool  res;
int  i, ndl_chg, ndl_rad, ndl_iso;
char ndl_atype[4], ndl_el[3];
  if ((ez_search || rs_search) && (ndl_n_heavyatoms > 3)) return false;
  if (((ndl_n_atoms < 1) || (n_atoms < 1)) ||
     ((ndl_n_atoms > n_atoms) || (ndl_n_bonds > n_bonds))) return false; // just to be sure...
  res = true;
  for (i =1;i<=ndl_n_atoms;i++) {
      //if atom^[i].atype <> ndl_atom^[i].atype then res := false;    // changed in
    if (strcmp(atom[i]->element,ndl_atom[i]->element)) res = false;  // v0.2k
    if ((opt_chg) && (atom[i]->formal_charge != ndl_atom[i]->formal_charge)) res = false;  // v0.3o, v0.3p
    if (opt_iso && (atom[i]->nucleon_number != ndl_atom[i]->nucleon_number)) res = false;  // v0.3p
    if (opt_rad && (atom[i]->radical_type != ndl_atom[i]->radical_type)) res = false;  // v0.3p
  }
  if (ndl_n_bonds > 0) {
      for (i=0;i< ndl_n_bonds;i++) {
          if ((ndl_bond[i]->a1 != bond[i]->a1) ||
             (ndl_bond[i]->a2 != bond[i]->a2) ||
             (ndl_bond[i]->btype != bond[i]->btype)) res = false;
      }
  }
  // added in v0.2d: special case: needle contains only one heavy atom; refined in v0.3l, v0.3o
  if (ndl_n_heavyatoms == 1) {
      res = false;  // v0.3p
      // first, find out the element and atom type of the only heavy atom      
      for (i= 1;i<=ndl_n_atoms;i++) {
        if (ndl_atom[i]->heavy) {
            strcpy(ndl_atype,ndl_atom[i]->atype);
            strcpy(ndl_el,ndl_atom[i]->element);  // v0.3l
            ndl_chg   = ndl_atom[i]->formal_charge;  // v0.3l
            ndl_iso   = ndl_atom[i]->nucleon_number; // v0.3p
            ndl_rad   = ndl_atom[i]->radical_type;   // v0.3p
        }
      }
      for (i = 1;i<=n_atoms;i++) {
          if (!strcmp(atom[i]->atype,ndl_atom[i]->atype) && !strcmp(atom[i]->element,ndl_el) &&
             (ndl_chg == atom[i]->formal_charge))  res = true;  // v0.3l, v0.3o
          if (res && opt_iso) {
            if (atom[i]->nucleon_number != ndl_iso) res = false;
          }// v0.3p
          if (res && opt_rad) {
            if (atom[i]->radical_type != ndl_rad) res = false; // v0.3p
          }
      }
  }
  if (res && use_gmm && valid_gmm && opt_matchnum1) {
      for (i = 1;i<=ndl_n_atoms;i++) *(*(gmm+i)+i)= true;
  }
  return res;
}

int hetatom_count(int a) {
int  i,hac,nb_count,j;
NEIGHBOR_REC  nb;
  hac = 0;
  if ((a > 0) && (a <= n_atoms)) {
      nb_count = 0;
      for (j= 0;j< n_bonds;j++ ){ //get ra2's neighbors!
        if ((bond[j]->a1 == a) && (nb_count < max_neighbors) && (atom[bond[j]->a2]->heavy)) {
          nb[nb_count]=bond[j]->a2;
          nb_count++;
        }
        if ((bond[j]->a2== a)&& (nb_count < max_neighbors)&& (atom[bond[j]->a1]->heavy)) {
          nb[nb_count]=bond[j]->a1;
          nb_count++;
        }
      }
      if (atom[a]->neighbor_count > 0) {
        for (i=0;i<atom[a]->neighbor_count;i++) {
          if (strcmp(atom[(nb[i])]->element,"C ") && strcmp(atom[(nb[i])]->element,"H ") &&
            strcmp(atom[(nb[i])]->element,"D ") && strcmp(atom[(nb[i])]->element,"LP") &&
            strcmp(atom[(nb[i])]->element,"DU")) hac++;  // added 'D ' in v0.3n
          }
      }
  }
  return hac;
}

void perform_match() {
int   i, j;  //ndl_ref_atom : integer;  // since v0.3j as a global variable
int  ndl_n_nb, ndl_n_hc, n_nb, n_hc ;
bool  qm;  // v0.3l
  // check for NoStruct (0 atoms);  v0.3l
  if (n_atoms <1 || ndl_n_atoms<1) {
      matchresult = false;
      return;
  }
  // if we perform an exact match, needle and haystack must have
  // the same number of atoms, bonds, and rings
  if (opt_exact && opt_iso) {
      if ((n_heavyatoms != ndl_n_heavyatoms) || (n_heavybonds != ndl_n_heavybonds)) {
          matchresult = false;
          return;
      }
  }
  if (opt_exact && (!opt_iso) && (n_trueheavyatoms != ndl_n_trueheavyatoms)) {
      matchresult = false;
      return;
  }    
  // have a quick look if needle and haystack are identical molfiles
  if ((use_gmm == false) || (valid_gmm == false) || (opt_matchnum == false)) {
      qm = quick_match();  // v0.3l
      if (qm) {
          matchresult = true; 
          clear_ndl_atom_tags();  // v0.3o
          return;
      }
  }
  // if we have only one heavy atom and quick_match fails, return "false";  v0.3l
  if ((qm == false) && (ndl_n_heavyatoms == 1)) {
      matchresult = false;
      return;
  }
  ndl_n_nb = ndl_atom[ndl_ref_atom]->neighbor_count;
  ndl_n_hc = ndl_hetatom_count(ndl_ref_atom);
  matchresult = false;
  for (j = 0;j< max_matchpath_length;j++) {
      ndl_matchpath[j] = 0;
      hay_matchpath[j] = 0;
  }
  ndl_matchpath[0] = ndl_ref_atom;
  i = 0;
  found_untagged = false;
  while ((i <= n_atoms) && (matchresult == false)) {
      recursion_depth = 0;  // v0.3p
      n_nb = atom[i]->neighbor_count;
      n_hc = hetatom_count(i);
      if (((n_nb >= ndl_n_nb) && (n_hc >= ndl_n_hc))
         && (!(use_gmm && valid_gmm && tmp_tag[i]))) {
          if (use_gmm && valid_gmm) found_untagged = true;   // v0.4a
          hay_matchpath[0] = i;
          matchresult = is_matching(ndl_matchpath, hay_matchpath);
          if (use_gmm && valid_gmm) {
             if (gmm_collision()==true) {
                 matchresult = false;
                 memset(gmm,0,sizeof(GLOBAL_MATCHMATRIX));
             }
          }
          if (matchresult) atom[i]->tag = true;  // v0.3o; mark this fragment as matched
          tmp_tag[i] = true;    // v0.4a
      }
      i++;
  }
}

int count_gmm_matches() {  // v0.4a
int  i, j, n ;
  n = 0;
  for (i= 1;i<=ndl_n_atoms;i++) {
      for (j= 1;j<=n_atoms;j++) if (*(*(gmm+i)+j)) n++;
  }
  return n;
}

void init_gmm() {
int i,j;
  for (i= 0;i<max_ndl_gmmsize;i++) {
      for (j= 0;j<max_hay_gmmsize;j++) *(*(gmm+i)+j)=0;
  }
}
 
void init_gmm_total() {
int i,j;
  for (i= 0;i<max_ndl_gmmsize;i++) {
      for (j= 0;j<max_hay_gmmsize;j++) *(*(gmm_total+i)+j)=0;
  }
}

void cleanup_gmm() {  // v0.4a
int  i, j, k, l;
int  old_n_matches, new_n_matches;
int  n_hits, hit_pos, nb_count;
NEIGHBOR_REC  ndl_nb, hay_nb;
int  ndl_n_nb, hay_n_nb;
int  na, ha;
bool  found_neighbor, found_all_neighbors;
  new_n_matches = count_gmm_matches();
  while (new_n_matches != old_n_matches) {
    old_n_matches = new_n_matches;
    // first step: check for unique matches and remove other substituents from here
    for (i= 1;i<= ndl_n_atoms;i++) {
        n_hits = 0;
        for (j= 1;j<=n_atoms;j++) {
            if (*(*(gmm+i)+j)) 
                n_hits++;
                hit_pos = j;
            }
        }
        if (n_hits == 1) {
            for (k = 1; k<=ndl_n_atoms;k++) {
                if (k != i) *(*(gmm+k)+hit_pos) = false;
            }
        }
    }
    // second step: check for dangling neighbor atoms
    for (i= 1;i<=ndl_n_atoms;i++) {
        memset(ndl_nb, 0, sizeof(NEIGHBOR_REC)); 
        nb_count = 0;
        for (j= 0;j< ndl_n_bonds;j++ ){ //get ra2's neighbors!
          if ((ndl_bond[j]->a1 == i) && (nb_count < max_neighbors) && 
            (ndl_atom[ndl_bond[j]->a2]->heavy)) {
            ndl_nb[nb_count]=ndl_bond[j]->a2;
            nb_count++;
          }
          if ((ndl_bond[j]->a2== i) && (nb_count < max_neighbors) &&
            (ndl_atom[bond[j]->a1]->heavy)) {
            ndl_nb[nb_count]=ndl_bond[j]->a1;
            nb_count++;
          }
        }
        ndl_n_nb = ndl_atom[i]->neighbor_count;
        for (j = 1;j<=n_atoms;j++) {
          if (*(*(gmm+i)+j)) {
            memset(hay_nb, 0, sizeof(NEIGHBOR_REC)); 
            nb_count = 0;
            for (k= 0;k< n_bonds;k++ ){ //get ra2's neighbors!
              if ((bond[k]->a1 == j) && (nb_count < max_neighbors) &&
                (atom[bond[k]->a2]->heavy)) {
                hay_nb[nb_count]=bond[k]->a2;
                nb_count++;
              }
              if ((bond[k]->a2== j)&& (nb_count < max_neighbors)&& 
                (atom[bond[k]->a1]->heavy)) {
                hay_nb[nb_count]=bond[k]->a1;
                nb_count++;
              }
            }
                hay_n_nb = atom[j]->neighbor_count;
                if (ndl_n_nb > 0) {
                    found_all_neighbors = true;
                    for (k =0;k<ndl_n_nb;k++) {
                        found_neighbor = false;
                        for (l =0;l<hay_n_nb;l++) {
                            na = ndl_nb[k];
                            ha = hay_nb[l];
                            if (*(*(gmm+na)+ha)) {
                                found_neighbor = true;
                            }
                        }
                        if (found_neighbor == false) found_all_neighbors = false;
                    }
                    if (found_all_neighbors == false) *(*(gmm+i)+j) = false;
                }
            }          
        }
    new_n_matches = count_gmm_matches();
  }
}

void copy_gmm_to_total() {  // v0.4a
int  i, j ;
  for (i= 1;i<= ndl_n_atoms;i++) {
    for (j = 1;j<= n_atoms;j++) {
      if (*(*(gmm+i)+j)) *(*(gmm_total+i)+j) = true;
    }
  }
}

int count_tagged_ndl_heavyatoms() {
int  i, n;
  n = 0;
  if (ndl_n_atoms > 0) {
      for (i= 1;i<=ndl_n_atoms;i++) {
          if (ndl_atom[i]->heavy && ndl_atom[i]->tag) n++;
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

void show_usage() {
      printf("molmatch 2014--for MOLBASE\n");
      printf("Usage: molmatch [options] <needle> <haystack>\n");
      printf("  where options can be:\n");
      printf("    -a  check charges strictly\n");
      printf("    -i  check isotopes strictly\n");
      printf("    -d  check radicals strictly\n");
      printf("    -x  exact match\n");
      printf("    -s  standard in mode\n");
      printf("    -g  check trans-, cis- double bonds, which are not found by -x\n");
      printf("    -G  check geometry of chiral centers, with 1 and 6, represent up or down stereo\n");
      printf("    -S  strict match, comparison of atom and bond types\n");
      printf("    -F  fingerprint mode (1 haystack, multiple needles) with decimal output\n");
      printf("    -M  accept metal atoms as ring members\n");
    // the "debug" option (-D) remains undocumented
}

void write_needle_mol() {
int  i, j;
RINGPATH_TYPE  testring;
int  ring_size, a_prev, a_ref;
bool aromatic;
  printf("Molecule name: %s\n",ndl_molname);
  printf("atoms: %d bonds: %d rings: %d\n",ndl_n_atoms,ndl_n_bonds,ndl_n_rings);
  if (n_atoms < 1) return;
  if (n_bonds < 1) return;
  for (i = 1; i<=ndl_n_atoms; i++) {
    if (i <   10) printf(" ");
    if (i <  100) printf(" ");
    if (i < 1000) printf(" ");
    printf("%d %s %s %9.4f %9.4f %9.4f",i,ndl_atom[i]->element,ndl_atom[i]->atype,ndl_atom[i]->x,ndl_atom[i]->y,ndl_atom[i]->z);
    printf("  (%d heavy-atom neighbors, Hexp: %d, Htot: %d)",ndl_atom[i]->neighbor_count,ndl_atom[i]->Hexp,ndl_atom[i]->Htot);
    if (ndl_atom[i]->formal_charge != 0) printf("  charge: %d",ndl_atom[i]->formal_charge);
    if (ndl_atom[i]->arom) printf(" aromatic");
    printf("\n");
  }
  for (i=0; i<ndl_n_bonds;i++) {
      if (i <   10) printf(" ");
      if (i <  100) printf(" ");
      if (i < 1000) printf(" ");
      printf("%d %d %d %c",i,ndl_bond[i]->a1,ndl_bond[i]->a2,ndl_bond[i]->btype);
      if (ndl_bond[i]->ring_count > 0) printf(", contained in %d rings",ndl_bond[i]->ring_count);
      if (ndl_bond[i]->arom) printf(" (aromatic) ");
      printf("\n");
  }
  if (ndl_n_rings > 0) {
    for (i =0; i<ndl_n_rings; i++) {
      aromatic = true;
      memset(testring,0,sizeof(RINGPATH_TYPE));
      for (j=0;j<max_ringsize;j++) if (*(*(ndl_ring+i)+j)>0) testring[j]=*(*(ndl_ring+i)+j);
      ring_size = path_pos(0,testring);
      printf("ring %d: ",i);
      a_prev = testring[ring_size-1];
      for (j=0;j<ring_size;j++) {
        printf("%d ",testring[j]);
        a_ref= testring[j];
        if (!ndl_bond[get_ndl_bond(a_prev,a_ref)]->arom) aromatic = false;  // v0.3k
        a_prev = a_ref;
      }
      if (aromatic) printf(" (aromatic)");
      printf("\n");
    }
  }
}

int main (int argc, char** argv) {
  int i,c,n;
  init_globals();
  while (1)
  {
    int option_index = 0;
    static struct option long_options[] =
    {
      {"add", 1, 0, 0},
      {"append", 0, 0, 0},
      {"delete", 1, 0, 0},
      {"verbose", 0, 0, 0},
      {"create", 1, 0, 'c'},
      {"file", 1, 0, 0},{0, 0, 0, 0},
      {0, 0, 0, 0},{0, 0, 0, 0},{0, 0, 0, 0},{0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "axrsmidgGMFS", long_options, &option_index);
    if (c == -1) break;

          switch (c)
            {
            case 0:
              printf ("option %s", long_options[option_index].name);
              if (optarg)
                printf (" with arg %s", optarg);
              printf ("\n");
              break;
            case 'x':
              opt_exact=true;
              break;
            case 'S':
              opt_strict=true;
              break;
            case 's':
              opt_stdin=true;
              break;
            case 'a':
              opt_chg= true;
              break;
            case 'r':
              opt_rs=rs_ssr;
              break;
            case 'm':
              opt_molout=true;
//              fglang=1;
              break;
            case 'i':
              opt_iso=true;
              break;
            case 'd':
              opt_rad=true;
              break;
            case 'g':
              opt_geom=true;
              break;
            case 'G':
              opt_chiral=true;
              break;
            case 'M':
              opt_metalrings=true;
              break;
            case 'F':
              opt_fp=true;
              fpformat = fpf_decimal;
              break;
            case '?':
              break;
            default:
              printf ("?? getopt returned character code 0%o ??\n", c);
            }
       }
/*if (optind < argc)
  {
   printf ("non-option ARGV-elements: ");
   while (optind < argc)
   printf ("%s ", argv[optind++]);
   printf ("\n");
  }*/
if(argc<2) {show_usage(); return 0; }
if (opt_nochirality) opt_chiral = false;  // v0.4e
if (opt_geom) ez_search = true;  // v0.3d
if (opt_chiral) rs_search = true;  // v0.3f
if (opt_chiral && opt_strict && (opt_exact || opt_fp)) rs_strict= true; // new in v0.3j, v0.3m
if (opt_strict) use_gmm = true;  // v0.4b
if (opt_fp) {    // v0.3m
  opt_molout = false;
  opt_exact  = false;
  opt_matchnum  = false;  // v0.4a
  opt_matchnum1 = false;  // v0.4a
  use_gmm       = false;  // v0.4a
}
ringsearch_mode = opt_rs;
mol_OK=true;
if(!opt_stdin) { 
    ndl_molfilename=argv[2]; 
    hay_molfilename=argv[3];
    readfile(ndl_molfilename); 
} else {
    readfile(ndl_molfilename); 
}
if (ringsearch_mode == rs_sar) {max_vringsize = max_ringsize;} else {
                                   max_vringsize = ssr_vringsize;}
read_molfile();//ndl_in
count_neighbors();
if (!found_arominfo) {
  chk_ringbonds();
  if (ringsearch_mode == rs_ssr) remove_redundant_rings();
  if (n_rings == max_rings) {
    if (opt_verbose) printf("warning: max. number of rings exceeded, reverting to SSR search");
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
} else {
  if (!opt_xmdlout) normalize_ionic_bonds();  // v0.3k
  update_atypes();  // added in v0.2f
  update_Htotal();  // end v0.2b snippet
}
if (use_gmm) {
  for(n=0;n<max_ndl_gmmsize;n++) gmm[n]=(bool *)malloc(max_hay_gmmsize*sizeof(bool));
  if ( gmm==NULL) {printf("Error allocating memory for ring!");exit(1);}
  for(n=0;n<max_ndl_gmmsize;n++) gmm_total[n]=(bool *)malloc(max_hay_gmmsize*sizeof(bool));
  if ( gmm_total==NULL) {printf("Error allocating memory for ring!");exit(1);}
}

      // now transfer all data to the "needle" set of variables, except for "fingerprint" mode
if (!opt_fp) {   // v0.3m
    copy_mol_to_needle();
    chk_wildcard_rings();  // v0.3p
    set_ndl_atom_tags();   // v0.3o
//    if (opt_verbose) write_needle_mol();
    if (ndl_n_atoms > max_ndl_gmmsize) {valid_gmm = false;} else {valid_gmm = true;}  // v0.4a
    if (rs_strict) {
        ndl_ref_atom = find_ndl_ref_atom_cv();
    } else {
      ndl_ref_atom = find_ndl_ref_atom();
    }
} else {
    copy_mol_to_tmp();  // v0.3m
    if (opt_verbose) printf("1st molecule stored in buffer: %s\n",tmp_molname);
    valid_gmm = false;  // v0.4a global match matrix is not needed in fingerprint mode
}
      // next, read the "haystack" file and process it
li = 1;
mol_count = 0;
fpdecimal = 0;  // v0.3m
fpindex   = 0;  // v0.3m
mol_in_queue = true;
fin=0;
while (mol_in_queue == true) {    
    // new in v0.3i: reset ringsearch_mode to its initial value
    // for each new molecule
    ringsearch_mode = opt_rs;
    if (ringsearch_mode == rs_sar) { max_vringsize = max_ringsize;} else {
                                     max_vringsize = ssr_vringsize;}  // v0.3n (was: 10)
    readfile(hay_molfilename);
if (ringsearch_mode == rs_sar) {max_vringsize = max_ringsize;} else {
                                   max_vringsize = ssr_vringsize;}
//    li = 1;
    found_arominfo = false;  // added in v0.2b
    mol_OK = true;           // added in v0.2i
    read_molfile();    
    mol_count++;
    fpindex++;
    count_neighbors();
 //if (not mol_OK) or (n_atoms < 1) then writeln(mol_count,':no valid structure found') else
    if (((! mol_OK) || (n_atoms < 1)) && (!(opt_fp && (fpformat == fpf_decimal)))) {
      printf("%d:F\n",mol_count); 
    } else {   // v0.3l
      if (opt_exact && !ndl_querymol && ((n_Ctot != ndl_n_Ctot) || (n_Otot != ndl_n_Otot) || (n_Ntot != ndl_n_Ntot))) {       // new in v0.3g
        if ((!opt_molout) && (!(opt_fp && (fpformat == fpf_decimal)))) { 
          printf("%d:F\n",mol_count);}
      } else {
        if ((!found_arominfo) || (opt_strict && tmfmismatch)) { n_rings=0;// temporary try, just to make sure n_rings is not tested yet, added in v0.3mmj
          chk_ringbonds();
          if (ringsearch_mode == rs_ssr) remove_redundant_rings();
          if (n_rings == max_rings) {
            if (opt_verbose) printf("Warning: max. number of rings reached, reverting to SSR search");
            ringsearch_mode = rs_ssr;
            clear_rings();
            max_vringsize = ssr_vringsize;  // v0.3n (was: 10)
            chk_ringbonds();
            remove_redundant_rings();
          }
          update_ringcount();
          update_atypes();
          update_Htotal();   // added in v0.3
          chk_arom();
          if (ringsearch_mode == rs_ssr){
            do{
              prev_n_ar = count_aromatic_rings();
              chk_arom();
              n_ar = count_aromatic_rings();
            } while((prev_n_ar - n_ar) != 0);
          } 
        } else { 
          if (opt_strict) update_atypes();  // added in v0.2f
          update_Htotal();
        }
//        init_molstat(ndl_molstat);
        if (normalize_ionic_bonds()) update_atypes();   // new in v0.3k, modified in v0.3m
      // if in "fingerprint mode", exchange needle and haystack
        if (opt_fp) {
          zap_needle(); 
          copy_mol_to_needle();
          chk_wildcard_rings(); // v0.3p
          zap_molecule();
          copy_tmp_to_mol();
//          if (opt_verbose) write_needle_mol();
          if (rs_strict) {      // v0.3j
            ndl_ref_atom= find_ndl_ref_atom_cv();
          } else { 
            ndl_ref_atom = find_ndl_ref_atom();}
//          
        }  //printf("%d\n",mol_count);  // v0.3m
          // now that we have both molecules, perform the comparison
          // v0.3o: takes care of disconnected fragment...
        update_atypes_quick();  // v0.4b
        if ((opt_verbose) && !opt_fp) write_mol();  // v0.4b
        clear_atom_tags();
        overall_match = false;  // v0.4a
//              fillchar(tmp_tag,sizeof(tmp_tag),false);   // v0.4a
        found_untagged = false;
//        if (use_gmm && valid_gmm) memset(gmm_total,0,max_ndl_gmmsize*max_hay_gmmsize*sizeof(bool));  // v0.4a
         n_matches = 0;  // v0.4b        // this repeat...until is new in v0.4a
         do {// v0.4b
           set_ndl_atom_tags();
           matchsummary = true; 
           if (use_gmm && valid_gmm) init_gmm();  // v0.4a
           perform_match();
           matchsummary = matchresult;
           if ((count_tagged_ndl_heavyatoms() > 0) && (matchsummary == true)) {
             do {  
               if (rs_strict) {
                 ndl_ref_atom = find_ndl_ref_atom_cv();
               } else {ndl_ref_atom = find_ndl_ref_atom();}
               perform_match();
               if (matchresult == false) matchsummary = false;
             } while ((count_tagged_ndl_heavyatoms() != 0) && (matchsummary == true)); 
           }                  // end of disconnected-fragment matching (v0.3o)
           if (matchsummary==true) {
             overall_match = true;
             n_matches++;
             if (use_gmm && valid_gmm) {
               cleanup_gmm();
               init_gmm_total();
               copy_gmm_to_total();
//               write_matches();   // v0.4b
             }
           }
       } while ((use_gmm == true) && (valid_gmm == true) && (n_matches < max_n_matches) &&  
           (opt_matchnum1 == false) && (found_untagged == true));
         if (overall_match == true) {  // v0.3o, v0.4a
           if (opt_molout) {
             for (i=0; i<molbufindex;i++) printf("%s\n",molbuf[i]);
//             if (opt_verbose && use_gmm) write_gmm();   // v0.4a
           } else {
             if (!opt_fp) {
               printf("%d:T\n",mol_count);   // v0.4a
             } else {
     //if (ndl_n_heavyatoms = n_heavyatoms) and (ndl_n_heavybonds = n_heavybonds) then 
            //   fp_exacthit := true else fp_exacthit := false;
               if ((ndl_n_atoms == n_atoms) && (ndl_n_bonds == n_bonds)) { 
// v0.4b; also explicit H must match!
                 fp_exacthit = true;} else {fp_exacthit = false; }
               if (fp_exacthit) fp_exactblock = true;
               if (fpformat == fpf_decimal) {
                 fpincrement = 1;
                 for (i =1; i<=fpindex;i++) fpincrement=fpincrement<<1;
                 fpdecimal=fpdecimal+fpincrement;                        
               }                                   
           }
         }
       } else {
         if (!(opt_molout || (opt_fp && (fpformat == fpf_decimal)))) {
           printf("%d:F\n",mol_count);} 
       }
       if ((opt_fp && (fpformat == fpf_decimal)) && (fpindex == fp_blocksize)) {
         if (fp_exactblock) fpdecimal = fpdecimal + 1;
         printf("%ld\n",fpdecimal);
         fpindex   = 0;
         fpdecimal = 0;
         fp_exactblock = false;
       }
       zap_molecule();
       molbufindex = 0;
     }
   }  // mol_OK
 }
 if ((opt_fp && (fpformat == fpf_decimal)) && (fpindex > 0)) {
   if (fp_exactblock) fpdecimal = fpdecimal + 1;
   printf("%ld\n",fpdecimal);
 }
 zap_needle();
 if (use_gmm) {
    if(gmm!=NULL){
          free(*gmm);
          *gmm=NULL;}
    if(gmm_total!=NULL){
          free(*gmm_total);
          *gmm_total=NULL;}
 }  
return 0;
}
