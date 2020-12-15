
#include <stdio.h>
#include "InfoStream.h"
#include "DumpBench.h"
#include "SimParameters.h"
#include "ComputeNonbondedUtil.h"
#include "LJTable.h"
#include "Molecule.h"
#include "Node.h"
#include "PatchMap.h"
#include "HomePatch.h"
#include "NamdState.h"
#include "ComputeMap.h"

inline void dump_param(FILE *file, const char *name, int value) {
  fprintf(file,"%s %d\n",name,value);
}

inline void dump_param(FILE *file, const char *name, double value) {
  fprintf(file,"%s %g\n",name,value);
}

inline void dump_param(FILE *file, const char *name, Vector value) {
  fprintf(file,"%s %f %f %f\n",name,value.x,value.y,value.z);
}

int dumpbench(FILE *file) {

  Node *node = Node::Object();

  fprintf(file,"SIMPARAMETERS_BEGIN\n");

  SimParameters *simParams = node->simParameters;

#define SIMPARAM(T,N,V) dump_param(file,#N,simParams->N)
#include "DumpBenchParams.h"
#undef SIMPARAM

  fprintf(file,"SIMPARAMETERS_END\n");

  fprintf(file,"LJTABLE_BEGIN\n");

  const LJTable *ljTable = ComputeNonbondedUtil::ljTable;

  int table_dim = ljTable->get_table_dim();
  fprintf(file,"%d\n",table_dim);

  const LJTable::TableEntry *table = ljTable->get_table();
  int i,j;
  for ( i=0; i < table_dim; ++i) {
    for ( j=i; j < table_dim; ++j)
    {
      const LJTable::TableEntry *curij = &(table[2*(i*table_dim+j)]);
      fprintf(file,"%g %g %g %g\n",curij->A,curij->B,
				(curij+1)->A,(curij+1)->B);
    }
  }

  fprintf(file,"LJTABLE_END\n");

  fprintf(file,"MOLECULE_BEGIN\n");

  const Molecule *mol = node->molecule;

  fprintf(file,"%d %d\n",mol->numAtoms,mol->numCalcExclusions);
 
  for ( i=0; i<mol->numAtoms; ++i) {
    int vdw = mol->atomvdwtype(i);
    #ifdef MEM_OPT_VERSION
    Index exclIdx = mol->getAtomExclSigId(i);
    const ExclusionCheck *excl = mol->get_excl_check_for_idx(exclIdx);
    //if(excl->flags==NULL || excl->flags == (char *)-1){
    //    fprintf(file,"%d: %d ====\n",i, vdw);
    //    continue;
    //}
    int min = i+excl->min;
    int max = i+excl->max;
    #else
    const ExclusionCheck *excl = mol->get_excl_check_for_atom(i);
    //if(excl->flags==NULL || excl->flags == (char *)-1){
    //    fprintf(file,"%d: %d ====\n",i, vdw);
    //    continue;
    //}
    int min = excl->min;
    int max = excl->max;
    #endif
    // fprintf(file,"%d: %d %d %d |",i,vdw,min,max);
    fprintf(file,"%d %d %d",vdw,min,max);
    if ( min <= max ) {
      int s = max - min + 1;
      const char *f = excl->flags;
      for ( int k=0; k<s; ++k ) {
        int fk = f[k];
        fprintf(file," %d",fk);
      }
    }
    fprintf(file,"\n");
  }

  fprintf(file,"MOLECULE_END\n");

#if 0
  fprintf(file, "BONDS_BEGIN\n");
  fprintf(file, "%d %d\n", mol->numBonds, mol->numCalcBonds);
#ifdef MEM_OPT_VERSION
  for(i=0; i<mol->numAtoms; i++){
      int sigId = node->molecule->getAtomSigId(i);
      AtomSignature *sig = &(mol->atomSigPool[sigId]);
      if(sig->bondCnt==0){
          fprintf(file, "%d: ===\n", i);
          continue;
      }
      fprintf(file, "%d:", i);
      for(j=0; j<sig->bondCnt; j++){
          fprintf(file, " (%d | %d)", (sig->bondSigs[j]).offset[0], sig->bondSigs[j].tupleParamType);
      }
      fprintf(file, "\n");
  }
#else
  for(i=0; i<mol->numAtoms; i++){      
      int *p = node->molecule->get_bonds_for_atom(i);      
      if(*p==-1){
          fprintf(file, "%d: ===\n", i);
          continue;
      }
      fprintf(file, "%d:", i);
      for(; *p!=-1;p++){
          Bond *t = mol->get_bond(*p);
          fprintf(file, " (%d | %d)", t->atom2-i, t->bond_type);
      }
      fprintf(file, "\n");
  }
#endif
  fprintf(file, "BONDS_END\n");

  fprintf(file, "ANGLES_BEGIN\n");
  fprintf(file, "%d %d\n", mol->numAngles, mol->numCalcAngles);
#ifdef MEM_OPT_VERSION
  for(i=0; i<mol->numAtoms; i++){
      int sigId = node->molecule->getAtomSigId(i);
      AtomSignature *sig = &(mol->atomSigPool[sigId]);
      if(sig->angleCnt==0){
          fprintf(file, "%d: ===\n", i);
          continue;
      }
      fprintf(file, "%d:", i);
      for(j=0; j<sig->angleCnt; j++){
          int offset0 = (sig->angleSigs[j]).offset[0];
          int offset1 = (sig->angleSigs[j]).offset[1];
          fprintf(file, " (%d, %d | %d)", offset0, offset1, sig->angleSigs[j].tupleParamType);
      }
      fprintf(file, "\n");
  }
#else
  for(i=0; i<mol->numAtoms; i++){      
      int *p = node->molecule->get_angles_for_atom(i);      
      if(*p==-1){
          fprintf(file, "%d: ===\n", i);
          continue;
      }
      fprintf(file, "%d:", i);
      for(; *p!=-1;p++){
          Angle *t = mol->get_angle(*p);
          int offset0 = t->atom2 - i;
          int offset1 = t->atom3 - i;
          fprintf(file, " (%d, %d | %d)", offset0, offset1, t->angle_type);
      }
      fprintf(file, "\n");
  }
#endif
  fprintf(file, "ANGLES_END\n");

  fprintf(file, "DIHEDRALS_BEGIN\n");
  fprintf(file, "%d %d\n", mol->numDihedrals, mol->numCalcDihedrals);
#ifdef MEM_OPT_VERSION
  for(i=0; i<mol->numAtoms; i++){
      int sigId = node->molecule->getAtomSigId(i);
      AtomSignature *sig = &(mol->atomSigPool[sigId]);
      if(sig->dihedralCnt==0){
          fprintf(file, "%d: ===\n", i);
          continue;
      }
      fprintf(file, "%d:", i);
      for(j=0; j<sig->dihedralCnt; j++){
          int offset0 = (sig->dihedralSigs[j]).offset[0];
          int offset1 = (sig->dihedralSigs[j]).offset[1];
          int offset2 = (sig->dihedralSigs[j]).offset[2];
          fprintf(file, " (%d, %d, %d | %d)", offset0, offset1, offset2, sig->dihedralSigs[j].tupleParamType);
      }
      fprintf(file, "\n");
  }
#else
  for(i=0; i<mol->numAtoms; i++){      
      int *p = node->molecule->get_dihedrals_for_atom(i);      
      if(*p==-1){
          fprintf(file, "%d: ===\n", i);
          continue;
      }
      fprintf(file, "%d:", i);
      for(; *p!=-1;p++){
          Dihedral *t = mol->get_dihedral(*p);
          int offset0 = t->atom2 - i;
          int offset1 = t->atom3 - i;
          int offset2 = t->atom4 - i;
          fprintf(file, " (%d, %d, %d | %d)", offset0, offset1, offset2, t->dihedral_type);
      }
      fprintf(file, "\n");
  }
#endif
  fprintf(file, "DIHEDRALS_END\n");

  fprintf(file, "IMPROPERS_BEGIN\n");
  fprintf(file, "%d %d\n", mol->numImpropers, mol->numCalcImpropers);
#ifdef MEM_OPT_VERSION
  for(i=0; i<mol->numAtoms; i++){
      int sigId = node->molecule->getAtomSigId(i);
      AtomSignature *sig = &(mol->atomSigPool[sigId]);
      if(sig->improperCnt==0){
          fprintf(file, "%d: ===\n", i);
          continue;
      }
      fprintf(file, "%d:", i);
      for(j=0; j<sig->improperCnt; j++){
          int offset0 = (sig->improperSigs[j]).offset[0];
          int offset1 = (sig->improperSigs[j]).offset[1];
          int offset2 = (sig->improperSigs[j]).offset[2];
          fprintf(file, " (%d, %d, %d | %d)", offset0, offset1, offset2, sig->improperSigs[j].tupleParamType);
      }
      fprintf(file, "\n");
  }
#else
  for(i=0; i<mol->numAtoms; i++){      
      int *p = node->molecule->get_impropers_for_atom(i);      
      if(*p==-1){
          fprintf(file, "%d: ===\n", i);
          continue;
      }
      fprintf(file, "%d:", i);
      for(; *p!=-1;p++){
          Improper *t = mol->get_improper(*p);
          int offset0 = t->atom2 - i;
          int offset1 = t->atom3 - i;
          int offset2 = t->atom4 - i;
          fprintf(file, " (%d, %d, %d | %d)", offset0, offset1, offset2, t->improper_type);
      }
      fprintf(file, "\n");
  }
#endif
  fprintf(file, "IMPROPERS_END\n");
#endif

  fprintf(file,"PATCHLIST_BEGIN\n");

  PatchMap *patchMap = PatchMap::Object();
  int numPatches = patchMap->numPatches();
  fprintf(file,"%d\n",numPatches);

  for ( i=0; i<numPatches; ++i) {
    HomePatch *patch = patchMap->homePatch(i);
    fprintf(file,"PATCH_BEGIN\n");
    int numAtoms = patch->getNumAtoms();
    fprintf(file,"%d\n",numAtoms);
    FullAtomList &atoms = patch->getAtomList();
    for ( j=0; j<numAtoms; ++j) {
      FullAtom &a = atoms[j];
      double x,y,z,q;
      int id,hgs,ngia,af,gf,part;
      x = a.position.x;
      y = a.position.y;
      z = a.position.z;
      q = a.charge;
      id = a.id;
      hgs = a.hydrogenGroupSize;
      ngia = ( hgs != 1 && a.nonbondedGroupSize == 1 ) ? 1 : 0;
      af = a.atomFixed;
      gf = a.groupFixed;
      part = a.partition;
      fprintf(file,"%f %f %f %f %d %d %d %d %d %d\n",
        x,y,z,q,id,hgs,ngia,af,gf,part);
    }
    fprintf(file,"PATCH_END\n");
  }

  fprintf(file,"PATCHLIST_END\n");

  fprintf(file,"COMPUTEPAIR_BEGIN\n");

  ComputeMap *computeMap = ComputeMap::Object();
  int numComputes = computeMap->numComputes();
  int numPairComputes = 0;
  for ( i=0; i<numComputes; ++i) {
    if ( computeMap->type(i) == computeNonbondedPairType
         && computeMap->partition(i) == 0 ) ++numPairComputes;
  }
  fprintf(file,"%d\n",numPairComputes);
  for ( i=0; i<numComputes; ++i) {
    if ( computeMap->type(i) == computeNonbondedPairType
         && computeMap->partition(i) == 0 ) {
      int pid1 = computeMap->pid(i,0);
      int trans1 = computeMap->trans(i,0);
      int pid2 = computeMap->pid(i,1);
      int trans2 = computeMap->trans(i,1);
      fprintf(file,"%d %d %d %d\n",pid1,trans1,pid2,trans2);
    }
  }

  fprintf(file,"COMPUTEPAIR_END\n");

  return 0;
}

