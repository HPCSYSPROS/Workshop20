/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "Node.h"
#include "Molecule.h"
#include "NamdTypes.h"

#include "SimParameters.h"
#include "GlobalMasterTMD.h"
#include "PDB.h"
#include "PDBData.h"
#include "fitrms.h"
//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"
//using namespace __gnu_cxx;
using namespace std;

class Matrix4TMD {
  BigReal mat[16];                               ///< the matrix itself
public:
  Matrix4TMD(void) { identity(); }
  Matrix4TMD(const BigReal *m)  { memcpy(mat, m, 16*sizeof(BigReal)); }
  void multpoint(BigReal point[3]) const {
    BigReal tmp[3];
    BigReal itmp3 = 1.0f / (point[0]*mat[3] + point[1]*mat[7] +
                            point[2]*mat[11] + mat[15]);
    tmp[0] = itmp3*point[0];
    tmp[1] = itmp3*point[1];
    tmp[2] = itmp3*point[2];
    point[0]=tmp[0]*mat[0] + tmp[1]*mat[4] + tmp[2]*mat[ 8] + itmp3*mat[12];
    point[1]=tmp[0]*mat[1] + tmp[1]*mat[5] + tmp[2]*mat[ 9] + itmp3*mat[13];
    point[2]=tmp[0]*mat[2] + tmp[1]*mat[6] + tmp[2]*mat[10] + itmp3*mat[14];
  }

  void identity() {
    memset(mat, 0, 16*sizeof(BigReal));
    mat[0]=1.0f;
    mat[5]=1.0f;
    mat[10]=1.0f;
    mat[15]=1.0f;
  }
  void transpose() {
    BigReal tmp[16];
    int i,j;
    for(i=0;i<4;i++) {
      for(j=0;j<4;j++) {
        tmp[4*i+j] = mat[i+4*j];
      }
    }
    for(i=0;i<16;i++) mat[i] = tmp[i];
  }
  /// premultiply the matrix by the given matrix, this->other * this
  void multmatrix(const Matrix4TMD &m) {
    BigReal tmp[4];
    for (int j=0; j<4; j++) {
      tmp[0] = mat[j];
      tmp[1] = mat[4+j];
      tmp[2] = mat[8+j]; 
      tmp[3] = mat[12+j];
      for (int i=0; i<4; i++) {
        mat[4*i+j] = m.mat[4*i]*tmp[0] + m.mat[4*i+1]*tmp[1] +
          m.mat[4*i+2]*tmp[2] + m.mat[4*i+3]*tmp[3]; 
      }
    } 
  }
  void translate(BigReal x, BigReal y, BigReal z) {
    Matrix4TMD m;		
    m.mat[12] = x;
    m.mat[13] = y;
    m.mat[14] = z;
    multmatrix(m);
  }
  void translate(BigReal d[3]) { translate(d[0], d[1], d[2]); }
};

GlobalMasterTMD::GlobalMasterTMD() {
  DebugM(3,"initialize called\n");
  SimParameters *params = Node::Object()->simParameters;
  outputFreq = params->TMDOutputFreq;
  K = params->TMDk;
/*  if (params->TMDInitialRMSD < 0){
    initialRMS = -1; // get from first coordinates
  }
  else */
    initialRMS = params->TMDInitialRMSD;
  finalRMS = params->TMDFinalRMSD;
  
  currentStep = params->firstTimestep;
  firstStep = params->TMDFirstStep;
  lastStep = params->TMDLastStep;
  qDiffRMSD=params->TMDDiffRMSD;
  altloc = 0;
  target = 0;
  target2 = 0;
  weight = 0;
  if (qDiffRMSD) parseAtoms(params->TMDFile2,Node::Object()->molecule->numAtoms, 1);
  parseAtoms(params->TMDFile,Node::Object()->molecule->numAtoms, 0);

  //iterate through all domains to see if altloc is used
  map <int, vector<int> >::iterator it;
  for (it = dmap.begin(); it != dmap.end(); ++it){  
    int refcount = 0;
    int biascount = 0;
    for(int i = 0; i<it->second.size(); i++){
      char aloc = altloc[it->second[i]];
      if ( aloc & 1 ) ++biascount;
      if ( aloc & 2 ) ++refcount;
    }
    altlocmap[it->first] = ( refcount ? 1 : 0 );
    if ( ! refcount ) refcount = it->second.size();
    iout << iINFO << "TMD domain " << it->first <<
      " has " << it->second.size() << " atoms " <<
      refcount << " fitted " << biascount << " biased\n" << endi;
  }

 // k /= numTMDatoms;
  iout << iINFO << numTMDatoms << " TMD ATOMS\n" << endi;
  DebugM(1,"done with initialize\n");
}

void GlobalMasterTMD::parseAtoms(const char *file, int numTotalAtoms, bool isTwo) {
  DebugM(3,"parseAtoms called\n");
  PDB tmdpdb(file);
  numatoms = tmdpdb.num_atoms();
  if (numatoms < 1) 
    NAMD_die("No atoms found in TMDFile\n");
  if (numatoms < numTotalAtoms)
    iout << iWARN << "The number of atoms in TMDFile is less than the total number of atoms in the structure.\n" << endi;
  if (numatoms > numTotalAtoms)
    NAMD_die("The number of atoms in TMDFile must not exceed the total number of atoms in the structure!");
  if ( modifyRequestedAtoms().size() )
    NAMD_bug("GlobalMasterTMD::parseAtoms() modifyRequestedAtoms() not empty");

  numTMDatoms = 0;
  if(isTwo){
    target2 = new BigReal[3*numatoms];
    atompos2 = new Vector[numatoms];
    tmdpdb.get_all_positions(atompos2);

  }
  else{
    target = new BigReal[3*numatoms];
    atompos = new Vector[numatoms];
    tmdpdb.get_all_positions(atompos);
  }
  if ( ! altloc ) altloc = new char[numatoms];
  int i;
  for (i=0; i<numatoms; i++) {
#ifdef MEM_OPT_VERSION
    PDBCoreData *atom = tmdpdb.atom(i);
    char aloc = tmdpdb.alternatelocation(i);
#else
    PDBAtom *atom = tmdpdb.atom(i); // get an atom from the file
    char aloc = atom->alternatelocation()[0];
#endif
    if ( aloc ) aloc -= '0';
    if ( aloc ) aloc = 2;  // 2 bit == reference
    if ( atom->occupancy() ) aloc |= 1;  // 1 bit == biased
    altloc[i] = aloc;
    if ( aloc ) {
      if(isTwo){
        target2[3*numTMDatoms  ] = atompos2[i].x;
        target2[3*numTMDatoms+1] = atompos2[i].y;
        target2[3*numTMDatoms+2] = atompos2[i].z;
      }
      else{
        target[3*numTMDatoms  ] = atompos[i].x;
        target[3*numTMDatoms+1] = atompos[i].y;
        target[3*numTMDatoms+2] = atompos[i].z;
 //       aidmap[i] = numTMDatoms++;
        numTMDatoms++;
        // add the atom to the list
        modifyRequestedAtoms().add(i);
        if(!K){ kmap[i] = atom->occupancy();}
        //check to see if domain is already in the map
        map <int, vector<int> >::iterator it = dmap.find((int)atom->temperaturefactor());
        if (it != dmap.end()){
          it->second.push_back(i); //add atomid to vector in proper existing domain
        }
        else{
           dmap[(int)atom->temperaturefactor()] = vector <int> (1,i); //create new domain with atomid
        }
      }
    }
  }
  
  DebugM(1,"done with parseAtoms\n");
}

//recreates target array for domain selected
void GlobalMasterTMD::NewTarget(int domain)
{
  map <int, vector<int> >::iterator it = dmap.find(domain);
  //target_aid = new int[it->second.size()];
  delete [] target;
  target = new BigReal [3*it->second.size()];
  for(int i = 0; i<it->second.size(); i++){
   target[3*i  ] = atompos[it->second[i]].x;
   target[3*i+1] = atompos[it->second[i]].y;
   target[3*i+2] = atompos[it->second[i]].z; 
  }
  if(qDiffRMSD){
   delete [] target2;
   target2 = new BigReal [3*it->second.size()];
   for(int i = 0; i<it->second.size(); i++){
     target2[3*i  ] = atompos2[it->second[i]].x;
     target2[3*i+1] = atompos2[it->second[i]].y;
     target2[3*i+2] = atompos2[it->second[i]].z; 
   }   
  }
  delete [] weight;  weight = 0;
  if ( altlocmap.find(domain)->second ) {
    weight = new BigReal [it->second.size()];
    for(int i = 0; i<it->second.size(); i++){
      weight[i] = ( (altloc[it->second[i]] & 2) ? 1.0 : 0.0 );
    }
  }
//   target_aid[i] = it->second[i];
//   aidmap[it->second[i]] = i; 

}
GlobalMasterTMD::~GlobalMasterTMD() { 
  delete [] target;
 // delete [] target_aid;
//  delete [] aidmap;
  delete [] atompos;
  delete [] atompos2;
  delete [] altloc;
  delete [] target2;
}
void GlobalMasterTMD::calculate() {
  // have to reset my forces every time.  
  modifyAppliedForces().resize(0);
  modifyForcedAtoms().resize(0);
  
  // see if TMD should be active
  if (currentStep < firstStep || currentStep >= lastStep) {
    currentStep++;
    return;
  }
  
  map<int, Position> positions;
  AtomIDList::const_iterator a_i = getAtomIdBegin();
  AtomIDList::const_iterator a_e = getAtomIdEnd();
  PositionList::const_iterator p_i = getAtomPositionBegin();

  //create mapping of positions now to avoid
  //going through these iterators for each domain
  for ( ; a_i != a_e; ++a_i, ++p_i ){
    positions[*a_i] = *p_i;
  }

  //iterate through all domains
  map <int, vector<int> >::iterator it;
  for (it = dmap.begin(); it != dmap.end(); ++it){  
  NewTarget(it->first);  //set new target
  // fetch the current coordinates
  BigReal *curpos = new BigReal[3*(it->second.size())];
  for (int i = 0; i < it->second.size(); i++){
    //int ind = 3*aidmap[it->second[i]];
    curpos[3*i  ] = positions[it->second[i]].x;
    curpos[3*i+1] = positions[it->second[i]].y;
    curpos[3*i+2] = positions[it->second[i]].z;
  }
  BigReal *curpos2;
if(qDiffRMSD){
  curpos2 = new BigReal[3*(it->second.size())];
  for (int i = 0; i < it->second.size(); i++){
    //int ind = 3*aidmap[it->second[i]];
    curpos2[3*i  ] = positions[it->second[i]].x;
    curpos2[3*i+1] = positions[it->second[i]].y;
    curpos2[3*i+2] = positions[it->second[i]].z;
  }
}
  // align target with current coordinates.   Uses same weight for all
  // atoms.  Maybe instead use weight from occupancy?
  BigReal ttt[16], pre[3], post[3];
  BigReal curRMS = MatrixFitRMS(it->second.size(), target, curpos, weight, ttt);
  // Compute targetRMS.
  if (initialRMS < 0) {
    initialRMS = curRMS;
  }

  BigReal curRMS0 = curRMS;
  BigReal curRMS1 =  1.; 
  BigReal ttt1[16]; 
  if(qDiffRMSD){
    curRMS1 = MatrixFitRMS(it->second.size(), target2, curpos2, weight, ttt1);
    curRMS = curRMS0 - curRMS1 ;
  }


  BigReal frac = (BigReal(currentStep-firstStep)) /
                            (lastStep-firstStep);
  
  BigReal targetRMS = initialRMS * (1-frac) + frac * finalRMS;

  
  BigReal maxforce2 = 0.;
//orig finalRMS < initialRMS...changed to <= when allowing initialRMS = 0
//qdiff part and the whole && section new finalRMS <=initialRMS
  if (((finalRMS <= initialRMS && targetRMS <= curRMS) ||
      (finalRMS >= initialRMS && targetRMS >= curRMS) ||
      qDiffRMSD) && (curRMS0 > 0. && curRMS1 > 0) ) {


    // compute transformation to align target structure with current structure
    // Might be more stable to instead align current positions with target,
    // although then we have to back-transform the forces.
    int j;
    for (j=0; j<3; j++) {
      post[j] = ttt[4*j+3];
      ttt[4*j+3] = 0;
      pre[j] = ttt[12+j];
      ttt[12+j] = 0;
    }
    Matrix4TMD result;
    result.translate(pre);
    result.multmatrix(Matrix4TMD(ttt));
    result.translate(post);
  
    // compute forces on each atom
    BigReal myrms = 0;
      for (int i=0; i<it->second.size(); i++) {
      BigReal k = 0.; 
      if(!K){
        k = kmap[it->second[i]];  
      } else if ( (! weight) || (altloc[it->second[i]] & 1) ) {
        k = K/it->second.size();  
      }
   //   BigReal prefac = k * (targetRMS / curRMS - 1); 
      BigReal prefac = k * (targetRMS - curRMS)/curRMS0;
      result.multpoint(target+3*i);
      BigReal dx = curpos[3*i  ] - target[3*i  ];
      BigReal dy = curpos[3*i+1] - target[3*i+1];
      BigReal dz = curpos[3*i+2] - target[3*i+2];
      
      BigReal fvec[3] = { dx, dy, dz };
      Vector force(fvec[0]*prefac, fvec[1]*prefac, fvec[2]*prefac);
      modifyForcedAtoms().add(it->second[i]);
      modifyAppliedForces().add(force);
      BigReal force2 = force.length2();
      if ( force2 > maxforce2 ) maxforce2 = force2;
    }

    if(qDiffRMSD){
       int j;
      for (j=0; j<3; j++) {
        post[j] = ttt1[4*j+3];
        ttt1[4*j+3] = 0;
        pre[j] = ttt1[12+j];
        ttt1[12+j] = 0;
      }
      Matrix4TMD result2;
      result2.identity();
      result2.translate(pre);
      result2.multmatrix(Matrix4TMD(ttt1));
      result2.translate(post);
    
      // compute forces on each atom
      BigReal myrms = 0;
        for (int i=0; i<it->second.size(); i++) {
        BigReal k = 0.;
        if(!K){
          k = kmap[it->second[i]];  
        } else if ( (! weight) || (altloc[it->second[i]] & 1) ) {
          k = K/it->second.size();  
        }
     //   BigReal prefac = k * (targetRMS / curRMS - 1); 
  //      BigReal prefac = k * (targetRMS - curRMS)/curRMS0;
          BigReal prefac1 = - k * (targetRMS  - curRMS)/curRMS1; // included with a negative term in the potential

        result2.multpoint(target2+3*i);
        BigReal dx = curpos2[3*i  ] - target2[3*i  ];
        BigReal dy = curpos2[3*i+1] - target2[3*i+1];
        BigReal dz = curpos2[3*i+2] - target2[3*i+2];
        BigReal fvec[3] = { dx, dy, dz };
        Vector force(fvec[0]*prefac1, fvec[1]*prefac1, fvec[2]*prefac1);
        modifyForcedAtoms().add(it->second[i]);
        modifyAppliedForces().add(force);
        BigReal force2 = force.length2();
        if ( force2 > maxforce2 ) maxforce2 = force2;
      }   
    }
  }
  //delete [] target_aid;
  delete [] curpos;
  if(qDiffRMSD){delete [] curpos2;}
// write output if needed
  if (currentStep % outputFreq == 0) {
    iout << "TMD  " << currentStep << " Domain: "<< it->first << " " << targetRMS << ' ' << curRMS << '\n' << endi; //*it
    // iout << "TMD  " << currentStep << " " << targetRMS << ' ' << curRMS << ' ' << sqrt(maxforce2) << '\n' << endi;
  }
  }
/*  // write output if needed
  if (currentStep % outputFreq == 0) {
    iout << "TMD  " << currentStep << " " << targetRMS << ' ' << curRMS << '\n' << endi;
    // iout << "TMD  " << currentStep << " " << targetRMS << ' ' << curRMS << ' ' << sqrt(maxforce2) << '\n' << endi;
  }*/
  currentStep++;
 // delete [] curpos;
}

