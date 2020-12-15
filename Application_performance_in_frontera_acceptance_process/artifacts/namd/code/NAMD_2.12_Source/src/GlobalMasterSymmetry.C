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
#include "GlobalMasterSymmetry.h"
#include "PDB.h"
#include "PDBData.h"
#include "fitrms.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <stdlib.h>
//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"
//using namespace __gnu_cxx;
using namespace std;
//Read in and parse matrix file
void GlobalMasterSymmetry::parseMatrix(int id, char fileName []){
  int count = 1;
  string line;
  ifstream matrixFile (fileName);
  if (matrixFile.is_open()){
    while (!matrixFile.eof()){
      vector <string> tmp;
      for (int i = 0; i < 4; i++){
      getline(matrixFile, line);
      istringstream iss(line);
      copy(istream_iterator<string>(iss),
            istream_iterator<string>(),
            back_inserter<vector <string>  >(tmp));
      }
      getline(matrixFile, line);
      Matrix4Symmetry tmpmatrix;
      if (tmp.size() < 16){NAMD_die("Error reading matrix file.  Please check layout of the matrix file(s).");}
      for(int j = 0; j < 16; j++){
        tmpmatrix.mat[j] = atof(tmp[j].c_str());
      }
      tmpmatrix.transpose();
      matrices.push_back(tmpmatrix);
      count++;
    }
    matrixFile.close();
  } 
}
GlobalMasterSymmetry::GlobalMasterSymmetry() {
  DebugM(3,"initialize called\n");
  SimParameters *params = Node::Object()->simParameters;
  currentStep = params->firstTimestep;
  firstStep = params->symmetryFirstStep;
  lastStep = params->symmetryLastStep;
  firstFullStep = params->symmetryFirstFullStep;
  lastFullStep = params->symmetryLastFullStep;
  K = params->symmetryk;

  StringList *klist = Node::Object()->configList->find("symmetrykfile");
  if (!K){
    //if (!params->symmetrykfile){NAMD_die("A pdb file containing per-atom force constants must be specified if symmetryk is not in configuration file!");}
    if (!klist){NAMD_die("A pdb file containing per-atom force constants must be specified if symmetryk is not in configuration file!");}
    //symmetrykfile = params->symmetrykfile;
  }
  scaleForces = params->symmetryScaleForces;
  if (scaleForces && lastStep == -1){
    NAMD_die("symmetryLastStep must be specified if symmetryScaleForces is enabled!");
  }
  StringList *matrixList = Node::Object()->configList->find("symmetryMatrixFile");
  StringList *symmetryList = Node::Object()->configList->find("symmetryFile");
  int symfileindex = 0;
  if (!K) {symmetrykfile = klist->data;}
  for (; symmetryList; symmetryList = symmetryList->next) {  
    parseAtoms(symmetryList->data, Node::Object()->molecule->numAtoms, symfileindex);
    if(!K){
      klist = klist->next;
     if (klist){ symmetrykfile = klist->data;}
    }
  }

  map<int, vector<int> >::iterator it = simmap.begin();
  if (!matrixList) {initialTransform();}
  for (; matrixList; matrixList = matrixList->next, ++it){
      parseMatrix(it->first, matrixList->data);
  }

  DebugM(1,"done with initialize\n");
}

//Aligns monomers based on transformation matrices
//found in matrix file
/*
void GlobalMasterSymmetry::alignMonomers(){
  //this is assuming the matrices are written
  //in order of monomer id designation (0, 1, 2,..etc)
  map<int, vector<int> >::iterator simit = simmap.begin();
  for (; simit!=simmap.end(); ++simit){
    for(int x = 0; x < simit->second.size(); x++){
    map<int, vector<int> >::iterator mit = dmap.find(simit->second[x]);
    for (int i = 0; i < mit->second.size(); i++){
      map<int, BigReal *>::iterator it = posmap.find(mit->second[i]);
      matrices[(mit->first)-1].multpoint(it->second);
    }
  }
  }
}
*/
bool GlobalMasterSymmetry::gluInvertMatrix(const BigReal m[16], BigReal invOut[16])
{
  BigReal inv[16], det;
  int i;

  inv[0] =   m[5]*m[10]*m[15] - m[5]*m[11]*m[14] - m[9]*m[6]*m[15]
  + m[9]*m[7]*m[14] + m[13]*m[6]*m[11] - m[13]*m[7]*m[10];
  inv[4] =  -m[4]*m[10]*m[15] + m[4]*m[11]*m[14] + m[8]*m[6]*m[15]
  - m[8]*m[7]*m[14] - m[12]*m[6]*m[11] + m[12]*m[7]*m[10];
  inv[8] =   m[4]*m[9]*m[15] - m[4]*m[11]*m[13] - m[8]*m[5]*m[15]
  + m[8]*m[7]*m[13] + m[12]*m[5]*m[11] - m[12]*m[7]*m[9];
  inv[12] = -m[4]*m[9]*m[14] + m[4]*m[10]*m[13] + m[8]*m[5]*m[14]
  - m[8]*m[6]*m[13] - m[12]*m[5]*m[10] + m[12]*m[6]*m[9];
  inv[1] =  -m[1]*m[10]*m[15] + m[1]*m[11]*m[14] + m[9]*m[2]*m[15]
  - m[9]*m[3]*m[14] - m[13]*m[2]*m[11] + m[13]*m[3]*m[10];
  inv[5] =   m[0]*m[10]*m[15] - m[0]*m[11]*m[14] - m[8]*m[2]*m[15]
  + m[8]*m[3]*m[14] + m[12]*m[2]*m[11] - m[12]*m[3]*m[10];
  inv[9] =  -m[0]*m[9]*m[15] + m[0]*m[11]*m[13] + m[8]*m[1]*m[15]
  - m[8]*m[3]*m[13] - m[12]*m[1]*m[11] + m[12]*m[3]*m[9];
  inv[13] =  m[0]*m[9]*m[14] - m[0]*m[10]*m[13] - m[8]*m[1]*m[14]
  + m[8]*m[2]*m[13] + m[12]*m[1]*m[10] - m[12]*m[2]*m[9];
  inv[2] =   m[1]*m[6]*m[15] - m[1]*m[7]*m[14] - m[5]*m[2]*m[15]
  + m[5]*m[3]*m[14] + m[13]*m[2]*m[7] - m[13]*m[3]*m[6];
  inv[6] =  -m[0]*m[6]*m[15] + m[0]*m[7]*m[14] + m[4]*m[2]*m[15]
  - m[4]*m[3]*m[14] - m[12]*m[2]*m[7] + m[12]*m[3]*m[6];
  inv[10] =  m[0]*m[5]*m[15] - m[0]*m[7]*m[13] - m[4]*m[1]*m[15]
  + m[4]*m[3]*m[13] + m[12]*m[1]*m[7] - m[12]*m[3]*m[5];
  inv[14] = -m[0]*m[5]*m[14] + m[0]*m[6]*m[13] + m[4]*m[1]*m[14]
  - m[4]*m[2]*m[13] - m[12]*m[1]*m[6] + m[12]*m[2]*m[5];
  inv[3] =  -m[1]*m[6]*m[11] + m[1]*m[7]*m[10] + m[5]*m[2]*m[11]
  - m[5]*m[3]*m[10] - m[9]*m[2]*m[7] + m[9]*m[3]*m[6];
  inv[7] =   m[0]*m[6]*m[11] - m[0]*m[7]*m[10] - m[4]*m[2]*m[11]
  + m[4]*m[3]*m[10] + m[8]*m[2]*m[7] - m[8]*m[3]*m[6];
  inv[11] = -m[0]*m[5]*m[11] + m[0]*m[7]*m[9] + m[4]*m[1]*m[11]
  - m[4]*m[3]*m[9] - m[8]*m[1]*m[7] + m[8]*m[3]*m[5];
  inv[15] =  m[0]*m[5]*m[10] - m[0]*m[6]*m[9] - m[4]*m[1]*m[10]
  + m[4]*m[2]*m[9] + m[8]*m[1]*m[6] - m[8]*m[2]*m[5];

  det = m[0]*inv[0] + m[1]*inv[4] + m[2]*inv[8] + m[3]*inv[12];
  if (det == 0)
            return false;

            det = 1.0 / det;

            for (i = 0; i < 16; i++)
                      invOut[i] = inv[i] * det;

                      return true;
}

void GlobalMasterSymmetry::initialTransform(){
  map<int, vector<int> >::iterator simit = simmap.begin();

  for (; simit!=simmap.end(); ++simit){
  map<int, vector<int> >::iterator fit = dmap.find(simit->second[0]);
  BigReal * first = new BigReal [3*fit->second.size()];
  for(int i = 0; i < fit->second.size(); i++){
    map<int, BigReal *>::iterator sit = posmap.find(fit->second[i]);
    first[3*i] = sit->second[0];
    first[3*i+1] = sit->second[1];
    first[3*i+2] = sit->second[2];  
  }

  map<int, vector<int> >::iterator it = dmap.begin();
  for(; it != dmap.end(); ++it){
  if (std::find(simit->second.begin(), simit->second.end(), it->first)!=simit->second.end()){
  BigReal * arr = new BigReal [3*it->second.size()];
    for (int i = 0; i<it->second.size(); i++){
      map<int, BigReal *>::iterator sit = posmap.find(it->second[i]);
      arr[3*i] = sit->second[0];
      arr[3*i+1] = sit->second[1];
      arr[3*i+2] = sit->second[2];  
    }  
      BigReal ttt[16], pre[3], post[3];
      BigReal curRMS = MatrixFitRMS(it->second.size(), arr, first, NULL, ttt);
      int j;
      for (j=0; j<3; j++) {
        post[j] = ttt[4*j+3];
        ttt[4*j+3] = 0;
        pre[j] = ttt[12+j];
        ttt[12+j] = 0;
      }
      Matrix4Symmetry result;
      result.translate(pre);
      result.multmatrix(Matrix4Symmetry(ttt));
      result.translate(post);
      matrices.push_back(result);

      delete [] arr;
  }
  }
  delete [] first;
  }

  for (int i = 0; i < matrices.size(); i++) {
    matrices[i].transpose();
    //iout << "Matrix: " << i << " " << matrices[i].mat[0] << " " << matrices[i].mat[1] << " " << matrices[i].mat[2] << " "<<matrices[i].mat[3] << "\n";
    //iout << "Matrix: " << i << " " << matrices[i].mat[4] << " " << matrices[i].mat[5] << " " << matrices[i].mat[6] << " "<<matrices[i].mat[7] << "\n";  
    //iout << "Matrix: " << i << " " << matrices[i].mat[8] << " " << matrices[i].mat[9] << " " << matrices[i].mat[10] << " "<<matrices[i].mat[11] << "\n";
    //iout << "Matrix: " << i << " " << matrices[i].mat[12] << " " << matrices[i].mat[13] << " " << matrices[i].mat[14] << " "<<matrices[i].mat[15] << "\n";
    //iout <<"\n"<<endi;
    matrices[i].transpose();
  }
}

void GlobalMasterSymmetry::parseAtoms(const char *file, int numTotalAtoms, int symfileindex) {
  DebugM(3,"parseAtoms called\n");
  PDB tmdpdb(file);
  int numatoms = tmdpdb.num_atoms();
  if (numatoms < 1) 
    NAMD_die("No atoms found in symmetryFile\n");
  if (numatoms < numTotalAtoms)
    iout << iWARN << "The number of atoms in symmetryFile is less than the total number of atoms in the structure.\n" << endi;
  if (numatoms > numTotalAtoms)
    NAMD_die("The number of atoms in symmetryFile must not exceed the total number of atoms in the structure!");
 // if ( modifyRequestedAtoms().size() )
 //   NAMD_bug("GlobalMasterSymmetry::parseAtoms() modifyRequestedAtoms() not empty");

  Vector * atompos = new Vector[numatoms];  
  tmdpdb.get_all_positions(atompos);
  if (K){symmetrykfile = file;} 
  PDB kpdb(symmetrykfile);

  int i;
  for (i=0; i<numatoms; i++) {
#ifdef MEM_OPT_VERSION
    PDBCoreData *atom = tmdpdb.atom(i);
#else
    PDBAtom *atom = tmdpdb.atom(i); // get an atom from the file
#endif

    if (atom->occupancy() && atom->temperaturefactor()) { // if occupancy and beta are not 0, then add it!
      // add the atom to the list
      modifyRequestedAtoms().add(i);
      if(!K){
        #ifdef MEM_OPT_VERSION
          PDBCoreData *atomk = kpdb.atom(i);
        #else
          PDBAtom *atomk = kpdb.atom(i); // get an atom from the file
        #endif 
        //kmap[i] = atomk->occupancy();
        kdmap[atom->temperaturefactor()][i] = atomk->occupancy();
      }
      BigReal *arr = new BigReal [3];
      arr[0] = atompos[i].x;
      arr[1] = atompos[i].y;
      arr[2] = atompos[i].z;
      posmap[i] = arr;

      bmap[atom->temperaturefactor()] = atom->occupancy();
        //check to see if monomer id is already in the map
      map <int, vector<int> >::iterator it = dmap.find((int)atom->temperaturefactor());
      map <int, vector<int> >::iterator simit = simmap.find((int)atom->occupancy());
      if (it != dmap.end()){
        it->second.push_back(i); //add atomid to vector in proper existing monomer id
      }
      else{
         dmap[(int)atom->temperaturefactor()] = vector <int> (1,i); //create new monomer id with atomid
      }
      if (simit != simmap.end()){
        if (std::find(simit->second.begin(), simit->second.end(), atom->temperaturefactor()) == simit->second.end()){
        simit->second.push_back(atom->temperaturefactor()); 
        }
      }
      else {
        simmap[(int)atom->occupancy()] = vector <int> (1, atom->temperaturefactor());
      }
    }
  }

   map <int, vector<int> >::iterator simit = simmap.begin();
   for (; simit != simmap.end(); ++simit){
    map <int, vector<int> >::iterator sit = dmap.find(simit->second[0]);
    int numatoms = sit->second.size();
    for (int i = 0; i<simit->second.size(); i++){
      map <int, vector<int> >::iterator fit = dmap.find(simit->second[i]);
      if (fit->second.size() != numatoms){
        NAMD_die("Every monomer must contain the same number of atoms!");
      }
    }
   }  
  delete [] atompos;
}

void GlobalMasterSymmetry::determineAverage() {

   map <int, vector<BigReal *> >::iterator delit = averagePos.begin();
   for (; delit != averagePos.end(); ++delit){
     for (int i = 0; i < delit->second.size(); i++){
       delete [] delit->second[i];
     }
     delit->second.erase(delit->second.begin(), delit->second.end());
   }
   //std::map <int, BigReal * > posmap;
   map <int, BigReal *>::iterator posit;
   map <int, vector<int> >::iterator simit = simmap.begin();
   for (; simit != simmap.end(); ++simit){     
    

    map <int, BigReal *>::iterator pit = posmap.begin();
    for (; pit != posmap.end(); ++pit){delete [] pit->second;}
    posmap.clear();

    map <int, vector<int> >::iterator dit = dmap.begin();
    for (; dit!=dmap.end(); ++dit){
      for (int i = 0; i < dit->second.size(); i++){
        if (std::find(simit->second.begin(), simit->second.end(), dit->first)!=simit->second.end()){

          BigReal* arr = new BigReal[3];
          arr[0] = positions[dit->second[i]].x;
          arr[1] = positions[dit->second[i]].y;
          arr[2] = positions[dit->second[i]].z;

          posmap[dit->second[i]] = arr;
        }
      } 
    }
    averagePos[simit->first] = vector <BigReal *> (); 
    map <int, vector<int> >::iterator it = dmap.begin();
    map <int, vector<int> >::iterator sit = dmap.find(simit->second[0]);
    int numatoms = sit->second.size();
    for (int i = 0; i < numatoms; i++){
       BigReal *arr = new BigReal [3];
      arr[0] = 0;
      arr[1] = 0;
      arr[2] = 0;
      for (; it!=dmap.end(); ++it){
        if (std::find(simit->second.begin(), simit->second.end(), it->first)!=simit->second.end()){
          posit = posmap.find(it->second[i]);
          matrices[(it->first)-1].multpoint(posit->second);
          arr[0] += posit->second[0];
          arr[1] += posit->second[1];
          arr[2] += posit->second[2];
        }
      }
      it = dmap.begin();
      BigReal *avg = new BigReal[3];
      avg[0] = arr[0]/(simit->second.size());
      avg[1] = arr[1]/(simit->second.size());
      avg[2] = arr[2]/(simit->second.size());
      averagePos[simit->first].push_back(avg);
      delete [] arr;
    }
    
   }

}

void GlobalMasterSymmetry::backTransform(){
  map <int, BigReal *>::iterator bit = backavg.begin();
  for (; bit != backavg.end(); ++bit){delete [] bit->second;}
  backavg.clear();

  map <int, vector<int> >::iterator it = dmap.begin();
  for (; it!=dmap.end(); ++it){
    map<int, int >::iterator bmit = bmap.find(it->first);
    int bm = bmit->second;
    map<int, vector <BigReal *> >::iterator avit = averagePos.find(bmit->second);
    int numatoms = it->second.size();
    BigReal *avg = new BigReal [3*numatoms];
    for (int i = 0; i < numatoms; i++){
      avg[3*i] = avit->second[i][0];
      avg[3*i+1] = avit->second[i][1];
      avg[3*i+2] = avit->second[i][2];
    }
    BigReal inverse[16];
    matrices[it->first-1].transpose();
    gluInvertMatrix(matrices[it->first-1].mat, inverse);
    Matrix4Symmetry inv(inverse);
    inv.transpose();
    matrices[it->first-1].transpose();
 
    for (int k = 0; k < numatoms; k++){
      inv.multpoint(avg+3*k);
    }
    backavg[it->first] = avg;
   }

}
GlobalMasterSymmetry::~GlobalMasterSymmetry() { 
  map <int, BigReal *>::iterator pit = posmap.begin();
  for (; pit != posmap.end(); ++pit){delete pit->second;}
}
void GlobalMasterSymmetry::calculate() {
  // have to reset my forces every time.  
  modifyAppliedForces().resize(0);
  modifyForcedAtoms().resize(0);

  // see if symmetry restraints should be active
  if (currentStep < firstStep || (currentStep >= lastStep && lastStep != -1)) {
    currentStep++;
    return;
  }

  //map<int, Position> positions;
  AtomIDList::const_iterator a_i = getAtomIdBegin();
  AtomIDList::const_iterator a_e = getAtomIdEnd();
  PositionList::const_iterator p_i = getAtomPositionBegin();

  //create mapping of positions now to avoid
  //going through these iterators for each domain
  for ( ; a_i != a_e; ++a_i, ++p_i ){
    positions[*a_i] = *p_i;
  }

  map <int, vector<int> >::iterator it;
  //map <int, BigReal *>::iterator pit = posmap.begin();
  //for (; pit != posmap.end(); ++pit){delete [] pit->second;}

 // posmap.clear();
  //for (it = dmap.begin(); it != dmap.end(); ++it){
    // fetch the current coordinates
    
    //for (int i = 0; i < it->second.size(); i++){

      //BigReal* arr = new BigReal[3];
      //arr[0] = positions[it->second[i]].x;
      //arr[1] = positions[it->second[i]].y;
      //arr[2] = positions[it->second[i]].z;

     // posmap[it->second[i]] = arr;
   // } 
//}

//  alignMonomers();
  determineAverage();
  backTransform();

  //iterate through all domains
  for (it = dmap.begin(); it != dmap.end(); ++it){

  // fetch the current coordinates
  BigReal *curpos = new BigReal[3*(it->second.size())];
  for (int i = 0; i < it->second.size(); i++){
    curpos[3*i  ] = positions[it->second[i]].x;
    curpos[3*i+1] = positions[it->second[i]].y;
    curpos[3*i+2] = positions[it->second[i]].z;  
  }


  BigReal *tmpavg = backavg[it->first];

  for (int i=0; i<it->second.size(); i++) {
    BigReal k; 
    if(!K){
     //k = kmap[it->second[i]];  
     k = kdmap[it->first][it->second[i]];
    }
    else{
     k = K/it->second.size();  
    }
    BigReal maxforce2 = 0.;
    BigReal frac = -k;
    if (scaleForces){

    if (currentStep < firstFullStep){
      BigReal linear_evolve = (BigReal(currentStep-firstStep)) /
                            (firstFullStep-firstStep);
      frac = frac*(linear_evolve);
    }
    if (currentStep > lastFullStep)
    {
      BigReal linear_evolve = (BigReal(currentStep-lastFullStep)) /
                            (lastStep-lastFullStep);
      frac = frac*(1-linear_evolve);
    }

    }
    BigReal dx = curpos[3*i  ] - tmpavg[3*i];
//    iout << "DX: " << dx << "\n" << endi;
    BigReal dy = curpos[3*i+1] - tmpavg[3*i+1];
    BigReal dz = curpos[3*i+2] - tmpavg[3*i+2];
    BigReal fvec[3] = { dx, dy, dz };
    Vector force(fvec[0]*frac, fvec[1]*frac, fvec[2]*frac);
    modifyForcedAtoms().add(it->second[i]);
    modifyAppliedForces().add(force);
    BigReal force2 = force.length2();
    if ( force2 > maxforce2 ) maxforce2 = force2;
  } 
  delete [] curpos;
 }
   
  currentStep++;
}
