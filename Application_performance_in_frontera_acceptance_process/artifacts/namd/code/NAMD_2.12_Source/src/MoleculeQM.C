/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifdef DEBUG_QM
  #define DEBUGM
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "ResizeArray.h"
#include "InfoStream.h"
#include "Molecule.h"
#include "strlib.h"
#include "MStream.h"
#include "Communicate.h"
#include "Node.h"
#include "ObjectArena.h"
#include "Parameters.h"
#include "PDB.h"
#include "SimParameters.h"
#include "Hydrogen.h"
#include "UniqueSetIter.h"
#include "charm++.h"

#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"
#include "CompressPsf.h"
#include "ParallelIOMgr.h"
#include <deque>
#include <algorithm>

#ifndef CODE_REDUNDANT
#define CODE_REDUNDANT 0
#endif

#include <string>
#include <sstream>

class qmSolvData {
public:
    char segName[5];
    int resID, begAtmID, size;
    std::vector<int> atmIDs ;
    Real qmGrpID ;
    
    qmSolvData() : resID(-1), begAtmID(-1), size(0) {}
    qmSolvData(int newResID, int newBegAtm, int newSize) {
        resID = newResID;
        begAtmID = newBegAtm;
        size = newSize;
    }
    qmSolvData(int newResID, int newBegAtm, int newSize, 
               const char* newSegName, Real newQMID) {
        resID = newResID;
        begAtmID = newBegAtm;
        size = newSize;
        strncpy(segName, newSegName,5);
        atmIDs.push_back(newBegAtm) ;
        qmGrpID = newQMID;
    }
    
    bool operator <(const qmSolvData& ref) {return begAtmID < ref.begAtmID;}
    bool operator ==(const qmSolvData& ref) {return begAtmID == ref.begAtmID;}
};

std::vector<std::string> split(const std::string &text, std::string delimiter) {
    
    std::vector<std::string> tokens;
    std::size_t start = 0, end = 0;
    
    while ((end = text.find(delimiter, start)) != std::string::npos) {
        
        std::string temp = text.substr(start, end - start);
        
        if (! temp.empty()) tokens.push_back(temp);
        
        start = end + delimiter.length();
    }
    
    // Gets last item
    std::string temp = text.substr(start);
    if (! temp.empty()) tokens.push_back(temp);
    
    return tokens;
}

struct refSelStr {
    std::string segid;
    int resid;
    
    refSelStr() {} ;
    refSelStr(std::string newSegID, int newResid) {
        segid = newSegID;
        resid = newResid;
    }
} ;

typedef std::vector<refSelStr> refSelStrVec ;
typedef std::map<Real, refSelStrVec> refSelStrMap ;
typedef std::pair<Real,refSelStrVec> refSelStrPair ;

void Molecule::prepare_qm(const char *pdbFileName, 
                               Parameters *params, ConfigList *cfgList) {
#ifdef MEM_OPT_VERSION
    NAMD_die("QMMM interface is not supported in memory optimized builds.");
#else
    
    PDB *pdbP;      //  Pointer to PDB object to use
    
    qmNumQMAtoms = 0;
    
    qmNumGrps = 0 ;
    
    iout << iINFO << "Using the following PDB file for QM parameters: " << 
    pdbFileName << "\n" << endi;
    if (pdbFileName)
        pdbP = new PDB(pdbFileName);
    else
        iout << "pdbFile name not defined!\n\n" << endi ;
    
    if (pdbP->num_atoms() != numAtoms)
    {
       NAMD_die("Number of atoms in QM parameter PDB file doesn't match coordinate PDB");
    }
    
    qmElementArray = new char*[numAtoms] ;
    
    qmAtomGroup = new Real[numAtoms] ;
    
    BigReal bondVal;
    Real atmGrp;
    
    std::set<Real> atmGrpSet;
    
    std::vector<Real> grpChrgVec;
    
    // Value in the bond column fro QM-MM bonds.
    std::vector<BigReal> qmBondValVec;
    // Identifiers of different QM regions
    std::vector<Real> qmGroupIDsVec;
    // This maps the qm group ID with the serail index of this group in 
    // various arrays.
    std::map<Real,int> qmGrpIDMap ;
    
    std::vector<int> qmAtmIndxVec, qmGrpSizeVec;
    
    // Used to parse user selection in different places.
    std::vector<std::string> strVec;
    
    qmNoPC = simParams->qmNoPC;
    
    // We set to zero by default, So there is no need for extra processing.
    qmPCFreq = 0;
    if (simParams->qmPCSelFreq > 1)
        qmPCFreq = simParams->qmPCSelFreq;
    
    
    ///////////////////////////////
    /// Prepares Live Solvent Selection
    
    
    qmLSSTotalNumAtms = 0;
    SortedArray< qmSolvData> lssGrpRes;
    std::map<Real,std::vector<int> > lssGrpRefIDs;
    refSelStrMap lssRefUsrSel;
    int totalNumRefAtms = 0;
    int lssClassicResIndx = 0 ;
    int lssCurrClassResID = -1 ;
    char lssCurrClassResSegID[5];
    if (simParams->qmLSSOn) {
        DebugM(4, "LSS is ON! Processing QM solvent.\n") ;
        
        if (resLookup == NULL)
            NAMD_die("We need residue data to conduct LSS.");
         
        if (simParams->qmLSSMode == QMLSSMODECOM) {
            
            StringList *current = cfgList->find("QMLSSRef");
            for ( ; current; current = current->next ) {
                
                strVec = split( std::string(current->data) , " ");
                
                if (strVec.size() != 3 ) {
                    iout << iERROR << "Format error in QM LSS size: " 
                    << current->data
                    << "\n" << endi;
                    NAMD_die("Error processing QM information.");
                }
                
                std::stringstream storConv ;
                
                storConv << strVec[0] ;
                Real grpID ;
                storConv >> grpID;
                if (storConv.fail()) {
                    iout << iERROR << "Error parsing QMLSSRef selection: " 
                    << current->data
                    << "\n" << endi;
                    NAMD_die("Error processing QM information.");
                }
                
                std::string segName = strVec[1].substr(0,4);
                
                storConv.clear() ;
                storConv << strVec[2];
                int resID ;
                storConv >> resID;
                if (storConv.fail()) {
                    iout << iERROR << "Error parsing QMLSSRef selection: " 
                    << current->data
                    << "\n" << endi;
                    NAMD_die("Error processing QM information.");
                }
                
                auto it = lssRefUsrSel.find(grpID) ;
                if (it == lssRefUsrSel.end())
                    lssRefUsrSel.insert(refSelStrPair(grpID,refSelStrVec()));
                
                lssRefUsrSel[grpID].push_back(refSelStr(segName, resID));
            }
            
            for (auto it=lssRefUsrSel.begin(); it!=lssRefUsrSel.end(); it++) {
                iout << iINFO << "LSS user selection for COM composition of group "
                << it->first << ":\n" << endi ;
                for (int i=0; i<it->second.size();i++) {
                    iout << iINFO << "Segment " << it->second[i].segid 
                    << " ; residue " << it->second[i].resid << "\n" << endi ;
                }
            }
        }
    }
    
    
    
    ///////////////////////////////
    /// Data gathering from PDB to find QM atom and bond info
    
    
    for (int atmInd = 0 ; atmInd < numAtoms; atmInd++) {
        
        // If we are looking for QM-MM bonds, then we need to store extra info.
        if (simParams->qmBondOn) {
            
            // We store both the qm group and the bond value
            if ( strcmp(simParams->qmColumn,"beta") == 0 ){
                atmGrp = pdbP->atom(atmInd)->temperaturefactor() ;
                
                bondVal = pdbP->atom(atmInd)->occupancy() ;
            }
            else {
                atmGrp = pdbP->atom(atmInd)->occupancy() ;
                
                bondVal = pdbP->atom(atmInd)->temperaturefactor() ;
            }
            
            qmBondValVec.push_back(bondVal);
        }
        else {
            
            if ( strcmp(simParams->qmColumn,"beta") == 0 ){
                atmGrp = pdbP->atom(atmInd)->temperaturefactor() ;
            }
            else {
                atmGrp = pdbP->atom(atmInd)->occupancy() ;
            }
        }
        
        // We store all atom QM Group IDs in an array for 
        // future transport to all PEs.
        qmAtomGroup[atmInd] = atmGrp;
        
        // We store all atom's elements for quick access in the QM code.
        // Elements are used to tell the QM code the atom type.
        qmElementArray[atmInd] = new char[3];
        strncpy(qmElementArray[atmInd],pdbP->atom(atmInd)->element(),3);
        
        // For QM atoms, we keep global atom indices and parse different QM Group
        // IDs, keeping track of each groups size
        if (atmGrp > 0){
            
            if (simParams->fixedAtomsOn) {
                if ( fixedAtomFlags[atmInd] == 1 ) {
                    iout << iERROR << "QM atom cannot be fixed in space!\n" 
                    << endi;
                    NAMD_die("Error processing QM information.");
                }
            }
            
            // Changes the VdW type of QM atoms.
            // This is may be used to counter balance the overpolarization 
            // that QM atoms sufer.
            if (simParams->qmVDW) {
                // Creating a new type string
                std::string qmType(qmElementArray[atmInd]) ;
                // Erases empty spaces
                qmType.erase(
                        std::remove_if(qmType.begin(), qmType.end(), isspace ),
                    qmType.end());
                // pre-pends a "q" to the element name.
                qmType = std::string("q") + qmType;
            
//                 iout << "QM VdW type: " << atoms[atmInd].vdw_type 
//                 << " atom type: " << atomNames[atmInd].atomtype 
//                 << " new type "<< qmType << "\n" << endi;
                
                /*  Look up the vdw constants for this atom    */
                // I am casting a non-const here because the function doesn't actually
                // change the string value, but it doesn't ask for a const char* as 
                // an argument.
                params->assign_vdw_index(const_cast<char*>( qmType.c_str()), 
                                         &(atoms[atmInd]));
                
//                 iout << "----> new VdW type: " << atoms[atmInd].vdw_type << "\n" << endi;
            }
            
            // Indexes all global indices of QM atoms.
            qmAtmIndxVec.push_back(atmInd);
            
            auto retVal = atmGrpSet.insert(atmGrp);
            
            // If a new QM group is found, store the reverse reference in a map
            // and increase the total count.
            if (retVal.second) {
                // This mak makes the reverse identification from the group ID
                // to the sequential order in which each group was first found.
                qmGrpIDMap.insert(std::pair<BigReal,int>(atmGrp, qmNumGrps));
                
                // This vector keeps track of the group ID for each group
                qmGroupIDsVec.push_back(atmGrp);
                
                // This counter is used to keep track of the sequential order in
                // which QM groups were first seen in the reference PDB file.
                qmNumGrps++ ;
                
                // If this is a new QM group, initialize a new variable in the 
                // vector to keep track of the number of atoms in this group.
                qmGrpSizeVec.push_back(1);
                
                // For each new QM group, a new entry in the total charge
                // vector is created
                grpChrgVec.push_back( atoms[atmInd].charge );
            }
            else {
                // If the QM group has already been seen, increase the count
                // of atoms in that group.
                qmGrpSizeVec[ qmGrpIDMap[atmGrp] ]++ ;
                
                // If the QM group exists, adds current atom charge to total charge.
                grpChrgVec[ qmGrpIDMap[atmGrp] ] += atoms[atmInd].charge;
            }
            
            // In case we are using live solvent selection
            if(simParams->qmLSSOn) {
                
                int resID = pdbP->atom(atmInd)->residueseq();
                char segName[5];
                strncpy(segName, pdbP->atom(atmInd)->segmentname(),5);
                
                int resSize = get_residue_size(segName,resID);
                
                int i =-1, end =-1;
                
                resLookup->lookup(segName,resID,&i, &end);
                
                if ( strncasecmp(pdbP->atom(atmInd)->residuename(),
                           simParams->qmLSSResname, 4) == 0) {
                    
                    // We try to insert every residue from every atom
                    qmSolvData *resP = lssGrpRes.find(qmSolvData(resID, i, resSize));
                    
                    if (resP != NULL) {
                        resP->atmIDs.push_back(atmInd) ;
//                         DebugM(3, "Existing residue " << resP->resID 
//                         << " from segID " << resP->segName
//                         << " received atom "
//                         << atmInd << "\n" );
                    }
                    else {
                        lssGrpRes.insert(qmSolvData(resID, i, resSize, segName, atmGrp));
//                         DebugM(3, lssGrpRes.size() << ") Inserted new residue: " 
//                         << resID << " from segID " << segName
//                         << " with atom "
//                         << i << "\n" ) ;
                    }
                }
                else {
                    // We store the QM atoms, per group, which are NOT part of
                    // solvent molecules.
                    
                    // Checks if we have a vector for this QM group.
                    auto itNewGrp = lssGrpRefIDs.find(atmGrp) ;
                    if (itNewGrp == lssGrpRefIDs.end()) {
                        lssGrpRefIDs.insert(std::pair<Real, std::vector<int> >(
                            atmGrp, std::vector<int>() ) );
                    }
                    
                    switch (simParams->qmLSSMode)
                    {
                    
                    case QMLSSMODECOM:
                    {
                        // If we are using COM selection, checks if the atom
                        // is part of the user selected residues
                        for(int i=0; i<lssRefUsrSel[atmGrp].size(); i++) {
                            if (lssRefUsrSel[atmGrp][i].resid == resID &&
                                strncmp(lssRefUsrSel[atmGrp][i].segid.c_str(),
                                        segName,5) == 0 ) {
                                
                                lssGrpRefIDs[atmGrp].push_back(atmInd);
                                totalNumRefAtms++;
                            }
                        }
                        
                    } break;
                    
                    case QMLSSMODEDIST:
                    {
                    
                        lssGrpRefIDs[atmGrp].push_back(atmInd);
                        totalNumRefAtms++;
                    
                    } break;
                    
                    }
                }
                
            }
            
        }
        else if (atmGrp == 0) {
            
            if(simParams->qmLSSOn) {
                
                if ( strncasecmp(pdbP->atom(atmInd)->residuename(),
                           simParams->qmLSSResname, 4) == 0) {
                    
                    int resID = pdbP->atom(atmInd)->residueseq();
                    char segName[5];
                    strncpy(segName, pdbP->atom(atmInd)->segmentname(),5);
                    
                    if (lssCurrClassResID < 0) {
                        lssCurrClassResID = resID ;
                        strncpy(lssCurrClassResSegID, segName,5);
                        lssClassicResIndx = 0;
                    }
                    else if (lssCurrClassResID != resID ||
                            strcmp(lssCurrClassResSegID, segName) != 0 ) {
                        lssCurrClassResID = resID ;
                        strncpy(lssCurrClassResSegID, segName,5);
                        lssClassicResIndx++;
                    }
                    
                    qmClassicSolv.insert(std::pair<int,int>(atmInd,lssClassicResIndx));
                    
                }
            }
        }
        else if(atmGrp < 0) {
            iout << iERROR << "QM group ID cannot be less than zero!\n" << endi;
            NAMD_die("Error processing QM information.");
        }
    }
    
    // Sanity check
    if (simParams->qmLSSOn) {
        if (lssGrpRes.size() == 0)
            NAMD_die("LSS was activated but there are no QM solvent molecules!\n");
        
        for (auto it=lssRefUsrSel.begin(); it!=lssRefUsrSel.end(); it++) {
            auto itTarget = qmGrpIDMap.find(it->first);
            if (itTarget == qmGrpIDMap.end()) {
                iout << iERROR << "QM group ID for LSS could not be found in input!"
                << " QM group ID: " << it->first << "\n" << endi;
                NAMD_die("Error processing QM information.");
            }
        }
        
        DebugM(3,"We found " << lssClassicResIndx << " classical solvent molecules totalling "
        << qmClassicSolv.size() << " atoms.\n" );
        
    }
    
    qmNumQMAtoms = qmAtmIndxVec.size();
    
    if (qmNumQMAtoms == 0)
        NAMD_die("No QM atoms were found in this QM simulation!") ;
    
    iout << iINFO << "Number of QM atoms (excluding Dummy atoms): " << 
    qmNumQMAtoms << "\n" << endi;
    
    qmAtmChrg = new Real[qmNumQMAtoms] ;
    qmAtmIndx = new int[qmNumQMAtoms] ;
    for (int i=0; i<qmNumQMAtoms; i++) {
        // qmAtmChrg gets initialized with the PSF charges at the end of this
        // function, but values may change as QM steps are taken.
        qmAtmChrg[i] = 0;  
        qmAtmIndx[i] = qmAtmIndxVec[i] ;
    }
    
    // This map relates the QM group index with a vector of pairs
    // of bonded MM-QM atoms (with the bonded atoms ins this order: 
    // MM first, QM second).
    std::map<int, std::vector<std::pair<int,int> > > qmGrpIDBonds ;
    int bondCounter = 0;
    if (simParams->qmBondOn) {
        
        // Checks all bonds for QM-MM pairs.
        // Makes sanity checks against QM-QM pairs and MM-MM pairs that
        // are flagged by the user to be bonds between QM and MM regions.
        // QM-QM bonds will be removed in another function.
        for (int bndIt = 0; bndIt < numBonds; bndIt++) {
            
            bond curr = bonds[bndIt] ;
            
            // In case either atom has a non-zero
            if ( qmBondValVec[curr.atom1] != 0 &&
                 qmBondValVec[curr.atom2] != 0 ) {
                
                // Sanity checks (1 of 2)
                if (qmAtomGroup[curr.atom1] != 0 &&
                    qmAtomGroup[curr.atom2] != 0) {
                    iout << iERROR << "Atoms " << curr.atom1 << " and " << 
                    curr.atom2 << " are assigned as QM atoms.\n" << endi;
                    NAMD_die("Error in QM-MM bond assignment.") ;
                }
                
                // Sanity checks (2 of 2)
                if (qmAtomGroup[curr.atom1] == 0 &&
                    qmAtomGroup[curr.atom2] == 0) {
                    iout << iERROR << "Atoms " << curr.atom1 << " and " << 
                    curr.atom2 << " are assigned as MM atoms.\n" << endi;
                    NAMD_die("Error in QM-MM bond assignment.") ;
                }
                
                int currGrpID = 0;
                std::pair<int,int> newPair(0,0);
                
                // We create a QM-MM pair with the MM atom first
                if (qmAtomGroup[curr.atom1] != 0) {
                    newPair = std::pair<int,int>(curr.atom2, curr.atom1) ;
                    currGrpID = qmAtomGroup[curr.atom1];
                } else {
                    newPair = std::pair<int,int>(curr.atom1, curr.atom2) ;
                    currGrpID = qmAtomGroup[curr.atom2];
                } 
                
                int grpIndx = qmGrpIDMap[currGrpID] ;
                
                // Stores bonds in vectors key-ed by the QM group ID of the QM atom.
                auto retIter = qmGrpIDBonds.find(grpIndx) ;
                // In case thi QM-MM bonds belong to a QM group we have not seen 
                // yet, stores a new vector in the map.
                if (retIter == qmGrpIDBonds.end()) {
                    qmGrpIDBonds.insert(std::pair<BigReal,std::vector<std::pair<int,int> > >(
                    grpIndx, std::vector<std::pair<int,int> >() ) );
                }
                
                qmGrpIDBonds[grpIndx].push_back( newPair );
                
                bondCounter++ ;
            }
            
            
        }
        
//         iout << "Finished processing QM-MM bonds.\n" << endi ;
        
        if (bondCounter == 0)
            iout << iWARN << "We found " << bondCounter << " QM-MM bonds.\n" << endi ;
        else
            iout << iINFO << "We found " << bondCounter << " QM-MM bonds.\n" << endi ;
        
    }
    
    // Initializes several arrays used to setup the QM simulation.
    qmNumBonds = bondCounter ;
    
    qmMMBond = new int*[qmNumBonds];
    
    qmDummyBondVal = new BigReal[qmNumBonds];
    
    qmMMChargeTarget = new int*[qmNumBonds] ;
    qmMMNumTargs = new int[qmNumBonds] ;
    
    qmDummyElement = new char*[qmNumBonds] ;
    
    
    qmNumGrps = atmGrpSet.size();
    
    qmGrpSizes = new int[qmNumGrps];
    
    qmGrpID = new Real[qmNumGrps];
    
    qmGrpChrg = new Real[qmNumGrps];
    
    qmGrpMult = new Real[qmNumGrps] ;
    
    qmGrpNumBonds = new int[qmNumGrps];
    
    qmGrpBonds = new int*[qmNumGrps];
    qmMMBondedIndx = new int*[qmNumGrps];
    
    
    ///////////////////////////////
    /// Multiplicity of each QM group
    
    
    // We first initialize the multiplicity vector with 
    // default values, then populate it with user defined values.
    for (int grpIter = 0; grpIter < qmNumGrps; grpIter++) {
        qmGrpMult[grpIter] = 0;
    }
    
    int multCount = 0;
    StringList *current = cfgList->find("QMMult");
    for ( ; current; current = current->next ) {
        
        auto strVec = split( std::string(current->data) , " ");
        
        if (strVec.size() != 2 ) {
            iout << iERROR << "Format error in QM Multiplicity string: " 
            << current->data
            << "\n" << endi;
            NAMD_die("Error processing QM information.");
        }
        
        std::stringstream storConv ;
        
        storConv << strVec[0] ;
        Real grpID ;
        storConv >> grpID;
        if (storConv.fail()) {
            iout << iERROR << "Error parsing QMMult selection: " 
            << current->data
            << "\n" << endi;
            NAMD_die("Error processing QM information.");
        }
        
        storConv.clear() ;
        storConv << strVec[1];
        Real multiplicity ;
        storConv >> multiplicity;
        if (storConv.fail()) {
            iout << iERROR << "Error parsing QMMult selection: " 
            << current->data
            << "\n" << endi;
            NAMD_die("Error processing QM information.");
        }
        
        auto it = qmGrpIDMap.find(grpID);
        
        if (it == qmGrpIDMap.end()) {
            iout << iERROR << "Error parsing QMMult selection. Could not find QM group ID: " 
            << grpID
            << "\n" << endi;
            NAMD_die("Error processing QM information.");
        }
        else {
            iout << iINFO << "Applying user defined multiplicity "
            << multiplicity << " to QM group ID " << grpID << "\n" << endi;
            qmGrpMult[it->second] = multiplicity;
        }
        
        multCount++;
    }
    
    if (multCount != qmNumGrps && simParams->qmFormat == QMFormatORCA) {
        NAMD_die("ORCA neds multiplicity values for all QM regions.");
    }
    
    
    ///////////////////////////////
    /// Charge of each QM group
    
    
    for (size_t grpIndx = 0; grpIndx < qmGrpSizeVec.size(); grpIndx++) {
        
        bool nonInteger = true;
        if ((fabsf(roundf(grpChrgVec[grpIndx]) - grpChrgVec[grpIndx]) <= 0.001f) ) {
            grpChrgVec[grpIndx] = roundf(grpChrgVec[grpIndx]) ;
            nonInteger = false;
        }
        
        iout << iINFO << grpIndx + 1 << ") Group ID: " << qmGroupIDsVec[grpIndx]
        << " ; Group size: " << qmGrpSizeVec[grpIndx] << " atoms"
        << " ; Total charge: " << grpChrgVec[grpIndx] << "\n" << endi ;
        
        if (nonInteger && simParams->PMEOn)
            NAMD_die("QM atoms do not add up to a whole charge, which is needed for PME.") ;
    }
    
    current = cfgList->find("QMCharge");
    for ( ; current; current = current->next ) {
        
        auto strVec = split( std::string(current->data) , " ");
        
        if (strVec.size() != 2 ) {
            iout << iERROR << "Format error in QM Charge string: " 
            << current->data
            << "\n" << endi;
            NAMD_die("Error processing QM information.");
        }
        
        std::stringstream storConv ;
        
        storConv << strVec[0] ;
        Real grpID ;
        storConv >> grpID;
        if (storConv.fail()) {
            iout << iERROR << "Error parsing QMCharge selection: " 
            << current->data
            << "\n" << endi;
            NAMD_die("Error processing QM information.");
        }
        
        storConv.clear() ;
        storConv << strVec[1];
        Real charge ;
        storConv >> charge;
        if (storConv.fail()) {
            iout << iERROR << "Error parsing QMCharge selection: " 
            << current->data
            << "\n" << endi;
            NAMD_die("Error processing QM information.");
        }
        
        auto it = qmGrpIDMap.find(grpID);
        
        if (it == qmGrpIDMap.end()) {
            iout << iERROR << "Error parsing QMCharge selection. Could not find QM group ID: " 
            << grpID
            << "\n" << endi;
            NAMD_die("Error processing QM information.");
        }
        else {
            iout << iINFO << "Applying user defined charge "
            << charge << " to QM group ID " << grpID << "\n" << endi;
            grpChrgVec[it->second] = charge;
        }
    }
    
    
    // If mechanichal embedding was requested but we have QM-MM bonds, we need 
    // to send extra info to ComputeQM to preserve calculation speed.
    if (qmNumBonds > 0 && qmNoPC) {
        qmMeNumBonds = qmNumBonds;
        qmMeMMindx = new int[qmMeNumBonds] ;
        qmMeQMGrp = new Real[qmMeNumBonds] ;
    }
    else {
        qmMeNumBonds = 0 ;
        qmMeMMindx = NULL;
        qmMeQMGrp = NULL;
    }
    
    
    ///////////////////////////////
    /// Populate arrays that are used throughout the the calculations.
    
    
    bondCounter = 0;
    for (int grpIter = 0; grpIter < qmNumGrps; grpIter++) {
        
        qmGrpID[grpIter] = qmGroupIDsVec[grpIter] ;
        
        qmGrpSizes[grpIter] = qmGrpSizeVec[grpIter] ;
        
        qmGrpChrg[grpIter] = grpChrgVec[grpIter];
        
//         iout << "Loaded " << qmGrpSizes[grpIter] << " atoms to this group.\n" << endi;
        
        int currNumbBonds = qmGrpIDBonds[grpIter].size() ;
        
        // Assigns the number of bonds that the current QM group has.
        qmGrpNumBonds[grpIter] = currNumbBonds;
        
        if (currNumbBonds > 0) {
            
            qmGrpBonds[grpIter] = new int[currNumbBonds];
            qmMMBondedIndx[grpIter] = new int[currNumbBonds];
            
            for (int bndIter = 0; bndIter<currNumbBonds; bndIter++) {
                
                // Adds the bonds to the overall sequential list.
                qmMMBond[bondCounter] = new int[2] ;
                qmMMBond[bondCounter][0] = qmGrpIDBonds[grpIter][bndIter].first ;
                qmMMBond[bondCounter][1] = qmGrpIDBonds[grpIter][bndIter].second ;
                
                // For the current QM group, and the current bond, gets the bond index.
                qmGrpBonds[grpIter][bndIter] =  bondCounter;
                
                // For the current QM group, and the current bond, gets the MM atom.
                qmMMBondedIndx[grpIter][bndIter] = qmGrpIDBonds[grpIter][bndIter].first ;
                
                // Assign the default value of dummy element
                qmDummyElement[bondCounter] = new char[3];
                strcpy(qmDummyElement[bondCounter],"H\0");
                
                // Determines the distance that will separate the new Dummy atom
                // and the Qm atom to which it will be bound.
                bondVal = 0;
                if (simParams->qmBondDist) {
                    if (qmBondValVec[qmMMBond[bondCounter][0]] != 
                        qmBondValVec[qmMMBond[bondCounter][1]]
                    ) {
                        iout << iERROR << "qmBondDist is ON but the values in the bond column are different for atom "
                        << qmMMBond[bondCounter][0] << " and " << qmMMBond[bondCounter][1]
                        << ".\n" << endi ;
                        NAMD_die("Error in QM-MM bond processing.");
                    }
                    
                    bondVal = qmBondValVec[qmMMBond[bondCounter][0]] ;
                } 
                else {
                    
                    if (strcmp(qmElementArray[qmMMBond[bondCounter][1]],"C") == 0 ) {
                        bondVal = 1.09 ;
                    }
                    else if (strcmp(qmElementArray[qmMMBond[bondCounter][1]],"O") == 0 ) {
                        bondVal = 0.98 ;
                    }
                    else if (strcmp(qmElementArray[qmMMBond[bondCounter][1]],"N") == 0 ) {
                        bondVal = 0.99 ;
                    }
                    else {
                        iout << iERROR << "qmBondDist is OFF but the QM atom has no default distance value. Atom "
                        << qmMMBond[bondCounter][1] << ", with element: " 
                        << qmElementArray[qmMMBond[bondCounter][1]] 
                        <<  ".\n" << endi ;
                        NAMD_die("Error in QM-MM bond processing.");
                    }
                    
                }
                
                iout << iINFO << "MM-QM pair: " << qmMMBond[bondCounter][0] << ":"
                << qmMMBond[bondCounter][1] 
                << " -> Value (distance or ratio): " << bondVal
                << " (QM Group " << grpIter << " ID " << qmGrpID[grpIter] << ")"
                << "\n" << endi ;
                
                qmDummyBondVal[bondCounter] = bondVal;
                
                // In case we are preparing for a mechanical embedding simulation
                // with no point charges, populate the following vectors
                if (qmMeNumBonds > 0) {
                    qmMeMMindx[bondCounter] = qmMMBond[bondCounter][0];
                    qmMeQMGrp[bondCounter] = qmGrpID[grpIter];
                }
                
                bondCounter++ ;
            }
        }
        
    }
    
    
    ///////////////////////////////
    /// Overides Link Atom element with user selection.
    
    current = NULL;
    if (qmNumBonds > 0)
        current = cfgList->find("QMLinkElement");
    
    int numParsedLinkElem = 0;
    for ( ; current != NULL; current = current->next ) {
        
        DebugM(3,"Parsing link atom element: " << current->data << "\n" );
        
        strVec = split( std::string(current->data) , " ");
        
        // We need the two atoms that compose the QM-MM bonds and 
        // then the element.
        if (strVec.size() != 3 && qmNumBonds > 1) {
            iout << iERROR << "Format error in QM link atom element selection: " 
            << current->data
            << "\n" << endi;
            NAMD_die("Error processing QM information.");
        }
        
        // If there is only one QM-MM bond, we can accept the element only.
        if (strVec.size() != 1 && strVec.size() != 3 && qmNumBonds == 1) {
            iout << iERROR << "Format error in QM link atom element selection: " 
            << current->data
            << "\n" << endi;
            NAMD_die("Error processing QM information.");
        }
        
        std::stringstream storConv ;
        int bondAtom1, bondAtom2;
        std::string linkElement ;
        
        if (strVec.size() == 1) {
            linkElement = strVec[0].substr(0,2);
        }
        else if (strVec.size() == 3) {
            
            storConv << strVec[0] ;
            storConv >> bondAtom1;
            if (storConv.fail()) {
                iout << iERROR << "Error parsing link atom element selection: " 
                << current->data
                << "\n" << endi;
                NAMD_die("Error processing QM information.");
            }
            
            storConv.clear() ;
            storConv << strVec[1];
            storConv >> bondAtom2;
            if (storConv.fail()) {
                iout << iERROR << "Error parsing link atom element selection: " 
                << current->data
                << "\n" << endi;
                NAMD_die("Error processing QM information.");
            }
            
            linkElement = strVec[2].substr(0,2);
        }
        
        numParsedLinkElem++;
        
        if (numParsedLinkElem>qmNumBonds) {
            iout << iERROR << "More link atom elements were set than there"
            " are QM-MM bonds.\n" << endi;
            NAMD_die("Error processing QM information.");
        }
        
        int bondIter;
        
        if (strVec.size() == 1) {
            bondIter = 0;
        }
        else if (strVec.size() == 3) {
            
            Bool foundBond = false;
            
            for (bondIter=0; bondIter<qmNumBonds; bondIter++) {
                
                if ( (qmMMBond[bondIter][0] == bondAtom1 &&
                    qmMMBond[bondIter][1] == bondAtom2 ) ||
                    (qmMMBond[bondIter][0] == bondAtom2 &&
                    qmMMBond[bondIter][1] == bondAtom1 ) ) {
                    
                    foundBond = true;
                    break;
                }
            }
            
            if (! foundBond) {
                iout << iERROR << "Error parsing link atom element selection: " 
                << current->data
                << "\n" << endi;
                iout << iERROR << "Atom pair was not found to be a QM-MM bond.\n"
                << endi;
                NAMD_die("Error processing QM information.");
            }
        }
        
        strcpy(qmDummyElement[bondIter],linkElement.c_str());
        qmDummyElement[bondIter][2] = '\0';
        
        iout << iINFO << "Applying user defined link atom element "
        << qmDummyElement[bondIter] << " to QM-MM bond "
        << bondIter << ": " << qmMMBond[bondIter][0]
        << " - " << qmMMBond[bondIter][1]
        << "\n" << endi;
    }
    
    
    
    ///////////////////////////////
    /// Bond Schemes. Prepares for treatment of QM-MM bonds in ComputeQM.C
    
    
    int32 **bondsWithAtomLocal = NULL ;
    int32 *byAtomSizeLocal = NULL;
    ObjectArena <int32 >* tmpArenaLocal = NULL;
    if (qmNumBonds > 0) {
        
        bondsWithAtomLocal = new int32 *[numAtoms];
        byAtomSizeLocal = new int32[numAtoms];
        tmpArenaLocal = new ObjectArena<int32>;
        
        //  Build the bond lists
       for (int i=0; i<numAtoms; i++)
       {
         byAtomSizeLocal[i] = 0;
       }
       for (int i=0; i<numRealBonds; i++)
       {
         byAtomSizeLocal[bonds[i].atom1]++;
         byAtomSizeLocal[bonds[i].atom2]++;
       }
       for (int i=0; i<numAtoms; i++)
       {
         bondsWithAtomLocal[i] = tmpArenaLocal->getNewArray(byAtomSizeLocal[i]+1);
         bondsWithAtomLocal[i][byAtomSizeLocal[i]] = -1;
         byAtomSizeLocal[i] = 0;
       }
       for (int i=0; i<numRealBonds; i++)
       {
         int a1 = bonds[i].atom1;
         int a2 = bonds[i].atom2;
         bondsWithAtomLocal[a1][byAtomSizeLocal[a1]++] = i;
         bondsWithAtomLocal[a2][byAtomSizeLocal[a2]++] = i;
       }
    }
    
    // In this loops we try to find other bonds in which the MM atoms from
    // QM-MM bonds may be involved. The other MM atoms (which we will call M2 and M3)
    // will be involved in charge manipulation. See ComputeQM.C for details.
    for (int qmBndIt = 0; qmBndIt < qmNumBonds; qmBndIt++) {
        
        // The charge targets are accumulated in a temporary vector and then 
        // transfered to an array that will be transmited to the ComputeQMMgr object.
        std::vector<int> chrgTrgt ;
        
        int MM1 = qmMMBond[qmBndIt][0], MM2, MM2BondIndx, MM3, MM3BondIndx;
        
        switch (simParams->qmBondScheme) {
            
            case QMSCHEMERCD:
            
            case QMSCHEMECS:
            {
                // Selects ALL MM2 atoms.
                for (int i=0; i<byAtomSizeLocal[MM1]; i++) {
                    MM2BondIndx = bondsWithAtomLocal[MM1][i] ;
                    
                    // Checks which of the atoms in the bond structure is the
                    // MM2 atom.
                    if (bonds[MM2BondIndx].atom1 == MM1)
                        MM2 = bonds[MM2BondIndx].atom2;
                    else
                        MM2 = bonds[MM2BondIndx].atom1;
                    
                    // In case the bonded atom is a QM atom,
                    // skips the index.
                    if (qmAtomGroup[MM2] > 0)
                        continue;
                    
                    chrgTrgt.push_back(MM2);
                }
                
            } break;
            
            case QMSCHEMEZ3:
            {
                // Selects all MM3 atoms
                for (int i=0; i<byAtomSizeLocal[MM1]; i++) {
                    MM2BondIndx = bondsWithAtomLocal[MM1][i] ;
                    
                    // Checks which of the atoms in the bond structure is the
                    // MM2 atom.
                    if (bonds[MM2BondIndx].atom1 == MM1)
                        MM2 = bonds[MM2BondIndx].atom2;
                    else
                        MM2 = bonds[MM2BondIndx].atom1;
                    
                    // In case the bonded atom is a QM atom,
                    // skips the index.
                    if (qmAtomGroup[MM2] > 0)
                        continue;
                    
                    for (int i=0; i<byAtomSizeLocal[MM2]; i++) {
                        MM3BondIndx = bondsWithAtomLocal[MM2][i] ;
                        
                        // Checks which of the atoms in the bond structure is the
                        // MM3 atom.
                        if (bonds[MM3BondIndx].atom1 == MM2)
                            MM3 = bonds[MM3BondIndx].atom2;
                        else
                            MM3 = bonds[MM3BondIndx].atom1;
                        
                        // In case the bonded atom is a QM atom,
                        // skips the index.
                        // We also keep the search from going back to MM1.
                        if (qmAtomGroup[MM3] > 0 || MM3 == MM1)
                            continue;
                        
                        chrgTrgt.push_back(MM3);
                    }
                    
                }
                
            };
            
            case QMSCHEMEZ2:
            {
                // Selects all MM2 atoms
                for (int i=0; i<byAtomSizeLocal[MM1]; i++) {
                    MM2BondIndx = bondsWithAtomLocal[MM1][i] ;
                    
                    // Checks which of the atoms in the bond structure is the
                    // MM2 atom.
                    if (bonds[MM2BondIndx].atom1 == MM1)
                        MM2 = bonds[MM2BondIndx].atom2;
                    else
                        MM2 = bonds[MM2BondIndx].atom1;
                    
                    // In case the bonded atom is a QM atom,
                    // skips the index.
                    if (qmAtomGroup[MM2] > 0)
                        continue;
                    
                    chrgTrgt.push_back(MM2);
                }
                
            };
            
            case QMSCHEMEZ1:
            {
                // Selects all MM1 atoms
                chrgTrgt.push_back(MM1);
            } break;
        }
        
        
        qmMMChargeTarget[qmBndIt] = new int[chrgTrgt.size()] ;
        qmMMNumTargs[qmBndIt] =  chrgTrgt.size();
        
        DebugM(3, "MM-QM bond " << qmBndIt << "; MM atom " 
        << qmMMBond[qmBndIt][0] << " conections: \n" );
        
        for (size_t i=0; i < chrgTrgt.size(); i++) {
            qmMMChargeTarget[qmBndIt][i] = chrgTrgt[i];
            DebugM(3,"MM Bonded to: " << chrgTrgt[i] << "\n" );
        }
        
        chrgTrgt.clear();
    }
    
    if (bondsWithAtomLocal != NULL)
        delete [] bondsWithAtomLocal;  bondsWithAtomLocal = 0;
    if (byAtomSizeLocal != NULL)
        delete [] byAtomSizeLocal;  byAtomSizeLocal = 0;
    if (tmpArenaLocal != NULL)
        delete tmpArenaLocal;  tmpArenaLocal = 0;
    
    
    ///////////////////////////////
    /// Live Solvent Selection
    
    
    if(simParams->qmLSSOn) {
        
        std::map<Real,int> grpLSSSize ;
        std::map<Real,int>::iterator itGrpSize;
        
        qmLSSTotalNumAtms = 0;
        qmLSSResidueSize = 0;
        
        if (simParams->qmLSSFreq == 0)
            qmLSSFreq = simParams->stepsPerCycle ;
        else
            qmLSSFreq = simParams->qmLSSFreq;
            
        #ifdef DEBUG_QM
        int resSize = -1;
        #endif
        
        std::map<Real, int> grpNumLSSRes;
        std::map<Real, int>::iterator itGrpNumRes;
        
        for( auto it = lssGrpRes.begin(); it != lssGrpRes.end();it++ ) {
            
            if (it->atmIDs.size() != it->size) {
                iout << iERROR << "The number of atoms loaded for residue "
                << it->resID << " does not match the expected for this residue type.\n" 
                << endi;
                NAMD_die("Error parsing data for LSS.");
            }
            
            qmLSSTotalNumAtms += it->size;
            
            #ifdef DEBUG_QM
            if (resSize < 0) resSize = it->size ;
            if (resSize > 0 and resSize != it->size) {
                iout << iERROR << "The number of atoms loaded for residue "
                << it->resID << " does not match previously loaded residues.\n" 
                << endi;
                NAMD_die("Error parsing data for LSS.");
            }
                
//             DebugM(3,"Residue " << it->resID << ": " << it->segName
//             << " - from " << it->begAtmID << " with size "
//             << it->size << " (QM ID: " << it->qmGrpID
//             << ") has " << it->atmIDs.size() << " atoms: \n" ) ;
//             for (int i=0; i<it->atmIDs.size(); i++) 
//                 DebugM(3, it->atmIDs[i] << "\n" );
            #endif
            
            // Calculating total number of atoms per group
            itGrpSize = grpLSSSize.find(it->qmGrpID) ;
            if (itGrpSize != grpLSSSize.end())
                itGrpSize->second += it->size;
            else
                grpLSSSize.insert(std::pair<Real,int>(it->qmGrpID, it->size));
            
            // Calculating total number of solvent residues per group
            itGrpNumRes = grpNumLSSRes.find(it->qmGrpID) ;
            if (itGrpNumRes != grpNumLSSRes.end())
                itGrpNumRes->second += 1;
            else
                grpNumLSSRes.insert(std::pair<Real,int>(it->qmGrpID, 1));
        }
        
        qmLSSResidueSize = lssGrpRes.begin()->size ;
        
        qmLSSSize = new int[qmNumGrps];
        
        qmLSSIdxs = new int[qmLSSTotalNumAtms];
        int lssAtmIndx = 0;
        
        switch (simParams->qmLSSMode) {
        
        case QMLSSMODECOM:
        {
            
            qmLSSRefSize = new int[qmNumGrps];
            
            qmLSSMass = new Mass[qmLSSTotalNumAtms];
            
            qmLSSRefIDs = new int[totalNumRefAtms];
            qmLSSRefMass = new Mass[totalNumRefAtms];
            int lssRefIndx = 0;
            
            for (int grpIndx=0; grpIndx<qmNumGrps; grpIndx++) {
                
                itGrpSize = grpNumLSSRes.find(qmGrpID[grpIndx]) ;
                
                if (itGrpSize != grpNumLSSRes.end())
                    qmLSSSize[grpIndx] =  itGrpSize->second;
                else
                    qmLSSSize[grpIndx] =  0;
                
                for( auto it = lssGrpRes.begin(); it != lssGrpRes.end();it++ ) {
                    
                    if (it->qmGrpID == qmGrpID[grpIndx]) {
                        for (int i=0; i<it->atmIDs.size(); i++) {
                            qmLSSIdxs[lssAtmIndx] = it->atmIDs[i];
                            qmLSSMass[lssAtmIndx] = atoms[it->atmIDs[i]].mass;
                            lssAtmIndx++;
                        }
                    }
                }
                
                DebugM(4, "QM group " << qmGrpID[grpIndx]
                << " has " << qmLSSSize[grpIndx] << " solvent molecules, "
                << "totalling " << grpLSSSize[qmGrpID[grpIndx]]
                << " atoms (check " << lssAtmIndx << ").\n" );
                
                qmLSSRefSize[grpIndx] = lssGrpRefIDs[qmGrpID[grpIndx]].size();
                for(int i=0; i < qmLSSRefSize[grpIndx]; i++) {
                    qmLSSRefIDs[lssRefIndx] = lssGrpRefIDs[qmGrpID[grpIndx]][i];
                    qmLSSRefMass[lssRefIndx] =  atoms[qmLSSRefIDs[lssRefIndx]].mass;
                    lssRefIndx++;
                }
                
                DebugM(4, "QM group " << qmGrpID[grpIndx]
                << " has " << qmLSSRefSize[grpIndx] << " non-solvent atoms for LSS.\n" );
            }
            
        } break ;
        
        case QMLSSMODEDIST:
        {
            for (int grpIndx=0; grpIndx<qmNumGrps; grpIndx++) {
                
                itGrpSize = grpNumLSSRes.find(qmGrpID[grpIndx]) ;
                
                if (itGrpSize != grpNumLSSRes.end())
                    qmLSSSize[grpIndx] =  itGrpSize->second;
                else
                    qmLSSSize[grpIndx] =  0;
                
                for( auto it = lssGrpRes.begin(); it != lssGrpRes.end();it++ ) {
                    
                    if (it->qmGrpID == qmGrpID[grpIndx]) {
                        for (int i=0; i<it->atmIDs.size(); i++) {
                            qmLSSIdxs[lssAtmIndx] = it->atmIDs[i];
                            lssAtmIndx++;
                        }
                    }
                }
                
                DebugM(4, "QM group " << qmGrpID[grpIndx]
                << " has " << qmLSSSize[grpIndx] << " solvent molecules, "
                << "totalling " << grpLSSSize[qmGrpID[grpIndx]]
                << " atoms (check " << lssAtmIndx << ").\n" );
            }
            
        } break ;
            
        }
    }
    
    
    ///////////////////////////////
    /// Custom Point Charge selection
    
    
    PDB *customPCPDB;
    
    // In case we have a custom and fixed set of point charges for each QM group,
    // we process the files containing information.
    current = NULL;
    if (simParams->qmCustomPCSel) {
        current = cfgList->find("QMCustomPCFile");
    }
    
    std::map<Real,std::vector<int> > qmPCVecMap ;
    
    int numParsedPBDs = 0;
    for ( ; current != NULL; current = current->next ) {
        
        iout << iINFO << "Parsing QM Custom PC file " << current->data << "\n" << endi;
        customPCPDB = new PDB(current->data);
        
        if (customPCPDB->num_atoms() != numAtoms)
           NAMD_die("Number of atoms in QM Custom PC PDB file doesn't match coordinate PDB");
        
        std::vector< int > currPCSel ;
        Real currQMID = 0 ;
        int currGrpSize = 0 ;
        
        for (int atmInd = 0 ; atmInd < numAtoms; atmInd++) {
            
            BigReal beta = customPCPDB->atom(atmInd)->temperaturefactor() ;
            BigReal occ = customPCPDB->atom(atmInd)->occupancy() ;
            
            if ( beta != 0 && occ != 0)
                NAMD_die("An atom cannot be marked as QM and poitn charge simultaneously!");
            
            // If this is not a QM atom and 
            if (occ != 0) {
                currPCSel.push_back(atmInd) ;
            }
            
            if (beta != 0) {
                if (pdbP->atom(atmInd)->temperaturefactor() != beta)
                    NAMD_die("QM Group selection is different in reference PDB and Custom-PC PDB!");
                
                if (currQMID == 0) {
                    // If this is the first QM atom we find, keep the QM Group ID.
                    currQMID = beta;
                }
                else {
                    // For all other atoms, check if it is the same group. It must be!!
                    if (currQMID != beta)
                        NAMD_die("Found two different QM group IDs in this file!");
                }
                
                currGrpSize++;
            }
            
        }
        
        if (currGrpSize != qmGrpSizeVec[ qmGrpIDMap[currQMID] ])
            NAMD_die("Reference PDB and Custom-PC PDB have different QM group sizes!") ;
        
        qmPCVecMap.insert(std::pair<Real,std::vector<int>>(
            currQMID, currPCSel ));
        
        numParsedPBDs++;
        delete customPCPDB;
    }
    
    delete pdbP;
    
    if (numParsedPBDs != qmNumGrps && simParams->qmCustomPCSel) {
        iout << iWARN << "The number of files provided for custom point "
        "charges is not equal to the number of QM groups!\n" << endi;
    }
    
    // Initializes an array with the number of Custom Point Charges per 
    // QM group.
    qmCustPCSizes = new int[qmNumGrps];
    for (int i=0; i<qmNumGrps; i++)
        qmCustPCSizes[i] = 0;
    
    qmTotCustPCs = 0;
    
    // Stores the size of each Custom PC vector in the array.
    // We may not have PCs for all QM groups.
    for (auto it = qmPCVecMap.begin(); it != qmPCVecMap.end(); it++) {
        qmTotCustPCs += (*it).second.size();
        int qmIndx = qmGrpIDMap[(*it).first];
        qmCustPCSizes[qmIndx] = (*it).second.size();
    }
    
    qmCustomPCIdxs = new int[qmTotCustPCs];
    
    if (simParams->qmCustomPCSel) {
        
        int auxIter = 0;
        for (int grpIndx=0; grpIndx<qmNumGrps; grpIndx++) {
            
            Real currQMID = qmGrpID[grpIndx];
            
            iout << iINFO << "Loading " << qmPCVecMap[currQMID].size()
            << " custom point charges to QM Group " << grpIndx
            << " (ID: " << currQMID << ")\n" << endi;
            
            for (int pcIndx=0; pcIndx<qmPCVecMap[currQMID].size(); pcIndx++) {
                qmCustomPCIdxs[auxIter] = qmPCVecMap[currQMID][pcIndx] ;
                auxIter++;
            }
        }
    } 
    
    
    ///////////////////////////////
    /// Topology preparation
    
    // Modifies Atom charges.
    // QM atoms cannot have charges in the standard location, to keep
    // NAMD from calculating electrostatic interactions between QM and MM atoms.
    // We handle electrostatics ourselves in ComputeQM.C and in special
    // modifications for PME.
    for (int i=0; i<qmNumQMAtoms; i++) {
        qmAtmChrg[i] = atoms[qmAtmIndx[i]].charge;
        atoms[qmAtmIndx[i]].charge = 0;
    }
    
    return;
    
#endif // MEM_OPT_VERSION
}

// Adapted from Molecule::delete_alch_bonded
void Molecule::delete_qm_bonded(void)  {
    
#ifdef MEM_OPT_VERSION
    NAMD_die("QMMM interface is not supported in memory optimized builds.");
#else

    DebugM(3,"Cleaning QM bonds, angles, dihedrals, impropers and crossterms.\n")
    
    // Bonds
    Bond *nonQMBonds;
    nonQMBonds = new Bond[numBonds] ;
    int nonQMBondCount = 0;
    qmDroppedBonds = 0;
    for (int i = 0; i < numBonds; i++) {
        int part1 = qmAtomGroup[bonds[i].atom1];
        int part2 = qmAtomGroup[bonds[i].atom2];
        if (part1 > 0 && part2 > 0 ) {
            qmDroppedBonds++;
        } else {
            // Just a simple sanity check.
            // If there are no QM bonds defined by the user but we find a
            // bond between a QM and an MM atom, error out.
            if ((part1 > 0 || part2 > 0) && qmNumBonds == 0) {
                iout << iERROR << " Atoms " << bonds[i].atom1
                << " and " << bonds[i].atom2 << " form a QM-MM bond!\n" << endi ;
                NAMD_die("No QM-MM bonds were defined by the user, but one was found.");
            }
            nonQMBonds[nonQMBondCount++] = bonds[i];
        }
    }
    numBonds = nonQMBondCount;
    delete [] bonds;
    bonds = new Bond[numBonds];
    for (int i = 0; i < nonQMBondCount; i++) {
        bonds[i]=nonQMBonds[i];
    }
    delete [] nonQMBonds;
  numRealBonds = numBonds;
    
  // Angles
  Angle *nonQMAngles;
  nonQMAngles = new Angle[numAngles];
  int nonQMAngleCount = 0;
  qmDroppedAngles = 0;
  for (int i = 0; i < numAngles; i++) {
    int part1 = qmAtomGroup[angles[i].atom1];
    int part2 = qmAtomGroup[angles[i].atom2];
    int part3 = qmAtomGroup[angles[i].atom3];
    if (part1 > 0 && part2 > 0 && part3 > 0) {
      //CkPrintf("-----ANGLE ATOMS %i %i %i partitions %i %i %i\n",angles[i].atom1, angles[i].atom2, angles[i].atom3, part1, part2, part3);
      qmDroppedAngles++;
    }
    else {
      nonQMAngles[nonQMAngleCount++] = angles[i];
    }
  }
  numAngles = nonQMAngleCount;
  delete [] angles;
  angles = new Angle[numAngles];
  for (int i = 0; i < nonQMAngleCount; i++) {
    angles[i]=nonQMAngles[i];
  }
  delete [] nonQMAngles;


  // Dihedrals
  Dihedral *nonQMDihedrals;
  nonQMDihedrals = new Dihedral[numDihedrals];
  int nonQMDihedralCount = 0;
  qmDroppedDihedrals = 0;
  for (int i = 0; i < numDihedrals; i++) {
    int part1 = qmAtomGroup[dihedrals[i].atom1];
    int part2 = qmAtomGroup[dihedrals[i].atom2];
    int part3 = qmAtomGroup[dihedrals[i].atom3];
    int part4 = qmAtomGroup[dihedrals[i].atom4];
    if (part1 > 0 && part2 > 0 && part3 > 0 && part4 > 0) {
      //CkPrintf("-----i %i DIHEDRAL ATOMS %i %i %i %i partitions %i %i %i %i\n",i,dihedrals[i].atom1, dihedrals[i].atom2, dihedrals[i].atom3, dihedrals[i].atom4, part1, part2, part3,part4);
      qmDroppedDihedrals++;
    }
    else {
      nonQMDihedrals[nonQMDihedralCount++] = dihedrals[i];
    }
  }
  numDihedrals = nonQMDihedralCount;
  delete [] dihedrals;
  dihedrals = new Dihedral[numDihedrals];
  for (int i = 0; i < numDihedrals; i++) {
    dihedrals[i]=nonQMDihedrals[i];
  }
  delete [] nonQMDihedrals;

  // Impropers
  Improper *nonQMImpropers;
  nonQMImpropers = new Improper[numImpropers];
  int nonQMImproperCount = 0;
  qmDroppedImpropers = 0;
  for (int i = 0; i < numImpropers; i++) {
    int part1 = qmAtomGroup[impropers[i].atom1];
    int part2 = qmAtomGroup[impropers[i].atom2];
    int part3 = qmAtomGroup[impropers[i].atom3];
    int part4 = qmAtomGroup[impropers[i].atom4];
    if (part1 > 0 && part2 > 0 && part3 > 0 && part4 > 0) {
      //CkPrintf("-----i %i IMPROPER ATOMS %i %i %i %i partitions %i %i %i %i\n",i,impropers[i].atom1, impropers[i].atom2, impropers[i].atom3, impropers[i].atom4, part1, part2, part3,part4);
      qmDroppedImpropers++;
    }
    else {
      nonQMImpropers[nonQMImproperCount++] = impropers[i];
    }
  }
  numImpropers = nonQMImproperCount;
  delete [] impropers;
  impropers = new Improper[numImpropers];
  for (int i = 0; i < numImpropers; i++) {
    impropers[i]=nonQMImpropers[i];
  }
  delete [] nonQMImpropers;
  
  // Crossterms 
  Crossterm *nonQMCrossterms;
  nonQMCrossterms = new Crossterm[numCrossterms];
  int nonQMCrosstermCount = 0;
  qmDroppedCrossterms = 0 ;
  for (int i = 0; i < numCrossterms; i++) {
    int part1 = qmAtomGroup[crossterms[i].atom1];
    int part2 = qmAtomGroup[crossterms[i].atom2];
    int part3 = qmAtomGroup[crossterms[i].atom3];
    int part4 = qmAtomGroup[crossterms[i].atom4];
    int part5 = qmAtomGroup[crossterms[i].atom5];
    int part6 = qmAtomGroup[crossterms[i].atom6];
    int part7 = qmAtomGroup[crossterms[i].atom7];
    int part8 = qmAtomGroup[crossterms[i].atom8];
    if (part1 > 0 && part2 > 0 && part3 > 0 && part4 > 0 &&
        part5 > 0 && part6 > 0 && part7 > 0 && part8 > 0) {
      //CkPrintf("-----i %i IMPROPER ATOMS %i %i %i %i partitions %i %i %i %i\n",i,impropers[i].atom1, impropers[i].atom2, impropers[i].atom3, impropers[i].atom4, part1, part2, part3,part4);
      qmDroppedCrossterms++;
    }
    else {
      nonQMCrossterms[nonQMCrosstermCount++] = crossterms[i];
    }
  }
  numCrossterms = nonQMCrosstermCount;
  delete [] crossterms;
  crossterms = new Crossterm[numCrossterms] ;
  for (int i=0; i < numCrossterms; i++){
      crossterms[i] = nonQMCrossterms[i] ;
  }
  delete [] nonQMCrossterms ;
  
  if (!CkMyPe()) {
      iout << iINFO << "The QM region will remove " << qmDroppedBonds << " bonds, " << 
          qmDroppedAngles << " angles, " << qmDroppedDihedrals << " dihedrals, "
          << qmDroppedImpropers << " impropers and " << qmDroppedCrossterms 
          << " crossterms.\n" << endi ;
  }
  
#endif // MEM_OPT_VERSION
}

