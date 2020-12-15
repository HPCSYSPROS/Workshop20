/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "vmdsock.h"
#include "Node.h"
#include "IMDOutput.h"
#include "imd.h"
#include "SimParameters.h"
#include "UniqueSortedArray.h"
#include "GlobalMaster.h"
#include "GlobalMasterIMD.h"
#include "Vector.h"

#include <errno.h>

//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"

struct vmdforce {
  int index;
  Vector force;
  int operator <(const vmdforce& v) {return index < v.index;}
  // XXX the following is an abuse of overloading!
  int operator ==(const vmdforce& v) {return index == v.index;}
  vmdforce& operator=(const vmdforce& v) {
    index=v.index;
    force=v.force; 
    return *this;
  }  
};

//
// XXX static and global variables are unsafe for shared memory builds.
// The use of global and static vars should be eliminated.
//
static UniqueSortedArray<vmdforce> vmdforces;

// Search for a free port in the range 1025-4096; return the successful port,
// or -1 if we failed.

static int find_free_port(void *sock, int defport) {
  if (vmdsock_bind(sock, defport)==0) return defport; // success
  for (int port=1025; port < 4096; port++) 
    if (vmdsock_bind(sock, port)==0) return port;
  return -1;
}
 
GlobalMasterIMD::GlobalMasterIMD() {
  DebugM(3,"Constructing\n");
  SimParameters *simparams = Node::Object()->simParameters;
  int port = simparams->IMDport;
  IMDwait = simparams->IMDwait;
  IMDignore = simparams->IMDignore;
  coordtmp = NULL;
  coordtmpsize = 0;

  if ( vmdsock_init() ) {
    NAMD_die("Unable to initialize socket interface for IMD.\n");
  }
  sock = vmdsock_create();
  int newport = find_free_port(sock, port);
  if (newport != port) {
    iout << iWARN << "Interactive MD failed to bind to port "
                  << port << ".\n" << endi;
  }
  if (newport < 0) {
    vmdsock_destroy(sock);
    NAMD_die("Interactive MD failed to find free port.\n");
  }
  vmdsock_listen(sock); 
  iout << iINFO << "Interactive MD listening on port "
                  << newport << ".\n" << endi;
  DebugM(2,"Done constructing ("<<requestedGroups().size()<<" initial groups)\n");

  Node::Object()->imd->use_imd(this);
}

GlobalMasterIMD::~GlobalMasterIMD() {
  if (sock) 
    vmdsock_destroy(sock);
  for (int i=0; i<clients.size(); i++)
    vmdsock_destroy(clients[i]);
  delete [] coordtmp;
}

static int my_imd_connect(void *s) {
  if (imd_handshake(s)) {
    iout << iWARN << "IMD handshake failed\n" << endi;
    return 0;
  }

  // Wait a second, then see if VMD has responded.
  int32 length;
  if (vmdsock_selread(s,1) != 1 || imd_recv_header(s, &length) != IMD_GO) {
    iout << iWARN << "Incompatible Interactive MD, use VMD v1.4b2 or higher\n"
         << endi;
    return 0;
  }
  return 1;
}

void GlobalMasterIMD::calculate() {
  /* clear out the requested forces first! */
  if (!IMDignore) {
    modifyAppliedForces().resize(0);
    modifyForcedAtoms().resize(0);
    modifyGroupForces().resize(0);
  }

  // check for incoming connection
  do {
  int rc;
  if (IMDwait && !clients.size()) {
    iout << iINFO << "INTERACTIVE MD AWAITING CONNECTION\n" << endi;
    do { rc = vmdsock_selread(sock, 3600); } while (rc <= 0);
  } else {
    rc = vmdsock_selread(sock, 0);
  } 
  if (rc > 0) {
    void *clientsock = vmdsock_accept(sock);
    if (!clientsock) {
      iout << iWARN << "IMD socket accept failed\n" << endi;
    } else {
      if (!my_imd_connect(clientsock)) {
        iout << iWARN << "IMD connection failed\n" << endi;
        vmdsock_destroy(clientsock);
      } else {
        iout << iINFO << "IMD connection opened\n" <<endi;	
        clients.add(clientsock);
      }
    }
  }
  } while (IMDwait && !clients.size());

  // Assume for now that the only thing we get from VMD is a set of forces.
  // Later we'll want to look for and implement more sophisticated control
  // parameters. ie have a specified protocol

  // Check/get new forces from VMD
  get_vmd_forces();

  // Right now I don't check to see if any new forces were obtained.
  // An optimization would be cache the results message.  However, there
  // would still be copying since it looks like the messages get deleted
  // by the receiver.

  /* set our arrays to be big enough to hold all of the forces */
  int num = vmdforces.size();

  DebugM(2,"Setting " << num << " forces.\n");
  
  if (!IMDignore) {
    modifyForcedAtoms().resize(num);
    modifyAppliedForces().resize(num);
  
    int i;
    UniqueSortedArray<vmdforce>::iterator v_i = vmdforces.begin();
    for ( i = 0; i < num; ++i, ++v_i) {
      modifyForcedAtoms().item(i) = v_i->index;
      modifyAppliedForces().item(i) = v_i->force;
    }
  }
}

void GlobalMasterIMD::get_vmd_forces() {
  IMDType type;
  int32 length;
  int32 *vmd_atoms;
  float *vmd_forces;
  int paused = 0;
  int warned = 0;
  vmdforce *vtest, vnew;

  // Loop through each socket one at a time.  By doing this, rather than 
  // polling all sockets at once, NAMD only has to keep up with one IMD
  // connection; if it tried to read from all of them, it could more easily
  // fall behind and never finish draining all the sockets.
  // It would be better to have a system where VMD couldn't DOS NAMD by 
  // spamming it with messages, but in practice NAMD is able to keep up with 
  // VMD's send rate.
  for (int i_client=0; i_client<clients.size(); i_client++) {
    void *clientsock = clients[i_client];
    while (vmdsock_selread(clientsock,0) > 0 || paused) {  // Drain the socket
      type = imd_recv_header(clientsock, &length);
      switch (type) {
        case IMD_MDCOMM:
          // Expect the msglength to give number of indicies, and the data
          // message to consist of first the indicies, then the coordinates
          // in xyz1 xyz2... format.
          vmd_atoms = new int32[length];
          vmd_forces = new float[3*length];
          if (imd_recv_mdcomm(clientsock, length, vmd_atoms, vmd_forces)) {
            iout << iWARN <<
              "Error reading IMD forces, killing connection\n" << endi;
            goto vmdDestroySocket;
          } 
          if (IMDignore) {
            if ( ! warned ) {
              warned = 1;
              iout << iWARN << "Ignoring IMD forces due to IMDignore\n" << endi;
            }
          } else {
            for (int i=0; i<length; i++) {
              vnew.index = vmd_atoms[i];
              if ( (vtest=vmdforces.find(vnew)) != NULL) {
                // find was successful, so overwrite the old force values
                if (vmd_forces[3*i] != 0.0f || vmd_forces[3*i+1] != 0.0f
                    || vmd_forces[3*i+2] != 0.0f) {
                  vtest->force.x = vmd_forces[3*i];
                  vtest->force.y = vmd_forces[3*i+1];
                  vtest->force.z = vmd_forces[3*i+2];
                } else {
                  // or delete it from the list if the new force is ZERO
                  vmdforces.del(vnew);
                }
              }
              else {
                // Create a new entry in the table if the new force isn't ZERO
                if (vmd_forces[3*i] != 0.0f || vmd_forces[3*i+1] != 0.0f
                    || vmd_forces[3*i+2] != 0.0f) {
                  vnew.force.x = vmd_forces[3*i];
                  vnew.force.y = vmd_forces[3*i+1];
                  vnew.force.z = vmd_forces[3*i+2];
                  vmdforces.add(vnew);
                }
              }
            } 
          }
          delete [] vmd_atoms;
          delete [] vmd_forces;
          break;
        case IMD_TRATE:
          iout << iINFO << "Setting transfer rate to " << length<<'\n'<<endi;	
          Node::Object()->imd->set_transrate(length);
          break;
        case IMD_PAUSE:
          if (IMDignore) {
            iout << iWARN << "Ignoring IMD pause due to IMDignore\n" << endi;
            break;
          }
          if ( paused ) {
            iout << iINFO << "Resuming IMD\n" << endi;
            IMDwait = Node::Object()->simParameters->IMDwait;
          }
          paused = ! paused;
          if ( paused ) {
            iout << iINFO << "Pausing IMD\n" << endi;
            IMDwait = 1;
          }
          break;
        case IMD_IOERROR:
          iout << iWARN << "IMD connection lost\n" << endi;
        case IMD_DISCONNECT:
          iout << iINFO << "IMD connection detached\n" << endi;
          vmdDestroySocket:
          vmdsock_destroy(clientsock);
          clients.del(i_client);
          goto vmdEnd;
        case IMD_KILL:
          if (IMDignore) {
            iout << iWARN << "Ignoring IMD kill due to IMDignore\n" << endi;
            break;
          }
          NAMD_quit("Received IMD kill from client\n");
          break;
        case IMD_ENERGIES:
          IMDEnergies junk;
          imd_recv_energies(clientsock, &junk);
          break;
        case IMD_FCOORDS:
          vmd_forces = new float[3*length];
          imd_recv_fcoords(clientsock, length, vmd_forces);
          delete [] vmd_forces;
          break;
        default: ;
      }
    }
  vmdEnd: ;
  }
}

void GlobalMasterIMD::send_energies(IMDEnergies *energies) {
  for (int i=0; i<clients.size(); i++) {
    void *clientsock = clients[i];
    if (!clientsock || !vmdsock_selwrite(clientsock,0)) continue;
    imd_send_energies(clientsock, energies);
  }
}

void GlobalMasterIMD::send_fcoords(int N, FloatVector *coords) {
  for (int i=0; i<clients.size(); i++) {
    void *clientsock = clients[i];
    if (!clientsock || !vmdsock_selwrite(clientsock,0)) continue;
    if (sizeof(FloatVector) == 3*sizeof(float)) {
      imd_send_fcoords(clientsock, N, (float *)coords);
    } else {
      if (coordtmpsize < N) {
        delete [] coordtmp;
        coordtmp = new float[3*N];
        coordtmpsize = N;
      }
      for (int i=0; i<N; i++) {
        coordtmp[3*i] = coords[i].x; 
        coordtmp[3*i+1] = coords[i].y; 
        coordtmp[3*i+2] = coords[i].z; 
      } 
      imd_send_fcoords(clientsock, N, coordtmp);
    }
  }
}
