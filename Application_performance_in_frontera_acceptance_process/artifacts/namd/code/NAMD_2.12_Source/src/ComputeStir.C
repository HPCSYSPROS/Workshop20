/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000, 2001 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/* Barry Isralewitz, July 2001 */

#include "InfoStream.h" 
#include "ComputeStir.h"
#include "Node.h"
#include "SimParameters.h"
#include "HomePatch.h"
#include "DataStream.h"
#include "Molecule.h"
#include "CollectionMgr.h" //for error corr of buffer problem


ComputeStir::ComputeStir(ComputeID c, PatchID pid)
  : ComputeHomePatch(c,pid)
{

  

  
  SimParameters* simParams = Node::Object()->simParameters;
  axisUnit = simParams->stirAxis.unit();
  pivot = simParams->stirPivot;
  
  
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  
  //nitVars(p, axisDir);
  initVars();
}
/*			END OF FUNCTION ComputeStir		*/


ComputeStir::~ComputeStir()

{
delete reduction;
}
/*			END OF FUNCTION ~ComputeStir		*/


Tensor ComputeStir::matMult (Tensor a, Tensor b) {
  // should really use operator*(Tensor, Tensor)(, but this should be rarely
  // used, maybe inline function.  Or extendor Tensor class.
  // Also abusing Tensor math concept to mere 3x3 matrix
  // just so I can borrow a method from Tensor class. 
  Tensor tmp;

  //multiply two 3x3 matrices...

  tmp.xx = a.xx*b.xx + a.xy*b.yx + a.xz*b.zx;
  tmp.xy = a.xx*b.xy + a.xy*b.yy + a.xz*b.zy;
  tmp.xz = a.xx*b.xz + a.xy*b.yz + a.xz*b.zz;

  tmp.yx = a.yx*b.xx + a.yy*b.yx + a.yz*b.zx;
  tmp.yy = a.yx*b.xy + a.yy*b.yy + a.yz*b.zy;
  tmp.yz = a.yx*b.xz + a.yy*b.yz + a.yz*b.zz;

  tmp.zx = a.zx*b.xx + a.zy*b.yx + a.zz*b.zx;
  tmp.zy = a.zx*b.xy + a.zy*b.yy + a.zz*b.zy;
  tmp.zz = a.zx*b.xz + a.zy*b.yz + a.zz*b.zz;
  
  return tmp;
}


void ComputeStir::initVars () {
  //find leftMat and rightMat, using 3x3 arrrays for speed,
  //so stil must provide translation before/after matrix multuplies
  BigReal a,b,c,d;
  Tensor testA, testB;
  Tensor rotXMat, invRotXMat, rotYMat, invRotYMat;
  Tensor temp;
  Vector refPos;
  SimParameters *simParams = Node::Object()->simParameters;
  Molecule *molecule = Node::Object()->molecule;
  char statbuf[1024];
  const int globalNumAtoms = molecule->numAtoms;



 

  //iout << "DEBUG: Starting Init ComputeStir Vars...\n" << endi;
  omega = ( (simParams->stirVel) /  360) * (2 * PI);
  startingTheta  = ( (simParams->stirStartingTheta) /  360) * (2 * PI);
  //iout <<"DEBUG: omega(deg/ts,stirVel)= "<< simParams->stirVel <<"  omega(rad/ts)= "<<omega << "  startingTheta(deg,stirStartingTheta)= "<< simParams->stirStartingTheta <<"  startingTheta(rad)= "<< startingTheta << "\n"<<endi;

  //rotation speed in radians/timestep

  a = axisUnit.x; b = axisUnit.y; c = axisUnit.z;
  d = sqrt( (b * b) + (c * c) );
  //iout <<"DEBUG: a="<<a<<" b= "<<b<< "c= "<<c<<" d= "<<d<<"\n"<<endi;

  testA.xx = 1.1; testA.xy = 3; testA.xz = 3;
  //printTensor(testA);  
  //iout << "DEBUG:   axis unit vector = ("  << endi;
  //iout  << axisUnit.x << "," <<    axisUnit.y << "," << axisUnit.z << ")\n" << endi;
  ////findRotMats
  //yes, could pre-mupltiply these matrices,
  //but makes code easier to check
  //Besides, only done at startup, so speed not an issue
  rotXMat.xx = 1;  rotXMat.xy =   0;  rotXMat.xz =    0;
  rotXMat.yx = 0;  rotXMat.yy = c/d;  rotXMat.yz = -b/d;
  rotXMat.zx = 0;  rotXMat.zy = b/d;  rotXMat.zz =  c/d;

  
  invRotXMat.xx =  1; invRotXMat.xy =     0;  invRotXMat.xz =   0;
  invRotXMat.yx =  0; invRotXMat.yy =   c/d;  invRotXMat.yz = b/d;
  invRotXMat.zx =  0; invRotXMat.zy =  -b/d;  invRotXMat.zz = c/d;


  rotYMat.xx = d; rotYMat.xy = 0; rotYMat.xz = -a;
  rotYMat.yx = 0; rotYMat.yy = 1;  rotYMat.yz = 0;
  rotYMat.zx = a; rotYMat.zy = 0;  rotYMat.zz = d;


  invRotYMat.xx = d;   invRotYMat.xy =  0;      invRotYMat.xz = a;
  invRotYMat.yx = 0;   invRotYMat.yy =  1;      invRotYMat.yz = 0;
  invRotYMat.zx =-a;   invRotYMat.zy =  0;      invRotYMat.zz = d;

  temp = rotYMat;
  //iout << "DEBUG: Mat check rotYMat\n" << endi;
  //printTensor (temp);
 
  temp = invRotYMat;
  //iout << "DEBUG: Mat check invRotYMat\n" << endi;
  //printTensor (temp);

  temp = matMult (invRotYMat, rotYMat);
  //iout << "DEBUG: Mat check Y^-1 Y\n" << endi;
  //printTensor (temp);

    temp = matMult (invRotXMat, rotXMat);
  //iout << "DEBUG: Mat check X^-1 X\n" << endi;
  //printTensor (temp);
  

  //tmp.xx =  ; tmp.xy =   ;  tmp.xz = 
  //tmp.yx = ;  tmp.yy =  ;  tmp.yz =
  //tmp.zx =  ; tmp.zy =  ;  tmp.zz = 
 
  leftMat = matMult (invRotXMat, invRotYMat);
  //translate after use
  //e.g. :  matVec(leftMat,p) ) - offset
  rightMat = matMult ( rotYMat, rotXMat) ;
  //translate before use
  //e.g. :  matVec(rightMat,p+offset))

  temp =  leftMat;
  //iout << "DEBUG: Mat check leftMat \n" << endi;
  //printTensor (temp); 

     temp = rightMat;
  //iout << "DEBUG: Mat check rightMat \n" << endi;
  //printTensor (temp); 

  temp = matMult (leftMat,rightMat);
  //iout << "DEBUG: Mat check matMult (leftMat,rightMat) \n" << endi;
  //printTensor (temp); 
  
  temp = arbRotMat (0.34);
  //iout << "DEBUG: Mat check arbRotMat (0.34) \n" << endi;
  //printTensor (temp);   

  //now make table of initial thetas
  // (perhaps could do this back in Molecule.C)
  for (int i=0; i<globalNumAtoms; i++) {
      if (molecule->is_atom_stirred(i))
      {
	//CkPrintf ("DEBUG: now to get stir params");	
	
	//use same index as params
	molecule->get_stir_refPos (refPos, i);

	molecule->put_stir_startTheta (findTheta(refPos), i); 
	//CkPrintf ("DEBUG: now to get sprintf x y z's");
	//iout << "DEBUG: For atom i="<<i<<" refPos.x= " << refPos.x << " refPos.y= "<<refPos.y<<" refPos.z "<< refPos.z  <<  " theta= " << molecule->get_stir_startTheta(i)<<"\n"<<endi; 
	  // CollectionMgr::Object()->sendDataStream(statbuf);	
      
      }

  }
}
  

BigReal ComputeStir::findTheta (Vector pos) {
  BigReal theta;
  //translate to origin and align to z axis
  Vector aligned_pos = rightMat * (pos - pivot);
  
  if ( ! ((aligned_pos.x == 0 ) && (aligned_pos.y  == 0) )) {
    //we have our rotation counter-clockwise around z-axis (right hand rule
    //along increasing z
    theta = atan (aligned_pos.y/aligned_pos.x);
    if (aligned_pos.x < 0) {
      theta = theta + PI;
      }
  } else { 
    //iout << iWARN << "stirred atom exactly along stir axis, setting theta of this atom to 0.0" << endi;
    theta = 0.0;
  }
  
  
  return theta;
}

Tensor ComputeStir::arbRotMat (BigReal theta) {
  //needed only to show graphics in output
  //not called by main doForce
  //must use +pivot, -pivot after/before 
  Tensor rotZMat;
  Tensor temp;
  rotZMat.xx = cos (theta) ; rotZMat.xy = -sin (theta);  rotZMat.xz = 0;
  rotZMat.yx = sin (theta);  rotZMat.yy =  cos(theta);   rotZMat.yz = 0;
  rotZMat.zx =  0;           rotZMat.zy =      0;        rotZMat.zz = 1;
  //iout <<"DEBUG:: in arbRotMat, print Tensor rotZMat\n"<<endi;
  //printTensor (rotZMat);
  //iout <<"DEBUG:: in arbRotMat, print Tensor rotZMat*rightMat\n"<<endi;
  temp =matMult (rotZMat, rightMat);
  //printTensor (temp);

  temp = matMult (leftMat,rightMat);
  //iout << "DEBUG: Mat check matMult (leftMat,rightMat) \n" << endi;
  //printTensor (temp); 

  temp = matMult (leftMat, (matMult (rotZMat, rightMat) ) ) ;
;
 //iout <<"DEBUG:: in arbRotMat, print Tensor temp(should be whats returned)\n"<<endi;
 //printTensor (temp); 



  return matMult (leftMat, (matMult (rotZMat, rightMat) ) );
}

void ComputeStir::printTensor (Tensor &t) {
  //just for debugging
  iout << "DEBUG: Tensor =" << "\n" << endi;
  iout <<"DEBUG:   " << t.xx << " " <<t.xy  << " " << t.xz << \
    "\nDEBUG:   " << t.yx << " " << t.yy << " " << t.yz << \
    "\nDEBUG:   " << t.zx << " " << t.zy << " " <<t.zz << "\n"  << endi;
}

BigReal ComputeStir::distanceToRay(Vector dir,Vector origin,Vector x) {
          return (x - projectionOnRay(dir, origin,x)).length();
}
BigReal ComputeStir::findHeight(Vector x) {
  //rotate so axisUnitVec is along z-axis
  //we have our rotation theta counter-clockwise around z (right hand rule along increasing z)
  Vector pos =rightMat * (x - pivot);  
  //so, just read height as z
  return pos.z;
};

Vector ComputeStir::projectionOnRay (Vector dir, Vector origin, Vector x) {
  return origin + ( ( ((x - origin ).dot(dir) )/    \
		     (dir.dot(dir))                )
		               * dir);
};

Vector ComputeStir::placeThetaRadius (BigReal theta, BigReal height, BigReal radius) {
  Vector a; 
  Tensor rotZMat;
  ///iout <<"DEBUG:: in placeThetaRadius,  theta = " << theta << " height= "<<height<<" radius= \n"<< radius<<endi;
  rotZMat.xx = cos (theta) ; rotZMat.xy = -sin (theta);  rotZMat.xz = 0;
  rotZMat.yx = sin (theta);  rotZMat.yy =  cos(theta);   rotZMat.yz = 0;
  rotZMat.zx =  0;           rotZMat.zy =      0;        rotZMat.zz = 1;
  ///iout <<"DEBUG:: in placeThetaRadius, print Tensor rotZMat\n"<<endi;
  ///printTensor (rotZMat);
  ///iout <<"DEBUG:: in placeThetaRadius, pivot = "<<pivot<<"\n"<<endi;
  //set position to (radius,0,height)  
  a.x = radius; a.y =0; a.z = height;
		  ///  //move position along z
		  ///  a += Vector(0,0,height);
  ///iout <<"DEBUG:: in placeThetaRadius, prerot,, a.x="<<a.x<<" a.y= "<<a.y<<" a.z= "<<a.z<<endi;
  //rotate it by theta, then rotate it to axis system, then transl. to axis sytem
  a = (matMult(leftMat, rotZMat) * a) +pivot;
  
   ///iout <<"DEBUG:: in placeThetaRadius, a.x="<<a.x<<" a.y= "<<a.y<<" a.z= "<<a.z<<endi;
  return a;
}
void ComputeStir::doForce (FullAtom* p, Results* r) {
  
  
  Molecule *molecule = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  Lattice &lattice = homePatch->lattice;
  char statbuf[1024];
  Vector curPos; //current atom position
  Vector proj;  //projection on the axis ray
  Vector forceDir;  //direction to apply force in for torque
  Vector theForce; // the force to apply
  Vector targPos; //target Position
  BigReal rayDist, height,targTheta,forceMag;
  //intermediate calc. results for readability

  //convert from degrees/timestep to radians/timestep
  int currentStep = patch->flags.step;

  //Vector  = simParams->eField;
  ////// Get list, soemthing like List = simParams->StirIndexes

  int GlobalId;
  Force *forces = r->f[Results::normal];
  BigReal energy = 0;
  Force extForce = 0.;
  Tensor extVirial;
  
  //CkPrintf("DEBUG: In ComputeStir::doForce");
  //  Loop through and check each atom
  //CkPrintf ("DEBUG: now to loop atoms, numAtoms = %d\n",numAtoms);	
    for (int i=0; i<numAtoms; i++) {
      if (molecule->is_atom_stirred(p[i].id))
      {
	//CkPrintf ("DEBUG: now to get stir params");	


	molecule->get_stir_startTheta (p[i].id);
	
	curPos = lattice.reverse_transform(p[i].position,p[i].transform);

	rayDist = distanceToRay(axisUnit,pivot,curPos);
	proj = projectionOnRay (axisUnit,pivot,curPos);
	height = findHeight (proj);
	targTheta = (omega * currentStep) + \
	             startingTheta + \
	              molecule->get_stir_startTheta (p[i].id);

	targPos = placeThetaRadius (targTheta,height,rayDist);

	//project the actual displacement ontoh a unit vector along ( (axis) cross (moment arm) ),
	//so we are applying only  a torque
      
	forceDir = (cross(axisUnit, (curPos - proj))).unit();
	forceMag = forceDir.dot(targPos - curPos); 
	
		       //unitDisp = disp.unit();
	
	theForce = (simParams->stirK) * forceMag * forceDir;
	
	forces[i]  += theForce;
	extForce += theForce;
	extVirial += outer(theForce,curPos);

			     
          //CkPrintf ("DEBUG: now to get sprintf x y z's");
	  sprintf (statbuf, "DEBUG: step= %d atom i= %d globID= %d T= %g %g %g C= %g %g %g F= %g %g %g startTheta= %g targTheta= %g forceMag= %g F.len= %g FC.len= %g\n",currentStep,i,p[i].id, targPos.x, targPos.y, targPos.z, curPos.x, curPos.y, curPos.z,theForce.x, theForce.y, theForce.z, molecule->get_stir_startTheta (p[i].id),  targTheta,  forceMag, theForce.length(), forces[i].length());

	
	  CollectionMgr::Object()->sendDataStream(statbuf);	
	  
          //////add force as Restraints does

      }

      //CkPrintf ("DEBUG: got past check for stirred\n");	
      //GlobalId = p[i].id;
      //sprintf (statbuf, "DEBUG: i= %d Global= %d\n",i,GlobalId);
      //dout << "DEBUG: i= " << i << " Global=" << GlobalId <<"\n"<<endd;
      //CollectionMgr::Object()->sendDataStream(statbuf);
    
    }

  reduction->item(REDUCTION_MISC_ENERGY) += energy;
  ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,extForce);
  ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,extVirial);
  reduction->submit();

}
/*			END OF FUNCTION force				*/
 
