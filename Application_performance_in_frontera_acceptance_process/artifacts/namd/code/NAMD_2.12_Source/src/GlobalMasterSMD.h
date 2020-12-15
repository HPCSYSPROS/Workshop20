/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef GLOBALMASTERSMD_H
#define GLOBALMASTERSMD_H

class GlobalMasterSMD : public GlobalMaster {
public:
  /* initializes this to be a SMD master maintaing the group specified
     in <filename> bound by a spring with k=<spring_constant> to a
     point moving from its initial CM in the direction <direction>
     with speed <velocity>.  Must also be passed the desired rate of
     output <output_frequency> and the inital timestep # of the
     simulation <first_timestep> */
  GlobalMasterSMD(BigReal spring_constant,
		  BigReal transverse_spring_constant,
		  BigReal velocity,
		  const Vector direction, int output_frequency,
		  int first_timestep, const char *filename, int);
  ~GlobalMasterSMD();

private:
  virtual void calculate();

  /* prints out a message for timestep <t> reporting that the forced
     atoms are at CM position <p> and are feeling force <f> */
  void output(int t, Position p, Force f);

  /* read the atoms from the file into group 0 (only call this once!) */
  void parseAtoms(const char *file, int);
 
  BigReal k;
  BigReal k2;
  BigReal moveVel;   // A/timestep
  Vector moveDir;
  int outputFreq;
  Position cm;       // Initial center of mass
  int currentTime;   // Keep track of elapsed time steps for yourself!
};
#endif

