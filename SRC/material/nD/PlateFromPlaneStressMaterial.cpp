/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.0 $
// $Date: 2012-05-28 22:03:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlateFromPlaneStressMaterial.cpp,v $

//
// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//
// Generic Plate Material from Plane Stress Material
//
/* Ref: Lu X, Lu XZ, Guan H, Ye LP, Collapse simulation of reinforced 
concrete high-rise building induced by extreme earthquakes, 
Earthquake Engineering & Structural Dynamics, 2013, 42(5): 705-723*/


#include <PlateFromPlaneStressMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>
#include <elementAPI.h>
#include <cmath>
#include <cstring>

void* OPS_PlateFromPlaneStressMaterial()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 3) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: nDMaterial PlateFromPlaneStress tag? matTag? gmod?" << endln;
	return 0;
    }

    int tag[2];
    numdata = 2;
    if (OPS_GetIntInput(&numdata,tag)<0) {
	opserr << "WARNING invalid nDMaterial PlateFromPlaneStress tag and matTag" << endln;
	return 0;
    }

    NDMaterial *threeDMaterial = OPS_getNDMaterial(tag[1]);
    if (threeDMaterial == 0) {
	opserr << "WARNING nD material does not exist\n";
	opserr << "nD material: " << tag[1];
	opserr << "\nPlateFromplanestress nDMaterial: " << tag[0] << endln;
	return 0;
    }

    double gmod;
    numdata = 1;
    if (OPS_GetDoubleInput(&numdata,&gmod)<0) {
	opserr << "WARNING invalid gmod" << endln;
	return 0;
    }
      
    NDMaterial* mat = new PlateFromPlaneStressMaterial( tag[0], *threeDMaterial, gmod);

    if (mat == 0) {
	opserr << "WARNING: failed to create PlateFromplanestress material\n";
	return 0;
    }

    return mat;
}

//null constructor
PlateFromPlaneStressMaterial::PlateFromPlaneStressMaterial( ) : 
NDMaterial(0, ND_TAG_PlateFromPlaneStressMaterial ), 
theMat(nullptr), gmod(0.0), strain(5), planeStressStrain(3),
planeStressRate(3), stress(5), tangent(5,5)
{ }


//full constructor
PlateFromPlaneStressMaterial::PlateFromPlaneStressMaterial(
				   int tag, NDMaterial &ndMat, double g ) :
NDMaterial( tag, ND_TAG_PlateFromPlaneStressMaterial ),
theMat(nullptr), gmod(g), strain(5), planeStressStrain(3),
planeStressRate(3), stress(5), tangent(5,5)
{
  // a material that is already plane stress is copied as is
  const char *materialType = ndMat.getType();
  if (strcmp(materialType, "PlaneStress") == 0 ||
      strcmp(materialType, "PlaneStress2D") == 0)
    theMat = ndMat.getCopy() ;
  else
    theMat = ndMat.getCopy("PlaneStress") ;
}


// re-synchronise theMat with the in-plane part of strain when they differ
// (after revert, copy or restart of a wrapper whose trial was not committed)
static int
syncPlaneStressStrain(NDMaterial *theMat, const Vector &strain, Vector &planeStressStrain)
{
  const Vector &current = theMat->getStrain();
  bool consistent = current.Size() == 3;
  for (int i = 0; i < 3 && consistent; i++)
    consistent = current(i) == strain(i);
  if (consistent)
    return 0;
  for (int i = 0; i < 3; i++)
    planeStressStrain(i) = strain(i);
  return theMat->setTrialStrain(planeStressStrain);
}


//destructor
PlateFromPlaneStressMaterial::~PlateFromPlaneStressMaterial( ) 
{ 
  if (theMat != 0) delete theMat ;
} 



//make a clone of this material
NDMaterial*
PlateFromPlaneStressMaterial::getCopy( ) 
{
  // copy the wrapped material together with its state
  PlateFromPlaneStressMaterial *clone = new PlateFromPlaneStressMaterial();
  clone->setTag(this->getTag());
  clone->theMat = theMat->getCopy();
  if (clone->theMat == nullptr) {
    delete clone;
    return nullptr;
  }
  clone->gmod = gmod;
  clone->strain = strain;
  for (int i = 0; i < 5; i++)
    clone->committedStrain[i] = committedStrain[i];
  syncPlaneStressStrain(clone->theMat, clone->strain, clone->planeStressStrain);

  return clone ;
}


//make a clone of this material
NDMaterial* 
PlateFromPlaneStressMaterial::getCopy( const char *type ) 
{
  if (strcmp(type, this->getType()) == 0)
    return this->getCopy( ) ;
  else
    return 0;
}


//send back order of strain in vector form
int 
PlateFromPlaneStressMaterial::getOrder( ) const
{
  return 5 ;
}


const char*
PlateFromPlaneStressMaterial::getType( ) const 
{
  return "PlateFiber" ; 
}



//swap history variables
int 
PlateFromPlaneStressMaterial::commitState( ) 
{
  const int result = theMat->commitState( ) ;
  if (result == 0)
    for (int i = 0; i < 5; i++)
      committedStrain[i] = strain(i);
  return result ;
}



//revert to last saved state
int
PlateFromPlaneStressMaterial::revertToLastCommit( )
{
  for (int i = 0; i < 5; i++)
    strain(i) = committedStrain[i];
  int result = theMat->revertToLastCommit( ) ;
  if (result == 0)
    result = syncPlaneStressStrain(theMat, strain, planeStressStrain);
  return result ;
}


//revert to start
int
PlateFromPlaneStressMaterial::revertToStart( )
{
  for (int i = 0; i < 5; i++)
    committedStrain[i] = 0.0;

  strain.Zero();

  return theMat->revertToStart( ) ;
}


//mass per unit volume
double
PlateFromPlaneStressMaterial::getRho( )
{
  return theMat->getRho( ) ;
}


//receive the strain
int 
PlateFromPlaneStressMaterial::setTrialStrain( const Vector &strainFromElement )
{
  if (strainFromElement.Size() != 5)
    return -1;
  for (int i = 0; i < 5; i++)
    if (!std::isfinite(strainFromElement(i)))
      return -1;

  for (int i = 0; i < 3; i++)
    planeStressStrain(i) = strainFromElement(i);

  // the wrapper strain changes only if the wrapped material accepts the trial
  const int result = theMat->setTrialStrain(planeStressStrain);
  if (result == 0)
    strain = strainFromElement;
  return result;
}


int
PlateFromPlaneStressMaterial::setTrialStrain( const Vector &strainFromElement,
                                              const Vector &rate )
{
  if (strainFromElement.Size() != 5 || rate.Size() != 5)
    return -1;
  for (int i = 0; i < 5; i++)
    if (!std::isfinite(strainFromElement(i)) || !std::isfinite(rate(i)))
      return -1;

  for (int i = 0; i < 3; i++) {
    planeStressStrain(i) = strainFromElement(i);
    planeStressRate(i) = rate(i);
  }

  const int result = theMat->setTrialStrain(planeStressStrain, planeStressRate);
  if (result == 0)
    strain = strainFromElement;
  return result;
}


int
PlateFromPlaneStressMaterial::setTrialStrainIncr( const Vector &strainIncrement )
{
  if (strainIncrement.Size() != 5)
    return -1;
  const Vector next = strain + strainIncrement;
  return this->setTrialStrain(next);
}


int
PlateFromPlaneStressMaterial::setTrialStrainIncr( const Vector &strainIncrement,
                                                  const Vector &rate )
{
  if (strainIncrement.Size() != 5)
    return -1;
  const Vector next = strain + strainIncrement;
  return this->setTrialStrain(next, rate);
}


//send back the strain
const Vector& 
PlateFromPlaneStressMaterial::getStrain( )
{
  return strain ;
}


//send back the stress 
const Vector&  
PlateFromPlaneStressMaterial::getStress( )
{
  //three dimensional stress
  const Vector &PSStress = theMat->getStress();

  stress(0) = PSStress(0);
  stress(1) = PSStress(1);
  stress(2) = PSStress(2);

  stress(3) = gmod * strain(3);
  stress(4) = gmod * strain(4);

  return stress ;
}


//send back the tangent 
const Matrix&  
PlateFromPlaneStressMaterial::getTangent( )
{
  const Matrix &PSTangent = theMat->getTangent( ) ;

  tangent.Zero();
  
  tangent(0,0) = PSTangent(0,0);
  tangent(0,1) = PSTangent(0,1);
  tangent(0,2) = PSTangent(0,2);
  tangent(1,0) = PSTangent(1,0);
  tangent(1,1) = PSTangent(1,1);
  tangent(1,2) = PSTangent(1,2);
  tangent(2,0) = PSTangent(2,0);
  tangent(2,1) = PSTangent(2,1);
  tangent(2,2) = PSTangent(2,2);

  tangent(3,3) = gmod;
  tangent(4,4) = gmod;

  return tangent ;
}

//send back the tangent 
const Matrix&  
PlateFromPlaneStressMaterial::getInitialTangent
( )
{
  const Matrix &PSTangent = theMat->getInitialTangent( ) ;

  tangent.Zero();
  
  tangent(0,0) = PSTangent(0,0);
  tangent(0,1) = PSTangent(0,1);
  tangent(0,2) = PSTangent(0,2);
  tangent(1,0) = PSTangent(1,0);
  tangent(1,1) = PSTangent(1,1);
  tangent(1,2) = PSTangent(1,2);
  tangent(2,0) = PSTangent(2,0);
  tangent(2,1) = PSTangent(2,1);
  tangent(2,2) = PSTangent(2,2);
  tangent(3,3) = gmod;
  tangent(4,4) = gmod;

  return tangent ;
}


//print out data
void  
PlateFromPlaneStressMaterial::Print( OPS_Stream &s, int flag )
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "PlateFromPlaneStress Material tag: " << this->getTag() << "" << endln;
        s << "G: " << gmod << endln;
        s << "using PlaneStress material: " << endln;
        theMat->Print(s, flag);
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"PlateFromPlaneStressMaterial\", ";
        s << "\"G\": " << gmod << ", ";
        s << "\"material\": \"" << theMat->getTag() << "\"}";
    }
}


int 
PlateFromPlaneStressMaterial::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0;

  int dataTag = this->getDbTag();

  int matDbTag;
  
  // idData(1) < 0 marks a packet that carries the wrapper strain history;
  // idData(3) is reserved (0) to keep the layout of OpenSeesCUDA checkpoints
  static ID idData(4);
  idData(0) = this->getTag();
  idData(1) = -theMat->getClassTag();
  matDbTag = theMat->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMat->setDbTag(matDbTag);
  }
  idData(2) = matDbTag;
  idData(3) = 0;

  res = theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "PlateFromPlaneStressMaterial::sendSelf() - failed to send data" << endln;
    return res;
  }

  // gmod, trial strain(5), committed strain(5)
  Vector vecData(11);
  vecData(0) = gmod;
  for (int i = 0; i < 5; i++) {
    vecData(i + 1) = strain(i);
    vecData(i + 6) = committedStrain[i];
  }

  res = theChannel.sendVector(dataTag, commitTag, vecData);
  if (res < 0) {
    opserr << "PlateFromPlaneStressMaterial::sendSelf() - failed to send data" << endln;
    return res;
  }

  // now send the materials data
  res += theMat->sendSelf(commitTag, theChannel);
  if (res < 0) 
    opserr << "PlateFromPlaneStressMaterial::sendSelf() - failed to send material1" << endln;

  return res;
}

int 
PlateFromPlaneStressMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  int dataTag = this->getDbTag();

  // recv an id containing the tag and associated materials class and db tags
  static ID idData(4);
  res = theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "PlateFromPlaneStressMaterial::sendSelf() - failed to receive id data" << endln;
    return res;
  }

  this->setTag(idData(0));
  // negative class tag: packet with committed strain (11 values);
  // positive: earlier packet with trial strain only (6 values)
  const bool hasHistory = idData(1) < 0;
  int matClassTag = hasHistory ? -idData(1) : idData(1);
  if (theMat == nullptr || theMat->getClassTag() != matClassTag) {
    if (theMat != 0) delete theMat;
    theMat = theBroker.getNewNDMaterial(matClassTag);
    if (theMat == 0) {
      opserr << "PlateFromPlaneStressMaterial::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMat->setDbTag(idData(2));

  Vector vecData(hasHistory ? 11 : 6);
  res = theChannel.recvVector(dataTag, commitTag, vecData);
  if (res < 0) {
    opserr << "PlateFromPlaneStressMaterial::sendSelf() - failed to receive vector data" << endln;
    return res;
  }
  gmod = vecData(0);
  for (int i = 0; i < 5; i++) {
    strain(i) = vecData(i + 1);
    committedStrain[i] = hasHistory ? vecData(i + 6) : strain(i);
  }

  // now receive the materials data
  res = theMat->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0)
    opserr << "PlateFromPlaneStressMaterial::sendSelf() - failed to receive material1" << endln;
  else
    res = syncPlaneStressStrain(theMat, strain, planeStressStrain);

  return res;
}

//setResponse - added by V.K. Papanikolaou [AUTh] - start
Response*
PlateFromPlaneStressMaterial::setResponse(const char** argv, int argc, OPS_Stream& output)
{
    if (strcmp(argv[0], "Tangent") == 0 || strcmp(argv[0], "tangent") == 0 ||
        strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0 ||
        strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0) {
        
        return NDMaterial::setResponse(argv, argc, output);  // for stresses/strains, get response from NDMaterial
    }

    // for material-specific output

    Response *theResponse = 0;
    theResponse = theMat->setResponse(argv, argc, output);

    if (theResponse == 0) 
        return NDMaterial::setResponse(argv, argc, output);  // not implemented, get default (zero) response from NDMaterial

    return theResponse;  // implemented, get damage response from theMat
}

//setResponse - added by V.K. Papanikolaou [AUTh] - end

int PlateFromPlaneStressMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
    return theMat->setParameter(argv, argc, param);
}
