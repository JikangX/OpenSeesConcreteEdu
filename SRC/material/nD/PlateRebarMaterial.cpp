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
// $Date: 2012-05-26 22:03:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlateRebarMaterial.cpp,v $

//
// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//
// Generic Plate Rebar Material
//
/* Ref: Lu X, Lu XZ, Guan H, Ye LP, Collapse simulation of reinforced 
concrete high-rise building induced by extreme earthquakes, 
Earthquake Engineering & Structural Dynamics, 2013, 42(5): 705-723*/

#include <PlateRebarMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>   //Antonios Vytiniotis used for the recorder
#include <math.h>
#include <elementAPI.h>
#include <cmath>
#include <cstring>

//null constructor
PlateRebarMaterial::PlateRebarMaterial( ) :
NDMaterial(0, ND_TAG_PlateRebarMaterial ),
theMat(nullptr), angle(0.0), c(0.0), s(0.0), strain(5),
stress(5), tangent(5,5)
{ }

void* OPS_PlateRebarMaterial()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 3) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: nDMaterial PlateRebar tag? matTag? angle?" << endln;
	return 0;
    }

    int tag[2];
    numdata = 2;
    if (OPS_GetIntInput(&numdata,tag)<0) {
	opserr << "WARNING invalid nDMaterial PlateRebar tag or matTag" << endln;
	return 0;
    }

    UniaxialMaterial *theMaterial = OPS_getUniaxialMaterial(tag[1]);
    if (theMaterial == 0) {
	opserr << "WARNING uniaxialmaterial does not exist\n";
	opserr << "UniaxialMaterial: " << tag[1];
	opserr << "\nPlateRebar nDMaterial: " << tag[0] << endln;
	return 0;
    }

    double angle;
    numdata = 1;
    if (OPS_GetDoubleInput(&numdata,&angle)<0) {
	opserr << "WARNING invalid angle" << endln;
	return 0;
    }
      
    NDMaterial* mat = new PlateRebarMaterial( tag[0], *theMaterial, angle);

    if (mat == 0) {
	opserr << "WARNING: failed to create PlateRebar material\n";
	return 0;
    }

    return mat;
}


//full constructor
PlateRebarMaterial::PlateRebarMaterial(int tag,
                                       UniaxialMaterial &uniMat,
                                       double ang) :
NDMaterial( tag, ND_TAG_PlateRebarMaterial ),
theMat(nullptr), angle(ang), c(0.0), s(0.0), strain(5),
stress(5), tangent(5,5)
{
  theMat = uniMat.getCopy() ;
  double rang = ang * 4.0 * asin(1.0)/360.0;
  c = cos(rang);
  s = sin(rang);
}


double
PlateRebarMaterial::axialStrain( const Vector &plateStrain ) const
{
  if (angle == 0)
    return plateStrain(0);
  else if (angle == 90)
    return plateStrain(1);
  return plateStrain(0) * c * c + plateStrain(1) * s * s + plateStrain(2) * c * s;
}


//destructor
PlateRebarMaterial::~PlateRebarMaterial( ) 
{ 
  if (theMat != 0) delete theMat ;
} 


//make a clone of this material
NDMaterial*
PlateRebarMaterial::getCopy( ) 
{
  // copy the uniaxial material together with its state
  PlateRebarMaterial *clone = new PlateRebarMaterial( this->getTag(),
                                                      *theMat,
                                                      angle ) ;
  clone->strain = strain;
  for (int i = 0; i < 5; i++)
    clone->committedStrain[i] = committedStrain[i];
  const double axial = this->axialStrain(strain);
  if (clone->theMat->getStrain() != axial)
    clone->theMat->setTrialStrain(axial);

  return clone ;
}


//make a clone of this material
NDMaterial* 
PlateRebarMaterial::getCopy( const char *type ) 
{
  if (strcmp(type,this->getType()) == 0)
    return this->getCopy( ) ;
  else
    return 0;
}


//send back order of strain in vector form
int 
PlateRebarMaterial::getOrder( ) const
{
  return 5 ;
}


const char*
PlateRebarMaterial::getType( ) const 
{
  return "PlateFiber" ; 
}



//swap history variables
int 
PlateRebarMaterial::commitState( ) 
{
  const int result = theMat->commitState( ) ;
  if (result == 0)
    for (int i = 0; i < 5; i++)
      committedStrain[i] = strain(i);
  return result ;
}



//revert to last saved state
int
PlateRebarMaterial::revertToLastCommit( )
{
  for (int i = 0; i < 5; i++)
    strain(i) = committedStrain[i];
  int result = theMat->revertToLastCommit( ) ;
  const double axial = this->axialStrain(strain);
  if (result == 0 && theMat->getStrain() != axial)
    result = theMat->setTrialStrain(axial);
  return result ;
}


//revert to start
int
PlateRebarMaterial::revertToStart( )
{
  for (int i = 0; i < 5; i++)
    committedStrain[i] = 0.0;
  strain.Zero();
  return theMat->revertToStart( ) ;
}


//mass per unit volume
double
PlateRebarMaterial::getRho( )
{
  return theMat->getRho( ) ;
}


//receive the strain
int 
PlateRebarMaterial::setTrialStrain( const Vector &strainFromElement )
{
  if (strainFromElement.Size() != 5)
    return -1;
  for (int i = 0; i < 5; i++)
    if (!std::isfinite(strainFromElement(i)))
      return -1;

  // the wrapper strain changes only if the uniaxial material accepts the trial
  const int result = theMat->setTrialStrain(this->axialStrain(strainFromElement));
  if (result == 0)
    strain = strainFromElement;
  return result;
}


int
PlateRebarMaterial::setTrialStrain( const Vector &strainFromElement,
                                    const Vector &rate )
{
  if (strainFromElement.Size() != 5 || rate.Size() != 5)
    return -1;
  for (int i = 0; i < 5; i++)
    if (!std::isfinite(strainFromElement(i)) || !std::isfinite(rate(i)))
      return -1;

  const int result = theMat->setTrialStrain(this->axialStrain(strainFromElement),
                                            this->axialStrain(rate));
  if (result == 0)
    strain = strainFromElement;
  return result;
}


int
PlateRebarMaterial::setTrialStrainIncr( const Vector &strainIncrement )
{
  if (strainIncrement.Size() != 5)
    return -1;
  const Vector next = strain + strainIncrement;
  return this->setTrialStrain(next);
}


int
PlateRebarMaterial::setTrialStrainIncr( const Vector &strainIncrement,
                                        const Vector &rate )
{
  if (strainIncrement.Size() != 5)
    return -1;
  const Vector next = strain + strainIncrement;
  return this->setTrialStrain(next, rate);
}


//send back the strain
const Vector& 
PlateRebarMaterial::getStrain( )
{
  return strain ;
}


//send back the stress 
const Vector&  
PlateRebarMaterial::getStress( )
{
  double sig = theMat->getStress();

  stress.Zero();
  if (angle == 0) 
    stress(0) = sig; 
  else if (angle == 90)
    stress(1)= sig;
  else {
    stress(0) = sig * c * c;
    stress(1) = sig * s * s;
    stress(2) = sig * c * s;
  }

  return stress ;
}


//send back the tangent 
const Matrix&  
PlateRebarMaterial::getTangent( )
{
  double tan = theMat->getTangent( ) ;

  tangent.Zero();
  if (angle == 0) 
    tangent(0,0) = tan;
  else if (angle == 90)
    tangent(1,1) = tan;
  else {
    tangent(0,0) = tan * c * c * c * c ;
    tangent(0,1) = tan * c * c * s * s ;
    tangent(0,2) = tan * c * c * c * s ;
    tangent(1,0) = tangent(0,1) ;
    tangent(1,1) = tan * s * s * s * s ;
    tangent(1,2) = tan * c * s * s * s ;
    tangent(2,0) = tangent(0,2) ;
    tangent(2,1) = tangent(1,2) ;
    tangent(2,2) = tangent(0,1) ;
  }

  return tangent ;
}

const Matrix&  
PlateRebarMaterial::getInitialTangent
( )
{
  double tan = theMat->getInitialTangent( ) ;
  tangent.Zero();
  if (angle == 0) 
    tangent(0,0) = tan;
  else if (angle == 90)
    tangent(1,1) = tan;
  else {
    tangent(0,0) = tan * c * c * c * c ;
    tangent(0,1) = tan * c * c * s * s ;
    tangent(0,2) = tan * c * c * c * s ;
    tangent(1,0) = tangent(0,1) ;
    tangent(1,1) = tan * s * s * s * s ;
    tangent(1,2) = tan * c * s * s * s ;
    tangent(2,0) = tangent(0,2) ;
    tangent(2,1) = tangent(1,2) ;
    tangent(2,2) = tangent(0,1) ;
  }

  return tangent ;
}


//print out data
void  
PlateRebarMaterial::Print( OPS_Stream &s, int flag )
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "PlateRebar Material tag: " << this->getTag() << endln;
        s << "angle: " << angle << endln;
        s << "using uniaxial material: " << endln;
        theMat->Print(s, flag);
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"PlateRebarMaterial\", ";
        s << "\"angle\": " << angle << ", ";
        s << "\"material\": \"" << theMat->getTag() << "\"}";
    }
}


int 
PlateRebarMaterial::sendSelf(int commitTag, Channel &theChannel) 
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
    opserr << "PlateRebarMaterial::sendSelf() - failed to send ID" << endln;
    return res;
  }

  // angle, trial strain(5), committed strain(5)
  Vector vecData(11);
  vecData(0) = angle;
  for (int i = 0; i < 5; i++) {
    vecData(i + 1) = strain(i);
    vecData(i + 6) = committedStrain[i];
  }

  res = theChannel.sendVector(dataTag, commitTag, vecData);
  if (res < 0) {
    opserr << "PlateRebarMaterial::sendSelf() - failed to send Vector" << endln;
    return res;
  }

  // now send the materials data
  res += theMat->sendSelf(commitTag, theChannel);
  if (res < 0) 
    opserr << "PlateRebarMaterial::sendSelf() - failed to send material" << endln;

  return res;
}

int 
PlateRebarMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  int dataTag = this->getDbTag();

  // recv an id containing the tag and associated materials class and db tags
  static ID idData(4);
  res = theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "PlateRebarMaterial::sendSelf() - failed to receive ID" << endln;
    return res;
  }

  this->setTag(idData(0));
  // negative class tag: packet with committed strain (11 values);
  // positive: earlier packet with trial strain only (6 values)
  const bool hasHistory = idData(1) < 0;
  int matClassTag = hasHistory ? -idData(1) : idData(1);
  if (theMat == nullptr || theMat->getClassTag() != matClassTag) {
    if (theMat != 0) delete theMat;
    theMat = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMat == 0) {
      opserr << "PlateRebarMaterial::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMat->setDbTag(idData(2));

  Vector vecData(hasHistory ? 11 : 6);
  res = theChannel.recvVector(dataTag, commitTag, vecData);
  if (res < 0) {
    opserr << "PlateRebarMaterial::sendSelf() - failed to receive Vector" << endln;
    return res;
  }
  angle = vecData(0);
  for (int i = 0; i < 5; i++) {
    strain(i) = vecData(i + 1);
    committedStrain[i] = hasHistory ? vecData(i + 6) : strain(i);
  }
  // same conversion as the constructor (the former truncated factor
  // 0.0174532925 changed non-axis rebar results after a restart)
  double rang = angle * 4.0 * asin(1.0)/360.0;
  c = cos(rang);
  s = sin(rang);

  // now receive the materials data
  res = theMat->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0)
    opserr << "PlateRebarMaterial::sendSelf() - failed to receive material1" << endln;
  else {
    const double axial = this->axialStrain(strain);
    if (theMat->getStrain() != axial)
      res = theMat->setTrialStrain(axial);
  }

  return res;
}

int PlateRebarMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
    return theMat->setParameter(argv, argc, param);
}
 
