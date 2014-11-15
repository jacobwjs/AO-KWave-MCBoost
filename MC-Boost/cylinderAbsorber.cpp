//
//  CylinderAbsorber.cpp
//  Xcode
//
//  Created by jacob on 6/20/11.
//
#include "absorber.h"
#include "cylinderAbsorber.h"


// Cylinder constructor takes the radius and the top and bottom coordinates
// of the cylinder.
CylinderAbsorber::CylinderAbsorber(const double radius, const boost::shared_ptr<Vector3d> cap_A,
                                                        const boost::shared_ptr<Vector3d> cap_B)
:Absorber(center)
{
    // Calculate the length of the cylinder.
    this->length = sqrt(pow(cap_A->location.x - cap_B->location.x, 2) +
                        pow(cap_A->location.y - cap_B->location.y, 2) +
                        pow(cap_A->location.z - cap_B->location.z, 2));
    
    this->radius = radius;
    this->cap_A  = (*cap_A);
    this->cap_B  = (*cap_B);
}


CylinderAbsorber::~CylinderAbsorber()
{
    // STUB
}


bool CylinderAbsorber::crossedAbsorber(const boost::shared_ptr<Vector3d> photonVector)
{
    // STUB
}


bool CylinderAbsorber::hitAbsorberBoundary(const boost::shared_ptr<Vector3d> photonVector)
{
    // STUB
}

// Need to draw a line from top and bottom disc of the cylinder.
// If the length of the line from the photon to that line is 
// larger than the cylinder's radius, then return false.  Else
// return true.
// Also need to check if the photon had passed through the cylinder.
bool CylinderAbsorber::inAbsorber(const boost::shared_ptr<Vector3d> photonVector)
{
    if (inCylinderVolume(photonVector))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool CylinderAbsorber::inCylinderVolume(const boost::shared_ptr<Vector3d> photonLocation)
{
    // Subtract the previous location (B) of the photon from the current location (A)
    // to form a new vector.
    boost::shared_ptr<Vector3d> AB = (*cap_A) - (*cap_B);
    
    // Subtract the current location of the photon (A) from the center location of the
    // absorber (center) to yield a new vector.
    boost::shared_ptr<Vector3d> AP = (*cap_A) - (*photonLocation);
    
    // Take the cross-product of AB and AC to get the area of the parallelogram formed
    // from A-B and A-C
    double c1 = VectorMath::dotProduct(AB, AP);
    double c2 = VectorMath::dotProduct(AB, AB);
    
    double t = c1 / c2;
    
    /// 't' is the portion along the line from 'cap_A' to 'cap_B' that represents the point along
    /// that line that is the minimum distance from the line to point 'P'. If 't' is less than zero
    /// or greater than 1 we know we don't lie between the caps of the cylinder and we're done. If 't'
    /// is between 0 and 1 we lie along the line, then we need to check the distance from the point
    /// along the line, i.e. L(t), to the location of the photon ('P'). If this distance is less than
    /// the radius we are in the volume of the cylinder and we should return true.
    if ((t < 0.0) || (t > 1.0))
    {
        return false;
    }
    
    /// If we make it here the photon is between the two caps of the cylinder somewhere in 3D space,
    /// now we need to check if it is within the radius of the cylinder.
    
    /// Point on the line between the caps of the cylinder that is closest to the photon.
    Pt = (*cap_A) + t*((*cap_B)-(*cap_A));
    
    double distance_to_photon = VectorMath::Distance(abs((*Pt) - (*photonLocation)));
    
    if (distance_to_photon < radius)
    {
        return true;
    }
    else
    {
        return false;
    }
}


void CylinderAbsorber::cartesianToCylindrical(void)
{
    // STUB
}



