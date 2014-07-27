//
//  detector.cpp
//  Xcode
//
//  Created by jacob on 7/18/11.
//
#include <iostream>
using std::cout;

#include "detector.h"
#include "photon.h"


Detector::Detector(void)
{
    detector_center.location.x = 0.0f;
    detector_center.location.y = 0.0f;
    detector_center.location.z = 0.0f;
    
    m_logger = new Logger();
}



Detector::Detector(const double x, const double y, const double z)
{
    detector_center.location.x = x;
    detector_center.location.y = y;
    detector_center.location.z = z;
    
    m_logger = new Logger();
}


Detector::Detector(const Vector3d &centerPoint)
{
    detector_center.location.x = centerPoint.location.x;
    detector_center.location.y = centerPoint.location.y;
    detector_center.location.z = centerPoint.location.z;
    
    m_logger = new Logger();
}


Detector::Detector(const boost::shared_ptr<Vector3d> centerPoint)
{
    detector_center.location.x = centerPoint->location.x;
    detector_center.location.y = centerPoint->location.y;
    detector_center.location.z = centerPoint->location.z;
    
    m_logger = new Logger();
}


Detector::~Detector()
{
    if (m_logger != NULL)
        delete m_logger;
}
