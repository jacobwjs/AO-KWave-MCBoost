//
//  logger.cpp
//  Xcode
//
//  Created by jacob on 7/18/11.
//  Copyright 2011 BMPI. All rights reserved.
//

#include "multikey.h"
#include "vector3D.h"
#include "photon.h"
#include "logger.h"
#include <cmath>
using std::cos;







Logger::Logger()
:LoggerBase()
{
   
}


void Logger::Destroy()
{
    this->~Logger();
}

Logger::~Logger()
{
    //LoggerBase::Destroy();
}





