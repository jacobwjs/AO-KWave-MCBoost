//
//  LoggerBase.cpp
//  Xcode
//
//  Created by Jacob Staley on 7/18/14.
//  Copyright 2014 BMPI, University of Twente. All rights reserved.
//

#include "multikey.h"
#include "vector3D.h"
#include "photon.h"
#include "loggerSingleton.h"
#include <cmath>
using std::cos;


LoggerSingleton * LoggerSingleton::pInstance = NULL;




LoggerSingleton * LoggerSingleton::getInstance()
{
    if (!pInstance)
    {
        pInstance = new LoggerSingleton;
    }
    
    return pInstance;
}



void LoggerSingleton::Destroy()
{
    this->~LoggerSingleton();
}

LoggerSingleton::~LoggerSingleton()
{
    
}







