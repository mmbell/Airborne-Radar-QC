/*
 *  QCscript.h
 *
 *  Copyright 2011 Michael Bell and Cory Wolff
 *  All rights reserved.
 *
 */

#ifndef QCSCRIPT_H
#define QCSCRIPT_H

#include "AirborneRadarQC.h"

class QCscript //: public AirborneRadarQC
{

public:
     QCscript();
     ~QCscript();

     // Ruby wrappers
     bool rb_setInputPath(const std::string& in);
     bool rb_setOutputPath(const std::string& out);

private:

};

#endif
