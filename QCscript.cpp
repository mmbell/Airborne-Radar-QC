/*
 *  QCscript.cpp
 *
 *  Copyright 2011 Michael Bell and Cory Wolff
 *  All rights reserved.
 *
 */

#include "rice/Data_Type.hpp"
#include "rice/Constructor.hpp"
#include "QCscript.h"
#include "AirborneRadarQC.h"
#include "QString.h"

using namespace Rice;

extern "C"
void Init_QCscript(void)
{
  RUBY_TRY
  {

    define_class<AirborneRadarQC>("AirborneRadarQC");
//      .define_constructor(Constructor<AirborneRadarQC>());

    define_class<QCscript, AirborneRadarQC>("QCscript")
      .define_constructor(Constructor<QCscript>())
      .define_method("processSweeps", &QCscript::processSweeps)
      .define_method("setInputPath", &QCscript::rb_setInputPath)
      .define_method("setOutputPath", &QCscript::rb_setOutputPath);
  }
  RUBY_CATCH
} 

QCscript::QCscript()
   : AirborneRadarQC()
{

}

QCscript::~QCscript()
{
}


bool QCscript::rb_setInputPath(const std::string& in)
{
	QString path = QString::fromStdString(in);
	return setInputPath(path);
}

bool QCscript::rb_setOutputPath(const std::string& out)
{
        QString path = QString::fromStdString(out);
        return setOutputPath(path);
}

