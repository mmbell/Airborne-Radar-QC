TEMPLATE = app
TARGET = 
DEPENDPATH += .
INCLUDEPATH += .

# Input
HEADERS += AirborneRadarQC.h \
           Dorade.h \
           precision.h \
           read_dorade.h \
           RecursiveFilter.h
SOURCES += AirborneRadarQC.cpp \
	   Dorade.cpp \
	   RecursiveFilter.cpp
LIBS += -lgeotiff -ltiff -lgeographic

