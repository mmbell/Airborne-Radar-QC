require 'mkmf-rice'
$DEFINES  = " -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED"
$CFLAGS   += " -pipe -g -gdwarf-2 -Wall -W #$DEFINES"
$CXXFLAGS += " -pipe -g -gdwarf-2 -Wall -W #$DEFINES"
$INCFLAGS += " -I/usr/local/Qt4.7/mkspecs/macx-g++ -I. -I/Library/Frameworks/QtCore.framework/Versions/4/Headers -I/usr/include/QtCore -I/Library/Frameworks/QtGui.framework/Versions/4/Headers -I/usr/include/QtGui -I/usr/include -I. -I. -F/Library/Frameworks"
$DLDFLAGS += " -headerpad_max_install_names"
$LIBS     += " -F/Library/Frameworks -L/Library/Frameworks -framework QtGui -framework QtCore"

create_makefile('QCscript')
