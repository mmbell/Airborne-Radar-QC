
SHELL = /bin/sh

#### Start of system configuration section. ####

srcdir = .
topdir = /Users/mbell/.rvm/rubies/ruby-1.9.2-head/include/ruby-1.9.1
hdrdir = /Users/mbell/.rvm/rubies/ruby-1.9.2-head/include/ruby-1.9.1
arch_hdrdir = /Users/mbell/.rvm/rubies/ruby-1.9.2-head/include/ruby-1.9.1/$(arch)
VPATH = $(srcdir):$(arch_hdrdir)/ruby:$(hdrdir)/ruby

prefix = $(DESTDIR)/Users/mbell/.rvm/rubies/ruby-1.9.2-head

exec_prefix = $(prefix)

rubylibprefix = $(libdir)/$(RUBY_BASE_NAME)

bindir = $(exec_prefix)/bin

sbindir = $(exec_prefix)/sbin

libexecdir = $(exec_prefix)/libexec

datarootdir = $(prefix)/share

datadir = $(datarootdir)

sysconfdir = $(prefix)/etc

sharedstatedir = $(prefix)/com

localstatedir = $(prefix)/var

includedir = $(prefix)/include

oldincludedir = $(DESTDIR)/usr/include

docdir = $(datarootdir)/doc/$(PACKAGE)

infodir = $(datarootdir)/info

htmldir = $(docdir)

dvidir = $(docdir)

pdfdir = $(docdir)

psdir = $(docdir)

libdir = $(exec_prefix)/lib

localedir = $(datarootdir)/locale

mandir = $(datarootdir)/man

ridir = $(datarootdir)/$(RI_BASE_NAME)

sitedir = $(rubylibprefix)/site_ruby

vendordir = $(rubylibprefix)/vendor_ruby

rubyhdrdir = $(includedir)/$(RUBY_BASE_NAME)-$(ruby_version)

sitehdrdir = $(rubyhdrdir)/site_ruby

vendorhdrdir = $(rubyhdrdir)/vendor_ruby

rubylibdir = $(rubylibprefix)/$(ruby_version)

archdir = $(rubylibdir)/$(arch)

sitelibdir = $(sitedir)/$(ruby_version)

sitearchdir = $(sitelibdir)/$(sitearch)

vendorlibdir = $(vendordir)/$(ruby_version)

vendorarchdir = $(vendorlibdir)/$(sitearch)


CC = gcc
CXX = g++
LIBRUBY = $(LIBRUBY_SO)
LIBRUBY_A = lib$(RUBY_SO_NAME)-static.a
LIBRUBYARG_SHARED = -l$(RUBY_SO_NAME)
LIBRUBYARG_STATIC = -lruby.1.9.1-static
OUTFLAG = -o 
COUTFLAG = -o 

RUBY_EXTCONF_H = 
cflags   =  $(optflags) $(debugflags) $(warnflags)
optflags = -O3
debugflags = -ggdb
warnflags = -Wextra -Wno-unused-parameter -Wno-parentheses -Wpointer-arith -Wwrite-strings -Wno-missing-field-initializers -Wshorten-64-to-32 -Wno-long-long
CFLAGS   = -fno-common $(cflags)  -fno-common -pipe -pipe -g -gdwarf-2 -Wall -W  -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED 
INCFLAGS = -I. -I$(arch_hdrdir) -I$(hdrdir)/ruby/backward -I$(hdrdir) -I$(srcdir) -I/usr/local/Qt4.7/mkspecs/macx-g++ -I. -I/Library/Frameworks/QtCore.framework/Versions/4/Headers -I/usr/include/QtCore -I/Library/Frameworks/QtGui.framework/Versions/4/Headers -I/usr/include/QtGui -I/usr/include -I. -I. -F/Library/Frameworks
DEFS     = 
CPPFLAGS =  -D_XOPEN_SOURCE -D_DARWIN_C_SOURCE $(DEFS) $(cppflags) -arch x86_64 -I/Users/mbell/.rvm/gems/ruby-1.9.2-head/gems/rice-1.4.2/ruby/lib/include
CXXFLAGS = $(CFLAGS)  -Wall -g -pipe -g -gdwarf-2 -Wall -W  -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED
ldflags  = -L.   -L/Users/mbell/.rvm/gems/ruby-1.9.2-head/gems/rice-1.4.2/ruby/lib/lib
dldflags = -Wl,-undefined,dynamic_lookup -Wl,-multiply_defined,suppress -Wl,-flat_namespace -headerpad_max_install_names
ARCH_FLAG = 
DLDFLAGS = $(ldflags) $(dldflags)
LDSHARED = g++ -dynamic -bundle
LDSHAREDXX = $(CXX) -dynamic -bundle
AR = ar
EXEEXT = 

RUBY_BASE_NAME = ruby
RUBY_INSTALL_NAME = ruby
RUBY_SO_NAME = ruby.1.9.1
arch = x86_64-darwin10.6.0
sitearch = $(arch)
ruby_version = 1.9.1
ruby = /Users/mbell/.rvm/rubies/ruby-1.9.2-head/bin/ruby
RUBY = $(ruby)
RM = rm -f
RM_RF = $(RUBY) -run -e rm -- -rf
RMDIRS = $(RUBY) -run -e rmdir -- -p
MAKEDIRS = mkdir -p
INSTALL = /usr/bin/install -c
INSTALL_PROG = $(INSTALL) -m 0755
INSTALL_DATA = $(INSTALL) -m 644
COPY = cp

#### End of system configuration section. ####

preload = 


CXX = g++

libpath = . $(libdir)
LIBPATH =  -L. -L$(libdir)
DEFFILE = 

CLEANFILES = mkmf.log
DISTCLEANFILES = 
DISTCLEANDIRS = 

extout = 
extout_prefix = 
target_prefix = 
LOCAL_LIBS = 
LIBS = -lrice -lruby.1.9.1 -lpthread -ldl -lobjc  -F/Library/Frameworks -L/Library/Frameworks -framework QtGui -framework QtCore
SRCS = AirborneRadarQC.cpp Dorade.cpp main.cpp QCscript.cpp RecursiveFilter.cpp
OBJS = AirborneRadarQC.o Dorade.o main.o QCscript.o RecursiveFilter.o
TARGET = QCscript
DLLIB = $(TARGET).bundle
EXTSTATIC = 
STATIC_LIB = 

BINDIR        = $(bindir)
RUBYCOMMONDIR = $(sitedir)$(target_prefix)
RUBYLIBDIR    = $(sitelibdir)$(target_prefix)
RUBYARCHDIR   = $(sitearchdir)$(target_prefix)
HDRDIR        = $(rubyhdrdir)/ruby$(target_prefix)
ARCHHDRDIR    = $(rubyhdrdir)/$(arch)/ruby$(target_prefix)

TARGET_SO     = $(DLLIB)
CLEANLIBS     = $(TARGET).bundle 
CLEANOBJS     = *.o  *.bak

all:    $(DLLIB)
static: $(STATIC_LIB)
.PHONY: all install static install-so install-rb
.PHONY: clean clean-so clean-rb

clean-rb-default::
clean-rb::
clean-so::
clean: clean-so clean-rb-default clean-rb
		@-$(RM) $(CLEANLIBS) $(CLEANOBJS) $(CLEANFILES)

distclean-rb-default::
distclean-rb::
distclean-so::
distclean: clean distclean-so distclean-rb-default distclean-rb
		@-$(RM) Makefile $(RUBY_EXTCONF_H) conftest.* mkmf.log
		@-$(RM) core ruby$(EXEEXT) *~ $(DISTCLEANFILES)
		@-$(RMDIRS) $(DISTCLEANDIRS)

realclean: distclean
install: install-so install-rb

install-so: $(RUBYARCHDIR)
install-so: $(RUBYARCHDIR)/$(DLLIB)
$(RUBYARCHDIR)/$(DLLIB): $(DLLIB)
	@-$(MAKEDIRS) $(@D)
	$(INSTALL_PROG) $(DLLIB) $(@D)
install-rb: pre-install-rb install-rb-default
install-rb-default: pre-install-rb-default
pre-install-rb: Makefile
pre-install-rb-default: Makefile
$(RUBYARCHDIR):
	$(MAKEDIRS) $@

site-install: site-install-so site-install-rb
site-install-so: install-so
site-install-rb: install-rb

.SUFFIXES: .c .m .cc .cxx .cpp .C .o

.cc.o:
	$(CXX) $(INCFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(COUTFLAG)$@ -c $<

.cxx.o:
	$(CXX) $(INCFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(COUTFLAG)$@ -c $<

.cpp.o:
	$(CXX) $(INCFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(COUTFLAG)$@ -c $<

.C.o:
	$(CXX) $(INCFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(COUTFLAG)$@ -c $<

.c.o:
	$(CC) $(INCFLAGS) $(CPPFLAGS) $(CFLAGS) $(COUTFLAG)$@ -c $<

$(DLLIB): $(OBJS) Makefile
	@-$(RM) $(@)
	$(LDSHAREDXX) -o $@ $(OBJS) $(LIBPATH) $(DLDFLAGS) $(LOCAL_LIBS) $(LIBS)



$(OBJS): $(hdrdir)/ruby.h $(hdrdir)/ruby/defines.h $(arch_hdrdir)/ruby/config.h
