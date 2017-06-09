SHELL=/bin/sh
LIB=$(TECADDONDEVDIR)/lib/$(PLATFORM)/$(OBJSUBDIR)/lib$(LIBBASENAME).$(LIBEXT)
STATLIB=$(ADDONLIBDIR)/lib$(LIBBASENAME).a
OBJDIR=$(OBJBUILDDIR)/obj/$(PLATFORM)/$(OBJSUBDIR)
EXTLIBDIR=$(EXTBASEDIR)/$(PLATFORM)
ADKUTILLIB=$(EXTLIBDIR)/libadkutil.a
IOUTILLIB=$(EXTLIBDIR)/libioutil.a
 
shlib: $(LIB)
statlib: $(STATLIB)
#
# list of obj files....
#
OBJS= \
	$(OBJDIR)/arrlist.o \
	$(OBJDIR)/bundles.o \
	$(OBJDIR)/critpoints.o \
	$(OBJDIR)/engine.o \
	$(OBJDIR)/guicb.o \
	$(OBJDIR)/guidefs.o \
	$(OBJDIR)/linearalgebra.o \
	$(OBJDIR)/main.o

#
# The next two dependencies probably can remain unchanged...
#

$(LIB) : $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(LIB)

$(STATLIB) : $(OBJS)
	$(AR) $(ARFLAGS) lib$(LIBBASENAME).a $(OBJS) ; mv lib$(LIBBASENAME).a  $(STATLIB)

$(OBJDIR)/arrlist.o  : arrlist.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/arrlist.cpp 

$(OBJDIR)/bundles.o  : bundles.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/bundles.cpp 

$(OBJDIR)/critpoints.o  : critpoints.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/critpoints.cpp 

$(OBJDIR)/engine.o  : engine.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/engine.cpp 

$(OBJDIR)/guicb.o  : guicb.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/guicb.cpp 

$(OBJDIR)/guidefs.o  : guidefs.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/guidefs.cpp 

$(OBJDIR)/linearalgebra.o  : linearalgebra.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/linearalgebra.cpp 

$(OBJDIR)/main.o  : main.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/main.cpp 



# DO NOT DELETE

$(OBJDIR)/arrlist.o: arrlist.cpp ALLOC.h ARRLIST.h
$(OBJDIR)/bundles.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ALLOC.h ARRLIST.h BUNDLES.h
$(OBJDIR)/bundles.o: CRITPOINTS.h ENGINE.h
$(OBJDIR)/critpoints.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ALLOC.h ARRLIST.h
$(OBJDIR)/critpoints.o: LINEARALGEBRA.h CRITPOINTS.h ENGINE.h BUNDLES.h
$(OBJDIR)/engine.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ALLOC.h ARRLIST.h CRITPOINTS.h
$(OBJDIR)/engine.o: BUNDLES.h ENGINE.h
$(OBJDIR)/guicb.o: ADDGLBL.h GUIDEFS.h ALLOC.h ARRLIST.h CRITPOINTS.h ENGINE.h
$(OBJDIR)/guicb.o: BUNDLES.h guibld.cpp
$(OBJDIR)/guidefs.o: GUIDEFS.h
$(OBJDIR)/linearalgebra.o: ALLOC.h LINEARALGEBRA.h
$(OBJDIR)/main.o: ADDGLBL.h ADDONVER.h GUIDEFS.h
