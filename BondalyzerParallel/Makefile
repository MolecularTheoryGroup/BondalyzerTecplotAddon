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
	$(OBJDIR)/arrlisttools.o \
	$(OBJDIR)/atomcagelist.o \
	$(OBJDIR)/bondatompairs.o \
	$(OBJDIR)/bondbundle.o \
	$(OBJDIR)/bundles.o \
	$(OBJDIR)/createvolumezone.o \
	$(OBJDIR)/critpoints.o \
	$(OBJDIR)/elemshapefunc.o \
	$(OBJDIR)/engine.o \
	$(OBJDIR)/geomtools.o \
	$(OBJDIR)/gradpath.o \
	$(OBJDIR)/guicb.o \
	$(OBJDIR)/guidefs.o \
	$(OBJDIR)/integrate.o \
	$(OBJDIR)/isosurfgradpath.o \
	$(OBJDIR)/isotopo.o \
	$(OBJDIR)/linearalgebra.o \
	$(OBJDIR)/main.o \
	$(OBJDIR)/midplanegradpath.o \
	$(OBJDIR)/normals.o \
	$(OBJDIR)/oderungekutta.o \
	$(OBJDIR)/optimize.o \
	$(OBJDIR)/orderedpairs.o \
	$(OBJDIR)/spheregradpath.o \
	$(OBJDIR)/surfacefit.o \
	$(OBJDIR)/surfelemmap.o \
	$(OBJDIR)/surftoposeg.o \
	$(OBJDIR)/svd.o \
	$(OBJDIR)/zonevarinfo.o


#
# The next two dependencies probably can remain unchanged...
#

$(LIB) : $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) $(ADKUTILLIB) -o $(LIB)

$(STATLIB) : $(OBJS)
	$(AR) $(ARFLAGS) lib$(LIBBASENAME).a $(OBJS) ; mv lib$(LIBBASENAME).a  $(STATLIB)


$(OBJDIR)/arrlist.o  : arrlist.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/arrlist.cpp 

$(OBJDIR)/arrlisttools.o  : arrlisttools.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/arrlisttools.cpp 

$(OBJDIR)/atomcagelist.o  : atomcagelist.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/atomcagelist.cpp 

$(OBJDIR)/bondatompairs.o  : bondatompairs.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/bondatompairs.cpp 

$(OBJDIR)/bondbundle.o  : bondbundle.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/bondbundle.cpp 

$(OBJDIR)/bundles.o  : bundles.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/bundles.cpp 

$(OBJDIR)/createvolumezone.o  : createvolumezone.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/createvolumezone.cpp 

$(OBJDIR)/critpoints.o  : critpoints.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/critpoints.cpp 

$(OBJDIR)/elemshapefunc.o  : elemshapefunc.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/elemshapefunc.cpp 

$(OBJDIR)/engine.o  : engine.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/engine.cpp 

$(OBJDIR)/geomtools.o  : geomtools.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/geomtools.cpp 

$(OBJDIR)/gradpath.o  : gradpath.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/gradpath.cpp 

$(OBJDIR)/guicb.o  : guicb.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/guicb.cpp 

$(OBJDIR)/guidefs.o  : guidefs.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/guidefs.cpp 

$(OBJDIR)/integrate.o  : integrate.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/integrate.cpp 

$(OBJDIR)/isosurfgradpath.o  : isosurfgradpath.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/isosurfgradpath.cpp 

$(OBJDIR)/isotopo.o  : isotopo.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/isotopo.cpp 

$(OBJDIR)/linearalgebra.o  : linearalgebra.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/linearalgebra.cpp 

$(OBJDIR)/main.o  : main.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/main.cpp 

$(OBJDIR)/midplanegradpath.o  : midplanegradpath.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/midplanegradpath.cpp 

$(OBJDIR)/normals.o  : normals.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/normals.cpp 

$(OBJDIR)/oderungekutta.o  : oderungekutta.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/oderungekutta.cpp 

$(OBJDIR)/optimize.o  : optimize.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/optimize.cpp 

$(OBJDIR)/orderedpairs.o  : orderedpairs.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/orderedpairs.cpp 

$(OBJDIR)/spheregradpath.o  : spheregradpath.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/spheregradpath.cpp 

$(OBJDIR)/surfacefit.o  : surfacefit.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/surfacefit.cpp 

$(OBJDIR)/surfelemmap.o  : surfelemmap.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/surfelemmap.cpp 

$(OBJDIR)/surftoposeg.o  : surftoposeg.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/surftoposeg.cpp 

$(OBJDIR)/svd.o  : svd.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/svd.cpp 

$(OBJDIR)/zonevarinfo.o  : zonevarinfo.cpp
	cd $(OBJDIR); $(CPPCOMP) $(CPPFLAGS) $(SOURCEDIR)/zonevarinfo.cpp 


# DO NOT DELETE

arrlist.o: ALLOC.h ARRLIST.h
arrlisttools.o: ADDGLBL.h ADDONVER.h ALLOC.h ARRLIST.h ARRLISTTOOLS.h
atomcagelist.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ENGINE.h ALLOC.h ARRLIST.h
atomcagelist.o: ATOMCAGELIST.h
bondatompairs.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ENGINE.h ALLOC.h ARRLIST.h
bondatompairs.o: BONDATOMPAIRS.h
bondbundle.o: ADDGLBL.h ALLOC.h ARRLIST.h GEOMTOOLS.h SURFELEMMAP.h NORMALS.h
bondbundle.o: ZONEVARINFO.h CRITPOINTS.h BUNDLES.h BONDBUNDLE.h GRADPATH.h
bundles.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ENGINE.h ALLOC.h ARRLIST.h
bundles.o: BUNDLES.h
createvolumezone.o: ADDGLBL.h ALLOC.h ARRLIST.h GEOMTOOLS.h SURFELEMMAP.h
createvolumezone.o: NORMALS.h ZONEVARINFO.h CRITPOINTS.h BUNDLES.h
createvolumezone.o: BONDBUNDLE.h GRADPATH.h CREATEVOLUMEZONE.h
critpoints.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ENGINE.h ALLOC.h ARRLIST.h
critpoints.o: ARRLISTTOOLS.h LINEARALGEBRA.h ODERUNGEKUTTA.h SURFELEMMAP.h
critpoints.o: NORMALS.h SURFACEFIT.h ZONEVARINFO.h CRITPOINTS.h SVD.h
critpoints.o: ELEMSHAPEFUNC.h GRADPATH.h
elemshapefunc.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ENGINE.h ALLOC.h
elemshapefunc.o: LINEARALGEBRA.h ELEMSHAPEFUNC.h
engine.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ENGINE.h ALLOC.h ARRLIST.h
engine.o: LINEARALGEBRA.h ODERUNGEKUTTA.h NORMALS.h SURFELEMMAP.h
engine.o: ZONEVARINFO.h CRITPOINTS.h GRADPATH.h BUNDLES.h MIDPLANEGRADPATH.h
engine.o: SPHEREGRADPATH.h GEOMTOOLS.h SURFTOPOSEG.h ISOSURFGRADPATH.h
engine.o: ORDEREDPAIRS.h BONDBUNDLE.h INTEGRATE.h CREATEVOLUMEZONE.h
geomtools.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ENGINE.h GEOMTOOLS.h
gradpath.o: ADDGLBL.h ALLOC.h ARRLIST.h ARRLISTTOOLS.h ODERUNGEKUTTA.h
gradpath.o: NORMALS.h SURFELEMMAP.h ZONEVARINFO.h CRITPOINTS.h GRADPATH.h
gradpath.o: ENGINE.h GEOMTOOLS.h ELEMSHAPEFUNC.h
guicb.o: ADDGLBL.h GUIDEFS.h ENGINE.h guibld.cpp
guidefs.o: GUIDEFS.h
integrate.o: INTEGRATE.h
isosurfgradpath.o: ADDGLBL.h ALLOC.h ARRLIST.h ODERUNGEKUTTA.h NORMALS.h
isosurfgradpath.o: SURFELEMMAP.h ZONEVARINFO.h CRITPOINTS.h BUNDLES.h
isosurfgradpath.o: GRADPATH.h GEOMTOOLS.h SURFTOPOSEG.h ISOSURFGRADPATH.h
isosurfgradpath.o: SURFACEFIT.h ENGINE.h
isotopo.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ENGINE.h ALLOC.h ARRLIST.h
isotopo.o: ISOTOPO.h SURFELEMMAP.h
linearalgebra.o: ALLOC.h LINEARALGEBRA.h
main.o: ADDGLBL.h ADDONVER.h GUIDEFS.h
midplanegradpath.o: ADDGLBL.h ALLOC.h ARRLIST.h ODERUNGEKUTTA.h NORMALS.h
midplanegradpath.o: SURFELEMMAP.h ZONEVARINFO.h CRITPOINTS.h BUNDLES.h
midplanegradpath.o: GRADPATH.h MIDPLANEGRADPATH.h SURFACEFIT.h ENGINE.h
normals.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ENGINE.h ALLOC.h ARRLIST.h
normals.o: SURFELEMMAP.h NORMALS.h
oderungekutta.o: ODERUNGEKUTTA.h
optimize.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ENGINE.h ALLOC.h ARRLIST.h
optimize.o: OPTIMIZE.h
orderedpairs.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ENGINE.h ALLOC.h ARRLIST.h
orderedpairs.o: ORDEREDPAIRS.h
spheregradpath.o: ADDGLBL.h ALLOC.h ARRLIST.h ODERUNGEKUTTA.h NORMALS.h
spheregradpath.o: SURFELEMMAP.h ZONEVARINFO.h CRITPOINTS.h GRADPATH.h
spheregradpath.o: SPHEREGRADPATH.h ENGINE.h
surfacefit.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ENGINE.h ALLOC.h ARRLIST.h
surfacefit.o: SURFACEFIT.h SVD.h
surfelemmap.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ENGINE.h ALLOC.h ARRLIST.h
surfelemmap.o: SURFELEMMAP.h
surftoposeg.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ENGINE.h ALLOC.h ARRLIST.h
surftoposeg.o: ARRLISTTOOLS.h SURFTOPOSEG.h
svd.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ENGINE.h ALLOC.h ARRLIST.h SVD.h
zonevarinfo.o: ADDGLBL.h ADDONVER.h GUIDEFS.h ENGINE.h ALLOC.h ARRLIST.h
zonevarinfo.o: LINEARALGEBRA.h ODERUNGEKUTTA.h NORMALS.h SURFELEMMAP.h
zonevarinfo.o: ZONEVARINFO.h
