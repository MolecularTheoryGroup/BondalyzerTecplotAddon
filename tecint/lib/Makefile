SHELL=/bin/sh
LIB=$(TECADDONDEVDIR)/lib/$(PLATFORM)/$(OBJSUBDIR)/lib$(LIBBASENAME).$(LIBEXT)
STATLIB=$(ADDONLIBDIR)/lib$(LIBBASENAME).a
OBJDIR=$(OBJBUILDDIR)/obj/$(PLATFORM)/$(OBJSUBDIR)
EXTLIBDIR=$(EXTBASEDIR)/$(PLATFORM)
IOUTILLIB=$(EXTLIBDIR)/libioutil.a
GUILIB=$(EXTLIBDIR)/libgui.a
 
shlib: $(LIB)
statlib: $(STATLIB)
#
# list of obj files....
#
OBJS= \
	$(OBJDIR)/integrate.o \
	$(OBJDIR)/iresult.o \
	$(OBJDIR)/pctdone.o \
	$(OBJDIR)/zoneinfo.o 
#
# The next two dependencies probably can remain unchanged...
#

$(LIB) : $(OBJS)
	cd $(OBJDIR); $(LD) $(LDFLAGS) *.o $(GUILIB) -o $(LIB)

$(STATLIB) : $(OBJS)
	cd $(OBJDIR); $(AR) $(ARFLAGS) lib$(LIBBASENAME).a *.o ; mv lib$(LIBBASENAME).a  $(STATLIB)



$(OBJDIR)/integrate.o  : integrate.c  \
                            iresult.h \
                            pctdone.h \
                            tecint.h \
                            zoneinfo.h 
	cd $(OBJDIR); $(CCOMP) $(CFLAGS) $(SOURCEDIR)/integrate.c 


$(OBJDIR)/iresult.o  : iresult.c  \
                            iresult.h \
                            tecint.h 
	cd $(OBJDIR); $(CCOMP) $(CFLAGS) $(SOURCEDIR)/iresult.c 


$(OBJDIR)/pctdone.o  : pctdone.c  \
                            pctdone.h 
	cd $(OBJDIR); $(CCOMP) $(CFLAGS) $(SOURCEDIR)/pctdone.c 


$(OBJDIR)/zoneinfo.o  : zoneinfo.c  \
                            zoneinfo.h 
	cd $(OBJDIR); $(CCOMP) $(CFLAGS) $(SOURCEDIR)/zoneinfo.c 
