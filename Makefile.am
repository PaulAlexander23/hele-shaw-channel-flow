# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS=hele_shaw_channel_flow

#---------------------------------------------------------------------------

# Local sources that Jack's Poisson code depends on:
hele_shaw_channel_flow_SOURCES = hele_shaw_channel_flow.cc

# Required libraries: 
hele_shaw_channel_flow_LDADD = -L@libdir@ -lmeshes -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------
# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
AM_CPPFLAGS += -I@includedir@
