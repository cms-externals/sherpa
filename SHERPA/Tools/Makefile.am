pkglib_LTLIBRARIES       =	libSherpaTools.la

GITTAG = SHERPA/Tools
include ../../git.make

libSherpaTools_la_SOURCES = \
    Git_Info.C \
    Output_Base.C \
    Event_Reader_Base.C \
    HepEvt_Interface.C \
    HepMC2_Interface.C \
    HepMC3_Interface.C \
    Pythia_HepEvt_Translator.C \
    Output_HepEvt.C \
    Output_LHEF.C \
    Analysis_Interface.C \
    Userhook_Base.C \
    Variations.C

localinc_HEADERS = \
    Definitions.H \
    Output_Base.H \
    Event_Reader_Base.H \
    HepEvt_Interface.H \
    HepMC2_Interface.H \
    HepMC3_Interface.H \
    Pythia_HepEvt_Translator.H \
    Output_HepEvt.H \
    Output_LHEF.H \
    Analysis_Interface.H \
    Userhook_Base.H \
    Variations.H

libSherpaTools_la_LIBADD = @CONDITIONAL_HEPMC3LIBS@ @CONDITIONAL_HEPMC2LIBS@  @CONDITIONAL_LHAPDFLIBS@

libSherpaTools_la_CXXFLAGS = $(AM_CXXFLAGS)

libSherpaTools_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_GZIPINCS@ @CONDITIONAL_HEPMC3INCS@ @CONDITIONAL_HEPMC2INCS@ @CONDITIONAL_LHAPDFINCS@

