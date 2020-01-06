SAGE_SPKG_CONFIGURE([nauty], [
  AC_CHECK_PROG(DIRECTG_CHECK,directg,yes)
  AC_CHECK_PROG(GENTOURNG_CHECK,gentourng,yes)
  AC_CHECK_PROG(GENG_CHECK,geng,yes)
  AC_CHECK_PROG(GENBG_CHECK,genbg,yes)
  AS_IF([test x"$DIRECTG_CHECK" != "xyes"],[sage_spkg_install_nauty=yes])
  AS_IF([test x"$GENTOURNG_CHECK" != "xyes"],[sage_spkg_install_nauty=yes])
  AS_IF([test x"$GENG_CHECK" != "xyes"],[sage_spkg_install_nauty=yes])
  AS_IF([test x"$GENBG_CHECK" != "xyes"],[sage_spkg_install_nauty=yes])
])
