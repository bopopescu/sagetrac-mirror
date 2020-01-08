SAGE_SPKG_CONFIGURE([nauty], [
  dnl At least Debian and Ubuntu prefix these executables with "nauty-",
  dnl so we have to check for both. Due to some packaging concerns, we
  dnl eventually store the paths in SageMath environment variables
  dnl rather than substituting them directly into the source code.
  AC_PATH_PROGS([DIRECTG_PATH],[directg nauty-directg])
  AS_IF([test x"$DIRECTG_PATH" = "x"], [sage_spkg_install_nauty=yes])

  AC_PATH_PROGS([GENTOURNG_PATH],[gentourng nauty-gentourng])
  AS_IF([test x"$GENTOURNG_PATH" = "x"], [sage_spkg_install_nauty=yes])

  AC_PATH_PROGS([GENG_PATH],[geng nauty-geng])
  AS_IF([test x"$GENG_PATH" = "x"], [sage_spkg_install_nauty=yes])

  AC_PATH_PROGS([GENBG_PATH],[genbg nauty-genbg])
  AS_IF([test x"$GENBG_PATH" = "x"], [sage_spkg_install_nauty=yes])
],[],[],[
  dnl These commands are always run AFTER the check, even
  dnl if --with-system-nauty=no was passed to ./configure.
  dnl If we're installing the SageMath spkg, we don't want
  dnl to use any of the paths that we may have found earlier;
  dnl instead we use the bare executable names because that's
  dnl what we used to do and it works with the spkg.
  AS_IF([test x"$sage_spkg_install_nauty" = "xyes"], [
           DIRECTG_PATH=directg
           GENTOURNG_PATH=gentourng
           GENG_PATH=geng
           GENBG_PATH=genbg
        ])
  AC_SUBST(DIRECTG_PATH)
  AC_SUBST(GENTOURNG_PATH)
  AC_SUBST(GENG_PATH)
  AC_SUBST(GENBG_PATH)
])
