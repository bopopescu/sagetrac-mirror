"""
Utility functions for libGAP
"""

#*****************************************************************************
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function, absolute_import

from cpython.exc cimport PyErr_SetObject
from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_NE, Py_GT, Py_GE
from cysignals.signals cimport sig_on, sig_off, sig_error

from .gap_includes cimport *
from .element cimport *
from sage.cpython.string import FS_ENCODING
from sage.cpython.string cimport str_to_bytes, char_to_str
from sage.interfaces.gap_workspace import prepare_workspace_dir
from sage.env import SAGE_LOCAL, GAP_ROOT_DIR

# these do nothing now
cdef void libgap_enter() {}
cdef void libgap_exit() {}

############################################################################
### Hooking into the GAP memory management #################################
############################################################################

cdef class ObjWrapper(object):
    """
    Wrapper for GAP master pointers

    EXAMPLES::

        sage: from sage.libs.gap.util import ObjWrapper
        sage: x = ObjWrapper()
        sage: y = ObjWrapper()
        sage: x == y
        True
    """

    def __richcmp__(ObjWrapper self, ObjWrapper other, int op):
        r"""
        Comparison wrapped Obj.

        INPUT:

        - ``lhs``, ``rhs`` -- :class:`ObjWrapper`.

        - ``op`` -- integer. The comparison operation to be performed.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.gap.util import ObjWrapper
            sage: x = ObjWrapper()
            sage: y = ObjWrapper()
            sage: x == y
            True
        """
        cdef result
        cdef Obj self_value = self.value
        cdef Obj other_value = other.value
        if op == Py_LT:
            return self_value < other_value
        elif op == Py_LE:
            return self_value <= other_value
        elif op == Py_EQ:
            return self_value == other_value
        elif op == Py_GT:
            return self_value > other_value
        elif op == Py_GE:
            return self_value >= other_value
        elif op == Py_NE:
            return self_value != other_value
        else:
            assert False  # unreachable

    def __hash__(self):
        """
        Return a hash value

        EXAMPLES::

            sage: from sage.libs.gap.util import ObjWrapper
            sage: x = ObjWrapper()
            sage: hash(x)
            0
        """
        return <int>(self.value)


cdef ObjWrapper wrap_obj(Obj obj):
    """
    Constructor function for :class:`ObjWrapper`
    """
    cdef ObjWrapper result = ObjWrapper.__new__(ObjWrapper)
    result.value = obj
    return result


# a dictionary to keep all GAP elements
# needed for GASMAN callbacks
#
cdef dict owned_objects_refcount = dict()

#
# used in Sage's libgap.Gap.count_GAP_objects
#
cpdef get_owned_objects():
    """
    Helper to access the refcount dictionary from Python code
    """
    return owned_objects_refcount


cdef void reference_obj(Obj obj):
    """
    Reference ``obj``
    """
    cdef ObjWrapper wrapped = wrap_obj(obj)
    global owned_objects_refcount
    if wrapped in owned_objects_refcount:
        owned_objects_refcount[wrapped] += 1
    else:
        owned_objects_refcount[wrapped] = 1


cdef void dereference_obj(Obj obj):
    """
    Reference ``obj``
    """
    cdef ObjWrapper wrapped = wrap_obj(obj)
    global owned_objects_refcount
    refcount = owned_objects_refcount.pop(wrapped)
    if refcount > 1:
        owned_objects_refcount[wrapped] = refcount - 1


cdef void gasman_callback():
    """
    Callback before each GAP garbage collection
    """
    global owned_objects_refcount
    for obj in owned_objects_refcount:
        MarkBag((<ObjWrapper>obj).value)





############################################################################
### Initialization of libGAP ###############################################
############################################################################

def gap_root():
    """
    Find the location of the GAP root install which is stored in the gap
    startup script.

    EXAMPLES::

        sage: from sage.libs.gap.util import gap_root
        sage: gap_root()   # random output
        '/home/vbraun/opt/sage-5.3.rc0/local/gap/latest'
    """
    import os.path
    if os.path.exists(GAP_ROOT_DIR):
        return GAP_ROOT_DIR
    print('The gap-4.5.5.spkg (or later) seems to be not installed!')
    gap_sh = open(os.path.join(SAGE_LOCAL, 'bin', 'gap')).read().splitlines()
    gapdir = filter(lambda dir:dir.strip().startswith('GAP_DIR'), gap_sh)[0]
    gapdir = gapdir.split('"')[1]
    gapdir = gapdir.replace('$SAGE_LOCAL', SAGE_LOCAL)
    return gapdir


# To ensure that we call initialize_libgap only once.
cdef bint _gap_is_initialized = False
cdef extern char **environ

cdef initialize():
    """
    Initialize the GAP library, if it hasn't already been
    initialized.  It is safe to call this multiple times.

    TESTS::

        sage: libgap(123)   # indirect doctest
        123
    """
    global _gap_is_initialized, environ
    if _gap_is_initialized: return

    # Define argv and environ variables, which we will pass in to
    # initialize GAP. Note that we must pass define the memory pool
    # size!
    cdef char* argv[14]
    argv[0] = "sage"
    argv[1] = "-l"
    s = str_to_bytes(gap_root(), FS_ENCODING, "surrogateescape")
    argv[2] = s

    from sage.interfaces.gap import _get_gap_memory_pool_size_MB
    memory_pool = str_to_bytes(_get_gap_memory_pool_size_MB())
    argv[3] = "-o"
    argv[4] = memory_pool
    argv[5] = "-s"
    argv[6] = memory_pool

    argv[7] = "-m"
    argv[8] = "64m"

    argv[9] = "-q"    # no prompt!
    argv[10] = "-T"    # no debug loop
    argv[11] = NULL
    cdef int argc = 11   # argv[argc] must be NULL

    from .saved_workspace import workspace
    workspace, workspace_is_up_to_date = workspace()
    ws = str_to_bytes(workspace, FS_ENCODING, "surrogateescape")
    if workspace_is_up_to_date:
        argv[11] = "-L"
        argv[12] = ws
        argv[13] = NULL
        argc = 13

    # Initialize GAP and capture any error messages
    # The initialization just prints error and does not use the error handler
    # libgap_initialize(argc, argv)
    # gap_error_msg = char_to_str(libgap_get_output())
    if GAP_Initialize(argc, argv, environ, &gasman_callback, &error_handler):
        raise RuntimeError('libGAP initialization failed\n' + gap_error_msg)

    # The error handler is called if a GAP evaluation fails, e.g. 1/0
    # libgap_set_error_handler(&error_handler)

    # Prepare global GAP variable to hold temporary GAP objects
    global reference_holder
    libgap_enter()
    reference_holder = GVarName("$SAGE_libgap_reference_holder")
    libgap_exit()

    # Finished!
    _gap_is_initialized = True

    # Save a new workspace if necessary
    if not workspace_is_up_to_date:
        prepare_workspace_dir()
        from sage.misc.temporary_file import atomic_write
        with atomic_write(workspace) as f:
            f.close()
            gap_eval('SaveWorkspace("{0}")'.format(f.name))


############################################################################
### Evaluate string in GAP #################################################
############################################################################

cdef Obj gap_eval(str gap_string) except? NULL:
    r"""
    Evaluate a string in GAP.

    INPUT:

    - ``gap_string`` -- string. A valid statement in GAP.

    OUTPUT:

    The resulting GAP object or NULL+Python Exception in case of error.

    EXAMPLES::

        sage: libgap.eval('if 4>3 then\nPrint("hi");\nfi')
        NULL
        sage: libgap.eval('1+1')   # testing that we have successfully recovered
        2

        sage: libgap.eval('if 4>3 thenPrint("hi");\nfi')
        Traceback (most recent call last):
        ...
        ValueError: libGAP: Syntax error: then expected
        if 4>3 thenPrint("hi");
        fi;
                       ^
        sage: libgap.eval('1+1')   # testing that we have successfully recovered
        2

    TESTS:

    Check that we fail gracefully if this is called within
    ``libgap_enter()``::

        sage: cython('''
        ....: # distutils: libraries = gap
        ....: from sage.libs.gap.gap_includes cimport libgap_enter
        ....: libgap_enter()
        ....: ''')
        sage: libgap.eval('1+1')
        Traceback (most recent call last):
        ...
        ValueError: libGAP: Entered a critical block twice
    """
    initialize()
    cdef ExecStatus status
    cdef Obj *result

    # Careful: We need to keep a reference to the bytes object here
    # so that Cython doesn't dereference it before libGAP is done with
    # its contents.
    cmd = str_to_bytes(gap_string + ';\n')
    try:
        try:
            sig_on()
            libgap_enter()
            status = ReadEvalCommand(BottomLVars, result, NULL)
            if status != STATUS_END:
                # libgap_call_error_handler()
                if error_func==NULL:
                        print("An error occurred, but libGAP has no handler set")
                        return # needs work
            sig_off()
        except RuntimeError as msg:
            raise ValueError('libGAP: '+str(msg).strip())

        if Symbol != S_SEMICOLON:
            raise ValueError('did not end with semicolon')
        GetSymbol()
        if Symbol != S_EOF:
            raise ValueError('can only evaluate a single statement')

    finally:
        libgap_exit()

    if ReadEvalResult != NULL:
        libgap_enter()
        AssGVar(Last3, VAL_GVAR(Last2))
        AssGVar(Last2, VAL_GVAR(Last))
        AssGVar(Last, ReadEvalResult)
        libgap_exit()

    return result # may be NULL, thats ok


############################################################################
### Helper to protect temporary objects from deletion ######################
############################################################################

# Hold a reference (inside the GAP kernel) to obj so that it doesn't
# get deleted this works by assigning it to a global variable. This is
# very simple, but you can't use it to keep two objects alive. Be
# careful.
cdef UInt reference_holder

cdef void hold_reference(Obj obj):
    """
    Hold a reference (inside the GAP kernel) to obj

    This ensures that the GAP garbage collector does not delete
    ``obj``. This works by assigning it to a global variable. This is
    very simple, but you can't use it to keep two objects alive. Be
    careful.
    """
    libgap_enter()
    global reference_holder
    AssGVar(reference_holder, obj)
    libgap_exit()


############################################################################
### Error handler ##########################################################
############################################################################

cdef void error_handler(char* msg):
    """
    The libgap error handler

    We call ``sig_error()`` which causes us to jump back to the Sage
    signal handler. Since we wrap libGAP C calls in ``sig_on`` /
    ``sig_off`` blocks, this then jumps back to the ``sig_on`` where
    the ``RuntimeError`` we raise here will be seen.
    """
    msg_py = char_to_str(msg)
    msg_py = msg_py.replace('For debugging hints type ?Recovery from NoMethodFound\n', '')
    PyErr_SetObject(RuntimeError, msg_py)
    sig_error()


############################################################################
### Debug functions ########################################################
############################################################################

cdef inline void DEBUG_CHECK(Obj obj):
    """
    Check that ``obj`` is valid.

    This function is only useful for debugging.
    """
    libgap_enter()
    CheckMasterPointers()
    libgap_exit()
    if obj == NULL:
        print('DEBUG_CHECK: Null pointer!')




cpdef memory_usage():
    """
    Return information about the memory usage.

    See :meth:`~sage.libs.gap.libgap.Gap.mem` for details.
    """
    cdef size_t SizeMptrsArea = OldBags - MptrBags
    cdef size_t SizeOldBagsArea = YoungBags - OldBags
    cdef size_t SizeYoungBagsArea = AllocBags - YoungBags
    cdef size_t SizeAllocationArea = StopBags - AllocBags
    cdef size_t SizeUnavailableArea = EndBags - StopBags
    return (SizeMptrsArea, SizeOldBagsArea, SizeYoungBagsArea, SizeAllocationArea, SizeUnavailableArea)


cpdef error_enter_libgap_block_twice():
    """
    Demonstrate that we catch errors from entering a block twice.

    EXAMPLES::

        sage: from sage.libs.gap.util import error_enter_libgap_block_twice
        sage: error_enter_libgap_block_twice()
        Traceback (most recent call last):
        ...
        RuntimeError: Entered a critical block twice
    """
    from sage.libs.gap.libgap import libgap
    try:
        # The exception will be seen by this sig_on() after being
        # raised by the second libgap_enter().
        sig_on()
        libgap_enter()
        libgap_enter()
        sig_off()
    finally:
        libgap_exit()


cpdef error_exit_libgap_block_without_enter():
    """
    Demonstrate that we catch errors from omitting libgap_enter.

    EXAMPLES::

        sage: from sage.libs.gap.util import error_exit_libgap_block_without_enter
        sage: error_exit_libgap_block_without_enter()
        Traceback (most recent call last):
        ...
        RuntimeError: Called libgap_exit without previous libgap_enter
    """
    from sage.libs.gap.libgap import libgap
    sig_on()
    libgap_exit()
    sig_off()

############################################################################
### Auxilliary functions ###################################################
############################################################################


def command(command_string):
    """
    Playground for accessing Gap via libGap.

    You should not use this function in your own programs. This is
    just here for convenience if you want to play with the libgap
    libray code.

    EXAMPLES::

        sage: from sage.libs.gap.util import command
        sage: command('1')
        Output follows...
        1

        sage: command('1/0')
        Traceback (most recent call last):
        ...
        ValueError: libGAP: Error, Rational operations: <divisor> must not be zero

        sage: command('NormalSubgroups')
        Output follows...
        <Attribute "NormalSubgroups">

        sage: command('rec(a:=1, b:=2)')
        Output follows...
        rec( a := 1, b := 2 )
    """
    initialize()
    cdef ExecStatus status
    cdef Obj *result
    cmd = str_to_bytes(command_string + ';\n')
    try:
        libgap_enter()
        try:
            sig_on()
            status = ReadEvalCommand(BottomLVars, result)
            if status != STATUS_END:
                print("oops - not ready...")  # needs work
            sig_off()
        except RuntimeError as msg:
            raise ValueError('libGAP: '+str(msg).strip())

        assert Symbol == S_SEMICOLON, 'Did not end with semicolon?'
        GetSymbol()
        if Symbol != S_EOF:
            raise ValueError('command() expects a single statement.')

        if result:
            print('Output follows...')
            ViewObjHandler(result)
        else:
            print('No output.')

    finally:
        libgap_exit()

    DEBUG_CHECK(result)

    if result != NULL:
        libgap_enter()
        AssGVar(Last3, VAL_GVAR(Last2))
        AssGVar(Last2, VAL_GVAR(Last))
        AssGVar(Last, ReadEvalResult)
        libgap_exit()

