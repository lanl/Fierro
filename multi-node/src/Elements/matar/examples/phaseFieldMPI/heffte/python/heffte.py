'''
    Python Interface
    -- heFFTe --
    Univ. of Tennessee, Knoxville
'''

from heffte_config import *
import heffte_config
__version__ = heffte_config.__version__
from ctypes import sizeof, Structure, c_char_p, c_int, c_double, c_float, c_void_p, POINTER, cast, cdll
import sys, traceback
from numpy.ctypeslib import ndpointer

import numpy as np
from mpi4py import MPI as mpi

class heffte_input_error(Exception):
    '''
    Exception raided by incorrect input to heFFTe

    Initialized with a message that is saved in self.message
    and can be printed with the print_what() method.
    '''
    def __init__(self, message):
        '''
        Constructor, saves a copy of the message.
        '''
        self.message = message

    def print_what(self):
        '''
        Prints the message set by the constructor.
        '''
        print(message)

class heffte_plan(Structure):
    _fields_ = [ ("backend_type", c_int), ("using_r2c", c_int), ("fft", c_void_p) ]
LP_plan = POINTER(heffte_plan)

class plan_options(Structure):
    _fields_ = [ ("use_reorder", c_int), ("use_alltoall", c_int), ("use_pencils", c_int) ]

# double-check this!
MPI_Comm = c_int if mpi._sizeof(mpi.COMM_WORLD) == sizeof(c_int) else c_void_p

libheffte = cdll.LoadLibrary(libheffte_path)

# create and destroy
libheffte.heffte_plan_create.argtypes = [c_int, ndpointer(c_int, flags="C_CONTIGUOUS"), ndpointer(c_int, flags="C_CONTIGUOUS"), \
                                         ndpointer(c_int, flags="C_CONTIGUOUS"), ndpointer(c_int, flags="C_CONTIGUOUS"), \
                                         ndpointer(c_int, flags="C_CONTIGUOUS"), ndpointer(c_int, flags="C_CONTIGUOUS"), \
                                         MPI_Comm, POINTER(plan_options), POINTER(LP_plan)]
libheffte.heffte_plan_create.restype = c_int
libheffte.heffte_plan_create_r2c.argtypes = [c_int, ndpointer(c_int, flags="C_CONTIGUOUS"), ndpointer(c_int, flags="C_CONTIGUOUS"), \
                                             ndpointer(c_int, flags="C_CONTIGUOUS"), ndpointer(c_int, flags="C_CONTIGUOUS"), \
                                             ndpointer(c_int, flags="C_CONTIGUOUS"), ndpointer(c_int, flags="C_CONTIGUOUS"), \
                                             c_int, MPI_Comm, POINTER(plan_options), POINTER(LP_plan)]
libheffte.heffte_plan_create_r2c.restype = c_int

libheffte.heffte_plan_destroy.argtypes = [LP_plan]
libheffte.heffte_plan_destroy.restype = c_int

# forward and backward (unbuffered)
libheffte.heffte_forward_s2c.argtypes = [LP_plan, POINTER(c_float), c_void_p, c_int]
libheffte.heffte_forward_s2c.restype = None
libheffte.heffte_forward_c2c.argtypes = [LP_plan, c_void_p, c_void_p, c_int]
libheffte.heffte_forward_c2c.restype = None
libheffte.heffte_forward_d2z.argtypes = [LP_plan, POINTER(c_double), c_void_p, c_int]
libheffte.heffte_forward_d2z.restype = None
libheffte.heffte_forward_z2z.argtypes = [LP_plan, c_void_p, c_void_p, c_int]
libheffte.heffte_forward_z2z.restype = None

libheffte.heffte_backward_c2s.argtypes = [LP_plan, c_void_p, POINTER(c_float), c_int]
libheffte.heffte_backward_c2s.restype = None
libheffte.heffte_backward_c2c.argtypes = [LP_plan, c_void_p, c_void_p, c_int]
libheffte.heffte_backward_c2c.restype = None
libheffte.heffte_backward_z2d.argtypes = [LP_plan, c_void_p, POINTER(c_double), c_int]
libheffte.heffte_backward_z2d.restype = None
libheffte.heffte_backward_z2z.argtypes = [LP_plan, c_void_p, c_void_p, c_int]
libheffte.heffte_backward_z2z.restype = None

# forward and backward (buffered)
libheffte.heffte_forward_s2c_buffered.argtypes = [LP_plan, POINTER(c_float), c_void_p, c_void_p, c_int]
libheffte.heffte_forward_s2c_buffered.restype = None
libheffte.heffte_forward_c2c_buffered.argtypes = [LP_plan, c_void_p, c_void_p, c_void_p, c_int]
libheffte.heffte_forward_c2c_buffered.restype = None
libheffte.heffte_forward_d2z_buffered.argtypes = [LP_plan, POINTER(c_double), c_void_p, c_void_p, c_int]
libheffte.heffte_forward_d2z_buffered.restype = None
libheffte.heffte_forward_z2z_buffered.argtypes = [LP_plan, c_void_p, c_void_p, c_void_p, c_int]
libheffte.heffte_forward_z2z_buffered.restype = None

libheffte.heffte_backward_c2s_buffered.argtypes = [LP_plan, c_void_p, POINTER(c_float), c_void_p, c_int]
libheffte.heffte_backward_c2s_buffered.restype = None
libheffte.heffte_backward_c2c_buffered.argtypes = [LP_plan, c_void_p, c_void_p, c_void_p, c_int]
libheffte.heffte_backward_c2c_buffered.restype = None
libheffte.heffte_backward_z2d_buffered.argtypes = [LP_plan, c_void_p, POINTER(c_double), c_void_p, c_int]
libheffte.heffte_backward_z2d_buffered.restype = None
libheffte.heffte_backward_z2z_buffered.argtypes = [LP_plan, c_void_p, c_void_p, c_void_p, c_int]
libheffte.heffte_backward_z2z_buffered.restype = None

# read meta-data
libheffte.heffte_size_inbox.argtypes = [LP_plan]
libheffte.heffte_size_inbox.restype = c_int
libheffte.heffte_size_outbox.argtypes = [LP_plan]
libheffte.heffte_size_outbox.restype = c_int
libheffte.heffte_size_workspace.argtypes = [LP_plan]
libheffte.heffte_size_workspace.restype = c_int


class backend:
    stock = 0
    fftw = 1
    mkl = 2
    cufft = 10
    rocm = 11
    valid = [0, 1, 2, 10, 11]

class scale:
    none = 0
    full = 1
    symmetric = 2

class box3d:
    def __init__(self, clow, chigh, corder = np.array([0,1,2])):
        self.low = np.array(clow, dtype = np.int32)
        self.high = np.array(chigh, dtype = np.int32)
        self.order = np.array(corder, dtype = np.int32)
        assert self.low.size == 3
        assert self.high.size == 3
        assert self.order.size == 3
        for i in range(3):
            assert self.order[i] in [0, 1, 2]

def fft3d(backend_tag, inbox, outbox, comm):
    '''
    Initialize a fft3d operation, the syntax is near identical to C++.

    backend_tag replaces the template type-tag, use of the heffte.backend
                constants, e.g., heffte.backend.fftw
    '''
    if backend_tag not in backend.valid:
        raise heffte_input_error("Invalid backend, use one of the entries in heffte.backend")

    plan = heffte_fft_plan()
    plan.use_r2c = False

    # Define ctypes API for each library method
    comm_value = MPI_Comm.from_address( mpi._addressof(comm) )

    # Initialize
    plan.fft_comm = comm
    plan.plan = LP_plan()
    options = plan_options(0,1,1)

    herr = libheffte.heffte_plan_create(backend_tag, inbox.low, inbox.high, inbox.order,
                                        outbox.low, outbox.high, outbox.order,
                                        comm_value, options, plan.plan)

    if herr != 0:
        raise heffte_input_error("heFFTe encountered internal error with code: {0:1d}".format(herr))
    return plan

def fft3d_r2c(backend_tag, inbox, outbox, r2c_direction, comm):
    '''
    Initialize a fft3d_r2c operation, the syntax is near identical to C++.

    backend_tag replaces the template type-tag, use of the heffte.backend
                constants, e.g., heffte.backend.fftw
    '''
    if backend_tag not in backend.valid:
        raise heffte_input_error("Invalid backend, use one of the entries in heffte.backend")

    if r2c_direction not in [0, 1, 2]:
        raise heffte_input_error("fft3d_r2c() called with invalid r2c_direction, must use 0, 1, or 2")

    plan = heffte_fft_plan()
    plan.use_r2c = False

    # Define ctypes API for each library method
    comm_value = MPI_Comm.from_address( mpi._addressof(comm) )

    # Initialize
    plan.fft_comm = comm
    plan.plan = LP_plan()
    options = plan_options(0,1,1)

    herr = libheffte.heffte_plan_create_r2c(backend_tag, inbox.low, inbox.high, inbox.order,
                                            outbox.low, outbox.high, outbox.order,
                                            r2c_direction, comm_value, options, plan.plan)

    if herr != 0:
        raise heffte_input_error("heFFTe encountered internal error with code: {0:1d}".format(herr))
    return plan

class heffte_fft_plan:
    def __init__(self):
        '''
        Empty constructor, the object cannot be used after this call.

        Use one of the factory methods, either fft3d() or fft3d_r2c().
        '''
        pass

    def __del__(self):
        libheffte.heffte_plan_destroy(self.plan)

    def _get_ctypes_data(self, x, xtype):
        '''
        Internal use, extracts ctypes array from x
        '''
        if type(x) is np.ndarray:
            return x.ctypes.data_as(xtype)
        else:
            return cast(x.device_ctypes_pointer.value, xtype)

    def forward(self, inarray, outarray, scaling = scale.none):

        if scaling not in (0, 1, 2):
            raise heffte_input_error("forward() called with invalid scaling")

        if inarray.size != self.size_inbox() or outarray.size != self.size_outbox():
            raise heffte_input_error("forward() called with invalid array size")

        if self.use_r2c:
            if inarray.dtype == np.complex64 or inarray.dtype == np.complex128:
                raise heffte_input_error("forward() called with r2c can use only real inarray types")

        if inarray.dtype == np.float32:
            if outarray.dtype != np.complex64:
                raise heffte_input_error("forward() with inarray.dtype == float32 needs outarray.dtype == complex64")
            libheffte.heffte_forward_s2c(self.plan,
                                         self._get_ctypes_data(inarray, POINTER(c_float)),
                                         self._get_ctypes_data(outarray, c_void_p),
                                         scaling)
        elif inarray.dtype == np.complex64:
            if outarray.dtype != np.complex64:
                raise heffte_input_error("forward() with inarray.dtype == complex64 needs outarray.dtype == complex64")
            libheffte.heffte_forward_c2c(self.plan,
                                         self._get_ctypes_data(inarray, c_void_p),
                                         self._get_ctypes_data(outarray, c_void_p),
                                         scaling)
        elif inarray.dtype == np.float64:
            if outarray.dtype != np.complex128:
                raise heffte_input_error("forward() with inarray.dtype == float64 needs outarray.dtype == complex128")
            libheffte.heffte_forward_d2z(self.plan,
                                         self._get_ctypes_data(inarray, POINTER(c_double)),
                                         self._get_ctypes_data(outarray, c_void_p),
                                         scaling)
        elif inarray.dtype == np.complex128:
            if outarray.dtype != np.complex128:
                raise heffte_input_error("forward() with inarray.dtype == complex128 needs outarray.dtype == complex128")
            libheffte.heffte_forward_z2z(self.plan,
                                         self._get_ctypes_data(inarray, c_void_p),
                                         self._get_ctypes_data(outarray, c_void_p),
                                         scaling)
        else:
            raise heffte_input_error("forward() called with wrong inarray.dtype, use float32, float64, complex64, or complex128")

    def forward_buffered(self, inarray, outarray, workspace, scaling = scale.none):

        if scaling not in (0, 1, 2):
            raise heffte_input_error("forward() called with invalid scaling")

        if inarray.size != self.size_inbox() or outarray.size != self.size_outbox():
            raise heffte_input_error("forward() called with invalid array size")

        if workspace.size < self.size_workspace():
            raise heffte_input_error("forward() called with insufficient workspace size")

        if self.use_r2c:
            if inarray.dtype == np.complex64 or inarray.dtype == np.complex128:
                raise heffte_input_error("forward() called with r2c can use only real inarray types")

        if (inarray.dtype == np.float32 or inarray.dtype == np.complex64) and workspace.dtype != np.complex64:
            raise heffte_input_error("forward() called with incorrect workspace type, must be complex64")
        if (inarray.dtype == np.float64 or inarray.dtype == np.complex128) and workspace.dtype != np.complex128:
            raise heffte_input_error("forward() called with incorrect workspace type, must be complex128")

        if inarray.dtype == np.float32:
            if outarray.dtype != np.complex64:
                raise heffte_input_error("forward() with inarray.dtype == float32 needs outarray.dtype == complex64")
            libheffte.heffte_forward_s2c_buffered(self.plan,
                                                  self._get_ctypes_data(inarray, POINTER(c_float)),
                                                  self._get_ctypes_data(outarray, c_void_p),
                                                  self._get_ctypes_data(workspace, c_void_p),
                                                  scaling)
        elif inarray.dtype == np.complex64:
            if outarray.dtype != np.complex64:
                raise heffte_input_error("forward() with inarray.dtype == complex64 needs outarray.dtype == complex64")
            libheffte.heffte_forward_c2c_buffered(self.plan,
                                                  self._get_ctypes_data(inarray, c_void_p),
                                                  self._get_ctypes_data(outarray, c_void_p),
                                                  self._get_ctypes_data(workspace, c_void_p),
                                                  scaling)
        elif inarray.dtype == np.float64:
            if outarray.dtype != np.complex128:
                raise heffte_input_error("forward() with inarray.dtype == float64 needs outarray.dtype == complex128")
            libheffte.heffte_forward_d2z_buffered(self.plan,
                                                  self._get_ctypes_data(inarray, POINTER(c_double)),
                                                  self._get_ctypes_data(outarray, c_void_p),
                                                  self._get_ctypes_data(workspace, c_void_p),
                                                  scaling)
        elif inarray.dtype == np.complex128:
            if outarray.dtype != np.complex128:
                raise heffte_input_error("forward() with inarray.dtype == complex128 needs outarray.dtype == complex128")
            libheffte.heffte_forward_z2z_buffered(self.plan,
                                                  self._get_ctypes_data(inarray, c_void_p),
                                                  self._get_ctypes_data(outarray, c_void_p),
                                                  self._get_ctypes_data(workspace, c_void_p),
                                                  scaling)
        else:
            raise heffte_input_error("forward() called with wrong inarray.dtype, use float32, float64, complex64, or complex128")

    def backward(self, inarray, outarray, scaling = scale.none):

        if scaling not in (0, 1, 2):
            raise heffte_input_error("backward() called with invalid scaling")

        if inarray.size != self.size_inbox() or outarray.size != self.size_outbox():
            raise heffte_input_error("backward() called with invalid array size")

        if self.use_r2c:
            if outarray.dtype == np.complex64 or outarray.dtype == np.complex128:
                raise heffte_input_error("backward() called with r2c can use only real outarray types")

        if outarray.dtype == np.float32:
            if inarray.dtype != np.complex64:
                raise heffte_input_error("backward() with outarray.dtype == float32 needs inarray.dtype == complex64")
            libheffte.heffte_backward_c2s(self.plan,
                                          self._get_ctypes_data(inarray, c_void_p),
                                          self._get_ctypes_data(outarray, POINTER(c_float)),
                                          scaling)
        elif outarray.dtype == np.complex64:
            if inarray.dtype != np.complex64:
                raise heffte_input_error("backward() with outarray.dtype == complex64 needs inarray.dtype == complex64")
            libheffte.heffte_backward_c2c(self.plan,
                                          self._get_ctypes_data(inarray, c_void_p),
                                          self._get_ctypes_data(outarray, c_void_p),
                                          scaling)
        elif outarray.dtype == np.float64:
            if inarray.dtype != np.complex128:
                raise heffte_input_error("backward() with outarray.dtype == float64 needs inarray.dtype == complex128")
            libheffte.heffte_backward_z2d(self.plan,
                                          self._get_ctypes_data(inarray, c_void_p),
                                          self._get_ctypes_data(outarray, POINTER(c_double)),
                                          scaling)
        elif outarray.dtype == np.complex128:
            if inarray.dtype != np.complex128:
                raise heffte_input_error("backward() with outarray.dtype == complex128 needs inarray.dtype == complex128")
            libheffte.heffte_backward_z2z(self.plan,
                                          self._get_ctypes_data(inarray, c_void_p),
                                          self._get_ctypes_data(outarray, c_void_p),
                                          scaling)
        else:
            raise heffte_input_error("backward() called with wrong outarray.dtype, use float32, float64, complex64, or complex128")

    def backward_buffered(self, inarray, outarray, workspace, scaling = scale.none):

        if scaling not in (0, 1, 2):
            raise heffte_input_error("backward() called with invalid scaling")

        if inarray.size != self.size_inbox() or outarray.size != self.size_outbox():
            raise heffte_input_error("backward() called with invalid array size")

        if self.use_r2c:
            if outarray.dtype == np.complex64 or outarray.dtype == np.complex128:
                raise heffte_input_error("backward() called with r2c can use only real outarray types")

        if (outarray.dtype == np.float32 or outarray.dtype == np.complex64) and workspace.dtype != np.complex64:
            raise heffte_input_error("backward() called with incorrect workspace type, must be complex64")
        if (outarray.dtype == np.float64 or outarray.dtype == np.complex128) and workspace.dtype != np.complex128:
            raise heffte_input_error("backward() called with incorrect workspace type, must be complex128")

        if outarray.dtype == np.float32:
            if inarray.dtype != np.complex64:
                raise heffte_input_error("backward() with outarray.dtype == float32 needs inarray.dtype == complex64")
            libheffte.heffte_backward_c2s_buffered(self.plan,
                                                   self._get_ctypes_data(inarray, c_void_p),
                                                   self._get_ctypes_data(outarray, POINTER(c_float)),
                                                   self._get_ctypes_data(workspace, c_void_p),
                                                   scaling)
        elif outarray.dtype == np.complex64:
            if inarray.dtype != np.complex64:
                raise heffte_input_error("backward() with outarray.dtype == complex64 needs inarray.dtype == complex64")
            libheffte.heffte_backward_c2c_buffered(self.plan,
                                                   self._get_ctypes_data(inarray, c_void_p),
                                                   self._get_ctypes_data(outarray, c_void_p),
                                                   self._get_ctypes_data(workspace, c_void_p),
                                                   scaling)
        elif outarray.dtype == np.float64:
            if inarray.dtype != np.complex128:
                raise heffte_input_error("backward() with outarray.dtype == float64 needs inarray.dtype == complex128")
            libheffte.heffte_backward_z2d_buffered(self.plan,
                                                   self._get_ctypes_data(inarray, c_void_p),
                                                   self._get_ctypes_data(outarray, POINTER(c_double)),
                                                   self._get_ctypes_data(workspace, c_void_p),
                                                   scaling)
        elif outarray.dtype == np.complex128:
            if inarray.dtype != np.complex128:
                raise heffte_input_error("backward() with outarray.dtype == complex128 needs inarray.dtype == complex128")
            libheffte.heffte_backward_z2z_buffered(self.plan,
                                                   self._get_ctypes_data(inarray, c_void_p),
                                                   self._get_ctypes_data(outarray, c_void_p),
                                                   self._get_ctypes_data(workspace, c_void_p),
                                                   scaling)
        else:
            raise heffte_input_error("backward() called with wrong outarray.dtype, use float32, float64, complex64, or complex128")

    def size_inbox(self):
        '''
        Returns the size of the inbox used in the constructor.
        '''
        return int(libheffte.heffte_size_inbox(self.plan))

    def size_outbox(self):
        '''
        Returns the size of the outbox used in the constructor.
        '''
        return int(libheffte.heffte_size_outbox(self.plan))

    def size_workspace(self):
        '''
        Returns the size of the workspace that will be used in computation.
        '''
        return int(libheffte.heffte_size_workspace(self.plan))


# Create a processor grid using the minimum surface algorithm
def proc_setup(world, num_procs):
    assert(world.count() > 0)
    all_indexes = world.size

    best_grid = [1, 1, num_procs]
    surface = lambda x: np.sum( all_indexes/x * np.roll(all_indexes/x,-1) )
    best_surface = surface([1, 1, num_procs])

    for i in np.arange(1, num_procs+1):
        if (num_procs % i == 0):
            remainder = num_procs / i
            for j in np.arange(1, remainder+1):
                if (remainder % j == 0):
                    candidate_grid = [i, j, remainder / j]
                    candidate_surface = surface(candidate_grid)
                    if (candidate_surface < best_surface):
                        best_surface = candidate_surface
                        best_grid    = candidate_grid

    best_grid = np.array(best_grid, dtype=np.int32)
    assert(np.prod(best_grid) == num_procs)
    return best_grid

def split_world(world, proc_grid):
    fast = lambda i : world.low[0] + i * (world.size[0] / proc_grid[0]) + min(i, (world.size[0] % proc_grid[0]))
    mid  = lambda i : world.low[1] + i * (world.size[1] / proc_grid[1]) + min(i, (world.size[1] % proc_grid[1]))
    slow = lambda i : world.low[2] + i * (world.size[2] / proc_grid[2]) + min(i, (world.size[2] % proc_grid[2]))

    result = []
    for k in np.arange(proc_grid[2]):
        for j in np.arange(proc_grid[1]):
            for i in np.arange(proc_grid[0]):
                result.append(box3d([fast(i), mid(j), slow(k)], [fast(i+1)-1, mid(j+1)-1, slow(k+1)-1], world.order))
    return result
