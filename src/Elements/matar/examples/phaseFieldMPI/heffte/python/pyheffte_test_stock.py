
import heffte
import numpy as np
import mpi4py

def make_reference(num_entries, dtype, scale):
    reference = np.zeros((num_entries,), dtype)
    reference[0] = -512.0
    if scale == heffte.scale.symmetric:
        reference /= np.sqrt(float(2 * num_entries))
    elif scale == heffte.scale.full:
        reference /= float(2 * num_entries)
    return reference

comm = mpi4py.MPI.COMM_WORLD
me = comm.Get_rank()

assert comm.Get_size() == 2

box = (heffte.box3d([0, 0, 0], [3, 3, 1])
       if me == 0 else
       heffte.box3d([0, 0, 2], [3, 3, 3]))

fft = heffte.fft3d(heffte.backend.stock,
                   box, box, comm)

assert fft.size_inbox() == 32
assert fft.size_outbox() == 32

test_types = [[np.float32,    np.complex64,  heffte.scale.none, heffte.scale.full, 1.E-3],
              [np.complex64,  np.complex64,  heffte.scale.symmetric, heffte.scale.symmetric, 1.E-3],
              [np.float64,    np.complex128, heffte.scale.full, heffte.scale.none, 1.E-11],
              [np.complex128, np.complex128, heffte.scale.none, heffte.scale.full, 1.E-11],
              ]

for tt in test_types:

    in_array = np.array(range(fft.size_inbox()), tt[0])
    out_array = np.empty(fft.size_outbox(), tt[1])
    workspace = np.empty(fft.size_workspace(), tt[1])

    fft.forward_buffered(in_array, out_array, workspace, tt[2])

    reference = make_reference(fft.size_outbox(), tt[1], tt[2]) # num_entries, type, scale
    if me == 1:
        assert np.max(np.abs(reference - out_array)) < tt[4]

    in_array = np.zeros((fft.size_inbox(),), tt[0])
    fft.backward_buffered(out_array, in_array, workspace, tt[3])

    assert np.max(np.abs(in_array - np.array(range(fft.size_inbox()), tt[0]))) < tt[4]
    if me == 0:
        print("  pass   ",tt[0],"  ->  ",tt[1])

