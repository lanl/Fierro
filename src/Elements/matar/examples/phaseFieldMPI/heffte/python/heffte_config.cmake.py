
libheffte_path = "@heffte_library_path@"

enable_fftw = ("@Heffte_ENABLE_FFTW@" == "ON")
enable_mkl  = ("@Heffte_ENABLE_MKL@" == "ON")
enable_cuda = ("@Heffte_ENABLE_CUDA@" == "ON")
enable_rocm = ("@Heffte_ENABLE_ROCM@" == "ON")

__version__ = "@PROJECT_VERSION@"
