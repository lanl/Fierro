# Heffte forwards a link to sysroot libpthread from fftw 
# This is annoying for anaconda, so we just patch it out.
sed -i.backup "s|;/home/[^;]*/sysroot/usr/lib/libpthread.so||g" ${PREFIX}/lib/cmake/Heffte/HeffteConfig.cmake
rm ${PREFIX}/lib/cmake/Heffte/HeffteConfig.cmake.backup