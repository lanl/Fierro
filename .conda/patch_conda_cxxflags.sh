arch_flags=(-fno-tree-vectorize)
if [ "$PLATFORM" == "linux-64" ]; then
    arch_flags+=(-march=x86-64 -mtune=generic)
    echo "Compiling with generic x86_64 instruction set."
    # VECTOR_ARCH_FLAGS is used by FindVector.cmake
    export VECTOR_ARCH_FLAGS=" -march=x86-64 -mtune=generic "
fi
if [ "$PLATFORM" == "linux-ppc64le" ]; then
    arch_flags+=(-mtune=powerpc64)
    echo "Tuning for powerpc64 CPUs"
    export VECTOR_ARCH_FLAGS=" -mcpu=powerpc64le "
fi
if [ "$PLATFORM" == "linux-aarch64" ]; then
    arch_flags+=(-mtune=generic)
    echo "Tuning for generic aarch64 CPUs"
    export VECTOR_ARCH_FLAGS=" -mtune=generic "
fi

# Some apple architectures. 
# Only applies for clang.
if [ "$PLATFORM" == "osx-arm64" ]; then
    arch_flags+=(--target=arm64-apple-darwin)
    echo "Tuning for generic aarch64 CPUs"
    export VECTOR_ARCH_FLAGS=" -mtune=generic "
fi
if [ "$PLATFORM" == "osx-64" ]; then
    arch_flags+=(--target=x86_64-apple-darwin)
    echo "Tuning for generic aarch64 CPUs"
    export VECTOR_ARCH_FLAGS=" -mtune=generic "
fi

PATCHED_CXXFLAGS=()
for arg in $CXXFLAGS
do
    case $arg in 
        -march* | -mtune* | -mcpu*)
            echo "Removing architecture CXXFLAG: $arg"
        ;;
        -fvisibility-inlines-hidden | -ftree-vectorize)
            echo "Removing CXXFLAG $arg"
        ;;
        *)
            PATCHED_CXXFLAGS+=($arg)
        ;;
    esac
done

PATCHED_CXXFLAGS+=(${arch_flags[@]})
PATCHED_CXXFLAGS=${PATCHED_CXXFLAGS[@]}
export PATCHED_CXXFLAGS