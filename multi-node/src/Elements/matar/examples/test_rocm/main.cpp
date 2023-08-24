#include "matar.h"
#include "SomeClass.h"

int main(int argc, char* argv[])
{
    Kokkos::initialize(argc, argv);
    {            
        FOR_ALL(i, 0, 10, {
            SomeClass s;
            s.some_func();
        });
    }
    Kokkos::finalize();
    return 0;
}
