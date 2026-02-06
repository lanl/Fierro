#pragma once

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif


EXTERNC
void vumat_(
        int* nblock, int* ndir, int* nshr, int* nstatev, int* nfieldv, int* nprops, int* lanneal,
        double* stepTime, double* totalTime, double* dt, const char* cmname, double* coordMp, double* charLength,
        double* props, double* density, double* strainInc, double* relSpinInc,
        double* tempOld, double* stretchOld, double* defgradOld, double* fieldOld,
        double* stressOld, double* stateOld, double* enerInternOld, double* enerInelasOld,
        double* tempNew, double* stretchNew, double* defgradNew, double* fieldNew,
        double* stressNew, double* stateNew, double* enerInternNew, double* enerInelasNew );