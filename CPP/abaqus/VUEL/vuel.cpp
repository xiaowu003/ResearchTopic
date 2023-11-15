#include <aba_for_c.h>              //abaqus fortran 到c接口的转换宏

extern "C"

void FOR_NAME(vuel) {
    nblock,

    /* variables to be defined */
    rhs,
    amass,
    dtimeStable,
    svars,
    nsvars,
    energy,

    /*  */
    nnode,
    ndofel,
    props,
    nprops,
    jprops,
    njprops,
    coords,
    mcrd,
    u,
    du,
    v,
    a,
    jtype,
    jElem,
    time,
    eriod,
    dtimeCur,
    dtimePrev,
    kstep,
    kinc,
    lflags,
    dMassScaleFactor, 
    predef,
    npredef,
    jdltyp, 
    adlmag
}