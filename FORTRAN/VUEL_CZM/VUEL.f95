subroutine ComputeMassMatrix (                  &
    nblock,rho,coords,ndofel,amass              &
    )

    !   Function:
    !   calculate the mass matrix of cohesive element
    !
    !   Notices:
    !   The area of a triangle can be calculated by the outer product of two vectors 
    !
    !   Parameters:
    !   nblock  : the total number of units for which vuel is called this time
    !   rho : ρ density
    !   coords  : an array containing node coordinates, coords(a,b,c) : the c coordinate of the b node of the a element
    !


end subroutine ComputeMassMatrix

! ===================================== Obtain parameters from ABAQUS main program ========================
subroutine GetParameter (                       &
    props,constant,                             &
    rkn,rkt,                                    &
    strengthN,strengthT,strength,               &               
    Gnc,Gtc,                                    &
    rho,exp,rknt,                               & 
    delta10,delta20,delta1f,delta2f,gaussPara   &                          
    )
    
    !! 
    dimension props(8),rknt(2),strength(2),gaussPara(2,4)



end subroutine GetParameter

! ====================================== VUEL subroutine ====================================================
subroutine vuel (                               &
    nblock,rhs,amass,dtimeStable,svars,nsvars,  &
    energy,                                     &
    nnode,ndofel,props,nprops,jprops,njprops,   &
    coords,mcrd,u,du,v,a,                       &
    jtype,jElem,                                &
    time,period,dtimeCur,dtimePrev,kstep,kinc,  &
    lflags,                                     &
    dMassScaleFactor,                           &
    predef,npredef,                             &
    jdltyp, adlmag                              &
    )

    ! nblock,rhs,amass,dtimeStable,svars,nsvars,  
    ! energy,                                     
    ! nnode,ndofel,props,nprops,jprops,njprops,   
    ! coords,mcrd,u,du,v,a,                       
    ! jtype,jElem,                                
    ! time,period,dtimeCur,dtimePrev,kstep,kinc,  
    ! lflags,                                     
    ! dMassScaleFactor,                           
    ! predef,npredef,                             
    ! jdltyp, adlmag 

    !! operational code keys
    parameter (                                 &
        jMassCalc            = 1,               &   ! define the mass matrix 'amass' in the beginning of the analysis            
        jIntForceAndDtStable = 2,               &   ! define the element internal force and the stable time increment
        jExternForce         = 3                &   ! define the distributed load eddect on the ecternal force associated with the element
    )

    !! flag indices
    parameter (                                 &
        iProcedure = 1,                         &
        iNlgeom    = 2,                         &
        iOpCode    = 3,                         &
        nFlags     = 3                          &
    )

    !! procedure flags
    parameter (                                 &
        jDynExplicit      =   17,               &   ! direct integration explicit dynamic analysis
        jThermoMechanical =   74                &   ! transient fully coupled thermal-stress analysis
    )

    !! energy array indices
    parameter (                                 &
        iElPd               = 1,                &
        iElCd               = 2,                &
        iElIe               = 3,                &
        iElTs               = 4,                &
        iElDd               = 5,                &
        iElBv               = 6,                &
        iElDe               = 7,                &
        iElHe               = 8,                &
        iUnused             = 9,                &
        iElTh               = 10,               &
        iElDmd              = 11,               &
        iElDc               = 12,               &
        nElEnergy           = 12                &
    )

    !! predefined variables indices
    parameter (                                 &
        iPredValueNew = 1,                      &
        iPredValueOld = 2,                      &
        nPred         = 2                       &
    )

    !! time indices
    parameter (                                 &
        iStepTime  = 1,                         &
        iTotalTime = 2,                         &
        nTime      = 2                          &
    )

! #BUG 这里的数字定义不完整，或者有些用不到
    dimension                                   &
        rhs(nblock,ndofel),                     &   
        amass(nblock,ndofel,ndofel),            &
        dtimeStable(nblock),                    &
        svars(nblock,nsvars),                   &
        energy(nblock,nElEnergy),               &
        props(nprops),                          &
        jprops(njprops),                        &
        jElem(nblock),                          &
        time(nTime),                            &
        lflags(nFlags),                         &
        coords(nblock,nnode,mcrd),              &
        u(nblock,ndofel),                       &
        du(nblock,ndofel),                      &
        v(nblock,ndofel),                       &
        a(nblock, ndofel),                      &
        dMassScaleFactor(nblock),               &
        predef(nblock,nnode,npredef,nPred),     &
        adlmag(nblock)                          


!#BUG 这里的jtyoe可能不是2，师兄的单元中是2，暂用2代替
    if ( jtype == 2 .and. lflags(iProcedure) == jDynExplicit ) then
        !! jtype == 2                               : determines if the element type matches
        !! lflags(iProcedure) == jDynExplicit       : direct integration explicit dynamic analysis

        !! parameter  assignment
!#TODO 完成参数赋值操作
        call GetParameter(                              &
            props,constant,rkn,rkt,strengthN,           &
            trengthT,strength,Gnc,Gtc,rho,exp,rknt,     &
            delta10,delta20,delta1f,delta2f,gaussPara   &
            )

        !! write mass matrix
        if ( lflags(iOpCode) == jMassCalc ) then
            !! compute mass matrix of the cohesive element

        else if ( lflags(iOpCode) == jIntForceAndDtStable ) then
            !! write the internal force vector, define the upper limit of the time increment and solve the dependent state variable

        end if




    else if ( jtype == 2 .and. lflags(iProcedure) == jThermoMechanical ) then
        !! jtype == 2                               : determines if the element type matches
        !! lflags(iProcedure) == jThermoMechanical  : transient fully coupled thermal-stress analysis
    
    end if

end subroutine vuel




