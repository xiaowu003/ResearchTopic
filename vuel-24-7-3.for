!============================================Some system constant====================================================
module system_constant
    implicit none
    integer, parameter :: RK = selected_real_kind(p = 15)
end module system_constant

!============================================Some math functions=========================================================
module math_mod
    use system_constant
    contains
    !============================================vector cross================================================
    subroutine vector_cross(vec1, vec2, cross)
        implicit none
        real(kind=RK), intent(in) :: vec1(3)
        real(kind=RK), intent(in) :: vec2(3)
        real(kind=RK), intent(out) :: cross(3)

        cross(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
        cross(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
        cross(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)
    end subroutine vector_cross


    !============================================vector mod================================================
    function vector_mod(vec) result(vec_mod)
        implicit none
        real(kind=RK), intent(in) :: vec(3)
        real(kind=RK) :: vec_mod

        vec_mod = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
    end function vector_mod

    !============================================dot product================================================
    function dot_product_1(vec1, vec2) result(dot)
        implicit none
        real(kind=RK), intent(in) :: vec1(3)
        real(kind=RK), intent(in) :: vec2(3)
        real(kind=RK) :: dot

        dot = vec1(1) * vec2(1) + vec1(2) * vec2(2) + vec1(3) * vec2(3)
    end function dot_product_1

    function get_unit_vector(vec) result(unit_vec)
        implicit none
        real(kind=RK), intent(in) :: vec(3)
        real(kind=RK) :: unit_vec(3)

        unit_vec = vec / vector_mod(vec)
        !write(*,*) "vecmod= ", vector_mod(vec)
    end function get_unit_vector

end module math_mod

!============================================Some parameters=========================================================
module parameter_mod
    use system_constant
    implicit none
    integer, save :: num_of_ele
    integer, save :: node_num_of_an_ele
    integer, save :: half_node_num_of_an_ele
    integer, save :: dof_num_of_a_node
    integer, save :: dof_num_of_an_ele
    integer, save :: dimmension
    real(kind=RK), save :: rho

    real(kind=RK), save :: pi = 3.14159265358979323846d0
    real(kind=RK), save :: thickness = 1.0d0
    !real(kind=RK), save :: Tet_gaussPoint(3, 2) = ( (0.0d0, 0.0d0), (1.0d0, 0.0d0), (0.0d0, 1.0d0) )
    real(kind=RK), save :: Tet_gaussPoint(2, 3) = reshape((/ 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 1.0d0 /), shape(Tet_gaussPoint))
end module parameter_mod


!============================================Some Cohesive functions=========================================================
module cohesive_element_mod
    use system_constant
    use math_mod
    use parameter_mod
    contains

    !============================================shape function q3================================================
    subroutine shape_function_q3(Tet_gaussPoint1, Tet_gaussPoint2, shape)
        implicit none
        real(kind=RK), intent(in) :: Tet_gaussPoint1, Tet_gaussPoint2
        real(kind=RK), intent(out) :: shape(3)

        shape(1) = 1.0d0 - Tet_gaussPoint1 - Tet_gaussPoint2
        shape(2) = Tet_gaussPoint1
        shape(3) = Tet_gaussPoint2

    end subroutine shape_function_q3

    !============================================shape function derivate================================================
    subroutine shape_function_derivate(Tet_gaussPoint, dshape1, dshape2)
        implicit none
        real(kind=RK), intent(in) :: Tet_gaussPoint(2)
        real(kind=RK), intent(out) :: dshape1(3)
        real(kind=RK), intent(out) :: dshape2(3)

        dshape1(1) = -1.0d0
        dshape1(2) = 1.0d0
        dshape1(3) = 0.0d0

        dshape2(1) = -1.0d0
        dshape2(2) = 0.0d0
        dshape2(3) = 1.0d0

    end subroutine shape_function_derivate

    !===========================================distribution of force=========================================
    subroutine distribute_force(i, j, force_vec, rhs_temp, Tet_gaussPoint)
        implicit none
        integer, intent(in) :: i, j
        real(kind=RK), intent(in) :: force_vec(3)
        real(kind=RK), intent(out) :: rhs_temp(dof_num_of_an_ele)
        real(kind=RK), intent(in) :: Tet_gaussPoint(2)

        integer :: k
        real(kind=RK) :: shape(3)
        real(kind=RK) :: rhs_temp_reshape(dof_num_of_a_node, node_num_of_an_ele)        !! here may be (3, 6) or (4, 6)
        !!write(*,*) "here 11"

        !rhs_temp_reshape = reshape(rhs, (/dof_num_of_a_node, node_num_of_an_ele/))     !! here may be (3, 6) or (4, 6)

        rhs_temp_reshape = 0.0d0        ! initialize the rhs_temp_reshape

        call shape_function_q3(Tet_gaussPoint(1), Tet_gaussPoint(2), shape)

        !write(*,*) "Tet_gaussPoint= ", Tet_gaussPoint, " shape= ", shape, " force_vec= ", force_vec
        
        !!write(*,*) "here 12"
        
        do k = 1, half_node_num_of_an_ele
            !!write(*,*) "k= ", k, " half_node_num_of_an_ele", half_node_num_of_an_ele, " force_vec= ", force_vec, " shape(k)= ", shape(k)
            rhs_temp_reshape(1:3, (k + half_node_num_of_an_ele)) = rhs_temp_reshape(1:3, (k + half_node_num_of_an_ele)) + shape(k) * force_vec
        end do
        
        !!write(*,*) "here 13"
        
        !write(*,*) "rhs_temp_reshape= ", rhs_temp_reshape
        !write(*,*) "before rhs_temp= ", rhs_temp
        rhs_temp_reshape(1:3, 1:half_node_num_of_an_ele) = -rhs_temp_reshape(1:3, (half_node_num_of_an_ele + 1) : node_num_of_an_ele) 
        rhs_temp = rhs_temp + reshape(rhs_temp_reshape, (/dof_num_of_an_ele/))
        !write(*,*) "after rhs_temp= ", rhs_temp
        !!write(*,*) "here 14"
    end subroutine distribute_force

    !============================================Friction_COH================================================
    subroutine Friction_COH( shear_seperation_old, sigm_shear_old, delt_n, pk_n, pk_s, FMU_c, sigm_n, sigm_shear, unit_vec_t, unit_vec_s, shear_s, shear_t )
        implicit none
        real(kind=RK), intent(in) :: shear_seperation_old(3)
        real(kind=RK), intent(in) :: sigm_shear_old(3)
        real(kind=RK), intent(in) :: delt_n
        !real(kind=RK), intent(in) :: delt_shear_inc
        real(kind=RK), intent(in) :: pk_n
        real(kind=RK), intent(in) :: pk_s
        real(kind=RK), intent(in) :: FMU_c
        !real(kind=RK), intent(in) :: delt_shear_inc_vec(3)
        real(kind=RK), intent(out) :: sigm_n
        real(kind=RK), intent(out) :: sigm_shear
        real(kind=RK), intent(in) :: unit_vec_t(3), unit_vec_s(3), shear_s, shear_t
        

        real(kind=RK) :: friction_max
        real(kind=RK) :: friction_temp(2), friction_temp_mod
        real(kind=RK) :: delt_shear_s_t(2)

        real(kind=RK) :: shear_s_old, shear_t_old, sigm_shear_s_old, sigm_shear_t_old, sigm_old_temp(2)
        real(kind=RK) :: delt_shear_s, delt_shear_t

        shear_s_old = dot_product_1(shear_seperation_old, unit_vec_s)
        shear_t_old = dot_product_1(shear_seperation_old, unit_vec_t)
        delt_shear_s = shear_s - shear_s_old
        delt_shear_t = shear_t - shear_t_old
        sigm_shear_s_old = dot_product_1(sigm_shear_old, unit_vec_s)
        sigm_shear_t_old = dot_product_1(sigm_shear_old, unit_vec_t)

        sigm_n = pk_n * delt_n
        friction_max = abs(FMU_c * sigm_n)
        
        sigm_old_temp = (/ sigm_shear_s_old, sigm_shear_t_old /)
        delt_shear_s_t = (/ delt_shear_s, delt_shear_t /)
        friction_temp = sigm_old_temp + delt_shear_s_t * pk_s
        friction_temp_mod = norm2(friction_temp)

        if (friction_temp_mod <= friction_max) then
            sigm_shear = friction_temp_mod
        else
            sigm_shear = friction_max
        end if
    end subroutine Friction_COH

    !============================================BK_M_C================================================
    subroutine BK_M_C(delt_n0, delt_s0, delt_t0,                                    &         !三个方向的初始断裂位移
                      delt_n, delt_shear, penaltyK_all, GIC, GIIC, ETA, delt_sc,    &         !位移、罚刚度、断裂能及系数
                      delt_sf, beta, pk_n, penaltyK_shear, SIGM_NMAX, FMU_i, C,     &        
                      SDEG)
        real(kind = RK) :: delt_n0,delt_s0,delt_t0,                                &       !定义变量为双精度浮点数   
                        delt_n, delt_shear,penaltyK_all,GIC,GIIC,ETA,delt_sc,   &  
                        delt_sf,beta, pk_n,penaltyK_shear,SIGM_NMAX,FMU_i,C,    &
                        SDEG,SDEG_NEW,delt_m,delt_m0, mix_rate,delt_mf,beta_0
        !C-------------------------------------------------------------
        delt_m = SQRT((DIM(delt_n, 0.0d0))**2 + delt_shear**2)  !当前分离量 

        if (delt_n > 0.0d0) then                                                   !初始断裂分离量
           beta    = ABS(delt_shear / delt_n )                                    !混合比 beta>=0
           beta_0  = pk_n * (C - SIGM_NMAX * FMU_i) / (penaltyK_shear * SIGM_NMAX)       !计算delt_m0的分界点
        
           !C--------------------------        
           if ( beta <= beta_0 ) then
                delt_m0 = SIGM_NMAX * SQRT(1 + beta**2) / pk_n
           else
                delt_m0 = C * SQRT(1 + beta**2) / (penaltyK_shear * beta + pk_n * FMU_i)
           end if

           !C--------------------------
           mix_rate = penaltyK_shear * beta**2 / (pk_n + penaltyK_shear * beta**2)
           delt_mf = 2.0d0 / penaltyK_all / delt_m0 * (GIC + (GIIC - GIC) * (mix_rate)**ETA)
        else                                                  !不然则为压剪模式
           delt_m0 = delt_s0
           delt_mf = delt_sf
        end if
        !C-------------------------------------------------------------
        if (delt_m > delt_m0) then
            SDEG_NEW =( (delt_mf - delt_sc) * (delt_m - delt_m0) ) /( (delt_mf - delt_m0) * (delt_m - delt_sc) )    
        else
             SDEG_NEW = 0.0d0                                    !取零，不影响实际sdeg取值
        end if   

        SDEG  = MAX( SDEG ,SDEG_NEW ) 
        return 
    end subroutine BK_M_C

    !============================================vuel_vumat_sub================================================
    subroutine vuel_vumat_constitutive( i, j, time, props, nprops, svars, nsvars, normal_seperation, tangential_seperation, delt_n_inc, sigm_n, sigm_shear, shear_s, shear_t, unit_vec_t, unit_vec_s )
        implicit none
        integer, intent(in) :: i, j, nprops, nsvars
        real(kind=RK), intent(in) :: time(2)
        real(kind=RK), intent(in) :: props(nprops)
        real(kind=RK), intent(in) :: normal_seperation
        real(kind=RK), intent(in) :: tangential_seperation
        real(kind=RK), intent(in) :: delt_n_inc
        !real(kind=RK), intent(in) :: delt_shear_inc
        !real(kind=RK), intent(in) :: delt_shear_inc_vec(3)
        !real(kind=RK), intent(in) :: unit_tangential_vec(3)

        real(kind=RK), intent(in) :: shear_t, shear_s

        real(kind=RK), intent(inout) :: svars(nsvars)
        real(kind=RK), intent(out) :: sigm_n
        real(kind=RK), intent(out) :: sigm_shear
        real(kind=RK), intent(in) :: unit_vec_t(3), unit_vec_s(3)

        real(kind=RK) :: pk_n, pk_s, pk_t
        real(kind=RK) :: epsil_n, epsil_s, epsil_t
        real(kind=RK) :: delt_n, delt_s, delt_t, delt_shear
        real(kind=RK) :: SIGM_fr, SIGM_fi, SIGM_TMAX, SIGM_SMAX
        real(kind=RK) :: delt_n0, delt_s0, delt_t0, delt_sc
        real(kind=RK) :: epsil_n0, epsil_s0, epsil_t0
        real(kind=RK) :: delt_nf, delt_sf, delt_tf
        real(kind=RK) :: DMICRT, SDEG, fracture_mode
        real(kind=RK) :: penaltyK_all, penaltyK_shear
        real(kind=RK) :: pk_nd, pk_sd
        real(kind=RK) :: delt_shearInc
        real(kind=RK) :: beta
        real(kind=RK) :: delete_flag
        real(kind=RK) :: delt_shear_D
        real(kind=RK) :: pc
        real(kind=RK) :: EN, ES, SIGM_NMAX, C, PHI, GIC, GIIC, ETA, D_MODE, FMU_c, ET, FMU_i, FMU_r
        real(kind=RK) :: DMICRT_n, DMICRT_shear

        !real(kind=RK) :: delt_shear_old
        real(kind=RK) :: sigm_shear_old(3), shear_seperation_old(3)
        !real(kind=RK) :: node_normal_vec(3)
        !real(kind=RK) :: node_tangential_vec(3)
        !real(kind=RK) :: unit_shear_vec_old(3)

        !real(kind=RK) :: delt_shear_s, delt_shear_t
        !real(kind=RK) :: sigm_shear_s, sigm_shear_t
        !real(kind=RK) :: sigm_shear_s_old, sigm_shear_t_old

        real(kind=RK) :: shear_vec_old(3)

        !real(kind=RK), parameter :: pi = 3.14159265358979323D0
        real(kind=RK), parameter :: crit_lens = 1.18D0 
        real(kind=RK), parameter :: crit_lenn = 1.1D0

        ! ----------------------------提取材料参数
        EN      = props(1)      !杨氏模量
        ES      = props(2)      !切变模量
        SIGM_NMAX = props(3)    !界面抗拉强度
        C       = props(4)      !内聚力(当PHI=0时为界面抗剪强度）
        PHI     = props(5)      !内摩擦角/deg
        GIC     = props(6)      !I型断裂能
        GIIC    = props(7)      !II型断裂能
        ETA     = props(8)      !混合模式系数
        D_MODE  = props(9)      !初始破坏准则标志 M-C=0  Quads=2 
        FMU_c   = props(10)     !完全破坏后的 摩擦系数

        ET      = ES                !切变模量
        FMU_i   = tan(PHI/180*PI)   !摩擦角等效摩擦系数--初始破坏处
        FMU_r   = FMU_i             !摩擦角等效摩擦系数--完全破坏处

        pk_n   =    EN                              !罚刚度
        pk_s   =    ES   
        pk_t   =    pk_s

        !write(*,*) "pk_n= ", pk_n, " pk_s= ", pk_s, " c = ", C, " SIGM_NMAX= ", SIGM_NMAX, " FMU_i= ", FMU_i, " FMU_c= ", FMU_c

        delt_n = normal_seperation              !位移
        delt_shear = tangential_seperation     

        !write(*,*) "delt_n= ", delt_n, " delt_shear= ", delt_shear
   
        !计算初始断裂位移和应变、完全破坏位移
        SIGM_fr   = pk_n * DIM(-delt_n, 0.0d0) * FMU_i                !残余摩擦应力--初始破坏处
        SIGM_fi   = pk_n * DIM(-delt_n, 0.0d0) * FMU_i    !!!初步修改,需进一步调试!!!!
        SIGM_TMAX = C + SIGM_fi                                !仅在压剪和纯剪模式下有意义   
        SIGM_SMAX = SIGM_TMAX                                  !剪切强度
        !C     
        delt_n0 = SIGM_NMAX / (pk_n)                            !初始断裂位移，不适用于拉剪混合模式
        delt_s0 = SIGM_SMAX / (pk_s)
        delt_t0 = SIGM_TMAX / (pk_t)
        delt_sc = SIGM_fr   / (pk_s)
        !C        
        epsil_n0 = SIGM_NMAX / EN                               !初始断裂应变
        epsil_s0 = SIGM_SMAX / ES
        epsil_t0 = SIGM_TMAX / ET
        !C      
        delt_nf = 2.d0 * GIC / SIGM_NMAX                           !完全破坏位移
        delt_sf = 2.d0 * GIIC / C + delt_sc
        delt_tf = delt_sf
        !C      
        DMICRT = svars((j - 1) * 12 + 2)                            !继承上一步变量
        SDEG   = svars((j - 1) * 12 + 3)                            
        fracture_mode = -1.0d0                                 !赋初值  -1表示未完全损伤

        !if (time(2) <= 8e-8) then
        !    write(*,*) "j = ", j, "(j - 1) * 11 + 2 = ", (j - 1) * 11 + 2, "svars() = ", svars((j - 1) * 12 + 2)
        !    write(*,*)"1 time: ",time(2), " delt_n: ", delt_n, " delt_shear: ",delt_shear," D: ",DMICRT," SDEG: ",SDEG
        !end if

        if ( SDEG >= 1.0d0 ) then
            !write(*,*) "1 time: ", time(2), " SDEG: ", SDEG, " DMICRT: ", DMICRT, " sigm_n: ", sigm_n, " delt_n: ", delt_n," sigm_shear: ", sigm_shear, " delt_shear: ", delt_shear
        end if

        !!write(*,*) "here 2"      
        !C-------------------------罚刚度混合--------------------------
        if (delt_shear == 0) then         !避免纯压缩时 分母为零 该值为NaN
            penaltyK_all = pk_n
        else
            penaltyK_all =(pk_n * (DIM(delt_n, 0.0d0))**2 + pk_s * delt_shear**2) / ((DIM(delt_n, 0.0d0))**2 + delt_shear**2)
        end if
    
        penaltyK_shear = pk_s                                 !先不混合，避免当纯压缩时分母为0

        !!write(*,*) "here 3"
        !C--------------------------------------------------------------------------            
        if ( SDEG >= 1.0D0 ) then
            SDEG = 1.0D0 
        !初始损伤判断   
        else if ( DMICRT < 1.0d0 .and. SDEG < 1.0d0) then                 !判断单元是否已经损伤
            !若未损伤，根据计算上一增量步的DMICRT
            !选着初始损伤判断准则 M-C 纯摩尔库仑
            SIGM_TMAX    = C - pk_n * delt_n * FMU_i                    !判断损伤时，coh必处于弹性阶段
            DMICRT_n     = pk_n * DIM(delt_n, 0.0d0) / SIGM_NMAX          !未损伤时，按弹性计算。
            DMICRT_shear = (penaltyK_shear * delt_shear) / SIGM_TMAX 
            DMICRT = MAX(DMICRT_n, DMICRT_shear)

            !判断上一增量步是否损伤  
            if (DMICRT >= 1.0d0) then                                !若损伤，根据混合模式计算sdeg值
                DMICRT = 1.0d0                                       !将dmicrt置为1，保证其在[0,1]范围内

                call BK_M_C(delt_n0,delt_s0,delt_t0,delt_n, delt_shear,penaltyK_all,GIC,GIIC,ETA,delt_sc, &
                            delt_sf, beta, pk_n,penaltyK_shear,SIGM_NMAX, FMU_i, C, SDEG )   
            else                                                    !若未损伤，SDEG=0
                SDEG = 0.0d0   
            end if   
            
            if ( SDEG >= 1.0d0 ) then
                write(*,*) "time: ", time(2), " SDEG > 1.0d0"
            end if
        !若损伤，根据混合模式计算SDEG值              
        else if ( DMICRT >= 1.0d0 .and.SDEG < 1.0d0) then
            !!write(*,*) "here 4"
            DMICRT = 1.0d0
            call  BK_M_C(delt_n0,delt_s0,delt_t0,delt_n, delt_shear,penaltyK_all,GIC,GIIC,ETA,delt_sc,   &
                         delt_sf, beta, pk_n,penaltyK_shear,SIGM_NMAX,FMU_i,C, SDEG)  

                         
            if ( SDEG >= 1.0d0 ) then
                write(*,*) "time: ", time(2), " SDEG > 1.0d0"
            end if


        end if
        !C-------------------------------------------------------------
        if (SDEG >= 1.0d0 )then  !如果完全损伤，则执行摩擦
            
            !!write(*,*) "here 5"

            svars((j - 1) * 12 + 1) = 0
            delete_flag = 1

            !计算断裂时的断裂类型，仅在完全破坏后的下一增量步计算，更新一次后不再计算
            if ( svars((j - 1) * 12 + 5) == -1.0d0) then     !fracture_mode = -1 表示未计算
                if ( delt_n <= 0  ) then 
                    fracture_mode = 90.0d0
                else if( delt_n > 0  ) then    
                    fracture_mode = atan(delt_shear / delt_n) / pi * 180
                end if
            else
                fracture_mode = svars((j - 1) * 12 + 5)
            end if
            !-----------------------
            if(delt_n < 0 ) then    
                if ( delt_shear >= crit_lens ) then 
                    sigm_n = 0.0d0        
                    sigm_shear = 0.0d0

                    !sigm_shear_s = 0.0d0
                    !sigm_shear_t = 0.0d0

                    delete_flag = 0

                else
                    !unit_shear_vec_old(1) = svars((j - 1) * 12 + 9)
                    !unit_shear_vec_old(2) = svars((j - 1) * 12 + 10)
                    !unit_shear_vec_old(3) = svars((j - 1) * 12 + 11)

                    !sigm_shear_old = svars((j - 1) * 12 + 8)                !继承上一步切向应力 #bug
                    !sigm_shear_s_old = svars((j - 1) * 12 + 9)
                    !sigm_shear_t_old = svars((j - 1) * 12 + 10)

                    !sigm_shear_old = (/ sigm_shear_s_old, sigm_shear_t_old /)

                    !delt_shear_s = shear_s - svars((j - 1) * 12 + 11)
                    !delt_shear_t = shear_t - svars((j - 1) * 12 + 12)

                    sigm_shear_old(1) = svars((j - 1) * 12 + 7)
                    sigm_shear_old(2) = svars((j - 1) * 12 + 8)
                    sigm_shear_old(3) = svars((j - 1) * 12 + 9)

                    shear_seperation_old(1) = svars((j - 1) * 12 + 10)
                    shear_seperation_old(2) = svars((j - 1) * 12 + 11)
                    shear_seperation_old(3) = svars((j - 1) * 12 + 12)

                    call Friction_COH( shear_seperation_old, sigm_shear_old, delt_n, pk_n, pk_s, FMU_c, sigm_n, sigm_shear, unit_vec_t, unit_vec_s, shear_s, shear_t )                     
                end if
            else
                !delt_n >= 0
                if (delt_n >= crit_lenn .OR. delt_shear >= crit_lens) then
                    delete_flag = 0
                end if
                sigm_n = 0.0d0        
                sigm_shear = 0.0d0

                !sigm_shear_s = 0.0d0
                !sigm_shear_t = 0.0d0

            end if

        else if ( SDEG < 1.0d0 ) then

            !!write(*,*) "here 6"

            svars((j - 1) * 12 + 1) = 1
            delete_flag = 1
            !----------------------------------------
            if (delt_n  >= 0.0d0 ) then                  !更新应力 对损伤后的卸载有效
              pk_nd = (1.0d0 - SDEG) * pk_n
              pc    = 0.0d0
            else if (delt_n  < 0.0d0 ) then 
              pk_nd = pk_n
              delt_shear_D = (delt_s0 - delt_sc) * (delt_sf - delt_sc) / ( (delt_sf - delt_sc) - SDEG * (delt_sf - delt_s0) ) + delt_sc
              pc = SDEG * delt_sc / delt_shear_D
            end if
            !----------------------------------------

            pk_sd = (1.0d0 - SDEG + pc) * pk_s

            !---------------------------------------
            !调整压力变化阶段切应力的变化
            shear_vec_old(1) = svars((j - 1) * 12 + 10)
            shear_vec_old(2) = svars((j - 1) * 12 + 11)
            shear_vec_old(3) = svars((j - 1) * 12 + 12)

            delt_shearInc   = tangential_seperation - vector_mod(shear_vec_old)              !切向总位移增量 这一步的总切向位移-上一步的总切向位移

            if ( delt_shearInc <= 0 .and. delt_n < 0 .and. delt_n_inc <= 0 ) then
                pk_sd = svars((j - 1) * 12 + 4)
            end if

            !---------------------------------------
            sigm_n  = pk_nd * delt_n
            sigm_shear  = pk_sd * delt_shear
            
            !sigm_shear_s = pk_sd * shear_s
            !sigm_shear_t = pk_sd * shear_t
        
        end if
        !!write(*,*) "here 7"

        svars((j - 1) * 12 + 2) = DMICRT
        svars((j - 1) * 12 + 3) = SDEG
        svars((j - 1) * 12 + 4) = pk_sd
        svars((j - 1) * 12 + 5) = fracture_mode         !记录断裂类型
        svars((j - 1) * 12 + 6) = delete_flag


        !svars((j - 1) * 12 + 7) = delt_shear 
        !svars((j - 1) * 12 + 8) = sigm_shear

        !svars((j - 1) * 12 + 9) = unit_tangential_vec(1)
        !svars((j - 1) * 12 + 10) = unit_tangential_vec(2)
        !svars((j - 1) * 12 + 11) = unit_tangential_vec(3)
        !svars((j - 1) * 12 + 9) = sigm_shear_s
        !svars((j - 1) * 12 + 10) = sigm_shear_t

        !svars((j - 1) * 12 + 11) = shear_s
        !svars((j - 1) * 12 + 12) = shear_t


        !if ( SDEG >= 1.0d0 ) then
        !    !write(*,*) "2 time: ", time(2), " SDEG: ", SDEG, " DMICRT: ", DMICRT, " sigm_n: ", sigm_n, " delt_n: ", delt_n," sigm_shear: ", sigm_shear, " delt_shear: ", delt_shear
        !    write(*,*) "time: ", time(2), " SDEG > 1.0d0"
        !end if

        !write(*,*) "time: ", time(1), " SDEG: ", SDEG, " DMICRT: ", DMICRT, " sigm_n: ", sigm_n, " sigm_shear: ", sigm_shear
        !!write(*,*) "here 8"

        !if (time(2) <= 8e-8) then
        !    write(*,*) "2 time: ", time(2), " delt_n: ", delt_n," delt_shear: ",delt_shear, " sigm_n: ", sigm_n, " sigm_shear: ",sigm_shear, " D: ",DMICRT," SDEG: ",SDEG
        !end if

    end subroutine vuel_vumat_constitutive

    !============================================seperation_component_and_vector_ele================================================
    subroutine seperation_component_and_vector_ele(mid_surf_coords, gap_of_nodes, unit_vec_t, unit_vec_s, unit_normal_vector, normal_seperation, tangential_seperation, unit_tangential_vec_ele, shear_s, shear_t, unit_tangential)
        implicit none
        real(kind=RK), intent(in) :: mid_surf_coords(3, 3), gap_of_nodes(3, 3)
        real(kind=RK), intent(in) :: unit_vec_t(3), unit_vec_s(3), unit_normal_vector(3)
        real(kind=RK), intent(inout) :: unit_tangential_vec_ele(3, 3), normal_seperation(3), tangential_seperation(3), shear_t(3), shear_s(3)

        integer :: i, k
        real(kind=RK) :: temp_tangential_seperation(3)
        real(kind=RK) :: edge(3, 2), unit_edge(3, 2)
        real(kind=RK) :: temp_s, temp_t, test_vec(3), edge_seperation(3)
        real(kind=RK), intent(inout) :: unit_tangential(3, 3)

        do i = 1, 3
            normal_seperation(i) = dot_product_1(gap_of_nodes(:, i), unit_normal_vector)
            unit_tangential_vec_ele(:, i) = 0.75d0 * ( gap_of_nodes(:, i) - normal_seperation(i) * unit_normal_vector )
            tangential_seperation(i) = vector_mod(unit_tangential_vec_ele(:, i))

            if (tangential_seperation(i) < 1.0e-18) then       !!#BUG
                unit_tangential(:, i) = 0.0d0
                tangential_seperation(i) = 0.0d0
            else
                unit_tangential(:, i) = unit_tangential_vec_ele(:, i) / tangential_seperation(i)
            end if

            shear_s(i) = dot_product_1( unit_tangential_vec_ele(:, i), unit_vec_s )
            shear_t(i) = dot_product_1( unit_tangential_vec_ele(:, i), unit_vec_t )

        end do

        !write(*,*) "normal_seperation= ", normal_seperation(1), normal_seperation(2), normal_seperation(3)
        !write(*,*) "tangential_seperation= ", tangential_seperation(1), tangential_seperation(2), tangential_seperation(3)

        !shear_s = 0.0d0
        !shear_t = 0.0d0

        !do i = 1, 3

            !do k = 0, 1
                !tangential_seperation(MOD(k + i, 3) + 1) = tangential_seperation(MOD(k + i, 3) + 1) + dot_product_1( unit_tangential_vec_ele(:, i), unit_tangential(:, MOD(k + i, 3) + 1) )
            !end do
            
            !do k = 0, 1
            !    该方法也不行，大变形错误
            !    edge(:, k + 1) = mid_surf_coords(:, i) - mid_surf_coords(:, MOD(k + i, 3) + 1)
            !    unit_edge(:, k + 1) = get_unit_vector(edge(:, k + 1))
            !    test_vec = unit_edge(:, k + 1) * dot_product_1( unit_tangential_vec_ele(:, i), unit_edge(:, k + 1) )
            !    temp_t = dot_product_1( test_vec , unit_vec_t )
            !    temp_s = dot_product_1( test_vec , unit_vec_s )
            !    shear_s(MOD(k + i, 3) + 1) = shear_s(MOD(k + i, 3) + 1) + temp_s
            !    shear_t(MOD(k + i, 3) + 1) = shear_t(MOD(k + i, 3) + 1) + temp_t
            !    !!shear_t(MOD(k + i, 3) + 1) = shear_t(MOD(k + i, 3) + 1) + dot_product_1( dot_product_1( unit_tangential_vec_ele(:, i), unit_edge(:, k + 1) ) *, unit_vec_t )
            !    !!shear_s(MOD(k + i, 3) + 1) = shear_s(MOD(k + i, 3) + 1) + dot_product_1( dot_product_1( unit_tangential_vec_ele(:, i), unit_edge(:, k + 1) ), unit_vec_s )
            !    !if (temp_s < 0.0d0 .or. temp_t < 0.0d0) then
            !    !    write(*,*) "temp_s= ", temp_s, " temp_t= ", temp_t
            !    !end if
            !    !shear_s(i) = ( unit_tangential_vec_ele(:, i) + unit_tangential_vec_ele(:, MOD(k + i, 3) + 1) ) * unit_edge(:, k + 1) 
            !end do

            !do k = 0, 0
                !该方法不行，还是提前断裂
            !   edge(:, k + 1) = mid_surf_coords(:, i) - mid_surf_coords(:, MOD(k + i, 3) + 1)
            !   unit_edge(:, k + 1) = get_unit_vector(edge(:, k +  
            !   edge_seperation = ( dot_product_1( unit_tangential_vec_ele(:, i), unit_edge(:, k + 1) ) + dot_product_1( unit_tangential_vec_ele(:, MOD(k + i, 3) + 1), unit_edge(:, k + 1) ) ) * 0.5d0 * unit_edge(:, k + 1)
            !   temp_t = dot_product_1( edge_seperation, unit_vec_t )
            !   temp_s = dot_product_1( edge_seperation, unit_vec_s )
            !   
            !   shear_s(i) = shear_s(i) + temp_s
            !   shear_s(MOD(k + i, 3) + 1) = shear_s(MOD(k + i, 3) + 1) + temp_s
            !   shear_t(i) = shear_t(i) + temp_t
            !   shear_t(MOD(k + i, 3) + 1) = shear_t(MOD(k + i, 3) + 1) + temp_t
            !end do


            !shear_t(i) = shear_t(i) +  dot_product_1( unit_tangential_vec_ele(:, i), unit_vec_t )
            !shear_s(i) = shear_s(i) +  dot_product_1( unit_tangential_vec_ele(:, i), unit_vec_s )
        !end do

        !write(*,*) "shear_s= ", shear_s(1), shear_s(2), shear_s(3)
        !write(*,*) "shear_t= ", shear_t(1), shear_t(2), shear_t(3)
        !write(*,*) "before tangentail_seperation= ", tangential_seperation(1), tangential_seperation(2), tangential_seperation(3)
        
        !do i = 1, 3
        !    unit_tangential_vec_ele(:, i) = unit_tangential(:, i) * tangential_seperation(i)
        !    tangential_seperation(i) = vector_mod(unit_tangential_vec_ele(:, i))
        !    if (tangential_seperation(i) < 1.0e-18) then       !!#BUG
        !        unit_tangential(:, i) = 0.0d0
        !        tangential_seperation(i) = 0.0d0
        !    else
        !        unit_tangential(:, i) = unit_tangential_vec_ele(:, i) / tangential_seperation(i)
        !    end if
        !end do


        
        !do i = 1, 3
        !    unit_tangential_vec_ele(:, i) = shear_t(i) * unit_vec_t + shear_s(i) * unit_vec_s
        !    tangential_seperation(i) = vector_mod( unit_tangential_vec_ele(:, i) )
        !    
        !    if (tangential_seperation(i) < 1.0e-18) then       !!#BUG
        !        unit_tangential_vec_ele(:, i) = 0.0d0
        !        tangential_seperation(i) = 0.0d0
        !    else
        !        unit_tangential_vec_ele(:, i) = unit_tangential_vec_ele(:, i) / tangential_seperation(i)
        !    end if
        !end do

        !write(*,*) "after tangentail_seperation= ", tangential_seperation(1), tangential_seperation(2), tangential_seperation(3)

    end subroutine seperation_component_and_vector_ele

    !============================================seperation component and vector================================================
    subroutine seperation_component_and_vector(j, gap_of_intergration_points, unit_normal_vector, normal_seperation, tangential_seperation, unit_tangential_vec)
        implicit none
        integer, intent(in) :: j
        real(kind=RK), intent(in) :: gap_of_intergration_points(3)
        real(kind=RK), intent(in) :: unit_normal_vector(3)
        real(kind=RK), intent(out) :: normal_seperation
        real(kind=RK), intent(out) :: tangential_seperation
        real(kind=RK), intent(out) :: unit_tangential_vec(3)

        real(kind=RK) :: node_normal_vec(3)
        real(kind=RK) :: node_tangential_vec(3)

        normal_seperation = dot_product_1(gap_of_intergration_points, unit_normal_vector)

        !write(*,*) "gap_of_intergration_points= "
        !write(*,*) gap_of_intergration_points
        !write(*,*) "unit_normal_vector= "
        !write(*,*) unit_normal_vector
        !write(*,*) "normal_seperation= ", normal_seperation

        node_normal_vec = normal_seperation * unit_normal_vector
        node_tangential_vec = gap_of_intergration_points - node_normal_vec
        tangential_seperation = vector_mod(node_tangential_vec)

        !write(*,*) "tangential_seperation= ", tangential_seperation

        if (tangential_seperation < 1.0e-20) then       !!#BUG
            unit_tangential_vec = 0.0d0
            tangential_seperation = 0.0d0
        else
            unit_tangential_vec = node_tangential_vec / tangential_seperation
        end if
    end subroutine seperation_component_and_vector

    !============================================gap of integral points================================================
    subroutine gap_of_integral_points(j, gap_of_nodes, gap_of_intergration_points, Tet_gaussPoint)
        implicit none
        integer, intent(in) :: j
        real(kind=RK), intent(in) :: gap_of_nodes(3,3)
        real(kind=RK), intent(out) :: gap_of_intergration_points(3)
        real(kind=RK), intent(in) :: Tet_gaussPoint(2)

        integer :: k
        real(kind=RK) :: shape(3)

        call shape_function_q3(Tet_gaussPoint(1), Tet_gaussPoint(2), shape)

        do k = 1, half_node_num_of_an_ele
            gap_of_intergration_points = gap_of_intergration_points + gap_of_nodes(:, k) * shape(k)
        end do

    end subroutine gap_of_integral_points

    !============================================calculate the area================================================
    subroutine calculate_area(j, mid_surf_coords, area)
        implicit none
        integer, intent(in) :: j
        real(kind=RK), intent(in) :: mid_surf_coords(3,3)
        real(kind=RK), intent(out) :: area

        real(kind=RK) :: center_of_gravity(3)
        real(kind=RK) :: cross_temp1(3), cross_temp2(3)
        real(kind=RK) :: temp_vec(3, 3)
        integer :: k

        center_of_gravity = 0.0d0

        do k = 1, half_node_num_of_an_ele
            center_of_gravity = center_of_gravity + mid_surf_coords(:, k)
        end do

        !write(*,*) "before center_of_gravity= ", center_of_gravity

        center_of_gravity = center_of_gravity / 3.0d0

        !write(*,*) "after center_of_gravity= ", center_of_gravity


        do k = 0, 1
            temp_vec(:, k + 1) = 0.5 * ( mid_surf_coords(:, MOD(k + j, 3) + 1) - mid_surf_coords(:, j) )     !!#bug, 第一列为指向下一节点的向量，第二列为指向上一节点的向量
        end do
        temp_vec(:, 3) = center_of_gravity - mid_surf_coords(:, j)   !! 第三列为指向重心的向量

        !write(*,*) "temp_vec1= ", temp_vec(:, 1)
        !write(*,*) "temp_vec2= ", temp_vec(:, 2)
        !write(*,*) "temp_vec3= ", temp_vec(:, 3)

        call vector_cross(temp_vec(:, 1), temp_vec(:, 3), cross_temp1)
        call vector_cross(temp_vec(:, 2), temp_vec(:, 3), cross_temp2)

        area = 0.5d0 * ( vector_mod( cross_temp1 )  + vector_mod( cross_temp2 ) ) 

    end subroutine calculate_area

    !=======================================calculate the local tangential basis vector======================================
    subroutine get_local_tangential_basis_vector(mid_surf_coords, unit_normal_vector, unit_vec_t, unit_vec_s)
        implicit none
        real(kind=RK), intent(in) :: mid_surf_coords(3,3)
        real(kind=RK), intent(in) :: unit_normal_vector(3)
        real(kind=RK), intent(out) :: unit_vec_t(3)
        real(kind=RK), intent(out) :: unit_vec_s(3)

        real(kind=RK) :: x_axis(3), z_axis(3), reference_vec(3)
        real(kind=RK) :: cos_theta, angle_radians, angle_degrees
        real(kind=RK) :: cross_temp(3)
        !real(kind=RK) :: temp_x_dot,
        

        x_axis = (/1.0d0, 0.0d0, 0.0d0/)
        z_axis = (/0.0d0, 0.0d0, 1.0d0/)
        reference_vec = x_axis

        !temp_x_dot = dot_product_1(unit_normal_vector, x_axis)
        cos_theta = dot_product_1(unit_normal_vector, x_axis)
        !cos_theta = temp_x_dot
        angle_radians = acos(cos_theta)
        angle_degrees = angle_radians * 180.0d0 / pi
        
        !call vector_cross(unit_normal_vector, x_axis, cross_temp)

        !if (vector_mod(cross_temp) < 0.0d0) then       #bug
        !    angle_degrees = 360.0d0 - angle_degrees
        !end if

        if (angle_degrees <= 0.1d0) then
            reference_vec = z_axis
        end if

        unit_vec_s = reference_vec - dot_product_1(unit_normal_vector, reference_vec) * unit_normal_vector
        unit_vec_s = unit_vec_s / vector_mod(unit_vec_s)

        call vector_cross(unit_normal_vector, unit_vec_s, unit_vec_t)

    end subroutine get_local_tangential_basis_vector

    !============================================calculate the unit normal vector===========================================
    subroutine get_unit_normal_vector(mid_surf_coords, unit_normal_vector)
        implicit none
        real(kind=RK), intent(in) :: mid_surf_coords(3,3)
        !real(kind=RK), intent(in) :: Tet_gaussPoint(2)
        real(kind=RK), intent(out) :: unit_normal_vector(3)

        integer :: k
        !real(kind=RK) :: dshape1(3), dshape2(3)
        real(kind=RK) :: temp_vec1(3), temp_vec2(3)
        real(kind=RK) :: mod_vector

        !call shape_function_derivate(Tet_gaussPoint, dshape1, dshape2)

        temp_vec1 = 0.0d0
        temp_vec2 = 0.0d0

        !do k = 1, half_node_num_of_an_ele
        !    temp_vec1 = temp_vec1 + mid_surf_coords(:, k) * dshape1(k)
        !    temp_vec2 = temp_vec2 + mid_surf_coords(:, k) * dshape2(k)
        !end do

        temp_vec1 = mid_surf_coords(:, 2) - mid_surf_coords(:, 1)
        temp_vec2 = mid_surf_coords(:, 3) - mid_surf_coords(:, 1)

        call vector_cross(temp_vec1, temp_vec2, unit_normal_vector)

        mod_vector = vector_mod(unit_normal_vector)
        if (abs(mod_vector) < 1.0d-15) then
            write(6,*) "The area of the integration point is negative"
            call xplb_exit
        end if

        unit_normal_vector = unit_normal_vector / mod_vector

    end subroutine get_unit_normal_vector

    !============================================calculate the shear increment================================================
    subroutine get_shear_inc(i, du, shear_inc)
        implicit none
        integer, intent(in) :: i
        real(kind=RK), intent(in) :: du(num_of_ele, dof_num_of_an_ele)
        real(kind=RK), intent(out) :: shear_inc(3, 3)

        integer :: j
        real(kind=RK) :: reshape_du(dof_num_of_a_node, node_num_of_an_ele)

        reshape_du = reshape(du(i,:), (/dof_num_of_a_node, node_num_of_an_ele/))

        do j = 1, half_node_num_of_an_ele
            shear_inc(:, j) = reshape_du(1:3, j+half_node_num_of_an_ele) - reshape_du(1:3, j)
        end do
    end subroutine get_shear_inc

    !============================================calculate the seperation and middle surface==============================
    subroutine seperation_and_middle_surface(i, u, coords, gap_of_nodes, mid_surf_coords)
        implicit none
        integer, intent(in) :: i
        real(kind=RK), intent(in) :: u(num_of_ele, dof_num_of_an_ele)
        real(kind=RK), intent(in) :: coords(num_of_ele, node_num_of_an_ele, dimmension)
        real(kind=RK), intent(out) :: gap_of_nodes(3,3)
        real(kind=RK), intent(out) :: mid_surf_coords(3,3)

        integer :: j
        real(kind=RK) :: reshape_u(dof_num_of_a_node, node_num_of_an_ele)

        reshape_u = reshape(u(i,:), (/dof_num_of_a_node, node_num_of_an_ele/))

        do j = 1, half_node_num_of_an_ele
            gap_of_nodes(:, j) = coords(i, j + half_node_num_of_an_ele,:) - coords(i, j, :) + reshape_u(1:3, j + half_node_num_of_an_ele) - reshape_u(1:3, j)
            mid_surf_coords(:, j) = 0.5d0 * ((coords(i, j + half_node_num_of_an_ele, :) + coords(i, j, :)) +  (reshape_u(1:3, j + half_node_num_of_an_ele) + reshape_u(1:3, j)))
        end do

    end subroutine seperation_and_middle_surface


    subroutine k_local_coordinates(co_de, coord_l,coords,u,ndofel,nnode, mcrd)
        implicit none

        integer, intent(in) :: ndofel, nnode, mcrd
        real(kind=RK), intent(in) :: co_de(mcrd,nnode), coords(mcrd,nnode), u(ndofel)
        real(kind=RK), intent(inout) :: coord_l(mcrd,nnode)
        real(kind=RK) :: aJacob_M(2,3), aLen, a_Jacob, dum, dum1, dum2, dum3, Rn1
        real(kind=RK) :: co_de_m(3,3), SFD(2,4), R(mcrd,mcrd), Transformation_M(ndofel,ndofel), Transformation_M_T(ndofel,ndofel)
        integer :: i, j, k, num


        !real(kind) R(mcrd,mcrd),coord_l(mcrd,nnode),aJacob_M(2,3),
        !& Transformation_M(ndofel,ndofel),coords(mcrd,nnode),
        !& Transformation_M_T(ndofel,ndofel),u(ndofel),
        !& co_de(mcrd,nnode), co_de_m(3,3),SFD(2,4)
   
        co_de_m = 0.0d0
        !计算中间截面坐标
        do i = 1, 3
           co_de_m(i,1)=(co_de(i,1)+co_de(i,4))*0.5
           co_de_m(i,2)=(co_de(i,2)+co_de(i,5))*0.5
           co_de_m(i,3)=(co_de(i,3)+co_de(i,6))*0.5
        end do
         
        !形函数导数
        SFD(1,1) =-1
        SFD(1,2) = 1
        SFD(1,3) = 0
        SFD(2,1) =-1
        SFD(2,2) = 0
        SFD(2,3) = 1
        
        aJacob_M = 0.0d0

        do i = 1,2
           do j = 1,3
              do k = 1, 3
                 aJacob_M(i,j) = aJacob_M(i,j) + SFD(i,k)*co_de_m(j,k)
              end do
           end do
        end do
   
        dum1 = aJacob_M(1,2)*aJacob_M(2,3) - aJacob_M(1,3)*aJacob_M(2,2)
        dum2 = aJacob_M(1,3)*aJacob_M(2,1) - aJacob_M(1,1)*aJacob_M(2,3)
        dum3 = aJacob_M(1,1)*aJacob_M(2,2) - aJacob_M(1,2)*aJacob_M(2,1)
   
        a_Jacob = sqrt(dum1**2 + dum2**2 + dum3**2) / 2.0d0
        Rn1 = sqrt(dum1**2 + dum2**2 + dum3**2)

        R = 0.0d0
        R(3,1) = dum1/Rn1
        R(3,2) = dum2/Rn1
        R(3,3) = dum3/Rn1
   
        aLen=sqrt(aJacob_M(1,1)**2.0 + aJacob_M(1,2)**2.0 + aJacob_M(1,3)**2.0)
        R(1,1)=aJacob_M(1,1)/aLen
        R(1,2)=aJacob_M(1,2)/aLen
        R(1,3)=aJacob_M(1,3)/aLen

        R(2,1)=R(3,2)*R(1,3)-R(3,3)*R(1,2)
        R(2,2)=R(3,3)*R(1,1)-R(3,1)*R(1,3)
        R(2,3)=R(3,1)*R(1,2)-R(3,2)*R(1,1)
   

        num=nnode
        Transformation_M = 0.0d0
        do i = 1, num
           dum=3.0*(i-1.0)
           Transformation_M(dum+1,dum+1)=R(1,1)
           Transformation_M(dum+1,dum+2)=R(1,2) 
           Transformation_M(dum+1,dum+3)=R(1,3)
           Transformation_M(dum+2,dum+1)=R(2,1)
           Transformation_M(dum+2,dum+2)=R(2,2)
           Transformation_M(dum+2,dum+3)=R(2,3)
           Transformation_M(dum+3,dum+1)=R(3,1)
           Transformation_M(dum+3,dum+2)=R(3,2)
           Transformation_M(dum+3,dum+3)=R(3,3)
        end do

        Transformation_M_T = 0.0d0
        call k_matrix_transpose(Transformation_M,Transformation_M_T, ndofel,ndofel)

        do i = 1, nnode
           coord_l(1,i)=(R(1,1)*co_de(1,i)+R(1,2)*co_de(2,i) + R(1,3)*co_de(3,i))
           coord_l(2,i)=(R(2,1)*co_de(1,i)+R(2,2)*co_de(2,i) + R(2,3)*co_de(3,i))
           coord_l(3,i)=(R(3,1)*co_de(1,i)+R(3,2)*co_de(2,i) + R(3,3)*co_de(3,i))
        end do

        return
    end subroutine k_local_coordinates

    subroutine k_matrix_transpose(A,B,n,m)
        implicit none
        integer, intent(in) :: n,m
        real(kind=RK), intent(in) :: A(n,m)
        real(kind=RK), intent(out) ::B(m,n)
        integer :: i,j
            do i = 1, n
               do j = 1, m
                  B(j,i)=A(i,j)
               end do
            end do
        return
    end subroutine k_matrix_transpose


    !============================================get now coords================================================
    subroutine get_now_coords(i, u, coords, now_coords)
        implicit none
        integer, intent(in) :: i
        real(kind=RK), intent(in) :: u(num_of_ele, dof_num_of_an_ele)
        real(kind=RK), intent(in) :: coords(num_of_ele, node_num_of_an_ele, dimmension)
        real(kind=RK), intent(out) :: now_coords(3,3)

        integer :: j
        real(kind=RK) :: reshape_u(dof_num_of_a_node, node_num_of_an_ele)

        reshape_u = reshape(u(i,:), (/dof_num_of_a_node, node_num_of_an_ele/))

        do j = 1, node_num_of_an_ele
            now_coords(:, j) = coords(i, j, :) + reshape_u(1:3, j)
        end do
    end subroutine get_now_coords


    !============================================vuel_vumat================================================
    subroutine vuel_vumat(i, u, du, coords, rhs_temp, svars, nsvars, time, props, nprops)
        integer, intent(in) :: i, nsvars
        integer, intent(in) :: nprops
        real(kind=RK), intent(in) :: time(2)
        real(kind=RK), intent(in) :: props(nprops)
        real(kind=RK), intent(in) :: u(num_of_ele, dof_num_of_an_ele)
        real(kind=RK), intent(in) :: du(num_of_ele, dof_num_of_an_ele)
        real(kind=RK), intent(in) :: coords(num_of_ele, node_num_of_an_ele, dimmension)
        real(kind=RK), intent(inout) :: rhs_temp(dof_num_of_an_ele)
        real(kind=RK), intent(inout) :: svars(num_of_ele, nsvars)

        real(kind=RK) :: gap_of_nodes(3,3)                  !! the gap of the nodes of the integration point
        real(kind=RK) :: mid_surf_coords(3,3)               !! the middle surface coordinates
        real(kind=RK) :: shear_inc(3, 3)                    !! the shear increment of the integration point
        real(kind=RK) :: unit_normal_vector(3)           !! the unit normal vector of the integration point
        real(kind=RK) :: area(3)                            !! the area of the integration point
        real(kind=RK) :: gap_of_intergration_points(3,3)    !! the gap of the integration points
        real(kind=RK) :: normal_seperation(3)               !! the normal seperation
        real(kind=RK) :: tangential_seperation(3)           !! the tangential seperation
        real(kind=RK) :: node_tangential_vec(3)             !! the tangential vector of the node
        real(kind=RK) :: unit_tangential_vec(3)             !! the unit tangential vector
        real(kind=RK) :: delt_n_inc                         !! the normal displacement increment in this incrace time step
        !real(kind=RK) :: delt_shear_inc                     !! the shear displacement increment in this incrace time step
        real(kind=RK) :: sigm_n(3)                          !! the normal stress
        real(kind=RK) :: sigm_shear(3)                      !! the shear stress
        !real(kind=RK) :: delt_shear_inc_vec(3)              !! the shear displacement increment vector
        real(kind=RK) :: stress_vec(3)                      !! the stress vector
        real(kind=RK) :: force_vec(3)                       !! the force vector
        real(kind=RK) :: unit_vec_t(3)                      !! the unit tangential vector
        real(kind=RK) :: unit_vec_s(3)                      !! the unit tangential vector
        real(kind=RK) :: delt_s, delt_t
        integer :: j

        real(kind=RK) :: unit_tangential_vec_ele(3, 3), shear_s(3), shear_t(3), unit_tangential(3, 3), sigm_shear_vec(3)

        !real(kind=RK) :: now_coords(3,3), coord_l(3,3), ds1(3), ds2(3), dn(3), ds_shear(3)
        !integer :: x
        !now_coords = 0.0d0
        !coord_l = 0.0d0
        !call get_now_coords(i, u, coords, now_coords) !! calculate the current coordinates
        !call k_local_coordinates(now_coords, coord_l, coords, u, dof_num_of_an_ele, node_num_of_an_ele, dimmension) !! calculate the local coordinates
        !do x = 1, 3
        !    ds1(x)=coord_l(1,x+3)-coord_l(1,x)
        !    ds2(x)=coord_l(2,x+3)-coord_l(2,x)
        !    dn(x) =coord_l(3,x+3)-coord_l(3,x)
        !end do
        !do x = 1, 3
        !    ds_shear(x) = sqrt(ds1(x)**2 + ds2(x)**2)
        !end do
        !write(*,*) "ds1= ", ds1, " ds2= ", ds2, " dn= ", dn, " ds_shear= ", ds_shear

        gap_of_intergration_points = 0.0d0

        call seperation_and_middle_surface(i, u, coords, gap_of_nodes, mid_surf_coords) !! calculate the seperation and middle surface

        call get_shear_inc(i, du, shear_inc)                                            !! calculate the shear increment

        call get_unit_normal_vector( mid_surf_coords, unit_normal_vector )         !!这里可以优化，只需要计算一次，因为三个点只有一个平面，所以只有一个法向量

        call get_local_tangential_basis_vector(mid_surf_coords, unit_normal_vector, unit_vec_t, unit_vec_s)

        call seperation_component_and_vector_ele(mid_surf_coords, gap_of_nodes, unit_vec_t, unit_vec_s, unit_normal_vector, normal_seperation, tangential_seperation, unit_tangential_vec_ele, shear_s, shear_t, unit_tangential)

        do j = 1, half_node_num_of_an_ele

            call calculate_area(j, mid_surf_coords, area(j))

            !call gap_of_integral_points(j, gap_of_nodes, gap_of_intergration_points(:, j), Tet_gaussPoint(:, j))

            !call seperation_component_and_vector(j, gap_of_intergration_points(:, j), unit_normal_vector, normal_seperation(j), tangential_seperation(j), unit_tangential_vec)

            !delt_s = dot_product_1(gap_of_intergration_points(:, j), unit_vec_s)
            !delt_t = dot_product_1(gap_of_intergration_points(:, j), unit_vec_t)
            !write(*,*) "j = ", j, "delt_s= ", delt_s, " delt_t= ", delt_t

            delt_n_inc = dot_product_1(shear_inc(:, j), unit_normal_vector)
            !delt_shear_inc_vec = shear_inc(:, j) - delt_n_inc * unit_normal_vector
            !delt_shear_inc = vector_mod(delt_shear_inc_vec)

            !write(*,*) "inc: ", shear_inc(1, j), shear_inc(2, j), shear_inc(3, j)

            call vuel_vumat_constitutive(i, j, time, props, nprops, svars(i, :), nsvars,                    &
                                         normal_seperation(j), tangential_seperation(j), delt_n_inc,        &
                                         sigm_n(j), sigm_shear(j), shear_s(j), shear_t(j), unit_vec_t, unit_vec_s)

            !!write(*,*) "here 9"
            !stress_vec = sigm_n(j) * unit_normal_vector + sigm_shear(j) * unit_tangential_vec_ele(:, j)
            sigm_shear_vec = sigm_shear(j) * unit_tangential(:, j)
            stress_vec = sigm_n(j) * unit_normal_vector + sigm_shear_vec
            force_vec = stress_vec * area(j)
            !!write(*,*) "here 10"

            svars(i, ((j - 1) * 12 + 7) ) = sigm_shear_vec(1)
            svars(i, ((j - 1) * 12 + 8) )= sigm_shear_vec(2)
            svars(i, ((j - 1) * 12 + 9) )= sigm_shear_vec(3)
            svars(i, ((j - 1) * 12 + 10)) = unit_tangential_vec_ele(1, j)
            svars(i, ((j - 1) * 12 + 11)) = unit_tangential_vec_ele(2, j)
            svars(i, ((j - 1) * 12 + 12)) = unit_tangential_vec_ele(3, j)

            call distribute_force(i, j, force_vec, rhs_temp, Tet_gaussPoint(:, j) )

        end do   
    end subroutine vuel_vumat


    !============================================calculate the mass matrix================================================
    subroutine mass_matrix(nblock,nnode,mcrd,ndofel,coords,amass)
        implicit none
        integer, intent(in) :: nblock,nnode,mcrd,ndofel
        real(kind=RK), intent(in) :: coords(nblock,nnode,mcrd)
        real(kind=RK), intent(out) :: amass(nblock,ndofel,ndofel)
        integer :: i, j
        real(kind=RK) :: area, am0, S1, S2, elemMass
        real(kind=RK) :: edge(3,3), V1(3), V2(3)

        do i = 1, nblock
            edge(1,:) = coords(i,2,:) - coords(i,1,:)
            edge(2,:) = coords(i,3,:) - coords(i,1,:)
            edge(3,:) = coords(i,4,:) - coords(i,1,:)

            call vector_cross( edge(1,:),edge(2,:),V1 )

            S1 = vector_mod( V1 )

            call vector_cross(edge(2,:),edge(3,:),V2)

            S2      = vector_mod( V2 )
            area    = 0.5d0 * ( S1 + S2 )
            elemMass = thickness * area * rho
            am0     = elemMass / real(node_num_of_an_ele,kind=RK)

            !write(*,*) "am0 = ", am0

            do j = 1, dof_num_of_an_ele
                ! When the mass matrix is calculated, lflags(iProcedure) = 17, that is, there are only 18 degrees of freedom are involved, and
                ! the contribution of the temperature degree of freedom to the mass matrix is zero.
                amass(i,j,j) = am0
            end do
        end do
    end subroutine mass_matrix

    !============================================initialize the parameters================================================
    subroutine set_parameter(nblock,nnode,mcrd,ndofel,nprops,props)
        implicit none
        integer, intent(in) :: nblock,nnode,mcrd,ndofel,nprops
        real(kind=RK), intent(in) :: props(nprops)

        num_of_ele = nblock
        node_num_of_an_ele = nnode
        half_node_num_of_an_ele = nnode / 2
        dof_num_of_a_node = ndofel / nnode
        dof_num_of_an_ele = ndofel
        dimmension = mcrd
        rho = props(11)
    end subroutine set_parameter

end module cohesive_element_mod

!================================================VUEL=========================================================
subroutine vuel(nblock,rhs,amass,dtimeStable,svars,nsvars,energy,nnode,ndofel,  &
                props,nprops,jprops,njprops,coords,mcrd,u,du,v,a,               &
                jtype,jElem,time,period,dtimeCur,dtimePrev,kstep,kinc,lflags,   &
                dMassScaleFactor,predef,npredef,jdltyp,adlmag)
    !     used module
    use system_constant
    use math_mod
    use parameter_mod
    use cohesive_element_mod
    !     if you want to use `implicit none`, you must add the folowing three lines of code
    implicit none
    integer, parameter :: j_sys_Dimension   = 2
    integer, parameter :: maxblk            = 512
    
    !     operation code keys
    integer, parameter :: jMassCalc             = 1
    integer, parameter :: jIntForceAndDtStable  = 2
    integer, parameter :: jExternForce          = 3
    !     flags indices
    integer, parameter :: iProcedure = 1
    integer, parameter :: iNlgeom    = 2
    integer, parameter :: iOpCode    = 3
    integer, parameter :: nFlags     = 3
    !     energy array indices
    integer, parameter :: iElPd     = 1
    integer, parameter :: iElCd     = 2
    integer, parameter :: iElIe     = 3
    integer, parameter :: iElTs     = 4
    integer, parameter :: iElDd     = 5
    integer, parameter :: iElBv     = 6
    integer, parameter :: iElDe     = 7
    integer, parameter :: iElHe     = 8
    integer, parameter :: iUnused   = 9
    integer, parameter :: iElTh     = 10
    integer, parameter :: iElDmd    = 11
    integer, parameter :: iElDc     = 12
    integer, parameter :: nElEnergy = 12
    !     time indices
    integer, parameter :: iStepTime  = 1
    integer, parameter :: iTotalTime = 2
    integer, parameter :: nTime      = 2
    !     procedure flags
    integer, parameter :: jDynExplicit     = 17   ! Direct-Integration Explicit Dynamic Analysis
    integer, parameter :: jDynExplicitTM   = 74   ! Transient Fullly Coupled Thermal-Stress Analysis
    !     predefined variables indices
    integer, parameter :: iPredValueNew = 1
    integer, parameter :: iPredValueOld = 2
    integer, parameter :: nPred         = 2 

    integer, intent(in) :: jType    !! VUn, n = jType
    integer, intent(in) :: jdlTyp   !! integer identifying the load type
    integer, intent(in) :: kStep    !! Current step number
    integer, intent(in) :: kinc     !! Current increment number
    integer, intent(in) :: npredef  !! field variable number   not used now  
    integer, intent(in) :: nblock   !! total solid element number
    integer, intent(in) :: mcrd     !! n coordinates
    integer, intent(in) :: ndofel   !! total dof number of one element = 32;  right
    integer, intent(in) :: njprops  !! user defined properties (integer)
    integer, intent(in) :: nnode    !! total node number of one element  (real) =8
    integer, intent(in) :: nprops   !! The number of user defined properties 
    integer, intent(in) :: nsvars   !! The number of SDV of one element
    real(kind=RK), intent(in) :: dTimeCur   !! Current time increment
    real(kind=RK), intent(in) :: dTimePrev  !! Previous time increment.
    real(kind=RK), intent(in) :: period     !! Time period of the current step
    real(kind=RK), intent(out) :: rhs(nblock,ndofel) !! the contributions of each element to the right-hand-side vector of the overall system of equations (e.g.internal force)
    real(kind=RK), intent(out) :: amass(nblock,ndofel,ndofel)  !!  mass matrix of each element
    real(kind=RK), intent(inout) :: dtimeStable(nblock) !! stable time increment of each element
    real(kind=RK), intent(out) :: svars(nblock,nsvars) !! the solution-dependent state variables of all elements
    real(kind=RK), intent(in)  :: energy(nblock,nElEnergy)  !! not computed
    real(kind=RK), intent(in)  :: props(nprops) !! User-defined number of real property values associated with the elements processed
    integer, intent(in)        :: jProps(njProps)  !! User-defined number of integer property values 
    integer, intent(in)        :: jElem(nblock)    !! the element ID
    real(kind=RK), intent(in)  :: time(nTime)    !! Current value of total time
    integer, intent(in)        :: lFlags(nFlags)
    real(kind=RK), intent(in)  :: coords(nblock,nnode,mcrd) !! An array containing the original coordinates of the nodes of the element
    real(kind=RK), intent(in)  :: u(nblock,ndofel) !! total displacements
    real(kind=RK), intent(in)  :: du(nblock,ndofel) !! Incremental displacements in the current increment
    real(kind=RK), intent(in)  :: v(nblock,ndofel)  !! velocities at the midpoint of the increment
    real(kind=RK), intent(in)  :: a(nblock, ndofel) !! Accelerations at the end of the current increment
    real(kind=RK), intent(in)  :: dMassScaleFactor(nblock)  !! the mass scale factors for each element
    real(kind=RK), intent(in)  :: preDef(nblock, nnode, nPreDef, nPred)  !! predefined field variable
    real(kind=RK), intent(in)  :: adlmag(nblock) !! the load magnitude of the load type jdltyp for each elem 

    integer :: i
    real(kind=RK) :: rhs_temp(ndofel)

    !     Debug
    integer :: tempRead   !! debug flag
    
    if( (jType == 2 .and. lflags(iProcedure) == jDynExplicit) .or. (jType == 2 .and. lFlags(iProcedure) == jDynExplicitTM) ) then 
    !if( jType == 2 .and. lFlags(iProcedure) == jDynExplicitTm ) then    !热力耦合分析
    !if( jType==2 .and. lflags(iProcedure)==jDynExplicit ) then          !直接积分的显示动力学分析
        
        call set_parameter(nblock, nnode, mcrd, ndofel, nprops, props)
        
        if( lflags(iOpCode) == jMassCalc )then  
            !!!write(*,*) "calculate the mass matrix"
            call mass_matrix(nblock, nnode, mcrd, ndofel, coords, amass) 
        
        else if ( lflags(iOpCode) == jIntForceAndDtStable) then
                
              !if (time(2) > 0.0d0) then
              !        write(*,*) "Please input an integer:"
              !        read(*,*) tempRead
              !end if
                
            do i = 1, nblock
                rhs_temp = 0.0d0        !initialize the rhs_temp

                call vuel_vumat(i, u, du, coords, rhs_temp, svars, nsvars, time, props, nprops)
                
                rhs(i, :) = rhs_temp             
            end do
            !open(2019, file = "G:\ABAQUS\VUEL\TM-test-10-1-18\C3D8T-vuel-test\data.csv", position = "append")
            !write(2019,'(*(G0.5, :, ",", X))') time(1), rhs, u
            !close(2019) 
        end if
    end if
    return
end subroutine vuel
