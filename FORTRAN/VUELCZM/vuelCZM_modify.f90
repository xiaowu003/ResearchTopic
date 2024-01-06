module system_constant
    implicit none
    !integer, parameter :: RK4 = selected_real_kind(p = 6)
    !integer, parameter :: RK8 = selected_real_kind(p = 15)
    integer, parameter :: RK = selected_real_kind(p = 15)
    integer, parameter :: MNL = 256 ! max_name_length
    !real(kind = RK), private :: temp
    !integer, parameter :: IHUGE = huge(IK1)
    !real(kind = RK), parameter :: RHUGE = huge(temp)
    !real(kind = RK), parameter :: RTINY = tiny(temp)
end module system_constant

module parameter_mod
    use system_constant
    implicit none
    integer, parameter :: nInteg1D = 2                          !!  
    integer, parameter :: nInteg2D = 4                          !!平面积分点数量 4
    integer, save :: nCoh                                       !!内聚力单元数？ #todo
    integer, save :: nNdCoh = 8                                 !!一个内聚力单元节点数 8
    integer, save :: nHalfNd = 4                                !!一个内聚力单元一半节点数 4
    integer, save :: nDim = 3                                   !!自由度数 3
    integer, save :: nDofCoh = 24                               !!一个内聚力单元总自由度数 24
    real(kind=RK), save :: m_thickness=1.0d0                    !!内聚力单元厚度 1
    real(kind=RK), save :: m_rho                                !!密度
    real(kind=RK), save :: m_nStrength                          !!拉伸强度
    real(kind=RK), save :: m_tStrength                          !!剪切强度
    real(kind=RK), save :: m_Gnc                                !! 2型临界能量释放率，即剪切模式
    real(kind=RK), save :: m_Gtc                                !! 1型临界能量释放率，即拉伸模式
    real(kind=RK), save :: m_alpha                              !! = props(8)  即 0.3
    real(kind=RK), save :: m_deltaNo                            !!
    real(kind=RK), save :: m_deltaSo                            !!
    real(kind=RK), save :: m_stiff                              !!法向罚刚度
    real(kind=RK), save :: xiGaussPoint(2,30)                   !!（不是高斯点）积分点, one element consists of four Gaussian integral points
    real(kind=RK), save :: weight(30)                           !!权重

end module parameter_mod

module math_mod
    use system_constant
    implicit none
    private
    public :: vector_cross, vector_dot
    public :: vector_mod, vector_unit
    public :: shap_function_derivative
    public :: shape_q4
    contains

    subroutine vector_cross(vec1,vec2,cross) !!向量叉乘
        use system_constant
        implicit none
            !
            ! Function:
            !
            !  Calculate the cross product of two vectors.
            !
            ! Remarks:
            !
            !   RESULT = VEC1    VEC2
            !
            ! Parameters:
            !
            !  Vec1, input, an array of rank one with 3 elements.
            !  Vec2, input, an array of rank one with 3 elements.
            !  result, output, an array of rank one with 3 elements.

        real(kind=RK),intent(in):: vec1(3),vec2(3)
        real(kind=RK),intent(out):: cross(3)
        cross(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
        cross(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
        cross(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
        return
    end subroutine vector_cross
     
    subroutine vector_dot( vec1, vec2, dot )
        implicit none 
        !
        ! Function:
        !
        !  Calculate the dot product of two vectors.
        !
        ! Remarks:
        !
        !   DOT = VEC1 * VEC2
        !
        ! Parameters:
        !
        !  Vec1, input, an array of rank one with 3 elements.
        !  Vec2, input, an array of rank one with 3 elements.
        !  DOT, output, an array of rank one with 3 elements.

        real(kind=RK),intent(in) :: vec1(3), vec2(3)
        real(kind=RK),intent(out):: dot

        dot=vec1(1)*vec2(1) + vec1(2)*vec2(2) + vec1(3)*vec2(3)
        return
    end subroutine vector_dot
      
    function vector_mod( vec ) result(rmod)
        implicit none

        real(kind=RK) :: vec(3)
        real(kind=RK) :: rmod

        rmod = sqrt( vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3) )
        return
    end function vector_mod

    subroutine vector_unit( vec )
        implicit none
        real(kind=RK), intent(inout) :: vec(3)
        real(kind=RK), parameter :: TINY = 1.0e-25
        real(kind=RK) :: rmod

        rmod = sqrt( vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3) )

        if( rmod < TINY) then
            write(6,*)'Error in vector_unit:the rmod of the vector is 0.'
            call xplb_exit
        end if

        vec = vec/rmod
        return
    end subroutine vector_unit

    subroutine shap_function_derivative(gpos,dshap1,dshap2)  !!形函数的导数
        implicit none
        real(kind=RK),intent(in)::gpos(2)
        real(kind=RK),intent(out)::dshap1(4),dshap2(4)
        
        ! Take the partial derivative of r, and r is gpos(1)
        dshap1(1)=-(1.0d0-gpos(2))
        dshap1(2)=1.0d0-gpos(2)
        dshap1(3)=1.0d0+gpos(2)
        dshap1(4)=-(1.0d0+gpos(2))
        dshap1 = dshap1 / 4.0d0
        
        ! Take the partial derivative of s, and s is gpos(2)
        dshap2(1)=-(1.0d0-gpos(1))
        dshap2(2)=-(1.0d0+gpos(1))
        dshap2(3)=1.0d0+gpos(1)
        dshap2(4)=1.0d0-gpos(1)
        dshap2 = dshap2 / 4.0d0
        
        return
    end subroutine shap_function_derivative
       
    subroutine shape_q4 ( r, s, t )
        implicit none
        !! SHAPE_Q4 evaluates shape functions for a 4 node reference quadrilateral.
        !  Reference Element Q4:
        !    |
        !    1  4-----3
        !    |  |     |
        !    |  |     |
        !    S  |     |
        !    |  |     |
        !    |  |     |
        !   -1  1-----2
        !    |
        !    +-(-1)--R--(1)-->
        !
        !  Parameters:
        !    Input : R, S, the reference coordinates of a point.
        !    Output: T(4), the basis functions at the point.
        real(kind=RK), intent(in) :: r, s
        real(kind=RK), intent(out) :: t(4)

        t(1) = 0.25d0 * ( 1.0d0 - r ) * ( 1.0d0 - s )
        t(2) = 0.25d0 * ( 1.0d0 + r ) * ( 1.0d0 - s )
        t(3) = 0.25d0 * ( 1.0d0 + r ) * ( 1.0d0 + s )
        t(4) = 0.25d0 * ( 1.0d0 - r ) * ( 1.0d0 + s )     
        return
    end subroutine shape_q4

end module math_mod

module cohesive_law_mod
    use system_constant
    use math_mod
    implicit none
    public :: compute_input_coh_displacement
    public :: cohesive_stress_at_one_point
    public :: calculate_effective_displacement_point
    public :: calculate_tangent_vector
    character(len=MNL) :: msg
    contains
    
    subroutine compute_input_coh_displacement(stiff,nStreng,tStreng,GN,GS,alfa,deltaNo,deltaSo,tdelta,ndelta,deltai,deltaf) !! 计算单元开始损伤时的分离量和完全损伤时的分离量
        
        implicit none
        real(kind=RK),intent(in) :: stiff !!法向罚刚度
        real(kind=RK),intent(in) :: nStreng !!拉伸强度
        real(kind=RK),intent(in) :: tStreng !!剪切强度
        real(kind=RK),intent(in) :: GN !!一型断裂能  GIC
        real(kind=RK),intent(in) :: GS !!二型断裂能  GIIC
        real(kind=RK),intent(in) :: alfa !! 0.3
        real(kind=RK),intent(in) :: deltaNo  !!拉伸分离量？ deltaI0
        real(kind=RK),intent(in) :: deltaSo  !!剪切分离量？ deltaII0
        real(kind=RK),intent(in) :: tdelta  !! =tGap 上下两节点的切向分离量
        real(kind=RK),intent(in) :: ndelta  !! =nGap 上下两节点的法向分离量
        real(kind=RK),intent(out) :: deltai !!输出 开始损伤时的分离量
        real(kind=RK),intent(out) :: deltaf !!输出 完全损伤时的分离量
        real(kind=RK)::beita,nsig,tsig
        real(kind=RK)::gna,gta,beita2,alfam,rnsig

        deltai = 0.0d0
        deltaf = 0.0d0
        
        call calculate_mixed_mode_rate(tdelta,ndelta,beita)
        
        beita2 = beita*beita  ! #bug 计算混合比的平方？
        
        !beita2=0.7
        !for reference: Pinho S T, Iannucci L, Robinson P. Formulation and implementation
        !of decohesion elements in an explicit finite element code[J].
        !Composites Part A: Applied science and manufacturing, 2006, 37(5): 778-789.
        !call obtain_transition_value_in_intrinsic_model(stiff,nStrengh,tStrength,delta_no,delta_so)
        
        if(ndelta<=0.0d0) then
            deltai=deltaSo !effective displacement at transition point.   法向分离量小于等于0时，总分离量=切向分离量

            if(deltai<=0.0d0) goto 222
            if(tStreng<=0.0d0) goto 111

            deltaf=2.0d0 * Gs / tStreng !final displacement of the mixed mode.  !! Gs: G2c; Gn: G1c; #bug 与论文中公式有所不同,  
        else
            !effective displacement at transition point of the mixed triangle.
            deltai = deltaNo * deltaSo * sqrt( (1+beita2) / ( beita2*deltaNo*deltaNo + deltaSo*deltaSo) )  !混合模式 开始损伤时分离量 计算公式

            if(deltai<=0.0d0) goto 222
            if(stiff<=0.0d0) goto 333

            deltaf=2.0d0*(1+beita2)/(stiff*deltai)

            if(GN<=0.0d0.or.GS<=0.0d0) goto 444

            gna=(1.0d0/GN)**alfa; gta=(beita2/GS)**alfa
            alfam=-1.0d0/alfa
            deltaf=deltaf * (gna+gta)**alfam !final displacement of the mixed mode.

        end if
        return

        111 continue
        msg='Negative tstrength in compute_input_coh_displacement.'
        write(6,*) msg
        call xplb_exit
        222 continue
        msg='Negative deltai in compute_input_coh_displacement.'
        write(6,*) msg
        call xplb_exit
        333 continue
        msg='Negative stiffness in compute_input_coh_displacement.'
        write(6,*) msg
        call xplb_exit
        444 continue
        msg='Negative GN or GS in compute_input_coh_displacement.'
        write(6,*) msg
        call xplb_exit

    end subroutine compute_input_coh_displacement

    subroutine compress_force(stiff, depth,normal,area,pfac,facc)
    
        implicit none
        real(kind=RK),intent(in):: stiff
        real(kind=RK),intent(in):: depth,normal(3),area
        real(kind=RK),intent(out):: facc(3)
        real(kind=RK),intent(in):: pfac
        real(kind=RK)::det,norm(3)
        real(kind=RK)::ustiff=0.0d0
        
        facc=0.0d0
        !ustiff=intrin_stiff
        ustiff=pfac*stiff
        facc=abs(depth)*ustiff*normal
        return
    end subroutine compress_force

    subroutine cohesive_stress_at_one_point(stiff,deltai,deltaf,delta,det_ndelta,det_tdelta,eff_delta,ghis,nvect,tvect,fail,sig,damage)
        !该子程序的参数解释为GPT生成，可能有错误，请参考
        implicit none
        real(kind=RK),intent(in):: stiff            !!刚度
        real(kind=RK),intent(in):: delta(3)         !!表示相对位移的向量
        real(kind=RK),intent(in):: det_ndelta       !!下方节点指向上方节点的向量 与 积分点单位法向量 的 点乘，所以可能会出现 负值
        real(kind=RK),intent(in):: det_tdelta       !!切向分离量
        real(kind=RK),intent(in):: ghis             !!实际分离量
        real(kind=RK),intent(in):: eff_delta        !!合成的总分离量
        real(kind=RK),intent(in):: nvect(3)         !!单位法向量
        real(kind=RK),intent(in):: tvect(3)         !!计算得出单位切向量
        real(kind=RK),intent(in):: deltai           !!开始损伤分离量
        real(kind=RK),intent(in):: deltaf           !!完全失效分离量
        !integer,intent(in)::co_eid, gith
        real(kind=RK),intent(out):: fail            !!失效标志
        real(kind=RK),intent(out):: sig(3)          !!计算得到的应力向量， 法向应力向量 + 切向应力向量
        real(kind=RK):: nsig(3)                     !!法向应力向量， 损伤后的法向应力 * 单位法向量
        real(kind=RK):: tsig(3)                     !!切向应力向量， 损伤后的切向应力 * 单位切向量
        real(kind=RK):: dnsig                       !!开始损伤后的法向应力 = 刚度 * 法向分离量 * (1 - d)
        real(kind=RK):: dtsig                       !!开始损伤后的切向应力 = 刚度 * 切向分离量 * (1 - d)
        real(kind=RK):: det
        real(kind=RK):: temp
        real(kind=RK):: eta
        real(kind=RK):: modu     !!由上下两节点间的向量计算出的模，即上下两节点的距离
        real(kind=RK):: bvect(3)
        real(kind=RK):: delta_tc
        real(kind=RK):: rnvect(3)                   !!单位法向量
        real(kind=RK):: rtvect(3)                   !!单位切向量
        real(kind=RK), parameter::tol=1.0e-15
        real(kind=RK):: nsig0, tsig0, effsig0 !normal and tangential stress components at the transition point in a TSL.
        real(kind=RK):: coef=0.0d0                      !!损伤变量D，#bug 这里为什么将它初始化为0？没有确保该值始终不会变小
        real(kind=RK),intent(inout) :: damage       !!传入的旧损伤变量和新计算的损伤变量比较进行更新
        
        fail    =   0.0d0
        sig     =   0.0d0
        dnsig   =   0.0d0
        dtsig   =   0.0d0
        coef    =   0.0d0

        !if(interface_failure/=1) then
        if(ghis <= 0.0d0) goto 111
        !end if

        !for reference: Wu L, Sket F, Molina-Aldareguia J M, et al. A study 
        ! of composite laminates failure using an anisotropic gradient-enhanced 
        !damage mean-field homogenization model[J]. Composite Structures, 2015, 126: 246-264.
        !for reference: Hu N, Zemba Y, Okabe T, et al. A new cohesive model for simulating delamination propagation
        !in composite laminates under transverse loads[J]. Mechanics of Materials, 2008, 40(11): 920-935.
        !if(interface_failure==1)then   !intrinsic cohesive modeling
        
        if(ghis <= deltai) then
            if(det_ndelta > 0.0d0) then
                !!弹性阶段
                !nsig = stiff * det_ndelta * nvect    !#bug 这里两个实数相乘得到了数组，错误, 再乘单位法向量也错，这里需要的是发向应力值而不是数组
                dnsig = stiff * det_ndelta !Right!
            end if

            dtsig = stiff * det_tdelta    !#bug 这里的刚度和计算法向应力的shift刚度一样吗？

            damage = 0.0d0   !分离量小于开始损伤分离量，故积分点还未损伤

        else if(ghis > deltai .and. ghis <= deltaf) then
            !!损伤阶段
            fail = 1.01d0
            coef = deltaf / ghis
            coef = coef * (ghis - deltai) / (deltaf - deltai)  ! Damage variable calculation formula

            if(coef <= damage) then
                coef = damage
            else
                damage = coef   ! 更新损伤变量，上一步的损伤变量小于当前的，则更新
            end if

            if(det_ndelta > 0.0d0) then
                dnsig = (1.0d0 - coef) * stiff * det_ndelta   ! #bug 这里从1改成了1.0d0
            end if

            dtsig = (1.0d0 - coef) * stiff * det_tdelta
        end if
        !else
        !goto 555
        !end if

        modu = sqrt(delta(1)*delta(1) + delta(2)*delta(2) + delta(3)*delta(3))

        if(modu <= 0.0d0) goto 222

        bvect = delta
        bvect = bvect/modu  !!得到上下两节点间的单位向量
        rnvect = nvect
        rtvect = tvect
        
        ! #bug 下面代码完全没有作用
        !if(det_ndelta > 0.0d0) then
        !    call vector_dot(nvect,bvect,det)
        !end if
        !call vector_dot(tvect,bvect,det)

        if(det_ndelta > 0.0d0) then
            nsig = dnsig * rnvect
        else
            nsig = 0.0d0
        end if

        tsig = dtsig * rtvect
        sig  = nsig + tsig
        return

        111 continue
        msg='Negative ghis in cohesive_stress_at_one_point.'
        write(6,*) msg
        call xplb_exit
        222 continue
        msg='Zero mod in cohesive_stress_at_one_point.'
        write(6,*) msg
        call xplb_exit
    end subroutine cohesive_stress_at_one_point

    subroutine calculate_tangent_vector(tdelta,tvect)
        implicit none
        real(kind=RK),intent(in)::tdelta(3) !切向量
        real(kind=RK),intent(out)::tvect(3) !切向单位向量
        real(kind=RK)::det
        real(kind=RK),parameter::tol=1.0e-16

        det = sqrt(tdelta(1)*tdelta(1) + tdelta(2)*tdelta(2) + tdelta(3)*tdelta(3))
        if(det < tol) then  !若分离量 极小 直接当做 分离量为0 进行计算？
            tvect=(/0.0d0 ,0.0d0, 0.0d0/)
            return
        end if
        
        tvect = tdelta/det  !向量除以自身模长
        return
    end subroutine calculate_tangent_vector

    subroutine obtain_transition_value(stiff,nstreng,tstreng,delta_no,delta_so)

        !!该子程序可以直接获得输入参数对应的开始损伤分离量和完全损伤分离量

        implicit none
        real(kind=RK),intent(in) :: stiff
        real(kind=RK),intent(in) :: nstreng, tstreng
        real(kind=RK),intent(out):: delta_no,delta_so
        real(kind=RK),save::no=0.0d0,so=0.0d0
        integer,save::count=0
        if(count>0) then
            delta_no=no
            delta_so=so
            return
        end if

        count = count+1

        if(stiff<=0.0) goto 111

        no = nstreng/stiff !opening displacement at transition point in pure mode I (intrinsic cohesive model)
        so = tstreng/stiff !opening displacement at transition point in pure mode II (intrinsic cohesive model)
        delta_no = no
        delta_so = so
        return

        111 continue
        msg='Nagetive stiffness in obtain_transition_value.'
        write(6,*) msg
        call xplb_exit
    end subroutine obtain_transition_value

    subroutine calculate_effective_displacement_point(det_ndelta,det_tdelta,eff_delta)
        implicit none
        real(kind=RK),INTENT(in) :: det_ndelta           !!下方节点指向上方节点的向量 与 积分点单位法向量 的 点乘，所以可能会出现 负值
        real(kind=RK),INTENT(in) :: det_tdelta           !!上下节点在该（不是高斯点）积分点切向量方向的分量，即切向距离
        real(kind=RK),intent(out) :: eff_delta           !!输出 合成的总分离量
        real(kind=RK)::ndel,tdel

        ndel = (det_ndelta+abs(det_ndelta))*0.5d0  !#bug 这里为什么要限制法向分离量为 非负 ？
        tdel = det_tdelta
        eff_delta = sqrt(ndel*ndel+tdel*tdel)      !合成的总分离量
        return
    end subroutine calculate_effective_displacement_point

    subroutine calculate_mixed_mode_rate(tdelta,ndelta,beita)  !!计算混合比
        implicit none
        real(kind=RK),intent(in)::tdelta,ndelta
        real(kind=RK),intent(out)::beita        !!混合比

        beita = 0.0d0
        if(ndelta<=0.0d0) then
            beita=1.0e20            !如果法向分离量为负，则混合比为无穷大 
        else
            beita=tdelta/ndelta     !混合比 = 切向分离量 / 法向分离量  
        end if
        return
    end subroutine

end module cohesive_law_mod

module cohesive_element_mod
    use math_mod
    use parameter_mod
    use system_constant
    use cohesive_law_mod
    implicit none
    private
    public :: mass_matrix, set_parameter, one_element_cohesive_force
    character(len=MNL) :: msg
    contains

    subroutine set_parameter(nBlock,nNode,nCrd,nDofEl,nProps,props) !!参数传递
        implicit none
        integer, intent(in) :: nBlock,nNode,nCrd,nDofEl,nProps
        real(kind=RK), intent(in) :: props(nProps)
        real(kind=RK), parameter :: constant = 1.0d0 !0.577350269189626d0

        nCoh    = nBlock
        nNdCoh  = nNode
        !nInteg2D = nNode/2
        nHalfNd = nNdCoh/2
        nDim    = nCrd
        nDofCoh = nDofEl
        m_stiff = props(1)
        m_nStrength = props(3)
        m_tStrength = props(4)
        m_Gnc   = props(5) !! m_GIC
        m_Gtc   = props(6)
        m_rho   = props(7)
        m_alpha = props(8)
        m_deltaNo = m_nStrength/m_stiff
        m_deltaSo = m_tStrength/m_stiff

        ! Four Gaussian integral points
        xiGaussPoint(1:2,1) = (/-constant,  -constant /)
        xiGaussPoint(1:2,2) = (/ constant,  -constant /)
        xiGaussPoint(1:2,3) = (/ constant,   constant /)
        xiGaussPoint(1:2,4) = (/-constant,   constant /)
        weight(1:4) = 1.0d0
        return
    end subroutine set_parameter

    subroutine mass_matrix(nBlock,nNode,nCrd,nDofEl,coords,amass)
        implicit none
        integer, intent(in) :: nBlock,nNode,nCrd,nDofEl
        real(kind=RK), intent(in) :: coords(nBlock,nNode,nCrd)
        real(kind=RK), intent(out) :: amass(nBlock,nDofEl,nDofEl)
        integer :: i, j
        real(kind=RK) :: area, am0, S1, S2, elemMass
        real(kind=RK) :: edge(3,3), V1(3), V2(3)

        do i = 1, nCoh
            edge(1,:) = coords(i,2,:)-coords(i,1,:)
            edge(2,:) = coords(i,3,:)-coords(i,1,:)
            edge(3,:) = coords(i,4,:)-coords(i,1,:)

            call vector_cross( edge(1,:),edge(2,:),V1 )
            
            S1 = vector_mod( V1 )
            
            call vector_cross(edge(2,:),edge(3,:),V2)
            
            S2      = vector_mod( V2 )
            area    = 0.5d0*(S1+S2)
            elemMass = m_thickness*area*m_rho
            am0     = elemMass/real(nNdCoh,kind=RK)

            do j=1,nDofCoh
                amass(i,j,j) = am0
            end do
        end do
    end subroutine mass_matrix

    subroutine one_element_cohesive_force(iEleM,nVar,u,coords,var,fcNd)     !!计算一个内聚力单元的 force 这里的force是指力 N
        implicit none
        integer, intent(in) :: iElem                    !!iElem是  第几个element? #todo
        integer, intent(in) :: nVar                     !!nvar是var数组的个数，即svars中值个数
        real(kind=RK),intent(in) :: u(nDofCoh)          !!单元各节点，个自由度方向的变化，如x y z 三方向的位移，若是温度，可能是温度差值? #todo
        real(kind=RK),intent(in) :: coords(nNdCoh,nDim) !!坐标(8,3)
        real(kind=RK),intent(out) :: var(nvar)          !!svars值传给var,即自定义接变量sdv值,  var(1)为element失效标识，为2时表示完全失效

        real(kind=RK),intent(out) :: fcNd(nDim,nNdCoh)  !! =rhs array #bug (3,8) may be not allowed, value range error， 这里传入的是(1,24),但是定义成(3,8)

        real(kind=RK):: normal(3)   !!法向量, The unit normal vector at Gaussian point
        real(kind=RK):: det         !! The length of normal vector， it is also an area of the middle plane
        integer :: fCount           !! #todo 计数删除的积分点，为4时表示所有积分点失效，单元完全失效
        integer :: i                !!计数变量
        integer :: j                !!计数变量  此子程序中暂未用到
        integer :: ibgn             !!计数变量  各积分点是否损伤的变量的下标位置 2, 5, 8, 11
        integer :: iEffGap          !!计数变量  各积分点实际分离量的下标位置 3, 6, 9, 12
        integer :: damage_num       !!计数变量  各积分点 损伤变量值 所在下标位置 4， 7， 10， 13
        real(kind=RK) :: fail
        real(kind=RK) :: failElem  !!失效单元？#todo
        real(kind=RK) :: maxEffGap
        real(kind=RK) :: xi(2)
        real(kind=RK) :: gapPt(3)
        real(kind=RK) :: xyzMid(3,4)    !! The coordinates of the intermediate  points between the upper and lower nodes of the hexadecimal element
        real(kind=RK) :: gapNd(3,4)     !! The vector from the node below the cell to the node above it
        real(kind=RK) :: fcAdd(3,4)
        real(kind=RK) :: fcPt(3)         !!计算传出的 力向量
        real(kind=RK) :: damage_old      !!新计算的损伤变量, #bug 这里的单元损伤变量是整个单元的吗，还是每个节点都有一个损伤变量值

        if( nVar < nInteg2D * 2 + 1 ) goto 111

        failElem = var(1)
        
        if( failElem >= 2 ) return  ! element complete failure

        call seperation_and_middle_surface(u,coords,gapNd,xyzMid)
        fcNd    = 0.0d0
        fcount  = 0
        ibgn    = 2

        do i=1,nInteg2D

            fail = var(ibgn) ! Here ibgn can be 2, 5, 8, 11
            
            !First, determine whether the integration point is in failure.
            if( fail>=2 ) then
                fcount = fcount + 1 !means this gauss point is in failure.
                cycle   !end this loop, begin next loop
            end if

            xi = xiGaussPoint(:,i) !第 i 列

            call get_unit_normal_vector(xi,xyzMid,det,normal)
            call calculate_gap_at_integration_point(xi,gapNd,gapPt)
            !call obtain_failure_values_of_newton_integration_points(
            !co_eid,ith,fail=fail)
            
            iEffGap = ibgn+1  ! Here iEffGap can be 3, 5, 7, 9
            maxEffGap = var(iEffGap)  !The maxEffGap here is the data from the previous step

            damage_num = iEffGap + 1  !4， 7， 10， 13
            damage_old = var(damage_num)

            
            call one_point_cohesive_force(var,nvar,det,fCount,maxEffGap,gapPt,normal,fail,fcPt,damage_old)
            
            var(damage_num) = damage_old !更新积分点 损伤变量的值
            var(iEffGap) = maxEffGap   !The maxEffGap here is the updated data.
            var(ibgn) = fail !更新积分点是否失效
            
            call distribute_force_to_nodes(xi,fcPt,fcAdd)
            
            fcNd(:,1:nHalfNd) = fcNd(:,1:nHalfNd) + fcAdd
            ibgn = ibgn + 3
        end do

        fcNd(:,nHalfNd+1:nNdCoh) = fcNd(:,1:nHalfNd)
        fcNd(:,1:nHalfNd)  = -fcNd(:,1:nHalfNd)

        if( fCount == nInteg2D ) var(1) = 2.01d0          !单元以完全失效，设置失效标识 var(1) = 2.01d0
        return

        111 continue
        msg='Wrong length of svars in one_element_cohesive_force.'
        write(6,*) msg
        call xplb_exit
        return
    end subroutine one_element_cohesive_force

    subroutine distribute_force_to_nodes(xi,fcPt,fcNd)
        implicit none
        real(kind=RK),intent(in)::xi(2), fcPt(3)
        real(kind=RK),intent(out)::fcNd(3,4)
        real(kind=RK)::shap(4)
        integer::nNum, i
        
        fcNd = 0.0d0
        call shape_q4(xi(1),xi(2),shap)
        nNum = nInteg2D
        do i=1,nNum
            fcNd(:,i) = shap(i)*fcPt
        end do
        return
    end subroutine distribute_force_to_nodes

    subroutine calculate_gap_at_integration_point(xi,gapNd,gapPoint)
        implicit none
        real(kind=RK), intent(in) :: xi(2) !!接收 （（不是高斯点）积分）积分点的两个参数
        real(kind=RK), intent(in) :: gapNd(3,4) !! The coordinate difference between the upper and lower nodes at the incoming Gaussian point, &
                                                !! that is, the vector from the lower node to the upper node
        real(kind=RK), intent(out) :: gapPoint(3) !! Output, but it is acutally gapNd, which is the vector from the bottom node at the Gaussian point to the upper onde. 
        real(kind=RK) :: t(4)
        integer :: i

        call shape_q4(xi(1),xi(2),t) !由（不是高斯点）积分点参数计算形函数

        gapPoint = 0.0d0
        do i = 1, nHalfNd
            gapPoint = gapPoint + t(i) * gapNd(:,i) !#bug 这里是获得加权平均的分离量？实际上形函数只在该（不是高斯点）积分点位置为1，其他位置为0,所以加权平均就是该（不是高斯点）积分点位置的分离量        
        end do
        return
    end subroutine calculate_gap_at_integration_point

    subroutine get_unit_normal_vector(xi,xyz,det,normal)  !!计算单位法向向量
        implicit none
        real(kind=RK),intent(in)::xi(2)          !!接收 （不是高斯点）积分点数组的第 i 列
        real(kind=RK),intent(in)::xyz(3,4)       !!接收中间截面节点坐标
        real(kind=RK),intent(out)::normal(3)     !!法向量，最后传出为单位法向量,中间截面处（不是高斯点）积分点的法向量
        real(kind=RK),intent(out)::det           !!传出 法向量的长度
        real(kind=RK):: tvect1(3)
        real(kind=RK):: tvect2(3)
        real(kind=RK):: tmod !!法向量长度
        real(kind=RK),parameter:: tol=1.0e-10

        call calculate_tangential_vector(xyz,xi,tVect1,tVect2)
        call vector_cross(tVect1,tVect2,normal)   !两切向量叉乘即为法向量

        tmod    = vector_mod(normal)            !取模，即计算法向量长度
        if(abs(tmod)<tol) goto 333              !说明计算错误
        det     = tMod
        normal  = normal/tMod                   !将法向量转化为单位法向量
        return

        333 continue
        msg='Negative tmod in get_unit_normal_vector'
        write(6,*) msg
        call xplb_exit
    end subroutine get_unit_normal_vector

    subroutine calculate_tangential_vector(xyz,gpos,tvect1,tvect2)  !!计算切向量
        implicit none
        real(kind=RK),intent(in) :: xyz(3,4) !!接收 中间截面节点坐标
        real(kind=RK),intent(in) :: gpos(2)  !!接收 （不是高斯点）积分点数组 第 i 列
        real(kind=RK),intent(out) :: tvect1(3) !! The vector of this Gaussian point pointing to the adjacent Gaussian point
        real(kind=RK),intent(out) :: tvect2(3) !! Same as above,. pointing to another adjacent Gaussian point
        real(kind=RK)::dshap1(4),dshap2(4)
        integer::id
        
        tvect1=0.0d0; tvect2=0.0d0  !The array is initialized to zero
        call shap_function_derivative(gpos,dshap1,dshap2)

        ! #todo 有限元求解切向量的方法还不理解; 这里计算了该（不是高斯点）积分点的切向量，分别为该点指向相邻两点的的向量
        do id=1,4
          tvect1=tvect1+xyz(:,id)*dshap1(id)
          tvect2=tvect2+xyz(:,id)*dshap2(id)
        end do

        return
    end subroutine calculate_tangential_vector

    subroutine seperation_and_middle_surface(u,coords,gap,xyzMidsurf)     !!计算中间接截面节点坐标，以及上下节点间距向量
        implicit none
        real(kind=RK), intent(in) :: u(nDofCoh)             !!u(24)
        real(kind=RK), intent(in) :: coords(nNdCoh,nDim)    !!coords(8,3)
        real(kind=RK), intent(out) :: gap(nDim,nInteg2D)    !!(3,4)二维数组， 由下方节点指向上方节点的向量
        real(kind=RK), intent(out) :: xyzMidsurf(:,:)       !!6面体单元变形后上下对应两节点之间 中间点 的坐标
        real(kind=RK) :: disNode(nDim,nNdCoh)               !!(3,8)二维数组
        integer :: k

        disNode = reshape( u, (/nDim,nCoh/) )  ! nDim = 3, so nCoh = 8, (3,8)
        !iBgn = 1
        !do k=1,nnode
        !   disNode(:,k)=u(iBgn:iBgn+2)
        !   iBgn = iBgn + 3
        !end do
        do k = 1, nHalfNd  ! Here the value of k ranges from 1 to 4, but coords
           xyzMidsurf(:,k) = ( coords(k,:)+coords(k+nInteg2D,:)+ disNode(:,k)+disNode(:,k+nInteg2D))/2.0d0  !原来相对应两个节点的坐标值加位移量 再 除以2，得到两节点中间点坐标
           gap(:,k)=disNode(:,k+nInteg2D)-disNode(:,k)  !上下对应节点坐标值相减，为上下节点间距向量，为下方节点指向上方节点的向量
        end do
        return
    end subroutine seperation_and_middle_surface

    subroutine one_point_cohesive_force(var,nvar,area,fCount,ghis,gapVect,normal,fail,fc,damage)  !! gpos, Hpos, fcount
        implicit none
        real(kind=RK),intent(in) :: var(nvar)
        integer, intent(in) :: nvar
        real(kind=RK),intent(in):: area         !! the area of middle surface
        real(kind=RK),intent(in):: gapVect(3)   !! It is acutally gapNd, which is the vector from the bottom node at the Gaussian point to the upper onde.
        real(kind=RK),intent(in):: normal(3)    !! The unit normal vector at the Gaussian point.
        integer, intent(inout) :: fCount 
        real(kind=RK),intent(inout) :: ghis     !!实际分离量
        !real(kind=RK),intent(in) :: gpos(2)     !!
        !real(kind=RK),intent(in) :: Hpos(2)     !!
        real(kind=RK),intent(out) :: fail
        real(kind=RK),intent(out) :: fc(3)       !!最后传出的力向量
        real(kind=RK),intent(inout) :: damage    !!传入的旧损伤变量和新计算的损伤变量比较进行更新
        !integer,intent(inout)::fcount
        real(kind=RK) :: tmod
        real(kind=RK) :: nGap               !! 下方节点指向上方节点的向量 与 积分点单位法向量 的 点乘，所以可能会出现 负值
        real(kind=RK) :: tGap               !!上下两节点的切向量
        real(kind=RK) :: co_sig(3)          !!计算得到的应力向量， 法向应力向量 + 切向应力向量
        real(kind=RK) :: nGapVect(3)        !!normal separation vector
        real(kind=RK) :: tGapVect(3)        !!通过向量首尾相连的理论，获得切向量
        real(kind=RK) :: effDelta           !!合成的总分离量
        real(kind=RK) :: tfc(3)             !!临时储存力的值
        real(kind=RK) :: tvect(3)           !!计算得出单位切向量
        real(kind=RK) :: deltai             !!开始损伤分离量
        real(kind=RK) :: deltaf             !!完全失效分离量
        real(kind=RK) :: gVect(3)           !! =gapVect, the vector from the lower node to the upper node at the Gaussian point.

        fc      = 0.0d0
        gVect   = gapVect
        nGap    = gVect(1)*normal(1) + gVect(2)*normal(2) + gVect(3)*normal(3) ! 法向分离量（间隔量？） The component of gVect in the direction of the unit normal vector of the Gaussian point
        nGapVect = nGap * normal          !normal separation vector
        tGapVect = gVect - nGapVect     !通过向量首尾相连的理论，获得切向量
        tGap = sqrt( tGapVect(1)*tGapVect(1) + tGapVect(2)*tGapVect(2) + tGapVect(3)*tGapVect(3) )   !得到上下两个节点间的切向分离量

        if(nGap < 0.0) then               !单位法向量与下节点指向上节点的向量方向不同，角度超过90度， 发生了穿透所以用 罚函数法 限制穿透
            tfc = -abs(nGap) * m_stiff * area * normal  ! 这里计算出来的是什么？ --罚函数法 得到的法向力
            fc  = fc + tfc
        end if

        call calculate_effective_displacement_point(nGap,tGap,effDelta)  !计算总分离量，不是法向或切向，是合成的总分离量
        if(effDelta <= 0.0d0) return

        call compute_input_coh_displacement(m_stiff,m_nStrength,m_tStrength,m_Gnc,m_Gtc,m_alpha,m_deltaNo,m_deltaSo,tGap,nGap,deltai,deltaf)

        if(effDelta >= deltaf) then !failure of the gauss point.
            fCount = fCount + 1
            ghis   = effDelta
            !m_BindDem(iDem)%status(iFacet) = 2
            fc     = 0.0d0      !该积分点已完全失效，故不承力
            fail   = 2.01d0
            damage = 1.0d0      !该积分点损伤变量值为1，即完全损伤
        else
            !m_BindDem(iDem)%status(iFacet) = 1
            if(ghis <= effDelta)then
                ghis = effDelta  !更新实际分离量
                !m_bindDem(iDem)%maxEffGap(iFacet) = ghis
            end if

            call calculate_tangent_vector(tGapVect,tvect)  !!计算出单位切向量
            call cohesive_stress_at_one_point(m_stiff,deltai,deltaf,gVect,nGap,tGap,effDelta,ghis,normal,tvect,fail,co_sig,damage)
            tfc = co_sig * area
            fc = fc + tfc
        end if
        return

        111 continue
        msg='mistake of the value of effDelta. in one_point_bind_force'
        write(6,*) msg
        call xplb_exit
        222 continue
        msg='mistake of the value of ghis. in one_point_bind_force'
        write(6,*) msg
        call xplb_exit
    end subroutine one_point_cohesive_force
end module cohesive_element_mod

!================================================VUEL=========================================================
    subroutine vuel(nblock,rhs,amass,dTimeStable,svars,nsvars,energy,nnode,ndofel,props,nprops,jprops,njprops,coords,ncrd,u,du,v,a,jtype,jelem,time,period,dtimeCur,dtimePrev,kstep,kinc,lflags,dMassScaleFactor,predef,npredef,jdltyp,adlmag)

!       used module
        use system_constant
        use math_mod
        use parameter_mod
        use cohesive_element_mod

!       if you want to use `implicit none`, you must add the folowing three lines of code
        implicit none
        integer, parameter :: j_sys_Dimension = 2
        integer, parameter :: maxblk = 512
        
!       operation code
        integer, parameter :: jMassCalc = 1
        integer, parameter :: jIntForceAndDtStable = 2
        integer, parameter :: jExternForce         = 3

!       flags
        integer, parameter :: iProcedure = 1
        integer, parameter ::  iNlgeom    = 2
        integer, parameter :: iOpCode    = 3
        integer, parameter :: nFlags     = 3

!       time
        integer, parameter :: iStepTime  = 1
        integer, parameter :: iTotalTime = 2
        integer, parameter :: nTime      = 2

!       procedure flags
        integer, parameter ::  jDynExplicit = 17

!       predefined variables
        integer, parameter :: iPredValueNew = 1
        integer, parameter :: iPredValueOld = 2
        integer, parameter :: nPred         = 2 

        integer, parameter :: nElEnergy = 12

        !integer :: tempRead   !! debug flag

        integer, intent(in) :: jType
        integer, intent(in) :: jdlTyp   !! integer identifying the load type
        integer, intent(in) :: kStep    !! Current step number
        integer, intent(in) :: kinc     !! Current increment number
        integer, intent(in) :: nPreDef  !! field variable number   not used now  
        integer, intent(in) :: nBlock   !! total solid element number
        integer, intent(in) :: nCrd     !! n deminsion
        integer, intent(in) :: nDofEl   !! total dof number of one element = 24
        integer, intent(in) :: njprops  !! user defined properties (integer)
        integer, intent(in) :: nNode    !! total node number of one element  (real)
        integer, intent(in) :: nProps   !! user defined properties 
        integer, intent(in) :: nsvars   !! state variable number
        integer, intent(in) :: jElem(nBlock)    !! the element ID
        integer, intent(in) :: jProps(njProps)  !! User-defined number of integer property values  
        integer, intent(in) :: lFlags(nFlags)
        integer, intent(in) :: dTimeStable(nFlags) !! stable time increment of each element
        real(kind=RK), intent(in) :: dTimeCur   !! Current time increment
        real(kind=RK), intent(in) :: dTimePrev  !! Previous time increment.
        real(kind=RK), intent(in) :: period     !! Time period of the current step
        real(kind=RK), intent(in) :: props(nProps) !! User-defined number of real property values associated with the elements processed
        real(kind=RK), intent(in) :: coords(nBlock,nNode,nCrd) !! inital coordinates     
        real(kind=RK), intent(in) :: u(nBlock,nDofEl) !! total displacements
        real(kind=RK), intent(in) :: du(nBlock,nDofEl) !! Incremental displacements in the current increment
        real(kind=RK), intent(in) :: v(nBlock,nDofEl)  !! velocities at the midpoint of the increment
        real(kind=RK), intent(in) :: a(nBlock, nDofEl) !! Accelerations at the end of the current increment
        real(kind=RK), intent(in) :: adlmag(nBlock) !! the load magnitude of the load type jdltyp for each elem 
        real(kind=RK), intent(in) :: time(iTotalTime)    !! Current value of total time
        real(kind=RK), intent(in) :: preDef(nBlock, nNode, nPreDef, nPred)  !! predefined field variable
        real(kind=RK), intent(in) :: dMassScaleFactor(nBlock)  !! the mass scale factors for each element
        real(kind=RK), intent(out) :: rhs(nBlock,nDofEl) !! the contributions of each element to the right-hand-side vector of the overall system of equations (e.g.internal force)
        real(kind=RK), intent(out) :: amass(nBlock,nDofEl,nDofEl)  !!  mass matrix of each element
        !real(kind=RK), intent(out) :: dTimeStable(nBlock) !! stable time increment of each element
        real(kind=RK), intent(out) :: svars(nBlock,nsvars) !! the solution-dependent state variables
        real(kind=RK), intent(in) ::  energy(nBlock,nElEnergy)  !! not computed

        integer :: i

        if( jtype==2 .and. lflags(iProcedure)==jDynExplicit ) then 

            !write(*,*) "Please input an integer:"
            !read(*,*) tempRead
            !call sleep(5)

            call set_parameter(nBlock,nNode,nCrd,nDofEl,nProps,props)

            if( lflags(iOpCode) == jMassCalc )then  
                call mass_matrix(nBlock,nNode,nCrd,nDofEl,coords,amass) 
            
            else if ( lflags(iOpCode) == jIntForceAndDtStable) then  
                do i = 1, nblock
                    call one_element_cohesive_force(i,nsvars,u(i,:),coords(i,:,:),svars(i,:),rhs(i,:))
                    write(6,*) 'node force'
                    write(6,*) rhs(i,1:12)
                    write(6,*) rhs(i,13:24)
                    if( svars(1,2)>0.0 )then
                        write(6,*) 'time', time(1)
                        write(6,*) 'svars'
                        write(6,*) svars(1,:)
                        write(6,*) 'fc'
                        write(6,*) rhs(1,1:12)
                        !call xplb_exit
                    end if
                end do
            end if
        end if
        return
    end subroutine vuel