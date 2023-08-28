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
      end module

      module math_mod                                                     !! 基本的数学函数和形函数的运算等
        include 'vaba_param.inc'
        private 
        public :: vector_cross, vector_dot, vector_mod, vector_unit
        public :: shape_q4
      contains
       subroutine vector_cross( vec1, vec2 , cross) 
         !
         ! Function:
         !
         !  Calculate the cross product of two vectors.
         !
         ! Remarks:
         !
         !   RESULT = VEC1 × VEC2
         !
         ! Parameters:
         !
         !  Vec1, input, an array of rank one with 3 elements.
         !  Vec2, input, an array of rank one with 3 elements.
         !  result, output, an array of rank one with 3 elements.
         !
         !! use system_constant
         include 'vaba_param.inc'
         dimension :: vec1(3), vec2(3)
         dimension :: cross(3)
         cross(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
         cross(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
         cross(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)!叉积结果向量三个值simon 2011-11-1
         return
     
        end subroutine vector_cross
       
        subroutine vector_dot( vec1, vec2, dot) 
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

          include 'vaba_param.inc'
          dimension :: vec1(3), vec2(3)
          dot = vec1(1)*vec2(1) + vec1(2)*vec2(2) + vec1(3)*vec2(3)
          
          return
         end subroutine vector_dot !点积结果
        
         function vector_mod( vec ) result( rmod )
         include 'vaba_param.inc'
         dimension :: vec(3)
         rmod = sqrt( vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3) )
         return
         end Function vector_mod 

         subroutine vector_unit( vec )
         include 'vaba_param.inc'
         dimension :: vec(3)
         parameter( TINY = 1.0e-25 )

         rmod = sqrt( vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3) )
         if( rmod<TINY)then
           write(6,*)'Error in vector_unit:the rmod of the vector is 0.'
           stop
         end if
         vec = vec/rmod
         return
         end subroutine vector_unit 
         
      subroutine shape_q4 ( r, s, t, dtdr, dtds )
c*****************************************************************************80
      !
      !! SHAPE_Q4 evaluates shape functions for a 4 node reference quadrilateral.
      !
      !  Reference Element Q4:
      !
      !    |
      !    1  4-----3
      !    |  |     |
      !    |  |     |
      !    S  |     |
      !    |  |     |
      !    |  |     |
      !   -1  1-----2
      !    |
      !    +--1--R--1-->
      !
      !
      !  Parameters:
      !
      !    Input : R, S, the reference coordinates of a point.
      !
      !    Output: T(4), the basis functions at the point.
      !
      !    Output: DTDR(4), the R basis derivatives at the point.
      !
      !    Output: DTDS(4), the S basis derivatives at the point.
      !
      !! use system_constant   !use system_constant必须在implicit none前面
      include 'vaba_param.inc'
      dimension ::  dtdr(4), dtds(4), t(4)
      parameter ( quarter = 0.25d0 )
      parameter ( one = 1.d0 )
      t(1) = quarter * ( one - r ) * ( one - s )
      t(2) = quarter * ( one + r ) * ( one - s )
      t(3) = quarter * ( one + r ) * ( one + s )
      t(4) = quarter * ( one - r ) * ( one + s )     

      dtdr(1) = quarter * (- one + s )
      dtdr(2) = quarter * (  one - s )
      dtdr(3) = quarter * (  one + s )
      dtdr(4) = quarter * (- one - s )

      dtds(1) = quarter * (- one + r )
      dtds(2) = quarter * (- one - r )
      dtds(3) = quarter * (  one + r )
      dtds(4) = quarter * (  one - r )

      return
      end subroutine shape_q4
      end module

      module cohesive_law_mod                                                                                                                                                                                                                                         
      use system_constant                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
      implicit none                                                                                                                                                                                                                                                  
      private                                                                                                                                                                                                                                                        
      public :: compute_input_coh_displacement                                                                                                                                                                                                                
      public :: cohesive_stress_at_one_point
      public :: calculate_effective_displacement_point                                                                                                                                                                                 
      public :: calculate_tangent_vector
      character(len=MNL) :: msg                                                                                                                                                                                                                                      
      contains                                                                                                                                                                                                                                                       
      subroutine one_point_cohesive_force(iDem,iFacet,gapVect,
     *           normal,area,fc)  !! gpos, Hpos, fcount                                                                                                                                                                
      implicit none                                                                                                                                                                                                                                               
      integer,intent(in)::iDem, iFacet                                                                                                                                                                                                                            
      real(kind=RK),intent(in):: gapVect(3)   !! gapVect是该高斯积分点上的分离量向量                                                                                                                                                                              
      real(kind=RK),intent(in):: normal(3)    !! normal是中间面上该高斯积分点上的单位法向量                                                                                                                                                                       
      !real(kind=RK),intent(in):: gpos(2)     !! gpos 高斯积分点的自然坐标                                                                                                                                                                                        
      real(kind=RK),intent(in):: area         !! 高斯积分点对应的雅可比行列式值，也就是面积                                                                                                                                                                       
      !real(kind=RK),intent(in):: Hpos(2)     !! Hpos是高斯积分点的权系数                                                                                                                                                                                         
      real(kind=RK),intent(out)::fc(3)      !!  fc 该积分点算出来全局坐标系下的cohesive force                                                                                                                                                                     
      !integer,intent(inout)::fcount                                                                                                                                                                                                                              
      real(kind=RK)::tmod,nGap,ghis,tGap                                                                                                                                                                                                                          
      real(kind=RK)::co_sig(3),nGapVect(3),tGapVect(3),effDelta                                                                                                                                                                                                   
      real(kind=RK)::tfc(3),tvect(3)                                                                                                                                                                                                                              
      real(kind=RK)::deltai,deltaf,idis(3), gVect(3)                                                                                                                                                                                                              
      fc=0.0                                                                                                                                                                                                                                                      
      gVect = gapVect                                                                                                                                                                                                                                             
      nGap=gVect(1)*normal(1)+gVect(2)*normal(2)+gVect(3)*normal(3) !! 法向分离值                                                                                                                                                                               
      nGapVect = nGap*normal          !! 法向分离向量 normal separation vector                                                                                                                                                                                    
      tGapVect = gVect - nGapVect     !! 切向分离向量                                                                                                                                                                                                             
      tGap = sqrt( tGapVect(1)*tGapVect(1)+tGapVect(2)*tGapVect(2)+tGapVect(3)*tGapVect(3) ) !! 切向分离值                                                                                                                                                        
      if(nGap<0.0)then !受压状态                                                                                                                                                                                                                                  
         tfc=-abs(nGap)*m_stiff*area*normal                                                                                                                                                                                                                      
         fc=fc+tfc                                                                                                                                                                                                                                               
      end if                                                                                                                                                                                                                                                      
      call calculate_effective_displacement_point(nGap,tGap,effDelta)  !! 计算有效张开量                                                                                                                                                                          
      if(effDelta==0.0) return                                                                                                                                                                                                                                    
      call compute_input_coh_displacement(m_stiff, 
     *     m_nStrength,m_tStrength,m_GIc,m_GIIc,m_alpha,m_deltaNo,
     *     m_deltaSo,tGap,nGap,deltai,deltaf) !! 计算混合三角形的两个位移                                                                                                                                                     
      if(effDelta>=deltaf)then !failure of the gauss point.                                                                                                                                                                                                       
        !fcount=fcount+1                                                                                                                                                                                                                                        
        !call update_failure_value_of_gauss_integration_points(co_eid,ith)                                                                                                                                                                                      
        m_BindDem(iDem)%status(iFacet) = 2                                                                                                                                                                                                                      
      else                                                                                                                                                                                                                                                        
        !call obtain_failure_values_of_gauss_integration_points(co_eid,ith,ghis=ghis) !! ghis有效张开量的历史最大值                                                                                                                                             
        ghis = m_bindDem(iDem)%maxEffGap(iFacet)                                                                                                                                                                                                                
        m_BindDem(iDem)%status(iFacet) = 1                                                                                                                                                                                                                      
        if(ghis<=effDelta)then                                                                                                                                                                                                                                  
          ghis = effDelta                                                                                                                                                                                                                                     
          m_bindDem(iDem)%maxEffGap(iFacet) = ghis                                                                                                                                                                                                            
          !call update_failure_delta_of_gauss_integration_points(co_eid,ith,ghis=ghis)                                                                                                                                                                        
        end if                                                                                                                                                                                                                                                  
        call calculate_tangent_vector(tGapVect,tvect)  !! 计算切向单位向量                                                                                                                                                                                      
        call cohesive_stress_at_one_point(m_stiff,deltai,deltaf,gVect,
     *       nGap,tGap,effDelta,ghis,normal,tvect,co_sig)                                                                                                                                                                      
        tfc = co_sig*area                                                                                                                                                                                                                                       
        fc = fc + tfc                                                                                                                                                                                                                                         
      end if                                                                                                                                                                                                                                                      
      return                                                                                                                                                                                                                                                      
      111 continue                                                                                                                                                                                                                                                    
        msg='mistake of the value of effDelta. in one_point_bind_force'                                                                                                                                                                                             
            !call errormsg(msg)                                                                                                                                                                                                                                          
        stop                                                                                                                                                                                                                                                      
      222 continue                                                                                                                                                                                                                                                    
        msg='mistake of the value of ghis. in one_point_bind_force'                                                                                                                                                                                                 
        !call errormsg(msg)                                                                                                                                                                                                                                          
        stop                                                                                                                                                                                                                                                      
      end subroutine
                                                                                                                                                                                                                                                                        
      subroutine compute_input_coh_displacement(stiff,nStreng,
     * tStreng,GN,GS,alfa,deltaNo,deltaSo,tdelta,ndelta,deltai,deltaf)                                                                                                                                                                                           
      implicit none                                                                                                                                                                                                                                               
      real(kind=RK),intent(in)::  stiff, nStreng, tStreng                                                                                                                                                                                                         
      real(kind=RK),intent(in):: GN, GS !! GIC, GIIC                                                                                                                                                                                                              
      real(kind=RK),intent(in):: alfa                                                                                                                                                                                                                             
      real(kind=RK),intent(in):: deltaNo,deltaSo !! deltaI0, deltaII0                                                                                                                                                                                             
      real(kind=RK),intent(in):: tdelta,ndelta                                                                                                                                                                                                                    
      real(kind=RK),intent(out):: deltai,deltaf                                                                                                                                                                                                                   
      real(kind=RK)::beita,nsig,tsig                                                                                                                                                                                                                              
      real(kind=RK)::gna,gta,beita2,alfam,rnsig                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                        
      deltai=0.0; deltaf=0.0                                                                                                                                                                                                                                      
      call calculate_mixed_mode_rate(tdelta,ndelta,beita)                                                                                                                                                                                         
      beita2=beita*beita                                                                                                                                                                                                                                          
      !beita2=0.7                                                                                                                                                                                                                                                 
      !for reference: Pinho S T, Iannucci L, Robinson P. Formulation and implementation 
      !of decohesion elements in an explicit finite element code[J].                                                                                                             
      !Composites Part A: Applied science and manufacturing, 2006, 37(5): 778-789.                                                                                                                                                                                
      !call obtain_transition_value_in_intrinsic_model(stiff,nStrengh,tStrength,delta_no,delta_so)                                                                                                                                                            
      if(ndelta<=0.0)then                                                                                                                                                                                                                                     
         deltai=deltaSo !effective displacement at transition point.                                                                                                                                                                                         
         if(deltai<=0.0) goto 222                                                                                                                                                                                                                            
         if(tStreng<=0.0) goto 111                                                                                                                                                                                                                           
         deltaf=2*Gs/tStreng !final displacement of the mixed mode.  !! Gs: G2c; Gn: G1c                                                                                                                                                                     
      else                                                                                                                                                                                                                                                    
         !effective displacement at transition point of the mixed triangle.                                                                                                                                                                                  
         deltai=deltaNo*deltaSo*sqrt( (1+beita2)/
     *          ( beita2*deltaNo*deltaNo+deltaSo*deltaSo) )                                                                                                                                                                    
         if(deltai<=0.0) goto 222                                                                                                                                                                                                                            
         if(stiff<=0.0) goto 333                                                                                                                                                                                                                             
         deltaf=2*(1+beita2)/(stiff*deltai)                                                                                                                                                                                                                  
         if(GN<=0.0.or.GS<=0.0) goto 444                                                                                                                                                                                                                     
            gna=(1/GN)**alfa; gta=(beita2/GS)**alfa                                                                                                                                                                                                             
            alfam=-1/alfa                                                                                                                                                                                                                                       
            deltaf=deltaf*(gna+gta)**alfam !final displacement of the mixed mode.                                                                                                                                                                               
      end if                                                                                                                                                                                                                                                  
      return                                                                                                                                                                                                                                                      
      111 continue                                                                                                                                                                                                                                                    
        msg='Negative tstrength in compute_input_coh_displacement.'                                                                                                                                                                             
            !call errormsg(msg)                                                                                                                                                                                                                                          
        stop                                                                                                                                                                                                                                                      
      222 continue                                                                                                                                                                                                                                                    
        msg='Negative deltai in compute_input_coh_displacement.'                                                                                                                                                                              
            !call errormsg(msg)                                                                                                                                                                                                                                          
        stop                                                                                                                                                                                                                                                      
      333 continue                                                                                                                                                                                                                                                    
        msg='Negative stiffness in compute_input_coh_displacement.'                                                                                                                                                                        
        !call errormsg(msg)                                                                                                                                                                                                                                          
         stop                                                                                                                                                                                                                                                      
      444 continue                                                                                                                                                                                                                                                    
        msg='Negative GN or GS in compute_input_coh_displacement.'                                                                                                                                                                            
        !call errormsg(msg)                                                                                                                                                                                                                                          
        stop                                                                                                                                                                                                                                                      
      end subroutine                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                        
      subroutine compress_force(stiff, depth,normal,area,pfac,facc)                                                                                                                                                                                                  
      implicit none                                                                                                                                                                                                                                               
      real(kind=RK),intent(in):: stiff                                                                                                                                                                                                                            
      real(kind=RK),intent(in):: depth,normal(3),area                                                                                                                                                                                                             
      real(kind=RK),intent(out):: facc(3)                                                                                                                                                                                                                         
      real(kind=RK),intent(in):: pfac                                                                                                                                                                                                                             
      real(kind=RK)::det,norm(3)                                                                                                                                                                                                                                  
      real(kind=RK)::ustiff=0.0                                                                                                                                                                                                                                   
       facc=0.0                                                                                                                                                                                                                                                    
            !ustiff=intrin_stiff                                                                                                                                                                                                                                        
       ustiff=pfac*stiff                                                                                                                                                                                                                                          
       facc=abs(depth)*ustiff*normal                                                                                                                                                                                                                              
      return                                                                                                                                                                                                                                                      
      end subroutine                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                        
      subroutine cohesive_stress_at_one_point(stiff,deltai,deltaf,
     * delta,det_ndelta,det_tdelta,eff_delta,ghis,nvect,tvect,sig)                                                                                                                                                                                
      implicit none                                                                                                                                                                                                                                                
      real(kind=RK),intent(in):: stiff                                                                                                                                                                                                                             
      real(kind=RK),intent(in)::delta(3),det_ndelta,det_tdelta,ghis,eff_delta                                                                                                                                                                                      
      real(kind=RK),intent(in)::nvect(3),tvect(3),deltai,deltaf                                                                                                                                                                                                    
      !integer,intent(in)::co_eid, gith                                                                                                                                                                                                                            
      real(kind=RK),intent(out)::sig(3)                                                                                                                                                                                                                            
      real(kind=RK)::nsig(3),tsig(3)                                                                                                                                                                                                                               
      real(kind=RK)::dnsig,dtsig,det,temp                                                                                                                                                                                                                          
      real(kind=RK)::eta,modu,bvect(3),delta_tc                                                                                                                                                                                                                    
      real(kind=RK)::rnvect(3),rtvect(3)                                                                                                                                                                                                                           
      real(kind=RK),parameter::tol=1.0e-15                                                                                                                                                                                                                         
      real(kind=RK)::nsig0, tsig0, effsig0 !normal and tangential stress components at the transition point in a TSL.                                                                                                                                              
      real(kind=RK)::coef=0.0                                                                                                                                                                                                                                      
                                                                                                                                                                                                                                                                        
      sig=0.0; dnsig=0.0; dtsig=0.0; coef=0.0                                                                                                                                                                                                                      
      !if(interface_failure/=1) then                                                                                                                                                                                                                               
      if(ghis<=0.0) goto 111  !! !! ghis有效张开量的历史最大值                                                                                                                                                                                                 
      !end if                                                                                                                                                                                                                                                      
                                                                                                                                                                                                                                                                        
      !for reference: Wu L, Sket F, Molina-Aldareguia J M, et al. A study 
      ! of composite laminates failure using an anisotropic gradient-enhanced 
      !damage mean-field homogenization model[J]. Composite Structures, 2015, 126: 246-264.                                                                                                                                       
      !for reference: Hu N, Zemba Y, Okabe T, et al. A new cohesive model for simulating delamination propagation                                                                                                                                                  
      !in composite laminates under transverse loads[J]. Mechanics of Materials, 2008, 40(11): 920-935.                                                                                                                                                            
      !if(interface_failure==1)then   !intrinsic cohesive modeling                                                                                                                                                                                                 
      if(ghis<=deltai)then                                                                                                                                                                                                                                     
        if(det_ndelta>0.0) dnsig=stiff*det_ndelta                                                                                                                                                                                                            
        dtsig=stiff*det_tdelta                                                                                                                                                                                                                               
      else if(ghis>deltai.and.ghis<=deltaf)then                                                                                                                                                                                                                
        coef=deltaf/ghis; coef=coef*(ghis-deltai)/(deltaf-deltai)                                                                                                                                                                                            
        if(det_ndelta>0.0) dnsig=(1-coef)*stiff*det_ndelta                                                                                                                                                                                                   
        dtsig=(1.0-coef)*stiff*det_tdelta                                                                                                                                                                                                                    
      end if                                                                                                                                                                                                                                                   
      !else                                                                                                                                                                                                                                                        
      !    goto 555                                                                                                                                                                                                                                                
      !end if                                                                                                                                                                                                                                                      
      modu=sqrt(delta(1)*delta(1)+delta(2)*delta(2)+delta(3)*delta(3))                                                                                                                                                                                            
      if(modu==0.0) goto 222                                                                                                                                                                                                                                      
         bvect=delta                                                                                                                                                                                                                                                 
         bvect=bvect/modu                                                                                                                                                                                                                                            
         rnvect=nvect; rtvect=tvect                                                                                                                                                                                                                                  
         if(det_ndelta>0.0) then                                                                                                                                                                                                                                     
           det=vector_dot(nvect,bvect)                                                                                                                                                                                                                             
         end if                                                                                                                                                                                                                                                      
         det=vector_dot(tvect,bvect)                                                                                                                                                                                                                                 
         if(det_ndelta>0.0) then                                                                                                                                                                                                                                     
            nsig=dnsig*rnvect                                                                                                                                                                                                                                       
         else                                                                                                                                                                                                                                                        
            nsig=0.0                                                                                                                                                                                                                                                
         end if                                                                                                                                                                                                                                                      
         tsig=dtsig*rtvect                                                                                                                                                                                                                                           
         sig=nsig+tsig                                                                                                                                                                                                                                               
      return                                                                                                                                                                                                                                                      
      111 continue                                                                                                                                                                                                                                                    
          msg='Negative ghis in cohesive_stress_at_one_point.'                                                                                                                                                                                         
          !call errormsg(msg)                                                                                                                                                                                                                                          
          stop                                                                                                                                                                                                                                                      
      222 continue                                                                                                                                                                                                                                                    
          msg='Zero mod in cohesive_stress_at_one_point.'                                                                                                                                                                                          
          !call errormsg(msg)                                                                                                                                                                                                                                          
          stop                                                                                                                                                                                                                                                      
      333 continue                                                                                                                                                                                                                                                    
          msg='mistake of the value of det. in cohesive_stress_at_one_point'                                                                                                                                                                                          
          !call errormsg(msg)                                                                                                                                                                                                                                          
          stop                                                                                                                                                                                                                                                      
      444 continue                                                                                                                                                                                                                                                    
          msg='mistake of the value of e_strength. in cohesive_stress_at_one_point'                                                                                                                                                                                   
          !call errormsg(msg)                                                                                                                                                                                                                                          
          stop                                                                                                                                                                                                                                                      
      555 continue                                                                                                                                                                                                                                                    
          msg='mistake of the value. in cohesive_stress_at_one_point'                                                                                                                                                                                                 
          !call errormsg(msg)                                                                                                                                                                                                                                          
          stop                                                                                                                                                                                                                                                      
      end subroutine                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                        
      subroutine calculate_tangent_vector(tdelta,tvect)                                                                                                                                                                                                              
      implicit none                                                                                                                                                                                                                                               
      real(kind=RK),intent(in)::tdelta(3)                                                                                                                                                                                                                         
      real(kind=RK),intent(out)::tvect(3)                                                                                                                                                                                                                         
      real(kind=RK)::det                                                                                                                                                                                                                                          
      real(kind=RK),parameter::tol=1.0e-16                                                                                                                                                                                                                        
       det=sqrt(tdelta(1)*tdelta(1)+tdelta(2)*tdelta(2)+
     *     tdelta(3)*tdelta(3))                                                                                                                                                                                       
       if(det<tol)then                                                                                                                                                                                                                                             
         tvect=(/0,0,0/)                                                                                                                                                                                                                                         
         return                                                                                                                                                                                                                                                  
       end if                                                                                                                                                                                                                                                      
       tvect=tdelta/det                                                                                                                                                                                                                                            
       return                                                                                                                                                                                                                                                      
      end subroutine                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                        
      subroutine obtain_transition_value(stiff,nstreng,tstreng, 
     *   delta_no,delta_so)                                                                                                                                                              
      implicit none                                                                                                                                                                                                                                               
      real(kind=RK),intent(in) :: stiff                                                                                                                                                                                                                           
      real(kind=RK),intent(in) :: nstreng, tstreng                                                                                                                                                                                                                
      real(kind=RK),intent(out)::delta_no,delta_so                                                                                                                                                                                                                
      real(kind=RK),save::no=0.0,so=0.0                                                                                                                                                                                                                           
      integer,save::count=0                                                                                                                                                                                                                                       
       if(count>0)then                                                                                                                                                                                                                                             
         delta_no=no; delta_so=so                                                                                                                                                                                                                                
         return                                                                                                                                                                                                                                                  
       end if                                                                                                                                                                                                                                                      
       count=count+1                                                                                                                                                                                                                                               
       if(stiff<=0.0) goto 111                                                                                                                                                                                                                                     
       no=nstreng/stiff !opening displacement at transition point in pure mode I (intrinsic cohesive model)                                                                                                                                                        
       so=tstreng/stiff !opening displacement at transition point in pure mode II (intrinsic cohesive model)                                                                                                                                                       
       delta_no=no; delta_so=so                                                                                                                                                                                                                                    
       return                                                                                                                                                                                                                                                      
       111 continue                                                                                                                                                                                                                                                    
       msg='Nagetive stiffness in obtain_transition_value.'                                                                                                                                                                   
       !call errormsg(msg)                                                                                                                                                                                                                                          
       stop                                                                                                                                                                                                                                                      
      end subroutine                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                        
      subroutine calculate_effective_displacement_point(det_ndelta,
     *  det_tdelta,eff_delta)                                                                                                                                                                              
      implicit none                                                                                                                                                                                                                                               
      real(kind=RK),INTENT(in)::det_ndelta,det_tdelta                                                                                                                                                                                                             
      real(kind=RK),intent(out)::eff_delta                                                                                                                                                                                                                        
      real(kind=RK)::ndel,tdel                                                                                                                                                                                                                                    
       ndel=(det_ndelta+abs(det_ndelta))*0.5                                                                                                                                                                                                                       
       tdel=det_tdelta                                                                                                                                                                                                                                             
       eff_delta=sqrt(ndel*ndel+tdel*tdel)                                                                                                                                                                                                                         
       return                                                                                                                                                                                                                                                      
      end subroutine                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                        
      subroutine calculate_mixed_mode_rate(tdelta,ndelta,beita)                                                                                                                                                                                       
      implicit none                                                                                                                                                                                                                                               
      real(kind=RK),intent(in)::tdelta,ndelta                                                                                                                                                                                                                     
      real(kind=RK),intent(out)::beita                                                                                                                                                                                                                            
       beita=0.0                                                                                                                                                                                                                                                   
       if(ndelta<=0.0)then                                                                                                                                                                                                                                         
         beita=1.0e20                                                                                                                                                                                                                                            
       else                                                                                                                                                                                                                                                        
         beita=tdelta/ndelta                                                                                                                                                                                                                                     
       end if                                                                                                                                                                                                                                                      
       return                                                                                                                                                                                                                                                      
      end subroutine                                                                                                                                                                                                                                                  
      end module                                                                                                                                                                                                                                                      

      subroutine mass_matrix(nblock,rho,coords,ndofel,amass)          !! 计算质量矩阵
         use math_mod
         include 'vaba_param.inc'
         parameter ( eighth=0.125d0, thickness=1.d0)
         parameter ( nnode=8, ncrd=3)
         dimension amass( nblock,ndofel,ndofel )
         dimension coords( nblock,nnode,ncrd )
         dimension edge(3,3),V1(3),V2(3)
 
         do kblock = 1, nblock
             edge(1,:)=coords(kblock,2,:)-coords(kblock,1,:)
             edge(2,:)=coords(kblock,3,:)-coords(kblock,1,:)
             edged(3,:)=coords(kblock,4,:)-coords(kblock,1,:)
             call vector_cross( edge(1,:),edge(2,:),V1 )  
             S1= vector_mod( V1 )
             call vector_cross(edge(2,:),edge(3,:),V2)
             S2= vector_mod( V2 )
             area=0.5d0*(S1+S2)
             !thickness=one
             !Volume=thickness*S
             elementMass=thickness*area*rho
             am0=eighth*elementMass
             do i=1,ndofel
                 amass(kblock,i,i) = am0
             end do
         end do
      end subroutine  
      
      subroutine key_parameter(props,constant,rkn,rkt,strenthN,       !! 输入需要的参数
     * strenthT,strength,Gnc,Gtc,rho,exp,rknt, delta10,delta20,
     * delta1f,delta2f,gaussPara)
       include 'vaba_param.inc'
       !intent(out) :: props,constant,rkn,rkt,strenthN 
       parameter ( two=2.d0 )
       dimension props(8),rknt(2),strength(2),gaussPara(2,4)
         
         rkn = props(1)
         rkt = props(2)
         strenthN = props(3)
         strenthT = props(4)
         Gnc = props(5)
         Gtc = props(6)
         rho = props(7)
         exp = props(8)
         
         rknt=(/rkn,rkt/)  
         strength=(/strenthN,strenthT/)
         delta10=strenthN/rkn                         !! 法向初始损伤分离量
         delta20=strenthT/rkt                         !! 切向初始损伤分离量  
         delta1f=two*Gnc/strenthN                     !! 法向完全损伤分离量
         delta2f=two*Gtc/strenthT                     !! 切向完全损伤分离量
         
         gaussPara(1:2,1)= (/-constant,  -constant /)
         gaussPara(1:2,2)= (/ constant,  -constant /)
         gaussPara(1:2,3)= (/ constant,   constant /)
         gaussPara(1:2,4)= (/-constant,   constant /)
         return
       end subroutine      

    

c================================================VUEL子程序=========================================================
      subroutine vuel(
     *     nblock,
c          to be defined
     *     rhs,amass,dtimeStable,
     *     svars,nsvars,
     *     energy,
c          
     *     nnode,ndofel,
     *     props,nprops,
     *     jprops,njprops,
     *     coords,ncrd,
     *     u,du,v,a,
     *     jtype,jelem,
     *     time,period,dtimeCur,dtimePrev,kstep,kinc,lflags,
     *     dMassScaleFactor,
     *     predef,npredef,
     *     jdltyp,adlmag)
c   
      use math_mod
      include 'vaba_param.inc'
c     
      parameter ( zero = 0.d0, half = 0.5d0, one = 1.d0, two=2.d0,
     *           three=3.d0, quarder=0.25d0, eighth=0.125d0 )

c     operation code
      parameter ( jMassCalc            = 1,
     *            jIntForceAndDtStable = 2,
     *            jExternForce         = 3)

c     flags
      parameter (iProcedure = 1,
     *           iNlgeom    = 2,
     *           iOpCode    = 3,
     *           nFlags     = 3)

c     time
      parameter (iStepTime  = 1,
     *           iTotalTime = 2,
     *           nTime      = 2)

c     procedure flags
      parameter ( jDynExplicit = 17 )

c     energies 
      parameter ( iElPd = 1,
     *            iElCd = 2,
     *            iElIe = 3,
     *            iElTs = 4,
     *            iElDd = 5,
     *            iElBv = 6,
     *            iElDe = 7,
     *            iElHe = 8,
     *            iElKe = 9,
     *            iElTh = 10,
     *            iElDmd = 11,
     *            iElDc = 12,
     *            nElEnergy = 12)

c     predefined variables
      parameter ( iPredValueNew = 1,
     *            iPredValueOld = 2,
     *            nPred         = 2)    
    
      parameter ( nGP = 4 )
      parameter ( intervalU = nGP*3 )
      parameter ( TimeIncrementMax=2.0e-8 )
      parameter ( cos=0.9999984769132877d0 )
      parameter ( constant=one )  
      parameter ( rMin=1.0e-16 )
**-------------------------------------------------------------------
c
      dimension rhs(nblock,ndofel), amass(nblock,ndofel,ndofel),
     *     dtimeStable(nblock),
     *     svars(nblock,nsvars), energy(nblock,nElEnergy),
     *     props(nprops), jprops(njprops),
     *     jelem(nblock), time(nTime), lflags(nFlags),
     *     coords(nblock,nnode,ncrd), u(nblock,ndofel),
     *     du(nblock,ndofel), v(nblock,ndofel), a(nblock, ndofel),    
     *     dMassScaleFactor(nblock),
     *     predef(nblock, nnode, npredef, nPred), adlmag(nblock),
     *     disNode(3,nnode),uRelative(3,nGP),xyzMidsurf(3,nGP),
     *     gaussPara(2,nGP),tangs(3), tangt(3), unit(3),
     *     rShape(nGP), drds(nGP),drdt(nGP),uRelaGP(3), 
     *     shear(3),vNormal(3),rknt(2),
     *     Damage(2),strength(2),uRelaLocal(2),rJudge(3),tau(2),
     *     rNodeFc(3,nGP),TransCoords(3,3),rNodeLocal(2),rNodeGlobal(3),
     *     tauGlobal(3),uRelaGlobal(3),shearUnit(3),EleCoordVector(3,3),
     *     V1(3),V2(3),rNormaTauVect(3),rShearTauVect(3),rTauVect(3)
      
      if (jtype .eq. 2 .and.
     *    lflags(iProcedure).eq.jDynExplicit) then  
          
          call key_parameter(props,constant,rkn,rkt,strenthN,          !! 参数赋值
     *                      strengthT,strength,Gnc,Gtc,rho,exp,rknt,
     *                  delta10,delta20,delta1f,delta2f,gaussPara)

         if ( lflags(iOpCode).eq.jMassCalc ) then  
            call mass_matrix(nblock,rho,coords,ndofel,amass)           !! compute mass matrix of the cohesive element
            
         else if ( lflags(iOpCode) .eq.
     *             jIntForceAndDtStable) then  
             
           do kblock = 1, nblock
               call seperation_and_middle_suface(u,coords,disNode,nGP, !! 计算中性面点的坐标以及每对节点的分开量
     *            kblock,nblock,ndofel,ncrd,nnode,uRelative,xyzMidsurf)
                rNodeFc=zero
                iNumber = 1
                do iGp=1,nGP
                    DMICRT=svars(kblock,iNumber)
                    OldSDEG=svars(kblock,iNumber+1)
                    tau(1)=svars(kblock,iNumber+2)
                    tau(2) =svars(kblock,iNumber+3)

               call shape_q4(gaussPara(1,iGp), gaussPara(2,iGp), !! 求形函数及其偏导
     *                           rShape, drds, drdt)
               call normal_seperation_and_shear_seperation_and_unitVect
     *                   (drds,drdt,xyzMidsurf,rShape,uRelative,
     *                   vNormalMod,vNormal,shearUnit,uRelaLocal)      !! 求单位法向量和单位切向量，切法向分离量 l

               call Constitutive( DMICRT,delta10,delta20,uRelaLocal,Gnc,
     *            Gtc,exp,delta2f,delta1f,rkn,OldSDEG,SDEG,tau,strength)
                       
               if ( SDEG >= 1.0d0 ) then 
                   write(6,*) 'time(iTotalTime)',time(iTotalTime)
                   stop
               end if
                    rNormaTauVect=tau(1)*vNormal                       !! 法向应力向量
                    rShearTauVect=tau(2)*shearUnit                     !! 切向应力向量
                    rTauVect=rNormaTauVect+rShearTauVect               !! 总应力向量
                    !!  由积分点应力插值得到节点力 
                    do i=1,nGP
                        rNodeGlobal=vNormalMod*rShape(i)*rTauVect
                        rNodeFc(1:3,i)=rNodeGlobal+rNodeFc(1:3,i)
                    end do               
                    svars(kblock,iNumber)=DMICRT
                    svars(kblock,iNumber+1)=SDEG
                    svars(kblock,iNumber+2)=tau(1)
                    svars(kblock,iNumber+3)=tau(2)
                    iNumber = iNumber + 4
                end do
              !! 稳定时间增量
              dtimeStable(kblock) = TimeIncrementMax
              !! write(17,*) time(iTotalTime)
              !! 算出单元8个节点的节点力
              iBgn = 1
              jBgn = iBgn + intervalU 
              do i=1,nGP
                 rhs(kblock,iBgn:iBgn+2)=-rNodeFc(1:3,i)
                 rhs(kblock,jBgn:jBgn+2)=rNodeFc(1:3,i)
                 iBgn=iBgn+3
                 jBgn=jBgn+3
              end do  
      end do
      end if
      end if
      return
      end subroutine 