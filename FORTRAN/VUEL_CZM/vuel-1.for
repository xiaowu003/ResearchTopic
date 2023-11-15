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

      module parameter_mod
       use system_constant                                                                                                                                                                                                                                          
      implicit none
       integer, parameter :: nInteg1D = 2
       integer, parameter :: nInteg2D = 4                                                                                                                                                                                                                                                 
       integer, save :: nCoh                                                                                                                                                                                                                                                                                                                                                                                                                                                        
       integer, save :: nNdCoh = 8
       integer, save :: nHalfNd = 4
       integer, save :: nDim = 3
       integer, save :: nDofCoh = 24
       real(kind=RK), save :: m_thickness=1.0d0
       real(kind=RK), save :: m_rho
       real(kind=RK), save :: m_nStrength, m_tStrength
       real(kind=RK), save :: m_Gnc, m_Gtc
       real(kind=RK), save :: m_alpha,m_deltaNo, m_deltaSo
       real(kind=RK), save :: m_stiff
       real(kind=RK), save :: xiGaussPoint(2,30)
       real(kind=RK), save :: weight(30)
                                                                                                                                                                                                                                                                                                                                                                                                                         
      end module

      module math_mod
      use system_constant
      implicit none
      private 
      public :: vector_cross, vector_dot
      public :: vector_mod, vector_unit
      public :: shap_function_derivative
      public :: shape_q4
      contains

      subroutine vector_cross(vec1,vec2,cross) 
      implicit none
         !
         ! Function:
         !
         !  Calculate the cross product of two vectors.
         !
         ! Remarks:
         !
         !   RESULT = VEC1 �� VEC2
         !
         ! Parameters:
         !
         !  Vec1, input, an array of rank one with 3 elements.
         !  Vec2, input, an array of rank one with 3 elements.
         !  result, output, an array of rank one with 3 elements.
         !
         !! use system_constant

      real(kind=RK),intent(in):: vec1(3),vec2(3)
      real(kind=RK),intent(out):: cross(3)
      cross(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
      cross(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
      cross(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
      return
      end subroutine
       
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

          real(kind=RK),intent(in):: vec1(3), vec2(3)
          real(kind=RK),intent(out):: dot
          dot=vec1(1)*vec2(1) + vec1(2)*vec2(2)+vec1(3)*vec2(3)
          return
         end subroutine
        
      function vector_mod( vec ) result(rmod)
      !取模运算
      implicit none
      real(kind=RK) :: vec(3)
      real(kind=RK) :: rmod
      rmod = sqrt( vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3) )
      return
      end function 

      subroutine vector_unit( vec )
      !计算单位向量
      implicit none
      real(kind=RK),intent(inout):: vec(3)
      real(kind=RK), parameter :: TINY = 1.0e-25
      real(kind=RK) :: rmod
      rmod = sqrt( vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3) )
      if( rmod<TINY)then
        write(6,*)'Error in vector_unit:the rmod of the vector is 0.'
        stop
      end if
      vec = vec/rmod
      return
      end subroutine
 
      subroutine shap_function_derivative(gpos,dshap1,dshap2)
      !形函数的导数  
      implicit none
      real(kind=RK),intent(in)::gpos(2)
      real(kind=RK),intent(out)::dshap1(4),dshap2(4)
   
      dshap1(1)=-(1.0d0-gpos(2))
      dshap1(2)=1.0d0-gpos(2)
      dshap1(3)=1.0d0+gpos(2)
      dshap1(4)=-(1.0d0+gpos(2))
   
      dshap1=dshap1/4.0d0
   
      dshap2(1)=-(1.0d0-gpos(1))
      dshap2(2)=-(1.0d0+gpos(1))
      dshap2(3)=1.0d0+gpos(1)
      dshap2(4)=1.0d0-gpos(1)
   
      dshap2=dshap2/4.0d0
   
      return
      end subroutine
         
      subroutine shape_q4 ( r, s, t )
      !计算4节点参考4边形的形状函数
      implicit none
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
      real(kind=RK), intent(in) :: r, s
      real(kind=RK), intent(out) :: t(4)
      t(1) = 0.25d0 * ( 1.0d0 - r ) * ( 1.0d0 - s )
      t(2) = 0.25d0 * ( 1.0d0 + r ) * ( 1.0d0 - s )
      t(3) = 0.25d0 * ( 1.0d0 + r ) * ( 1.0d0 + s )
      t(4) = 0.25d0 * ( 1.0d0 - r ) * ( 1.0d0 + s )     
      return
      end subroutine shape_q4
      end module

      module cohesive_law_mod                                                                                                                                                                                                                                         
      use system_constant
      use math_mod                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
      implicit none                                                                                                                                                                                                                                                  
      private                                                                                                                                                                                                                                                        
      public :: compute_input_coh_displacement                                                                                                                                                                                                                
      public :: cohesive_stress_at_one_point
      public :: calculate_effective_displacement_point                                                                                                                                                                                 
      public :: calculate_tangent_vector
      character(len=MNL) :: msg                                                                                                                                                                                                                                      
      contains  

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
                                                                                                                                                                                                                                                                        
      deltai=0.0d0; deltaf=0.0d0                                                                                                                                                                                                                                      
      call calculate_mixed_mode_rate(tdelta,ndelta,beita)                                                                                                                                                                                         
      beita2=beita*beita                                                                                                                                                                                                                                          
      !beita2=0.7                                                                                                                                                                                                                                                 
      !for reference: Pinho S T, Iannucci L, Robinson P. Formulation and implementation 
      !of decohesion elements in an explicit finite element code[J].                                                                                                             
      !Composites Part A: Applied science and manufacturing, 2006, 37(5): 778-789.                                                                                                                                                                                
      !call obtain_transition_value_in_intrinsic_model(stiff,nStrengh,tStrength,delta_no,delta_so)                                                                                                                                                            
      if(ndelta<=0.0d0)then                                                                                                                                                                                                                                     
         deltai=deltaSo !effective displacement at transition point.                                                                                                                                                                                         
         if(deltai<=0.0d0) goto 222                                                                                                                                                                                                                            
         if(tStreng<=0.0d0) goto 111                                                                                                                                                                                                                           
         deltaf=2.0d0*Gs/tStreng !final displacement of the mixed mode.  !! Gs: G2c; Gn: G1c                                                                                                                                                                     
      else                                                                                                                                                                                                                                                    
         !effective displacement at transition point of the mixed triangle.                                                                                                                                                                                  
         deltai=deltaNo*deltaSo*sqrt( (1+beita2)/
     *          ( beita2*deltaNo*deltaNo+deltaSo*deltaSo) )                                                                                                                                                                    
         if(deltai<=0.0d0) goto 222                                                                                                                                                                                                                            
         if(stiff<=0.0d0) goto 333                                                                                                                                                                                                                             
         deltaf=2.0d0*(1+beita2)/(stiff*deltai)                                                                                                                                                                                                                  
         if(GN<=0.0d0.or.GS<=0.0d0) goto 444                                                                                                                                                                                                                     
            gna=(1.0d0/GN)**alfa; gta=(beita2/GS)**alfa                                                                                                                                                                                                             
            alfam=-1.0d0/alfa                                                                                                                                                                                                                                       
            deltaf=deltaf*(gna+gta)**alfam !final displacement of the mixed mode.                                                                                                                                                                               
      end if                                                                                                                                                                                                                                                  
      return                                                                                                                                                                                                                                                      
111   continue                                                                                                                                                                                                                                                    
        msg='Negative tstrength in compute_input_coh_displacement.'                                                                                                                                                                             
        write(6,*) msg
        stop                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
222   continue                                                                                                                                                                                                                                                    
        msg='Negative deltai in compute_input_coh_displacement.'                                                                                                                                                                              
        write(6,*) msg
        stop                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
333    continue                                                                                                                                                                                                                                                    
        msg='Negative stiffness in compute_input_coh_displacement.'                                                                                                                                                                        
        write(6,*) msg
        stop                                                                                                                                                                                                                                                      
444    continue                                                                                                                                                                                                                                                    
        msg='Negative GN or GS in compute_input_coh_displacement.'                                                                                                                                                                            
        write(6,*) msg
        stop                                                                                                                                                                                                                                                     
      end subroutine                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                        
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
      end subroutine                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                        
      subroutine cohesive_stress_at_one_point(stiff,deltai,deltaf,
     * delta,det_ndelta,det_tdelta,eff_delta,ghis,nvect,tvect,fail,sig)                                                                                                                                                                                
      implicit none                                                                                                                                                                                                                                                
      real(kind=RK),intent(in):: stiff                                                                                                                                                                                                                             
      real(kind=RK),intent(in)::delta(3),det_ndelta,det_tdelta,ghis
      real(kind=RK),intent(in)::eff_delta                                                                                                                                                                                      
      real(kind=RK),intent(in)::nvect(3),tvect(3),deltai,deltaf                                                                                                                                                                                                    
      !integer,intent(in)::co_eid, gith
      real(kind=RK),intent(out):: fail                                                                                                                                                                                                                            
      real(kind=RK),intent(out)::sig(3)                                                                                                                                                                                                                            
      real(kind=RK)::nsig(3),tsig(3)                                                                                                                                                                                                                               
      real(kind=RK)::dnsig,dtsig,det,temp                                                                                                                                                                                                                          
      real(kind=RK)::eta,modu,bvect(3),delta_tc                                                                                                                                                                                                                    
      real(kind=RK)::rnvect(3),rtvect(3)                                                                                                                                                                                                                           
      real(kind=RK),parameter::tol=1.0e-15                                                                                                                                                                                                                         
      real(kind=RK)::nsig0, tsig0, effsig0 !normal and tangential stress components at the transition point in a TSL.                                                                                                                                              
      real(kind=RK)::coef=0.0d0                                                                                                                                                                                                                                      
      
      fail = 0.0d0                                                                                                                                                                                                                                                                  
      sig=0.0d0; dnsig=0.0d0; dtsig=0.0d0; coef=0.0d0                                                                                                                                                                                                                      
      !if(interface_failure/=1) then                                                                                                                                                                                                                               
      if(ghis<=0.0d0) goto 111  !! !! ghis��Ч�ſ�������ʷ���ֵ                                                                                                                                                                                                 
      !end if                                                                                                                                                                                                                                                      
                                                                                                                                                                                                                                                                        
      !for reference: Wu L, Sket F, Molina-Aldareguia J M, et al. A study 
      ! of composite laminates failure using an anisotropic gradient-enhanced 
      !damage mean-field homogenization model[J]. Composite Structures, 2015, 126: 246-264.                                                                                                                                       
      !for reference: Hu N, Zemba Y, Okabe T, et al. A new cohesive model for simulating delamination propagation                                                                                                                                                  
      !in composite laminates under transverse loads[J]. Mechanics of Materials, 2008, 40(11): 920-935.                                                                                                                                                            
      !if(interface_failure==1)then   !intrinsic cohesive modeling                                                                                                                                                                                                 
      if(ghis<=deltai)then                                                                                                                                                                                                                                     
        if(det_ndelta>0.0d0) dnsig=stiff*det_ndelta                                                                                                                                                                                                            
        dtsig=stiff*det_tdelta                                                                                                                                                                                                                               
      else if(ghis>deltai.and.ghis<=deltaf)then
        fail = 1.01d0                                                                                                                                                                                                                
        coef=deltaf/ghis; coef=coef*(ghis-deltai)/(deltaf-deltai)                                                                                                                                                                                            
        if(det_ndelta>0.0d0) dnsig=(1-coef)*stiff*det_ndelta                                                                                                                                                                                                   
        dtsig=(1.0d0-coef)*stiff*det_tdelta                                                                                                                                                                                                                    
      end if                                                                                                                                                                                                                                                   
      !else                                                                                                                                                                                                                                                        
      !    goto 555                                                                                                                                                                                                                                                
      !end if                                                                                                                                                                                                                                                      
      modu=sqrt(delta(1)*delta(1)+delta(2)*delta(2)+delta(3)*delta(3))                                                                                                                                                                                            
      if(modu<=0.0d0) goto 222                                                                                                                                                                                                                                      
         bvect=delta                                                                                                                                                                                                                                                 
         bvect=bvect/modu                                                                                                                                                                                                                                            
         rnvect=nvect; rtvect=tvect                                                                                                                                                                                                                                  
         if(det_ndelta>0.0d0) then                                                                                                                                                                                                                                     
            call vector_dot(nvect,bvect,det)                                                                                                                                                                                                                             
         end if                                                                                                                                                                                                                                                      
         call vector_dot(tvect,bvect,det)                                                                                                                                                                                                                                 
         if(det_ndelta>0.0d0) then                                                                                                                                                                                                                                     
            nsig=dnsig*rnvect                                                                                                                                                                                                                                       
         else                                                                                                                                                                                                                                                        
            nsig=0.0d0                                                                                                                                                                                                                                                
         end if                                                                                                                                                                                                                                                      
         tsig=dtsig*rtvect                                                                                                                                                                                                                                           
         sig=nsig+tsig                                                                                                                                                                                                                                               
      return                                                                                                                                                                                                                                                      
111   continue                                                                                                                                                                                                                                                    
      msg='Negative ghis in cohesive_stress_at_one_point.'                                                                                                                                                                                         
      write(6,*) msg
      stop                                                                                                                                                                                                                                                       
222   continue                                                                                                                                                                                                                                                    
      msg='Zero mod in cohesive_stress_at_one_point.'                                                                                                                                                                                          
      write(6,*) msg
      stop                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
      end subroutine                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                        
      subroutine calculate_tangent_vector(tdelta,tvect)
      !计算切向量                                                                                                                                                                                                              
      implicit none                                                                                                                                                                                                                                               
      real(kind=RK),intent(in)::tdelta(3)                                                                                                                                                                                                                         
      real(kind=RK),intent(out)::tvect(3)                                                                                                                                                                                                                         
      real(kind=RK)::det                                                                                                                                                                                                                                          
      real(kind=RK),parameter::tol=1.0e-16                                                                                                                                                                                                                        
       det=sqrt(tdelta(1)*tdelta(1)+tdelta(2)*tdelta(2)+
     *     tdelta(3)*tdelta(3))                                                                                                                                                                                       
       if(det<tol)then                                                                                                                                                                                                                                             
         tvect=(/0.0d0 ,0.0d0, 0.0d0/)                                                                                                                                                                                                                                         
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
      real(kind=RK),save::no=0.0d0,so=0.0d0                                                                                                                                                                                                                           
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
111   continue                                                                                                                                                                                                                                                    
      msg='Nagetive stiffness in obtain_transition_value.'                                                                                                                                                                   
      write(6,*) msg
      stop                                                                                                                                                                                                                                                      
      end subroutine                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                        
      subroutine calculate_effective_displacement_point(det_ndelta,
     *  det_tdelta,eff_delta)                                                                                                                                                                              
      implicit none                                                                                                                                                                                                                                               
      real(kind=RK),INTENT(in)::det_ndelta,det_tdelta                                                                                                                                                                                                             
      real(kind=RK),intent(out)::eff_delta                                                                                                                                                                                                                        
      real(kind=RK)::ndel,tdel                                                                                                                                                                                                                                    
       ndel=(det_ndelta+abs(det_ndelta))*0.5d0                                                                                                                                                                                                                       
       tdel=det_tdelta                                                                                                                                                                                                                                             
       eff_delta=sqrt(ndel*ndel+tdel*tdel)                                                                                                                                                                                                                         
       return                                                                                                                                                                                                                                                      
      end subroutine                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                        
      subroutine calculate_mixed_mode_rate(tdelta,ndelta,beita)                                                                                                                                                                                       
      implicit none                                                                                                                                                                                                                                               
      real(kind=RK),intent(in)::tdelta,ndelta                                                                                                                                                                                                                     
      real(kind=RK),intent(out)::beita                                                                                                                                                                                                                            
       beita=0.0d0                                                                                                                                                                                                                                                   
       if(ndelta<=0.0d0)then                                                                                                                                                                                                                                         
         beita=1.0e20                                                                                                                                                                                                                                            
       else                                                                                                                                                                                                                                                        
         beita=tdelta/ndelta                                                                                                                                                                                                                                     
       end if                                                                                                                                                                                                                                                      
       return                                                                                                                                                                                                                                                      
      end subroutine                                                                                                                                                                                                                                                  
      end module                                                                   

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

      subroutine set_parameter(nBlock,nNode,nCrd,nDofEl,nProps,props)
       implicit none
       integer, intent(in) :: nBlock,nNode,nCrd,nDofEl,nProps
       real(kind=RK), intent(in) :: props(nProps)
       real(kind=RK), parameter :: constant = 1.0d0 !0.577350269189626d0
       nCoh = nBlock
       nNdCoh = nNode
       !nInteg2D = nNode/2
       nHalfNd = nNdCoh/2
       nDim = nCrd
       nDofCoh = nDofEl                                                                                                                                                                                                                                                                                                                                                                                                                                                          
       m_stiff = props(1)
       m_nStrength = props(3)
       m_tStrength = props(4)
       m_Gnc = props(5) !! m_GIC
       m_Gtc = props(6)
       m_rho = props(7)
       m_alpha = props(8)
       m_deltaNo = m_nStrength/m_stiff                                                                                                                                                                                                                            
       m_deltaSo = m_tStrength/m_stiff
       xiGaussPoint(1:2,1) = (/-constant,  -constant /)
       xiGaussPoint(1:2,2) = (/ constant,  -constant /)
       xiGaussPoint(1:2,3) = (/ constant,   constant /)
       xiGaussPoint(1:2,4) = (/-constant,   constant /)
       weight(1:4) = 1.0d0      
       return
      end subroutine      

      subroutine mass_matrix(nBlock,nNode,nCrd,nDofEl,coords,amass) 
      implicit none
      integer, intent(in) :: nBlock,nNode,nCrd,nDofEl
      real(kind=RK), intent(in) :: coords(nBlock,nNode,nCrd)
      real(kind=RK), intent(out) :: amass(nBlock,nDofEl,nDofEl)
      integer :: i, j
      real(kind=RK) :: area, am0, S1, S2, elemMass
      real(kind=RK) :: edge(3,3), V1(3), V2(3)
 
      do i = 1, nCoh
         edge(1,:)=coords(i,2,:)-coords(i,1,:)
         edge(2,:)=coords(i,3,:)-coords(i,1,:)
         edge(3,:)=coords(i,4,:)-coords(i,1,:)
         call vector_cross( edge(1,:),edge(2,:),V1 )  
         S1= vector_mod( V1 )
         call vector_cross(edge(2,:),edge(3,:),V2)
         S2= vector_mod( V2 )
         area=0.5d0*(S1+S2)
         elemMass=m_thickness*area*m_rho
         am0=elemMass/real(nNdCoh,kind=RK)
          do j=1,nDofCoh
            amass(i,j,j) = am0
          end do
         end do
      end subroutine  

      subroutine one_element_cohesive_force(iEleM,nVar,u,coords,
     * var,fcNd) 
      implicit none
      integer, intent(in) :: iElem, nVar
      real(kind=RK),intent(in) :: u(nDofCoh)
      real(kind=RK),intent(in) :: coords(nNdCoh,nDim)
      real(kind=RK),intent(inout) :: var(nvar)
      real(kind=RK),intent(out) :: fcNd(nDim,nNdCoh)  !! may be not allowed
  
      real(kind=RK)::normal(3),det
      integer::fCount
      integer::i,j, iEffGap, ibgn
      real(kind=RK) :: fail, failElem
      real(kind=RK) :: maxEffGap
      real(kind=RK) :: xi(2), gapPt(3)
      real(kind=RK) ::  xyzMid(3,4), gapNd(3,4)
      real(kind=RK):: fcAdd(3,4),fcPt(3)


      if( nVar<nInteg2D*2+1 ) goto 111
      failElem = var(1)
      if( failElem>=2 ) return !! element complete failure
 
      call seperation_and_middle_suface(u,coords,gapNd,xyzMid)
      fcNd=0.0d0; fcount=0
      ibgn = 2
      do i=1,nInteg2D
         fail = var(ibgn)
         if( fail>=2 )then
           fcount=fcount+1; cycle  !means this gauss point is in failure. 
         end if
         xi = xiGaussPoint(:,i)
         call get_unit_normal_vector(xi,xyzMid,det,normal)
         call calculate_gap_at_integration_point(xi,gapNd,gapPt)
         !call obtain_failure_values_of_newton_integration_points(
     *   !     co_eid,ith,fail=fail)
        iEffGap = ibgn+1 
        maxEffGap = var(iEffGap)             
        call one_point_cohesive_force(det,fCount,maxEffGap,gapPt,
     *       normal,fail,fcPt)
        var(iEffGap) = maxEffGap
        var(ibgn) = fail
        call distribute_force_to_nodes(xi,fcPt,fcAdd)
        fcNd(:,1:nHalfNd)=fcNd(:,1:nHalfNd)+fcAdd
        ibgn = ibgn + 2
      end do
      fcNd(:,nHalfNd+1:nNdCoh) = fcNd(:,1:nHalfNd)
      fcNd(:,1:nHalfNd)  = -fcNd(:,1:nHalfNd)
      if( fCount==nInteg2D ) var(1) = 2.01d0
      return
111   continue
      msg='Wrong length of svars in one_element_cohesive_force.'
      write(6,*) msg
      stop 
      return
      end subroutine

      subroutine distribute_force_to_nodes(xi,fcPt,fcNd)    
      implicit none
      real(kind=RK),intent(in)::xi(2), fcPt(3)
      real(kind=RK),intent(out)::fcNd(3,4)
      real(kind=RK)::shap(4)
      integer::nNum, i
   
      fcNd=0.0d0
      call shape_q4(xi(1),xi(2),shap)
      nNum = nInteg2D
      do i=1,nNum
        fcNd(:,i)=shap(i)*fcPt
      end do
             
      return
      end subroutine

      subroutine calculate_gap_at_integration_point(xi,gapNd,gapPoint)
      implicit none
      real(kind=RK), intent(in) :: xi(2) 
      real(kind=RK),intent(in) :: gapNd(3,4)
      real(kind=RK),intent(out) :: gapPoint(3)
      real(kind=RK) :: t(4)
      integer :: i
      call shape_q4(xi(1),xi(2),t)
      gapPoint = 0.0d0
      do i = 1, nHalfNd
         gapPoint = gapPoint+t(i)*gapNd(:,i)
      end do
      return
      end subroutine

      subroutine get_unit_normal_vector(xi,xyz,det,normal)
      implicit none
      real(kind=RK),intent(in)::xi(2)  
      real(kind=RK),intent(in)::xyz(3,4)
      real(kind=RK),intent(out)::normal(3),det
      real(kind=RK)::tvect1(3),tvect2(3)
      real(kind=RK)::tmod
      real(kind=RK),parameter::tol=1.0e-10
      call calculate_tangential_vector(xyz,xi,tVect1,tVect2)
      call vector_cross(tVect1,tVect2,normal)
      tmod=vector_mod(normal)
      if(abs(tmod)<tol) goto 333
      det=tMod
      normal=normal/tMod    
      return
333   continue
      msg='Negative tmod in get_unit_normal_vector'
      write(6,*) msg
      stop 
      end subroutine

      subroutine calculate_tangential_vector(xyz,gpos,tvect1,tvect2)          
      implicit none
      real(kind=RK),intent(in)::xyz(3,4),gpos(2)
      real(kind=RK),intent(out)::tvect1(3),tvect2(3)
      real(kind=RK)::dshap1(4),dshap2(4)
      integer::id
    
      tvect1=0.0d0; tvect2=0.0d0
      call shap_function_derivative(gpos,dshap1,dshap2)

      do id=1,4
        tvect1=tvect1+xyz(:,id)*dshap1(id)
        tvect2=tvect2+xyz(:,id)*dshap2(id)
      end do

      return
      end subroutine

      subroutine seperation_and_middle_suface(u,coords, gap, 
     *  xyzMidsurf)     
      implicit none
      real(kind=RK), intent(in) :: u(nDofCoh), coords(nNdCoh,nDim)
      real(kind=RK), intent(out) :: gap(nDim,nInteg2D)
      real(kind=RK), intent(out) :: xyzMidsurf(:,:) 
      real(kind=RK) :: disNode(nDim,nNdCoh)
      integer :: k
       disNode=reshape( u, (/nDim,nCoh/) )
       !iBgn = 1
       !do k=1,nnode
       !   disNode(:,k)=u(iBgn:iBgn+2)
       !   iBgn = iBgn + 3
       !end do
                
       do k = 1, nHalfNd
          xyzMidsurf(:,k) = ( coords(k,:)+coords(k+nInteg2D,:)+
     *                      disNode(:,k)+disNode(:,k+nInteg2D))/2.0d0   
          gap(:,k)=disNode(:,k+nInteg2D)-disNode(:,k)  !! ÿ�Խڵ�ķ�����
       end do
       return
      end subroutine

      subroutine one_point_cohesive_force(area,fCount,gHis,gapVect,
     * normal,fail,fc)  !! gpos, Hpos, fcount                                                                                                                                                                
      implicit none
      real(kind=RK),intent(in):: area         !! ��˹���ֵ��Ӧ���ſɱ�����ʽֵ��Ҳ�������                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
      real(kind=RK),intent(in):: gapVect(3)   !! gapVect�Ǹø�˹���ֵ��ϵķ���������                                                                                                                                                                              
      real(kind=RK),intent(in):: normal(3)    !! normal���м����ϸø�˹���ֵ��ϵĵ�λ������
      integer, intent(inout) :: fCount 
      real(kind=RK),intent(inout) :: gHis      !! gHis��Ч�ſ�������ʷ���ֵ                                                                                                                                                                 
      !real(kind=RK),intent(in):: gpos(2)     !! gpos ��˹���ֵ����Ȼ����                                                                                                                                                                                                                                                                                                                                                              
      !real(kind=RK),intent(in):: Hpos(2)     !! Hpos�Ǹ�˹���ֵ��Ȩϵ��
      real(kind=RK),intent(out):: fail                                                                                                                                                                                
      real(kind=RK),intent(out)::fc(3)      !!  fc �û��ֵ������ȫ������ϵ�µ�cohesive force                                                                                                                                                                     
      !integer,intent(inout)::fcount                                                                                                                                                                                                                              
      real(kind=RK)::tmod,nGap,tGap                                                                                                                                                                                                                          
      real(kind=RK)::co_sig(3),nGapVect(3),tGapVect(3),effDelta                                                                                                                                                                                                   
      real(kind=RK)::tfc(3),tvect(3)                                                                                                                                                                                                                              
      real(kind=RK)::deltai,deltaf, gVect(3)                                                                                                                                                                                                              
      fc=0.0d0                                                                                                                                                                                                                                                      
      gVect = gapVect                                                                                                                                                                                                                                             
      nGap=gVect(1)*normal(1)+gVect(2)*normal(2)+gVect(3)*normal(3) !! �������ֵ                                                                                                                                                                               
      nGapVect = nGap*normal          !! ����������� normal separation vector                                                                                                                                                                                    
      tGapVect = gVect - nGapVect     !! �����������                                                                                                                                                                                                             
      tGap = sqrt( tGapVect(1)*tGapVect(1)+tGapVect(2)*tGapVect(2)+
     *       tGapVect(3)*tGapVect(3) ) !! �������ֵ                                                                                                                                                        
      if(nGap<0.0)then !��ѹ״̬                                                                                                                                                                                                                                  
         tfc=-abs(nGap)*m_stiff*area*normal                                                                                                                                                                                                                      
         fc=fc+tfc                                                                                                                                                                                                                                               
      end if                                                                                                                                                                                                                                                      
      call calculate_effective_displacement_point(nGap,tGap,effDelta)  !! ������Ч�ſ���                                                                                                                                                                          
      if(effDelta<=0.0d0) return                                                                                                                                                                                                                                    
      call compute_input_coh_displacement(m_stiff, 
     *     m_nStrength,m_tStrength,m_Gnc,m_Gtc,m_alpha,m_deltaNo,
     *     m_deltaSo,tGap,nGap,deltai,deltaf) !! �����������ε�����λ��                                                                                                                                                     
      if(effDelta>=deltaf)then !failure of the gauss point.                                                                                                                                                                                                         
        fCount=fCount+1
        ghis = effDelta                                                                                                                                                                                                                                        
        !m_BindDem(iDem)%status(iFacet) = 2
        fc = 0.0d0
        fail = 2.01d0                                                                                                                                                                                                                      
      else                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
        !m_BindDem(iDem)%status(iFacet) = 1                                                                                                                                                                                                                     
        if(ghis<=effDelta)then                                                                                                                                                                                                                                  
          ghis = effDelta                                                                                                                                                                                                                                     
          !m_bindDem(iDem)%maxEffGap(iFacet) = ghis                                                                                                                                                                                                                                                                                                                                                                                  
        end if                                                                                                                                                                                                                                                  
        call calculate_tangent_vector(tGapVect,tvect)  !! ��������λ����                                                                                                                                                                                      
        call cohesive_stress_at_one_point(m_stiff,deltai,deltaf,gVect,
     *       nGap,tGap,effDelta,ghis,normal,tvect,fail,co_sig)                                                                                                                                                                      
        tfc = co_sig*area                                                                                                                                                                                                                                       
        fc = fc + tfc                                                                                                                                                                                                                                         
      end if                                                                                                                                                                                                                                                      
      return                                                                                                                                                                                                                                                      
111   continue                                                                                                                                                                                                                                                    
        msg='mistake of the value of effDelta. in one_point_bind_force'                                                                                                                                                                                             
        write(6,*) msg
        stop                                                                                                                                                                                                                                                      
222   continue                                                                                                                                                                                                                                                    
        msg='mistake of the value of ghis. in one_point_bind_force'                                                                                                                                                                                                 
        write(6,*) msg
        stop                                                                                                                                                                                                                                                       
      end subroutine
      end module

c================================================VUEL�ӳ���=========================================================
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
      use system_constant
      use math_mod
      use parameter_mod
      use cohesive_element_mod
      implicit none
c     operation code
      integer, parameter :: jMassCalc = 1
      integer, parameter :: jIntForceAndDtStable = 2
      integer, parameter :: jExternForce         = 3

c     flags
      integer, parameter :: iProcedure = 1
      integer, parameter ::  iNlgeom    = 2
      integer, parameter :: iOpCode    = 3
      integer, parameter :: nFlags     = 3

c     time
      integer, parameter :: iStepTime  = 1
      integer, parameter :: iTotalTime = 2
      integer, parameter :: nTime      = 2

c     procedure flags
      integer, parameter ::  jDynExplicit = 17

c     predefined variables
      integer, parameter :: iPredValueNew = 1
      integer, parameter :: iPredValueOld = 2
      integer, parameter :: nPred         = 2 

      integer, parameter :: nElEnergy = 12

      integer, intent(in) :: jType
      integer, intent(in) :: jdlTyp   !! integer identifying the load type
      integer, intent(in) :: kStep    !! Current step number
      integer, intent(in) :: kinc     !! Current increment number
      integer, intent(in) :: nPreDef  !! field variable number ��not used now��
      integer, intent(in) :: nBlock   !! total solid element number
      integer, intent(in) :: nCrd     !! n deminsion
      integer, intent(in) :: nDofEl   !! total dof number of one element
      integer, intent(in) :: njprops  !! user defined properties (integer)
      integer, intent(in) :: nNode    !! total node number of one element  (real)
      integer, intent(in) :: nProps   !! user defined properties 
      integer, intent(in) :: nsvars   !! state variable number
      integer, intent(in) :: jElem(nBlock)    !! the element ID
      integer, intent(in) :: jProps(njProps)  !! User-defined number of integer property values  
      integer, intent(in) :: lFlags(nFlags)
      integer, intent(in) :: dTimeStable(nFlags)
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
      real(kind=RK), intent(out) :: rhs(nBlock,nDofEl) ! the contributions of each element to the 
                          ! right-hand-side vector of the overall system of equations (e.g.internal force)
      real(kind=RK), intent(out) :: aMass(nBlock,nDofEl,nDofEl)  !!  mass matrix of each element
!      real(kind=RK), intent(out) :: dTimeStable(nBlock) !! stable time increment of each element
      real(kind=RK), intent(inout) :: sVars(nBlock,nsvars) !! the solution-dependent state variables
      real(kind=RK), intent(in) ::  energy(nBlock,nElEnergy)  !! not computed

      integer :: i

      if( jtype==2 .and. lflags(iProcedure)==jDynExplicit )then 
 
         call set_parameter(nBlock,nNode,nCrd,nDofEl,nProps,props)

         if( lflags(iOpCode)==jMassCalc )then  
            call mass_matrix(nBlock,nNode,nCrd,nDofEl,coords,amass) 
         else if ( lflags(iOpCode)==jIntForceAndDtStable)then  
            do i = 1, nblock
                call one_element_cohesive_force( i,nsvars,u(i,:),
     *          coords(i,:,:),svars(i,:),rhs(i,:) )
                !write(6,*) 'node force'
                !write(6,*)rhs(i,1:12)
                !write(6,*)rhs(i,13:24)
                if( svars(1,2)>0.0 )then
                   write(6,*) 'time', time(1)
                   write(6,*) 'svars', svars(1,:)
                   write(6,*) 'fc', rhs(1,1:12)
                   !stop
                end if
            end do
         end if
      end if
      return
      end subroutine                                                                                                                        