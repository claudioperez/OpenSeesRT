***************************************************************************************************
! Software Name: PeriFEM Software for Structural Failure Simulation Based on Quadrilateral Mesh for ABAQUS [PeriFEM-ABAQUS-Qua2D]
!
! PeriFEM-ABAQUS-Qua2D is a PeriFEM software for structural failure simulation
! based on quadrilateral mesh for ABAQUS, released by the Near-Field Dynamics
! and Multiscale Methods team led by Professor Han Fei of Dalian University of
! Technology on January 1, 2022. It was reviewed by the China Copyright
! Protection Center and registered with the National Copyright Administration of
! the People's Republic of China on August 25, 2022. 
!
! Dalian University of Technology owns all rights to this software.
!
***************************************************************************************************     
!
! Function: This program is the basic program for compiling ABAQUS, assembling
! the element stiffness matrix by determining whether the near-field dynamics
! bond is broken.
!      
***************************************************************************************************     
! Developers: Han Fei (hanfei@dlut.edu.cn)
!             Zhang Jianyu
!             Yao Chen
!     
***************************************************************************************************
!   
! Version:   ABAQUS 2019
!            Visual Studio 2015
!            Intel Parallel Studio XE 2017    
!
***************************************************************************************************
*                                                                                                 *
*                                        Main Program Section                                     *     
*                                                                                                 *
***************************************************************************************************
      module kvisual
      Implicit None
      real*8 ele_bond_sum_state_fem(4,1,100000)    ! FEM element damage, calculation limit is 10e5 elements
      integer ele_bond_state_pd(1,1,100,100000)    ! Used to determine whether it is a related element, calculation limit is 10e2 related elements, total number of elements is 10e5
      integer ele_bond_sum_state_pd(4,1,100,100000)! Number of broken bonds in PD elements
      save
      end module
      
      ! Call UEL
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
  
      Use kvisual
      INCLUDE 'ABA_PARAM.INC'

! ******************** State Variables ***************************************
! svars(1) = Number of new broken bonds per element
! svars(2) = Total bond damage per element
! svars(3) = subroutine entry count
! svars(4-19) = Bond status 0/1
! svars(20) = Total number of unbroken bonds per element
! svars(21) = Total number of broken bonds per element
! **********************************************************************
      
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(*),ENERGY(7),PROPS(*),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(*),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(4),
     5     JPROPS(*)
      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

                
      INTEGER          loop_number,new_bre_bond,i,j
     &                 para_build_bond
      DOUBLE PRECISION x_Gspoint(4),y_Gspoint(4),block_stiff(16,16),  
     &                 curr_gauss(2,8),ref_gauss(2,8),
     &                 ref_bond,curr_bond,curr_coord(2,8),
     &                 rhsmatrx(NDOFEL,1)             
      CHARACTER*4      label(3)
      parameter(ntens=4)

      ! Define Gauss points
      x_Gspoint(1)=-0.577350269189626d0
      x_Gspoint(2)=0.577350269189626d0
      x_Gspoint(3)=0.577350269189626d0
      x_Gspoint(4)=-0.577350269189626d0
    
      y_Gspoint(1)=-0.577350269189626d0
      y_Gspoint(2)=-0.577350269189626d0
      y_Gspoint(3)=0.577350269189626d0
      y_Gspoint(4)=0.577350269189626d0
      
      ! Initialize global variables
      if (kinc.eq.1) then
          ele_bond_state_pd(1,1,Int(COORDS(2,9)),Int(COORDS(1,9)))=0
          do i=1,NSVARS
              svars(i)=0.0d0
          enddo
          do i=1,4
            ele_bond_sum_state_fem(i,1,int(COORDS(1,9))) = 0.0d0
            ele_bond_sum_state_pd(i,1,int(COORDS(2,9)),int(COORDS(1,9)))
     &      =0
          enddo
      endif
      
      ! Initialize matrices  
      call mat_zero(curr_gauss,2,8)      
      call mat_zero(ref_gauss,2,8) 
      call mat_zero(curr_coord,2,8) 
      call mat_zero(AMATRX,NDOFEL,NDOFEL) 
      call mat_zero(block_stiff,16,16)
      call mat_zero(rhsmatrx,NDOFEL,NRHS)
      
      ! Calculate current node coordinates
      do i=1,8  
            curr_coord(1,i)=COORDS(1,i)+U(i*2-1)
            curr_coord(2,i)=COORDS(2,i)+U(i*2)
      enddo
      
      ! Calculate Gauss coordinates  
      Call gauss_loc(curr_coord,curr_gauss)
      Call gauss_loc(COORDS,ref_gauss)
      
      ! Determine bond breakage                
      loop_number=1           ! Initial loop variable count 
      new_bre_bond=0          ! Current number of broken bonds
      do i=1,4
          do j=1,4
      ref_bond=Sqrt((ref_gauss(1,i)-ref_gauss(1,j+4))*    ! Initial bond length
     &               (ref_gauss(1,i)-ref_gauss(1,j+4))+
     &               (ref_gauss(2,i)-ref_gauss(2,j+4))*
     &                (ref_gauss(2,i)-ref_gauss(2,j+4))) 
      curr_bond=Sqrt((curr_gauss(1,i)-curr_gauss(1,j+4))*  ! Current bond length
     &              (curr_gauss(1,i)-curr_gauss(1,j+4))+
     &              (curr_gauss(2,i)-curr_gauss(2,j+4))*
     &              (curr_gauss(2,i)-curr_gauss(2,j+4)))
          If(svars(3+loop_number).eq.0.0d0) then
              If(ref_bond.ne.0.0d0) then  
                  If (((curr_bond-ref_bond)/ref_bond).ge.PROPS(2)) then ! Determine if the current bond is broken based on damage value
                      new_bre_bond = new_bre_bond+1 
                      svars(3+loop_number) = 1.0d0  
                  endif
              endif            
          endif
          loop_number=loop_number+1
          enddo
      enddo   
      ! Assemble element stiffness matrix
      loop_number=1           ! Initial loop variable count
      para_build_bond=0       ! Number of unbroken bonds in the element
      do i=1,4
          do j=1,4
              If(svars(3+loop_number).ne.1.0d0) then
                  call calc_gauss_intg(x_Gspoint(i),y_Gspoint(i),COORDS,
     &                        x_Gspoint(j),y_Gspoint(j),block_stiff,
     &                        PROPS(1),PROPS(3),PROPS(4))
                   AMATRX=AMATRX+block_stiff     ! Element stiffness matrix
                   para_build_bond=para_build_bond+1
              endif
            loop_number=loop_number+1
          enddo
      enddo 
      
      call mat_mult(AMATRX,U,rhsmatrx,NDOFEL,NDOFEL,NDOFEL,1)
      RHS(1:NDOFEL,1) = -rhsmatrx(1:NDOFEL,1)

      ! Update state variables
      svars(1)=new_bre_bond            ! Current number of broken bonds          
      svars(20)=para_build_bond        ! Number of unbroken bonds per element
      svars(21)=16-para_build_bond     ! Total number of broken bonds per element
      svars(2)=(16-para_build_bond)/16 ! Total bond damage per element
      ! Store the number of broken bonds
      ele_bond_state_pd(1,1,Int(COORDS(2,9)),Int(COORDS(1,9)))=1
      ele_bond_sum_state_pd(1,1,Int(COORDS(2,9)),Int(COORDS(1,9)))=svars(21)
      ele_bond_sum_state_pd(2,1,Int(COORDS(2,9)),Int(COORDS(1,9)))=svars(21)
      ele_bond_sum_state_pd(3,1,Int(COORDS(2,9)),Int(COORDS(1,9)))=svars(21)
      ele_bond_sum_state_pd(4,1,Int(COORDS(2,9)),Int(COORDS(1,9)))=svars(21)
 
      RETURN
      END

***************************************************************************************************
*                                                                                                 *
*                                           subroutine                                            *    
*                                                                                                 *
***************************************************************************************************             
!------------------------------------------------------------------------------
!                            Calculate Jacobi Matrix  
!------------------------------------------------------------------------------ 

        subroutine calc_Jdet(xRef,i_xi,i_eta,j_xi,
     &                       j_eta,i_Jdet,j_Jdet)
        Double Precision xref(2,8),i_xi,i_eta,j_xi,j_eta,
     &                   i_J(2,2),i_N_diff(2,4),j_J(2,2),
     &                   i_Jdet,j_Jdet,j_N_diff(2,4)
        Integer i
        
        Call mat_zero(i_J,2,2)
        Call mat_zero(j_J,2,2)

        ! Partial derivatives of shape functions with respect to parametric coordinates in the main element
        i_N_diff(1,1) = (i_eta - 1.0d0) / 4.0d0
        i_N_diff(2,1) = (i_xi - 1.0d0) / 4.0d0
        i_N_diff(1,2) = (1.0d0 - i_eta) / 4.0d0
        i_N_diff(2,2) = -1.0d0 * (i_xi + 1.0d0) / 4.0d0
        i_N_diff(1,3) = (1.0d0 + i_eta) / 4.0d0
        i_N_diff(2,3) = (1.0d0 + i_xi) / 4.0d0
        i_N_diff(1,4) = -1.0d0 * (1.0d0 + i_eta) / 4.0d0
        i_N_diff(2,4) = (1.0d0 - i_xi) / 4.0d0
        
        ! Partial derivatives of shape functions with respect to parametric coordinates in the auxiliary element
        j_N_diff(1,1) = (j_eta - 1.0d0) / 4.0d0
        j_N_diff(2,1) = (j_xi - 1.0d0) / 4.0d0
        j_N_diff(1,2) = (1.0d0 - j_eta) / 4.0d0
        j_N_diff(2,2) = -1.0d0 * (j_xi + 1.0d0) / 4.0d0
        j_N_diff(1,3) = (1.0d0 + j_eta) / 4.0d0
        j_N_diff(2,3) = (1.0d0 + j_xi) / 4.0d0
        j_N_diff(1,4) = -1.0d0 * (1.0d0 + j_eta) / 4.0d0
        j_N_diff(2,4) = (1.0d0 - j_xi) / 4.0d0
        
        ! Calculate Jacobi matrix for the main element
        do i=1,4
            i_J(1,1)=i_J(1,1)+i_N_diff(1,i)*xref(1,i)
            i_J(1,2)=i_J(1,2)+i_N_diff(1,i)*xref(2,i)
            i_J(2,1)=i_J(2,1)+i_N_diff(2,i)*xref(1,i)
            i_J(2,2)=i_J(2,2)+i_N_diff(2,i)*xref(2,i)
        enddo
        
        ! Calculate Jacobi matrix for the auxiliary element
        do i=1,4
            j_J(1,1)=j_J(1,1)+j_N_diff(1,i)*xref(1,i+4)
            j_J(1,2)=j_J(1,2)+j_N_diff(1,i)*xref(2,i+4)
            j_J(2,1)=j_J(2,1)+j_N_diff(2,i)*xref(1,i+4)
            j_J(2,2)=j_J(2,2)+j_N_diff(2,i)*xref(2,i+4)
        enddo
        
        ! Calculate the determinant of the Jacobi matrix
        i_Jdet=i_J(1,1)*i_J(2,2)-i_J(1,2)*i_J(2,1)
        j_Jdet=j_J(1,1)*j_J(2,2)-j_J(1,2)*j_J(2,1)
        
        end subroutine calc_Jdet

!-----------------------------------------------------------------------------
!                       Calculate Gauss Integration Matrix
!-----------------------------------------------------------------------------    
       subroutine calc_gauss_intg(i_xi,i_eta,x_ref,j_xi,
     &                    j_eta,block_stiff,ele_esize,param_1,param_2)
        implicit none
    
        double precision peri_B(2,16),i_xi,i_eta,x_ref(2,8),i_N(4),j_xi,
     &                j_eta,j_N(4),peri_D(2,2),i_realcoord(2),
     &                j_realcoord(2),peri_DB(2,16),peri_Jdet,i_Jdet,
     &                j_Jdet,tran_peri_B(16,2),block_stiff(16,16),
     &                kexi,L,ele_esize,param_1,param_2
    
        integer i
        
        
        call calc_Jdet(x_Ref,i_xi,i_eta,j_xi,
     &        i_eta,i_Jdet,j_Jdet)

        ! Shape functions
        i_N(1)=1.0d0/4.0d0*(1.0d0-i_xi)*(1.0d0-i_eta)
        i_N(2)=1.0d0/4.0d0*(1.0d0+i_xi)*(1.0d0-i_eta)
        i_N(3)=1.0d0/4.0d0*(1.0d0+i_xi)*(1.0d0+i_eta)
        i_N(4)=1.0d0/4.0d0*(1.0d0-i_xi)*(1.0d0+i_eta)
    
        j_N(1)=1.0d0/4.0d0*(1.0d0-j_xi)*(1.0d0-j_eta)*-1.0d0
        j_N(2)=1.0d0/4.0d0*(1.0d0+j_xi)*(1.0d0-j_eta)*-1.0d0
        j_N(3)=1.0d0/4.0d0*(1.0d0+j_xi)*(1.0d0+j_eta)*-1.0d0
        j_N(4)=1.0d0/4.0d0*(1.0d0-j_xi)*(1.0d0+j_eta)*-1.0d0
        
        call mat_zero(peri_B,2,16)
       
        do i=0,3
            peri_B(1,i*2+1)=i_N(i+1)
            peri_B(1,i*2+9)=j_N(i+1)
            peri_B(2,i*2+2)=i_N(i+1)
            peri_B(2,i*2+10)=j_N(i+1)
        enddo
    
        ! Calculate current particle position coordinates
        i_realcoord(1)=0.0d0
        i_realcoord(2)=0.0d0
        j_realcoord(1)=0.0d0
        j_realcoord(2)=0.0d0
    
        DO i=1,4
            i_realcoord(1)=i_realcoord(1)+i_N(i)*x_ref(1,i)
            i_realcoord(2)=i_realcoord(2)+i_N(i)*x_ref(2,i)
        
            j_realcoord(1)=j_realcoord(1)+j_N(i)*x_ref(1,i+4)
            j_realcoord(2)=j_realcoord(2)+j_N(i)*x_ref(2,i+4)
        ENDDO
        
        ! Distance between interaction points
        kexi=sqrt((i_realcoord(1)+j_realcoord(1))
     &  *(i_realcoord(1)+j_realcoord(1))
     &  +(i_realcoord(2)+j_realcoord(2))
     &  * (i_realcoord(2)+j_realcoord(2)))

        L=ele_esize*param_2
        peri_Jdet=i_Jdet*j_Jdet*DEXP(-1.0d0*kexi/(L))*param_1
        
        peri_D(1,1)=(i_realcoord(1)+j_realcoord(1))
     &               *(i_realcoord(1)+j_realcoord(1))*peri_Jdet
        peri_D(1,2)=(i_realcoord(1)+j_realcoord(1))
     &               *(i_realcoord(2)+j_realcoord(2))*peri_Jdet
        peri_D(2,1)=(i_realcoord(2)+j_realcoord(2))
     &               *(i_realcoord(1)+j_realcoord(1))*peri_Jdet
        peri_D(2,2)=(i_realcoord(2)+j_realcoord(2))
     &               *(i_realcoord(2)+j_realcoord(2))*peri_Jdet

        call mat_zero(peri_DB,2,16)
        call mat_zero(block_stiff,16,16)
        
        peri_DB=matmul(peri_D,peri_B)
        call mat_tran(peri_B,tran_peri_B,2,16)
        block_stiff=matmul(tran_peri_B,peri_DB)
      end subroutine calc_gauss_intg  

!-----------------------------------------------------------------------------
!                     Matrix Multiplication
!-----------------------------------------------------------------------------    
      subroutine mat_mult(mat_1,mat_2,mat_3,row_1,col_1,row_2,col_2)
        Integer row_1,col_1,row_2,col_2,i,j,k
        Double Precision mat_1(row_1,col_1),mat_2(row_2,col_2),
     &                       mat_3(row_1,col_2)
    
        do i=1,row_1
            do j=1,col_2
                mat_3(i,j)=0.0d0
            enddo
        enddo
    
        do i=1,row_1
            do j=1,col_2
                do k=1,col_1
                    mat_3(i,j)=mat_3(i,j)+mat_1(i,k)*mat_2(k,j)
                enddo
            enddo
        enddo
    
        end subroutine mat_mult
    
!-----------------------------------------------------------------------------
!                    Matrix Transposition
!-----------------------------------------------------------------------------     
      subroutine mat_tran(mat_1,mat_2,row_1,col_1)
        Integer row_1,col_1,i,j
        Double Precision mat_1(row_1,col_1),mat_2(col_1,row_1)
    
        do i=1,row_1
            do j=1,col_1
                mat_2(j,i)=mat_1(i,j)
            enddo
        enddo
    
        end subroutine mat_tran

!-----------------------------------------------------------------------------
!                        Matrix Initialization
!-----------------------------------------------------------------------------      
      subroutine mat_zero(mat_1,row_1,col_1)
        integer row_1,col_1,i,j
        double precision mat_1(row_1,col_1)
    
        do i=1,row_1
            do j=1,col_1
                mat_1(i,j)=0.0d0
            enddo
        enddo
    
      end subroutine mat_zero
       
!------------------------------------------------------------------------------
!                           Function to Calculate Gauss Point Coordinates
!------------------------------------------------------------------------------   
       subroutine gauss_loc(xCur,loc_gauss)
       integer i,j
       double precision xCur(2,8),loc_gauss(2,8),x_Gspoint(4),
     &                  y_Gspoint(4),u_N(4)

        call mat_zero(loc_gauss,2,8)
        x_Gspoint(1)=-0.577350269189626d0
        x_Gspoint(2)=0.577350269189626d0
        x_Gspoint(3)=0.577350269189626d0
        x_Gspoint(4)=-0.577350269189626d0
    
        y_Gspoint(1)=-0.577350269189626d0
        y_Gspoint(2)=-0.577350269189626d0
        y_Gspoint(3)=0.577350269189626d0
        y_Gspoint(4)=0.577350269189626d0
       
       do i=1,4
            u_N(1)=1.0d0/4.0d0*(1.0d0-x_Gspoint(i))*(1.0d0-y_Gspoint(i))
            u_N(2)=1.0d0/4.0d0*(1.0d0+x_Gspoint(i))*(1.0d0-y_Gspoint(i))
            u_N(3)=1.0d0/4.0d0*(1.0d0+x_Gspoint(i))*(1.0d0+y_Gspoint(i))
            u_N(4)=1.0d0/4.0d0*(1.0d0-x_Gspoint(i))*(1.0d0+y_Gspoint(i))

            do j=1,4
                loc_gauss(1,i)=loc_gauss(1,i)+u_N(j)*xCur(1,j)
                loc_gauss(2,i)=loc_gauss(2,i)+u_N(j)*xCur(2,j)
                loc_gauss(1,i+4)=loc_gauss(1,i+4)+u_N(j)*xCur(1,j+4)
                loc_gauss(2,i+4)=loc_gauss(2,i+4)+u_N(j)*xCur(2,j+4)
            enddo
       enddo
       
       end subroutine gauss_loc
