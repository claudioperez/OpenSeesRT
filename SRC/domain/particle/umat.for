c*****************************************************************
      ! Call UMAT
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     1 drplde,drpldt,stran,dstran,time,dtime,temp2,dtemp,predef,dpred,
     2 cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      use kvisual
      include 'aba_param.inc' 

      character*8 cmname
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1 ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),
     2 time(2),predef(1),dpred(1),props(nprops),coords(3),drot(3,3),
     3 dfgrd0(3,3),dfgrd1(3,3),jstep(4)

      INTEGER          i,j,cor_is
   
      
      ddsdde=0.0d0

      cor_is=0 ! Initialize related element variable
      do j=1,4
          ele_bond_sum_state_fem(j,1,noel)=0.0d0 
      enddo

      ! Calculate the number of related elements
      do i=1,100
          if(ele_bond_state_pd(1,1,i,noel).eq.1)then
              cor_is = cor_is+1
          endif
          do j=1,4
          ele_bond_sum_state_fem(j,1,noel)=
     &    ele_bond_sum_state_fem(j,1,noel)
     &     +ele_bond_sum_state_pd(j,1,i,noel)
          enddo
      enddo

      ! Total number of broken bonds in FEM element
      statev(1)=ele_bond_sum_state_fem(npt,1,noel)
      ! Damage in FEM element
      statev(2)=ele_bond_sum_state_fem(npt,1,noel)/cor_is/16
      return
      end


