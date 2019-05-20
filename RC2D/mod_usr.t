module mod_usr
  use mod_mhd
  implicit none
  double precision :: q_e, parb,unit_currentdensity

  double precision, allocatable :: pbc(:),rbc(:)
  double precision :: usr_grav
  double precision :: heatunit,gzone,SRadius,bQ0,dya
  double precision, allocatable :: pa(:),ra(:),ya(:),Ta(:)
  double precision, allocatable :: bQa(:)
  integer, parameter :: jmax=20000

  integer :: numxQ^D,numFh,numFL,numLP
  integer :: iFshock,iFSmin,iFSmax
  double precision :: xFhmin,xFhmax,dFh
  double precision, allocatable :: xFh(:,:)
  double precision :: xQmin^D,xQmax^D,dxQ^D
  double precision, allocatable :: xQ(:^D&,:),Qe(:^D&)
  integer,allocatable :: numRL(:),numRR(:)

  integer :: numN
  double precision :: Nmin,Nmax,dNlog

  double precision :: t_update_Qe=-1.0e-7,dt_update_Qe=5.0e-4

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    call set_coordinate_system("Cartesian_2.5D")

    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3

    usr_set_parameters      => initglobaldata_usr
    usr_init_one_grid       => initonegrid_usr
    usr_special_bc          => specialbound_usr
    usr_aux_output          => specialvar_output
    usr_add_aux_names       => specialvarnames_output 
    usr_set_B0              => specialset_B0
    usr_set_J0              => specialset_J0
    usr_special_convert     => usrspecial_convert
    usr_special_resistivity => special_eta
    usr_var_for_errest      => p_for_errest
    usr_process_global      => special_global
    usr_gravity             => gravity
    usr_source              => special_source
    usr_refine_grid         => special_refine_grid



    call mhd_activate()
    parb=20.d0/3.d0
    ! unit of current density
    unit_currentdensity=unit_magneticfield/unit_length/4.d0/dpi
    ! unit of charge
    q_e=unit_currentdensity/unit_numberdensity/unit_velocity
    if(mype==0) print*,'unit of charge',q_e
    ! dimensionless charge of electron
    q_e=1.60217653d-19/q_e
    if(mype==0) print*,'dimensionless e',q_e

  end subroutine usr_init

  subroutine initglobaldata_usr()
    use mod_global_parameters

    heatunit=unit_pressure/unit_time          ! 3.697693390805347E-003 erg*cm^-3/s

    usr_grav=-2.74d4*unit_length/unit_velocity**2 ! solar gravity
    !bQ0=1.d3/heatunit ! background heating power density
    bQ0=1.0d-3/heatunit
    gzone=0.2d0 ! thickness of a ghostzone below the bottom boundary
    dya=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax) ! cells size of high-resolution 1D solar atmosphere
    Busr=Busr/unit_magneticfield ! magnetic field strength at the bottom
    SRadius=69.61d0 ! Solar radius
    ! hydrostatic vertical stratification of density, temperature, pressure
    call inithdstatic
    call get_bQa()

    ! for fast electron heating
    call init_Bfield()

  end subroutine initglobaldata_usr

  subroutine inithdstatic
  !! initialize the table in a vertical line through the global domain
    use mod_global_parameters

    integer :: j,na,nb,ibc
    double precision, allocatable :: gg(:)
    double precision:: rpho,Ttop,Tpho,wtra,res,rhob,pb,htra,Ttr,Fc,invT,kappa

    integer :: n_val=49, i
    double precision :: h_val(49),t_val(49)

    double precision :: ai(49),bi(49),ci(49),di(49)
    double precision :: hi


    !rpho=0.71d15/unit_numberdensity ! number density at the bottom relaxla
    rpho=0.3d15/unit_numberdensity ! number density at the bottom relaxla
    Tpho=1.d4/unit_temperature ! temperature of chromosphere
    Ttop=2.d6/unit_temperature ! estimated temperature in the top
    htra=0.2543d0 ! height of initial transition region
    wtra=0.01d0 ! width of initial transition region 
    !Ttr=1.d5/unit_temperature ! lowest temperature of upper profile
    Ttr=4.47d5/unit_temperature ! lowest temperature of upper profile
    Fc=6*2.d5/heatunit/unit_length ! constant thermal conduction flux
    kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3

   !VAL-C
    data    h_val   / 0., 50, 100, 150, 250, &
                      350, 450, 515, 555, 605, &
                      655, 705, 755, 855, 905, &
                      980, 1065, 1180, 1280, 1380, &
                      1515, 1605, 1785, 1925, 1990, &
                      2016, 2050, 2070, 2080, 2090, &
                      2104, 2107, 2109, 2113, 2115, &
                      2120, 2129, 2160, 2200, 2230, &
                      2255, 2263, 2267, 2271, 2274, &
                      2280, 2290, 2298, 2543 /

    data    t_val   / 6420, 5840, 5455, 5180, 4780, &
                      4465, 4220, 4170, 4230, 4420, &
                      4730, 5030, 5280, 5650, 5755, &
                      5925, 6040, 6150, 6220, 6280, &
                      6370, 6440, 6630, 6940, 7160, &
                      7360, 7660, 7940, 8180, 8440, &
                      9500, 10700, 12300, 18500, 21000, &
                      22500, 23000, 23500, 24000, 24200, &
                      24500, 25500, 28000, 32000, 37000, &
                      50000, 89100, 141000, 447000 /

    h_val(1:n_val)=h_val(1:n_val)/1.d4
    t_val(1:n_val)=t_val(1:n_val)/1.d6

    allocate(ya(jmax),Ta(jmax),gg(jmax),pa(jmax),ra(jmax))


    !call cubic_interp(h_val,t_val,ai,bi,ci,di,n_val)

    !do j=1,jmax
    !  ya(j)=(dble(j)-0.5d0)*dya-gzone
    !  if (ya(j)<htra .and. ya(j)>=0.d0) then
    !    do i=1,n_val-1
    !      if (ya(j)>=h_val(i) .and. ya(j)<h_val(i+1)) then
    !        !Ta(j)=t_val(i)+(ya(j)-h_val(i))*(t_val(i+1)-t_val(i))/(h_val(i+1)-h_val(i))
    !        hi=ya(j)-h_val(i)
    !        Ta(j)=ai(i)+bi(i)*hi+ci(i)*hi**2+di(i)*hi**3
    !        exit
    !      endif
    !    enddo
    !  else if (ya(j)<0.d0) then
    !    Ta(j)=t_val(1)+(ya(j)-h_val(1))*(t_val(2)-t_val(1))/(h_val(2)-h_val(1))
    !  else
    !    Ta(j)=(3.5d0*Fc/kappa*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
    !  endif
    !  gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
    !enddo


    do j=1,jmax
      ya(j)=(dble(j)-0.5d0)*dya-gzone
      if(ya(j)>=htra) then
        Ta(j)=(3.5d0*Fc/kappa*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
      else
        do i=1,n_val
          if (ya(j)<h_val(i+1)) then
            Ta(j)=t_val(i)+(ya(j)-h_val(i))*(t_val(i+1)-t_val(i))/(h_val(i+1)-h_val(i))
            exit
          endif
        enddo
      endif
      gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
    enddo

    !! solution of hydrostatic equation 
    nb=int(gzone/dya)
    ra(1)=rpho
    pa(1)=rpho*Ta(1)
    invT=gg(1)/Ta(1)
    invT=0.d0
    do j=2,jmax
       invT=invT+(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
       pa(j)=pa(1)*dexp(invT*dya)
       ra(j)=pa(j)/Ta(j)
    end do
    !! initialized rho and p in the fixed bottom boundary
    na=floor(gzone/dya+0.5d0)
    res=gzone-(dble(na)-0.5d0)*dya
    rhob=ra(na)+res/dya*(ra(na+1)-ra(na))
    pb=pa(na)+res/dya*(pa(na+1)-pa(na))
    allocate(rbc(nghostcells))
    allocate(pbc(nghostcells))
    do ibc=nghostcells,1,-1
      na=floor((gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0))/dya+0.5d0)
      res=gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0)-(dble(na)-0.5d0)*dya
      rbc(ibc)=ra(na)+res/dya*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+res/dya*(pa(na+1)-pa(na))
    end do

    if (mype==0) then
     print*,'minra',minval(ra)
     print*,'rhob',rhob
     print*,'pb',pb
    endif

  end subroutine inithdstatic

  subroutine cubic_interp(xi,yi,ai,bi,ci,di,ni)
    use mod_global_parameters

    integer :: ni,i
    double precision :: xi(ni),yi(ni)
    double precision :: ai(ni),bi(ni),ci(ni),di(ni)

    double precision :: xnew(ni+2),ynew(ni+2)
    double precision :: matrix1(ni+2,ni+2),ri1(ni+2)
    double precision :: matrix2(ni+2,ni+2),ri2(ni+2)
    double precision :: mi(ni+2),hi(ni+2)


    xnew(1:ni)=xi
    ynew(1:ni)=yi
    xnew(ni+1)=0.3
    ynew(ni+1)=0.7
    xnew(ni+2)=0.4
    ynew(ni+2)=0.95


    ! build equations
    matrix1=0.d0
    matrix2=0.d0

    do i=1,ni+1
      hi(i)=xnew(i+1)-xnew(i)
    enddo

    matrix1(1,1)=1.0
    matrix1(ni+2,ni+2)=1.0
    
    !matrix1(1,1)=2*hi(1)
    !matrix1(2,1)=hi(1)
    !matrix1(ni-1,ni)=hi(ni-1)
    !matrix1(ni,ni)=2*hi(ni-1)

    ri1(1)=0.d0
    ri1(ni+2)=0.d0
    do i=2,ni+1
      matrix1(i-1,i)=hi(i-1)
      matrix1(i,i)=2*(hi(i-1)+hi(i))
      matrix1(i+1,i)=hi(i)
      ri1(i)=6.d0*((ynew(i+1)-ynew(i))/hi(i)-(ynew(i)-ynew(i-1))/hi(i-1))
    enddo


    ! solve equations
    do i=1,ni+2
      matrix2(i,i)=1.d0
    enddo

    matrix2(2,1)=matrix1(2,1)/matrix1(1,1)
    ri2(1)=ri1(1)/matrix1(1,1)

    do i=2,ni+1
      matrix2(i+1,i)=matrix1(i+1,i)/(matrix1(i,i)-matrix1(i-1,i)*matrix2(i,i-1))
    enddo

    do i=2,ni+2
      ri2(i)=ri1(i)-matrix1(i-1,i)*ri2(i-1)
      ri2(i)=ri2(i)/(matrix1(i,i)-matrix1(i-1,i)*matrix2(i,i-1))
    enddo

    mi(ni+2)=ri2(ni+2)
    do i=ni+1,1,-1
      mi(i)=ri2(i)-matrix2(i+1,i)*mi(i+1)
    enddo


    ! get parameters for interpolation
    ai=yi
    ci=mi(1:ni)
    do i=1,ni-1
      bi(i)=(yi(i+1)-yi(i))/hi(i)-hi(i)*mi(i)/2.0-hi(i)*(mi(i+1)-mi(i))/6.0
      di(i)=(mi(i+1)-mi(i))/(6.0*hi(i))
    enddo


    !if (mype==0) print *, di

  end subroutine cubic_interp

  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    double precision                :: ggrid(ixI^S)

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,2)=ggrid(ixO^S)

  end subroutine gravity

  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*(SRadius/(SRadius+x(ixO^S,2)))**2
  end subroutine

  subroutine init_Bfield()
    use mod_global_parameters

    integer :: ix^D,refine_factor
    double precision :: lengthFL

    refine_factor=2**(refine_max_level-1)
    dFh=(xprobmax2-xprobmin2)/(domain_nx2*refine_factor)

    dt_update_Qe=dFh/5.0

    xFhmin=1.0
    xFhmax=5.0
    numFh=floor((xFhmax-xFhmin)/dFh)+1
    numFL=numFh
    lengthFL=10.0
    numLP=floor(lengthFL/dFh)
        
    xQmin1=-3.0
    xQmax1=3.0
    xQmin2=0.0
    xQmax2=6.0
    dxQ1=dFh
    dxQ2=dFh
    numXQ1=floor((xQmax1-xQmin1)/dxQ1)+1
    numXQ2=floor((xQmax2-xQmin2)/dxQ2)+1

    allocate(xFh(numFh,ndim))
    allocate(xQ(numXQ1,numXQ2,ndim),Qe(numXQ1,numXQ2))
    allocate(numRL(numFL),numRR(numFL))

    do ix2=1,numFh
      xFh(ix2,1)=0.0
      xFh(ix2,2)=xFhmin+(ix2-1)*dFh
    enddo

    do ix1=1,numXQ1
      do ix2=1,numXQ2
        xQ(ix1,ix2,1)=xQmin1+(ix1-1)*dxQ1
        xQ(ix1,ix2,2)=xQmin2+(ix2-1)*dxQ2
      enddo
    enddo

  end subroutine init_Bfield

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: res,Bf(ixI^S,1:ndir)
    integer :: ix^D,na
    logical, save :: first=.true.

    if(first)then
      if(mype==0) then
        write(*,*)'Simulating 2.5D solar atmosphere'
      endif
      first=.false.
    endif
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
        na=floor((x(ix^D,2)-xprobmin2+gzone)/dya+0.5d0)
        res=x(ix^D,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dya
        w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
        w(ix^D,p_)  =pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
    {end do\}
    w(ixO^S,mom(:))=zero
    if(B0field) then
      w(ixO^S,mag(:))=0.d0
    else
      call specialset_B0(ixI^L,ixO^L,x,Bf)
      w(ixO^S,mag(1:ndir))=Bf(ixO^S,1:ndir)
    endif
    if(mhd_glm) w(ixO^S,psi_)=0.d0
    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: lQgrid(ixI^S),bQgrid(ixI^S)

    ! add global background heating bQ
    call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)

    !! add localized external heating lQ concentrated at feet of loops
    if(iprob > 2)then
      call getlQ(lQgrid,ixI^L,ixO^L,qtC,wCT,x)
      w(ixO^S,e_)=w(ixO^S,e_)+qdt*lQgrid(ixO^S)
    endif

  end subroutine special_source

  subroutine getbQ(bQgrid,ixI^L,ixO^L,qt,w,x)
  ! calculate background heating bQ
    use mod_global_parameters
    use mod_radiative_cooling

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)

    double precision :: bQgrid(ixI^S),ens(ixI^S)
    double precision :: res
    integer :: na,ix^D

    !bQgrid(ixO^S)=bQ0*dexp(-x(ixO^S,2)/6.d0)
   
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      if (x(ix^D,2)>0.0) then
        bQgrid(ix^D)=bQ0*(abs(x(ix^D,2)-0.1)/0.5)**(-2.7)
        bQgrid(ix^D)=bQgrid(ix^D)/(exp(0.3/(x(ix^D,2)-0.2))-1.0)
        if (bQgrid(ix^D)<0.0) bQgrid(ix^D)=0.0
      endif
    {end do\}

    !{do ix^DB=ixOmin^DB,ixOmax^DB\}
    !    na=floor((x(ix^D,2)-xprobmin2+gzone)/dya+0.5d0)
    !    res=x(ix^D,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dya
    !    bQgrid(ix^D)=bQa(na)+(one-cos(dpi*res/dya))/two*(bQa(na+1)-bQa(na))
    !{end do\}

  end subroutine getbQ

  subroutine get_bQa()
    use mod_global_parameters
    use mod_radiative_cooling

    integer :: j
    double precision :: l1,gradt
    double precision :: tc_k_para,ke,dcdh(jmax)
    double precision :: lQa(jmax),rad(jmax),cond(jmax)
    character(20) :: fname
   
    if(si_unit) then
      ! Spitzer thermal conductivity with SI units
      tc_k_para=8.d-12*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3 
    else
      ! Spitzer thermal conductivity with cgs units
      tc_k_para=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3 
    end if


    allocate(bQa(jmax))


    ! radiative cooling
    do j=1,jmax
      if (Ta(j)<=tcoolmin .or. ya(j)<=0.0) then
        l1=zero
      else if (Ta(j)>=tcoolmax) then
        l1=lcool(4000)*sqrt(Ta(j)/tcoolmax)
        l1=l1*(ra(j)**2)
      else
        call findl(Ta(j),l1)
        l1=l1*(ra(j)**2)+smalldouble
      endif
      rad(j)=l1
      bQa(j)=l1
    enddo


    ! themal conduction
    do j=1,jmax
      gradt=0.0
      if (j==1) then
        gradt=(Ta(j+1)-Ta(j))/(ya(j+1)-ya(j))
      else if (j==jmax) then
        gradt=(Ta(j)-Ta(j-1))/(ya(j)-ya(j-1))
      else
        gradt=(Ta(j+1)-Ta(j-1))/(ya(j+1)-ya(j-1))
      endif

      ke=tc_k_para*dsqrt(Ta(j))**5
      cond(j)=ke*gradt
    enddo


    ! heating owing to thermal conduction
    do j=1,jmax-1
      dcdh(j)=(cond(j+1)-cond(j))/(ya(j+1)-ya(j))
    enddo

    if (mype==0) then
      fname='dEdt.txt'
      open(1,file=fname)
      write(1,*) jmax
      write(1,*) 'y rho Te cooling conduction dcdh'
      do j=1,jmax
          write(1,'(e15.7, e15.7, e15.7, e15.7, e15.7, e15.7)') ya(j), &
                ra(j),Ta(j),rad(j),cond(j),dcdh(j)
      enddo
      close(1)
    endif


  end subroutine get_bQa

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    use mod_global_parameters

    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pth(ixI^S),tmp(ixI^S),ggrid(ixI^S),invT(ixI^S)
    double precision :: delydelx, Bf(ixI^S,1:ndir)
    integer :: ix^D,idir,ixInt^L

    select case(iB)
    case(3)
      ! fixed zero velocity
      do idir=1,ndir
        w(ixO^S,mom(idir))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(idir))&
                   /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
      end do
      ! fixed b1 b2 b3
      if(B0field) then
        w(ixO^S,mag(:))=0.d0
      else
        call specialset_B0(ixI^L,ixO^L,x,Bf)
        w(ixO^S,mag(1:ndir))=Bf(ixO^S,1:ndir)
      endif
      ! fixed gravity stratification of density and pressure pre-determined in
      ! initial condition
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,rho_)=rbc(ix2)
        w(ixOmin1:ixOmax1,ix2,p_)=pbc(ix2)
      enddo
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      ixInt^L=ixO^L;
      ixIntmin2=ixOmin2-1;ixIntmax2=ixOmin2-1;
      call mhd_get_pthermal(w,x,ixI^L,ixInt^L,pth)
      ixIntmin2=ixOmin2-1;ixIntmax2=ixOmax2;
      call getggrav(ggrid,ixI^L,ixInt^L,x)
      ! fill pth, rho ghost layers according to gravity stratification
      invT(ixOmin2-1^%2ixO^S)=w(ixOmin2-1^%2ixO^S,rho_)/pth(ixOmin2-1^%2ixO^S)
      tmp=0.d0
      do ix2=ixOmin2,ixOmax2
        tmp(ixOmin2-1^%2ixO^S)=tmp(ixOmin2-1^%2ixO^S)+0.5d0*&
            (ggrid(ix2^%2ixO^S)+ggrid(ix2-1^%2ixO^S))*invT(ixOmin2-1^%2ixO^S)
        w(ix2^%2ixO^S,p_)=pth(ixOmin2-1^%2ixO^S)*dexp(tmp(ixOmin2-1^%2ixO^S)*dxlevel(2))
        w(ix2^%2ixO^S,rho_)=w(ix2^%2ixO^S,p_)*invT(ixOmin2-1^%2ixO^S)
      enddo
      ! fixed zero velocity
      do idir=1,ndir
        w(ixO^S,mom(idir))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,mom(idir))&
                     /w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,rho_)
      end do
      ! zero normal gradient extrapolation
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,mag(:))=(1.0d0/3.0d0)* &
                    (-w(ixOmin1:ixOmax1,ix2-2,mag(:))&
               +4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(:)))
      enddo
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

  subroutine p_for_errest(ixI^L,ixO^L,iflag,w,x,var)
    integer, intent(in)           :: ixI^L,ixO^L,iflag
    double precision, intent(in)  :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(out) :: var(ixI^S)

    call mhd_get_pthermal(w,x,ixI^L,ixO^L,var)
    
  end subroutine p_for_errest

  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
  ! Enforce additional refinement or coarsening
  ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! fix the bottom layer to the 4 level
    if (iprob==1) then
      if (any(x(ixO^S,2)<=xprobmin2+0.3d0)) then
        if (level<4) then
          refine=1
          coarsen=-1
        else
          refine=-1
          coarsen=-1
        endif
      else
        refine=-1
        coarsen=1
      endif
    endif

  end subroutine special_refine_grid

  subroutine special_global(iit,qt)
    !! run at the begining of a step
    use mod_global_parameters

    integer, intent(in) :: iit
    double precision, intent(in) :: qt
    character (20) :: fname
    integer :: ixFL,iyLH
    double precision :: t1,t2

    if (iprob==3 .and. convert) then
      call get_flare_eflux()
    endif

    if (iprob==3) then
      if (qt>t_update_Qe+10*dt_update_Qe) then
        t_update_Qe=qt-1.0e-7
      endif

      if (qt>t_update_Qe) then
        t1=mpi_wtime()
        call get_flare_eflux()
        t_update_Qe=t_update_Qe+dt_update_Qe
        t2=mpi_wtime()
        if (mype==0) print *, iit,t2-t1
      endif

    endif

  end subroutine  

  subroutine get_flare_eflux()
    ! calculate electron deposit energy flux
    use mod_global_parameters

    integer :: ix^D
    double precision, allocatable :: xFL(:,:,:),xFR(:,:,:)
    double precision, allocatable :: NpL(:,:),NpR(:,:)
    double precision, allocatable :: QeL(:,:),QeR(:,:)
    double precision, allocatable :: BvL(:,:,:),BvR(:,:,:)

    double precision :: ratio

    double precision, allocatable :: QeN(:)

    allocate(xFL(numFL,numLP,ndim),xFR(numFL,numLP,ndim))
    allocate(NpL(numFL,numLP),NpR(numFL,numLP))
    allocate(QeL(numFL,numLP),QeR(numFL,numLP))
    allocate(BvL(numFL,numLP,ndim),BvR(numFL,numLP,ndim))

    xFL(:,1,1)=xFh(:,1)
    xFL(:,1,2)=xFh(:,2)
    xFR(:,1,1)=xFh(:,1)
    xFR(:,1,2)=xFh(:,2)


    ! trace Bfield to prapare for heating
    call update_Bfield(xFL,xFR,NpL,NpR,BvL,BvR)


    ! find the shock to deside the location of fast electrons
    call locate_shock(xFL(:,1,:),NpL(:,1),BvL(:,1,:),ratio)


    ! update heating table based of compression ratio of shock
    numN=1001
    allocate(QeN(numN))
    Nmin=1.0e18
    Nmax=1.0e23
    dNlog=log10(Nmax/Nmin)/(numN-1.0)    
    call update_heating_table(QeN,ratio)


    ! get heating rate for each field line
    call get_heating_rate(xFL,NpL,BvL,QeL,QeN,numRL)
    call get_heating_rate(xFR,NpR,BvR,QeR,QeN,numRR)


    ! get local heating rate via interpolation
    call get_Qe(xFL,xFR,QeL,QeR)


    deallocate(QeN)
    deallocate(xFL,xFR,NpL,NpR,BvL,BvR,QeL,QeR)

  end subroutine get_flare_eflux

  subroutine update_Bfield(xFL,xFR,NpL,NpR,BvL,BvR)
    ! trace B field and get column depth
    use mod_global_parameters
    use mod_usr_methods
    use mod_trace_Bfield

    double precision :: xFL(numFL,numLP,ndim),xFR(numFL,numLP,ndim)
    double precision :: NpL(numFL,numLP),NpR(numFL,numLP)
    double precision :: BvL(numFL,numLP,ndim),BvR(numFL,numLP,ndim)

    integer :: ix^D,dirct
    logical :: forward

    ! left foot
    forward=.false.
    dirct=0
    do ix1=1,numFL
      call trace_Bfield(xFL(ix1,:,:),NpL(ix1,:),BvL(ix1,:,:),dFh,dirct,numLP,numRL(ix1),forward)
    enddo

    !right foot
    forward=.true.
    dirct=0
    do ix1=1,numFL
      call trace_Bfield(xFR(ix1,:,:),NpR(ix1,:),BvR(ix1,:,:),dFh,dirct,numLP,numRR(ix1),forward)
    enddo    

  end subroutine update_Bfield

  subroutine locate_shock(xFh,Nph,Bvh,ratio)
    ! find the location of shock and return index and compression ratio
    use mod_global_parameters
    use mod_usr_methods

    double precision :: xFh(numFL,ndim),Nph(numFL),Bvh(numFL,ndim)
    double precision :: ratio

    integer :: iFL,numCR
    double precision :: deltaB1,widthCR,rhoUp,rhoDw


    ! find shock front based on the change of magnetic field
    deltaB1=0.0
    do iFL=1,numFL-1
      if (Bvh(iFL,1)-Bvh(iFL+1,1)>deltaB1) then
        deltaB1=Bvh(iFL,1)-Bvh(iFL+1,1)
        iFshock=iFL
      endif
    enddo


    widthCR=0.2
    numCR=floor(widthCR/dFh)
    rhoUp=0.0
    rhoDw=0.0
    if (iFshock>numCR .and. iFshock<numFL-numCR) then

      do iFL=iFshock-numCR,iFshock-1
        rhoDw=rhoDw+Nph(iFL)
      enddo

      do iFL=iFshock+1,iFshock+numCR      
        rhoUp=rhoUp+Nph(iFL)
      enddo

      ratio=rhoDw/rhoUp
    else
      print *, 'unable to calculate compression ratio'
      print *, iFshock,numCR,numFL
    endif

  end subroutine locate_shock

  subroutine update_heating_table(QeN,ratio)
    ! calculate heating table based on shock compression ratio
    use mod_global_parameters
    use mod_usr_methods

    double precision :: QeN(numN),Ncol(numN)
    double precision :: ratio   

    integer :: j,iBx,numBx,iBc,iNcol
    double precision :: Atom,temp
    double precision :: Cllog,gammaN,beta,K,Nc,Nb,Efactor
    double precision :: dtx,tx,au,al,b,maxh
    double precision, allocatable :: Bxc(:)
    double precision :: Ec,delta,Ec_erg,mu0


    ! electron spectrum [electrons cm^-2 keV-1]
    delta=(ratio+2)/(ratio-1)/2.0 ! spectral index
    Ec=20.0  ! cutoff energy [keV]

    if (mype==0) print *,'r',ratio,'delta',delta

    ! parameters for electron deposition
    Cllog=25.0  ! Coulomb logarithns
    gammaN=Cllog  ! full ionized plasma
    beta=2.0    ! full ionized plasma
    mu0=0.732   ! cos pitch angle
    K=2.0*dpi*const_e**4
    Ec_erg=Ec*1.0e3*const_ev
    Nc=mu0*Ec_erg**2/((2.0+beta/2.0)*gammaN*K)
    Efactor=K*gammaN*(delta-2.0)/(2*mu0*Ec_erg**2)


    ! beta function
    numBx=100001
    allocate(Bxc(numBx))
    dtx=1.0/(numBx-1)
    Bxc=0
    al=delta/2.0
    b=1.0/3
    do iBx=2,numBx-1
      tx=dtx*(iBx-1.0)
      Bxc(iBx)=Bxc(iBx-1) + (tx**(al-1.0))*((1.0-tx)**(b-1.0))*dtx
    enddo
    Bxc(numBx)=Bxc(numBx-1)


    do iNcol=1,numN
      Ncol(iNcol)=Nmin*10**(dNlog*(iNcol-1.0))
      iBc=int((numBx-1)*Ncol(iNcol)/Nc+0.5)+1
      if (iBc>numBx)  then
        iBc=numBx
      endif
      QeN(iNcol)=Efactor*Bxc(iBc)*((Ncol(iNcol)/Nc)**(-delta/2.0))
    enddo
    deallocate(Bxc)

    QeN=QeN*1.0d10

  end subroutine update_heating_table

  subroutine get_heating_rate(xF,Np,Bv,QeLR,QeN,numR)
    ! calculate local heating rate
    use mod_global_parameters
    use mod_usr_methods

    integer :: numR(numFL)
    double precision :: xF(numFL,numLP,ndim),QeLR(numFL,numLP)
    double precision :: Bv(numFL,numLP,ndim)
    double precision :: Np(numFL,numLP),Ncol(numFL,numLP)
    double precision :: QeN(numN)

    integer :: ix^D,iNlog,iFscale
    double precision :: dl^D,dl,Bratio,Bxy0,Bxy,mu0,mu
    double precision :: Hscale,Qratio,const


    mu0=0.732 ! initial pitch angle

    ! deside which field line has fast electrons
    Hscale=0.5
    iFscale=floor(Hscale/dFh/2.0)
    iFSmin=iFshock-iFscale*2
    iFSmax=iFshock+iFscale*2
    if (iFSmin<1) print *, 'wrong iFmin', iFSmin
    if (iFSmax>numFL) print *, 'wrong iFmax', iFSmax


    ! calculate heating rate for each field line
    Ncol=0
    QeLR=0

    do ix1=iFSmin,iFSmax
      Loopline: do ix2=2,numR(ix1)

        ! special distribution of fast electrons
        Qratio=exp(-(ix1-iFshock)**2/(1.0*iFscale)**2)

        ! column depth
        dl1=xF(ix1,ix2,1)-xF(ix1,ix2-1,1)
        dl2=xF(ix1,ix2,2)-xF(ix1,ix2-1,2)
        dl=dsqrt(dl1**2+dl2**2)*unit_length
        Ncol(ix1,ix2)=Ncol(ix1,ix2-1)+dl*(Np(ix1,ix2)+Np(ix1,ix2+1))/2.0

        ! loop expension
        Bxy0=dsqrt(Bv(ix1,1,1)**2+Bv(ix1,1,2)**2)
        Bxy=dsqrt(Bv(ix1,ix2,1)**2+Bv(ix1,ix2,2)**2)
        Bratio=Bxy/Bxy0
        const=(1-mu0**2)/Bxy0
        if (const*Bxy>=1) then
          exit Loopline
        endif

        ! local heating rate
        iNlog=int(log10(Ncol(ix^D)/Nmin)/dNlog+0.5)+1
        if (iNlog<1) then
          iNlog=1
        endif
        QeLR(ix^D)=QeN(iNlog)*Np(ix^D)*Bratio*Qratio

      enddo Loopline
    enddo

  end subroutine get_heating_rate

  subroutine get_Qe(xFL,xFR,QeL,QeR)
    ! get local heating rate via interpolation
    use mod_global_parameters
    use mod_usr_methods

    double precision :: xFL(numFL,numLP,ndim),xFR(numFL,numLP,ndim)
    double precision :: QeL(numFL,numLP),QeR(numFL,numLP)

    integer :: ix^D,ipe,ihmin
    double precision :: hmin
    double precision,allocatable :: xb1(:,:),xb2(:,:),xb3(:,:),xb4(:,:)
    logical :: heating
    integer :: numP1,numP2,numP3,numP4


    Qe=0
    hmin=0.05 ! lower boundary for heating [10 Mm]
    ihmin=floor((hmin-xQ(1,1,2))/dxQ2)+1


    ! find boundaries for heating region
    numP1=numRL(iFSmax)
    numP2=numRL(iFSmin)
    numP3=numRR(iFSmax)
    numP4=numRR(iFSmin)
    allocate(xb1(numP1,ndim),xb2(numP2,ndim))
    allocate(xb3(numP3,ndim),xb4(numP4,ndim))
    call get_boundary(xb1,xFL(iFSmax,1:numP1,:),numP1)
    call get_boundary(xb2,xFL(iFSmin,1:numP2,:),numP2)
    call get_boundary(xb3,xFR(iFSmax,1:numP3,:),numP3)
    call get_boundary(xb4,xFR(iFSmin,1:numP4,:),numP4)


    ! check heating region and calculate heating rate
    do ix2=ihmin,numXQ2
      ipe=mod(ix2,npe)
      if (mype==ipe) then
        do ix1=1,numXQ1
          call check_heating(xQ(ix1,ix2,:),xb1,xb2,xb3,xb4, &
                             numP1,numP2,numP3,numP4,heating)
          if (heating) then
            if (xQ(ix1,ix2,1)<=0) then
              call interp_Qe(xQ(ix1,ix2,:),Qe(ix1,ix2),xFL,QeL,numRL)
            else
              call interp_Qe(xQ(ix1,ix2,:),Qe(ix1,ix2),xFR,QeR,numRR)
            endif
          endif
        enddo
      endif
    enddo


    !do ix2=ihmin,numXQ2    
    !  ipe=mod(ix2,npe)
    !  if (mype==ipe) then
    !    do ix1=1,numXQ1
    !      if (xQ(ix1,ix2,1)<=0) then
    !        call interp_Qe(xQ(ix1,ix2,:),Qe(ix1,ix2),xFL,QeL,numRL)
    !      else
    !        !call interp_Qe(xQ(ix1,ix2,:),Qe(ix1,ix2),xFR,QeR,numRR)
    !      endif
    !    enddo
    !  endif
    !enddo


    do ix2=ihmin,numXQ2
      ipe=mod(ix2,npe)
      call MPI_BCAST(Qe(:,ix2),numXQ1,MPI_DOUBLE_PRECISION,ipe,icomm,ierrmpi)
    enddo

    deallocate(xb1,xb2,xb3,xb4)

  end subroutine get_Qe

  subroutine get_boundary(xb,xf,numP)
    ! find boundary
    use mod_global_parameters
    use mod_usr_methods

    integer :: numP
    double precision :: xb(numP,ndim),xf(numP,ndim)

    integer :: ix^D,j

    do j=1,numP
      if (xf(j,1)>=xQ(1,1,1) .and. xf(j,1)<=xQ(numXQ1,1,1) .and. &
          xf(j,2)>=xQ(1,1,2) .and. xf(j,2)<=xQ(1,numXQ2,2)) then
        ix1=floor((xf(j,1)-xQ(1,1,1))/dxQ1+0.5)+1
        ix2=floor((xf(j,2)-xQ(1,1,2))/dxQ2+0.5)+1
        xb(j,1)=xQ(ix1,ix2,1)
        xb(j,2)=xQ(ix1,ix2,2)
      else
        if (mype==0) print *, 'Field line go out of xQe table!'
        if (mype==0) print *, xf(j,:)
        if (mype==0) print *, xQ(1,1,1),xQ(numXQ1,1,1),xQ(1,numXQ2,2),xQ(1,numXQ2,2)
      endif
    enddo

  end subroutine get_boundary

  subroutine check_heating(xc,xb1,xb2,xb3,xb4,numP1,numP2,numP3,numP4,heating)
    ! check whether or not the point is in heating region
    use mod_global_parameters
    use mod_usr_methods

    integer :: numP1,numP2,numP3,numP4
    double precision :: xc(ndim)
    double precision :: xb1(numP1,ndim),xb2(numP2,ndim)
    double precision :: xb3(numP3,ndim),xb4(numP4,ndim)
    logical :: heating

    integer :: timesL,timesR,t1,t2,t3,t4

    ! left across times
    call get_across_time(xc,xb1,numP1,0,t1)
    call get_across_time(xc,xb2,numP2,0,t2)
    call get_across_time(xc,xb3,numP3,0,t3)
    call get_across_time(xc,xb4,numP4,0,t4)
    timesL=t1+t2+t3+t4
    ! right across times
    call get_across_time(xc,xb1,numP1,1,t1)
    call get_across_time(xc,xb2,numP2,1,t2)
    call get_across_time(xc,xb3,numP3,1,t3)
    call get_across_time(xc,xb4,numP4,1,t4)
    timesR=t1+t2+t3+t4

    if (mod(timesL,2)==1 .and. mod(timesR,2)==1) then
      heating=.true.
    else
      heating=.false.
    endif

  end subroutine check_heating

  subroutine get_across_time(xc,curve,numc,LR,times)
    ! calculate that how many times that a parallel or a perpendicular
    ! line intersect a curve on the left of the given point
    ! given point xc is in the line.
    ! LR==0: left across
    ! LR==1: right across
    use mod_global_parameters

    integer :: numc,times,j,LR
    double precision :: xc(ndim),curve(numc,ndim)
    double precision :: dl,xit

    times=0

    do j=1,numc-1
      if (curve(j,2)==xc(2) .and. curve(j+1,2)/=xc(2)) then
        if (curve(j,1)<=xc(1) .and. LR==0) times=times+1
        if (curve(j,1)>=xc(1) .and. LR==1) times=times+1
      endif
    enddo


    !do j=1,numc-1
    !  if (curve(j,2)<=xc(2) .and. curve(j+1,2)>xc(2)) then
    !    dl=(xc(2)-curve(j,2))/(curve(j+1,2)-curve(j,2))
    !    xit=curve(j,1)*(1.0-dl)+curve(j+1,1)*dl
    !    if (xit<=xc(1) .and. LR==0) times=times+1
    !    if (xit>=xc(1) .and. LR==1) times=times+1
    !  endif

    !  if (curve(j,2)>=xc(2) .and. curve(j+1,2)<xc(2)) then
    !    dl=(xc(2)-curve(j+1,2))/(curve(j,2)-curve(j+1,2))
    !    xit=curve(j+1,1)*(1.0-dl)+curve(j,1)*dl
    !    if (xit<=xc(1) .and. LR==0) times=times+1
    !    if (xit>=xc(1) .and. LR==1) times=times+1
    !  endif
    !enddo

  end subroutine get_across_time


  subroutine interp_Qe(xc,lQc,xi,Qi,numR)
    ! get Qe for each height
    use mod_global_parameters
    use mod_usr_methods

    double precision :: xc(ndim),lQc
    double precision :: xi(numFL,numLP,ndim),Qi(numFL,numLP)
    integer :: numR(numFL)

    integer :: ix^D
    double precision :: xd3,x1n(2,2),x2n(2,2),xdn(2,2),Qn(2,2),dn
    double precision :: dyl,dyr,dxm
    double precision :: xil(2),xir(2),Ql,Qr
    logical :: heating
    integer :: timesL,timesR,t0,t1,t2
    integer :: factor
    double precision :: Dii(2,2),sumDii,logQe

    lQc=0

    factor=-3

    !looking for nearby points to do interpolation
    x1n(:,:)=xi(1,1,1)
    x2n(:,:)=xi(1,1,2)
    xdn=8.0*dFh
    Qn=0

    do ix1=1,numFL
      do ix2=1,numR(ix1)
        dn=sqrt((xi(ix1,ix2,1)-xc(1))**2 + (xi(ix1,ix2,2)-xc(2))**2)
        ! lower left
        if (xi(ix1,ix2,1)<=xc(1) .and. xi(ix1,ix2,2)<=xc(2) .and. &
            dn<xdn(1,1)) then
          x1n(1,1)=xi(ix1,ix2,1)
          x2n(1,1)=xi(ix1,ix2,2)
          Qn(1,1)=Qi(ix1,ix2)
          xdn(1,1)=dn
        endif
        ! lower right
        if (xi(ix1,ix2,1)>=xc(1) .and. xi(ix1,ix2,2)<=xc(2) .and. &
            dn<xdn(2,1)) then
          x1n(2,1)=xi(ix1,ix2,1)
          x2n(2,1)=xi(ix1,ix2,2)
          Qn(2,1)=Qi(ix1,ix2)
          xdn(2,1)=dn
        endif
        ! upper left
        if (xi(ix1,ix2,1)<=xc(1) .and. xi(ix1,ix2,2)>=xc(2) .and. &
            dn<xdn(1,2)) then
          x1n(1,2)=xi(ix1,ix2,1)
          x2n(1,2)=xi(ix1,ix2,2)
          Qn(1,2)=Qi(ix1,ix2)
          xdn(1,2)=dn
        endif
        ! upper left
        if (xi(ix1,ix2,1)>=xc(1) .and. xi(ix1,ix2,2)>=xc(2) .and. &
            dn<xdn(2,2)) then
          x1n(2,2)=xi(ix1,ix2,1)
          x2n(2,2)=xi(ix1,ix2,2)
          Qn(2,2)=Qi(ix1,ix2)
          xdn(2,2)=dn
        endif
      enddo
    enddo

    ! interpolation
    if (sum(Qn)>0) then
      Dii=0

      do ix1=1,2
        do ix2=1,2
          if (Qn(ix1,ix2)==0) then
            Dii(ix1,ix2)=0
          else if (xdn(ix1,ix2)==0) then
            Dii(ix1,ix2)=1.0e18
          else
            Dii(ix1,ix2)=(xdn(ix1,ix2))**factor
          endif
        enddo
      enddo

      sumDii=sum(Dii)
      Dii=Dii/sumDii
      logQe=0

      do ix1=1,2
        do ix2=1,2
          lQc=lQc+Qn(ix1,ix2)*Dii(ix1,ix2)
        enddo
      enddo
    endif

  end subroutine interp_Qe

  subroutine getlQ(lQgrid,ixI^L,ixO^L,qt,w,x)
    !!calculate localized heating at chromosphere
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision :: lQgrid(ixI^S)

    integer :: ix^D,ixO^D


    select case(iprob)    

    case(3)
    
    lQgrid(ixO^S)=0
    {do ixO^DB=ixOmin^DB,ixOmax^DB\}
      ix1=floor((x(ixO^D,1)-xQ(1,1,1))/dxQ1+0.5)+1
      ix2=floor((x(ixO^D,2)-xQ(1,1,2))/dxQ2+0.5)+1
      if (ix1>=1 .and. ix1<=numXQ1 .and. ix2>=1 .and. ix2<=numXQ2) then
        lQgrid(ixO^D)=Qe(ix1,ix2)
      endif
    {enddo\}

    end select

  end subroutine getlQ

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_global_parameters
    use mod_radiative_cooling

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision :: pth(ixI^S),B2(ixI^S),divb(ixI^S)
    double precision :: Btotal(ixI^S,1:ndir),current_o(ixI^S,3)
    integer :: idir,idirmin
    double precision :: lQgrid(ixI^S),ens(ixI^S)
    double precision :: t

    ! output temperature
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=pth(ixO^S)/w(ixO^S,rho_)
    if(B0field) then
      Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))+block%B0(ixI^S,1:ndir,0)
    else
      Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))
    endif
    ! B^2
    B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)
    ! output Alfven wave speed B/sqrt(rho)
    w(ixO^S,nw+2)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))
    ! output divB1
    call divvector(Btotal,ixI^L,ixO^L,divb)
    w(ixO^S,nw+3)=0.5d0*divb(ixO^S)/dsqrt(B2(ixO^S))/(^D&1.0d0/dxlevel(^D)+)
    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+4)=pth(ixO^S)*two/B2(ixO^S)
    ! output current
    call get_current(w,ixI^L,ixO^L,idirmin,current_o)
    w(ixO^S,nw+5)=current_o(ixO^S,1)
    w(ixO^S,nw+6)=current_o(ixO^S,2)
    w(ixO^S,nw+7)=current_o(ixO^S,3)
    ! output special resistivity eta
    call special_eta(w,ixI^L,ixO^L,idirmin,x,current_o,divb)
    w(ixO^S,nw+8)=divb(ixO^S)

    ! output heating rate
    t=global_time
    !call getlQ(lQgrid,ixI^L,ixO^L,t,w,x)
    !w(ixO^S,nw+9)=lQgrid(ixO^S)

    ! background heating
    call getbQ(ens,ixI^L,ixO^L,t,w,x)
    w(ixO^S,nw+9)=ens(ixO^S)

    ! store the cooling rate 
    if(mhd_radiative_cooling)call getvar_cooling(ixI^L,ixO^L,w,x,ens)
    w(ixO^S,nw+10)=ens(ixO^S)


  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te Alfv divB beta j1 j2 j3 eta bQ rad'
  end subroutine specialvarnames_output

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a time-independent background magnetic field
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    if (iprob==1) then
      wB0(ixO^S,1)=zero
      wB0(ixO^S,2)=Busr
      wB0(ixO^S,3)=zero
    else if (iprob>=2) then
      wB0(ixO^S,1)=zero
      wB0(ixO^S,2)=-Busr*dtanh(parb*x(ixO^S,1))
      wB0(ixO^S,3)=Busr/dcosh(parb*x(ixO^S,1))
    endif

  end subroutine specialset_B0

  subroutine specialset_J0(ixI^L,ixO^L,x,wJ0)
  ! Here add a time-independent background current density 
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wJ0(ixI^S,7-2*ndir:ndir)

    if (iprob==1) then
      wJ0(ixO^S,1)=zero
      wJ0(ixO^S,2)=zero
      wJ0(ixO^S,3)=zero
    else if (iprob>=2) then
      wJ0(ixO^S,1)=zero
      wJ0(ixO^S,2)=parb*Busr*dtanh(parb*x(ixO^S,1))/dcosh(parb*x(ixO^S,1))
      wJ0(ixO^S,3)=-parb*Busr/dcosh(parb*x(ixO^S,1))**2
    endif

  end subroutine specialset_J0

  subroutine special_eta(w,ixI^L,ixO^L,idirmin,x,current,eta)
    ! Set the common "eta" array for resistive MHD based on w or the
    ! "current" variable which has components between idirmin and 3.
    integer, intent(in) :: ixI^L, ixO^L, idirmin
    double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision :: current(ixI^S,7-2*ndir:3), eta(ixI^S)
    double precision :: rad(ixI^S),heta,reta,eta1,eta2,etam,vc,tar
    double precision :: jc,jabs(ixI^S)

    if (iprob==1) then
      eta(ixO^S)=0.d0
    else if (iprob>=2) then
      heta = 6.
      reta = 0.8d0 * 0.3d0
      eta1 = 0.002d0
      tar= 0.4d0
      !tar= 0.0d0
      vc=1.d-4
      eta2=2.d-4
      etam=2.d-2
      if (global_time<tar) then
        rad(ixO^S)=dsqrt(x(ixO^S,1)**2+(x(ixO^S,2)-heta)**2)
        where (rad(ixO^S) .lt. reta)
          eta(ixO^S)=eta1*(2.d0*(rad(ixO^S)/reta)**3-3.d0*(rad(ixO^S)/reta)**2+1.d0)
        elsewhere
          eta(ixO^S)=zero
        endwhere
      else
        rad(ixO^S)=dsqrt(sum(current(ixO^S,:)**2,dim=ndim+1))/w(ixO^S,rho_)/q_e
        where(rad(ixO^S)>vc)
          eta(ixO^S)=eta2*(rad(ixO^S)/vc-1.d0)
        elsewhere
          eta(ixO^S)=0.d0
        endwhere
        where(eta(ixO^S)>etam)
          eta(ixO^S)=etam
        endwhere
      end if 
    endif

    end subroutine special_eta

  subroutine usrspecial_convert(qunitconvert)
    integer, intent(in) :: qunitconvert
    character(len=20):: userconvert_type
  
    call spatial_integral_w
  end subroutine usrspecial_convert

  subroutine spatial_integral_w
    double precision :: dvolume(ixG^T), dsurface(ixG^T),timephy,dvone
    double precision, allocatable :: integral_ipe(:), integral_w(:)

    integer           :: nregions,ireg,ncellpe,ncell,idims,hxM^LL,nx^D
    integer           :: iigrid,igrid,status(MPI_STATUS_SIZE),ni
    character(len=100):: filename,region
    character(len=1024) :: line, datastr
    logical           :: patchwi(ixG^T),alive

    nregions=1
    ! number of integrals to perform
    ni=3
    allocate(integral_ipe(ni),integral_w(ni))
    integral_ipe=0.d0
    integral_w=0.d0
    nx^D=ixMhi^D-ixMlo^D+1;
    do ireg=1,nregions
      select case(ireg)
      case(1)
        region='fulldomain'
      case(2)
        region='cropped'
      end select
      ncellpe=0 
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        block=>ps(igrid)
        if(slab) then
          dvone={rnode(rpdx^D_,igrid)|*}
          dvolume(ixM^T)=dvone
          dsurface(ixM^T)=two*(^D&dvone/rnode(rpdx^D_,igrid)+)
        else
          dvolume(ixM^T)=ps(igrid)%dvolume(ixM^T)
          dsurface(ixM^T)= sum(ps(igrid)%surfaceC(ixM^T,:),dim=ndim+1)
          do idims=1,ndim
            hxM^LL=ixM^LL-kr(idims,^D);
            dsurface(ixM^T)=dsurface(ixM^T)+ps(igrid)%surfaceC(hxM^T,idims)
          end do
        end if
        ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
        patchwi(ixG^T)=.false.
        select case(region)
        case('cropped')
           call mask_grid(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,patchwi,ncellpe)
        case('fulldomain')
           patchwi(ixM^T)=.true.
           ncellpe=ncellpe+{nx^D*}
        case default
           call mpistop("region not defined")
        end select
        integral_ipe(1)=integral_ipe(1)+ &
                  integral_grid(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,dvolume,dsurface,1,patchwi)
        integral_ipe(2)=integral_ipe(2)+ &
                  integral_grid(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,dvolume,dsurface,2,patchwi)
        integral_ipe(3)=integral_ipe(3)+ &
                  integral_grid(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,dvolume,dsurface,3,patchwi)
      end do
      call MPI_ALLREDUCE(integral_ipe,integral_w,ni,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,icomm,ierrmpi)
      !call MPI_ALLREDUCE(ncellpe,ncell,1,MPI_INTEGER,MPI_SUM,icomm,ierrmpi)
      timephy=global_time
      if(mype==0) then
        write(filename,"(a,a,a)") TRIM(base_filename),TRIM(region),"mkc.csv"
        inquire(file=filename,exist=alive)
        if(alive) then
          open(unit=21,file=filename,form='formatted',status='old',access='append')
        else
          open(unit=21,file=filename,form='formatted',status='new')
          write(21,'(a)') 'time, emagnetic, einternal, current'
        endif
        write(datastr,'(es13.6, a)') timephy,','
        line=datastr
        write(datastr,"(es13.6, a)") integral_w(1),','
        line = trim(line)//trim(datastr)
        write(datastr,"(es13.6, a)") integral_w(2),','
        line = trim(line)//trim(datastr)
        write(datastr,"(es13.6)") integral_w(3)
        line = trim(line)//trim(datastr)
        write(21,'(a)') trim(line)
        close(21)
      endif
    enddo
    deallocate(integral_ipe,integral_w)
  end subroutine spatial_integral_w

  subroutine mask_grid(ixI^L,ixO^L,w,x,patchwi,cellcount)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    logical, intent(inout)             :: patchwi(ixG^T)

    double precision  ::  buff
    integer                            :: ix^D,cellcount

    buff=0.05d0*(xprobmax1-xprobmin1)
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
       if(x(ix^D,1)>xprobmin1+buff .and. x(ix^D,1)<xprobmax1-buff .and. &
          x(ix^D,2)>xprobmin2+buff .and. x(ix^D,2)<xprobmax2-buff) then
         patchwi(ix^D)=.true.
         cellcount=cellcount+1
       else
         patchwi(ix^D)=.false.
       endif
    {end do\}
    return

  end subroutine mask_grid

  function integral_grid(ixI^L,ixO^L,w,x,dvolume,dsurface,intval,patchwi)
    integer, intent(in)                :: ixI^L,ixO^L,intval
    double precision, intent(in)       :: x(ixI^S,1:ndim),dvolume(ixG^T),dsurface(ixG^T)
    double precision, intent(in)       :: w(ixI^S,nw)
    logical, intent(in) :: patchwi(ixG^T)
    
    double precision, dimension(ixG^T,1:ndir) :: bvec,qvec
    double precision :: current(ixG^T,7-2*ndir:3),tmp(ixG^T)
    double precision :: integral_grid,mcurrent
    integer :: ix^D,idirmin,idir,jdir,kdir

    integral_grid=0.d0
    select case(intval)
     case(1)
      ! magnetic energy
      if(B0field)then
        tmp(ixO^S)=0.5d0*sum((w(ixO^S,mag(:))+&
                      block%B0(ixO^S,:,0))**2,dim=ndim+1)
      else
        tmp(ixO^S)=0.5d0*sum(w(ixO^S,mag(:))**2,dim=ndim+1)
      endif
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         if(patchwi(ix^D)) integral_grid=integral_grid+tmp(ix^D)*dvolume(ix^D)
      {end do\}
     case(2)
      ! internal energy
      call mhd_get_pthermal(w,x,ixI^L,ixO^L,tmp)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         if(patchwi(ix^D))  integral_grid=integral_grid+tmp(ix^D)/(mhd_gamma-1.d0)*dvolume(ix^D)
      {end do\}
     case(3)
      ! current strength
      call get_current(w,ixI^L,ixO^L,idirmin,current)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         if(patchwi(ix^D)) integral_grid=integral_grid+dsqrt(sum(current(ix^D,:)**2))*dvolume(ix^D)
      {end do\}
     case default
         call mpistop("intval not defined")
    end select
    
    return
  end function integral_grid

end module mod_usr
