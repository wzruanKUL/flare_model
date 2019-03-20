module mod_usr
  use mod_mhd
  implicit none
  double precision, allocatable :: pbc(:),rbc(:)
  double precision :: usr_grav
  double precision :: heatunit,gzone,theta,SRadius,kx,ly,bQ0,dya
  double precision, allocatable :: pa(:),ra(:),ya(:)
  integer, parameter :: jmax=8000
  double precision, allocatable :: xr(:), hr(:),rhor(:,:),Tr(:,:),pr(:,:)

  integer :: numX(ndim),numX^D,ixF^L
  double precision, allocatable :: xFL(:^D&,:),xFR(:^D&,:)
  double precision :: dxL(ndim),dxR(ndim)
  double precision, allocatable :: QeL(:^D&),QeR(:^D&)
  double precision, allocatable :: NpL(:^D&),NpR(:^D&)
  double precision :: xW,dFh

  integer :: numXI^D,dxI^D
  double precision :: xmaxL^D,xminL^D,xmaxR^D,xminR^D
  double precision :: xLI(:^D&,:),xRI(:^D&,:),QeLI(:^D&),QeRI(:^D&)
  double precision, allocatable :: maxhL(:,:),minhL(:,:),maxhR(:,:),minhR(:,:)

  double precision :: Ec,Eb,delta_l,delta_u,Ec_erg,Eb_erg,dE,mu0
  double precision :: tmax,tw1,tw2,Fleft,Fright
  integer :: numN,numE
  double precision, allocatable :: QeN(:),Ncol(:)
  double precision :: Nmin,Nmax,dNlog
  integer :: ccount=0
  double precision :: time_save=0
  integer :: count_save=0
  double precision :: t_update_Qe=-1.0d-7,dt_update_Qe=5.0d-4
  double precision, allocatable :: Ee(:,:),spectra(:,:),muN(:,:)

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    call set_coordinate_system("Cartesian_3D")

    unit_length        = 1.d9                                         ! cm
    unit_temperature   = 1.d6                                         ! K
    unit_numberdensity = 1.d9                                         ! cm^-3

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_source          => special_source
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_set_B0          => specialset_B0
    usr_set_J0          => specialset_J0
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output 
    usr_transform_w     => transform_w_usr
    usr_process_global  => special_global
    usr_special_convert => special_output

    call mhd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr()
    use mod_global_parameters

    heatunit=unit_pressure/unit_time          ! 3.697693390805347E-003 erg*cm^-3/s

    usr_grav=-2.74d4*unit_length/unit_velocity**2 ! solar gravity
    bQ0=1.d-4/heatunit ! background heating power density
    gzone=0.2d0 ! thickness of a ghostzone below the bottom boundary
    dya=(2.d0*gzone+xprobmax3-xprobmin3)/dble(jmax) ! cells size of high-resolution 1D solar atmosphere
    Busr=Busr/unit_magneticfield ! magnetic field strength at the bottom
    theta=30.d0*dpi/180.d0 ! the angle to the plane xy, 90-theta is the angle to the polarity inversion line of the arcade 
    kx=dpi/(xprobmax2-xprobmin2)
    ly=kx*dcos(theta)
    SRadius=69.61d0 ! Solar radius
    ! hydrostatic vertical stratification of density, temperature, pressure
    call inithdstatic

  end subroutine initglobaldata_usr

  subroutine inithdstatic
  !! initialize the table in a vertical line through the global domain
    use mod_global_parameters
    
    integer :: j,na,nb,ibc
    double precision, allocatable :: Ta(:),gg(:)
    double precision:: rpho,Ttop,Tpho,wtra,res,rhob,pb,htra,Ttr,Fc,invT,kappa

    integer :: n_val=49, i
    double precision :: h_val(49),t_val(49)
  
    rpho=0.71d15/unit_numberdensity ! number density at the bottom relaxla
    Tpho=1.d4/unit_temperature ! temperature of chromosphere
    Ttop=2.d6/unit_temperature ! estimated temperature in the top
    htra=0.2543d0 ! height of initial transition region
    wtra=0.01d0 ! width of initial transition region 
    !Ttr=1.d5/unit_temperature ! lowest temperature of upper profile
    Ttr=4.47d5/unit_temperature ! lowest temperature of upper profile
    Fc=2.d5/heatunit/unit_length ! constant thermal conduction flux
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

    !read rho T profiles from txt file
    if (iprob==2) then
      call read_init()
    endif


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
      na=floor((gzone-dx(3,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0))/dya+0.5d0)
      res=gzone-dx(3,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0)-(dble(na)-0.5d0)*dya
      rbc(ibc)=ra(na)+res/dya*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+res/dya*(pa(na+1)-pa(na))
    end do

    if (mype==0) then
     print*,'minra',minval(ra)
     print*,'rhob',rhob
     print*,'pb',pb
    endif

    call init_flare_heating()

  end subroutine inithdstatic

  subroutine init_flare_heating
    ! setting for flare heating
    use mod_global_parameters
    
    double precision :: maxh,xL^L,xR^L
    double precision :: Apt,temp,ks
    integer :: ix^D,refine_factor

    integer :: j,iBx,numBx,iBc,iBb,iNcol
    double precision :: Cllog,gammaN,beta,mu0,K,Nc,Nb,Efactor
    double precision :: dtx,tx,au,al,b
    double precision, allocatable :: Bu(:),Bl(:)

    character (20) :: fname

    ! heating parameter
    Fleft=0.5e11    ! energy flux at left footpoint
    Fright=0.6e11   ! energy flux at right footpoint
    tmax=60/unit_time   ! when flux reach a maximum value
    tw1=30/unit_time    ! time scale for flux increasing before tmax
    tw2=180/unit_time   ! time scale for flux decreasing after tmax


    ! initialize magnetic field lines
    xRmin2=2.0d0
    xRmax2=2.2d0
    xLmin2=-xRmax2
    xLmax2=-xRmin2
    refine_factor=2**(refine_max_level-1)
    do j=1,ndim
      dxL(j)=(xprobmax3-xprobmin3)/(domain_nx3*refine_factor)
      dxR(j)=dxL(j)
    enddo
    dFh=dxL(3)

    maxh=1.0  ! 5 Mm
    xW=0.4  ! width of loop in x direction
    numX(2)=floor((xRmax2-xRmin2)/dxL(2)+0.5)+1
    numX(3)=floor(maxh/dxL(3)+0.5)
    numX(1)=floor(xW/dxL(1)+0.5)+1
    numX1=numX(1)
    numX2=numX(2)
    numX3=numX(3)

    ixFmin1=1
    ixFmin2=1
    ixFmin3=1
    ixFmax1=numX1
    ixFmax2=numX2
    ixFmax3=numX3

    allocate(xFL(numX^D,ndim),xFR(numX^D,ndim))
    allocate(QeL(numX^D),QeR(numX^D),NpL(numX^D),NpR(numX^D))
    allocate(maxhL(numX3,2),minhL(numX3,2),maxhR(numX3,2),minhR(numX3,2))

    ks=tan(theta)
    do ix1=ixFmin1,ixFmax1
      do ix2=ixFmin2,ixFmax2
        xFL(ix1,ix2,1,2)=xLmin2+(ix2-1)*dxL(2)
        xFR(ix1,ix2,1,2)=xRmin2+(ix2-1)*dxR(2)
        xFL(ix1,ix2,1,1)=xFL(ix1,ix2,1,2)*ks-xW/2+ix1*dxL(1)
        xFR(ix1,ix2,1,1)=xFR(ix1,ix2,1,2)*ks-xW/2+ix1*dxR(1)
        xFL(ix1,ix2,1,3)=0.0
        xFR(ix1,ix2,1,3)=0.0
      enddo
    enddo

    QeL=0
    QeR=0


    ! electron's spectrum [electrons cm^-2 keV-1]
    delta_l=3.0 ! spectral index for E<Eb
    delta_u=4.0 ! spectral index for E>Eb
    Ec=20.0  ! cutoff energy [keV]
    Eb=100.0  ! break energy [keV]
    numE=280
    dE=1.0    ! [keV]

    ! parameters for electron deposition    
    Cllog=25.0  ! Coulomb logarithns
    gammaN=Cllog  ! full ionized plasma
    beta=2.0    ! full ionized plasma
    mu0=0.5   ! cos pitch angle
    K=2.0*dpi*const_e**4
    Ec_erg=Ec*1.0e3*const_ev
    Eb_erg=Eb*1.0e3*const_ev
    Nc=mu0*Ec_erg**2/((2.0+beta/2.0)*gammaN*K)
    Nb=mu0*Eb_erg**2/((2.0+beta/2.0)*gammaN*K)
    Efactor=K*gammaN*(delta_u-2.0)*(delta_l-2.0)/(2*mu0*Ec_erg**2)
    Efactor=Efactor/(((Eb/Ec)**(2.0-delta_l))*(delta_l-delta_u)+(delta_u-2))

    ! beta function    
    numBx=100001
    allocate(Bu(numBx),Bl(numBx))
    dtx=1.0/(numBx-1)
    Bu=0
    Bl=0
    au=delta_u/2.0
    al=delta_l/2.0
    b=1.0/3
    do iBx=2,numBx-1
      tx=dtx*(iBx-1.0)
      Bu(iBx)=Bu(iBx-1) + (tx**(au-1.0))*((1.0-tx)**(b-1.0))*dtx
      Bl(iBx)=Bl(iBx-1) + (tx**(al-1.0))*((1.0-tx)**(b-1.0))*dtx
    enddo
    Bu(numBx)=Bu(numBx-1)
    Bl(numBx)=Bl(numBx-1)

    numN=1001
    allocate(Ncol(numN),QeN(numN))
    Nmin=1.0e18
    Nmax=1.0e23
    dNlog=log10(Nmax/Nmin)/(numN-1.0)
    do iNcol=1,numN
      Ncol(iNcol)=Nmin*10**(dNlog*(iNcol-1.0))
      iBc=int((numBx-1)*Ncol(iNcol)/Nc+0.5)+1
      if (iBc>numBx)  then
        iBc=numBx
      endif
      iBb=int((numBx-1)*Ncol(iNcol)/Nb+0.5)+1
      if (iBb>numBx)  then
        iBb=numBx
      endif
      QeN(iNcol)=((Ncol(iNcol)/Nc)**(-delta_l/2.0))*(Bl(iBc)-Bl(iBb))
      QeN(iNcol)=QeN(iNcol)+((Ncol(iNcol)/Nc)**(-delta_u/2.0))*((Eb/Ec)**(delta_u-delta_l))*Bu(iBb)
    enddo
    QeN=QeN*Efactor
    deallocate(Bu,Bl)

    ! create table about electron distribution for different column depth 
    allocate(Ee(numN,numE),spectra(numN,numE),muN(numN,numE))
    call get_spectra()

  end subroutine init_flare_heating

  subroutine get_spectra()
    ! create table about electron distribution for different column depth
    use mod_global_parameters

    double precision :: K,gamma_,beta,dE_new
    double precision :: Fe0,Fe,qt,B0,const,const2,const3,keV_erg
    double precision :: Eph,delta,alpha,mec2,r0,E0,N0,c,temp
    double precision :: sigma0
    integer :: iNcol,iE
    character (20) :: fname

    mec2=511.0      ! energy of a static electron [keV]
    alpha=1.0/137.0 ! fine structure constant
    c=2.9979d10     ! light speed [cm/s]
    r0=2.8179d-13   ! radiu of electron [cm]
    keV_erg=1.0d3*const_ev  ! 1 keV = * erg
    sigma0=7.9e-25  ! [cm^2 keV]
    K=2*dpi*const_e**4
    gamma_=25.0
    beta=2.0

    ! fast electron spectra
    const=1.0*(delta_l-2.0)*(delta_u-2.0)/Ec_erg**2
    const=const/((delta_l-2.0)-((Eb/Ec)**(2.0-delta_l))*(delta_u-delta_l))
    const2=const*(Eb/Ec)**(delta_u-delta_l)
    do iE=1,numE
      Ee(1,iE)=Ec+(iE-1)*dE
      if (Ee(1,iE)<Eb) then
        spectra(1,iE)=const*(Ee(1,iE)/Ec)**(-delta_l)
      else
        spectra(1,iE)=const2*(Ee(1,iE)/Ec)**(-delta_u)
      endif
    enddo
    spectra=spectra*keV_erg   ! [electrons cm^-2 s^-1 keV^-1]

    ! change of electron's energy
    const=(2.0+beta/2.0)*gamma_*K/mu0
    const2=2.0/(4.0+beta)
    const3=beta/(4.0+beta)
    do iNcol=2,numN
      do iE=1,numE
        temp=const*Ncol(iNcol)/((Ee(1,iE)*keV_erg)**2)
        if (temp>0.99999) then
          temp=0.99999
        endif
        Ee(iNcol,iE)=Ee(1,iE)*((1-temp)**const2)
        muN(iNcol,iE)=mu0*((1-temp)**const3)
      enddo
    enddo

    ! renew spectra based on flux conservation
    do iNcol=2,numN
      do iE=1,numE
        if (iE<numE) then
          dE_new=Ee(iNcol,iE+1)-Ee(iNcol,iE)
        else
          dE_new=Ee(iNcol,iE)-Ee(iNcol,iE-1)
        endif
        spectra(iNcol,iE)=spectra(1,iE)*dE/dE_new
      enddo
    enddo

  end subroutine get_spectra

  subroutine read_init()
    !read initial rho and T from a txt file
    use mod_global_parameters

    integer :: ihr, numh, ixr, numx

    open(1, file='solar_init.txt')
    read(1,*)
    read(1,*) numx, numh
    allocate(xr(numx),hr(numh),rhor(numx,numh),Tr(numx,numh),pr(numx,numh))

    do ixr=1,numx
      read(1,*)
      read(1,*) xr(ixr)
      read(1,*)
      do ihr=1,numh
        read(1,*) hr(ihr),rhor(ixr,ihr),Tr(ixr,ihr)
        pr(ixr,ihr)=rhor(ixr,ihr)*Tr(ixr,ihr)
      end do
    end do
    close(1)
  end subroutine read_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: res,Bf(ixI^S,1:ndir),dhr,dxr
    integer :: ix^D,na,ihr,ixr
    logical, save :: first=.true.

    if(first)then
      if(mype==0) then
        write(*,*)'Simulating 3D solar atmosphere'
      endif
      first=.false.
    endif

    !if iprob==2, read initial rho and T from a txt file
    if (iprob==2) then
      dhr=hr(2)-hr(1)
      dxr=xr(2)-xr(1)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        na=floor((x(ix^D,3)-xprobmin3)/dhr+1.d0)
        res=x(ix^D,3)-xprobmin3-(dble(na)-1.d0)*dhr
        if (na < 2) then
          na=floor((x(ix^D,3)-xprobmin3+gzone)/dya+0.5d0)
          res=x(ix^D,3)-xprobmin3+gzone-(dble(na)-0.5d0)*dya
          w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
          w(ix^D,p_)  =pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
        else
          ixr=floor((x(ix^D,2)-xr(2))/dxr+1.0)
          w(ix^D,rho_)=rhor(ixr,na)+(res/dhr)*(rhor(ixr,na+1)-rhor(ixr,na))
          w(ix^D,p_)=pr(ixr,na)+(res/dhr)*(pr(ixr,na+1)-pr(ixr,na))
        endif
      {end do\}
    else
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        na=floor((x(ix^D,3)-xprobmin3+gzone)/dya+0.5d0)
        res=x(ix^D,3)-xprobmin3+gzone-(dble(na)-0.5d0)*dya
        w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
        w(ix^D,p_)  =pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
      {end do\}
    end if

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

  !> regenerate w output
  subroutine transform_w_usr(ixI^L,ixO^L,nw_in,w_in,x,w_out)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, nw_in
    double precision, intent(in)  :: w_in(ixI^S,1:nw_in)
    double precision, intent(in)  :: x(ixI^S, 1:ndim)
    double precision, intent(out) :: w_out(ixI^S,1:nw)

    integer :: ix^D
    double precision :: Atop, Alow, Apo
    double precision :: lQxmin,lQxmax

    Atop=Busr/kx*cos(kx*2.1d0)
    Alow=Busr/kx*cos(kx*2.2d0)
    lQxmin=-4.0d0
    lQxmax=4.0d0

    w_out(ixI^S,rho_) = w_in(ixI^S,rho_)
    w_out(ixI^S,e_) = w_in(ixI^S,e_)
    w_out(ixI^S,mom(:)) = w_in(ixI^S,mom(:))
    w_out(ixI^S,mag(:)) = w_in(ixI^S,mag(:))

    if (x(ix^D,2)<=0.d0) then
      w_out(ix^D,tracer(1))=10.d0
    endif
    if (x(ix^D,2)>=0.d0) then
      w_out(ix^D,tracer(2))=10.d0
    endif


  end subroutine transform_w_usr

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
    case(5)
      ! fixed zero velocity
      do idir=1,ndir
        w(ixO^S,mom(idir)) =-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(idir))&
                   /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)
      end do
      ! fixed b1 b2 b3
      if(B0field) then
        w(ixO^S,mag(:))=0.d0
      else
        call specialset_B0(ixI^L,ixO^L,x,Bf)
        w(ixO^S,mag(1:ndir))=Bf(ixO^S,1:ndir)
      endif
      ! fixed gravity stratification of density and pressure pre-determined in initial condition
      do ix3=ixOmin3,ixOmax3
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,rho_)=rbc(ix3)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,p_)=pbc(ix3)
      enddo
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(6)
      ixInt^L=ixO^L;
      ixIntmin3=ixOmin3-1;ixIntmax3=ixOmin3-1;
      call mhd_get_pthermal(w,x,ixI^L,ixInt^L,pth)
      ixIntmin3=ixOmin3-1;ixIntmax3=ixOmax3;
      call getggrav(ggrid,ixI^L,ixInt^L,x)
      ! fill pth, rho ghost layers according to gravity stratification
      invT(ixOmin3-1^%3ixO^S)=w(ixOmin3-1^%3ixO^S,rho_)/pth(ixOmin3-1^%3ixO^S)
      tmp=0.d0
      do ix3=ixOmin3,ixOmax3
        tmp(ixOmin3-1^%3ixO^S)=tmp(ixOmin3-1^%3ixO^S)+0.5d0*&
            (ggrid(ix3^%3ixO^S)+ggrid(ix3-1^%3ixO^S))*invT(ixOmin3-1^%3ixO^S)
        w(ix3^%3ixO^S,p_)=pth(ixOmin3-1^%3ixO^S)*dexp(tmp(ixOmin3-1^%3ixO^S)*dxlevel(3))
        w(ix3^%3ixO^S,rho_)=w(ix3^%3ixO^S,p_)*invT(ixOmin3-1^%3ixO^S)
      enddo
      ! fixed zero velocity
      do idir=1,ndir
        w(ixO^S,mom(idir)) =-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(idir))&
                     /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,rho_)
      end do
      ! zero normal gradient extrapolation
      do ix3=ixOmin3,ixOmax3
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(:))=(1.0d0/3.0d0)* &
                    (-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-2,mag(:))&
               +4.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-1,mag(:)))
      enddo
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select
    
  end subroutine specialbound_usr

  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    double precision                :: ggrid(ixI^S)

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,3)=ggrid(ixO^S)

  end subroutine gravity

  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*(SRadius/(SRadius+x(ixO^S,3)))**2
  end subroutine

  subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: lQgrid(ixI^S),bQgrid(ixI^S)

    integer :: ixO^D

    ! add global background heating bQ
    call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)

    !! add localized external heating lQ concentrated at feet of loops
    if (iprob > 2)then
      call getlQ(lQgrid,ixI^L,ixO^L,qtC,wCT,x)
      w(ixO^S,e_)=w(ixO^S,e_)+qdt*lQgrid(ixO^S)
    endif

  end subroutine special_source

  subroutine getbQ(bQgrid,ixI^L,ixO^L,qt,w,x)
  ! calculate background heating bQ
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)

    double precision :: bQgrid(ixI^S)

    bQgrid(ixO^S)=bQ0*dexp(-x(ixO^S,3)/5.d0)

  end subroutine getbQ

  subroutine update_block_para(igrid,dxb^D,xb^L)
    ! update parameters of the block
    use mod_global_parameters

    integer          :: ixO^L,ixO^D,j
    integer :: igrid
    double precision :: dxb^D,xb^L

    ^D&ixOmin^D=ixmlo^D\
    ^D&ixOmax^D=ixmhi^D\

    ^D&dxb^D=rnode(rpdx^D_,igrid)\

    ^D&xbmin^D=ps(igrid)%x(ixOmin^DD,^D)-dxb^D/2.0\
    ^D&xbmax^D=ps(igrid)%x(ixOmax^DD,^D)+dxb^D/2.0\

  end subroutine update_block_para

  subroutine special_global(iit,qt)
    !! run at the begining of a step
    use mod_global_parameters

    integer, intent(in) :: iit
    double precision, intent(in) :: qt
    character (20) :: fname
    integer :: ix^D
    double precision :: t1,t2,t3

    if (convert) then
      call update_Bfield()
      call get_flare_eflux()
      call update_region()
    endif

    if (iprob>=6) then
      if (qt>t_update_Qe+10*dt_update_Qe) then
        t_update_Qe=qt-1.0e-7
      endif

      if (qt>t_update_Qe) then
        call update_Bfield()
        call get_flare_eflux()
        call update_region()
        t_update_Qe=t_update_Qe+dt_update_Qe
      endif
    endif

  end subroutine special_global

  subroutine update_region()
    ! update parameters that are helpful to find heating region
    use mod_global_parameters

    integer :: ix^D
    double precision :: maxx1,maxx2,minx1,minx2

    !left foot
    xmaxL1=xFL(1,1,1,1)
    xminL1=xFL(1,1,1,1)
    xmaxL2=xFL(1,1,1,2)
    xminL2=xFL(1,1,1,2)

    do ix3=1,numX3
      maxx1=xFL(1,1,ix3,1)
      maxx2=xFL(1,1,ix3,2)
      minx1=xFL(1,1,ix3,1)
      minx2=xFL(1,1,ix3,2)

      do ix2=1,numX2
        do ix1=1,numX1          
          if (xFL(ix^D,1)>maxx1) maxx1=xFL(ix^D,1)
          if (xFL(ix^D,1)<minx1) minx1=xFL(ix^D,1)
          if (xFL(ix^D,2)>maxx2) maxx2=xFL(ix^D,2)
          if (xFL(ix^D,2)<minx2) minx2=xFL(ix^D,2)
        enddo
      enddo

      maxhL(ix3,1)=maxx1
      maxhL(ix3,2)=maxx2
      minhL(ix3,1)=minx1
      minhL(ix3,2)=minx2

      if (maxx1>xmaxL1) xmaxL1=maxx1
      if (maxx2>xmaxL2) xmaxL2=maxx2
      if (minx1<xminL1) xminL1=minx1
      if (minx2<xminL2) xminL2=minx2
    enddo    


    !right foot
    xmaxR1=xFR(1,1,1,1)
    xminR1=xFR(1,1,1,1)
    xmaxR2=xFR(1,1,1,2)
    xminR2=xFR(1,1,1,2)

    do ix3=1,numX3
      maxx1=xFR(1,1,ix3,1)
      maxx2=xFR(1,1,ix3,2)
      minx1=xFR(1,1,ix3,1)
      minx2=xFR(1,1,ix3,2)

      do ix2=1,numX2
        do ix1=1,numX1          
          if (xFR(ix^D,1)>maxx1) maxx1=xFR(ix^D,1)
          if (xFR(ix^D,1)<minx1) minx1=xFR(ix^D,1)
          if (xFR(ix^D,2)>maxx2) maxx2=xFR(ix^D,2)
          if (xFR(ix^D,2)<minx2) minx2=xFR(ix^D,2)
        enddo
      enddo

      maxhR(ix3,1)=maxx1
      maxhR(ix3,2)=maxx2
      minhR(ix3,1)=minx1
      minhR(ix3,2)=minx2

      if (maxx1>xmaxR1) xmaxR1=maxx1
      if (maxx2>xmaxR2) xmaxR2=maxx2
      if (minx1<xminR1) xminR1=minx1
      if (minx2<xminR2) xminR2=minx2
    enddo    

    xminL3=xFL(1,1,1,3)
    xminR3=xFR(1,1,1,3)
    xmaxL3=xFL(1,1,numX3,3)
    xmaxR3=xFR(1,1,numX3,3)

  end subroutine update_region

  subroutine get_flare_eflux()
    ! calculate electron deposit energy flux
    use mod_global_parameters

    integer          :: ixO^L,ixO^D,j

    integer :: ix^D
    double precision :: F0LT(numX1,numX2),F0RT(numX1,numX2)
    double precision :: ixHalf,numQuar,qt

    ! energy flux in the left/right footpoints [erg cm^-2 s^-1]
    F0LT=Fleft
    F0RT=Fright

    qt=global_time
    if (qt<=tmax)then
      F0LT=F0LT*exp(-(qt-tmax)**2/tw1**2)
      F0RT=F0RT*exp(-(qt-tmax)**2/tw1**2)
    else
      F0LT=F0LT*exp(-(qt-tmax)**2/tw2**2)
      F0RT=F0RT*exp(-(qt-tmax)**2/tw2**2)
    endif

    if (mype==0) then
      print *, global_time,exp(-(qt-tmax)**2/tw1**2),F0LT(1,1),F0RT(1,1)
    endif

    ! update distribution of the energy flux
    call update_e_distr(xFL,QeL,1)
    call update_e_distr(xFR,QeR,2)

    do ix1=1,numX1
      do ix2=1,numX2
        QeL(ix1,ix2,:)=F0LT(ix1,ix2)*QeL(ix1,ix2,:)/heatunit 
        QeR(ix1,ix2,:)=F0RT(ix1,ix2)*QeR(ix1,ix2,:)/heatunit 
      enddo
    enddo

  end subroutine get_flare_eflux

  subroutine update_e_distr(xLR,QeLR,iLR)
    ! calculate electron deposit energy flux
    use mod_global_parameters

    integer          :: ixO^L,ixO^D,j

    double precision :: Np(numX^D),Nh(numX^D)
    double precision :: xLR(numX^D,ndim)
    double precision :: QeLR(numX^D)
    double precision :: dl,dl^D,rx,ry
    integer :: ix^D
    double precision :: sumQe
    integer :: iLR,iNlog


    if (iLR==1) then
      Np=NpL
    else
      Np=NpR
    endif


    ! column depth along field lines
    Nh=0
    do ix1=1,numX1
      do ix2=1,numX2
        do ix3=numX3-1,1,-1
          ! distance between two points
          dl1=XLR(ix1,ix2,ix3+1,1)-XLR(ix1,ix2,ix3,1)
          dl2=XLR(ix1,ix2,ix3+1,2)-XLR(ix1,ix2,ix3,2)
          dl3=XLR(ix1,ix2,ix3+1,3)-XLR(ix1,ix2,ix3,3)
          dl=sqrt(dl1**2+dl2**2+dl3**2)*unit_length

          ! column depth
          Nh(ix^D)=Nh(ix1,ix2,ix3+1)+ &
                   dl*(Np(ix1,ix2,ix3)+Np(ix1,ix2,ix3+1))/2.0

          ! expansion of loop
          if (ix1<numX1) then
            rx=abs((XLR(ix1+1,ix2,ix3,1)-XLR(ix1,ix2,ix3,1))/ &
                   (XLR(ix1+1,ix2,1,1)-XLR(ix1,ix2,1,1)))
          else
            rx=abs((XLR(ix1,ix2,ix3,1)-XLR(ix1-1,ix2,ix3,1))/ &
                   (XLR(ix1,ix2,1,1)-XLR(ix1-1,ix2,1,1)))
          endif
          !
          if (ix2<numX2) then
            ry=abs((XLR(ix1,ix2+1,ix3,1)-XLR(ix1,ix2,ix3,1))/ &
                   (XLR(ix1,ix2+1,1,1)-XLR(ix1,ix2,1,1)))
          else
            ry=abs((XLR(ix1,ix2,ix3,1)-XLR(ix1,ix2-1,ix3,1))/ &
                   (XLR(ix1,ix2,1,1)-XLR(ix1,ix2-1,1,1)))
          endif

          iNlog=int(log10(Nh(ix^D)/Nmin)/dNlog+0.5)+1
          if (iNlog<1) then
            iNlog=1
          endif
          QeLR(ix^D)=QeN(iNlog)*Np(ix^D)/(rx*ry)
        enddo
      enddo
    enddo

  end subroutine update_e_distr

  subroutine update_Bfield()
    ! update the location of mangetic field lines
    use mod_usr_methods
    use mod_global_parameters

    integer :: ix^D,flag
    logical :: forward
    integer :: status(mpi_status_size)
    integer :: ids

    !forward=.FALSE.
    !call update_line(xFL(2,5,:,:),NpL(2,5,:),numX3,forward)

    !print *, mype,mod(mype+npe-1,npe),mod(mype+1,npe)


    forward=.FALSE.
    do ix1=ixFmin1,ixFmax1
      do ix2=ixFmin2,ixFmax2
        !call MPI_SEND(mype,1,MPI_INTEGER,mod(mype+1,npe),ix1+ix2,icomm,ierrmpi)
        !call MPI_RECV(flag,1,MPI_INTEGER,mod(mype+npe-1,npe),ix1+ix2,icomm,status,ierrmpi)
        call update_line(xFL(ix1,ix2,:,:),NpL(ix1,ix2,:),numX3,forward)
      enddo
    enddo

    forward=.TRUE.
    do ix1=ixFmin1,ixFmax1
      do ix2=ixFmin2,ixFmax2
        call update_line(xFR(ix1,ix2,:,:),NpR(ix1,ix2,:),numX3,forward)
      enddo
    enddo

  end subroutine update_Bfield

  subroutine update_line(xf,Npf,numP,forward)
    ! trace a field line
    use mod_usr_methods
    use mod_global_parameters

    integer :: numP
    double precision :: xf(numP,ndim),Npf(numP)
    double precision :: xtemp(ndim)
    logical :: forward

    double precision :: dxb^D,xb^L
    integer :: inblock,flag,flag_sum

    integer :: ixO^L,ixO^D,j
    integer :: iigrid,igrid
    integer :: igrid_prev,igrid_now,igrid_next
    integer :: ipe_prev,ipe_now,ipe_next
    integer :: igrid_rc,ipe_rc
    integer :: mainpe
    integer :: status(mpi_status_size)

    logical :: newpe

    ! set processor 0 as main processor
    mainpe=0


    ! find the grid and pe of the first point
    LOOP1: do iigrid=1,igridstail; igrid=igrids(iigrid);
      call update_block_para(igrid,dxb^D,xb^L)
      inblock=0
      {if (xf(1,^DB)>=xbmin^DB .and. xf(1,^DB)<xbmax^DB) inblock=inblock+1\}
      if (inblock==ndim) then
        igrid_now=igrid
        ipe_now=mype
        call MPI_SEND(igrid_now,1,MPI_INTEGER,mainpe,0,icomm,ierrmpi)
        call MPI_SEND(ipe_now,1,MPI_INTEGER,mainpe,1,icomm,ierrmpi)
        exit LOOP1
      endif
    enddo LOOP1

    if (mype==mainpe) then
      call MPI_RECV(igrid_rc,1,MPI_INTEGER,MPI_ANY_SOURCE,0,icomm,status,ierrmpi)
      call MPI_RECV(ipe_rc,1,MPI_INTEGER,MPI_ANY_SOURCE,1,icomm,status,ierrmpi)
      igrid_now=igrid_rc
      ipe_now=ipe_rc    
    endif
    call MPI_BCAST(igrid_now,1,MPI_INTEGER,mainpe,icomm,ierrmpi)
    call MPI_BCAST(ipe_now,1,MPI_INTEGER,mainpe,icomm,ierrmpi)


    ! other points in field line    
    j=1
    do while (j<numP)
      newpe=.TRUE.
      mainpe=ipe_now

      if (mype==ipe_now) then
        igrid=igrid_now
        ! find next point in this field line and get density of current point
        call find_next(igrid,xf(j,:),xf(j+1,:),Npf(j),forward)
        ! check the pe of next point
        call check_next(igrid,igrid_next,ipe_next,xf(j+1,:),newpe)
      endif

      call MPI_BCAST(newpe,1,MPI_LOGICAL,mainpe,icomm,ierrmpi)
      call MPI_BCAST(xf(j+1,:),ndim,MPI_DOUBLE_PRECISION,mainpe,icomm,ierrmpi)
      call MPI_BCAST(Npf(j),1,MPI_DOUBLE_PRECISION,mainpe,icomm,ierrmpi)
      
      ! next point is in another pe, find out grid and pe numbers
      if (newpe .eqv. .TRUE.) then
        call find_grid(ipe_now,igrid_now,ipe_next,igrid_next,xf(j+1,:),j)
      endif

      ! prepare for next point 
      if (newpe .eqv. .TRUE.) then
        ipe_now=ipe_next
        igrid_now=igrid_next
      else           
        if (mype==ipe_now) then
          igrid_now=igrid_next
        endif
      endif
      j=j+1
    enddo
 

    ! find the density for the last point
    if (j==numP) then
      mainpe=ipe_now
      if (mype==ipe_now) then
        call find_next(igrid_now,xf(j,:),xtemp,Npf(j),forward)
      endif
      call MPI_BCAST(Npf(j),1,MPI_DOUBLE_PRECISION,mainpe,icomm,ierrmpi)
    endif

  end subroutine update_line

  subroutine find_next(igrid,xf0,xf1,Np0,forward)
    !find next point
    use mod_usr_methods
    use mod_global_parameters

    integer :: igrid
    logical :: forward

    integer          :: ixO^L,ixO^D,j
    double precision :: xf0(ndim),xf1(ndim),Np0
    double precision :: dxf(ndim)
    double precision :: dxb^D,xb^L,xd^D
    integer          :: ixb^D,ix^D,ixbl^D
    double precision :: Bnear(0:1^D&,ndim),Nnear(0:1^D&)
    double precision :: Bx(ndim),factor(0:1^D&)

    ^D&ixOmin^D=ixmlo^D\
    ^D&ixOmax^D=ixmhi^D\

    call update_block_para(igrid,dxb^D,xb^L)

    ^D&ixbl^D=floor((xf0(^D)-ps(igrid)%x(ixOmin^DD,^D))/dxb^D)+ixOmin^D\
    ^D&xd^D=(xf0(^D)-ps(igrid)%x(ixbl^DD,^D))/dxb^D\


    ! interpolation
    if (B0field) then
      {do ix^DB=0,1\}
        do j=1,ndim
          Bnear(ix^D,j)=ps(igrid)%w(ixbl^D+ix^D,mag(j))+&
                        ps(igrid)%B0(ixbl^D+ix^D,j,0)
        enddo
        Nnear(ix^D)=ps(igrid)%w(ixbl^D+ix^D,rho_)
      {enddo\}
    else
      {do ix^DB=0,1\}
        do j=1,ndim
          Bnear(ix^D,j)=ps(igrid)%w(ixbl^D+ix^D,mag(j))
        enddo
        Nnear(ix^D)=ps(igrid)%w(ixbl^D+ix^D,rho_)
      {enddo\}
    endif
    
    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
    {enddo\}      

    Bx=0
    Np0=0
    {do ix^DB=0,1\}
      {Bx(^DB)=Bx(^DB)+Bnear(ix^DD,^DB)*factor(ix^DD)\}
      Np0=Np0+Nnear(ix^D)*factor(ix^D)
    {enddo\}      
    Np0=Np0*unit_numberdensity

    ! find next point based on magnetic field direction
    if (forward .eqv. .TRUE.) then
      do j=1,ndim
        dxf(j)=dFh*Bx(j)/abs(Bx(ndim))
      enddo
    else
      do j=1,ndim
        dxf(j)=-dFh*Bx(j)/abs(Bx(ndim))
      enddo
    endif

    do j=1,ndim
      xf1(j)=xf0(j)+dxf(j)
    enddo

  end subroutine find_next

  subroutine check_next(igrid,igrid_next,ipe_next,xf1,newpe)
    ! check the grid and pe of next point
    use mod_usr_methods
    use mod_global_parameters
    use mod_forest

    integer :: igrid,igrid_next,ipe_next
    double precision :: xf1(ndim)
    double precision :: dxb^D,xb^L

    integer :: inblock,inblock_n,inblock_nc
    integer :: ix^D
    logical :: newpe

    integer :: igrid_nb,ipe_nb

    call update_block_para(igrid,dxb^D,xb^L)
    inblock=0
    {if (xf1(^D)>=xbmin^D .and. xf1(^D)<xbmax^D) inblock=inblock+1\}
    if (inblock==ndim) then
      ! in the same grid with previous point
      igrid_next=igrid
      ipe_next=mype
      newpe=.FALSE.
    else
      ! not in the same grid
      {do ix^D=-1,1,1\}
        ! check neighbor
        igrid_nb=neighbor(1,ix^D,igrid)
        ipe_nb=neighbor(2,ix^D,igrid)
        if (mype==ipe_nb .and. igrid_inuse(igrid_nb,ipe_nb)) then
          call update_block_para(igrid_nb,dxb^D,xb^L)
          inblock_n=0
          {if (xf1(^D)>=xbmin^D .and. xf1(^D)<xbmax^D) inblock_n=inblock_n+1\}
          if (inblock_n==ndim) then
            ! in neighbor
            igrid_next=igrid_nb
            ipe_next=mype
            newpe=.FALSE.
          endif
        endif

        ! check neighbor_child
        igrid_nb=neighbor_child(1,ix^D,igrid)
        ipe_nb=neighbor_child(2,ix^D,igrid)
        if (mype==ipe_nb .and. igrid_inuse(igrid_nb,ipe_nb)) then
          call update_block_para(igrid_nb,dxb^D,xb^L)
          inblock_nc=0
          {if (xf1(^D)>=xbmin^D .and. xf1(^D)<xbmax^D) inblock_nc=inblock_nc+1\}
          if (inblock_nc==ndim) then
            ! in neighbor child
            igrid_next=igrid_nb
            ipe_next=mype
            newpe=.FALSE.
          endif
        endif
      {enddo\}
    endif    

  end subroutine check_next

  subroutine find_grid(ipe_now,igrid_now,ipe_next,igrid_next,xf1,nj)
    ! find for grid and pe numbers
    use mod_usr_methods
    use mod_global_parameters
    use mod_forest

    integer :: igrid,iigrid,igrid_now,ipe_now,igrid_next,ipe_next
    integer :: igrid_recv,ipe_recv
    double precision :: xf1(ndim)

    double precision :: dxb^D,xb^L
    integer :: inblock
    integer :: ix^D,i,j,nj
    integer :: found,found_recv,mainpe
    integer :: status(mpi_status_size)
    integer :: grid_nb(2*(3**ndim)),pe_nb(2*(3**ndim))
    integer :: numblock,flag

    double precision :: igrid_nb,ipe_nb

    mainpe=0

    numblock=2*(3**ndim)

    ! record neighbor blocks
    if (mype==ipe_now) then
      j=1
      {do ix^D=-1,1,1\}
        grid_nb(j)=neighbor(1,ix^D,igrid_now)
        pe_nb(j)=neighbor(2,ix^D,igrid_now)
        grid_nb(j+numblock/2)=neighbor_child(1,ix^D,igrid_now)
        pe_nb(j+numblock/2)=neighbor_child(2,ix^D,igrid_now)
        j=j+1
      {enddo\}
    endif    

    call MPI_BCAST(grid_nb,numblock,MPI_INTEGER,ipe_now,icomm,ierrmpi)
    call MPI_BCAST(pe_nb,numblock,MPI_INTEGER,ipe_now,icomm,ierrmpi)

    ! check neighbors
    LOOPNB: do j=1,numblock
      if (mype==pe_nb(j) .and. igrid_inuse(grid_nb(j),pe_nb(j))) then
        call update_block_para(grid_nb(j),dxb^D,xb^L)
        inblock=0
        {if (xf1(^D)>=xbmin^D .and. xf1(^D)<xbmax^D) inblock=inblock+1\}
        if (inblock==ndim) then
          igrid_next=grid_nb(j)
          ipe_next=mype
          found=1
          call MPI_SEND(j,1,MPI_INTEGER,ipe_now,nj,icomm,ierrmpi)
          exit LOOPNB
        endif
      endif
    enddo LOOPNB

    if (mype==ipe_now) then
      call MPI_RECV(i,1,MPI_INTEGER,MPI_ANY_SOURCE,nj,icomm,status,ierrmpi)
      j=i
    endif
    call MPI_BCAST(j,1,MPI_INTEGER,ipe_now,icomm,ierrmpi)
    igrid_next=grid_nb(j)
    ipe_next=pe_nb(j)


    if (mype==ipe_next) then
      call update_block_para(igrid_next,dxb^D,xb^L)
      if (xf1(1)<xbmin1 .or. xf1(1)>xbmax1 .or. &
          xf1(2)<xbmin2 .or. xf1(2)>xbmax2 .or. &
          xf1(3)<xbmin3 .or. xf1(3)>xbmax3) then
        if (found>0) then
          print *, 'wrong block is in neighbors in different pe'
        else
          print *, 'wrong block is not in neighbors in different pe'
        endif
        print *, xf1
        print *, xb^L
        print *, ipe_next,igrid_next,nj
        print *, igrid_inuse(igrid_next,ipe_next)
      endif
    endif

  end subroutine find_grid

  subroutine getlQ(lQgrid,ixI^L,ixO^L,qt,w,x)
    !!calculate localized heating at chromosphere
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    !double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in) :: w(ixI^S,1:nw)

    double precision :: lQgrid(ixI^S),asym,lQt,lQd,lQst,hp,flash,Apo,Atop,Alow
    double precision :: lQr,ks
    double precision :: sscale,tonset,tstop,tflare,tpeak,tscale,fstop,thstop
    double precision :: lQheat, lQt1, lQt2, lQtw
    integer          :: ixO^D
    integer :: ix,iy,iz

    integer :: ix^D
    integer :: flagL,flagR

    !-----------------------------------------------------------------------------
    ! 180 s of the flare heating
    lQt=5.d1/unit_time
    lQst=18.d1/unit_time
    lQt1=3.d1/unit_time
    lQt2=15.d1/unit_time
    lQd=0.05d0
    lQtw=60.d0/unit_time
    lQr=0.05
    asym=0.8d0
    hp=0.125d0
    lQgrid=0.d0
    Atop=Busr/kx*cos(kx*2.1d0)
    Alow=Busr/kx*cos(kx*2.2d0)
    lQheat=120.0d0/heatunit

    ks=tan(theta)

    select case(iprob)

    case(5)
    {do ixO^DB=ixOmin^DB,ixOmax^DB\}
      Apo=Busr/kx*cos(kx*x(ixO^D,2))*dexp(-ly*x(ixO^D,3))
      if(Apo>Alow .and. Apo<Atop .and. x(ixO^D,2)>zero) then
        lQgrid(ixO^D)=(lQheat/sqrt(3.14d0*lQd))*dexp(-(x(ixO^D,3)-hp)**2/lQd)*(1/(1+asym))
      endif

      if(Apo>Alow .and. Apo<Atop .and. x(ixO^D,2)<zero) then
        lQgrid(ixO^D)=(lQheat/sqrt(3.14d0*lQd))*dexp(-(x(ixO^D,3)-hp)**2/lQd)*(asym/(1+asym))
      endif
    {enddo\}
    lQgrid(ixO^S)=lQgrid(ixO^S)*dsqrt(1/(3.14*lQtw**2))*dexp(-(qt-2.0*lQtw)**2/(lQtw**2))

    case(6)
    !{do ixO^DB=ixOmin^DB,ixOmax^DB\}
    !  flagL=0
    !  flagR=0
    !  {if (x(ixO^DD,^D)>=xminL^D .and. x(ixO^DD,^D)<=xmaxL^D) flagL=flagL+1 \}
    !  {if (x(ixO^DD,^D)>=xminR^D .and. x(ixO^DD,^D)<=xmaxR^D) flagR=flagR+1 \}

    !  ! left foot
    !  if (flagL==ndim) then
    !    lQgrid(ixO^D)=0
    !    ix3=floor((x(ixO^D,3)-xFL(1^D&,3))/dFh)+1
    !    call interp_lQ(x(ixO^D,:),lQgrid(ixO^D),xFL(:,:,ix3:ix3+1,:),QeL(:,:,ix3:ix3+1))
    !  endif

    !  ! right foot
    !  if (flagR==ndim) then
    !    lQgrid(ixO^D)=0
    !    ix3=floor((x(ixO^D,3)-xFR(1^D&,3))/dFh)+1
    !    call interp_lQ(x(ixO^D,:),lQgrid(ixO^D),xFR(:,:,ix3:ix3+1,:),QeR(:,:,ix3:ix3+1))
    !  endif
    !{enddo\}

    end select



  end subroutine getlQ

  subroutine interp_lQ(xc,lQc,xLR,QeLR)
    ! get the heating of a cell via interpolation
    use mod_global_parameters

    double precision :: xc(ndim),lQc
    double precision :: xLR(numX1,numX2,2,ndim),QeLR(numX1,numX2,2)
    double precision :: xi(numX1,numX2,2),Qi(numX1,numX2)
    integer :: ix^D
    double precision :: xd3,x1n(2,2),x2n(2,2),xdn(2,2),Qn(2,2),dn
    double precision :: dyl,dyr,dxm
    double precision :: xil(2),xir(2),Ql,Qr

    ! interplot for z 
    xd3=(xc(3)-xLR(1,1,1,3))/dFh
    do ix1=1,numX1
      do ix2=1,numX2
        xi(ix1,ix2,1)=XLR(ix1,ix2,1,1)*(1-xd3)+XLR(ix1,ix2,2,1)*xd3
        xi(ix1,ix2,2)=XLR(ix1,ix2,1,2)*(1-xd3)+XLR(ix1,ix2,2,2)*xd3
        Qi(ix1,ix2)=QeLR(ix1,ix2,1)*(1-xd3)+QeLR(ix1,ix2,2)*xd3
      enddo
    enddo


    ! looking for nearby points to do interpolation
    x1n(:,:)=xi(1,1,1)
    x2n(:,:)=xi(1,1,2)
    xdn=8.0*dFh

    do ix1=1,numX1-1
      do ix2=1,numX2-1
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
    if (sum(xdn)==0) then
      lQc=Qn(1,1)
    elseif (sum(xdn)>0 .and. sum(xdn)<8.0*dFh) then
      dyl=(xc(2)-x2n(1,1))/(x2n(1,2)-x2n(1,1))
      dyr=(xc(2)-x2n(2,1))/(x2n(2,2)-x2n(2,1))
      xil(1)=x1n(1,1)*(1.0-dyl)+x1n(1,2)*dyl
      xil(2)=x2n(1,1)*(1.0-dyl)+x2n(1,2)*dyl
      Ql=Qn(1,1)*(1.0-dyl)+Qn(1,2)*dyl
      xir(1)=x1n(2,1)*(1.0-dyr)+x1n(2,2)*dyr
      xir(2)=x2n(2,1)*(1.0-dyr)+x2n(2,2)*dyr
      Qr=Qn(2,1)*(1.0-dyr)+Qn(2,2)*dyr

      dxm=(xc(1)-xil(1))/(xir(1)-xil(1))
      lQc=Ql*(1.0-dxm)+Qr*dxm
    endif

  end subroutine interp_lQ

  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
  ! Enforce additional refinement or coarsening
  ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen
    double precision :: Atop, Alow, ks, xW_max, l_scale

    Atop=Busr/kx*cos(kx*2.0d0)
    Alow=Busr/kx*cos(kx*2.4d0)
    ks=tan(theta)
    l_scale=1.5
    xW_max=0.5


    if (iprob<=2) then    
      ! fix the bottom layer to the highest level
      if (any(x(ixO^S,3)<=0.3d0)) then
        if (level<3) then
          refine=1
          coarsen=-1
        else
          refine=-1
          coarsen=-1
        endif
      endif
    endif


    if (iprob==5) then
      if (any(Busr/kx*cos(kx*x(ixO^S,2))*dexp(-ly*x(ixO^S,3))>Alow .and. &
              Busr/kx*cos(kx*x(ixO^S,2))*dexp(-ly*x(ixO^S,3))<Atop)) then
        refine=0
        coarsen=0
      else
        if (any(x(ixO^S,3)<=0.3d0)) then
          if (level<3) then
            refine=1
            coarsen=-1
          else
            refine=-1
            coarsen=-1
          endif
        endif
      endif
    endif


    if (iprob>5) then
      if (any(Busr/kx*cos(kx*x(ixO^S,2))*dexp(-ly*x(ixO^S,3))>Alow .and. &
              Busr/kx*cos(kx*x(ixO^S,2))*dexp(-ly*x(ixO^S,3))<Atop .and. &
              abs(x(ixO^S,1)-x(ixO^S,2)*ks)<=xW/2)) then
      !if (any(x(ixO^S,1)>x(ixO^S,2)*ks-xw/2 .and. &
      !        x(ixO^S,1)<x(ixO^S,2)*ks+xw/2)) then
          !near or inner loop
        refine=0
        coarsen=0
      !lower atmosphere
      else if (any(x(ixO^S,3)<=0.3d0)) then
        if (level<3) then
          refine=1
          coarsen=-1
        else
          refine=-1
          coarsen=-1
        endif
      else
        refine=-1
        coarsen=0
      endif
    endif

  end subroutine special_refine_grid

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_global_parameters
    use mod_radiative_cooling

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixI^S),B2(ixI^S),tmp2(ixI^S),dRdT(ixI^S)
    double precision :: ens(ixI^S),divb(ixI^S), lQgrid(ixI^S)
    double precision :: Btotal(ixI^S,1:ndir),curlvec(ixI^S,1:ndir)
    integer :: idirmin,idir,ix^D

    double precision :: t

    double precision :: I0,hMu,hMu_max,hMu_min,dhMu,kb,kev
    double precision :: kbT(ixI^S),gff(ixI^S),I_SXR(ixI^S),EM(ixI^S)
    integer          :: ixO^D,ihMu,num_hMu

    double precision :: aia_T(101),aia_131(101)
    double precision :: logT,logGT,LOS
    integer          :: iTl

    double precision :: dS(ixI^S)
    double precision :: e_k,e_b,e_t

    kb=1.38d-23
    !I0=1.07d-18
    I0=7.93d-24
    LOS=1.0d8
    keV=1.602d-16
    hMu_max=10
    hMu_min=4
    dhMu=0.1
    num_hMu=int((hMu_max-hMu_min)/dhMu)

   !aia T list
    data aia_T /4.        ,  4.05000019,  4.0999999 ,  4.1500001 ,  4.19999981, &
                4.25      ,  4.30000019,  4.3499999 ,  4.4000001 ,  4.44999981, &
                4.5       ,  4.55000019,  4.5999999 ,  4.6500001 ,  4.69999981, &
                4.75      ,  4.80000019,  4.8499999 ,  4.9000001 ,  4.94999981, &
                5.        ,  5.05000019,  5.0999999 ,  5.1500001 ,  5.19999981, &
                5.25      ,  5.30000019,  5.3499999 ,  5.4000001 ,  5.44999981, &
                5.5       ,  5.55000019,  5.5999999 ,  5.6500001 ,  5.69999981, &
                5.75      ,  5.80000019,  5.8499999 ,  5.9000001 ,  5.94999981, &
                6.        ,  6.05000019,  6.0999999 ,  6.1500001 ,  6.19999981, &
                6.25      ,  6.30000019,  6.3499999 ,  6.4000001 ,  6.44999981, &
                6.5       ,  6.55000019,  6.5999999 ,  6.6500001 ,  6.69999981, &
                6.75      ,  6.80000019,  6.8499999 ,  6.9000001 ,  6.94999981, &
                7.        ,  7.05000019,  7.0999999 ,  7.1500001 ,  7.19999981, &
                7.25      ,  7.30000019,  7.3499999 ,  7.4000001 ,  7.44999981, &
                7.5       ,  7.55000019,  7.5999999 ,  7.6500001 ,  7.69999981, &
                7.75      ,  7.80000019,  7.8499999 ,  7.9000001 ,  7.94999981, &
                8.        ,  8.05000019,  8.10000038,  8.14999962,  8.19999981, &
                8.25      ,  8.30000019,  8.35000038,  8.39999962,  8.44999981, &
                8.5       ,  8.55000019,  8.60000038,  8.64999962,  8.69999981, &
                8.75      ,  8.80000019,  8.85000038,  8.89999962,  8.94999981, 9./

    !unit = DN cm^5 s^-1 pix^-1
    data aia_131 /3.18403601e-37,   3.22254703e-36,   2.61657920e-35, &
                  1.59575286e-34,   6.65779556e-34,   2.07015132e-33, &
                  6.05768615e-33,   1.76074833e-32,   4.52633001e-32, &
                  8.57121883e-32,   1.09184271e-31,   1.10207963e-31, &
                  1.11371658e-31,   1.29105226e-31,   1.80385897e-31, &
                  3.27295431e-31,   8.92002136e-31,   3.15214579e-30, &
                  9.73440787e-30,   2.22709702e-29,   4.01788984e-29, &
                  6.27471832e-29,   8.91764995e-29,   1.18725647e-28, &
                  1.52888040e-28,   2.05082946e-28,   3.47651873e-28, &
                  8.80482184e-28,   2.66533063e-27,   7.05805149e-27, &
                  1.46072515e-26,   2.45282476e-26,   3.55303726e-26, &
                  4.59075911e-26,   5.36503515e-26,   5.68444094e-26, &
                  5.47222296e-26,   4.81119761e-26,   3.85959059e-26, &
                  2.80383406e-26,   1.83977650e-26,   1.11182849e-26, &
                  6.50748885e-27,   3.96843481e-27,   2.61876319e-27, &
                  1.85525324e-27,   1.39717024e-27,   1.11504283e-27, &
                  9.38169611e-28,   8.24801234e-28,   7.43331919e-28, &
                  6.74537063e-28,   6.14495760e-28,   5.70805277e-28, &
                  5.61219786e-28,   6.31981777e-28,   9.19747307e-28, &
                  1.76795732e-27,   3.77985446e-27,   7.43166191e-27, &
                  1.19785603e-26,   1.48234676e-26,   1.36673114e-26, &
                  9.61047146e-27,   5.61209353e-27,   3.04779780e-27, &
                  1.69378976e-27,   1.02113491e-27,   6.82223774e-28, &
                  5.02099099e-28,   3.99377760e-28,   3.36279037e-28, &
                  2.94767378e-28,   2.65740865e-28,   2.44396277e-28, &
                  2.28003967e-28,   2.14941419e-28,   2.04178995e-28, &
                  1.95031045e-28,   1.87011994e-28,   1.79777869e-28, &
                  1.73093957e-28,   1.66795789e-28,   1.60785455e-28, &
                  1.55002399e-28,   1.49418229e-28,   1.44022426e-28, &
                  1.38807103e-28,   1.33772767e-28,   1.28908404e-28, &
                  1.24196208e-28,   1.17437501e-28,   1.12854330e-28, &
                  1.08410498e-28,   1.04112003e-28,   9.99529904e-29, &
                  9.59358806e-29,   9.20512291e-29,   8.83009123e-29, &
                  8.46817043e-29,   8.11921928e-29/




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

   ! ! output divB1
   ! call divvector(Btotal,ixI^L,ixO^L,divb)
   ! w(ixO^S,nw+3)=0.5d0*divb(ixO^S)/dsqrt(B2(ixO^S))/(^D&1.0d0/dxlevel(^D)+)
   ! ! output the plasma beta p*2/B**2
   ! w(ixO^S,nw+4)=pth(ixO^S)*two/B2(ixO^S)
    ! output heating rate
    t=global_time
    !call getbQ(ens,ixI^L,ixO^L,t,w,x)
    call getlQ(lQgrid,ixI^L,ixO^L,t,w,x)
    !w(ixO^S,nw+3)=ens(ixO^S)
    w(ixO^S,nw+3)=lQgrid(ixO^S)
   ! ! store the cooling rate 
   ! if(mhd_radiative_cooling)call getvar_cooling(ixI^L,ixO^L,w,x,ens)
   ! w(ixO^S,nw+7)=ens(ixO^S)

   ! ! store current
   ! call curlvector(Btotal,ixI^L,ixO^L,curlvec,idirmin,1,ndir)
   ! do idir=1,ndir
   !   w(ixO^S,nw+7+idir)=curlvec(ixO^S,idir)
   ! end do

   ! !SXR
   ! w(ixO^S,nw+11)=0
   ! EM(ixO^S)=(w(ixO^S,rho_)*unit_numberdensity)**2
   ! kbT(ixO^S)=kb*w(ixO^S,nw+1)*unit_temperature/keV
   ! do ihMu=0,num_hMu
   !   hMu=ihMu*dhMu+hMu_min
   !   gff(ixO^S)=1
   !   {do ixO^DB=ixOmin^DB,ixOmax^DB\}
   !     if (kbT(ixO^D)<hMu) then
   !       gff(ixO^D)=(kbT(ixO^D)/hMu)**0.4
   !     endif
   !   {enddo\}
   !   I_SXR(ixO^S)=I0*EM(ixO^S)*gff(ixO^S)*exp(-hMu/(kbT(ixO^S)))/(hMu*sqrt(kbT(ixO^S)))
   !   w(ixO^S,nw+11)=w(ixO^S,nw+11)+I_SXR(ixO^S)*dhMu
   ! enddo

   ! !AIA 131
   ! {do ixO^DB=ixOmin^DB,ixOmax^DB\}
   !   logT=log10(w(ixO^D,nw+1)*unit_temperature)
   !   do iTl=1,100
   !     if (logT<aia_T(1)) then
   !       logGT=-50
   !     else  if (logT>=aia_T(iTl) .and. logT<aia_T(iTl+1)) then
   !       logGT=log10(aia_131(iTl))*(logT-aia_T(iTl+1))/(aia_T(iTl)-aia_T(iTl+1))+&
   !             log10(aia_131(iTl+1))*(logT-aia_T(iTl))/(aia_T(iTl+1)-aia_T(iTl))
   !     endif
   !   enddo
   !   w(ixO^D,nw+12)=EM(ixO^D)*(10**(logGT))
   ! {enddo\}


  
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='Te Alfv lQ'
    !varnames='Te Alfv divB beta bQ lQ rad j1 j2 j3 SXR AIA131'

  end subroutine specialvarnames_output

  subroutine special_output(qunit)
    use mod_global_parameters

    integer, intent(in) :: qunit

    call output_Bfield(qunit)

  end subroutine special_output

  subroutine output_Bfield(qunit)
    ! output magnetic field lines
    use mod_usr_methods
    use mod_global_parameters

    integer          :: ix^D

    integer, intent(in) :: qunit
    character(len=30)   :: filename
    integer :: iigrid,igrid
    integer :: filenr
    logical :: fileopen

    if (mype==0) then
      inquire(qunit,opened=fileopen)
      if(.not.fileopen)then
        ! generate filename 
        filenr=snapshotini
        if (autoconvert) filenr=snapshotnext
        write(filename,'(a,i4.4,a)') trim(base_filename),filenr,".txt"
      endif
      open(1,file=filename)
      write(1,*) global_time
      write(1,*) numX^D  
      write(1,*) 'xL yL zL xR yR zL NpL NpR QeL QeR'
      {do ix^D=1,numX^D \}
        write(1,'(e15.7, e15.7, e15.7, e15.7, e15.7, e15.7, e15.7, e15.7, e15.7, e15.7)') &
        {xFL(ix^DD,^D),},{xFR(ix^DD,^D),},NpL(ix^D),NpR(ix^D),QeL(ix^D),QeR(ix^D)
      {enddo \}
      close(1)
    endif

  end subroutine output_Bfield

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a steady (time-independent) potential or 
  ! linear force-free background field
    use mod_global_parameters

    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    wB0(ixO^S,2)=-Busr*dcos(kx*x(ixO^S,2))*dexp(-ly*x(ixO^S,3))*dcos(theta)
    wB0(ixO^S,3)= Busr*dsin(kx*x(ixO^S,2))*dexp(-ly*x(ixO^S,3))
    wB0(ixO^S,1)=-Busr*dcos(kx*x(ixO^S,2))*dexp(-ly*x(ixO^S,3))*dsin(theta)

  end subroutine specialset_B0

  subroutine specialset_J0(ixI^L,ixO^L,x,wJ0)
  ! Here add a time-independent background current density 
    use mod_global_parameters

    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wJ0(ixI^S,7-2*ndir:ndir)

    !wJ0(ixO^S,2)= ly*Busr*dcos(kx*x(ixO^S,2))*dexp(-ly*x(ixO^S,3))*dsin(theta)
    !wJ0(ixO^S,3)=-kx*Busr*dsin(kx*x(ixO^S,2))*dexp(-ly*x(ixO^S,3))*dsin(theta)
    !wJ0(ixO^S,1)= kx*Busr*dcos(kx*x(ixO^S,2))*dexp(-ly*x(ixO^S,3))&
    !             -ly*Busr*dcos(kx*x(ixO^S,2))*dexp(-ly*x(ixO^S,3))*dcos(theta)
    wJ0=0.d0

  end subroutine specialset_J0

end module mod_usr
