module mod_usr
  use mod_mhd
  implicit none
  double precision, allocatable :: pbc(:),rbc(:)
  double precision :: usr_grav
  double precision :: heatunit,gzone,theta,SRadius,kx,ly,bQ0,dya
  double precision, allocatable :: pa(:),ra(:),ya(:)
  integer, parameter :: jmax=8000
  integer :: numFL,numLH
  double precision :: Ec,Eb,delta_l,delta_u,Ec_erg,Eb_erg,dE,mu0
  double precision :: tmax,tw1,tw2,Fleft,Fright
  double precision, allocatable :: xL(:,:),yL(:,:),xR(:,:),yR(:,:)
  double precision, allocatable :: QeL(:,:),QeR(:,:),ApoL(:),ApoR(:)
  double precision, allocatable :: NpL(:,:),NpR(:,:)
  double precision, allocatable :: NcolL(:,:),NcolR(:,:)
  double precision :: dyL,dyR,dApo
  integer :: numN,numE
  double precision, allocatable :: QeN(:),Ncol(:)
  double precision :: Nmin,Nmax,dNlog
  integer :: ccount=0
  double precision :: time_save=0
  integer :: count_save=0
  double precision :: t_update_Qe=-1.0e-7,dt_update_Qe=5.0e-4
  double precision, allocatable :: Ee(:,:),spectra(:,:),muN(:,:)


contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    call set_coordinate_system("Cartesian_2.5D")

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
    dya=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax) ! cells size of high-resolution 1D solar atmosphere
    Busr=Busr/unit_magneticfield ! magnetic field strength at the bottom
    theta=30.d0*dpi/180.d0 ! the angle to the plane xy, 90-theta is the angle to the polarity inversion line of the arcade 
    kx=dpi/(xprobmax1-xprobmin1)
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

    call init_flare_heating()

  end subroutine inithdstatic

  subroutine init_flare_heating
    ! setting for flare heating
    use mod_global_parameters

    integer :: j,iBx,numBx,iBc,iBb,iNcol
    double precision :: Atom,temp
    double precision :: Cllog,gammaN,beta,K,Nc,Nb,Efactor
    double precision :: dtx,tx,au,al,b,maxh
    double precision, allocatable :: Bu(:),Bl(:)
    double precision :: xLmin,xLmax,xRmin,xRmax
    double precision :: dxL,dxR
    integer :: ix,iy,refine_factor

    ! heating parameter
    Fleft=1.1e11    ! energy flux at left footpoint
    Fright=1.3e11   ! energy flux at right footpoint
    tmax=60/unit_time   ! when flux reach a maximum value
    tw1=30/unit_time    ! time scale for flux increasing before tmax
    tw2=180/unit_time   ! time scale for flux decreasing after tmax

    ! electron's spectrum [electrons cm^-2 keV-1]
    delta_l=3.0 ! spectral index for E<Eb
    delta_u=4.0 ! spectral index for E>Eb
    Ec=20.0  ! cutoff energy [keV]
    Eb=100.0  ! break energy [keV]
    numE=280
    dE=1.0    ! [keV]

    ! looking for a magnetic field line
    ! initialize heating rate
    xRmin=2.1d0
    xRmax=2.2d0
    xLmin=-xRmax
    xLmax=-xRmin
    refine_factor=2**(refine_max_level-1)
    dyL=(xprobmax2-xprobmin2)/(domain_nx2*refine_factor)
    dyR=dyL
    maxh=1.0
    numFL=int(2*(xRmax-xRmin)/dyL)+1
    numLH=int(maxh/dyL)
    allocate(ApoL(numFl),ApoR(numFl))
    dApo=abs(Busr/kx*cos(kx*xRmin)-Busr/kx*cos(kx*xRmax))/(numFL-1)
    do j=1,numFL
      ApoL(j)=Busr/kx*cos(kx*xRmax)+dApo*(j-1)
      ApoR(j)=Busr/kx*cos(kx*xRmax)+dApo*(j-1)
    enddo

    allocate(xL(numFL,numLH),yL(numFL,numLH))    
    allocate(xR(numFL,numLH),yR(numFL,numLH))
    allocate(QeL(numFL,numLH),QeR(numFL,numLH))
    allocate(NpL(numFL,numLH),NpR(numFL,numLH))
    allocate(NcolL(numFL,numLH),NcolR(numFL,numLH))
    do j=1,numLH
      yL(:,j)=dyL*j
      yR(:,j)=dyR*j
    enddo
    do ix=1,numFL
      Atom=ApoL(ix)
      do iy=1,numLH
        temp=acos(Atom/(Busr/kx*exp(-(yL(ix,iy)*ly))))
        xL(ix,iy)=-temp/kx
        xR(ix,iy)=temp/kx
      enddo
    enddo
    QeL=0
    QeR=0    
    NcolL=1.0e18
    NcolR=1.0e18

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
      muN(1,iE)=mu0
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

    !! output e flux
    !if (mype==0) then
    !  write(fname,'(a17)') 'spectra_table.txt'
    !  open(1,file=fname)
    !  write(1,*) 'numN numE'
    !  write(1,*) numN,numE
    !  do iNcol=1,numN
    !    write(1,*) Ncol(iNcol)
    !    write(1,*) 'E spectra'
    !    do iE=1,numE
    !      write(1,'(e15.7, e15.7)') , Ee(iNcol,iE),spectra(iNcol,iE)
    !    enddo
    !  enddo
    !  close(1)
    !endif


    !! output pitch angle
    !if (mype==0) then
    !  write(fname,'(a12)') 'mu_table.txt'
    !  open(1,file=fname)
    !  write(1,*) 'numN numE'
    !  write(1,*) numN,numE
    !  do iNcol=1,numN
    !    write(1,*) Ncol(iNcol)
    !    write(1,*) 'E spectra'
    !    do iE=1,numE
    !      write(1,'(e15.7, e15.7)') , Ee(iNcol,iE),muN(iNcol,iE)
    !    enddo
    !  enddo
    !  close(1)
    !endif




  end subroutine get_spectra

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
        w(ixO^S,mom(idir)) =-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(idir))&
                   /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
      end do
      ! fixed b1 b2 b3
      if(B0field) then
        w(ixO^S,mag(:))=0.d0
      else
        call specialset_B0(ixI^L,ixO^L,x,Bf)
        w(ixO^S,mag(1:ndir))=Bf(ixO^S,1:ndir)
      endif
      ! fixed gravity stratification of density and pressure pre-determined in initial condition
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
        w(ixO^S,mom(idir)) =-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,mom(idir))&
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
    if(iprob > 1)then
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

    bQgrid(ixO^S)=bQ0*dexp(-x(ixO^S,2)/5.d0)

  end subroutine getbQ

  subroutine get_flare_eflux()
    ! calculate electron deposit energy flux
    use mod_global_parameters

    integer          :: ixO^L,ixO^D,j

    integer :: ixFL,iyLH
    double precision :: F0LT(numFL),F0RT(numFL)
    double precision :: ixHalf,numQuar,qt

    ixHalf=(numFL+1)/2.0
    numQuar=numFL/4.0
    ! energy flux in the left/right footpoints [erg cm^-2 s^-1]
    do ixFL=1,numFL
      F0LT(ixFL)=Fleft
      F0RT(ixFL)=Fright
    enddo

    qt=global_time
    if (qt<=tmax)then
      F0LT=F0LT*exp(-(qt-tmax)**2/tw1**2)
      F0RT=F0RT*exp(-(qt-tmax)**2/tw1**2)
    else
      F0LT=F0LT*exp(-(qt-tmax)**2/tw2**2)
      F0RT=F0RT*exp(-(qt-tmax)**2/tw2**2)
    endif

    ! update distribution of the energy flux
    call update_e_distr(xL,yL,QeL,1)
    call update_e_distr(xR,yR,QeR,2)

    do ixFL=1,numFL
      QeL(ixFL,:)=F0LT(ixFL)*QeL(ixFL,:)/heatunit    ! [erg cm^-2 s^-1] to defined unit
      QeR(ixFL,:)=F0RT(ixFL)*QeR(ixFL,:)/heatunit    ! [erg cm^-2 s^-1] to defined unit
    enddo


  end subroutine get_flare_eflux


  subroutine update_e_distr(xLR,yLR,QeLR,iLR)
    ! calculate electron deposit energy flux
    use mod_global_parameters

    integer          :: ixO^L,ixO^D,j

    double precision :: Np(numFL,numLH),Nh(numFL,numLH)
    double precision :: xLR(numFL,numLH),yLR(numFL,numLH)
    double precision :: QeLR(numFL,numLH)
    double precision :: dl,dlx,dly,dlz,dlxy,ry
    integer :: ixFL,iyLH
    double precision :: sumQe
    integer :: iLR,iNlog

    call get_Np_profile(xLR,yLR,Np)

    ! column depth along field lines
    Nh=0
    do ixFL=1,numFL
      do iyLH=numLH-1,1,-1
        dlx=xLR(ixFL,iyLH+1)-xLR(ixFL,iyLH)
        dly=yLR(ixFL,iyLH+1)-yLR(ixFL,iyLH)
        dlz=cos(theta)*dlx
        dl=sqrt(dlx**2+dly**2+dlz**2)
        !dl=sqrt((xL(ixFL,iyLH+1)-xL(ixFL,iyLH))**2 + &
        !        (yL(ixFL,iyLH+1)-yL(ixFL,iyLH))**2)
        dl=dl*unit_length
        Nh(ixFL,iyLH)=Nh(ixFL,iyLH+1) + & 
                      (Np(ixFL,iyLH)+Np(ixFL,iyLH+1))*dl/2.0
      enddo
    enddo

    ! energy flux as a function of height
    do ixFL=1,numFL
      do iyLH=1,numLH
        dlx=xLR(ixFL,iyLH+1)-xLR(ixFL,iyLH)
        dly=yLR(ixFL,iyLH+1)-yLR(ixFL,iyLH)
        dlz=cos(theta)*dlx
        dlxy=sqrt(dlx**2+dly**2)
        dl=sqrt(dlx**2+dly**2+dlz**2)
        ry=abs(xLR(numFL,iyLH)-xLR(1,iyLH))/abs(xLR(numFL,1)-xLR(1,1))
        iNlog=int(log10(Nh(ixFL,iyLH)/Nmin)/dNlog+0.5)+1
        if (iNlog<1) then
          iNlog=1
        endif
        QeLR(ixFL,iyLH)=QeN(iNlog)*Np(ixFL,iyLH)*(dl/dlxy)/ry
      enddo
    enddo

    if (iLR==1) then
      NpL=Np
      NcolL=Nh
    else 
      NpR=Np
      NcolR=Nh
    endif

  end subroutine update_e_distr

  subroutine get_Np_profile(xLR,yLR,Np)
    ! calculate electron deposit energy flux
    use mod_global_parameters

    integer          :: ixO^L,ixO^D,j

    integer :: iigrid,igrid,igrid_prev,igrid_now
    double precision :: Npsend(numFL,numLH),Nprecv(numFL,numLH)
    double precision :: Np(numFL,numLH),Nh(numFL,numLH)
    double precision :: xLR(numFL,numLH),yLR(numFL,numLH)
    double precision :: xLR_new(numFL,numLH),yLR_new(numFL,numLH)
    double precision :: QeLR(numFL,numLH)
    double precision :: dl,xbmin,xbmax,ybmin,ybmax,dxb,dyb
    double precision :: dlx,dly,dlz,dlxy,xd,yd
    double precision :: xLmin,xLmax,yLmin,yLmax
    integer :: iymin,iymax,ixb,iyb,iNlog,numS,ixFL,iyLH
    integer :: ixbl,ixbu,iybl,iybu
    double precision :: sumQe 
    integer :: iLR
    double precision :: Np00,Np10,NP01,Np11

    ixOmin1=ixmlo1
    ixOmin2=ixmlo2
    ixOmax1=ixmhi1
    ixOmax2=ixmhi2

    ! looking for density along magnetic field lines
    Npsend=0
    do ixFL=1,numFL
      igrid_prev=-1
      do iyLH=1,numLH
        if (igrid_prev==-1) then
          ! for iyLH==1
          LOOP1: do iigrid=1,igridstail; igrid=igrids(iigrid);
            call update_block_para(igrid,dxb,dyb,xbmin,xbmax,ybmin,ybmax)
            if (xLR(ixFL,iyLH)>=xbmin .and. xLR(ixFL,iyLH)<xbmax .and. &
                yLR(ixFL,iyLH)>=ybmin .and. yLR(ixFL,iyLH)<ybmax) then
              ixbl=floor((xLR(ixFL,iyLH)-ps(igrid)%x(ixOmin1,ixOmin2,1))/dxb)+ixOmin1 
              iybl=floor((yLR(ixFL,iyLH)-ps(igrid)%x(ixOmin1,ixOmin2,2))/dyb)+ixOmin2 
              ixbu=ixbl+1
              iybu=iybl+1
              xd=(xLR(ixFL,iyLH)-ps(igrid)%x(ixbl,iybl,1))/dxb
              yd=(yLR(ixFL,iyLH)-ps(igrid)%x(ixbl,iybl,2))/dyb
              Np00=ps(igrid)%w(ixbl,iybl,rho_)*(1-xd)*(1-yd)
              Np10=ps(igrid)%w(ixbu,iybl,rho_)*xd*(1-yd)
              Np01=ps(igrid)%w(ixbl,iybu,rho_)*(1-xd)*yd
              Np11=ps(igrid)%w(ixbu,iybu,rho_)*xd*yd
              Npsend(ixFL,iyLH)=Np00+Np10+Np01+Np11
              igrid_prev=igrid
              exit LOOP1
            endif
          end do LOOP1
        else
          ! for iyLH>1
          igrid=igrid_prev
          call update_block_para(igrid,dxb,dyb,xbmin,xbmax,ybmin,ybmax)
          if (xLR(ixFL,iyLH)>=xbmin .and. xLR(ixFL,iyLH)<xbmax .and. &
              yLR(ixFL,iyLH)>=ybmin .and. yLR(ixFL,iyLH)<ybmax) then
            ! in the same block with previous point
            ixbl=floor((xLR(ixFL,iyLH)-ps(igrid)%x(ixOmin1,ixOmin2,1))/dxb)+ixOmin1 
            iybl=floor((yLR(ixFL,iyLH)-ps(igrid)%x(ixOmin1,ixOmin2,2))/dyb)+ixOmin2 
            ixbu=ixbl+1
            iybu=iybl+1
            xd=(xLR(ixFL,iyLH)-ps(igrid)%x(ixbl,iybl,1))/dxb
            yd=(yLR(ixFL,iyLH)-ps(igrid)%x(ixbl,iybl,2))/dyb
            Np00=ps(igrid)%w(ixbl,iybl,rho_)*(1-xd)*(1-yd)
            Np10=ps(igrid)%w(ixbu,iybl,rho_)*xd*(1-yd)
            Np01=ps(igrid)%w(ixbl,iybu,rho_)*(1-xd)*yd
            Np11=ps(igrid)%w(ixbu,iybu,rho_)*xd*yd
            Npsend(ixFL,iyLH)=Np00+Np10+Np01+Np11
            igrid_prev=igrid
          else
            ! in different block
            LOOP2: do iigrid=1,igridstail; igrid=igrids(iigrid);
              call update_block_para(igrid,dxb,dyb,xbmin,xbmax,ybmin,ybmax)
              if (xLR(ixFL,iyLH)>=xbmin .and. xLR(ixFL,iyLH)<xbmax .and. &
                  yLR(ixFL,iyLH)>=ybmin .and. yLR(ixFL,iyLH)<ybmax) then
                ixbl=floor((xLR(ixFL,iyLH)-ps(igrid)%x(ixOmin1,ixOmin2,1))/dxb)+ixOmin1 
                iybl=floor((yLR(ixFL,iyLH)-ps(igrid)%x(ixOmin1,ixOmin2,2))/dyb)+ixOmin2 
                ixbu=ixbl+1
                iybu=iybl+1
                xd=(xLR(ixFL,iyLH)-ps(igrid)%x(ixbl,iybl,1))/dxb
                yd=(yLR(ixFL,iyLH)-ps(igrid)%x(ixbl,iybl,2))/dyb
                Np00=ps(igrid)%w(ixbl,iybl,rho_)*(1-xd)*(1-yd)
                Np10=ps(igrid)%w(ixbu,iybl,rho_)*xd*(1-yd)
                Np01=ps(igrid)%w(ixbl,iybu,rho_)*(1-xd)*yd
                Np11=ps(igrid)%w(ixbu,iybu,rho_)*xd*yd
                Npsend(ixFL,iyLH)=Np00+Np10+Np01+Np11
                igrid_prev=igrid
                exit LOOP2
              endif
            end do LOOP2
          endif
        endif        
      enddo
    enddo

    ! update density data
    numS=numFL*numLH
    call MPI_ALLREDUCE(Npsend,Nprecv,numS,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,icomm,ierrmpi)
    Np=Nprecv*unit_numberdensity



  end subroutine get_Np_profile


  subroutine update_block_para(igrid,dxb,dyb,xbmin,xbmax,ybmin,ybmax)
    ! update parameters of the block
    use mod_global_parameters

    integer          :: ixO^L,ixO^D,j
    integer :: igrid
    double precision :: dxb,dyb,xbmin,xbmax,ybmin,ybmax

    ixOmin1=ixmlo1
    ixOmin2=ixmlo2
    ixOmax1=ixmhi1
    ixOmax2=ixmhi2

    dxb=rnode(rpdx1_,igrid)
    dyb=rnode(rpdx2_,igrid)
    xbmin=ps(igrid)%x(ixOmin1,ixOmin2,1)-dxb/2.0
    xbmax=ps(igrid)%x(ixOmax1,ixOmax2,1)+dxb/2.0
    ybmin=ps(igrid)%x(ixOmin1,ixOmin2,2)-dyb/2.0
    ybmax=ps(igrid)%x(ixOmax1,ixOmax2,2)+dyb/2.0

  end subroutine update_block_para

  subroutine special_global(iit,qt)
    !! run at the begining of a step
    use mod_global_parameters

    integer, intent(in) :: iit
    double precision, intent(in) :: qt
    character (20) :: fname
    integer :: ixFL,iyLH
    
    if (iprob>=6 .and. convert) then
      call get_flare_eflux()
    endif

    if (iprob>=6) then
      if (qt>t_update_Qe+10*dt_update_Qe) then
        t_update_Qe=qt-1.0e-7
      endif

      if (qt>t_update_Qe) then
        call get_flare_eflux()
        t_update_Qe=t_update_Qe+dt_update_Qe
      endif

      !! output e flux
      !if (mype==0) then
      !  if (global_time>=time_save) then
      !    time_save=time_save+dtsave(2)
      !    write(fname,'(a4,i04.4,a4)') 'QehL',count_save,'.txt'
      !    open(1,file=fname)
      !    write(1,*) 'numFL','numLH'
      !    write(1,*) numFL,numLH
      !    do ixFL=1,numFL
      !      write(1,*) ixFL
      !      write(1,*) 'x y rho QeH'
      !      do iyLH=1,numLH
      !        write(1,'(e15.7, e15.7, e15.7, e15.7)') xL(ixFL,iyLH), &
      !              yL(ixFL,iyLH),NpL(ixFL,iyLH),QeL(ixFL,iyLH)
      !      enddo
      !    enddo
      !    close(1)
      !    count_save=count_save+1
      !  endif
      !endif
    endif

  end subroutine

  subroutine getlQ(lQgrid,ixI^L,ixO^L,qt,w,x)
    !!calculate localized heating at chromosphere
    use mod_global_parameters
    use mod_radiative_cooling

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)

    double precision :: lQgrid(ixI^S),ens(ixI^S)
    double precision :: asym,lQt,lQd,lQst,hp,flash,ly,Apo,Atop,Alow
    double precision :: sscale,tonset,tstop,tflare,tpeak,tscale,fstop,thstop
    double precision :: lQheat,lQt1,lQt2,lQtw
    integer          :: ixO^D,j,ixFL,iyLH

    !-----------------------------------------------------------------------------
    ! 180 s of the flare heating
    lQt=5.d1/unit_time
    lQst=18.d1/unit_time
    lQt1=3.d1/unit_time
    lQt2=15.d1/unit_time
    lQtw=60.d0/unit_time
    lQd=0.05d0
    asym=0.8d0
    hp=0.175d0
    lQgrid=0.d0
    ly=kx*dcos(theta)
    Atop=Busr/kx*cos(kx*2.2d0)
    Alow=Busr/kx*cos(kx*2.3d0)
    lQheat=180.0d0/heatunit
      
    select case(iprob)
    case(2)
    if(qt<lQst) then
      {do ixO^DB=ixOmin^DB,ixOmax^DB\}
        Apo=Busr/kx*cos(kx*x(ixO^D,1))*dexp(-ly*x(ixO^D,2))
        if(Apo>Alow .and. Apo<Atop) then
          lQgrid(ixO^D)=lQheat*dexp(-(x(ixO^D,2)-hp)**2/lQd)
        endif
      {enddo\}
      if(qt<lQt) then
        lQgrid(ixO^S)=lQgrid(ixO^S)*(qt/lQt)
      endif
      !lQgrid(ixO^S)=lQgrid(ixO^S)*dexp(-(qt-lQst/2)**2/0.15d0)/sqrt(0.15*dpi)
    endif

    case(3)
    if(qt<lQst) then
      {do ixO^DB=ixOmin^DB,ixOmax^DB\}
        Apo=Busr/kx*cos(kx*x(ixO^D,1))*dexp(-ly*x(ixO^D,2))
        if(Apo>Alow .and. Apo<Atop .and. x(ixO^D,1)>zero) then
          lQgrid(ixO^D)=lQheat*dexp(-(x(ixO^D,2)-hp)**2/lQd)
        endif
        if(Apo>Alow .and. Apo<Atop .and. x(ixO^D,1)<zero) then
          lQgrid(ixO^D)=lQheat*dexp(-(x(ixO^D,2)-hp)**2/lQd)*asym
        endif
      {enddo\}
      if(qt<lQt) then
        lQgrid(ixO^S)=lQgrid(ixO^S)*(qt/lQt)
      endif
      !lQgrid(ixO^S)=lQgrid(ixO^S)*dexp(-(qt-lQst/2)**2/0.15d0)/sqrt(0.15*dpi)
    endif

    case(4)
    if(qt<lQst) then
      {do ixO^DB=ixOmin^DB,ixOmax^DB\}
        Apo=Busr/kx*cos(kx*x(ixO^D,1))*dexp(-ly*x(ixO^D,2))
        if(Apo>Alow .and. Apo<Atop .and. x(ixO^D,1)>zero) then
          lQgrid(ixO^D)=lQheat*dexp(-(x(ixO^D,2)-hp)**2/lQd)
        endif
        if(Apo>Alow .and. Apo<Atop .and. x(ixO^D,1)<zero) then
          lQgrid(ixO^D)=lQheat*dexp(-(x(ixO^D,2)-hp)**2/lQd)*asym
        endif
      {enddo\}
      if (qt<lQt1) then
        lQgrid(ixO^S)=lQgrid(ixO^S)*(qt/lQt1)
      else if (qt<lQt2) then
        lQgrid(ixO^S)=lQgrid(ixO^S)
      else if (qt<lQst) then
        lQgrid(ixO^S)=lQgrid(ixO^S)*((lQst-qt)/lQt1)
      else
        lQgrid(ixO^S)=0
      endif
    endif

    case(5)
    {do ixO^DB=ixOmin^DB,ixOmax^DB\}
      Apo=Busr/kx*cos(kx*x(ixO^D,1))*dexp(-ly*x(ixO^D,2))
      if(Apo>Alow .and. Apo<Atop .and. x(ixO^D,1)>zero) then
        lQgrid(ixO^D)=(lQheat/sqrt(3.14d0*lQd))*dexp(-(x(ixO^D,2)-hp)**2/lQd)*(1/(1+asym))
      endif
      if(Apo>Alow .and. Apo<Atop .and. x(ixO^D,1)<zero) then
        lQgrid(ixO^D)=(lQheat/sqrt(3.14d0*lQd))*dexp(-(x(ixO^D,2)-hp)**2/lQd)*(asym/(1+asym))
      endif
    {enddo\}
    lQgrid(ixO^S)=lQgrid(ixO^S)*dsqrt(1/(3.14*lQtw**2))*dexp(-(qt-2.0*lQtw)**2/(lQtw**2))


    case(6)
    {do ixO^DB=ixOmin^DB,ixOmax^DB\}
      Apo=Busr/kx*cos(kx*x(ixO^D,1))*dexp(-ly*x(ixO^D,2))
      if(Apo>=ApoL(1) .and. Apo<=ApoL(numFL) .and. x(ixO^D,1)<zero .and. x(ixO^D,2)>zero) then
        ixFL=floor((Apo-ApoL(1))/dApo+0.5)+1
        iyLH=floor((x(ixO^D,2)-yL(ixFL,1))/dyL+0.5)+1
        if (iyLH>=1 .and. iyLH<=numLH) then
          lQgrid(ixO^D)=QeL(ixFL,iyLH)
        else
          lQgrid(ixO^D)=0.0
        endif
      endif

      if(Apo>=ApoL(1) .and. Apo<=ApoL(numFL) .and. x(ixO^D,1)>zero .and. x(ixO^D,2)>zero) then
        ixFL=floor((Apo-ApoR(1))/dApo+0.5)+1
        iyLH=floor((x(ixO^D,2)-yR(ixFL,1))/dyR+0.5)+1
        if (iyLH>=1 .and. iyLH<=numLH) then
          lQgrid(ixO^D)=QeR(ixFL,iyLH)
        else
          lQgrid(ixO^D)=0.0
        endif
      endif
    {enddo\}

    end select



  end subroutine getlQ

  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
  ! Enforce additional refinement or coarsening
  ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen
    double precision :: ly, Atop, Alow, psmall, rhosmall, Tesmall
    double precision :: pth(ixI^S), Te(ixI^S)

    ly=kx*dcos(theta)
    Atop=Busr/kx*cos(kx*1.9d0)
    Alow=Busr/kx*cos(kx*2.5d0)

    ! fix the bottom layer to the highest level
    if (iprob==1) then
      if (any(x(ixO^S,2)<=xprobmin2+0.05d0)) then
        refine=1
        coarsen=-1
      endif
    endif

    if (iprob>=3) then
      if (any(Busr/kx*cos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))>Alow .and. &
          Busr/kx*cos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))<Atop)) then
        if (any(x(ixO^S,1)>-1.0 .and. x(ixO^S,1)<1.0)) then
          if (level<6) then
            refine=1
            coarsen=-1
          else if (level==6) then
            refine=-1
            coarsen=-1
          else
            refine=-1
            coarsen=1
          endif
        else
          if (level<6) then
            refine=1
            coarsen=-1
          else if (level==6) then
            refine=-1
            coarsen=-1
          else
            refine=-1
            coarsen=1
          endif
        endif
      else if (any(x(ixO^S,2)<=xprobmin2+0.3d0)) then
        if (level<3) then
          refine=1
          coarsen=-1
        else if (level==3) then
          refine=-1
          coarsen=0
        else
          refine=-1
          coarsen=1
        endif
      else
        if (level==2) then
          refine=-1
          coarsen=0
        else if (level>2) then
          refine=-1
          coarsen=1
        else
          refine=0
          coarsen=0
        endif
      endif
    endif

  end subroutine special_refine_grid

  subroutine get_IRIS(EUV,ixI^L,ixO^L,w,x)
  ! synthesize IRIS 1354 emission line
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw)

    integer :: ix^D,ixO^D

    double precision :: EUV(ixI^S),EM(ixI^S)
    double precision :: Ttable(101),Gtable(101)
    double precision :: logT,logGT,LOS,Tu,Tl
    integer          :: iTt
    double precision :: depth,arcsec,AU

    depth=5.0e8     ! [cm]
    AU=1.496e13     ! 1 AU [cm]
    arcsec=7.35e7   ! 735 km

    data Ttable /4. ,  4.05, 4.1,  4.15, 4.2,  4.25, 4.3,  4.35, &
                 4.4,  4.45, 4.5,  4.55, 4.6,  4.65, 4.7,  4.75, &
                 4.8,  4.85, 4.9,  4.95, 5. ,  5.05, 5.1,  5.15, &
                 5.2,  5.25, 5.3,  5.35, 5.4,  5.45, 5.5,  5.55, &
                 5.6,  5.65, 5.7,  5.75, 5.8,  5.85, 5.9,  5.95, & 
                 6. ,  6.05, 6.1,  6.15, 6.2,  6.25, 6.3,  6.35, &
                 6.4,  6.45, 6.5,  6.55, 6.6,  6.65, 6.7,  6.75, &
                 6.8,  6.85, 6.9,  6.95, 7. ,  7.05, 7.1,  7.15, &
                 7.2,  7.25, 7.3,  7.35, 7.4,  7.45, 7.5,  7.55, &
                 7.6,  7.65, 7.7,  7.75, 7.8,  7.85, 7.9,  7.95, & 
                 8. ,  8.05, 8.1,  8.15, 8.2,  8.25, 8.3,  8.35, &
                 8.4,  8.45, 8.5,  8.55, 8.6,  8.65, 8.7,  8.75, &
                 8.8,  8.85, 8.9,  8.95, 9./

    data Gtable /0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, & 
                 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, & 
                 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, & 
                 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, & 
                 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, &
                 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, & 
                 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, & 
                 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, & 
                 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, & 
                 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, &
                 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 2.8306229e-43, & 
                 6.2249742e-40, 5.2502353e-37, 1.5532066e-34, 1.6648678e-32, & 
                 7.8517756e-31, 2.0530806e-29, 3.4992089e-28, 4.2820552e-27, & 
                 3.9942179e-26, 2.9558137e-25, 1.7773391e-24, 8.8169698e-24, & 
                 3.6290640e-23, 1.2338631e-22, 3.4156004e-22, 7.5194933e-22, &
                 1.2748505e-21, 1.6059518e-21, 1.4685663e-21, 9.9275103e-22, & 
                 5.3269078e-22, 2.4902077e-22, 1.0922949e-22, 4.7030177e-23, & 
                 2.0328035e-23, 8.9087640e-24, 3.9651675e-24, 1.7899798e-24, & 
                 8.1571408e-25, 3.7315855e-25, 1.7037092e-25, 7.7065973e-26, & 
                 3.4425015e-26, 1.5093498e-26, 6.4718403e-27, 2.7138098e-27, &
                 1.1132269e-27, 4.4855265e-28, 1.7830170e-28, 7.0251280e-29, & 
                 2.7573307e-29, 1.0824557e-29, 4.2623914e-30, 1.6892384e-30, & 
                 6.7457206e-31, 2.7204369e-31, 1.1077821e-31, 4.5605646e-32, & 
                 1.8963764e-32, 7.9683214e-33, 3.3808639e-33, 1.4484216e-33, & 
                 6.2624282e-34, 2.7303008e-34, 1.1997758e-34, 5.3099278e-35, &
                 2.3657412e-35/

    EUV(ixO^S)=0
    EM(ixO^S)=(w(ixO^S,rho_)*unit_numberdensity)**2
    {do ixO^DB=ixOmin^DB,ixOmax^DB\}
      logT=log10(w(ixO^D,nw+1)*unit_temperature)
      do iTt=1,100
        if (logT<6.1) then
          logGT=-50
        else if (logT>=Ttable(iTt) .and. logT<Ttable(iTt+1)) then
          Tu=Gtable(iTt+1)
          Tl=Ttable(iTt)
          logGT=log10(Gtable(iTt))*(logT-Tu)/(Tl-Tu)+&
                log10(Gtable(iTt+1))*(logT-Tl)/(Tu-Tl)
        endif
      enddo
      EUV(ixO^D)=EM(ixO^D)*(10**(logGT))
    {enddo\}

  end subroutine get_IRIS

  subroutine get_AIA(AIAgrid,ixI^L,ixO^L,w,x)
  ! synthesize AIA 131 emission line
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw)

    integer :: ix^D,ixO^D

    double precision :: AIAgrid(ixI^S),EM(ixI^S)
    double precision :: aia_T(101),aia_131(101)
    double precision :: logT,logGT,LOS
    integer          :: iTl

   !aia T list
    data aia_T /4.        ,  4.05000019,  4.0999999 ,  4.1500001 , & 
                4.19999981,  4.25      ,  4.30000019,  4.3499999 , &  
                4.4000001 ,  4.44999981,  4.5       ,  4.55000019, &  
                4.5999999 ,  4.6500001 ,  4.69999981,  4.75      , &  
                4.80000019,  4.8499999 ,  4.9000001 ,  4.94999981, &
                5.        ,  5.05000019,  5.0999999 ,  5.1500001 , &  
                5.19999981,  5.25      ,  5.30000019,  5.3499999 , &  
                5.4000001 ,  5.44999981,  5.5       ,  5.55000019, &  
                5.5999999 ,  5.6500001 ,  5.69999981,  5.75      , &  
                5.80000019,  5.8499999 ,  5.9000001 ,  5.94999981, &
                6.        ,  6.05000019,  6.0999999 ,  6.1500001 , &  
                6.19999981,  6.25      ,  6.30000019,  6.3499999 , &  
                6.4000001 ,  6.44999981,  6.5       ,  6.55000019, &  
                6.5999999 ,  6.6500001 ,  6.69999981,  6.75      , &  
                6.80000019,  6.8499999 ,  6.9000001 ,  6.94999981, &
                7.        ,  7.05000019,  7.0999999 ,  7.1500001 , &  
                7.19999981,  7.25      ,  7.30000019,  7.3499999 , &  
                7.4000001 ,  7.44999981,  7.5       ,  7.55000019, &  
                7.5999999 ,  7.6500001 ,  7.69999981,  7.75      , &  
                7.80000019,  7.8499999 ,  7.9000001 ,  7.94999981, &
                8.        ,  8.05000019,  8.10000038,  8.14999962, &  
                8.19999981,  8.25      ,  8.30000019,  8.35000038, &  
                8.39999962,  8.44999981,  8.5       ,  8.55000019, &  
                8.60000038,  8.64999962,  8.69999981,  8.75      , &  
                8.80000019,  8.85000038,  8.89999962,  8.94999981, &
                9./

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

    !AIA 131
    EM(ixO^S)=(w(ixO^S,rho_)*unit_numberdensity)**2
    {do ixO^DB=ixOmin^DB,ixOmax^DB\}
      logT=log10(w(ixO^D,nw+1)*unit_temperature)
      do iTl=1,100
        if (logT<aia_T(1)) then
          logGT=-50
        else  if (logT>=aia_T(iTl) .and. logT<aia_T(iTl+1)) then
          logGT=log10(aia_131(iTl))*(logT-aia_T(iTl+1))/(aia_T(iTl)-aia_T(iTl+1))+&
                log10(aia_131(iTl+1))*(logT-aia_T(iTl))/(aia_T(iTl+1)-aia_T(iTl))
        endif
      enddo
      AIAgrid(ixO^D)=EM(ixO^D)*(10**(logGT))
    {enddo\}

  end subroutine get_AIA

  subroutine get_SXR(SXR,ixI^L,ixO^L,w,x,Ephl,Ephu)
  !synthesize thermal SXR
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw)

    integer :: ix^D,ixO^D

    double precision :: I0,hMu,hMu_max,hMu_min,dhMu,kb,kev,LOS
    double precision :: kbT(ixI^S),gff(ixI^S),I_SXR(ixI^S)
    double precision :: pth(ixI^S),Te(ixI^S)
    double precision :: EM(ixI^S),SXR(ixI^S)
    integer          :: ihMu,num_hMu,Ephl,Ephu
    double precision :: depth,Rs,AU,arcsec,sigma0

    kb=1.38d-23
    !I0=1.07d-18
    I0=1.07d-42     ! [cm^-2 s^-1 keV^-1]
    LOS=5.0d8
    keV=1.602d-16
    depth=5.0e8     ! [cm]
    Rs=6.955e10     ! solar radiu [cm]
    AU=1.496e13     ! 1 AU [cm]
    arcsec=7.35e7   ! 735 km
    I0=I0*arcsec*arcsec*depth

    hMu_max=1.0*Ephu
    hMu_min=1.0*Ephl
    dhMu=0.1
    num_hMu=int((hMu_max-hMu_min)/dhMu)

    ! output temperature
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    Te(ixO^S)=pth(ixO^S)/w(ixO^S,rho_)

    !SXR
    SXR(ixO^S)=0
    EM(ixO^S)=(w(ixO^S,rho_)*unit_numberdensity)**2
    kbT(ixO^S)=kb*Te(ixO^S)*unit_temperature/keV
    do ihMu=0,num_hMu
      hMu=ihMu*dhMu+hMu_min
      gff(ixO^S)=1
      {do ixO^DB=ixOmin^DB,ixOmax^DB\}
        if (kbT(ixO^D)<hMu) then
          gff(ixO^D)=(kbT(ixO^D)/hMu)**0.4
        endif
      {enddo\}
      I_SXR(ixO^S)=I0*EM(ixO^S)*gff(ixO^S)*exp(-hMu/(kbT(ixO^S)))/(hMu*sqrt(kbT(ixO^S)))
      SXR(ixO^S)=SXR(ixO^S)+I_SXR(ixO^S)*dhMu
    enddo

  end subroutine get_SXR

  subroutine get_HXR_BR(HXR_BR,ixI^L,ixO^L,w,x,Ephl,Ephu)
  ! synthesize bremsstrahlung HXR based on equation (7) in Appendix A of 
  ! Krucker et al. (2008)
  ! The distribution of energetic electrons is the same as that in footpoint
  ! heating
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw)

    integer :: ix^D,ixO^D

    double precision :: Fe0,FL0,FR0,Fe
    double precision :: qt,Bg0,const,const2,keV_erg
    double precision :: Btotal(ixI^S,1:ndir),Bxy(ixI^S)
    double precision :: HXR_BR(ixI^S)
    double precision :: pth(ixI^S),Te(ixI^S),Ne(ixI^S)
    integer :: iE,iEph,numEph,Ephl,Ephu
    double precision :: dEph

    double precision :: Eph,delta,alpha,mec2,r0,E0,N0,c,temp
    double precision :: depth,Rs,AU,arcsec,sigma0
    double precision :: sigmaKr,sigmaBH,ivEE
    double precision :: Apo
    integer :: iNlog,ixFL,iyLH
    double precision :: ry

    mec2=511.0      ! energy of a static electron [keV]
    alpha=1.0/137.0 ! fine structure constant
    c=2.9979d10     ! light speed [cm/s]
    r0=2.8179d-13   ! radiu of electron [cm]
    keV_erg=1.0d3*const_ev  ! 1 keV = * erg
    sigma0=7.9e-25  ! [cm^2 keV]
    depth=5.0e8     ! [cm]
    Rs=6.955e10     ! solar radiu [cm]
    AU=1.496e13     ! 1 AU [cm]
    arcsec=7.35e7   ! 735 km


    ! energy flux carried by fast electrons
    qt=global_time
    if (qt<=tmax)then
      Fe0=((Fleft+Fright)/2.0)*exp(-(qt-tmax)**2/tw1**2)
      FL0=Fleft*exp(-(qt-tmax)**2/tw1**2)
      FR0=Fright*exp(-(qt-tmax)**2/tw1**2)
    else
      Fe0=((Fleft+Fright)/2.0)*exp(-(qt-tmax)**2/tw2**2)
      FL0=Fleft*exp(-(qt-tmax)**2/tw2**2)
      FR0=Fright*exp(-(qt-tmax)**2/tw2**2)
    endif


    ! magnetic field magnitude for when Fe=Fe0
    Bg0=Busr

    ! fast electron spectra
    const=1.0*(delta_l-2.0)*(delta_u-2.0)/Ec_erg**2
    const=const/((delta_l-2.0)-((Eb/Ec)**(2.0-delta_l))*(delta_u-delta_l))
    const2=const*(Eb/Ec)**(delta_u-delta_l)

    ! local magnetic field
    if(B0field) then
      Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))+block%B0(ixI^S,1:ndir,0)
    else
      Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))
    endif
    Bxy(ixO^S)=sqrt((Btotal(ixO^S,1))**2+(Btotal(ixO^S,2))**2+(Btotal(ixO^S,3))**2)

    ! local temperature and density
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    Te(ixO^S)=pth(ixO^S)/w(ixO^S,rho_)*unit_temperature
    Ne(ixO^S)=w(ixO^S,rho_)*unit_numberdensity

    dEph=1.0  ! [keV]
    numEph=floor((Ephu-Ephl)/dEph+0.5)
    HXR_BR(ixO^S)=0

    ! HXR flux
    !do iEph=1,numEph
    !  ! for different photon energy
    !Eph=Ephl+(iEph-1.0)*dEph
    Eph=Ephl
    {do ixO^DB=ixOmin^DB,ixOmax^DB\}
      Apo=Busr/kx*cos(kx*x(ixO^D,1))*dexp(-ly*x(ixO^D,2))
 
     ! left foot
      if (Apo>=ApoL(1) .and. Apo<=ApoL(numFL) .and. & 
          x(ixO^D,1)<zero .and. x(ixO^D,2)<=1.0) then
        ixFL=floor((Apo-ApoL(1))/dApo+0.5)+1
        iyLH=floor((x(ixO^D,2)-yL(ixFL,1))/dyL+0.5)+1
        ry=abs(xL(numFL,iyLH)-xL(1,iyLH))/abs(xL(numFL,1)-xL(1,1))
        if (iyLH>=1 .and. iyLH<=numLH) then
          do iEph=1,numEph
            ! for different photon energy
            Eph=Ephl+(iEph-1.0)*dEph
            iNlog=int(log10(NcolL(ixFL,iyLH)/Nmin)/dNlog+0.5)+1
            if (iNlog>numN) then
              iNlog=numN
            else if (iNlog<1) then
              iNlog=1
            endif
            Fe=FL0/ry
            do iE=1,numE-1
              if (Ee(iNlog,iE)>Eph) then
                sigmaKr=sigma0/(Eph*Ee(iNlog,iE))
                temp=sqrt(1-Eph/Ee(iNlog,iE))
                sigmaBH=sigmaKr*log((1+temp)/(1-temp))
                dE=Ee(iNlog,iE+1)-Ee(iNlog,iE)
                HXR_BR(ixO^D)=HXR_BR(ixO^D)+ &
                              Ne(ixO^D)*Fe*spectra(iNlog,iE)*dE*sigmaBH*dEph/muN(iNlog,iE)
              endif
            enddo
          enddo
        endif
      endif

      ! right foot
      if (Apo>=ApoL(1) .and. Apo<=ApoL(numFL) .and. & 
               x(ixO^D,1)>zero .and. x(ixO^D,2)<=1.0) then   
        ixFL=floor((Apo-ApoR(1))/dApo+0.5)+1
        iyLH=floor((x(ixO^D,2)-yR(ixFL,1))/dyR+0.5)+1
        ry=abs(xR(numFL,iyLH)-xR(1,iyLH))/abs(xR(numFL,1)-xR(1,1))
        if (iyLH>=1 .and. iyLH<=numLH) then
          do iEph=1,numEph
            ! for different photon energy
            Eph=Ephl+(iEph-1.0)*dEph
            iNlog=int(log10(NcolR(ixFL,iyLH)/Nmin)/dNlog+0.5)+1
            if (iNlog>numN) then
              iNlog=numN
            else if (iNlog<1) then
              iNlog=1
            endif
            Fe=FR0/ry
            do iE=1,numE-1
              if (Ee(iNlog,iE)>Eph) then
                sigmaKr=sigma0/(Eph*Ee(iNlog,iE))
                temp=sqrt(1-Eph/Ee(iNlog,iE))
                sigmaBH=sigmaKr*log((1+temp)/(1-temp))
                dE=Ee(iNlog,iE+1)-Ee(iNlog,iE)
                HXR_BR(ixO^D)=HXR_BR(ixO^D)+ &
                              Ne(ixO^D)*Fe*spectra(iNlog,iE)*dE*sigmaBH*dEph/muN(iNlog,iE)
              endif
            enddo
          enddo
        endif
      endif

      ! other place
      if (x(ixO^D,2)>1.0 .or. Apo<ApoL(1) .or. Apo>ApoL(numFL)) then 
        if (Te(ixO^D)>5.0e6) then
          do iEph=1,numEph
            ! for different photon energy
            Eph=Ephl+(iEph-1.0)*dEph
            Fe=(Fe0*Bxy(ixO^D)/Bg0)/mu0
            iNlog=1
            do iE=1,numE-1
              if (Ee(1,iE)>Eph) then
                sigmaKr=sigma0/(Eph*Ee(iNlog,iE))
                temp=sqrt(1-Eph/Ee(iNlog,iE))
                sigmaBH=sigmaKr*log((1+temp)/(1-temp))
                dE=Ee(iNlog,iE+1)-Ee(iNlog,iE)
                HXR_BR(ixO^D)=HXR_BR(ixO^D)+ &
                              Ne(ixO^D)*Fe*spectra(iNlog,iE)*dE*sigmaBH*dEph
              endif
            enddo
          enddo
        endif
      endif
    {enddo\}
    !enddo

    !near earth
    HXR_BR=HXR_BR*depth   ! [cm^-3 s^-1] -> [cm^-2 s^-1]
    HXR_BR=HXR_BR*arcsec*arcsec   ! [cm^-2 s^-1] -> [arcsec^-2 s^-1]
    HXR_BR=HXR_BR/(4.0*dpi*AU*AU)

  end subroutine get_HXR_BR

  subroutine get_HXR_BR_old(HXR_BR,ixI^L,ixO^L,w,x,Ephl,Ephu)
  ! synthesize bremsstrahlung HXR based on equation (7) in Appendix A of 
  ! Krucker et al. (2008)
  ! The distribution of energetic electrons is N(Ee) ~ N0 (Ee/E0)^{-delta}
  ! local fast electron flux is a fuction of local plasma density
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw)

    integer :: ix^D,ixO^D

    double precision :: Eph,delta,alpha,mec2,r0,E0,N0,c,temp
    double precision :: HXR_BR(ixI^S)
    double precision :: pth(ixI^S),Te(ixI^S),Ne(ixI^S)
    double precision :: Eth,Nnth,Enth,dEph
    double precision :: Nc,Ec,p,dE,keV_erg,sigma0,FE
    !double precision, allocatable :: sigmaKr(:),sigmaBH(:)
    double precision :: Ee_(numE),ve(numE)
    integer :: iE,numE,iEph,numEph,Ephl,Ephu
    double precision :: depth,Rs,AU,arcsec

    !Eph=20.0        ! energy of HXR photons [keV]
    !N0=1.0d8        ! local fast electron density [cm^-3]
    Enth=4.0        ! above this energy the electrons are non-thermal [keV]
    E0=20.0         ! minimum value of fast electron energy
    p=2.5           ! fast electron flux spectral index
    delta=p+0.5     ! spectrum index of fast electrons
    mec2=511.0      ! energy of a static electron [keV]
    alpha=1.0/137.0 ! fine structure constant
    c=2.9979d10     ! light speed [cm/s]
    r0=2.8179d-13   ! radiu of electron [cm]
    keV_erg=1.0d3*const_ev  ! 1 keV = * erg
    sigma0=7.9e-25  ! [cm^2 keV]
    depth=5.0e8     ! [cm]
    Rs=6.955e10     ! solar radiu [cm]
    AU=1.496e13     ! 1 AU [cm]
    arcsec=7.35e7   ! 735 km


    ! output temperature
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    Te(ixO^S)=pth(ixO^S)/w(ixO^S,rho_)*unit_temperature
    Ne(ixO^S)=w(ixO^S,rho_)*unit_numberdensity

    HXR_BR(ixO^S)=0

    dEph=1.0  ! [keV]
    numEph=floor((Ephu-Ephl)/dEph+0.5)


    do iEph=1,numEph
      ! for different photon energy
      Eph=Ephl+(iEph-1.0)*dEph
      do iE=1,numE
        Ee_(iE)=iE*dE+Eph
        ve(iE)=sqrt(Ee_(iE)*keV_erg/const_me)
        !sigmaKr(iE)=sigma0/(Eph*Ee_(iE))
        temp=sqrt(1-Eph/Ee_(iE))
        !sigmaBH(iE)=sigmaKr(iE)*log((1+temp)/(1-temp))
      enddo
  
      {do ixO^DB=ixOmin^DB,ixOmax^DB\}
        Eth=const_kB*Te(ixO^D)/keV_erg
        Nnth=(Ne(ixO^D)/sqrt(dpi*Eth))*exp(-Enth**2/Eth**2)
        do iE=1,numE
          FE=Nnth*ve(iE)*(Ee_(iE)/Enth)**(-delta)
          !HXR_BR(ixO^D)=HXR_BR(ixO^D)+Ne(ixO^D)*FE*dE*sigmaBH(iE)*dEph
        enddo
      {enddo\}
    enddo

    !near earth
    HXR_BR=HXR_BR*depth   ! [cm^-3 s^-1] -> [cm^-2 s^-1]
    HXR_BR=HXR_BR*arcsec*arcsec   ! [cm^-2 s^-1] -> [arcsec^-2 s^-1]
    HXR_BR=HXR_BR/(4.0*dpi*AU*AU)


  end subroutine get_HXR_BR_old

  subroutine get_HXR_IC(HXR_IC,ixI^L,ixO^L,w,x)
  ! synthesize inverse Compton HXR based on equation (5) in Appendix A of 
  ! Krucker et al. (2008)
  ! blackbody distribution from equation (2.58) of Blumenthal & Gould (1970)
  ! The distribution of energetic electrons is N(Ee) ~ N0 (Ee/E0)^{-delta}
    use mod_global_parameters
    use mod_constants

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw)

    integer :: ix^D,ixO^D

    double precision :: eph,r0,mec2,E0,N0,const,Q,delta,p
    double precision :: Tph,Fsun,Rs,Ibd,F0
    double precision :: HXR_IC(ixI^S)
    double precision :: nmu,ei,dei,temp1,temp2,temp3,kT,ei_keV
    double precision :: pth(ixI^S),Te(ixI^S),Ne(ixI^S)
    integer :: nume,iei
    double precision :: Eth,Nnth,Enth


    eph=20.0        ! energy of HXR photons [keV]
    !N0=1.0d8        ! local fast electron density [cm^-3]
    Enth=5.0        ! above this energy the electrons are non-thermal [keV]
    E0=20.0         ! minimum value of fast electron energy [keV]
    p=2.5          ! fast electron flux spectral index
    delta=p+0.5     ! spectrum index of fast electrons
    mec2=511.0      ! energy of a static electron [keV]
    r0=2.8179d-13   ! radiu of electron [cm]
    Q=2*(11+4*delta+delta**2)/((1+delta) * (3+delta)**2 * (5+delta))

    Tph=6000.0      ! temperature of photosphere [K]
    Fsun=3.828d33   ! solar luminosity [erg s^-1]
    Rs=6.96e10      ! solar radius [cm]
    Ibd=1.22        ! integration of the blackbody distribution
    F0=Fsun/(4*dpi*Rs**2) ! energy flux near solar surface
    temp1=1.0/(Ibd * dpi**2 * (const_h*const_c/dpi)**3)

    kT=const_kB*Tph
    dei=0.1*const_eV  ! delta photon energy [erg]
    nume=100
    temp2=8*dpi*(r0**2)*const_c/eph
    temp2=temp2*(delta-1)*(E0/mec2)**(delta-1)*Q

    HXR_IC(ixO^S)=0.0

    ! output temperature
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    Te(ixO^S)=pth(ixO^S)/w(ixO^S,rho_)*unit_temperature
    Ne(ixO^S)=w(ixO^S,rho_)*unit_numberdensity

    {do ixO^DB=ixOmin^DB,ixOmax^DB\}
      Eth=const_kB*Te(ixO^D)/const_eV/1000  ! [keV]
      Nnth=(Ne(ixO^D)/sqrt(dpi*Eth))*exp(-Enth**2/Eth**2)
      N0=Nnth*(E0/Enth)**(-delta) * (E0/(delta-1))
      do iei=1,nume
        ei=iei*dei  ! [erg]
        nmu=F0*temp1*(ei**2 / (exp(ei/kT)-1))*dei/const_c
        ei_keV=ei/const_eV/1000   ! [erg] -> [keV]
        temp3=N0*temp2*nmu*(eph/(4*ei_keV))**((1-delta)/2)
        HXR_IC(ixO^D)=HXR_IC(ixO^D)+temp3
      enddo
    {enddo\}

  end subroutine get_HXR_IC

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
    double precision :: EUV(ixI^S),SXR(ixI^S)
    double precision :: HXR_BR(ixI^S),HXR_IC(ixI^S)
    integer :: idirmin,idir,ix^D

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
    !call divvector(Btotal,ixI^L,ixO^L,divb)
    !w(ixO^S,nw+3)=0.5d0*divb(ixO^S)/dsqrt(B2(ixO^S))/(^D&1.0d0/dxlevel(^D)+)
    !! output the plasma beta p*2/B**2
    !w(ixO^S,nw+4)=pth(ixO^S)*two/B2(ixO^S)
    ! output heating rate
    t=global_time
    !call getbQ(ens,ixI^L,ixO^L,t,w,x)
    !w(ixO^S,nw+3)=ens(ixO^S)
    call getlQ(lQgrid,ixI^L,ixO^L,t,w,x)
    w(ixO^S,nw+3)=lQgrid(ixO^S)
    ! store the cooling rate 
    !if(mhd_radiative_cooling)call getvar_cooling(ixI^L,ixO^L,w,x,ens)
    !w(ixO^S,nw+7)=ens(ixO^S)

    !! store current
    !call curlvector(Btotal,ixI^L,ixO^L,curlvec,idirmin,1,ndir)
    !do idir=1,ndir
    !  w(ixO^S,nw+7+idir)=curlvec(ixO^S,idir)
    !end do
 
    !call get_AIA(AIAgrid,ixI^L,ixO^L,w,x)
    !w(ixO^S,nw+5)=AIAgrid(ixO^S)

    call get_SXR(SXR,ixI^L,ixO^L,w,x,6,12)
    w(ixO^S,nw+4)=SXR(ixO^S)

    !call get_HXR_BR(HXR_BR,ixI^L,ixO^L,w,x,12,25)
    !w(ixO^S,nw+7)=HXR_BR(ixO^S)

    !call get_HXR_BR(HXR_BR,ixI^L,ixO^L,w,x,25,50)
    !w(ixO^S,nw+8)=HXR_BR(ixO^S)

    !call get_HXR_IC(HXR_IC,ixI^L,ixO^L,w,x)
    !w(ixO^S,nw+9)=HXR_IC(ixO^S)

    call get_HXR_BR(HXR_BR,ixI^L,ixO^L,w,x,25,50)
    w(ixO^S,nw+5)=HXR_BR(ixO^S)

    !call get_IRIS(EUV,ixI^L,ixO^L,w,x)
    !w(ixO^S,nw+2)=EUV(ixO^S)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='Te Alfv lQ SXR_6_12 HXR_25_50'
    !varnames='Te Alfv lQ SXR AIA131 HXR_6_12 HXR_12_25 HXR_25_50 HXR_IC'
    !varnames='Te Alfv divB beta bQ lQ rad j1 j2 j3 SXR AIA131 HXR_BR HXR_IC'

  end subroutine specialvarnames_output

  subroutine special_output(qunit)
    use mod_global_parameters

    integer, intent(in) :: qunit

    call XR_flux_output(qunit)

    !call output_selected_region(qunit)

    !call output_spectra(qunit)

  end subroutine special_output

  subroutine output_spectra(qunit)
    ! output HXR spectra
    use mod_usr_methods

    integer          :: ixO^L,ixO^D,j

    integer, intent(in) :: qunit
    character(len=30)   :: filename
    integer :: iigrid,igrid
    integer :: filenr
    logical :: fileopen

    integer :: iSp,numSp,dSp,wSp,Emax,Emin
    integer, allocatable :: Esp(:)
    double precision, allocatable :: Fapex(:),Ffoot(:)
    double precision, allocatable :: Fapex_recv(:),Ffoot_recv(:)

    dSp=5   ! distance between two spectral points
    wSp=1   ! band width of each point
    Emin=20   ! [keV]
    Emax=100  ! [keV]
    numSp=(Emax-Emin)/dSp+1

    allocate(Esp(numSp),Fapex(numSp),Ffoot(numSp))
    allocate(Fapex_recv(numSp),Ffoot_recv(numSp))
    do iSp=1,numSp
      Esp(iSp)=Emin+(iSp-1)*dSp
    enddo
    Fapex=0   ! apex flux
    Ffoot=0   ! footpoints flux
    Fapex_recv=0
    Ffoot_recv=0

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      block=>ps(igrid)
      call get_HXR_spectra(ixg^LL,ixm^LL,ps(igrid)%w,ps(igrid)%x, &
                         Esp,Fapex,Ffoot,numSp,wSp)
    enddo    


    call MPI_ALLREDUCE(Fapex,Fapex_recv,numSp,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,icomm,ierrmpi)
    call MPI_ALLREDUCE(Ffoot,Ffoot_recv,numSp,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,icomm,ierrmpi)

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
      write(1,*) numSp
      write(1,*) 'E apex footpoint'
      do iSp=1,numSp
        write(1,'(I5, e15.7, e15.7)') Esp(iSp),Fapex_recv(iSp),Ffoot_recv(iSp)
      enddo    
      close(1)
    endif

    deallocate(Esp,Fapex,Ffoot,Fapex_recv,Ffoot_recv)

  end subroutine output_spectra

  subroutine get_HXR_spectra(ixI^L,ixO^L,w,x,Esp,Fapex,Ffoot,numSp,wSp)
    ! output total XR flux
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw)
    integer :: ix^D,ixO^D
    double precision :: SXR(ixI^S),HXR(ixI^S)
    double precision :: dxb,dyb
    double precision :: arcsec,area

    integer :: iSp,numSp,wSp
    integer :: Esp(numSp)
    double precision :: Fapex(numSp),Ffoot(numSp)

    arcsec=7.35d-2   ! 735 km
    dxb=x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1)
    dyb=x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2)
    area=(dxb/arcsec)*(dyb/arcsec)

    do iSp=1,numSp
      call get_HXR_BR(HXR,ixI^L,ixO^L,w,x,Esp(iSp),Esp(iSp)+wSp)
      {do ixO^DB=ixOmin^DB,ixOmax^DB\}
        if (x(ixO^D,2)>=1.0) then
          Fapex(iSp)=Fapex(iSp)+HXR(ixO^D)*area
        else if (x(ixO^D,2)<0.5) then
          Ffoot(iSp)=Ffoot(iSp)+HXR(ixO^D)*area
        endif
      {enddo\}
    enddo

  end subroutine get_HXR_spectra

  subroutine output_selected_region(qunit)  
    ! output the average value of selected region
    use mod_usr_methods

    integer          :: ixO^L,ixO^D,j

    integer, intent(in) :: qunit
    character(len=30)   :: filename
    integer :: iigrid,igrid
    integer :: filenr
    logical :: fileopen

    integer :: iP,numP
    double precision, allocatable :: xP(:),yP(:)
    double precision, allocatable :: SXRP(:),NeP(:),TeP(:),pthP(:)
    double precision, allocatable :: SXRP_recv(:),NeP_recv(:),pthP_recv(:)
    double precision, allocatable :: areaP(:),areaP_recv(:)
    double precision :: sideL

    numP=5
    allocate(xP(numP),yP(numP),SXRP(numP))
    allocate(NeP(numP),TeP(numP),pthP(numP))
    allocate(SXRP_recv(numP),NeP_recv(numP),pthP_recv(numP))
    allocate(areaP(numP),areaP_recv(numP))

    sideL=0.1
    xP=(/-0.64, -0.12, 0.76,  -1.2, 0.06/)
    yP=(/2.2,   2.3,   2.13,  1.6,  1.91/)

    SXRP=0
    SXRP_recv=0
    NeP=0
    NeP_recv=0
    pthP=0
    pthP_recv=0
    areaP=0
    areaP_recv=0


    do iigrid=1,igridstail; igrid=igrids(iigrid);
      block=>ps(igrid)
      call get_local_value(ixg^LL,ixm^LL,ps(igrid)%w,ps(igrid)%x,sideL, &
                             xP,yP,SXRP,NeP,pthP,areaP,numP)
    enddo

    call MPI_ALLREDUCE(SXRP,SXRP_recv,numP,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,icomm,ierrmpi)
    call MPI_ALLREDUCE(NeP,NeP_recv,numP,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,icomm,ierrmpi)
    call MPI_ALLREDUCE(pthP,pthP_recv,numP,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,icomm,ierrmpi)
    call MPI_ALLREDUCE(areaP,areaP_recv,numP,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,icomm,ierrmpi)

    do iP=1,numP
      SXRP_recv(iP)=SXRP_recv(iP)/areaP_recv(iP)
      NeP_recv(iP)=NeP_recv(iP)/areaP_recv(iP)
      pthP_recv(iP)=pthP_recv(iP)/areaP_recv(iP)
      TeP(iP)=pthP_recv(iP)/NeP_recv(iP)
    enddo

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
      write(1,*) numP
      write(1,*) 'x y SXR Ne Te'
      do iP=1,numP
        write(1,'(e15.7, e15.7, e15.7, e15.7, e15.7)') &
              xP(iP),yP(iP),SXRP_recv(iP),NeP_recv(iP),TeP(iP)
      enddo
      close(1)
    endif

    deallocate(xP,yP,SXRP,NeP,TeP,pthP,areaP)
    deallocate(SXRP_recv,NeP_recv,pthP_recv,areaP_recv)

  end subroutine output_selected_region

  subroutine get_local_value(ixI^L,ixO^L,w,x,sideL,xP,yP,SXRP,NeP,pthP,areaP,numP)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw)
    integer :: ix^D,ixO^D

    double precision :: SXR(ixI^S)
    double precision :: dxb,dyb
    integer :: numP,iP
    double precision :: areaP(numP)
    double precision :: sideL
    double precision :: xP(numP),yP(numP),SXRP(numP)
    double precision :: NeP(numP),pthP(numP)
    double precision :: pth(ixI^S)

    dxb=x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1)
    dyb=x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2)

    call get_SXR(SXR,ixI^L,ixO^L,w,x,6,12)
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)

    {do ixO^DB=ixOmin^DB,ixOmax^DB\}
      do iP=1,numP
        if (abs(x(ixO^D,1)-xP(iP))<=sideL .and. abs(x(ixO^D,2)-yP(iP))<=sideL) then
          SXRP(iP)=SXRP(iP)+SXR(ixO^D)*dxb*dyb
          NeP(iP)=NeP(iP)+w(ixO^D,rho_)*dxb*dyb
          pthP(iP)=pthP(iP)+pth(ixO^D)*dxb*dyb
          areaP(iP)=areaP(iP)+dxb*dyb
        endif
      enddo
    {enddo\}

  end subroutine get_local_value

  subroutine XR_flux_output(qunit)
    ! integrating XR emission and output
    use mod_usr_methods

    integer          :: ixO^L,ixO^D,j

    integer, intent(in) :: qunit
    character(len=30)   :: filename
    integer :: iigrid,igrid
    double precision :: SXRf,HXRf,HXRfA,HXRfF
    double precision :: SXRf_recv,HXRf_recv,HXRfA_recv,HXRfF_recv
    double precision :: heating,heating_recv
    integer :: filenr
    logical :: fileopen

    SXRf=0    ! SXR flux
    SXRf_recv=0
    HXRf=0    ! HXR flux
    HXRf_recv=0
    HXRfA=0   ! HXR flux in apex
    HXRfA_recv=0
    HXRfF=0   ! HXR flux in footpoints
    HXRfF_recv=0 
    heating=0
    heating_recv=0


    do iigrid=1,igridstail; igrid=igrids(iigrid);
      block=>ps(igrid)
      call get_XR_intens(ixg^LL,ixm^LL,ps(igrid)%w,ps(igrid)%x, &
                         SXRf,HXRf,HXRfA,HXRfF,heating)
    enddo    


    call MPI_ALLREDUCE(SXRf,SXRf_recv,1,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,icomm,ierrmpi)
    call MPI_ALLREDUCE(HXRf,HXRf_recv,1,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,icomm,ierrmpi)
    call MPI_ALLREDUCE(HXRfA,HXRfA_recv,1,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,icomm,ierrmpi)
    call MPI_ALLREDUCE(HXRfF,HXRfF_recv,1,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,icomm,ierrmpi)
    call MPI_ALLREDUCE(heating,heating_recv,1,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,icomm,ierrmpi)

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
      write(1,*) 'SXR_6_12_keV  HXR_25_50_keV  HXR_apex  HXR_footpoints'
      write(1,'(e15.7, e15.7, e15.7, e15.7)') SXRf_recv,HXRf_recv,HXRfA_recv,HXRfF_recv
      write(1,*) 'heating_rate'
      write(1,'(e15.7)') heating_recv
      close(1)
    endif

  end subroutine XR_flux_output

  subroutine get_XR_intens(ixI^L,ixO^L,w,x,SXRf,HXRf,HXRfA,HXRfF,lQ)
    ! output total XR flux
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw)
    integer :: ix^D,ixO^D
    double precision :: SXR(ixI^S),HXR(ixI^S),lQgrid(ixI^S)
    double precision :: dxb,dyb
    double precision :: SXRf,HXRf,HXRfA,HXRfF,lQ
    double precision :: arcsec,area

    arcsec=7.35d-2   ! 735 km
    dxb=x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1)
    dyb=x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2)
    area=(dxb/arcsec)*(dyb/arcsec)

    call get_SXR(SXR,ixI^L,ixO^L,w,x,6,12)
    call get_HXR_BR(HXR,ixI^L,ixO^L,w,x,25,50)
    call getlQ(lQgrid,ixI^L,ixO^L,global_time,w,x)

    {do ixO^DB=ixOmin^DB,ixOmax^DB\}
      SXRf=SXRf+SXR(ixO^D)*area
      HXRf=HXRf+HXR(ixO^D)*area
      lQ=lQ+lQgrid(ixO^D)*area
      if (x(ixO^D,2)>=1.0) then
        HXRfA=HXRfA+HXR(ixO^D)*area
      else if (x(ixO^D,2)<0.5) then
        HXRfF=HXRfF+HXR(ixO^D)*area
      endif
    {enddo\}

  end subroutine get_XR_intens

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a steady (time-independent) potential or 
  ! linear force-free background field
    use mod_global_parameters

    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    wB0(ixO^S,1)=-Busr*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)
    wB0(ixO^S,2)= Busr*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
    wB0(ixO^S,3)=-Busr*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)

  end subroutine specialset_B0

  subroutine specialset_J0(ixI^L,ixO^L,x,wJ0)
  ! Here add a time-independent background current density 
    use mod_global_parameters

    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wJ0(ixI^S,7-2*ndir:ndir)

    wJ0(ixO^S,1)= ly*Busr*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
    wJ0(ixO^S,2)=-kx*Busr*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
    wJ0(ixO^S,3)= kx*Busr*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))&
                 -ly*Busr*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)

  end subroutine specialset_J0

end module mod_usr
