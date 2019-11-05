module mod_usr
  use mod_mhd
  implicit none
  double precision :: q_e, parb,unit_currentdensity
  double precision :: unit_electricfield

  double precision, allocatable :: pbc(:),rbc(:)
  double precision :: usr_grav
  double precision :: heatunit,gzone,SRadius,bQ0,dya
  double precision, allocatable :: pa(:),ra(:),ya(:),Ta(:)
  double precision, allocatable :: bQa(:)
  integer, parameter :: jmax=20000

  integer :: numxQ^D,numFL,numLP
  double precision :: xQmin^D,xQmax^D,dxQ^D,dFh
  double precision, allocatable :: xQ(:^D&,:),Qe(:^D&)
  double precision, allocatable :: xFLb(:,:),xFRb(:,:)
  integer,allocatable :: numRL(:),numRR(:)
  integer,allocatable :: numTurnL(:),numTurnR(:)
  integer,allocatable :: BopenL(:),BopenR(:)
  double precision, allocatable :: EtotL(:),EtotR(:)
  double precision, allocatable :: muBL(:),muBR(:)
  double precision :: eta1,eta2,eta3,etam,tar,vc
  integer :: numValidL,numValidR,filenr

  integer :: numN
  double precision :: Nmin,Nmax,dNlog
  double precision, allocatable :: HXR(:^D&)

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
    unit_currentdensity=unit_magneticfield*const_c/unit_length/4.d0/dpi
    ! unit of electric field
    unit_electricfield=unit_magneticfield*unit_length/(unit_time*const_c)
    ! unit of charge
    q_e=unit_currentdensity/unit_numberdensity/unit_velocity
    if(mype==0) print*,'unit of charge',q_e
    ! dimensionless charge of electron
    q_e=1.60217653d-19/q_e
    if(mype==0) print*,'dimensionless e',q_e

    if (mype==0) then
      print *, 'rho', unit_density
      print *, 'B', unit_magneticfield
      print *, 'p', unit_pressure
      print *, 'v', unit_velocity
      print *, 't', unit_time
      print *, 'E', unit_electricfield
      print *, 'J', unit_currentdensity
      print *, 'eta', unit_electricfield/unit_currentdensity
    endif



  end subroutine usr_init

  subroutine initglobaldata_usr()
    use mod_global_parameters

    heatunit=unit_pressure/unit_time          ! 3.697693390805347E-003 erg*cm^-3/s

    usr_grav=-2.74d4*unit_length/unit_velocity**2 ! solar gravity
    !bQ0=1.d3/heatunit ! background heating power density
    bQ0=1.0d-2/heatunit
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
    rpho=0.9d15/unit_numberdensity ! number density at the bottom relaxla
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

    tar=0.4d0
    !eta1=1.d-1
    !vc=1.d3
    !eta2=1.d-3
    !etam=1.d0
    eta1=5.d-2
    eta3=1.d-2
    vc=1.d3
    eta2=1.d-3
    etam=1.d0


    refine_factor=2**(refine_max_level-1)
    dFh=(xprobmax2-xprobmin2)/(domain_nx2*refine_factor)

    dt_update_Qe=max(dFh/2.0,dt)

    xQmin1=-3.0d0
    xQmax1=3.0d0
    xQmin2=0.d0
    xQmax2=10.0d0
    dxQ1=dFh
    dxQ2=dFh
    numXQ1=floor((xQmax1-0.d0)/dxQ1)*2
    numXQ2=floor((xQmax2-xQmin2)/dxQ2)

    lengthFL=(xprobmax2-xprobmin2)*1.5d0
    numLP=floor(lengthFL/(dFh))
    numFL=numXQ1

    allocate(xFLb(numFL,ndim),xFRb(numFL,ndim))
    allocate(xQ(numXQ1,numXQ2,ndim),Qe(numXQ1,numXQ2))
    allocate(numRL(numFL),numRR(numFL))
    allocate(numTurnL(numFL),numTurnR(numFL))
    allocate(BopenL(numFL),BopenR(numFL))
    allocate(HXR(numXQ1,numXQ2))
    allocate(EtotL(numFL),EtotR(numFL))
    allocate(muBL(numFL),muBR(numFL))

    do ix1=1,numFL
      !xFLb(ix1,1)=0.d0-0.5d0*dxQ1-(2*ix1-1)*dxQ1
      !xFLb(ix1,1)=0.d0-0.5d0*dxQ1-(ix1-1)*dxQ1
      xFLb(ix1,1)=0.d0-0.5d0*dxQ1-(ix1-0.5)*dxQ1
      xFLb(ix1,2)=0.d0
      !xFRb(ix1,1)=0.d0+0.5d0*dxQ1+(2*ix1)*dxQ1
      !xFRb(ix1,1)=0.d0+0.5d0*dxQ1+(ix1-1)*dxQ1
      xFRb(ix1,1)=0.d0+0.5d0*dxQ1+(ix1)*dxQ1
      xFRb(ix1,2)=0.d0
    enddo


    do ix1=1,numXQ1
      do ix2=1,numXQ2
        xQ(ix1,ix2,1)=0.d0+(ix1-0.5d0*numxQ1-0.5)*dxQ1
        xQ(ix1,ix2,2)=0.d0+(ix2-0.5)*dxQ2
      enddo
    enddo

    xQmin1=xQ(1,1,1)-0.5d0*dxQ1
    xQmax1=xQ(numXQ1,1,1)+0.5d0*dxQ1
    xQmin2=xQ(1,1,2)-0.5d0*dxQ2
    xQmax2=xQ(1,numXQ2,2)+0.5d0*dxQ2


    ! initial fast electron energy
    EtotL=0.d0
    EtotR=0.d0
    Qe=0.d0
    muBL=0.d0
    muBR=0.d0

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

    double precision :: lQgrid(ixI^S),bQgrid(ixI^S),cQgrid(ixI^S)

    ! add global background heating bQ
    call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)

    !! add localized external heating lQ concentrated at feet of loops
    if(iprob > 2)then
      !energy loss owing to particle acceleration
      cQgrid=0.d0
      if (qt>tar) then
        call getcQ(cQgrid,ixI^L,ixO^L,qtC,wCT,x)
        w(ixO^S,e_)=w(ixO^S,e_)-qdt*cQgrid(ixO^S)
      endif

      lQgrid=0.d0
      !heating because of fast electron deposition
      if (qt>tar) then
        call getlQ(lQgrid,ixI^L,ixO^L,qtC,wCT,x)
        w(ixO^S,e_)=w(ixO^S,e_)+qdt*lQgrid(ixO^S)
      endif
    endif

  end subroutine special_source

  subroutine getcQ(cQgrid,ixI^L,ixO^L,qt,w,x)
  ! calculate background heating bQ
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)
    double precision :: cQgrid(ixI^S)

    integer :: idir,idirmin,ix^D
    double precision :: current(ixI^S,3),jpara(ixI^S),j2(ixI^S)
    double precision :: Btotal(ixI^S,ndir),B2(ixI^S),Alfv(ixI^S),Ma(ixI^S)
    double precision :: ED(ixI^S),pth(ixI^S),Te(ixI^S),eta(ixI^S)
    double precision :: Bn(ixI^S),v(ixI^S,ndir),Efield(ixI^S,ndir),jE3(ixI^S)
    

    if(B0field) then
      Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))+block%B0(ixI^S,1:ndir,0)
    else
      Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))
    endif
    ! B^2
    !B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)
    B2(ixO^S)=Btotal(ixO^S,1)**2+Btotal(ixO^S,2)**2+Btotal(ixO^S,3)**2
    Bn(ixO^S)=sqrt(Btotal(ixO^S,1)**2/(Btotal(ixO^S,1)**2+Btotal(ixO^S,2)**2))

    cQgrid=0.d0

    call get_current(w,ixI^L,ixO^L,idirmin,current)
    call special_eta(w,ixI^L,ixO^L,idirmin,x,current,eta)

    jpara(ixO^S)=abs(current(ixO^S,1)*Btotal(ixO^S,1)+current(ixO^S,2)*Btotal(ixO^S,2)+&
                     current(ixO^S,3)*Btotal(ixO^S,3))/sqrt(B2(ixO^S))
    j2(ixO^S)=current(ixO^S,1)**2+current(ixO^S,2)**2+current(ixO^S,3)**2
    !cQgrid(ixO^S)=Bn(ixO^S)*eta(ixO^S)*j2(ixO^S)

    Alfv(ixO^S)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))
    Ma(ixO^S)=(w(ixO^S,mom(2))/w(ixO^S,rho_))/Alfv(ixO^S)


    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    Te(ixO^S)=pth(ixO^S)/w(ixO^S,rho_)


    !do idir=1,ndir
    !  v(ixO^S,idir)=w(ixO^S,mom(idir))/w(ixO^S,rho_)
    !enddo
    !call cross_product(ixI^L,ixO^L,Btotal,v,Efield)
    !do idir=1,ndir
    !  Efield(ixO^S,idir)=Efield(ixO^S,idir)+eta(ixO^S)*current(ixO^S,idir)
    !enddo
    !jE3(ixO^S)=current(ixO^S,3)*Efield(ixO^S,3)
    !{do ix^DB=ixOmin^DB,ixOmax^DB\}
    !  if (eta(ix^D)>eta3) eta(ix^D)=eta3
    !  if (Te(ix^D)>5) cQgrid(ix^D)=eta(ix^D)*j2(ix^D)
    !{end do\}



    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      Bn(ix^D)=abs(Btotal(ix^D,1))/10
      if (Bn(ix^D)>1) Bn(ix^D)=1
      if (Te(ix^D)>5) cQgrid(ix^D)=Bn(ix^D)*eta(ix^D)*j2(ix^D)
    {end do\}


    !{do ix^DB=ixOmin^DB,ixOmax^DB\}
    !  if (Bn(ix^D)>0.1) Bn(ix^D)=0.1
    !  if (Te(ix^D)>5) cQgrid(ix^D)=Bn(ix^D)*eta(ix^D)*j2(ix^D)
    !{end do\}

    !{do ix^DB=ixOmin^DB,ixOmax^DB\}
    !  if (Te(ix^D)>5 .and. Bn(ix^D)>0.01) cQgrid(ix^D)=eta(ix^D)*j2(ix^D)
    !{end do\}

  end subroutine getcQ

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

    !bQgrid=0.d0
   
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

  end subroutine get_bQa

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    use mod_global_parameters

    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pth(ixI^S),tmp(ixI^S),ggrid(ixI^S),invT(ixI^S)
    double precision :: delydelx, Bf(ixI^S,1:ndir)
    integer :: ix^D,idir,ixInt^L,ixIntB^L
    integer :: ixA^L

    select case(iB)
    case(1)
      ixA^L=ixO^L;
      ixAmin1=ixOmax1+1;ixAmax1=ixOmax1+nghostcells;
      call mhd_get_pthermal(w,x,ixI^L,ixA^L,pth)
      !w(ixO^S,rho_)=w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,rho_)
      !w(ixO^S,p_)=pth(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2)
      do ix1=ixOmin1,ixOmax1
        w(ix1^%1ixO^S,mom(1))=w(ixOmax1+1^%1ixO^S,mom(1))/w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(2))=w(ixOmax1+1^%1ixO^S,mom(2))/w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(3))=w(ixOmax1+1^%1ixO^S,mom(3))/w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,rho_)=w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,p_)=pth(ixOmax1+1^%1ixO^S)
      enddo
      do ix1=ixOmax1,ixOmin1,-1
        w(ix1,ixOmin2:ixOmax2,mag(:))=(1.0d0/3.0d0)* &
                   (-w(ix1+2,ixOmin2:ixOmax2,mag(:)) &
              +4.0d0*w(ix1+1,ixOmin2:ixOmax2,mag(:)))
      !  w(ix1,ixOmin2:ixOmax2,mom(:))=(1.0d0/3.0d0)* &
      !             (-w(ix1+2,ixOmin2:ixOmax2,mom(:)) &
      !        +4.0d0*w(ix1+1,ixOmin2:ixOmax2,mom(:)))
      enddo
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(2)
      ixA^L=ixO^L;
      ixAmin1=ixOmin1-nghostcells;ixAmax1=ixOmin1-1;
      call mhd_get_pthermal(w,x,ixI^L,ixA^L,pth)
      !w(ixO^S,rho_)=w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,rho_)
      !w(ixO^S,p_)=pth(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2)
      do ix1=ixOmin1,ixOmax1
        w(ix1^%1ixO^S,mom(1))=w(ixOmin1-1^%1ixO^S,mom(1))/w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(2))=w(ixOmin1-1^%1ixO^S,mom(2))/w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(3))=w(ixOmin1-1^%1ixO^S,mom(3))/w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,rho_)=w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,p_)=pth(ixOmin1-1^%1ixO^S)
      enddo
      do ix1=ixOmin1,ixOmax1
        w(ix1,ixOmin2:ixOmax2,mag(:))=(1.0d0/3.0d0)* &
                   (-w(ix1-2,ixOmin2:ixOmax2,mag(:)) &
              +4.0d0*w(ix1-1,ixOmin2:ixOmax2,mag(:)))
        !w(ix1,ixOmin2:ixOmax2,mom(:))=(1.0d0/3.0d0)* &
        !           (-w(ix1-2,ixOmin2:ixOmax2,mom(:)) &
        !      +4.0d0*w(ix1-1,ixOmin2:ixOmax2,mom(:)))
      enddo
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(3)
      ixA^L=ixO^L;
      ixAmin2=ixOmax2+1;ixAmax2=ixOmax2+1;
      call mhd_get_pthermal(w,x,ixI^L,ixA^L,pth)
      do ix2=ixOmin2,ixOmax2
        w(ix2^%2ixO^S,rho_)=w(ixOmax2+1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(1))=w(ixOmax2+1^%2ixO^S,mom(1))/w(ixOmax2+1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mag(2))=w(ixOmax2+1^%2ixO^S,mag(2))
        w(ix2^%2ixO^S,mag(3))=w(ixOmax2+1^%2ixO^S,mag(3))
        w(ix2^%2ixO^S,p_)=pth(ixOmax2+1^%2ixO^S)
        if(mhd_glm) w(ix2^%2ixO^S,psi_)=w(ixOmax2+1^%2ixO^S,psi_)
      enddo
      w(ixO^S,mom(2:3))=zero
      w(ixO^S,mag(1))=zero
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      !!ixA^L=ixO^L;
      !!ixAmin2=ixOmin2-1;ixAmax2=ixOmin2-1;
      !!call mhd_get_pthermal(w,x,ixI^L,ixA^L,pth)
      !!do ix2=ixOmin2,ixOmax2
      !!  w(ix2^%2ixO^S,rho_)=w(ixOmin2-1^%2ixO^S,rho_)
      !!  w(ix2^%2ixO^S,mom(1))=w(ixOmin2-1^%2ixO^S,mom(1))/w(ixOmin2-1^%2ixO^S,rho_)
      !!  w(ix2^%2ixO^S,mom(2))=w(ixOmin2-1^%2ixO^S,mom(2))/w(ixOmin2-1^%2ixO^S,rho_)
      !!  w(ix2^%2ixO^S,mom(3))=w(ixOmin2-1^%2ixO^S,mom(3))/w(ixOmin2-1^%2ixO^S,rho_)
      !!  w(ix2^%2ixO^S,p_)=pth(ixOmin2-1^%2ixO^S)
      !!enddo
      !!do ix2=ixOmin2,ixOmax2
      !!  w(ixOmin1:ixOmax1,ix2,mag(:))=(1.0d0/3.0d0)* &
      !!              (-w(ixOmin1:ixOmax1,ix2-2,mag(:))&
      !!         +4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(:)))
      !!enddo

      !ixInt^L=ixO^L;
      !ixIntmin2=ixOmin2-1;ixIntmax2=ixOmin2-1;
      !call mhd_get_pthermal(w,x,ixI^L,ixInt^L,pth)
      !ixIntmin2=ixOmin2-1;ixIntmax2=ixOmax2;
      !call getggrav(ggrid,ixI^L,ixInt^L,x)
      !! fill pth, rho ghost layers according to gravity stratification
      !invT(ixOmin2-1^%2ixO^S)=w(ixOmin2-1^%2ixO^S,rho_)/pth(ixOmin2-1^%2ixO^S)
      !tmp=0.d0
      !do ix2=ixOmin2,ixOmax2
      !  tmp(ixOmin2-1^%2ixO^S)=tmp(ixOmin2-1^%2ixO^S)+0.5d0*&
      !      (ggrid(ix2^%2ixO^S)+ggrid(ix2-1^%2ixO^S))*invT(ixOmin2-1^%2ixO^S)
      !  w(ix2^%2ixO^S,p_)=pth(ixOmin2-1^%2ixO^S)*dexp(tmp(ixOmin2-1^%2ixO^S)*dxlevel(2))
      !  w(ix2^%2ixO^S,rho_)=w(ix2^%2ixO^S,p_)*invT(ixOmin2-1^%2ixO^S)
      !enddo
      !! fixed zero velocity
      !do idir=1,ndir
      !  w(ixO^S,mom(idir))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,mom(idir))&
      !               /w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,rho_)
      !end do
      !! zero normal gradient extrapolation
      !do ix2=ixOmin2,ixOmax2
      !  w(ixOmin1:ixOmax1,ix2,mag(:))=(1.0d0/3.0d0)* &
      !              (-w(ixOmin1:ixOmax1,ix2-2,mag(:))&
      !         +4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(:)))
      !enddo

      ixA^L=ixO^L;
      ixAmin2=ixOmin2-1;ixAmax2=ixOmin2-1;
      call mhd_get_pthermal(w,x,ixI^L,ixA^L,pth)
      do ix2=ixOmin2,ixOmax2
        w(ix2^%2ixO^S,rho_)=w(ixOmin2-1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,p_)=pth(ixOmin2-1^%2ixO^S)
        w(ix2^%2ixO^S,mom(1))=w(ixOmin2-1^%2ixO^S,mom(1))/w(ixOmin2-1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(2))=0
        w(ix2^%2ixO^S,mom(3))=0
        w(ix2^%2ixO^S,mag(1))=0
        w(ix2^%2ixO^S,mag(2))=w(ixOmin2-1^%2ixO^S,mag(2))
        w(ix2^%2ixO^S,mag(3))=w(ixOmin2-1^%2ixO^S,mag(3))
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
    character (30) :: fname,tempc
    integer :: ixFL,iyLH,ix1,ix2
    double precision :: t1,t2,t_output
    integer :: ifile,tempi

    ! for converting data
    if (iprob==3 .and. convert .and. qt>tar) then
      dt_update_Qe=dFh/5.0

      ! read Etot
      if (mype==0) then
        filenr=snapshotini
        write(fname, '(a,i4.4,a)') trim(base_filename),filenr,'_Etot.txt'
        open(1,file=fname,action='READ')
        read(1,*) tempc
        read(1,*) tempi
        read(1,*) tempc,tempc,tempc,tempc,tempc,tempc,tempc
        do ix1=1,numFL
          read(1,*) tempi,xFLb(ix1,1),xFRb(ix1,1),EtotL(ix1),EtotR(ix1),&
                    muBL(ix1),muBR(ix1)
        enddo
        close(1)
      endif
      call MPI_BCAST(xFLb(:,1),numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(xFRb(:,1),numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(EtotL,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(EtotR,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(muBL,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(muBR,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)

      !call get_flare_eflux()
    endif


    ! for restart
    if (iprob==3 .and. restart_from_file /= undefined)  then
      dt_update_Qe=dFh/5.0

      ! read Etot
      if (mype==0) then
        filenr=snapshotini
        write(fname, '(a,i4.4,a)') trim(base_filename),filenr,'_Etot.txt'
        open(1,file=fname,action='READ')
        read(1,*) tempc
        read(1,*) tempi
        read(1,*) tempc,tempc,tempc,tempc,tempc,tempc,tempc
        do ix1=1,numFL
          read(1,*) tempi,xFLb(ix1,1),xFRb(ix1,1),EtotL(ix1),EtotR(ix1),&
                    muBL(ix1),muBR(ix1)
        enddo
        close(1)
      endif
      call MPI_BCAST(xFLb(:,1),numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(xFRb(:,1),numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(EtotL,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(EtotR,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(muBL,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(muBR,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)

      call MPI_BCAST(filenr,1,MPI_INTEGER,0,icomm,ierrmpi)
      filenr=filenr+1
    endif


    ! calculate heating
    if (iprob==3 .and. qt>tar) then

      dt_update_Qe=max(dFh/5.0,dt)

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


    ! output Etot for each fieldline
    if (convert .eqv. .false.) then 
      if (iit==0) then
        filenr=0
      endif
      t_output=filenr*dtsave(2)

      if (qt>=t_output .and. qt<t_output+dt) then     
        if (mype==0) then
          write(fname, '(a,i4.4,a)') trim(base_filename),filenr,'_Etot.txt'
          open(1,file=fname)
          write(1,*) 'numFL'
          write(1,*) numFL
          write(1,*) 'iFL xFL(1) xFR(1) EtotL EtotR muBL muBR'
          do ix1=1,numFL
            write(1,'(i8, e15.7, e15.7, e15.7, e15.7, e15.7, e15.7)'), ix1,&
                  xFLb(ix1,1),xFRb(ix1,1),EtotL(ix1),EtotR(ix1),&
                  muBL(ix1),muBR(ix1)
          enddo
          close(1)
        endif
        filenr=filenr+1
      endif

    endif

  end subroutine special_global

  subroutine get_flare_eflux()
    ! calculate electron deposit energy flux
    use mod_global_parameters

    integer :: ix^D,j
    double precision :: xFL(numFL,numLP,ndim),xFR(numFL,numLP,ndim)
    double precision :: QeL(numFL,numLP),QeR(numFL,numLP)
    double precision :: wBL(numFL,numLP,nw+ndir),wBR(numFL,numLP,nw+ndir)
    double precision :: eFluxL(numFL),eFluxR(numFL)
    double precision :: HXRL(numFL,numLP),HXRR(numFL,numLP)

    double precision, allocatable :: QeN(:)

    double precision :: delta,Ec,dEe,keV_erg
    integer :: numEe,iEe
    double precision, allocatable :: Ee(:),spectra(:,:)
    character(20) :: fname

    ! initial fast electron spectral parameters
    delta=4.0
    Ec=20.0  ! cutoff energy [keV]
    dEe=1.0  
    numEe=280
    keV_erg=1.0d3*const_ev  ! 1 keV = * erg

    xFL=0
    xFR=0
    wBL=0
    wBR=0


    xFL(:,1,1)=xFLb(:,1)
    xFL(:,1,2)=xFLb(:,2)
    xFR(:,1,1)=xFRb(:,1)
    xFR(:,1,2)=xFRb(:,2)


    ! trace Bfield to prapare for heating
    call update_Bfield(xFL,xFR,wBL,wBR)


    ! find the heating region 
    call locate_heating(xFL,xFR)


    ! update heating table
    numN=1001
    allocate(QeN(numN))
    Nmin=1.0e18
    Nmax=1.0e23
    dNlog=log10(Nmax/Nmin)/(numN-1.0) 
    call update_heating_table(QeN,delta,Ec)
 
 
    ! get heating rate for each field line
    call get_heating_rate(xFL,wBL,QeL,QeN,numRL,numTurnL,numValidL,EtotL,muBL,eFluxL)
    call get_heating_rate(xFR,wBR,QeR,QeN,numRR,numTurnR,numValidR,EtotR,muBR,eFluxR)
    ! total number flux of fast electrons
    eFluxL=eFluxL*heatunit*unit_length*((delta-2)/(delta-1))/(Ec*keV_erg)
    eFluxR=eFluxR*heatunit*unit_length*((delta-2)/(delta-1))/(Ec*keV_erg)

  
    ! get local heating rate via interpolation
    call get_Qe(xFL,xFR,wBL,wBR,QeL,QeR)


    allocate(Ee(numEe),spectra(numN,numEe))


    if (convert) then
      ! get fast electron spectra table
      call get_spectra(delta,Ec,dEe,Ee,spectra,numEe)
  
      ! get HXR flux for each field line
      call get_HXR_line(Ee,spectra,numEe,dEe,xFL,wBL,numValidL,numTurnL,numRL,eFluxL,muBL,HXRL)
      call get_HXR_line(Ee,spectra,numEe,dEe,xFR,wBR,numValidR,numTurnR,numRR,eFluxR,muBR,HXRR)
  
      ! get local HXR flux via interpolation
      call interp_HXR(xFL,xFR,wBL,wBR,HXRL,HXRR)
  
    else 
      HXR=0.d0
    endif


    ! field line split
    call split_Bfield(xFL,xFLb,wBL,EtotL,numRL,numValidL)
    call split_Bfield(xFR,xFRb,wBR,EtotR,numRR,numValidR)


    deallocate(Ee,spectra)
    deallocate(QeN)

  end subroutine get_flare_eflux

  subroutine split_Bfield(xF,xFb,wB,Etot,numR,numValid)
    !
    use mod_global_parameters

    double precision :: xF(numFL,numLP,ndim),xFb(numFL,2)
    double precision :: Etot(numFL),wB(numFL,numLP,nw+ndir)
    integer :: numValid
    integer :: numR(numFL)
    
    integer :: ix^D,iFL
    double precision :: dx0,dxp,dxMax,widthMax,widthMin,B0,Bp
    double precision :: xFnew,Enew,dxL0,dxR0,dEtot
    logical :: splitB,mergeB


    widthMax=4.0*dFh
    widthMin=1.0*dFh

    ix1=3
    do while (ix1<numValid)
      dx0=abs(xF(ix1,1,1)-xF(ix1-1,1,1))
      dxMax=dx0
      B0=sqrt(wB(ix1,1,mag(1))**2+wB(ix1,1,mag(2))**2)

      splitB=.FALSE.
      mergeB=.FALSE.

      do ix2=1,numR(ix1)
        Bp=sqrt(wB(ix1,ix2,mag(1))**2+wB(ix1,ix2,mag(2))**2)
        if (Bp<1.d0) Bp=1.d0
        dxp=B0*dx0/Bp
        if (dxp>dxMax) dxMax=dxp
      enddo

      if (dxMax>widthMax) splitB=.TRUE.
      if (dxMax<widthMin) mergeB=.TRUE.


      ! for Bfield split
      if (splitB) then
        xFnew=0.5d0*(xF(ix1,1,1)+xF(ix1-1,1,1))
        Enew=0
        ! reduce energy from old field line and add
        ! energy to new field line
        dxL0=abs(xF(ix1+1,1,1)-xF(ix1,1,1))
        dEtot=Etot(ix1)*(0.25*dx0)/(0.5*(dxL0+dx0))
        Etot(ix1)=Etot(ix1)-dEtot
        Enew=Enew+dEtot
        dxR0=abs(xF(ix1-1,1,1)-xF(ix1-2,1,1))
        dEtot=Etot(ix1-1)*(0.25*dx0)/(0.5*(dxR0+dx0))
        Etot(ix1-1)=Etot(ix1-1)-dEtot
        Enew=Enew+dEtot

        do iFL=numFL,ix1+1,-1
          xF(iFL,1,1)=xF(iFL-1,1,1)
          xF(iFL,1,2)=0.d0
          Etot(iFL)=Etot(iFL-1)
        enddo
        xF(ix1,1,1)=xFnew
        Etot(ix1)=Enew
        if (numValid<numFL) numValid=numValid+1
      endif

      ! for Bfield merge
      if (mergeB) then
        ! add energy to nearby field lines
        dxL0=abs(xF(ix1+1,1,1)-xF(ix1,1,1))
        dEtot=Etot(ix1)*(0.5*dx0)/(0.5*dx0+0.5*dxL0)
        Etot(ix1+1)=Etot(ix1+1)+dEtot
        dEtot=Etot(ix1)*(0.5*dxL0)/(0.5*dx0+0.5*dxL0)
        Etot(ix1-1)=Etot(ix1-1)+dEtot

        do iFL=ix1,numFL-1
          xF(iFL,1,1)=xF(iFL+1,1,1)
          xF(iFL,1,2)=0.d0
          Etot(iFL)=Etot(iFL+1)
        enddo
        xF(numFL,1,1)=2*xF(numFL-1,1,1)-xF(numFL-2,1,1)
        Etot(numFL)=0.d0
        if (numValid>ix1) numValid=numValid-1
      endif

      ix1=ix1+1
    enddo


    xFb(:,1)=xF(:,1,1)
    xFb(:,2)=0.d0


  end subroutine split_Bfield

  subroutine update_Bfield(xFL,xFR,wBL,wBR)
    ! trace B field and get column depth
    use mod_global_parameters
    use mod_usr_methods
    use mod_trace_Bfield

    double precision :: xFL(numFL,numLP,ndim),xFR(numFL,numLP,ndim)
    double precision :: wBL(numFL,numLP,nw+ndir),wBR(numFL,numLP,nw+ndir)

    integer :: ix^D,dirct,j
    logical :: forward,interp
    double precision :: xmin
    logical :: wTrans(nw+ndir)
    double precision :: t0,t1
    
    t0=MPI_wtime()

    numValidL=numFL
    numValidR=numFL

    wTrans=.TRUE.
    interp=.TRUE.
    !wTrans(rho_)=.TRUE.
    !wTrans(mag(1))=.TRUE.
    !wTrans(mag(2))=.TRUE.
    !wTrans(mag(3))=.TRUE.


    ! left foot
    forward=.true.
    dirct=0
    TRACE1: do ix1=1,numFL

      call trace_Bfield(xFL(ix1,:,:),wBL(ix1,:,:),dFh,numLP,numRL(ix1),&
                        forward,wTrans,interp,ix1)

      !if field line is outside of the reconnection region, stop tracing
      xmin=abs(xFL(ix1,1,1))
      do ix2=1,numRL(ix1)
        if (abs(xFL(ix^D,1))<xmin) xmin=abs(xFL(ix^D,1))
      enddo
      if (xmin>0.5) then
        numValidL=ix1
        exit TRACE1
      endif
    enddo TRACE1


    !right foot
    forward=.false.
    dirct=0
    TRACE2: do ix1=1,numFL
      call trace_Bfield(xFR(ix1,:,:),wBR(ix1,:,:),dFh,numLP,numRR(ix1),&
                        forward,wTrans,interp,ix1)

      !if field line is outside of the reconnection region, stop tracing
      xmin=abs(xFR(ix1,1,1))
      do ix2=1,numRR(ix1)
        if (abs(xFR(ix^D,1))<xmin) xmin=abs(xFR(ix^D,1))
      enddo
      if (xmin>0.5) then
        numValidR=ix1
        exit TRACE2
      endif
    enddo TRACE2

    t1=MPI_wtime()


  end subroutine update_Bfield

  subroutine locate_heating(xFL,xFR)
    ! find the point that is most closed to x=0 for each field line
    ! fast electrons has a given pitch angle mu0 at this point
    use mod_global_parameters
    use mod_usr_methods

    double precision :: xFL(numFL,numLP,ndim),xFR(numFL,numLP,ndim)

    integer :: iFL,iLP,numRT
    double precision :: minx,absx

    numTurnL=1
    numTurnR=1

    ! for B field line start from left part
    do iFL=1,numValidL
      minx=abs(xFL(iFL,1,1))
      do iLP=1,numRL(iFL)
        absx=abs(xFL(iFL,iLP,1))
        if (absx<minx) then 
          minx=absx
          numTurnL(iFL)=iLP
        endif
      enddo
      if (numTurnL(iFL)<=2 .or. numTurnL(iFL)>=numRL(iFL)-1) then
        numTurnL(iFL)=(1+numRL(iFL))/2
      endif
    enddo


    ! for B field line start from right part
    do iFL=1,numValidR
      minx=abs(xFR(iFL,1,1))
      do iLP=1,numRR(iFL)
        absx=abs(xFR(iFL,iLP,1))
        if (absx<minx) then 
          minx=absx
          numTurnR(iFL)=iLP
        endif
      enddo
      if (numTurnR(iFL)<=2 .or. numTurnR(iFL)>=numRR(iFL)-1) then
        numTurnR(iFL)=(1+numRR(iFL))/2
      endif
    enddo


  end subroutine locate_heating

  subroutine update_heating_table(QeN,delta,Ec)
    ! calculate heating table based on shock compression ratio
    use mod_global_parameters
    use mod_usr_methods

    double precision :: QeN(numN),Ncol(numN)
    double precision :: delta,Ec   

    integer :: j,iBx,numBx,iBc,iNcol
    double precision :: Atom,temp
    double precision :: Cllog,gammaN,K,Nc,Nb,Efactor
    double precision :: dtx,tx,au,al,b,maxh
    double precision, allocatable :: Bxc(:)
    double precision :: Ec_erg


    ! parameters for electron deposition
    Cllog=25.0  ! Coulomb logarithns
    gammaN=Cllog  ! full ionized plasma
    K=2.0*dpi*const_e**4
    Ec_erg=Ec*1.0e3*const_ev
    Nc=Ec_erg**2/(2.0*gammaN*K)
    Efactor=K*gammaN*(delta-2.0)/(2*Ec_erg**2)


    ! beta function
    numBx=100001
    allocate(Bxc(numBx))
    dtx=1.0/(numBx-1)
    Bxc=0
    al=delta/2.0
    b=1.0/2
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

  end subroutine update_heating_table

  subroutine get_heating_rate(xF,wBLR,QeLR,QeN,numR,numTurn,numValid,Etot,muB,eFlux)
    ! calculate local heating rate
    use mod_global_parameters
    use mod_usr_methods

    integer :: numR(numFL),numTurn(numFL)
    double precision :: xF(numFL,numLP,ndim),QeLR(numFL,numLP)
    double precision :: wBLR(numFL,numLP,nw+ndir)
    double precision :: QeN(numN),Etot(numFL),muB(numFL),eFlux(numFL)
    integer :: numValid

    double precision :: Bv(numFL,numLP,ndim),Np(numFL,numLP)
    double precision :: Ncol(numLP,2)

    integer :: ix^D,iNlog,iFL,j
    double precision :: dl^D,dl,Bratio,Bxy0,BxyTurn,Bxy,mu0,mu,const

    logical :: Ereturn,ForReturn,BackReturn
    integer :: ireturn,iFreturn,iBreturn,iLPmin,iLPmax
    double precision :: ve,length,tau_e,dEdt

    double precision :: width,width0,area,el
    double precision :: jpara,Te,ED,B2,JdotB,Alfv,Ma,J2,Bn
    double precision :: E3,v1,v2,jE3

    double precision :: eta,heta,reta,heta2,rad,rad2,vd,muMin
    double precision :: mup,muT,muBsum,Esum


    mup=0.9
    muMin=0.1
    ve=1.d10/unit_velocity

    ! density and magnetic field
    Np(:,:)=wBLR(:,:,rho_)*unit_numberdensity
    do j=1,ndim
      Bv(:,:,j)=wBLR(:,:,mag(j))
    enddo


    ! calculate fast electron energy flux for each field line
    eFlux=0
    QeLR=0
    do ix1=2,numValid-1
      muBsum=0.d0
      Esum=0.d0

      ! initial width of the flux tube element
      Bxy0=dsqrt(Bv(ix1,1,1)**2+Bv(ix1,1,2)**2)
      width0=abs(xF(ix1+1,1,1)-xF(ix1-1,1,1))/2.d0
      BxyTurn=dsqrt(Bv(ix1,numTurn(ix1),1)**2+Bv(ix1,numTurn(ix1),2)**2)

      ! calculate particle acceleration rate
      do ix2=1,numR(ix1)
        Te=wBLR(ix^D,p_)/wBLR(ix^D,rho_)

        if (Te>5) then
          ! cosine pitch angle
          Bxy=dsqrt(Bv(ix1,ix2,1)**2+Bv(ix1,ix2,2)**2)
          const=(1-mup**2)/Bxy

          ! scan area
          dl1=xF(ix1,ix2,1)-xF(ix1,ix2-1,1)
          dl2=xF(ix1,ix2,2)-xF(ix1,ix2-1,2)
          dl=dsqrt(dl1**2+dl2**2)
          width=width0*Bxy0/Bxy
          area=dl*width

          ! local electron acceleration rate
          J2=0
          do j=1,ndir
            J2=J2+wBLR(ix^D,nw+j)**2
          enddo
          if (global_time<tar) then
            heta = 5.
            heta2 = 1.
            reta = 0.8d0 * 0.3d0
            rad=dsqrt(xF(ix^D,1)**2+(xF(ix^D,2)-heta)**2)
            rad2=exp(-((xF(ix^D,2)-heta)**2)/heta2**2)
            if (rad .lt. reta) then
              eta=eta1*(2.d0*(rad/reta)**3-3.d0*(rad/reta)**2+1.d0)*rad2
            else
              eta=zero
            endif
          else
            vd=sqrt(J2)/wBLR(ix^D,rho_)/q_e
            if (vd>vc) eta=eta2*(vd/vc-1)
            if (eta>etam) eta=etam
          end if 

          !Te=wBLR(ix^D,p_)/wBLR(ix^D,rho_)
          !if (eta>eta3) eta=eta3
          !if (Te>5) el=eta*J2
          !Esum=Esum+el*area*dt_update_Qe
          !muBsum=muBsum+el*area*dt_update_Qe*const
          

          ! energy got from local current/electric field
          Alfv=sqrt(wBLR(ix^D,mag(1))**2+wBLR(ix^D,mag(2))**2+&
                  wBLR(ix^D,mag(3))**2)/wBLR(ix^D,rho_)
          Ma=wBLR(ix^D,mom(2))/wBLR(ix^D,rho_)
          Te=wBLR(ix^D,p_)/wBLR(ix^D,rho_)
          B2=0
          do j=1,ndir
            B2=B2+wBLR(ix^D,mag(j))**2
          enddo
          Bn=abs(wBLR(ix^D,mag(1))/10)
          if (Bn>1) Bn=1
          if (Te>5) el=Bn*eta*J2
          Esum=Esum+el*area*dt_update_Qe
          muBsum=muBsum+el*area*dt_update_Qe*const
        endif
      enddo 

      if (Etot(ix1)+Esum>0) then
        muB(ix1)=(Etot(ix1)*muB(ix1)+muBsum)/(Etot(ix1)+Esum)
        Etot(ix1)=Etot(ix1)+Esum
      else
        Etot(ix1)=0
        muB(ix1)=0
      endif


      ! particle trapping
      length=0
      tau_e=0
      ForReturn=.false.
      BackReturn=.false.
      iFreturn=numR(ix1)
      iBreturn=1
      if (muB(ix1)*BxyTurn<1) then
        muT=sqrt(1-muB(ix1)*BxyTurn)
      else
        muT=0
      endif

      ! for the field lines have fast electron, calculate
      ! trapping length
      if (muT>0 .and. Etot(ix1)>0) then
        forward: do ix2=numTurn(ix1),numR(ix1)
          ! pitch angle
          Bxy=dsqrt(Bv(ix1,ix2,1)**2+Bv(ix1,ix2,2)**2)
          const=muB(ix1)
          ! return of electrons owing to loop expension
          if (const*Bxy>=1) then
            ForReturn=.true.
            iFreturn=ix2-1
            exit forward
          else
            mu=sqrt(1-const*Bxy)
          endif

          ! scan area
          dl1=xF(ix1,ix2,1)-xF(ix1,ix2-1,1)
          dl2=xF(ix1,ix2,2)-xF(ix1,ix2-1,2)
          dl=dsqrt(dl1**2+dl2**2)
          length=length+dl
        enddo forward

        backward: do ix2=numTurn(ix1),1,-1
          ! pitch angle
          Bxy=dsqrt(Bv(ix1,ix2,1)**2+Bv(ix1,ix2,2)**2)
          const=muB(ix1)
          ! return of electrons owing to loop expension
          if (const*Bxy>=1) then
            BackReturn=.true.
            iBreturn=ix2+1
            exit backward
          else
            mu=sqrt(1-const*Bxy)
          endif

          ! scan area
          dl1=xF(ix1,ix2+1,1)-xF(ix1,ix2,1)
          dl2=xF(ix1,ix2+1,2)-xF(ix1,ix2,2)
          dl=dsqrt(dl1**2+dl2**2)
          length=length+dl
        enddo backward
      endif


      ! time scale for fast electron moving
      tau_e=length/ve
      if (tau_e<dt_update_Qe) tau_e=dt_update_Qe
      
      ! energy flux of fast electron for rach field line
      eFlux(ix1)=Etot(ix1)/width0/tau_e


      ! region that has heating
      iLPmin=1
      iLPmax=numR(ix1)
      if (ForReturn) iLPmax=iFreturn
      if (BackReturn) iLPmin=iBreturn


      !calculate heating for field lines have fast electron
      if (muT>0 .and. Etot(ix1)>0) then
        ! calculate heating for one side
        Ncol=0
        dEdt=0
        if (numTurn(ix1)<iLPmax-1) then
          do ix2=numTurn(ix1),iLPmax-1
            Bxy=dsqrt(Bv(ix1,ix2,1)**2+Bv(ix1,ix2,2)**2)
            Bratio=Bxy/Bxy0
            const=muB(ix1)
            mu=sqrt(1-const*Bxy)
            if (mu<muMin) mu=muMin

            ! column depth
            dl1=xF(ix1,ix2+1,1)-xF(ix1,ix2,1)
            dl2=xF(ix1,ix2+1,2)-xF(ix1,ix2,2)
            dl=dsqrt(dl1**2+dl2**2)*unit_length/mu
            Ncol(ix2,1)=Ncol(ix2-1,1)+dl*(Np(ix1,ix2)+Np(ix1,ix2-1))/2.0
            width=width0*Bxy0/Bxy

            ! local heating rate
            iNlog=int(log10(Ncol(ix2,1)/Nmin)/dNlog+0.5)+1
            if (iNlog<1) iNlog=1
            if (iNlog>numN) iNlog=numN
            QeLR(ix^D)=QeN(iNlog)*Np(ix^D)*unit_length*(eFlux(ix1)*Bxy0/Bxy)*(muT/mu)/mu
            dEdt=dEdt+QeLR(ix^D)*width*dl
            !QeLR(ix^D)=eFlux(ix1)
          enddo
        endif

        ! reduce total energy
        if (ForReturn .eqv. .false.)  dEdt=eFlux(ix1)*width0
        if (dEdt>eFlux(ix1)*width0) dEdt=eFlux(ix1)*width0
        Etot(ix1)=Etot(ix1)-dEdt*dt_update_Qe


        ! calculate heating for the other side
        Ncol=0
        dEdt=0
        if (numTurn(ix1)>iLPmin+1) then
          do ix2=numTurn(ix1)-1,iLPmin,-1
            Bxy=dsqrt(Bv(ix1,ix2,1)**2+Bv(ix1,ix2,2)**2)
            Bratio=Bxy/Bxy0
            const=muB(ix1)
            mu=sqrt(1-const*Bxy)
            if (mu<muMin) mu=muMin

            ! column depth
            dl1=xF(ix1,ix2+1,1)-xF(ix1,ix2,1)
            dl2=xF(ix1,ix2+1,2)-xF(ix1,ix2,2)
            dl=dsqrt(dl1**2+dl2**2)*unit_length/mu
            Ncol(ix2,1)=Ncol(ix2+1,1)+dl*(Np(ix1,ix2)+Np(ix1,ix2+1))/2.0
            width=width0*Bxy0/Bxy

            ! local heating rate
            iNlog=int(log10(Ncol(ix2,1)/Nmin)/dNlog+0.5)+1
            if (iNlog<1) iNlog=1
            if (iNlog>numN) iNlog=numN
            QeLR(ix^D)=QeN(iNlog)*Np(ix^D)*unit_length*(eFlux(ix1)*Bxy0/Bxy)*(muT/mu)/mu
            dEdt=dEdt+QeLR(ix^D)*width*dl
            !QeLR(ix^D)=eFlux(ix1)
          enddo
        endif

        ! reduce total energy
        if (BackReturn .eqv. .false.)  dEdt=eFlux(ix1)*width0
        if (dEdt>eFlux(ix1)*width0) dEdt=eFlux(ix1)*width0
        Etot(ix1)=Etot(ix1)-dEdt*dt_update_Qe
      endif

      if (Etot(ix1)<0) Etot(ix1)=0

    enddo



  end subroutine get_heating_rate

  subroutine get_Qe(xFL,xFR,wBL,wBR,QeL,QeR)
    ! get local heating rate via interpolation
    use mod_global_parameters
    use mod_usr_methods

    double precision :: xFL(numFL,numLP,ndim),xFR(numFL,numLP,ndim)
    double precision :: wBL(numFL,numLP,nw+ndir),wBR(numFL,numLP,nw+ndir)
    double precision :: QeL(numFL,numLP),QeR(numFL,numLP)
 
    double precision :: sumWeights(numxQ1,numxQ2)
    double precision :: dxMax1,dxMax2,weight,dl,Bratio,Bratio1,Bratio2
    integer :: iFL,iLP,ixQ^L,ixQ^D,weightIndex

    integer :: ix^D,j
    character(20) :: fname

    Qe=0
    sumWeights=0
    dxMax1=1*dxQ1
    dxMax2=1*dxQ2
    weightIndex=1


    ! for the left foot
    do iFL=1,numValidL
      do iLP=1,numRL(iFL)
        weightIndex=4
        Bratio=Busr/sqrt(wBL(iFL,iLP,mag(1))**2+wBL(iFL,iLP,mag(2))**2)
        if (Bratio<3) then
          dxMax1=ceiling(Bratio)*dxQ1
          dxMax2=ceiling(Bratio)*dxQ2
        else
          dxMax1=3*dxQ1
          dxMax2=3*dxQ2
        endif

        ixQmin1=floor((xFL(iFL,iLP,1)-dxMax1-xQmin1)/dxQ1)+1
        ixQmin2=floor((xFL(iFL,iLP,2)-dxMax2-xQmin2)/dxQ2)+1
        ixQmax1=floor((xFL(iFL,iLP,1)+dxMax1-xQmin1)/dxQ1)+1
        ixQmax2=floor((xFL(iFL,iLP,2)+dxMax2-xQmin2)/dxQ2)+1
        if (ixQmin1<1) ixQmin1=1
        if (ixQmin2<1) ixQmin2=1
        if (ixQmax1>numxQ1) ixQmax1=numxQ1
        if (ixQmax2>numxQ2) ixQmax2=numxQ2

        !#
        do ixQ1=ixQmin1,ixQmax1
          do ixQ2=ixQmin2,ixQmax2
            dl=sqrt((xQ(ixQ1,ixQ2,1)-xFL(iFL,iLP,1))**2+&
                    (xQ(ixQ1,ixQ2,2)-xFL(iFL,iLP,2))**2)
            if (dl<1.0d-2*dxQ1) then
              weight=(1/(1.0d-2*dxQ1))**weightIndex
            else
              weight=(1/dl)**weightIndex
            endif
            sumWeights(ixQ1,ixQ2)=sumWeights(ixQ1,ixQ2)+weight
            Qe(ixQ1,ixQ2)=Qe(ixQ1,ixQ2)+weight*QeL(iFL,iLP)
          enddo
        enddo
        !#
      enddo
    enddo

    ! for the right foot
    do iFL=1,numValidR
      do iLP=1,numRR(iFL)
        weightIndex=4
        Bratio=Busr/sqrt(wBR(iFL,iLP,mag(1))**2+wBR(iFL,iLP,mag(2))**2)
        if (Bratio<3) then
          dxMax1=ceiling(Bratio)*dxQ1
          dxMax2=ceiling(Bratio)*dxQ2
        else
          dxMax1=3*dxQ1
          dxMax2=3*dxQ2
        endif

        ixQmin1=floor((xFR(iFL,iLP,1)-dxMax1-xQmin1)/dxQ1)+1
        ixQmin2=floor((xFR(iFL,iLP,2)-dxMax2-xQmin2)/dxQ2)+1
        ixQmax1=floor((xFR(iFL,iLP,1)+dxMax1-xQmin1)/dxQ1)+1
        ixQmax2=floor((xFR(iFL,iLP,2)+dxMax2-xQmin2)/dxQ2)+1
        if (ixQmin1<1) ixQmin1=1
        if (ixQmin2<1) ixQmin2=1
        if (ixQmax1>numxQ1) ixQmax1=numxQ1
        if (ixQmax2>numxQ2) ixQmax2=numxQ2

        !#
        do ixQ1=ixQmin1,ixQmax1
          do ixQ2=ixQmin2,ixQmax2
            dl=sqrt((xQ(ixQ1,ixQ2,1)-xFR(iFL,iLP,1))**2+&
                    (xQ(ixQ1,ixQ2,2)-xFR(iFL,iLP,2))**2)
            if (dl<1.0d-2*dxQ1) then
              weight=(1/(1.0d-2*dxQ1))**weightIndex
            else
              weight=(1/dl)**weightIndex
            endif
            sumWeights(ixQ1,ixQ2)=sumWeights(ixQ1,ixQ2)+weight
            Qe(ixQ1,ixQ2)=Qe(ixQ1,ixQ2)+weight*QeR(iFL,iLP)
          enddo
        enddo
        !#
      enddo
    enddo


    ! divide by total weights
    do ixQ1=1,numxQ1
      do ixQ2=1,numxQ2
        if (sumWeights(ixQ1,ixQ2)>0) then
          Qe(ixQ1,ixQ2)=Qe(ixQ1,ixQ2)/sumWeights(ixQ1,ixQ2)
        endif
      enddo
    enddo





  end subroutine get_Qe

  subroutine getlQ(lQgrid,ixI^L,ixO^L,qt,w,x)
    !!calculate localized heating at chromosphere
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision :: lQgrid(ixI^S)

    integer :: ix^D,ixO^D,ixQ^L,numQ
    double precision :: xc^L,sumQ


    select case(iprob)    

    case(3)
      lQgrid(ixO^S)=0
      {do ixO^DB=ixOmin^DB,ixOmax^DB\}
        xcmin^D=x(ixO^DD,^D)-0.5d0*dxlevel(^D);
        xcmax^D=x(ixO^DD,^D)+0.5d0*dxlevel(^D);
        ixQmin^D=floor((xcmin^D-xQmin^D)/dxQ^D)+1;
        ixQmax^D=floor((xcmax^D-xQmin^D)/dxQ^D);

        sumQ=0.d0
        numQ=0
        do ix1=ixQmin1,ixQmax1
          do ix2=ixQmin2,ixQmax2
            if (ix1>=1 .and. ix1<=numXQ1 .and. ix2>=1 .and. ix2<=numXQ2) then
              sumQ=sumQ+Qe(ix1,ix2)
              numQ=numQ+1
            endif
          enddo
        enddo

        if (numQ>0) lQgrid(ixO^D)=sumQ/numQ
      {enddo\}

    end select

  end subroutine getlQ

  subroutine get_spectra(delta,Ec,dEe,Ee,spectra,numEe)
    use mod_global_parameters

    double precision :: delta,Ec,dEe
    integer :: numEe
    double precision :: Ee(numEe),spectra(numN,numEe)

    integer :: iEe,iNcol
    double precision :: Ncol(numN)

    double precision :: K,lambda_,beta
    double precision :: Fe0,Fe,qt,B0,const,const2,keV_erg
    double precision :: Eph,alpha,mec2,r0,N0,c,temp
    double precision :: E0min,E0max,spec0min,spec0max,dEe0
    double precision :: sigma0
    character (20) :: fname

    mec2=511.0      ! energy of a static electron [keV]
    alpha=1.0/137.0 ! fine structure constant
    c=2.9979d10     ! light speed [cm/s]
    r0=2.8179d-13   ! radiu of electron [cm]
    keV_erg=1.0d3*const_ev  ! 1 keV = * erg
    sigma0=7.9e-25  ! [cm^2 keV]
    K=2*dpi*const_e**4
    lambda_=25.0
    beta=2.0


    ! initial fast electron spectra
    spectra=0.d0
    const=(delta-1)/Ec
    do iEe=1,numEe
      Ee(iEe)=Ec+(iEe-1)*dEe
      spectra(1,iEe)=const*(Ee(iEe)/Ec)**(-delta)
    enddo

    ! change of electron's energy
    Ncol=0.d0
    do iNcol=2,numN
      Ncol(iNcol)=Nmin*10**(dNlog*(iNcol-1.d0))
      temp=2.d0*lambda_*K*Ncol(iNcol)/keV_erg**2
      do iEe=1,numEe
        E0min=sqrt(Ee(iEe)**2+temp)
        E0max=sqrt((Ee(iEe)+dEe)**2+temp)
        dEe0=E0max-E0min
        spec0min=const*(E0min/Ec)**(-delta)
        spec0max=const*(E0max/Ec)**(-delta)
        spectra(iNcol,iEe)=0.5d0*(spec0min+spec0max)*dEe0/dEe
      enddo
    enddo

    spectra=spectra*keV_erg   ! [electrons cm^-2 s^-1 keV^-1]


    !! output e flux
    !if (mype==0) then
    !  write(fname,'(a17)') 'spectra_table.txt'
    !  open(1,file=fname)
    !  write(1,*) 'numN numE'
    !  write(1,*) numN,numEe
    !  do iNcol=1,numN
    !    write(1,*) Ncol(iNcol)
    !    write(1,*) 'E spectra'
    !    do iEe=1,numEe
    !      write(1,'(e15.7, e15.7)') , Ee(iEe),spectra(iNcol,iEe)
    !    enddo
    !  enddo
    !  close(1)
    !endif



  end subroutine get_spectra

  subroutine get_HXR_line(Ee,spectra,numEe,dEe,xF,wBLR,numValid,numTurn,numR,eFlux,muB,HXRLR)
    use mod_global_parameters

    double precision :: dEe
    integer :: numEe,numValid
    double precision :: Ee(numEe),spectra(numN,numEe)
    integer :: numR(numFL),numTurn(numFL)
    double precision :: xF(numFL,numLP,ndim),QeLR(numFL,numLP)
    double precision :: wBLR(numFL,numLP,nw+ndir)
    double precision :: eFlux(numFL),muB(numFL),HXRLR(numFL,numLP)

    integer :: j,ix^D,iNlog,ireturn
    double precision :: Ncol(numLP),Np(numFL,numLP),Bv(numFL,numLP,ndim)
    double precision :: mu0,mu,Bxy0,Bxy,BxyRT,Bratio,dl1,dl2,dl,const,muMin
    logical :: Ereturn
    double precision :: Ephmin,Ephmax,dEph,Eph
    integer :: numEph,iEph,iEe
    double precision :: sigma0,sigmaBH,keV_erg,temp,Fe,depth
    double precision :: HXRLRs(numFL,numLP)
    double precision :: BxyTurn

    HXRLR=0.d0
    HXRLRs=0.d0

    keV_erg=1.0d3*const_ev  ! 1 keV = * erg
    sigma0=7.9e-25  ! [cm^2 keV]
    depth=1.0d9
    
    mu0=0.95
    muMin=0.1
    Ephmin=25
    Ephmax=50
    dEph=1.d0
    numEph=floor((Ephmax-Ephmin)/dEph)


    ! density and magnetic field
    Np(:,:)=wBLR(:,:,rho_)*unit_numberdensity
    do j=1,ndim
      Bv(:,:,j)=wBLR(:,:,mag(j))
    enddo


    ! calculate fast electron energy flux for each field line
    do ix1=1,numValid

      ! initial width of the flux tube element
      Bxy0=dsqrt(Bv(ix1,1,1)**2+Bv(ix1,1,2)**2)
      BxyTurn=dsqrt(Bv(ix1,numTurn(ix1),1)**2+Bv(ix1,numTurn(ix1),2)**2)
      mu0=0.00001
      if (muB(ix1)*BxyTurn>=1) mu0=sqrt(1-muB(ix1)*BxyTurn)


      if (mype==mod(ix1,npe)) then
        Ncol=0
        FIELD1: do ix2=numTurn(ix1),numR(ix1)
          Bxy=dsqrt(Bv(ix1,ix2,1)**2+Bv(ix1,ix2,2)**2)
          Bratio=Bxy/Bxy0
          const=muB(ix1)
          if (const*Bxy>=1) exit FIELD1
          mu=sqrt(1-const*Bxy)
          if (mu<muMin) mu=muMin

          ! column depth
          dl1=xF(ix1,ix2+1,1)-xF(ix1,ix2,1)
          dl2=xF(ix1,ix2+1,2)-xF(ix1,ix2,2)
          dl=dsqrt(dl1**2+dl2**2)*unit_length/mu
          Ncol(ix2)=Ncol(ix2-1)+dl*(Np(ix1,ix2)+Np(ix1,ix2-1))/2.0

          ! local heating rate
          iNlog=int(log10(Ncol(ix2)/Nmin)/dNlog+0.5)+1
          if (iNlog<1) iNlog=1
          if (iNlog>numN) iNlog=numN

         !  calulate HXR flux for given energy range
          do iEph=1,numEph
            Eph=Ephmin+(iEph-1)*dEph
            do iEe=1,numEe
              if (Ee(iEe)>Eph) then
                temp=sqrt(1-Eph/Ee(iEe))
                sigmaBH=(sigma0/(Eph*Ee(iEe)))*log((1+temp)/(1-temp))
                HXRLRs(ix^D)=HXRLRs(ix^D)+spectra(iNlog,iEe)*Np(ix^D)*sigmaBH*&
                                          (eFlux(ix1)*Bxy0/Bxy)*(mu0/mu)/mu
              endif
            enddo
          enddo
        enddo FIELD1

        Ncol=0
        FIELD2: do ix2=numTurn(ix1)-1,1,-1
          Bxy=dsqrt(Bv(ix1,ix2,1)**2+Bv(ix1,ix2,2)**2)
          Bratio=Bxy/Bxy0
          const=muB(ix1)
          if (const*Bxy>=1) exit FIELD2
          mu=sqrt(1-const*Bxy)
          if (mu<muMin) mu=muMin

          ! column depth
          dl1=xF(ix1,ix2+1,1)-xF(ix1,ix2,1)
          dl2=xF(ix1,ix2+1,2)-xF(ix1,ix2,2)
          dl=dsqrt(dl1**2+dl2**2)*unit_length/mu
          Ncol(ix2)=Ncol(ix2+1)+dl*(Np(ix1,ix2)+Np(ix1,ix2+1))/2.0

          ! local heating rate
          iNlog=int(log10(Ncol(ix2)/Nmin)/dNlog+0.5)+1
          if (iNlog<1) iNlog=1
          if (iNlog>numN) iNlog=numN

          ! calulate HXR flux for given energy range
          do iEph=1,numEph
            Eph=Ephmin+(iEph-1)*dEph
            do iEe=1,numEe
              if (Ee(iEe)>Eph) then
                temp=sqrt(1-Eph/Ee(iEe))
                sigmaBH=(sigma0/(Eph*Ee(iEe)))*log((1+temp)/(1-temp))
                HXRLRs(ix^D)=HXRLRs(ix^D)+spectra(iNlog,iEe)*Np(ix^D)*sigmaBH*&
                                          (eFlux(ix1)*Bxy0/Bxy)*(mu0/mu)/mu
              endif
            enddo
          enddo
        enddo FIELD2

      endif
    enddo


    call MPI_ALLREDUCE(HXRLRs,HXRLR,numFL*numLP,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,icomm,ierrmpi)

    HXRLR=HXRLR*depth

  end subroutine get_HXR_line

  subroutine interp_HXR(xFL,xFR,wBL,wBR,HXRL,HXRR)
    ! get local heating rate via interpolation
    use mod_global_parameters
    use mod_usr_methods

    double precision :: xFL(numFL,numLP,ndim),xFR(numFL,numLP,ndim)
    double precision :: wBL(numFL,numLP,nw+ndir),wBR(numFL,numLP,nw+ndir)
    double precision :: HXRL(numFL,numLP),HXRR(numFL,numLP)
 
    double precision :: sumWeights(numxQ1,numxQ2)
    double precision :: dxMax1,dxMax2,weight,dl,Bratio
    integer :: iFL,iLP,ixQ^L,ixQ^D,weightIndex

    HXR=0
    sumWeights=0
    dxMax1=8*dxQ1
    dxMax2=8*dxQ2
    weightIndex=1


    ! for the left foot
    do iFL=1,numValidL
      do iLP=1,numRL(iFL)
        weightIndex=4
        Bratio=Busr/sqrt(wBL(iFL,iLP,mag(1))**2+wBL(iFL,iLP,mag(2))**2)
        if (Bratio<3) then
          dxMax1=ceiling(Bratio)*dxQ1
          dxMax2=ceiling(Bratio)*dxQ2
        else
          dxMax1=3*dxQ1
          dxMax2=3*dxQ2
        endif

        ixQmin1=floor((xFL(iFL,iLP,1)-dxMax1-xQmin1)/dxQ1)+1
        ixQmin2=floor((xFL(iFL,iLP,2)-dxMax2-xQmin2)/dxQ2)+1
        ixQmax1=floor((xFL(iFL,iLP,1)+dxMax1-xQmin1)/dxQ1)+1
        ixQmax2=floor((xFL(iFL,iLP,2)+dxMax2-xQmin2)/dxQ2)+1
        if (ixQmin1<1) ixQmin1=1
        if (ixQmin2<1) ixQmin2=1
        if (ixQmax1>numxQ1) ixQmax1=numxQ1
        if (ixQmax2>numxQ2) ixQmax2=numxQ2

        !#
        do ixQ1=ixQmin1,ixQmax1
          do ixQ2=ixQmin2,ixQmax2
            dl=sqrt((xQ(ixQ1,ixQ2,1)-xFL(iFL,iLP,1))**2+&
                    (xQ(ixQ1,ixQ2,2)-xFL(iFL,iLP,2))**2)
            if (dl<1.0d-2*dxQ1) then
              weight=(1/(1.0d-2*dxQ1))**weightIndex
            else
              weight=(1/dl)**weightIndex
            endif
            sumWeights(ixQ1,ixQ2)=sumWeights(ixQ1,ixQ2)+weight
            HXR(ixQ1,ixQ2)=HXR(ixQ1,ixQ2)+weight*HXRL(iFL,iLP)
          enddo
        enddo
        !#
      enddo
    enddo

    ! for the right foot
    do iFL=1,numValidR
      do iLP=1,numRR(iFL)
        weightIndex=4
        Bratio=Busr/sqrt(wBR(iFL,iLP,mag(1))**2+wBR(iFL,iLP,mag(2))**2)
        if (Bratio<3) then
          dxMax1=ceiling(Bratio)*dxQ1
          dxMax2=ceiling(Bratio)*dxQ2
        else
          dxMax1=3*dxQ1
          dxMax2=3*dxQ2
        endif

        ixQmin1=floor((xFR(iFL,iLP,1)-dxMax1-xQmin1)/dxQ1)+1
        ixQmin2=floor((xFR(iFL,iLP,2)-dxMax2-xQmin2)/dxQ2)+1
        ixQmax1=floor((xFR(iFL,iLP,1)+dxMax1-xQmin1)/dxQ1)+1
        ixQmax2=floor((xFR(iFL,iLP,2)+dxMax2-xQmin2)/dxQ2)+1
        if (ixQmin1<1) ixQmin1=1
        if (ixQmin2<1) ixQmin2=1
        if (ixQmax1>numxQ1) ixQmax1=numxQ1
        if (ixQmax2>numxQ2) ixQmax2=numxQ2

        !#
        do ixQ1=ixQmin1,ixQmax1
          do ixQ2=ixQmin2,ixQmax2
            dl=sqrt((xQ(ixQ1,ixQ2,1)-xFR(iFL,iLP,1))**2+&
                    (xQ(ixQ1,ixQ2,2)-xFR(iFL,iLP,2))**2)
            if (dl<1.0d-2*dxQ1) then
              weight=(1/(1.0d-2*dxQ1))**weightIndex
            else
              weight=(1/dl)**weightIndex
            endif
            sumWeights(ixQ1,ixQ2)=sumWeights(ixQ1,ixQ2)+weight
            HXR(ixQ1,ixQ2)=HXR(ixQ1,ixQ2)+weight*HXRR(iFL,iLP)
          enddo
        enddo
        !#
      enddo
    enddo


    ! divide by total weights
    do ixQ1=1,numxQ1
      do ixQ2=1,numxQ2
        if (sumWeights(ixQ1,ixQ2)>0) then
          HXR(ixQ1,ixQ2)=HXR(ixQ1,ixQ2)/sumWeights(ixQ1,ixQ2)
        endif
      enddo
    enddo

  end subroutine interp_HXR

  subroutine get_HXR(HXRgrid,ixI^L,ixO^L,qt,w,x)
    !!calculate localized heating at chromosphere
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision :: HXRgrid(ixI^S)

    integer :: ix^D,ixO^D,ixQ^L,numHXR
    double precision :: xc^L,sumHXR


    select case(iprob)    

    case(3)
      HXRgrid(ixO^S)=0
      {do ixO^DB=ixOmin^DB,ixOmax^DB\}
        xcmin^D=x(ixO^DD,^D)-0.5d0*dxlevel(^D);
        xcmax^D=x(ixO^DD,^D)+0.5d0*dxlevel(^D);
        ixQmin^D=floor((xcmin^D-xQmin^D)/dxQ^D)+1;
        ixQmax^D=floor((xcmax^D-xQmin^D)/dxQ^D);

        sumHXR=0.d0
        numHXR=0
        do ix1=ixQmin1,ixQmax1
          do ix2=ixQmin2,ixQmax2
            if (ix1>=1 .and. ix1<=numXQ1 .and. ix2>=1 .and. ix2<=numXQ2) then
              sumHXR=sumHXR+HXR(ix1,ix2)
              numHXR=numHXR+1
            endif
          enddo
        enddo

        if (numHXR>0) HXRgrid(ixO^D)=sumHXR/numHXR
      {enddo\}

    end select

  end subroutine get_HXR

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
    double precision :: lQgrid(ixI^S),ens(ixI^S),HXRgrid(ixI^S)
    double precision :: t
    double precision :: Efield(ixI^S,1:ndir),v(ixI^S,1:ndir)

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
    call getlQ(lQgrid,ixI^L,ixO^L,t,w,x)
    w(ixO^S,nw+9)=lQgrid(ixO^S)

    call getcQ(ens,ixI^L,ixO^L,t,w,x)
    w(ixO^S,nw+10)=ens(ixO^S)

    !! background heating
    call getbQ(ens,ixI^L,ixO^L,t,w,x)
    w(ixO^S,nw+11)=ens(ixO^S)

    ! store the cooling rate 
    if(mhd_radiative_cooling)call getvar_cooling(ixI^L,ixO^L,w,x,ens)
    w(ixO^S,nw+12)=ens(ixO^S)

    call get_HXR(HXRgrid,ixI^L,ixO^L,t,w,x)
    w(ixO^S,nw+13)=HXRgrid(ixO^S)

    ! electric field
    do idir=1,ndir
      v(ixO^S,idir)=w(ixO^S,mom(idir))/w(ixO^S,rho_)
    enddo
    call cross_product(ixI^L,ixO^L,Btotal,v,Efield)
    do idir=1,ndir
      Efield(ixO^S,idir)=Efield(ixO^S,idir)+divb(ixO^S)*current_o(ixO^S,idir)
      w(ixO^S,nw+13+idir)=Efield(ixO^S,idir)
    enddo

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te Alfv divB beta j1 j2 j3 eta lQ cQ bQ rad HXR_25_50 E1 E2 E3'
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
    double precision :: rad(ixI^S),heta,reta,reta2,heta2
    double precision :: jc,jabs(ixI^S),vd(ixI^S),rad2(ixI^S)

    if (iprob==1) then
      eta(ixO^S)=0.d0

    else if (iprob>=2) then
      heta = 5.
      heta2 = 1.
      reta = 0.8d0 * 0.3d0
      !eta1 = 0.01d0
      !tar= 0.0d0
      if (global_time<tar) then
        rad(ixO^S)=dsqrt(x(ixO^S,1)**2+(x(ixO^S,2)-heta)**2)
        where (rad(ixO^S) .lt. reta)
          eta(ixO^S)=eta1*(2.d0*(rad(ixO^S)/reta)**3-3.d0*(rad(ixO^S)/reta)**2+1.d0)
        elsewhere
          eta(ixO^S)=zero
        endwhere
      else
        vd(ixO^S)=dsqrt(sum(current(ixO^S,:)**2,dim=ndim+1))/w(ixO^S,rho_)/q_e
        rad2(ixO^S)=exp(-((x(ixO^S,2)-heta)**2)/heta2**2)
        where(vd(ixO^S)>vc)
          eta(ixO^S)=eta2*(vd(ixO^S)/vc-1.d0)*rad2(ixO^S)
        elsewhere
          eta(ixO^S)=0.d0
        endwhere
        where(eta(ixO^S)>etam)
          eta(ixO^S)=etam
        endwhere
        !eta(ixO^S)=eta0
      end if 
    endif


    !if (iprob==1) then
    !  eta(ixO^S)=0.d0

    !else if (iprob>=2) then
    !  heta = 5.
    !  reta = 0.8d0 * 0.3d0
    !  reta2 = 2.d0
    !  !eta1 = 0.01d0
    !  !tar= 0.0d0
    !  if (global_time<tar) then
    !    rad(ixO^S)=dsqrt(x(ixO^S,1)**2+(x(ixO^S,2)-heta)**2)
    !    where (rad(ixO^S) .lt. reta)
    !      eta(ixO^S)=eta1*(2.d0*(rad(ixO^S)/reta)**3-3.d0*(rad(ixO^S)/reta)**2+1.d0)
    !    elsewhere
    !      eta(ixO^S)=zero
    !    endwhere
    !  else
    !    rad(ixO^S)=dsqrt(x(ixO^S,1)**2+(x(ixO^S,2)-heta)**2)
    !    vd(ixO^S)=dsqrt(sum(current(ixO^S,:)**2,dim=ndim+1))/w(ixO^S,rho_)/q_e
    !    where(vd(ixO^S)>vc)
    !      eta(ixO^S)=eta2*(vd(ixO^S)/vc-1.d0)
    !    elsewhere
    !      eta(ixO^S)=0.d0
    !    endwhere
    !    where(eta(ixO^S)>etam)
    !      eta(ixO^S)=etam
    !    endwhere
    !    !eta(ixO^S)=eta0
    !  end if 
    !endif


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
