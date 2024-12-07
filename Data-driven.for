!module MyModule
!    implicit none
!    integer :: i1, i2, i3, num
!    real :: data1(112860, 2) ! 存储phi_ns数据的数组
!    real :: data2(112860, 3) ! 存储dphi_n数据的数组
!    real :: data3(338580, 3) ! 存储d2phi_n数据的数组
!
!contains
!
!    subroutine ReadData()
!        open(105, file="C:\Users\zxm\Desktop\jianhau\data\phi_ns.csv")
!        do i1 = 1, 112860
!            read(105, *) data1(i1, :)
!        end do
!        close(105)
!
!        open(105, file="C:\Users\zxm\Desktop\jianhau\data\dphi_n.csv")
!        do i2 = 1, 112860
!            read(105, *) data2(i2, :)
!        end do
!        close(105)
!
!        open(105, file="C:\Users\zxm\Desktop\jianhau\data\d2phi_n.csv")
!        do i3 = 1, 338580
!            read(105, *) data3(i3, :)
!        end do
!        close(105)
!
!        num = 1
!    end subroutine ReadData
!
!    ! subroutine InitData()
!    ! call ReadData()  ! 读取数据并存储到 data1 数组中
!    ! end subroutine InitData
!
!end module MyModule

module random_forest
    implicit none
    type Tree
        integer :: split_feature
        real :: split_value
        real :: leaf_value
        type(Tree), pointer :: tree_left => null()
        type(Tree), pointer :: tree_right => null()
    contains
        procedure :: calc_predict_value
    end type Tree

    type RandomForestRegression
        integer :: n_estimators
        integer :: max_depth
        integer :: min_samples_split
        integer :: min_samples_leaf
        real :: min_split_gain
        integer :: colsample_bytree
        real :: subsample
        integer :: random_state
        character(len=10) :: colsample_bytree_str
        !type(Tree), allocatable, target :: trees(:)
        type(Tree), pointer :: trees(:) 
        !integer, dimension(:), allocatable :: feature_importances
        
    contains
        procedure :: fit => random_forest_fit
        procedure :: predict => random_forest_predict
    end type RandomForestRegression

contains

    subroutine random_forest_fit(clf, file)
        implicit none
        class(RandomForestRegression), intent(inout) :: clf
        !real, dimension(:,:), intent(in) :: dataset
        !real, dimension(:), intent(in) :: targets
        character(len=100) :: file
        !file = 'C:\Users\zxm\Desktop\testf90\testf90\modelphi.dat'
        call load_model(clf, file)
    end subroutine random_forest_fit
    
    recursive subroutine load_tree(unit, trees)
        implicit none
        integer, intent(in) :: unit
        type(Tree), pointer :: trees
        integer :: node_exists

        read(unit, *) node_exists
        if (node_exists == 1) then
            allocate(trees)
            read(unit, *) trees%split_feature
            read(unit, *) trees%split_value
            read(unit, *) trees%leaf_value
            call load_tree(unit, trees%tree_left)
            call load_tree(unit, trees%tree_right)
        else
            nullify(trees)
        end if
    end subroutine load_tree

    subroutine load_model(clf, filename)
        implicit none
        class(RandomForestRegression), intent(inout) :: clf
        character(len=*), intent(in) :: filename
        integer :: i, unit
        type(Tree), pointer :: tree_ptr

        ! Open the file for reading
        open(newunit=unit, file=filename, status='old', action='read')

        ! Read the number of trees
        read(unit, *) clf%n_estimators
        read(unit, *) clf%max_depth
        read(unit, *) clf%min_samples_split
        read(unit, *) clf%min_samples_leaf
        read(unit, *) clf%min_split_gain
        read(unit, *) clf%colsample_bytree
        read(unit, *) clf%subsample
        read(unit, *) clf%random_state
        read(unit, '(A)') clf%colsample_bytree_str

        ! Allocate memory for trees
        allocate(clf%trees(clf%n_estimators))

        ! Read each tree's data
        do i = 1, clf%n_estimators
            !tree_ptr => clf%trees(i)
            call load_tree(unit, tree_ptr)
            clf%trees(i) = tree_ptr
        end do

        ! Close the file
        close(unit)
    end subroutine load_model

    recursive real function calc_predict_value(this, row)
        class(Tree), intent(in) :: this
        real, dimension(:), intent(in) :: row
        if (associated(this%tree_left) .and. associated(this%tree_right)) then
            if (row(this%split_feature) <= this%split_value) then
                calc_predict_value = this%tree_left%calc_predict_value(row)
            else
                calc_predict_value = this%tree_right%calc_predict_value(row)
            end if
        else
            calc_predict_value = this%leaf_value
        end if
    end function calc_predict_value

    subroutine random_forest_predict(self, dataset, res)
        class(RandomForestRegression), intent(in) :: self
        real, dimension(:,:), intent(in) :: dataset
        real, dimension(:), allocatable, intent(out) :: res

        integer :: i, j, num_trees, num_samples
        real :: pred_sum
        real, dimension(:), allocatable :: pred_list

        num_trees = size(self%trees)
        num_samples = size(dataset, 1)

        allocate(res(num_samples))
        allocate(pred_list(num_trees))

        do i = 1, num_samples
            pred_sum = 0.0
            do j = 1, num_trees
                pred_list(j) = self%trees(j)%calc_predict_value(dataset(i, :))
            end do
            res(i) = sum(pred_list) / num_trees
        end do

        deallocate(pred_list)
    end subroutine random_forest_predict

end module random_forest

Subroutine umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, props, nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, jstep, kinc)
!
  !INCLUDE 'ABA_PARAM.INC'
  !use MyModule
  use random_forest
  Implicit None
!
!
  Character *80 cmname
  Dimension stress(ntens), statev(nstatv), ddsdde(ntens, ntens), ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens), time(2), predef(1), dpred(1), props(nprops), coords(3), drot(3, 3), dfgrd0(3, 3), dfgrd1(3, 3), jstep(4)

  Integer ndi, nshr, ntens, nstatv, nprops, noel, npt, layer, kspt, jstep, kinc, kstep
  Real *8 stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, predef, dpred, props, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, syield

  Character *255 fndia, fnstr
  Integer nhrdc, nhrdp
  Parameter (nhrdc=4, nhrdp=1)

  Dimension cel(ntens, ntens), phi_ns(0:1), dphi_n(ntens), d2phi_n(ntens, ntens), stran_el_ns(0:1, ntens), stran_pl_ns(0:1, ntens), dstran_el(ntens), dstran_pl(ntens), spr(ntens), epr(ntens), hrdc(nhrdc), hrdp(nhrdp), stress_ns(0:1, ntens), dstress(ntens)
  Real *8 e, enu, g, ek, cel, dstress, deeq
!     predictor stress
  Real *8 spr, epr, stress_ns
  Real *8 hrdc, hrdp
  Dimension eeq_ns(0:1)
  Real *8 eeq_ns, flow_stress, dflow_stress
!     eeq_ns: eq plastic strain at steps n,n+1

  Real *8 phi_ns, dphi_n, d2phi_n
  Real *8 stran_el_ns, stran_pl_ns, dstran_el, dstran_pl
  Integer imsg, idia, i, istv, ihrd_law, iyld_law
  Logical idiaw, failnr
  Real *8 empa, gpa
  Parameter (empa=1D6, gpa=1D9)

!$$$  yld function parameters
  Integer nyldp, nyldc
  Parameter (nyldp=1, nyldc=9) ! this depends on iyld_law...          
  Dimension yldp_ns(0:1, nyldp), yldc(nyldc)
  Real *8 yldp_ns, yldc, deq_pr, deq

!$$$  Material data
  Dimension rvs(4)
  Real *8 rvs
! 以下代码为新增代码
  integer,save :: enterNumber = 0
  logical,save :: firstrun = .true.
  
  type(RandomForestRegression),save :: clf1, clf2, clf3, clf4
  character(len=100) :: file
  
  ! integer :: tempvar
  print *, "开始执行"
  if (firstrun) then
      !write(*,*)"please input an integer"
      file = 'G:\DataDriven\fourth\data\phi_1.dat'
      call clf1%fit(file)
      
      file = 'G:\DataDriven\fourth\data\dphi1_1.dat'
      call clf2%fit(file)
      
      file = 'G:\DataDriven\fourth\data\dphi2_1.dat'
      call clf3%fit(file)
      
      file = 'G:\DataDriven\fourth\data\dphi3_1.dat'
      call clf4%fit(file)
      
      !call ReadData()
      
      !read(*,*)tempvar
      firstrun = .false.
  end if
  
  enterNumber = enterNumber + 1

  !write(*,*) 'enterNumber', enterNumber, 'DSTRAN', DSTRAN(3)

  !***  hardwired material parameters
  !-----------------------------------------------------------------------
  !**   material constitutive laws
  ihrd_law = 0                ! Voce isotropic hardening
  !**   hardening parameters
  hrdc(1) = 747.852d0
  hrdc(2) = 0.01096d0
  hrdc(3) = 0.13768d0
  hrdc(4) = 70.87d0
  !-----------------------------------------------------------------------
  iyld_law = 2                ! yld2000-2d
  
  yldc(1) = 1.14315155941131d0
  yldc(2) = 0.888096416469737d0
  yldc(3) = 1.64578651879137d0
  yldc(4) = 0.994396992277993d0
  yldc(5) = 0.709496837009325d0
  yldc(6) = 0.140553274780613d0
  yldc(7) = 1.02562336441705d0
  yldc(8) = 0.642006782583844d0
  
  yldc(9) = 8.d0
  !-----------------------------------------------------------------------
  !**   elastic constants
  e = 71.967d0 * gpa
  enu = 0.30487d0
  !-----------------------------------------------------------------------
  imsg = 7
  idia = 0                    ! 0 (in case stdo is preferred)
  istv = 425 ! reserved for state variable output file
  stress_ns(0,:) = stress(:)
  
  idiaw = .false.
  
  ! restore stran_el, stran_pl, and yldp
  call restore_statev(statev, nstatv, eeq_ns(0), stran_el_ns(0,:), &
       stran_pl_ns(0,:), ntens, yldp_ns(0,:), nyldp, 0, .false., idia, &
       .false., kinc, noel, npt, time(1), stress)
  
  ! yld parameters pertaining to step n - need to find
  ! the updated yld parameters for step n+1 later...
  
  ! Moduluar pseudo code for stress integration
  !-----------------------------------------------------------------------
  ! i.   Check the current (given) variables
  call emod_iso_shell(e, enu, G, eK, Cel)
  !-----------------------------------------------------------------------
  ! ii.  Trial stress calculation
  ! Calculate predict elastic strain - assuming dstran = dstran_el
  call add_array2(stran_el_ns(0,:), dstran, epr, ntens)
  call mult_array(Cel, epr, ntens, stress_ns(1,:))
  spr(:) = stress_ns(1,:)
  
  !-----------------------------------------------------------------------
  ! iii. See if the pr stress (stress_ns(1,:)) calculation is in the plastic or
  !      elastic regime
  deq_pr = 0.0               ! assuming no plastic increment
  hrdp(1) = eeq_ns(0) + deq_pr
  call uhard(ihrd_law, hrdp, nhrdp, hrdc, nhrdc, &
       flow_stress, dflow_stress, empa)
  
  call update_yldp(iyld_law, yldp_ns, nyldp, deq_pr)
  
  ! predicting yield surface
   !call yld(iyld_law, yldp_ns(1,:), yldc, nyldp, nyldc, stress_ns(1,:), &
   !     phi_ns(0), dphi_n, d2phi_n, ntens)


  if (all(abs(stress_ns(1,:)) > 10 .and. abs(stress_ns(1,:)) < 1500000000)) then
      ! 缩放数据
      call process_stress(stress_ns(1,:), phi_ns(0), dphi_n, clf1, clf2, clf3, clf4)
  else
      call yld(iyld_law, yldp_ns(1,:), yldc, nyldp, nyldc, stress_ns(1,:), phi_ns(0), dphi_n, d2phi_n, ntens)
      print *, "调用了原函数"
  endif
          
  print *, "调用次数:", enterNumber
          
  !phi_ns = data1(enterNumber,:)
  !dphi_n = data2(enterNumber,:)
  !d2phi_n = data3(enterNumber*3-2:enterNumber*3,:) !这一行从未用过，可以直接赋零
          
  !open(105, file='C:\Users\zxm\Desktop\cae_test\rf_yl90_dphi_n.txt', position='append')
  !write(105,*) dphi_n
  !
  !open(106, file='C:\Users\zxm\Desktop\cae_test\rf_yl90_phi_ns.txt', position='append')
  !write(106,*) phi_ns
  !
  !open(107, file='C:\Users\zxm\Desktop\cae_test\rf_yl90_stress_ns1.txt', position='append')
  !write(107,*) stress_ns(1,:)
          
  ! 以下代码不要动
  !-----------------------------------------------------------------------
  if (phi_ns(0) < flow_stress) then ! elastic
     call update_elastic(idia, idiaw, iyld_law, ntens, nyldp, nstatv, &
          ddsdde, cel, stran, stran_el_ns, stran_pl_ns, dstran, stress, &
          eeq_ns, deq, yldp_ns, statev, kinc)
     dstran_pl(:) = 0d0
     dstran_el(:) = dstran(:)
     dstress(:) = stress(:) - stress_ns(0,:)
  else ! plastic
!-----------------------------------------------------------------------
! vi. Return mapping
!    Return mapping subroutine updates stress/statev
     call return_mapping(Cel, stress_ns(1,:), phi_ns(0), eeq_ns(0), &
          dphi_n, dstran, stran_el_ns(0,:), stran_pl_ns(0,:), &
          ntens, idiaw, idia, hrdp, nhrdp, hrdc, nhrdc, ihrd_law, &
          iyld_law, yldc, nyldc, yldp_ns, nyldp, &
          stress, deeq, dstran_pl, dstran_el, statev, nstatv, ddsdde, &
          failnr, kinc, noel, npt, time, clf1, clf2, clf3, clf4)
! new stress and stress increment
     stress_ns(1,:) = stress(:)
     dstress(:) = stress(:) - stress_ns(0,:)
     
     if (failnr) then
        ! reduce time step?
        pnewdt = 0.5d0
        
        return ! out of umat
     endif
     !-----------------------------------------------------------------------
  endif
  
  spd = 0d0
  sse = 0d0
  do i = 1, ntens
     spd = spd + (stress_ns(0,i) + 0.5 * dstress(i)) * dstran_pl(i)
     sse = sse + (stress_ns(0,i) + 0.5 * dstress(i)) * dstran_el(i)
  enddo
  
  !-----------------------------------------------------------------------
  
  return
end subroutine umat
    
subroutine process_stress(stress_ns, phi_ns, dphi_n, clf1, clf2, clf3, clf4)
    use random_forest
    type(RandomForestRegression) :: clf1, clf2, clf3, clf4
    real, dimension(:), allocatable :: min_val, max_val
    real, dimension(:,:), allocatable :: X_test, X_test_scaled
    real, allocatable :: y1(:), y2(:), y3(:), y4(:)
    integer :: jjj, n_features
    real, intent(out) :: phi_ns
    real, intent(out), dimension(3) :: dphi_n
    real, intent(in), dimension(3) :: stress_ns
    ! 输入参数
    !real, intent(in) :: stress_ns(1, 3)
    !real, intent(in) :: min_val(3)
    !real, intent(in) :: max_val(3)
    
    
    ! 分配数组的内存
    allocate(max_val(3))
    allocate(min_val(3))
    allocate(X_test(1,3))
    allocate(X_test_scaled(1,3))
    ! 逐个元素赋值
    max_val(1) = 9.999879E+08
    max_val(2) = 9.9997402E+08
    max_val(3) = 9.9999802E+08
  
    min_val(1) = -9.9999603E+09
    min_val(2) = 2001.004
    min_val(3) = 3001.06
    ! 处理stress_ns
    X_test(1,1) = real(stress_ns(1))
    X_test(1,2) = real(stress_ns(2))
    X_test(1,3) = real(stress_ns(3))

    ! 缩放数据
    n_features = size(X_test, 2)
    do jjj = 1, n_features
        X_test_scaled(1, jjj) = (X_test(1, jjj) - min_val(jjj)) / (max_val(jjj) - min_val(jjj))
    end do
    ! 调用分类器预测
    call clf1%predict(X_test_scaled, y1)
    call clf2%predict(X_test_scaled, y2)
    call clf3%predict(X_test_scaled, y3)
    call clf4%predict(X_test_scaled, y4)

    phi_ns = y1(1)
    dphi_n(1) = y2(1)
    dphi_n(2) = y3(1)
    dphi_n(3) = y4(1)
    
    deallocate(max_val)
    deallocate(min_val)
    deallocate(X_test)
    deallocate(X_test_scaled)

end subroutine process_stress
! 以上是新增子程序
    
!-----------------------------------------------------------------------
!     subroutine that stores/restores state variables.
!     This subroutine makes transactions between statev and other
!     variables passed as arguments.
Subroutine restore_statev(statev, nstatv, eqpl, stran_el, stran_pl, ntens, yldp, nyldp, iopt, verbose, iunit, & 
  iw, kinc, noel, npt, time, stress)
!-----------------------------------------------------------------------
!     Arguments
!     statev  : state variable array
!     nstatv  : len of statev
!     eqpl    : equivalent plastic strain
!     stran_el: cumulative elastic strain
!     stran_pl: cumulative plastic strain
!     ntens   : Len of stran_el and stran_pl
!     yldp    : parameters for yield surface
!     nyldp   : Len of yldp ! mostly static as for isotropic hardening...
!     iopt    : option to define the behavior of restore_statev
!                  0: read from statev
!                  1: save to statev
!     verbose : flag to be verbose or not
!     iunit   : file unit (if not zero) or std to which
!               inquiry stream will be written
!     kinc    : increment number
!     noel    : element number
!     npt     : integration point number
!     time    : time array passed to umat
!     stress  : stress is not stored to statev
!               (this is just to write to fn_stav)
!-----------------------------------------------------------------------
  Implicit None
  Integer nstatv, ntens
  Dimension statev(nstatv), stran_el(ntens), stran_pl(ntens), yldp(nyldp), stress(ntens)
  Real *8 statev, eqpl, stran_el, stran_pl, yldp, time, stress
  Character *80 fmt, fnstv, ncc
  Integer iopt, i, nyldp, iunit, kinc, noel, npt
  Logical verbose, iw
  If (iopt==0) Then
         ! read from statev
    Do i = 1, ntens
      stran_el(i) = statev(i)
      stran_pl(i) = statev(i+ntens)
    End Do
    eqpl = statev(2*ntens+1)
    Do i = 1, nyldp
      yldp(i) = statev(2*ntens+i+1)
    End Do
  Else If (iopt==1) Then
         ! save to statev
         ! read from statev
    Do i = 1, ntens
      statev(i) = stran_el(i)
      statev(i+ntens) = stran_pl(i)
    End Do
    statev(2*ntens+1) = eqpl
    Do i = 1, nyldp
      statev(2*ntens+i+1) = yldp(i)
    End Do

  End If

  Return
End Subroutine restore_statev
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!$$$  3D shell with S11,S22 and S12
Subroutine emod_iso_shell(e, enu, g, ekappa, cel)
!     intent(in) e, nu
!     intent(out) G,ekappa,c
!      e    : Young's modulus
!      enu  : Poisson ratio
!      G    : Shear modulus
!     ekappa: bulk modulus
!      cel  : elastic constants
  Implicit None
  Real *8 cel(3, 3)
  Real *8 enu, e, x, g, ekappa
  Integer i, j
  cel(:, :) = 0D0
!     Multiplier
  x = e/(1D0+enu)/(1D0-2D0*enu)
  Do i = 1, 2
    Do j = 1, 2
      cel(i, j) = x*enu
    End Do
    cel(i, i) = x*(1D0-enu)
  End Do
  cel(3, 3) = x*(1D0-2D0*enu)/2D0
  ekappa = e/3D0/(1D0-2D0*enu)
  g = e/2D0/(1D0+enu)
!      call w_mdim(0,cel,3,1d0)
  Return
End Subroutine emod_iso_shell
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Apply incremental update on array
!     ci = ai + d_ai
Subroutine add_array2(ai, d_ai, ci, ntens)
!     intent(in) ai, d_ai, ntens
!     intent(out) ci
!     ai   : old array
!     d_ai : increments
!     ci   : new array
!     ntens: len
  Implicit None
  Integer ntens, i
  Real *8 ai(ntens), d_ai(ntens), ci(ntens)
  Do i = 1, ntens
    ci(i) = ai(i) + d_ai(i)
  End Do
  Return
End Subroutine add_array2
!-----------------------------------------------------------------------
!     Apply tensor inner dot
!     ci = aij x bj
Subroutine mult_array(aij, bj, ntens, ci)
  Implicit None
  Integer, Intent (In) :: ntens
  Dimension aij(ntens, ntens), bj(ntens), ci(ntens)
  Real *8, Intent (In) :: aij, bj
  Real *8, Intent (Out) :: ci
  Integer i, j
!f2py intent(in) aij, bj, ntens
!f2py intent(out) ci
!f2py depend(ntens) aij, bj
  ci(:) = 0.D0
  Do i = 1, ntens
    Do j = 1, ntens
      ci(i) = ci(i) + aij(i, j)*bj(j)
    End Do
  End Do
  Return
End Subroutine mult_array
!-----------------------------------------------------------------------
!     custom hardening law manager...
Subroutine uhard(ihard_law, hrdp, nhrdp, hrdc, nhrdc, flow_stress, dflow_stress, fact)
!     intent(in) ihard_law,hrdp,nhrdp,hrdc,nhrdc, fact
!     intent(out) flow_stress,dflow_stress

!     ihard_law: hardening law
!     hrdp     : state parameter for hardening
!     hrdc     : hardening constants (invariant)
!      flow_stress: flow stress
!     dflow_stress: slope of flow stress (df/de^eq)
!     fact        : multiplicative factor applied to flow/dflow
  Implicit None
  Integer ihard_law, nhrdp, nhrdc
  Dimension hrdp(nhrdp), hrdc(nhrdc)
  Real *8 flow_stress, dflow_stress, hrdp, hrdc, fact
!-----------------------------------------------------------------------
  If (ihard_law==0) & !        hrdp(1) : equivalent plastic strain
    Then
! voce                                  
    Call voce(hrdp(1), hrdc(1), hrdc(2), hrdc(3), hrdc(4), flow_stress, dflow_stress)
    flow_stress = flow_stress*fact
    dflow_stress = dflow_stress*fact
  Else

         !call exit(-1)
  End If

  Return
End Subroutine uhard
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Voce - rate independent
Subroutine voce(e, a, b0, c, b1, sig, dsig)
  Implicit None
  Real *8 e, a, b0, c, b1, sig, dsig
!     Voce
  sig = a*(b0+e)**c
  dsig = a*c*(b0+e)**(c-1)
  Return
End Subroutine voce
!-----------------------------------------------------------------------
Subroutine update_yldp(iyld_law, yldp_ns, nyldp, deeq)
!     Arguments
!     iyld_law : yield function choice
!     yldp_ns  : yield parameters stored for two separate steps
!     nyldp    : len of yield parameters for each separate step
!     deeq     : incremental equivalent plastic strain
!     intent(in) iyld_law,yldp_ns,nyldp,deeq
!-----------------------------------------------------------------------
  Implicit None
  Integer iyld_law, nyldp
  Dimension yldp_ns(0:1, nyldp)
  Real *8 yldp_ns, deeq

  yldp_ns(1, 1) = deeq + yldp_ns(0, 1)

  Return
End Subroutine update_yldp
!-----------------------------------------------------------------------
Subroutine yld(iyld_law, yldp, yldc, nyldp, nyldc, stress, phi, dphi, d2phi, ntens)
!-----------------------------------------------------------------------
!***  Arguments
!     iyld_law  : choice of yield function
!     yldp      : state variables associated with yield function
!     yldc      : constants for yield function
!     nyldp     : Len of yldp
!     nyldc     : Len of yldc
!     stress    : stress tensor (cauchy stress is generally expected)
!     phi       : yield surface
!     dphi      : 1st derivative of yield surface w.r.t. stress
!     d2phi     : 2nd derivative of yield surface w.r.t. stress
!     ntens     : Len of stress
!-----------------------------------------------------------------------
  Implicit None
  Integer iyld_law, nyldp, nyldc, ntens
  Dimension yldp(nyldp), yldc(nyldc), dphi(ntens), d2phi(ntens, ntens)
  Real *8 yldp, yldc, phi, dphi, d2phi
!***  Local variables for better readibility
  Dimension stress(ntens)
  Real *8 stress
!***  Local - control
  Integer imsg
  Logical idiaw
!-----------------------------------------------------------------------
!f2py intent(in) iyld_law,yldp,yldc,nyldp,nyldc,stress,ntens
!f2py intent(out) phi,dphi,d2phi
!f2py depend(nyldp) yldp
!f2py depend(nyldc) yldc
!f2py depend(ntens) stress,dphi,d2phi
!***  Define phi,dphi,d2phi

  idiaw = .False.
  imsg = 0
  Call yld2000_2d(stress, phi, dphi, d2phi, yldc)

End Subroutine yld
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Subroutine yld2000_2d(cauchy, phi, dphi, d2phi, yldc)
  Implicit None
!     Arguments
!     cauchy: cauchy stress
!     phi   : yield surface
!     dphi  : yield surface 1st derivative w.r.t cauchy stress
!     d2phi : The 2nd derivative w.r.t. cauchy stress
!     yldc  : yield surface components
  Integer ntens
  Parameter (ntens=3)
  Dimension cauchy(ntens), sdev(ntens), dphi(ntens), d2phi(ntens, ntens), yldc(9)
  Real *8 cauchy, phi, dphi, d2phi, yldc, hydro, sdev
!f2py intent(in)  cauchy,yldc
!f2py intent(out) phi,dphi,d2phi
  Call deviat(ntens, cauchy, sdev, hydro)
!     call yld2000_2d_dev(sdev,phi,dphi,d2phi,yldc)
  Call yld2000_2d_cauchy(cauchy, phi, dphi, d2phi, yldc)
  Return
End Subroutine yld2000_2d
!     yld2000_2d_dev is implementation of the yld2000-2d model as a
!     function of 'deviatoric' stress and yield function constants
!     stored in yldc.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Calculate deviatoric part of cauchy tensor
Subroutine deviat(ntens, cauchy, sdev, p)
!     Arguments
!     ntens   : Len of <cauchy> and <sdev>
!     cauchy  : cauchy stress (can be other tensors of interest)
!     sdev    : deviatoric stress
!     p       : hydrostatic pressure
  Implicit None
  Integer, Intent (In) :: ntens
  Dimension cauchy(ntens), sdev(ntens)
  Real *8, Intent (In) :: cauchy
  Real *8, Intent (Out) :: sdev, p
  If (ntens==3) Then
    Call deviat3(cauchy, sdev, p)
  Else If (ntens==6) Then
    Call deviat6(cauchy, sdev, p)
  Else
         !write(*,*) 'unexpected case given to deviat'
         !call exit(-1)
  End If
  Return
End Subroutine deviat
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Calculate deviator and hydrostatic pressure
Subroutine deviat6(s, sd, p)
!     Arguments
!     s   : cauchy stress under plane stress space (s11,s22,s12 with s33=0)
!     sd  : stress deviator
!     p   : hydrostatic pressure
  Implicit None
  Real *8, Intent (In) :: s(6)
  Real *8, Intent (Out) :: sd(6), p
  p = (s(1)+s(2)+s(3))/3D0
  sd(1) = s(1) - p
  sd(2) = s(2) - p
  sd(3) = s(3) - p
  sd(4) = s(4)
  sd(5) = s(5)
  sd(6) = s(6)
  Return
End Subroutine deviat6
!-----------------------------------------------------------------------
!     Given stress s3; plane stress condition with s(3) = 0.
Subroutine deviat3(s, sd, p)
!     Arguments
!     s   : cauchy stress under plane stress space (s11,s22,s12 with s33=0)
!     sd  : stress deviator
!     p   : hydrostatic pressure
  Implicit None
  Real *8, Intent (In) :: s(3)
  Real *8, Intent (Out) :: sd(3), p
  p = (s(1)+s(2))/3D0
  sd(1) = s(1) - p
  sd(2) = s(2) - p
  sd(3) = s(3) !! shear component s12                  
  Return
End Subroutine deviat3
!----------------------------------------------------------------------c
Subroutine yld2000_2d_cauchy(cauchy, psi, dpsi, d2psi, yldc)
!     psi is defiend as ((phi`+phi``)/2)^(1/q)
!       where phi` and phi`` are analogous to Hershey isotropic functions;
!             and q is the yield function exponent
!       note that phi` and phi`` are functions of principal values of
!       'linearly' transformed stress tensors.
!     psi is a homogeneous function of degree 1.

!     Arguments
!     cauchy: cauchy stress
!     psi   : yield surface
!     dpsi  : yield surface 1st derivative w.r.t cauchy stress
!     d2psi : 2nd derivative - (not included yet)
!     yldc  : yield surface constants
!             yldc(1:8) - alpha values, yldc(9): yield function exponent
  Implicit None
  Integer ntens
  Parameter (ntens=3)
  Dimension cauchy(ntens), sdev(ntens), dpsi(ntens), d2psi(ntens, ntens), c(2, ntens, ntens), x1(ntens), x2(ntens), dx_dsig(2, ntens, ntens), phis(2), dphis(2, ntens), yldc(9), alpha(8), aux2(ntens, ntens), aux1(ntens), xs(2, ntens), xp1(2), xp2(2), chis(2, 2), dphis_dchi(2, 2), x(ntens), dchi_dx(2, 2, ntens), aux23(2, ntens), aux233(2, ntens, ntens), dphis_aux(2, ntens)
  Real *8 cauchy, sdev, psi, dpsi, d2psi, c, x1, x2, a, phis, dphis, yldc, alpha, dx_dsig, aux2, aux1, hydro, xs, xp1, xp2, chis, dphis_dchi, x, dchi_dx, aux23, aux233, dphis_aux
  Integer i
!     locals controlling
  Integer imsg, ind, j, k
  Logical idiaw
!f2py intent(in)  cauchy,yldc
!f2py intent(out) psi,dpsi,d2psi
  imsg = 0
  idiaw = .False.
  alpha(:) = yldc(1:8)
  a = yldc(9) ! yield surface exponent                
  dphis(:, :) = 0D0
  Call alpha2lc(alpha, c, dx_dsig)
  aux2(:, :) = c(1, :, :)
  Call deviat(ntens, cauchy, sdev, hydro)
  Call calc_x_dev(sdev, aux2, x1) ! x1: linearly transformed stress (p
  xs(1, :) = x1(:)
  aux2(:, :) = c(2, :, :)
  Call calc_x_dev(sdev, aux2, x2) ! x2: linearly transformed stress (d
  xs(2, :) = x2(:)
  Call calc_princ(x1, xp1)
  chis(1, :) = xp1(:)
  Call calc_princ(x2, xp2)
  chis(2, :) = xp2(:)
  Call hershey(xp1, xp2, a, phis(1), phis(2), dphis_dchi)
  psi = (0.5D0*(phis(1)+phis(2)))**(1D0/a)

  aux2(:, :) = c(1, :, :)
  Call calc_dphi_dev(sdev, aux2, a, aux1, 0)
  dphis(1, :) = aux1(:)

  aux2(:, :) = c(2, :, :)
  Call calc_dphi_dev(sdev, aux2, a, aux1, 1)
  dphis(2, :) = aux1(:)
  dpsi(:) = 0D0
  Do i = 1, ntens
    dpsi(i) = 1D0/2D0/a*((phis(1)+phis(2))/2D0)**(1D0/a-1D0)*(dphis(1,i)+dphis(2,i))
  End Do
!     Just to check if chain-rule based yield surface derivative is
!     equivalent to the derivations shown in the yld2000-2d paper
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Eq. 53 in Ref [2]
  Do ind = 1, 2
    x = xs(ind, :)
    Call calc_dchi_dx(x, aux23, aux233)
    dchi_dx(ind, :, :) = aux23(:, :)
  End Do
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     1st derivatives
!     See = Eq. 43 in Ref [43]
  dphis_aux(:, :) = 0D0
  Do ind = 1, 2
    Do i = 1, ntens
      Do j = 1, 2
        Do k = 1, 3
          dphis_aux(ind, i) = dphis_aux(ind, i) + dphis_dchi(ind, j)*dchi_dx(ind, j, k)*dx_dsig(ind, k, i)
        End Do
      End Do
    End Do
  End Do

!      do 20 ind=1,2
!      do 20 i=1,ntens
!         if (dabs(dphis(ind,i)-dphis_aux(ind,i)).gt.1e-3) then
!         !       write(*,*)
!         !       write(*,'(4a8)')'ind','i','dp_sol','dp_x'
!         !       write(*,'(i8,i8,f8.3,f8.3)') ind,i,
!         !$           dphis(ind,i),dphis_aux(ind,i)
!c            call exit(-1)
!         endif
! 20   continue

!     2nd derivatives - on my to-do list
  d2psi(:, :) = 0D0
  Call yld2000_2d_hessian(psi, dpsi, cauchy, chis, xs, dx_dsig, yldc(9), d2psi)

  Return
End Subroutine yld2000_2d_cauchy
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Convert alpha parameters to c and l matrices
Subroutine alpha2lc(alpha, c, l)
!     Arguments
!     alpha : the 8 parameters
!     C     : C matrix
!     L     : L matrix
  Implicit None
  Dimension alpha(8), c(2, 3, 3), l(2, 3, 3), aux266(2, 6, 6)
  Real *8, Intent (In) :: alpha
  Real *8, Intent (Out) :: c, l
  Real *8 aux266, aux66(6, 6), aux33(3, 3)
  Call alpha2l(alpha, aux266)
  aux66(:, :) = aux266(1, :, :)
  Call reduce_basis(aux66, aux33)
  l(1, :, :) = aux33
  aux66(:, :) = aux266(2, :, :)
  Call reduce_basis(aux66, aux33)
  l(2, :, :) = aux33
  Call alpha2c(alpha, aux266)
  aux66(:, :) = aux266(1, :, :)
  Call reduce_basis(aux66, aux33)
  c(1, :, :) = aux33
  aux66(:, :) = aux266(2, :, :)
  Call reduce_basis(aux66, aux33)
  c(2, :, :) = aux33
  Return
End Subroutine alpha2lc
!-----------------------------------------------------------------------
Subroutine alpha2l(a, l)
!     Arguments
!     a : the 8 parameters (alpha)
!     l : l matrix in the dimension of (2,3,3)
  Implicit None
  Dimension a(8), l(2, 6, 6)
  Real *8, Intent (In) :: a
  Real *8, Intent (Out) :: l
  l(:, :, :) = 0D0
!     l`
  l(1, 1, 1) = 2*a(1)/3
  l(1, 1, 2) = -a(1)/3
  l(1, 2, 1) = -a(2)/3
  l(1, 2, 2) = 2*a(2)/3
  l(1, 6, 6) = a(7)
  l(1, 3, 1) = -l(1, 1, 1) - l(1, 2, 1)
  l(1, 3, 2) = -l(1, 1, 2) - l(1, 2, 2)
  l(1, 4, 4) = 1D0
  l(1, 5, 5) = 1D0
  l(1, 1, 3) = -a(1)/3
  l(1, 2, 3) = -a(2)/3
  l(1, 3, 3) = -l(1, 1, 3) - l(1, 2, 3)
!     l``
  l(2, 1, 1) = (8*a(5)-2*a(3)-2*a(6)+2*a(4))/9
  l(2, 1, 2) = (4*a(6)-4*a(4)-4*a(5)+a(3))/9
  l(2, 2, 1) = (4*a(3)-4*a(5)-4*a(4)+a(6))/9
  l(2, 2, 2) = (8*a(4)-2*a(6)-2*a(3)+2*a(5))/9
  l(2, 6, 6) = a(8)
  l(2, 3, 1) = -l(2, 1, 1) - l(2, 2, 1)
  l(2, 3, 2) = -l(2, 1, 2) - l(2, 2, 2)
  l(2, 4, 4) = 1D0
  l(2, 5, 5) = 1D0
  l(2, 1, 3) = (a(3)-4*a(5)+2*a(4)-2*a(6))/9
  l(2, 2, 3) = (2*a(5)-2*a(3)+a(6)-4*a(4))/9
  l(2, 3, 3) = -l(2, 1, 3) - l(2, 2, 3)
  Return
End Subroutine alpha2l
!-----------------------------------------------------------------------
!     Reduce (6,6) matrix to (3,3) by ignoring 4 and 5 components
Subroutine reduce_basis(a, b)
!     Arguments
!     a (6,6) matrix
!     b (3,3) matrix
  Implicit None
  Dimension a(6, 6), b(3, 3)
  Real *8 a, b
  Integer i, j, ii, jj
  Do i = 1, 3
    Do j = 1, 3
      If (i==3) ii = 6
      If (i/=3) ii = i
      If (j==3) jj = 6
      If (j/=3) jj = j
      b(i, j) = a(ii, jj)
    End Do
  End Do
  Return
End Subroutine reduce_basis
!-----------------------------------------------------------------------
Subroutine alpha2c(a, c)
!     Arguments
!     a : the 8 parameters (alpha)
!     c : c matrix in the dimension of (2,3,3)
  Implicit None
  Dimension a(8), c(2, 6, 6)
  Real *8 a, c
  c(:, :, :) = 0D0
!     c`
  c(1, 1, 1) = a(1)
  c(1, 1, 2) = 0D0
  c(1, 2, 1) = 0D0
  c(1, 2, 2) = a(2)
  c(1, 6, 6) = a(7)
  c(1, 3, 1) = -a(1)
  c(1, 3, 2) = -a(2)
  c(1, 4, 4) = 1D0
  c(1, 5, 5) = 1D0
!     c``
  c(2, 1, 1) = (4*a(5)-a(3))/3
  c(2, 1, 2) = 2*(a(6)-a(4))/3
  c(2, 2, 1) = 2*(a(3)-a(5))/3
  c(2, 2, 2) = (4*a(4)-a(6))/3
  c(2, 6, 6) = a(8)
  c(2, 3, 1) = -c(2, 1, 1) - c(2, 2, 1)
  c(2, 3, 2) = -c(2, 1, 2) - c(2, 2, 2)
  c(2, 4, 4) = 1D0
  c(2, 5, 5) = 1D0
  Return
End Subroutine alpha2c
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Convert deviatoric stress to linearly transformed X (eq 14 )
Subroutine calc_x_dev(sdev, c, x)
!     Arguments
!     sdev : deviatoric stress
!     c    : c matrix
!     x    : x tensor
  Implicit None
  Dimension sdev(3), c(3, 3), x(3)
  Real *8 c, x, sdev
  Call mult_array(c, sdev, 3, x)
  Return
End Subroutine calc_x_dev
!-----------------------------------------------------------------------
Subroutine calc_princ(s, xp)
!     Arguments
!     s: stress tensor in the in-plane strss space (sxx,syy,sxy)
!     xp: the two principal components
  Implicit None
  Dimension xp(2), s(3)
  Real *8 xx, yy, xy, f, s, xp
  xx = s(1)
  yy = s(2)
  xy = s(3)
  f = dsqrt((xx-yy)**2+4D0*xy**2)
  xp(1) = 0.5D0*(xx+yy+f)
  xp(2) = 0.5D0*(xx+yy-f)
  Return
End Subroutine calc_princ
!-----------------------------------------------------------------------
Subroutine hershey(xp1, xp2, a, phi1, phi2, dphis_dchi)
!     Arguments
!     x1, x2: two linearly transformed stresses
!     phi1, phi2: the two Hershey components
!     a: exponent
  Implicit None
  Dimension xp1(2), xp2(2), dphis_dchi(2, 2), chis(2, 2)
  Real *8, Intent (In) :: xp1, xp2, a
  Real *8, Intent (Out) :: phi1, phi2, dphis_dchi
  Real *8 tempval1, tempval2, sgn1, sgn2, chis
  phi1 = dabs(xp1(1)-xp1(2))**a
  phi2 = dabs(2D0*xp2(2)+xp2(1))**a + dabs(2D0*xp2(1)+xp2(2))**a

!     dphis/dchi
  chis(1, :) = xp1(:)
  chis(2, :) = xp2(:)

!     dphi1/dchi
  dphis_dchi(1, 1) = a*(chis(1,1)-chis(1,2))**(a-1D0)
  dphis_dchi(1, 2) = -dphis_dchi(1, 1)
!     dphi2/dchi
  tempval1 = 2D0*chis(2, 2) + chis(2, 1)
  tempval2 = 2D0*chis(2, 1) + chis(2, 2)
  sgn1 = sign(1D0, tempval1)
  sgn2 = sign(1D0, tempval2)
  dphis_dchi(2, 1) = a*dabs(tempval1)**(a-1D0)*sgn1 + 2D0*a*dabs(tempval2)**(a-1D0)*sgn2
  dphis_dchi(2, 2) = 2D0*a*dabs(tempval1)**(a-1D0)*sgn1 + a*dabs(tempval2)**(a-1D0)*sgn2
End Subroutine hershey
!-----------------------------------------------------------------------
Subroutine calc_dphi_dev(sdev, c, a, dphi, iopt)
!     Arguments
!     sdev  : deviatoric stress
!     c     : c matrix
!     a     : yield exponent
!     dphi  : dphi  (dphi` or dphi``)
!     iopt  : iopt (0: dphi`; 1: dphi``)
  Implicit None
  Dimension sdev(3), c(3, 3), dphi(3), x(3), xp(2), dphi_dp(2), dp_dx(2, 3), dphi_dx(3), l(3, 3)
  Real *8 c, delta, x, dphi, xp, dp_dx, dphi_dp, a, dphi_dx, l, sdev
  Integer iopt, i, j
  Call calc_x_dev(sdev, c, x)
  Call calc_delta(x, delta)
  Call calc_princ(x, xp)
  dphi_dx(:) = 0D0
  If (iopt==0) Call calc_a14(xp, a, dphi_dp) ! eq A1.4             
  If (iopt==1) Call calc_a18(xp, a, dphi_dp) 
! eq A1.8             
  If (delta/=0) Then
    Call calc_a15(x, delta, dp_dx) !        eq A1.3
! eq A1.9 = eq A1.5               
    Do j = 1, 3
      Do i = 1, 2
        dphi_dx(j) = dphi_dx(j) + dphi_dp(i)*dp_dx(i, j)
      End Do
    End Do
  Else
    If (iopt==0) Then
!        eq A1.6
      dphi_dx(:) = 0D0
    Else If (iopt==1) Then
!        eq A1.10
      dphi_dx(1) = dphi_dp(1)
      dphi_dx(2) = dphi_dp(2)
      dphi_dx(3) = 0D0
    End If
  End If
      ! dphi_dx
  dphi(:) = 0D0
  Call calc_l(c, l)
  Do j = 1, 3
    Do i = 1, 3
      dphi(j) = dphi(j) + dphi_dx(i)*l(i, j)
    End Do
  End Do
  Return
End Subroutine calc_dphi_dev
!-----------------------------------------------------------------------
!     eq (A1.2) delta used in principal stress calculation
Subroutine calc_delta(x, delta)
!     Arguments
!     x: linearly transformed stress using a particular c matrix (intent: in)
!     delta: to be calculated (intent: out)
  Implicit None
  Dimension x(3)
  Real *8 x, delta
  delta = (x(1)-x(2))**2 + 4D0*x(3)**2D0
  Return
End Subroutine calc_delta
!-----------------------------------------------------------------------
Subroutine calc_a14(xp, a, dphi_dp)
!     Arguments
!     xp      : principal components of linearly transformed stress
!     a       : yield function exponent
!     dphi_dp : round(phi) / round(xp)
  Implicit None
  Dimension xp(2), dphi_dp(2)
  Real *8 xp, dphi_dp, a
  dphi_dp(1) = a*(xp(1)-xp(2))**(a-1D0)
  dphi_dp(2) = -a*(xp(1)-xp(2))**(a-1D0)
End Subroutine calc_a14
!-----------------------------------------------------------------------
Subroutine calc_a18(xp, a, dphi_dp)
!     Arguments
!     xp : principal components of linearly transformed stress
!     a  : yield function exponent
!     dphi_dp : round(phi) / round(xp)
  Implicit None
  Dimension xp(2), dphi_dp(2)
  Real *8 xp, dphi_dp, sgn1, sgn2, a
  sgn1 = sign(1D0, 2D0*xp(2)+xp(1))
  sgn2 = sign(1D0, 2D0*xp(1)+xp(2))

  dphi_dp(1) = a*dabs(2D0*xp(2)+xp(1))**(a-1D0)*sgn1 + 2D0*a*dabs(2D0*xp(1)+xp(2))**(a-1D0)*sgn2
  dphi_dp(2) = 2D0*a*dabs(2D0*xp(2)+xp(1))**(a-1D0)*sgn1 + a*dabs(2D0*xp(1)+xp(2))**(a-1D0)*sgn2
End Subroutine calc_a18
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Subroutine calc_a15(x, delta, dp_dx)
!     Arguments
!     x     : linearly transformed stress
!     delta : delta used in the principal stress calculation
!     dp_dx : round (principal stress) / round (linearly transformed stress)
  Implicit None
  Dimension x(3), dp_dx(2, 3)
  Real *8 x, delta, dp_dx, deno
  deno = dsqrt(delta)
  dp_dx(1, 1) = 0.5D0*(1D0+(x(1)-x(2))/deno)
  dp_dx(1, 2) = 0.5D0*(1D0-(x(1)-x(2))/deno)
  dp_dx(1, 3) = 2.0D0*x(3)/deno

  dp_dx(2, 1) = 0.5D0*(1D0-(x(1)-x(2))/deno)
  dp_dx(2, 2) = 0.5D0*(1D0+(x(1)-x(2))/deno)
  dp_dx(2, 3) = -2.0D0*x(3)/deno
  Return
End Subroutine calc_a15
!-----------------------------------------------------------------------
Subroutine calc_l(c, l)
!     Arguments
!     c: c matrix
!     l: l matrix
!     intent(in) c
!     intent(out) l
  Implicit None
  Dimension c(3, 3), l(3, 3), t(3, 3)
  Real *8 c, l, t
  t(:, :) = 0D0
  t(1, 1) = 2D0/3D0
  t(1, 2) = -1D0/3D0
  t(2, 2) = 2D0/3D0
  t(2, 1) = -1D0/3D0
  t(3, 3) = 1D0
  Call mult_array2(c, t, 3, l)
  Return
End Subroutine calc_l
!-----------------------------------------------------------------------
!     Apply tensor inner dot for 2nd x 2nd
!     cij = aik x bkj
Subroutine mult_array2(aik, bkj, ntens, cij)
!     intent(int) aik,bkj,ntens
!     intent(out) cij
  Implicit None
  Integer i, j, k, ntens
  Real *8 aik(ntens, ntens), bkj(ntens, ntens), cij(ntens, ntens)
  cij(:, :) = 0.D0
  Do i = 1, ntens
    Do j = 1, ntens
      Do k = 1, ntens
        cij(i, j) = cij(i, j) + aik(i, k)*bkj(k, j)
      End Do
    End Do
  End Do
  Return
End Subroutine mult_array2
!------------------------------------------------------------------------
Subroutine calc_dchi_dx(x, dchi, d2chi)
  Dimension x(3), dchi(2, 3), d2chi(2, 3, 3)
  Real *8, Intent (In) :: x
  Real *8, Intent (Out) :: dchi, d2chi
  Real *8 delta, delta_sqrt, delta_sqrt3
  Dimension dd_dx(3), d2d_dx(3, 3)
  Real *8 dd_dx, d2d_dx
  Integer i, j
  delta = (x(1)-x(2))*(x(1)-x(2)) + 4D0*x(3)*x(3)
  delta_sqrt = delta**(-1D0/2D0)
  delta_sqrt3 = delta**(-3D0/2D0)

!     Eq 51 in Ref [2]
  dd_dx(1) = 2D0*(x(1)-x(2))
  dd_dx(2) = -dd_dx(1)
  dd_dx(3) = 8D0*x(3)

!     Eq 52 in Ref [2]
  d2d_dx(:, :) = 0D0
  d2d_dx(1, 1) = 2D0
  d2d_dx(1, 2) = -2D0
  d2d_dx(2, 1) = -2D0
  d2d_dx(2, 2) = 2D0
  d2d_dx(3, 3) = 8D0

!     Eq 53 in Ref [2]
  dchi(1, 1) = 1D0 + 0.5D0*dd_dx(1)*delta_sqrt
  dchi(1, 2) = 1D0 - 0.5D0*dd_dx(1)*delta_sqrt
  dchi(1, 3) = 0.5D0*dd_dx(3)*delta_sqrt
  dchi(1, 1:3) = dchi(1, 1:3)*0.5D0

  dchi(2, 1) = dchi(1, 2)
  dchi(2, 2) = dchi(1, 1)
  dchi(2, 3) = -dchi(1, 3)

!     Eq 54 in Ref [2]
  Do i = 1, 3
    Do j = 1, 3
      d2chi(1, i, j) = d2d_dx(i, j)*delta_sqrt - 0.5D0*dd_dx(i)*dd_dx(j)*delta_sqrt3
      d2chi(1, i, j) = d2chi(1, i, j)*0.25D0
    End Do
  End Do
!     Using Eq 53 in Ref[2], d2chi(2,:,:) is obtained as below

  d2chi(2, 1, :) = d2chi(1, 2, :)
      !d2chi(2,1,2) = d2chi(1,2,2)
      !d2chi(2,1,3) = d2chi(1,2,3)
  d2chi(2, 2, :) = d2chi(1, 1, :)
      !d2chi(2,2,2) = d2chi(1,1,2)
      !d2chi(2,2,3) = d2chi(1,1,3)
  d2chi(2, 3, :) = -d2chi(1, 3, :)
      !d2chi(2,3,2) =-d2chi(1,3,2)
      !d2chi(2,3,3) =-d2chi(1,3,3)

  Return
End Subroutine calc_dchi_dx
!-----------------------------------------------------------------------
Subroutine yld2000_2d_hessian(psi, dpsi, cauchy, chis, xs, dx_dsig, a, d2psi)
!     Arguments
!     psi    : yld2000-2d yield surface
!     dpsi   : yld2000-2d yield surface 1st derivatives
!     cauchy : cauchy stress
!     chis   : principal values of transformed stress
!     xs     : transformed stress
!     dx_dsig: the L matrix.
!     a      : yield surface exponent
!     d2psi  : the send derivative of yld2000-2d (to be calculated)
  Implicit None
  Dimension cauchy(3), chis(2, 2), xs(2, 3), dx_dsig(2, 3, 3), dphis_dchi(2, 2), d2phis_dchidchi(2, 2, 2), d2psi(3, 3), dpsi(3)
  Real *8, Intent (In) :: cauchy, chis, xs, dx_dsig, a, psi, dpsi
  Real *8, Intent (Out) :: d2psi
  Dimension dchi_dx(2, 2, 3), d2chi_dxdx(2, 2, 3, 3), x(3), aux23(2, 3), aux233(2, 3, 3), d2phis_dsigdsig(2, 3, 3)
  Real *8 tempval1, tempval2, tempval3, sgn1, sgn2, dchi_dx, d2chi_dxdx, aux23, aux233, d2phis_dchidchi, d2phis_dsigdsig, dphis_dchi, x
  Integer i, j, k, l, m, n, ind
!      write(3,*) cauchy,chis,xs,dx_dsig,a
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Eq. 46 in Ref [2] - calculate dphis_dchi
  dphis_dchi(1, 1) = a*(chis(1,1)-chis(1,2))**(a-1D0)
  dphis_dchi(1, 2) = -dphis_dchi(1, 1)
  tempval1 = 2D0*chis(2, 2) + chis(2, 1)
  tempval2 = 2D0*chis(2, 1) + chis(2, 2)
  sgn1 = sign(1D0, tempval1)
  sgn2 = sign(1D0, tempval2)
  dphis_dchi(2, 1) = a*dabs(tempval1)**(a-1D0)*sgn1 + 2D0*a*dabs(tempval2)**(a-1D0)*sgn2
  dphis_dchi(2, 2) = 2D0*a*dabs(tempval1)**(a-1D0)*sgn1 + a*dabs(tempval2)**(a-1D0)*sgn2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Eq. 47 in Ref [2] - calculate (d2phis)_(dchi,dchi)
  d2phis_dchidchi(1, 1, 1) = a*(a-1D0)*(chis(1,1)-chis(1,2))**(a-2D0)
  d2phis_dchidchi(1, 1, 2) = -d2phis_dchidchi(1, 1, 1)
  d2phis_dchidchi(1, 2, 1) = -d2phis_dchidchi(1, 1, 1)
  d2phis_dchidchi(1, 2, 2) = d2phis_dchidchi(1, 1, 1)
!     --
  tempval1 = dabs(2D0*chis(2,2)+chis(2,1))**(a-2D0)
  tempval2 = dabs(2D0*chis(2,1)+chis(2,2))**(a-2D0)
  tempval3 = a*(a-1D0)
  d2phis_dchidchi(2, 1, 1) = tempval1 + 4D0*tempval2
  d2phis_dchidchi(2, 1, 2) = 2D0*tempval1 + 2D0*tempval2
  d2phis_dchidchi(2, 2, 1) = d2phis_dchidchi(2, 1, 2)
  d2phis_dchidchi(2, 2, 2) = 4D0*tempval1 + tempval2
  d2phis_dchidchi(2, :, :) = tempval3*d2phis_dchidchi(2, :, :)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Eq. 53 in Ref [2]
  Do ind = 1, 2
    x(:) = xs(ind, :)
    Call calc_dchi_dx(x, aux23, aux233)
    dchi_dx(ind, :, :) = aux23(:, :)
    d2chi_dxdx(ind, :, :, :) = aux233(:, :, :)
  End Do
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Eq. 45 in Ref [2]
  d2phis_dsigdsig(:, :, :) = 0D0
  Do ind = 1, 2
    Do i = 1, 3
      Do j = 1, 3
        Do n = 1, 3
          Do m = 1, 2
            Do l = 1, 3
              Do k = 1, 2
                d2phis_dsigdsig(ind, i, j) = d2phis_dsigdsig(ind, i, j) + d2phis_dchidchi(ind, k, m)*dchi_dx(ind, k, l)*dx_dsig(ind, l, i)*dchi_dx(ind, m, n)*dx_dsig(ind, n, j)
              End Do
            End Do
          End Do
        End Do
      End Do
    End Do
    Do i = 1, 3
      Do j = 1, 3
        Do m = 1, 3
          Do l = 1, 3
            Do k = 1, 2
              d2phis_dsigdsig(ind, i, j) = d2phis_dsigdsig(ind, i, j) + dphis_dchi(ind, k)*d2chi_dxdx(ind, k, l, m)*dx_dsig(ind, l, i)*dx_dsig(ind, m, j)
            End Do
          End Do
        End Do
      End Do
    End Do
  End Do
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Eq. 44 in Ref [2]
  d2psi(:, :) = 0D0
  Do i = 1, 3
    Do j = 1, 3
      d2psi(i, j) = (1D0-a)/psi*dpsi(i)*dpsi(j) + psi**(1D0-a)/(2D0*a)*(d2phis_dsigdsig(1,i,j)+d2phis_dsigdsig(2,i,j))
    End Do
  End Do
  Return
End Subroutine yld2000_2d_hessian
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Subroutine update_elastic(idia, idiaw, iyld_law, ntens, nyldp, nstatv, ddsdde, cel, stran, stran_el_ns, stran_pl_ns, dstran, stress, eeq_ns, deq, yldp_ns, statev, kinc)
  Implicit None
  Integer idia, iyld_law, nstatv, ntens, nyldp, kinc
  Logical idiaw
  Dimension ddsdde(ntens, ntens), cel(ntens, ntens), stran(ntens), stran_el_ns(0:1, ntens), stran_pl_ns(0:1, ntens), dstran(ntens), stress(ntens), yldp_ns(0:1, nyldp), statev(nstatv), eeq_ns(0:1)
  Real *8 ddsdde, cel, stran, stran_el_ns, stran_pl_ns, dstran, stress, yldp_ns, statev, eeq_ns, deq
  Real *8 empa, gpa
  empa = 1D6
  gpa = 1D9
  deq = 0D0
!$$$  1. Save jacobian as elastic moduli
  ddsdde(:, :) = cel(:, :)
!$$$  2. Update strain.
  stran_el_ns(1, :) = stran_el_ns(0, :) + dstran(:)
!     call add_array(stran,dstran,ntens)
!$$$  3. Updates stress
  Call mult_array(ddsdde, stran_el_ns(1,:), ntens, stress)
!$$$  4. Update all other state varaiables
  eeq_ns(1) = eeq_ns(0) + deq
  stran_pl_ns(1, :) = stran_pl_ns(0, :)
  Call update_yldp(iyld_law, yldp_ns, nyldp, deq)
!$$   5. Store updated state variables to statev
  Call restore_statev(statev, nstatv, eeq_ns(1), stran_el_ns(1,:), stran_pl_ns(1,:), ntens, yldp_ns(1,:), nyldp, 1, .False., idia, .False., -1, -1, -1, -1D0, stress)
  Return
End Subroutine update_elastic
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Subroutine return_mapping(cel, spr, phi_n, eeq_n, dphi_n, dstran, stran_el, stran_pl, ntens, idiaw, idia, hrdp, nhrdp, hrdc, nhrdc, ihard_law, iyld_law, yldc, nyldc, yldp_ns, nyldp, & 
!          variables to be updated.
  snew, deeq, dstran_pl, dstran_el, statev, nstatv, ddsdde, failnr, & 
!          variables to tell about the current integration point and time
  kinc, noel, npt, time, clf1, clf2, clf3, clf4)
!-----------------------------------------------------------------------
!***  Arguments
!     Cel     : elastic moduli
!     spr     : predictor stress at k=0
!     phi_n   : yield surface at the given step n
!     eeq_n   : accumulative plastic equivalent strain at step n
!     dphi_n  : dphi at step n
!     dstran  : total incremental strain given between steps n and n+1
!     stran_el: total elastic strain at step n
!     stran_pl: total plastic strain at step n
!     ntens   : len of stress/strain tensor (also applied to Cel)
!     idiaw   : whether or not turn on the diagnostic streams
!     idia    : idia
!     hrdp    : hardening state variable arrays associated with hrdp
!               (e.g., equivalent plastic strain for isotropic hardening)
!     hrdc    : hardening constants (invariable)
!     nhrdp   : Len of hrdp
!     nhrdc   : Len of hrdc
!     ihard_law: hardening law (refer to uhard.f for more details)
!     iyld_law : yield surface choice
!     yldc     : yield surface constants
!     yldp_ns  : yield surface state parameters
!     nyldp    : The number of yield surface state parameters
!     snew     : stress at step n+1 estimated by return_mapping method
!     deeq     : delta equivalent plastic strain
!     dstran_pl: delta plastic strain between steps n and n+1
!     dstran_el: delta elastic strain between steps n and n+1
!     statev   : state variables
!     nstatv   : the number of state variables
!     ddsdde   : Jacobian matrix for step n -> step n+1
!     failnr   : flag to tell if return mapping method fails
!     kinc     : increment number
!     noel     : element number
!     npt      : integration point number
!     time     : time staps at steps n and n+1
!-----------------------------------------------------------------------
!***  Intents of Arguments
!     intent(in) Cel, spr, phi_n, eeq_n, dphi_n, dstran, stran_el,
!                stran_pl,ntens,idiaw, hrdp,nhrdp,hrdc,nhrdc,ihard_law
!     intent(out) -- states at n+1 step
!-----------------------------------------------------------------------
  use random_forest
  Implicit None
  type(RandomForestRegression) :: clf1, clf2, clf3, clf4
  Character *255 fndia
  Character *20 chr
  Integer ntens, mxnr, nhrdc, nhrdp, ihard_law, iyld_law, nyldc, nyldp, nstatv, kinc, noel, npt
  Parameter (mxnr=10)
!-----------------------------------------------------------------------
  Dimension spr(ntens), dphi_n(ntens), snew(ntens), spr_ks(mxnr, ntens), statev(nstatv), dstran(ntens), time(2), & 
    stran_el(ntens), dstran_el(ntens), dstran_el_ks(mxnr, ntens), stran_el_ks(mxnr, ntens), stran_pl(ntens), dstran_pl(ntens), dstran_pl_ks(mxnr, ntens), stran_pl_ks(mxnr, ntens), & 
    aux_n(ntens), cel(ntens, ntens), eeq_ks(mxnr), fo_ks(mxnr), fp_ks(mxnr), dlamb_ks(mxnr), dphi_ks(mxnr, ntens), d2phi_ks(mxnr, ntens, ntens), phi_ks(mxnr), dh_ks(mxnr), h_flow_ks(mxnr), & 
    hrdc(nhrdc), hrdp(nhrdp), yldc(nyldc), yldp_ns(0:1, nyldp), ddsdde(ntens, ntens)
!-----------------------------------------------------------------------
  Real *8 cel, spr, dphi_n, dstran, stran_el, dstran_el, dstran_el_ks, stran_el_k, stran_el_ks, stran_pl, dstran_pl, dstran_pl_ks, stran_pl_k, stran_pl_ks, yldc, yldp_ns, statev, snew, deeq, ddsdde, time

  Real *8 & ! eq stress at nr-step k, stress predic at nr-
    spr_ks
  Real *8 fo_ks, & ! Fobjective, Jacobian for NR           
    fp_ks
  Real *8 dlamb_ks, phi_n
  Real *8 dphi_ks, d2phi_ks
  Real *8 delta_eeq, eeq_n, aux_n, eeq_ks, empa, gpa
  Real *8 h_flow_ks, dh_ks, phi_ks, tolerance, tol_val
  Real *8 hrdc, hrdp
  Integer k, idia, imsg
  Parameter (tolerance=1D-4)
  Logical idiaw, ibreak, failnr
  failnr = .False. 
! in case NR fails                     
  If (ntens/=3) Then

         !call exit(-1)
  End If

  empa = 1D6
  gpa = 1D9

!***  non-zero initial might be very risky.
!     delta_eeq = (dstran(1)**2+dstran(2)**2+dstran(2)**2)/3.d0 !! initial guess
  delta_eeq = & ! initial guess on equivalent strain rat
    0D0
  dlamb_ks(1) = delta_eeq
  spr_ks(1, :) = spr(:) !! stress predictor                    
  dphi_ks(1, :) = dphi_n(:)
  phi_ks(1) = phi_n
  dstran_el_ks(1, :) = dstran(:)
  dstran_pl_ks(1, :) = & 
    0D0
!------------------------------------------------------------------------
!     iv. return mapping (loop over k)
!! or using delta_eeq ...               
  k = 1

!      idia=315 ! write to diagnostic file
!      idia=0   ! write to stdo
!      idia=7   ! write to Abaqus msg file

  ibreak = .False.
  Do While (k<=mxnr)

!        spr_ks(k,:)    ! predictor stress at current k
!        dphi_ks(k,:) ! yield normal at current k
    eeq_ks(k) = eeq_n + dlamb_ks(k) 
!***  Hardening state variable updates according to NR step k
! assumed plastic strain at cur
    hrdp(1) = eeq_ks(k)
!***  --------------------------------
    Call uhard(ihard_law, hrdp, nhrdp, hrdc, nhrdc, h_flow_ks(k), dh_ks(k), empa)

    If (k==1) tol_val = h_flow_ks(1)*tolerance

    stran_el_ks(k, :) = dstran_el_ks(k, :) + stran_el(:)
    stran_pl_ks(k, :) = dstran_pl_ks(k, :) + stran_pl(:)
!-----------------------------------------------------------------------
!        f   = yield - hardening             (objective function)
    fo_ks(k) = phi_ks(k) - h_flow_ks(k)
    If (abs(fo_ks(k))<=tol_val) Then
      Goto 100
    Else
!           Find Fp
!           ** Use values pertaining to n+1 step (assuming that current eeq_ks(k) is correct)
      Call calc_fp(dphi_ks(k,:), cel, dh_ks(k), ntens, fp_ks(k))
    End If
!------------------------------------------------------------------------
!         2.  Update the multiplier^(k+1)  (dlamb)
!             dlamb^(k+1) = dlamb^k - fo_ks(k)/fp_ks(k)
    dlamb_ks(k+1) = dlamb_ks(k) + fo_ks(k)/fp_ks(k)

!     new plastic strain increment
    dstran_pl_ks(k+1, :) = dlamb_ks(k+1)*dphi_ks(k, :)
!     new elastic strain increment?
    dstran_el_ks(k+1, :) = dstran(:) - dstran_pl_ks(k+1, :)
!     new plastic acc strain
    stran_pl_ks(k+1, :) = stran_pl(:) + dstran_pl_ks(k+1, :)
!        new elastic acc strain
!     Update dE^(el)^(k+1) and update the predictor stress.
    stran_el_ks(k+1, :) = stran_el(:) + dstran_el_ks(k+1, :)

!     find the new predictor stress for next NR step
!     Using  dE = dE^(el)^(k+1) + dlamb^(k+1),
    Call mult_array(cel, stran_el_ks(k+1,:), ntens, spr_ks(k+1,:))
!------------------------------------------------------------------------
!        3. Find normal of the updated predictor stress (s_(n+1)^(k+1))
    Call update_yldp(iyld_law, yldp_ns, nyldp, dlamb_ks(k+1))
    if (all(abs(spr_ks(k+1,:)) > 1000 .and. abs(spr_ks(k+1,:)) < 1000000000)) then
        ! 缩放数据
        call process_stress(spr_ks(k+1,:), phi_ks(k+1), dphi_ks(k+1,:), clf1, clf2, clf3, clf4)
    else
        call yld(iyld_law, yldp_ns(1,:), yldc, nyldp, nyldc, spr_ks(k+1,:), phi_ks(k+1), dphi_ks(k+1,:), d2phi_ks(k+1,:,:), ntens)
        print *, "调用了原函数"
    endif
    print *, "mapping"
    !Call yld(iyld_law, yldp_ns(1,:), yldc, nyldp, nyldc, spr_ks(k+1,:), phi_ks(k+1), dphi_ks(k+1,:), d2phi_ks(k+1,:,:), ntens)
    k = k + 1

  End Do
!     case when k exceeds mxnr
! end of do while loop for NR procedure                     
  failnr = .True.

  Return
!-----------------------------------------------------------------------
! to get out before updating state variables.              
!***  update state variables
  100 Continue
! successful NR run                                      
  Call restore_statev(statev, nstatv, eeq_n+dlamb_ks(k), stran_el_ks(k,:), stran_pl_ks(k,:), ntens, yldp_ns(1,:), nyldp, 1, .False., idia, .False., kinc, noel, npt, time(1), spr_ks(k+1,:))
!***  new stress
  snew(:) = spr_ks(k, :)
!***  Equivalent plastic strain increment
  deeq = dlamb_ks(k)
!***  Plastic strain increment
  dstran_pl(:) = dstran_pl_ks(k, :)
!***  Elastic strain increment
  dstran_el(:) = dstran_el_ks(k, :)
!***  new jacobian
  Call calc_epl_jacob(cel, dphi_ks(k,:), dh_ks(k), ntens, ddsdde, idia, idiaw)

!     if (idiaw) close(idia)
  Return
End Subroutine return_mapping
!-----------------------------------------------------------------------
!     Calculate fp using the below formula
!     fp  = r(s^eq_(n+1)^k)/r(s_(n+1)^k) : -C^el : r(s^eq_(n+1)^k / r(s_(n+1)^k) + H`)
!     fp = dphi_i C_ij dphi_j + H
Subroutine calc_fp(dphi, cel, dh, ntens, fp)
!     intent(in) dphi,Cel,dh,ntens
!     intent(out) fp
!     dphi: round(s^eq)/round(s)
!     Cel : elastic modulus
!     dh  : round(s^flow)/round(e^eq)
!     ntens: len of dphi (also applied to Cel)
!     fp   : slope
  Implicit None
  Integer ntens
  Dimension s(ntens), cel(ntens, ntens), dphi(ntens)
  Real *8 s, seq, cel, dphi, fp, dh
  Integer i, j
  fp = 0.D0
  Do i = 1, ntens
    Do j = 1, ntens
      fp = fp + dphi(i)*cel(i, j)*dphi(j) + dh
    End Do
  End Do
  Return
End Subroutine calc_fp
!-----------------------------------------------------------------------
!     Calculate elasto-plastic consistent tangent modulus
!     Following J. W. Yoon et. al., 1999, IJP 15, p35-67
Subroutine calc_epl_jacob(cel, dphi, dh, ntens, jacob, idia, idiaw)
!     Arguments
!-----------------------------------------------------------------------
!     Cel   : elastic constants
!     dphi  : 1st derivative of yield potential w.r.t. stress
!     dh    : dh/de^eq
!     ntens : Len of tensor
!     jacob : Jacobian
!     idia  : File unit
!     idiaw : Flag to write the stream
  Implicit None
  Integer ntens
  Dimension cel(ntens, ntens), dphi(ntens), jacob(ntens, ntens)
  Real *8 cel, dphi, jacob, dh

!     local varaiables
  Dimension a(ntens), b(ntens)
  Real *8 deno, a, b, gpa
  Integer i, j, k, idia
  Logical idiaw

  gpa = 1D9
  deno = 0D0
  a(:) = 0D0
  b(:) = 0D0
  jacob(:, :) = 0D0

  Do i = 1, ntens
    Do j = 1, ntens
      deno = deno + dphi(i)*cel(i, j)*dphi(j)
      a(i) = a(i) + cel(i, j)*dphi(j)
      b(i) = b(i) + cel(i, j)*dphi(j)
    End Do
  End Do
  deno = deno + dh
  Do i = 1, ntens
    Do j = 1, ntens
      jacob(i, j) = cel(i, j) - a(i)*b(j)/deno
    End Do
  End Do
!-----------------------------------------------------------------------
  Return
End Subroutine calc_epl_jacob
!-----------------------------------------------------------------------