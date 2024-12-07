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

    subroutine initialize(clf, n_estimators, max_depth, min_samples_split, &
        & min_samples_leaf, min_split_gain, colsample_bytree, subsample, random_state)
        implicit none
        class(RandomForestRegression), intent(inout) :: clf
        integer, intent(in) :: n_estimators, max_depth, min_samples_split, min_samples_leaf, random_state
        real, intent(in) :: min_split_gain, subsample
        character(len=*), intent(in) :: colsample_bytree

        clf%n_estimators = n_estimators
        clf%max_depth = max_depth
        clf%min_samples_split = min_samples_split
        clf%min_samples_leaf = min_samples_leaf
        clf%min_split_gain = min_split_gain
        !clf%colsample_bytree = 3
        clf%subsample = subsample
        clf%random_state = random_state
        clf%colsample_bytree_str = colsample_bytree
    end subroutine initialize

    subroutine random_forest_fit(clf, dataset, targets)
        implicit none
        !include 'omp_lib.h'  ! ���� OpenMP ģ��
        class(RandomForestRegression), intent(inout) :: clf
        real, dimension(:,:), intent(in) :: dataset
        real, dimension(:), intent(in) :: targets
        integer :: i, num_samples, num_features, j, temp
        real, dimension(:,:), allocatable :: dataset_stage
        real, dimension(:), allocatable :: targets_stage
        integer, dimension(:), allocatable :: indices
        real, dimension(:), allocatable :: random_state_stages
        type(Tree), pointer :: trees
        integer, dimension(1) :: seed
        real :: rand
        integer, dimension(10000) :: temp_array
        real :: start_time, end_time
    
        !! ͳ�ƿ�ʼʱ��
        !start_time = omp_get_wtime()
        
        ! Initialize trees
        allocate(clf%trees(clf%n_estimators))
        num_samples = size(dataset, 1)
        num_features = size(dataset, 2)
    
        ! Set colsample_bytree based on the number of features
        if (trim(adjustl(clf%colsample_bytree_str)) == "sqrt") then
            clf%colsample_bytree = int(sqrt(real(num_features)))
        elseif (trim(adjustl(clf%colsample_bytree_str)) == "log2") then
            clf%colsample_bytree = int(log(real(num_features)))
        else
            clf%colsample_bytree = num_features
        end if
        
        ! Set random seed
        if (clf%random_state /= 0) then
            seed = clf%random_state
            call random_seed(put=seed)
        end if
        allocate(random_state_stages(clf%n_estimators))
        !�������������Ϊ����˳�������
        !������ô������10000����Զ����������������������������ı�����ı�����˳�����������Ӱ�죬���������ʧЧ
        do i = 1, 10000    
            temp_array(i) = i
        end do
        ! ��������
        do i = 10000, 2, -1
            call random_number(rand)
            j = 1 + int(rand * i)  ! ���ѡ��һ������
            temp = temp_array(i)
            temp_array(i) = temp_array(j)
            temp_array(j) = temp
        end do

        ! ѡ��ǰn_estimators����
        do i = 1, clf%n_estimators
            random_state_stages(i) = temp_array(i)
        end do
        
        ! �ϰ���������ɷ�ʽ1
        !do i = 1, clf%n_estimators
        !    random_state_stages(i) = i
        !end do
        !! ��������
        !do i = clf%n_estimators, 2, -1
        !    call random_number(rand)
        !    j = 1 + int(rand * i)  ! ���ѡ��һ������
        !    temp = random_state_stages(i)
        !    random_state_stages(i) = random_state_stages(j)
        !    random_state_stages(j) = temp
        !end do
        
        !! ��ʽ2�������������,ȫ��0-1֮��ʵ��
        !do i = 1, clf%n_estimators
        !    call random_number(random_state_stages(i))
        !end do
        
        !! ���й�����������ʽһ������޷��̶�����Ч�ʸ��ߡ�
        !!$omp parallel private(i, dataset_stage, targets_stage, trees, seed, rand, temp_array, j, temp)
        !!$omp do
        !do i = 1, clf%n_estimators
        !    ! Set random seed for each thread
        !    seed = random_state_stages(i)
        !    call random_seed(put=seed)
        !    ! bootstrap�зŻس�������ѵ��������
        !    call generate_dataset_stage(dataset, targets, clf, random_state_stages(i), dataset_stage, targets_stage)
        !
        !    ! Build decision tree
        !    call build_single_tree(trees, dataset_stage, targets_stage, 0, clf)
        !    clf%trees(i) = trees
        !end do
        !!$omp end do
        !!$omp end parallel

        ! ���й�����������ʽ����������Թ̶���Ч���Եͣ����ǽ����Ȼ�뵥�˼����в��
        !!$omp parallel do private(i, dataset_stage, targets_stage, trees) shared(clf, dataset, targets, random_state_stages) num_threads(4)
        do i = 1, clf%n_estimators
            !! Set random seed for each thread
            !seed = random_state_stages(i)
            !call random_seed(put=seed)
            ! bootstrap�зŻس�������ѵ��������
            call generate_dataset_stage(dataset, targets, clf, random_state_stages(i), dataset_stage, targets_stage)
            
            ! Build decision tree
            call build_single_tree(trees, dataset_stage, targets_stage, 0, clf)
            clf%trees(i) = trees
        end do
        !!$omp end parallel do
        
        !! ͳ�ƽ���ʱ��
        !end_time = omp_get_wtime()
        !
        !! ��ӡѵ��ʱ��
        !print *, "ģ��ѵ��ǽ��ʱ��:", end_time - start_time, "��"
        !! Save the model to a file
        !call save_model(clf, 'model.dat')
        
        !! Load the model from the file
        !call load_model(clf, 'model.dat')
    end subroutine random_forest_fit
    
    recursive subroutine save_tree(unit, trees)
        implicit none
        integer, intent(in) :: unit
        type(Tree), pointer :: trees

        if (associated(trees)) then
            write(unit, *) 1  ! �������ڵ����
            write(unit, *) trees%split_feature
            write(unit, *) trees%split_value
            write(unit, *) trees%leaf_value
            call save_tree(unit, trees%tree_left)
            call save_tree(unit, trees%tree_right)
        else
            write(unit, *) 0  ! �������ڵ㲻����
        end if
    end subroutine save_tree

    subroutine save_model(clf, filename)
        implicit none
        class(RandomForestRegression), intent(in) :: clf
        character(len=*), intent(in) :: filename
        integer :: i, unit
        type(Tree), pointer :: tree_ptr

        ! Open the file for writing
        open(newunit=unit, file=filename, status='replace', action='write')

        ! Write the number of trees
        !write(unit, *) clf%n_estimators
        write(unit, *) clf%n_estimators
        write(unit, *) clf%max_depth
        write(unit, *) clf%min_samples_split
        write(unit, *) clf%min_samples_leaf
        write(unit, *) clf%min_split_gain
        write(unit, *) clf%colsample_bytree
        write(unit, *) clf%subsample
        write(unit, *) clf%random_state
        write(unit, '(A)') clf%colsample_bytree_str

        ! Write each tree's data
        do i = 1, clf%n_estimators
            tree_ptr => clf%trees(i)  ! ʹ����ʱָ�����
            call save_tree(unit, tree_ptr)
        end do

        ! Close the file
        close(unit)
    end subroutine save_model
    
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
    
    subroutine generate_dataset_stage(dataset, targets, clf, random_state, dataset_stage, targets_stage)
        implicit none
        real, intent(in) :: dataset(:,:)
        real, intent(in) :: targets(:)
        class(RandomForestRegression), intent(inout) :: clf
        real, intent(in) :: random_state
        real, allocatable, intent(out) :: dataset_stage(:,:)
        real, allocatable, intent(out) :: targets_stage(:)
        integer :: num_samples, num_features, sample_size
        integer, allocatable :: indices(:)
        integer, allocatable :: unique_indices(:)
        integer, allocatable :: subcol_index(:)
        integer :: i
        integer, dimension(1) :: seed
    
        num_samples = size(dataset, 1)
        num_features = size(dataset, 2)
        sample_size = int(clf%subsample * num_samples)
        seed = random_state
        call random_seed(put=seed)
        call sample_with_replacement(num_samples, sample_size, indices)
        ! �޷Żصĳ�����Ҫ�������ظ�ֵҲ��������������ԡ�����Ҫ�޳��ظ����ݣ�����Ӱ�쵽�����������
        allocate(dataset_stage(sample_size, num_features))
        allocate(targets_stage(sample_size))
        do i = 1, sample_size
            dataset_stage(i,:) = dataset(indices(i), :)
            targets_stage(i) = targets(indices(i))
        end do
        !�޳��ظ�����
        !call get_unique_indices(indices, unique_indices)
        !call random_subset(num_features, clf%colsample_bytree, subcol_index)
        !call sort(subcol_index)
        !allocate(dataset_stage(size(unique_indices), size(subcol_index)))
        !allocate(targets_stage(size(unique_indices)))
        !do i = 1, size(unique_indices)
        !    dataset_stage(i,:) = dataset(unique_indices(i), subcol_index)
        !    targets_stage(i) = targets(unique_indices(i))
        !end do
    end subroutine generate_dataset_stage
    
    ! �зŻصĳ�ȡ���������س����������������clf%subsampleΪ1.������������Żᱻ�鵽����֮��
    subroutine sample_with_replacement(n, size, sample)
        implicit none
        integer, intent(in) :: n, size
        integer, allocatable, intent(out) :: sample(:)
        integer :: i
        real(8) :: rnd
    
        allocate(sample(size))
        do i = 1, size
            call random_number(rnd)
            sample(i) = 1 + int(rnd * n)  ! int������ȡ�����ض�С��
        end do
    end subroutine sample_with_replacement
    
    ! ��ȡ������ΨһԪ�ص�������������������
    subroutine get_unique_indices(indices, unique_indices)
        implicit none
        integer, intent(in) :: indices(:)
        integer, allocatable, intent(out) :: unique_indices(:)
        integer :: i, count, max_val
        logical, allocatable :: mask(:)
    
        max_val = maxval(indices)
        allocate(mask(max_val))
        mask = .false.
    
        do i = 1, size(indices)
            if (.not. mask(indices(i))) then
                mask(indices(i)) = .true.
            end if
        end do
    
        ! �ֶ�ͳ��ΨһԪ�ص�����
        count = 0
        do i = 1, max_val
            if (mask(i)) then
                count = count + 1
            end if
        end do
        allocate(unique_indices(count))
    
        count = 0
        do i = 1, size(indices)
            if (mask(indices(i))) then
                count = count + 1
                unique_indices(count) = indices(i)
                mask(indices(i)) = .false.
            end if
        end do
    end subroutine get_unique_indices
    
    subroutine random_subset(n, k, subset)
        implicit none
        integer, intent(in) :: n, k
        integer, allocatable, intent(out) :: subset(:)
        integer, allocatable :: indices(:)
        real(8) :: rnd
        integer :: i, j, temp

        ! Allocate memory for subset and indices arrays
        allocate(subset(k))
        allocate(indices(n))

        ! Initialize the indices array with values 1 to n
        do i = 1, n
            indices(i) = i
        end do

        ! Shuffle the indices array using Fisher-Yates algorithm
        do i = n, 2, -1
            call random_number(rnd)
            j = 1 + int(rnd * i)
            temp = indices(i)
            indices(i) = indices(j)
            indices(j) = temp
        end do

        ! Select the first k elements from the shuffled indices array
        subset = indices(1:k)

        ! Sort the subset array
        call insertion_sort(subset)
        
        ! Deallocate the indices array
        deallocate(indices)
    end subroutine random_subset
    
    subroutine insertion_sort(array)
        implicit none
        integer, intent(inout) :: array(:)
        integer :: i, j, key

        do i = 2, size(array)
            key = array(i)
            j = i - 1
            ! �ֿ��������
            do while (j >= 1)
                if (array(j) > key) then
                    array(j + 1) = array(j)
                    j = j - 1
                else
                    exit
                end if
            end do
            array(j + 1) = key
        end do
    end subroutine insertion_sort
    
    recursive subroutine build_single_tree(trees, dataset, targets, depth, clf)
        implicit none
        type(Tree), pointer :: trees
        real, dimension(:,:), intent(in) :: dataset
        real, dimension(:), intent(in) :: targets
        integer, intent(in) :: depth
        class(RandomForestRegression), intent(inout) :: clf
        integer :: best_split_feature
        real :: best_split_value, best_split_gain
        real, dimension(:,:), allocatable :: left_dataset, right_dataset
        real, dimension(:), allocatable :: left_targets, right_targets

        ! Create tree node
        allocate(trees)

        ! ����ýڵ�����ȫ��һ��/����С�ڷ���������С������������ѡȡ���ִ������������ֹ����.
        !����������û��ʵ�� len(np.unique(targets)) <= 1
        if (size(targets) <= 1 .or. size(targets) <= clf%min_samples_split) then
            trees%leaf_value = calc_leaf_value(targets)
            return
        end if
        
        if (depth < clf%max_depth) then
            call choose_best_feature(dataset, targets, best_split_feature, best_split_value, best_split_gain, clf)
            
            if (best_split_feature == -1) then
                trees%leaf_value = calc_leaf_value(targets)
                return
            end if

            call split_dataset(dataset, targets, best_split_feature, &
                & best_split_value, left_dataset, right_dataset, left_targets, right_targets)
            
            ! Check if we should stop splitting
            if (size(left_dataset, 1) <= clf%min_samples_leaf .or. size(right_dataset, 1) <= & 
                &clf%min_samples_leaf .or. best_split_gain <= clf%min_split_gain) then
                trees%leaf_value = calc_leaf_value(targets)
                return
            else
                !! Update feature importances
                !clf%feature_importances(best_split_feature) = clf%feature_importances(best_split_feature) + 1
            
                trees%split_feature = best_split_feature
                trees%split_value = best_split_value
                call build_single_tree(trees%tree_left, left_dataset, left_targets, depth + 1, clf)
                call build_single_tree(trees%tree_right, right_dataset, right_targets, depth + 1, clf)
            end if
        else
            trees%leaf_value = calc_leaf_value(targets)
        end if
    end subroutine build_single_tree
    
    real function calc_leaf_value(targets)
        implicit none
        real, dimension(:), intent(in) :: targets

        ! ѡ�����������ľ�ֵ��ΪҶ�ӽڵ�ȡֵ
        calc_leaf_value = sum(targets) / size(targets)
    end function calc_leaf_value

    subroutine choose_best_feature(dataset, targets, best_split_feature, best_split_value, best_split_gain, clf)
        implicit none
        real, dimension(:,:), intent(in) :: dataset
        real, dimension(:), intent(in) :: targets
        integer, intent(out) :: best_split_feature
        real, intent(out) :: best_split_value, best_split_gain
        class(RandomForestRegression), intent(in) :: clf

        integer :: num_features, feature, unique_count, i, k, j
        real, dimension(:), allocatable :: feature_data, unique_values, percentiles, temp_values
        real :: split_value, split_gain, left_sum, right_sum, left_mean, right_mean
        real, dimension(:), allocatable :: left_targets, right_targets
        logical, dimension(:), allocatable :: is_right_target
        integer :: left_size, right_size
        real :: start_time1, end_time1

        best_split_gain = huge(0.0)
        best_split_feature = -1
        best_split_value = 0.0

        ! ��ȡ�е�����
        num_features = size(dataset, 2)

        ! Ԥ����left_targets��right_targets���ڴ�
        allocate(left_targets(size(targets)))
        allocate(right_targets(size(targets)))
        allocate(is_right_target(size(targets)))
    
        ! Ԥ�ȷ���percentiles��temp_values�ڴ�
        allocate(percentiles(1000))
        allocate(temp_values(1000))

        ! ѭ������ÿһ��
        do feature = 1, num_features
            feature_data = dataset(:, feature)

            ! �ҵ�Ψһֵ������
            call cpu_time(start_time1)
            call find_unique_sorted(feature_data, unique_values)

            ! ����Ψһֵ������
            unique_count = size(unique_values)

            !! ���Ψһֵ����������1000��ѡ��1000���ٷ�λֵ��Ϊ��ѡ������ֵ,�������ȫ�����ݣ���if���ע�͵�
            !if (unique_count > 1000) then
            !    do i = 1, 1000
            !        percentiles(i) = (i - 1) * 100.0 / 999.0
            !    end do
            !    call percentile(unique_values, percentiles, temp_values)
            !    call move_alloc(temp_values, unique_values)
            !end if
            call cpu_time(end_time1)
            print *, "find_unique_sorted�Լ��ٷֳ�����ʱ:", end_time1 - start_time1, "��"

            ! ��ʼ��left_sum��right_sum
            left_sum = 0.0
            right_sum = sum(targets)

            ! ��ʼ��left_size��right_size
            left_size = 0
            right_size = size(targets)

            ! ��ʼ��right_targets
            right_targets = targets
            is_right_target = .true.

            ! �Կ��ܵķ�����ֵ��������棬ѡȡ����������ֵ
            call cpu_time(start_time1)
            do i = 1, size(unique_values) - 1
                split_value = unique_values(i)

                ! ����left_sum��right_sum���Լ�left_targets��right_targets
                do k = 1, size(targets)
                    if (dataset(k, feature) == split_value) then
                        left_size = left_size + 1
                        right_size = right_size - 1
                        left_sum = left_sum + targets(k)
                        right_sum = right_sum - targets(k)
                        left_targets(left_size) = targets(k)
                        is_right_target(k) = .false.
                    end if
                end do

                ! ����right_targets
                j = 0
                do k = 1, size(targets)
                    if (is_right_target(k)) then
                        j = j + 1
                        right_targets(j) = targets(k)
                    end if
                end do

                ! ����left_targets��right_targets��ƽ��ֵ
                if (left_size > 0) left_mean = left_sum / left_size
                if (right_size > 0) right_mean = right_sum / right_size

                ! ����ƽ�����
                split_gain = 0.0
                do k = 1, left_size
                    split_gain = split_gain + (left_targets(k) - left_mean) ** 2
                end do
                do k = 1, right_size
                    split_gain = split_gain + (right_targets(k) - right_mean) ** 2
                end do

                if (split_gain < best_split_gain) then
                    best_split_feature = feature
                    best_split_value = split_value
                    best_split_gain = split_gain
                end if
            end do
            call cpu_time(end_time1)
            print *, "���������ʱ:", end_time1 - start_time1, "��"
        end do

        ! �ͷ��ڴ�
        deallocate(left_targets)
        deallocate(right_targets)
        if (allocated(percentiles)) deallocate(percentiles)
        if (allocated(temp_values)) deallocate(temp_values)
    end subroutine choose_best_feature

    subroutine find_unique_sorted(arr, unique_values)
        implicit none
        real, dimension(:), intent(in) :: arr
        real, dimension(:), allocatable, intent(out) :: unique_values
        real, dimension(:), allocatable :: temp_values
        integer :: i, j, n, unique_count

        n = size(arr)
        allocate(temp_values(n))
        unique_count = 0

        do i = 1, n
            do j = 1, unique_count
                if (arr(i) == temp_values(j)) exit
            end do
            if (j > unique_count) then
                unique_count = unique_count + 1
                temp_values(unique_count) = arr(i)
            end if
        end do

        call quicksort(temp_values, 1, unique_count)
        ! ����unique_values�Ĵ�С��unique_count
        allocate(unique_values(unique_count))
        unique_values = temp_values(:unique_count)
        deallocate(temp_values)
        !call move_alloc(temp_values, unique_values) !���ڽ����ڴ棬�ܷ��㡣

    end subroutine find_unique_sorted

    recursive subroutine quicksort(arr, left, right)
        implicit none
        real, dimension(:), intent(inout) :: arr
        integer, intent(in) :: left, right
        integer :: i, j
        real :: pivot, temp

        if (left < right) then
            pivot = arr(left)
            i = left
            j = right

            do while (i < j)
                do while (arr(j) >= pivot .and. i < j)
                    j = j - 1
                end do
                if (i < j) then
                    arr(i) = arr(j)
                    i = i + 1
                end if
                do while (arr(i) <= pivot .and. i < j)
                    i = i + 1
                end do
                if (i < j) then
                    arr(j) = arr(i)
                    j = j - 1
                end if
            end do
            arr(i) = pivot
            call quicksort(arr, left, i - 1)
            call quicksort(arr, i + 1, right)
        end if
    end subroutine quicksort

    subroutine percentile(data, percentiles, result)
        implicit none
        real, dimension(:), intent(in) :: data, percentiles
        real, dimension(:), allocatable, intent(out) :: result
        integer :: i, n

        n = size(percentiles)
        allocate(result(n))

        do i = 1, n
            result(i) = data(int(percentiles(i) / 100.0 * (size(data) - 1)) + 1)
        end do
    end subroutine percentile
    
    function calc_r2(left_targets, right_targets) result(r2)
        implicit none
        real, dimension(:), intent(in) :: left_targets, right_targets
        real :: r2
        real :: mean
        integer :: i
    
        r2 = 0.0
    
        ! �������Ӽ���ƽ�����
        mean = sum(left_targets) / size(left_targets)
        do i = 1, size(left_targets)
            r2 = r2 + (left_targets(i) - mean) ** 2
        end do
    
        ! �������Ӽ���ƽ�����
        mean = sum(right_targets) / size(right_targets)
        do i = 1, size(right_targets)
            r2 = r2 + (right_targets(i) - mean) ** 2
        end do
    
        ! �����Ҫ������ȡ��ע�����´���������ƽ��ƽ�����
        ! if (size(left_targets) > 1) r2 = r2 / size(left_targets)
        ! if (size(right_targets) > 1) r2 = r2 / size(right_targets)
    end function calc_r2
    !function calc_r2(left_targets, right_targets) result(r2)
    !    implicit none
    !    real, dimension(:), intent(in) :: left_targets, right_targets
    !    real :: r2
    !    real :: mean_left, mean_right
    !    integer :: i, n_left, n_right
    !    real :: sum_left, sum_right, sum_sq_left, sum_sq_right
    !
    !    ! ��ʼ��
    !    n_left = size(left_targets)
    !    n_right = size(right_targets)
    !    sum_left = 0.0
    !    sum_right = 0.0
    !    sum_sq_left = 0.0
    !    sum_sq_right = 0.0
    !
    !    ! �������Ӽ������Ӽ��ĺͼ�ƽ����
    !    do i = 1, n_left
    !        sum_left = sum_left + left_targets(i)
    !        sum_sq_left = sum_sq_left + left_targets(i) * left_targets(i)
    !    end do
    !
    !    do i = 1, n_right
    !        sum_right = sum_right + right_targets(i)
    !        sum_sq_right = sum_sq_right + right_targets(i) * right_targets(i)
    !    end do
    !
    !    ! �����ֵ
    !    mean_left = sum_left / n_left
    !    mean_right = sum_right / n_right
    !
    !    ! ����ƽ�����
    !    r2 = sum_sq_left - sum_left * mean_left + sum_sq_right - sum_right * mean_right
    !end function calc_r2
    
    subroutine split_dataset(dataset, targets, split_feature, split_value, left_dataset, right_dataset, left_targets, right_targets)
        implicit none
        real, dimension(:,:), intent(in) :: dataset
        real, dimension(:), intent(in) :: targets
        integer, intent(in) :: split_feature
        real, intent(in) :: split_value
        real, dimension(:,:), allocatable, intent(out) :: left_dataset, right_dataset
        real, dimension(:), allocatable, intent(out) :: left_targets, right_targets

        integer :: i, left_count, right_count, n, m

        !��ȡ���ݼ�������������
        n = size(dataset, 1)
        m = size(dataset, 2)

        ! ���������Ӽ��Ĵ�С
        left_count = count(dataset(:, split_feature) <= split_value)
        right_count = count(dataset(:, split_feature) > split_value)

        ! ���������Ӽ��Ĵ�С
        allocate(left_dataset(left_count, m))
        allocate(right_dataset(right_count, m))
        allocate(left_targets(left_count))
        allocate(right_targets(right_count))

        ! ��ʼ��������
        left_count = 0
        right_count = 0

        ! ���ݷ��������ͷ�����ֵ�������ݼ���Ŀ��
        do i = 1, n
            if (dataset(i, split_feature) <= split_value) then
                left_count = left_count + 1
                left_dataset(left_count, :) = dataset(i, :)
                left_targets(left_count) = targets(i)
            else
                right_count = right_count + 1
                right_dataset(right_count, :) = dataset(i, :)
                right_targets(right_count) = targets(i)
            end if
        end do
    end subroutine split_dataset
    
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

program main
    use random_forest
    !use omp_lib
    implicit none 
    ! ��������
    real, allocatable :: data1(:,:)
    real, allocatable :: X(:,:)
    real, allocatable :: Y(:)
    real, dimension(:,:), allocatable :: X_train, X_test, X_train_scaled, X_test_scaled, X_scaled
    real, dimension(:), allocatable :: y_train, y_test, min_val, max_val
    type(RandomForestRegression) :: clf
    real, allocatable :: y_pred(:)
    real :: start_time, end_time, mre
    integer :: n_data

    ! �����ӳ����ȡ����
    call ReadData(data1)
    n_data = 10000
    allocate(X(n_data, 3))
    allocate(Y(n_data))
    X = data1(1:n_data,1:3)
    Y = data1(1:n_data,4)
    
    !real :: Y(10000, 1)
    !! �����ӳ����ȡ����
    !call ReadData(data1)
    !X = data1(1:10000,1:3)
    !Y = reshape(data1(1:10000,4), shape(Y))

    ! �������ݼ�
    call custom_train_test_split(X, Y, 0.1, 42, X_train, X_test, y_train, y_test)

    ! ��һ������
    call min_max_scaler(X_train, X_train_scaled, min_val, max_val)
    
    call cpu_time(start_time)
    allocate(X_test_scaled(size(X_test, 1), size(X_test, 2)))
    call apply_scaler(X_test, X_test_scaled, min_val, max_val)
    call cpu_time(end_time)
    print*, '�ֶ�����ʱ��: ', end_time - start_time
    !print*, '�ֶ�����ʱ��: ', start_time
    !print*, '�ֶ�����ʱ��: ', end_time

    ! ʹ����ͬ����Сֵ�����ֵ�������������ݼ�
    allocate(X_scaled(size(X, 1), size(X, 2)))
    call apply_scaler(X, X_scaled, min_val, max_val)

    ! ��ʼ�����ɭ�ֻع�ģ��
    call initialize(clf, 10, 5, 2, 2, 0.0, "none", 0.9, 66)
    
    ! ѵ��ģ��
    call cpu_time(start_time)
    call clf%fit(X_train_scaled, y_train)
    call cpu_time(end_time)
    print *, "ģ��ѵ��ʱ��:", end_time - start_time, "��,", "CPUʱ��,�ڲ��м�����ʧЧ"
    
    ! Ԥ��
    call cpu_time(start_time)
    call clf%predict(X_test_scaled, y_pred)
    call cpu_time(end_time)
    print *, "Ԥ��ʱ��:", end_time - start_time, "��"
    
    ! ����ƽ��������
    mre = mean_relative_error(y_test, y_pred)
    print *, "ƽ��������:", mre
    print *, "ƽ��������:", mre
    
contains    
    
    subroutine ReadData(data1)
        implicit none
        real, allocatable, intent(out) :: data1(:,:) ! �洢���ݵĶ�ά����
        integer :: i
        character(len=256) :: filename

        ! ����CSV�ļ���·��
        filename = 'C:\Users\zxm\Desktop\jianhau\pythoncode\cleaned_data.csv'

        ! ���ļ�
        open(105, file=filename)
        ! ��ȡ��������
        allocate(data1(2018024, 7))
        ! ��ȡ����
        do i = 1, 2018024
            read(105, *) data1(i, :)
        end do
        ! �ر��ļ�
        close(105)
    end subroutine ReadData

    subroutine custom_train_test_split(X, Y, test_size, random_state, X_train, X_test, y_train, y_test)
        implicit none
        real, intent(in) :: X(:,:)
        real, intent(in) :: Y(:)
        real, intent(in) :: test_size
        integer, intent(in) :: random_state
        real, dimension(:,:), allocatable, intent(out) :: X_train, X_test
        real, dimension(:), allocatable, intent(out) :: y_train, y_test
    
        integer :: total_samples, test_amount, i
        integer, dimension(:), allocatable :: all_indices, test_indices, train_indices
        integer, dimension(1) :: seed
    
        ! ��ȡ��������
        total_samples = size(X, 1)
    
        ! ���������������
        test_amount = int(total_samples * test_size)
    
        ! �������������
        allocate(all_indices(total_samples))
        allocate(test_indices(test_amount))
        allocate(train_indices(total_samples - test_amount))
        allocate(X_train(total_samples - test_amount, size(X, 2)))
        allocate(X_test(test_amount, size(X, 2)))
        allocate(y_train(total_samples - test_amount))
        allocate(y_test(test_amount))
    
        ! �����������
        seed = random_state
        call random_seed(put=seed)
    
        ! ��ʼ������
        all_indices = [(i, i=1,total_samples)]
    
        ! ��������
        call shuffle(all_indices)
    
        ! �ָ�����Ϊ���Լ���ѵ����
        test_indices = all_indices(1:test_amount)
        train_indices = all_indices(test_amount+1:total_samples)
    
        ! ���������ָ�����
        do i = 1, size(train_indices)
            X_train(i, :) = X(train_indices(i), :)
            y_train(i) = Y(train_indices(i))
        end do
        do i = 1, size(test_indices)
            X_test(i, :) = X(test_indices(i), :)
            y_test(i) = Y(test_indices(i))
        end do
        !implicit none
        !real, dimension(:,:), intent(in) :: X
        !real, dimension(:), intent(in) :: Y
        !real, intent(in) :: test_size
        !integer, intent(in) :: random_state
        !real, allocatable, intent(out) :: X_train(:,:), X_test(:,:)
        !real, allocatable, intent(out) :: y_train(:), y_test(:)
        !integer :: n, n_train, n_test, i, idx
        !integer, allocatable :: indices(:)
        !real :: r
        !
        !! �����������
        !call random_seed(put=random_state)
        !
        !n = size(X, 1)
        !n_test = int(n * test_size)
        !n_train = n - n_test
        !
        !allocate(indices(n))
        !indices = [(i, i=1, n)]
        !call shuffle(indices)
        !
        !allocate(X_train(n_train, size(X, 2)))
        !allocate(X_test(n_test, size(X, 2)))
        !allocate(y_train(n_train))
        !allocate(y_test(n_test))
        !
        !X_train = X(indices(1:n_train), :)
        !y_train = Y(indices(1:n_train))
        !X_test = X(indices(n_train+1:n), :)
        !y_test = Y(indices(n_train+1:n))
    end subroutine custom_train_test_split

    subroutine shuffle(indices)
        implicit none
        integer, dimension(:), intent(inout) :: indices
        integer :: i, j, temp, n
        real :: r

        n = size(indices)
        do i = n, 2, -1
            call random_number(r)
            j = int(r * i) + 1
            temp = indices(i)
            indices(i) = indices(j)
            indices(j) = temp
        end do
    end subroutine shuffle

    subroutine min_max_scaler(X, X_scaled, min_val, max_val)
        implicit none
        real, intent(in) :: X(:,:)
        real, dimension(:,:), allocatable, intent(out) :: X_scaled
        real, dimension(:), allocatable, intent(out) :: min_val, max_val
        integer :: i, j, n_features

        n_features = size(X, 2)
        allocate(min_val(n_features))
        allocate(max_val(n_features))
        allocate(X_scaled(size(X, 1), n_features))

        ! ������Сֵ�����ֵ
        do j = 1, n_features
            min_val(j) = minval(X(:, j))
            max_val(j) = maxval(X(:, j))
        end do

        ! ��������
        do j = 1, n_features
            X_scaled(:, j) = (X(:, j) - min_val(j)) / (max_val(j) - min_val(j))
        end do
        !do i = 1, size(X, 1)
        !    do j = 1, n_features
        !        X_scaled(i, j) = (X(i, j) - min_val(j)) / (max_val(j) - min_val(j))
        !    end do
        !end do
    end subroutine min_max_scaler

    subroutine apply_scaler(X, X_scaled, min_val, max_val)
        implicit none
        real, intent(in) :: X(:,:)
        real, intent(in) :: min_val(:)
        real, intent(in) :: max_val(:)
        real, dimension(:,:), allocatable, intent(out) :: X_scaled
        integer :: i, j, n_features

        allocate(X_scaled(size(X, 1), size(X, 2)))
        n_features = size(X, 2)
        ! ��������
        do j = 1, n_features
            X_scaled(:, j) = (X(:, j) - min_val(j)) / (max_val(j) - min_val(j))
        end do
        !do i = 1, size(X, 1)
        !    do j = 1, size(X, 2)
        !        X_scaled(i, j) = (X(i, j) - min_val(j)) / (max_val(j) - min_val(j))
        !    end do
        !end do
    end subroutine apply_scaler
    
    function mean_relative_error(y_true, y_pred) result(mre)
        implicit none
        real, dimension(:), intent(in) :: y_true
        real, dimension(:), intent(in) :: y_pred
        real :: mre
        real, dimension(size(y_true)) :: relative_errors
        integer :: i

        ! ����������
        do i = 1, size(y_true)
            relative_errors(i) = abs((y_true(i) - y_pred(i)) / y_true(i)) * 100.0
        end do

        ! ����ƽ��������
        mre = sum(relative_errors) / size(y_true)
    end function mean_relative_error
end program main


    
