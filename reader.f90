module mod_lammps_reader
    use sorting, only: argsort
    use iso_fortran_env, only: int64, dp => real64

    implicit none

    type :: metadata
        integer(int64) :: step, num_atoms
        integer :: num_chunks, num_columns

        ! boundary information
        integer :: triclinic, boundary_conditions(2, 3)
        logical :: is_triclinic
        real(dp) :: boundary(2, 3), & ! ((xmin,xmax),(ymin,ymax),(zmin,zmax))^T
                    angles(3)         ! if triclinic

        contains
            procedure :: read_metadata
    end type

    type :: lammps_reader
        integer :: funit

        real(dp), allocatable :: values(:,:)

        logical :: has_next_step, has_opened_file = .false.
        type(metadata) :: header, next_header

        contains
            procedure :: open_file_explicit
            procedure :: open_file_asterisk
            generic :: open_file => open_file_asterisk, open_file_explicit
            procedure :: read_next_header
            procedure :: read_step
            procedure :: reallocate_values
            final :: close_reader
    end type

    contains
        function read_metadata(self, funit) result(eof)
            class(metadata), intent(inout) :: self
            integer, intent(in) :: funit
            logical :: eof

            integer :: fstat

            not_eof: block

                read(funit, iostat=fstat) self%step, self%num_atoms, self%triclinic, &
                                          self%boundary_conditions, self%boundary

                if (fstat /= 0) exit not_eof

                self%is_triclinic = self%triclinic /= 0
                if (self%is_triclinic) then
                    read(funit, iostat=fstat) self%angles
                    if (fstat /= 0) exit not_eof
                end if

                read(funit, iostat=fstat) self%num_columns, self%num_chunks
                if (fstat /= 0) exit not_eof

                eof = .false.
                return
            end block not_eof

            eof = .true.
        end function

        subroutine open_file_explicit(self, filename)
            class(lammps_reader), intent(inout) :: self
            character(len=*), intent(in) :: filename

            integer :: fstat

            if (self%has_opened_file) close(self%funit)

            open(newunit=self%funit, file=filename, access="stream", &
                 action="read", iostat=fstat)

             if (fstat /= 0) then
                self%has_next_step = .false.
                return
            end if

            self%has_opened_file = .true.

            call self%read_next_header()
        end subroutine

        subroutine open_file_asterisk(self, filename, step)
            class(lammps_reader), intent(inout) :: self
            character(len=*), intent(in) :: filename
            integer, intent(in) :: step

            character(len=:), allocatable :: new_filename, prefix, suffix
            integer :: new_length, asterisk_location

            asterisk_location = index(filename, "*")

            prefix = filename(1:asterisk_location-1)
            suffix = filename(asterisk_location+1:)

            if (step == 0) then
                new_length = 1
            else
                new_length = floor(log10(1.0d0*step)) + 1
            end if

            new_length = new_length + len(prefix) + len(suffix)

            allocate(character(len=new_length) :: new_filename)

            write(new_filename, '(a,i0,a)') prefix, step, suffix

            write(*, '(a)') new_filename
            call self%open_file(new_filename)
        end subroutine

        function read_step(self) result(success)
            class(lammps_reader), intent(inout) :: self
            logical :: success

            integer :: chunk_number, chunk_start, chunk_end, &
                       values_in_chunk, atoms_in_chunk
            integer :: fstat

            not_eof: block
                if (.not. self%has_next_step) exit not_eof

                self%header = self%next_header
                call self%reallocate_values()

                chunk_end = 0
                do chunk_number = 1, self%header%num_chunks
                    read(self%funit, iostat=fstat) values_in_chunk
                    if (fstat /= 0) exit not_eof

                    atoms_in_chunk = values_in_chunk / self%header%num_columns
                    chunk_start = chunk_end + 1
                    chunk_end = chunk_start + atoms_in_chunk - 1

                    read(self%funit, iostat=fstat) self%values(:, chunk_start:chunk_end)
                    if (fstat /= 0) exit not_eof
                end do

                success = .true.
                call self%read_next_header
                return
            end block not_eof

            success = .false.
            self%has_next_step = .false.
        end function

        subroutine read_next_header(self)
            class(lammps_reader), intent(inout) :: self
            self%has_next_step = .not. self%next_header%read_metadata(self%funit)
        end subroutine

        subroutine reallocate_values(self)
            class(lammps_reader), intent(inout) :: self

            associate(num_columns => self%header%num_columns, &
                      num_atoms => self%header%num_atoms)

                if (allocated(self%values)) then
                    if (all(shape(self%values) == [int(num_columns, kind=int64), 1*num_atoms])) then
                        return
                    else
                        deallocate(self%values)
                    end if
                end if

                allocate(self%values(num_columns, num_atoms))

            end associate
        end subroutine

        subroutine close_reader(self)
            type(lammps_reader), intent(inout) :: self

            if (self%has_opened_file) close(self%funit)

            if (allocated(self%values)) deallocate(self%values)
        end subroutine
end module
