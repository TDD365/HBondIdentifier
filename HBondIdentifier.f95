! A simple program to decide where the hydrogen bond goes
! Program by Diandong Tang 10/8/2017
!===============================================================================

! This subroutine is used to define user input:
! Starting Point and Atom number;
! Hydrogen Bond Radius(A);
! Maximum water bridge length(A);
! pdb file name.
module para

    ! Declare the parameters
    integer :: SP, Natom, EP, maxcyc, WBL
    real*8  :: HBR
    character*200 :: pdb, OType, WType
    logical :: alive

    contains

subroutine makepara()

    implicit none

    ! Internal control
    character*200 :: userinput

    ! User interface
70  continue
    write(*,*) "==========================================="
    write(*,*) "Setting Parameters."
    write(*,*) "Please choose from menu."
    write(*,*) "==========================================="
    write(*,"(' 1 - Start Point (Mandatory)               |', i6)") SP
    write(*,"(' 2 - End Point (Mandatory)                 |', i6)") EP
    write(*,"(' 3 - Atom Number (Mandatory)               |', i6)") Natom
    write(*,"(' 4 - PDB File Name (Mandatory)             |', a)") trim(adjustl(pdb))
    write(*,"(' 5 - Hydrogen bond Range in A              |', f16.9)") HBR
    write(*,"(' 6 - Water Bridge Length                   |', i3)") WBL
    write(*,"(' 7 - O Type in PDB                         |', a)") trim(adjustl(OType))
    write(*,"(' 8 - Water Type in PDB                     |', a)") trim(adjustl(WType))
    write(*,"(' 9 - Maximum Tries                         |', i10)") maxcyc
    write(*,*) "                                          |"
    write(*,*) "0 - Finished, Ready to Go.                |"
    write(*,*) "==========================================="
    read(*,*) userinput

    if (userinput=='1') then
        write(*,*) "Start Point:"
        read(*,*) SP
        call SYSTEM('cls')
        goto 70
    elseif (userinput=='2') then
        write(*,*) "End Point:"
        read(*,*) EP
        call SYSTEM('cls')
        goto 70
    elseif (userinput=='3') then
        write(*,*) "Atom Number:"
        read(*,*) Natom
        call SYSTEM('cls')
        goto 70
    elseif (userinput=='4') then
        write(*,*) "Drag in PDB file:"
        read(*,*) pdb
        pdb=trim(adjustl(pdb))
        call SYSTEM('cls')
        goto 70
    elseif (userinput=='5') then
        write(*,*) "Hydrogen bond Range in A:"
        read(*,*) HBR
        call SYSTEM('cls')
        goto 70
    elseif (userinput=='6') then
        write(*,*) "Water Bridge Length:"
        read(*,*) WBL
        call SYSTEM('cls')
        goto 70
    elseif (userinput=='7') then
        write(*,*) "O Type in PDB:"
        read(*,*) OType
        OType=trim(adjustl(OType))
        call SYSTEM('cls')
        goto 70
    elseif (userinput=='8') then
        write(*,*) "Water Type in PDB:"
        read(*,*) WType
        WType=trim(adjustl(WType))
        call SYSTEM('cls')
        goto 70
    elseif (userinput=='9') then
        write(*,*) "Maximum Tries :"
        read(*,*) maxcyc
        call SYSTEM('cls')
        goto 70
    elseif (userinput=='0') then
        continue
    else
        write(*,*) "Cannot read input. Please input a number."
        call SYSTEM("cls")
        goto 70
    end if

    ! Check parameters
    alive = .False.
    inquire(file=pdb,exist=alive)
    if (alive) then
        continue
    else
        write(*,*) "No PDB file found. Please check the input."
    end if

    write(*,*) "==========================================="
    write(*,*) "Parameters identified."
    write(*,*) "Press anykey to continue."
    write(*,*) "==========================================="
    read(*,*)

end subroutine

end module para

program HBondIdentifier

!   Declare Modules
    use para

!   Implicit the rubbish information
    implicit character*200 (r)


! Variables
!==================================

    ! Declare Internal Command
    integer :: can_we_go
    character*200 :: menuinput

    ! Declare Variables
    character*200,allocatable :: AtomTyp(:), SegID(:)
    real*8,allocatable :: coord(:,:)
    integer,allocatable :: AtomNum(:), SegNum(:)

    ! Declare Dummy Variables
    integer,allocatable :: route(:), possto(:), route_all(:,:)!, route_fin(:,:)
    integer :: h_from, h_go, possnum
    real*8 :: hlength, randx
    integer :: i,j,k,l,roughcount,precisecount
    logical :: judge,found

    ! The defaults
    SP = 175
    EP = 173
    Natom = 15399
    HBR = 3.7
    WBL = 3
    pdb = 'D:\Programs\HbondIdentifier\bin\Debug\test.pdb'
    OType = 'OH2'
    WType = 'TIP3'
    maxcyc = 1000


! Menu Control
!==================================

    ! Start Menu
    call SYSTEM("cls")
    call banner()

    ! Start internal control
    can_we_go = 0
    do while (can_we_go < 1)

    ! Read user input
        read(*,*) menuinput
        if (menuinput == "1") then
            call SYSTEM("cls")
            call makepara()
            call SYSTEM("cls")
            call banner()
            Write(*,*) "Note: parameter changed."
        elseif (menuinput == "2") then
            call SYSTEM("cls")
            call banner()
            can_we_go = can_we_go + 1
        elseif (menuinput == "3") then
            stop
        else
            call SYSTEM("cls")
            call banner()
            Write(*,*) "Man, Try Enter 1 or 2."
            Write(*,*)
        end if

    end do

    ! Initialize the Paras
    allocate(AtomTyp(Natom))
    allocate(SegID(Natom))
    allocate(AtomNum(Natom))
    allocate(SegNum(Natom))
    allocate(coord(Natom,3))
    allocate(route(WBL+1))
    allocate(possto(100))

    found=.false.
    route = 0
    possto = 0
    roughcount=0
    precisecount=0

    write(*,*) "==========================================="
    write(*,*) "Parameters initialized."

    ! Output warnings
    if (HBR > 5.0) write(*,*) "Warning! Large HBond range found."
    if (WBL > 10) write(*,*) "Warning! Long Water Bridge found."

    ! Check pdb again
    alive = .False.
    inquire(file=pdb,exist=alive)
    if (alive) then
        continue
    else
        write(*,*) ""
        write(*,*) "Error - No PDB file found."
        write(*,*) "==========================================="
        write(*,*) "----Hydrogen Bond Identifier End *Error----"
        write(*,*) "==========================================="
        stop
    end if
    write(*,*) "==========================================="
    write(*,*) ""
    write(*,*) "Start Process."
    write(*,*) "==========================================="

! Start Operation, everything should be settled now
!==================================

    ! Read PDB file
    open(10,file=pdb,status='old')

    ! Get the Heading
    call loclabel(10,"ATOM ")

    write(*,*) "Collecting Atom Information."
    ! Get info
    do i = 1, Natom
        read(10,*) r1,AtomNum(i),AtomTyp(i),SegID(i),SegNum(i),coord(i,1),coord(i,2),coord(i,3),r2,r3,r4
    end do
!    write(*,"(i6, ' Atoms Read.')") Natom
    write(*,*) "==========================================="

    ! Organize the info
    do i = 1, Natom
        AtomTyp(i)=trim(adjustl(AtomTyp(i)))
        SegID(i)=trim(adjustl(SegID(i)))
    end do

    ! Bye-bye
    close(10)

    ! Start the Loop, i is total tries, if tries expired, the loop ends

    ! Open the output
    open(20,file="route_raw.txt",status='unknown')
    write(*,*) "==========================================="
    write(*,*) "Calculating possible path. Output enabled."
    write(*,*) "==========================================="

    write(*,*) "==========================================="
    write(*,*) "Analyzing Trajectory, Please Wait..."
    do i = 1, maxcyc
        !write(*,"(2a,' Analyzing Trajectory ', i9)") char(13),i
        ! Start the route. j is number of route.
        ! Because route(1)=SP, j start from 2

        ! Set up the Route
        h_from = SP
        h_go = SP
        route=0
        route(1)=SegNum(SP)

        j=2
        do while (j <= WBL+1)
!            write(*,*) j
            ! Confirm how many possible accepters
            ! First, find waters for the next step, k is atom investigated, l is water number
            l=0
            do k = 1, Natom
                hlength=0
                call distance(h_from,k,coord,hlength)
                call dojudge(k,route,judge)

!                write(*,*) k,trim(adjustl(AtomTyp(k))),SegNum(k)
!                write(*,*) hlength,k,trim(adjustl(SegID(k)))
!                write(*,*) HBR,h_from,trim(adjustl(WType))
!                write(*,*) (.not. judge)
!                if (hlength < HBR .and. k /= h_from .and. (.not. judge) .and. &
!                &trim(adjustl(SegID(k))) == trim(adjustl(WType)) .and. trim(adjustl(AtomTyp(k))) == trim(adjustl(OType))) stop

                ! When the destination is reached, the k loop ends
                if (hlength < HBR .and. j>2 .and. AtomNum(k) == EP) then
                    l=Natom+1
                    possto(1)=k
                    !write(*,"(' Destination Found at Residue ', i6)") SegNum(k)
                    found=.true.
                    exit
                ! When proper water is found, count the waters
                elseif (hlength < HBR .and. k /= h_from .and. (.not. judge) .and. &
                &trim(adjustl(SegID(k))) == trim(adjustl(WType)) .and. trim(adjustl(AtomTyp(k))) == trim(adjustl(OType))) then
                    !write(*,"(' Water Found at Residue ', i6)") SegNum(k)
                    l=l+1
                    possto(l)=k
                ! If not proper, jump to next atom
                else
                    l=l
                end if
            end do

            ! Assign the waters
            possnum=l
!            write(*,*) possto

            ! Get the list of possible destinations
            ! When no waters found, this route ends
            if (possnum==0) then
                write(*,*) "==========================================="
                write(*,*) "==========================================="
                write(*,*) "No proper Path found --- Dead End."
                write(*,*) "==========================================="
                ! It's all over but the crying
                write(*,*) ""
                write(*,*) "==========================================="
                write(*,*) "----Hydrogen Bond Identifier End Normal----"
                write(*,*) "==========================================="
                close(20)
                stop
            ! When destination is reached, the j loop ends
            elseif (possnum==Natom+1) then
                h_go = possto(1)
                route(j)=SegNum(EP)
                ! Now we have the route, we will print the route
                write(20,*) (route(j),j=1,WBL+1)
            ! When there are options, use random number to decide which to go
            else
                call init_random_seed()
                call random_number(randx)
!                write(*,*) randx
                h_go = possto(int(randx*possnum)+1)
!                write(*,*) h_go
                ! Assign the destination
                route(j)=SegNum(h_go)
                h_from=h_go
                h_go=0
            end if
            j=j+1
!            write(*,*) j
        end do
    end do

    ! All done, close the output
    close(20)
    write(*,*) "==========================================="

    if (found) then
        write(*,*) "==========================================="
        write(*,*) "Path found. Organizing the output."
        write(*,*) "==========================================="
    else
        write(*,*) "==========================================="
        write(*,*) "No proper Path found --- End of Loop."
        write(*,*) "==========================================="
        ! It's all over but the crying
        write(*,*) ""
        write(*,*) "==========================================="
        write(*,*) "----Hydrogen Bond Identifier End Normal----"
        write(*,*) "==========================================="
        stop
    end if

    ! Start to organize the output
!    allocate(route_all(roughcount,WBL+1))
!    open(20,file='route_raw.txt',status='unknown')
!    do i=1,roughcount
!        read(20,*) (route_all(i,j),j=1,WBL+1)
!    end do

    ! Call bash to output cls path
!    call SYSTEM("awk '!seen[$0]++' route_raw > route")

    write(*,*) "==========================================="
    write(*,*) "Formatted. Please check the output."
    write(*,*) "==========================================="
    write(*,*) ""
    write(*,*) "==========================================="
    write(*,*) "++++Hydrogen Bond Identifier End Normal++++"
    write(*,*) "==========================================="
    read(*,*)

!   End
    stop

Contains
!===============================================================================

! Banner, typical in my programs.
subroutine banner()

    implicit none

    write(*,*) "==========================================="
    write(*,*) "                                           "
    write(*,*) "    Hydrogen Bond Identifier for CHARMM  DEV0.1  "
    write(*,*) "                                           "
    write(*,*) "          Code By Diandong Tang            "
    write(*,*) "        BEIJING NORMAL UNIVERSITY          "
    write(*,*) "         TangDD@mail.bnu.edu.cn            "
    write(*,*) "                                           "
    write(*,*) "          Waiting For USER INPUT           "
    write(*,*) "                                           "
    write(*,*) "         ========================          "
    write(*,*) "          1. Start Configuration           "
    write(*,*) "          2. Start Calculation             "
    write(*,*) "          3. Exit                          "
    write(*,*) "         ========================          "
    write(*,*) "                                           "
    write(*,*) "==========================================="
!    write(*,*) "P.S. If you don't know how to use the code,"
!    write(*,*) "     just ask me! Let's fuck the manual!   "
    write(*,*) ""

end subroutine banner

!===============================================================================

! These are tools used in this program.

! Locate-label, a function to locate certain keywords.
subroutine loclabel(fileid,label,ifound,irewind)
!find the line where the label first appears in fileid
!return ifound=1 if found the label, else return 0
!default is rewind, if irewind=0 will not rewind
!if current line already has the label, calling this subroutine will do nothing
    integer fileid,ierror
    integer,optional :: ifound,irewind
    character*200 c200
    character(len=*) label
    if ((.not.present(irewind)).or.(present(irewind).and.irewind==1)) rewind(fileid)
    do while(.true.)
        read(fileid,"(a)",iostat=ierror) c200
        if (index(c200,label)/=0) then
            backspace(fileid)
        if (present(ifound)) ifound=1 !found result
            return
        end if
    if (ierror/=0) exit
    end do
    if (present(ifound)) ifound=0
end subroutine loclabel


! Get the distance between atoms.
subroutine distance(f,t,coord,hlength)

    use para

    implicit none

    integer,intent(in) :: f,t
    real*8,intent(in) :: coord(Natom,3)
    real*8,intent(out) :: hlength

    hlength = ((coord(f,1)-coord(t,1))**2+(coord(f,2)-coord(t,2))**2+(coord(f,3)-coord(t,3))**2)**0.5

end subroutine distance


! Judge if a Subject is in Array
subroutine dojudge(sub,array,judge)

    use para

    implicit none

    integer,intent(in) :: sub
    integer,intent(in) :: array(WBL+1)
    logical,intent(out) :: judge
    integer :: i

!    write(*,*) Natom
    judge=.false.
!    write(*,*) array
!    write(*,*) sub
    do i = 1, WBL+1
        if (array(i)==sub) then
            judge=.true.
            exit
        else
            judge=.false.
        end if
    end do

end subroutine

subroutine init_random_seed()

      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)
end subroutine
!===============================================================================

end program HBondIdentifier
!===============================================================================
