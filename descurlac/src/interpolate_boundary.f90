program interpolate_boundary
   implicit none

   real, allocatable :: a1(:), b1(:), a2(:), b2(:)
   integer           :: len, len2, num_int, ier, i, j, ic(4)
   character(100)    :: filename



   write (6, *) "length boundary file ?"
   read (5,*) len
   write (6,*) "number of points to interpolate ?"
   read (5,*) num_int
   len2 = len + num_int + num_int

   write (6, *) "input boundary file name ?"
   read (5,"(a)") filename
   open (1, file=ADJUSTL(TRIM(filename)), action="read", &
            form="formatted", iostat=ier)
   if (ier /= 0) then
      write(6,*) "error in opening ", trim(filename), ": ", ier
      stop
   end if

   allocate (a1(len), stat = ic(1))
   allocate (b1(len), stat = ic(2))
   allocate (a2(len2), stat = ic(3))
   allocate (b2(len2), stat = ic(4))
   if (ANY(ic /= 0)) then
      write(6,*) "allocation error: ", ic
      stop
   end if

   read (1,*,iostat=ier) (a1(i), b1(i), i=1,len)
   if (ier /= 0) then
      write(6,*) "error in reading ", trim(filename), ": ", ier
      stop
   end if

   close(1)

   write (6, *) "output boundary file name ?"
   read (5,"(a)") filename
   open (2, file=ADJUSTL(TRIM(filename)), action="write", &
            form="formatted", iostat=ier)
   if (ier /= 0) then
      write(6,*) "error in opening ", trim(filename), ": ", ier
      stop
   end if

   j = 0
   do i = 1, num_int
      a2(i+j) = a1(i)
      b2(i+j) = b1(i)
      j = j + 1
      a2(i+j) = a1(i) + ((a1(i+1) - a1(i)) / 2)
      b2(i+j) = b1(i) + ((b1(i+1) - b1(i)) / 2)
   end do
   do i = num_int+1, len-num_int
      a2(i+j) = a1(i)
      b2(i+j) = b1(i)
   end do
   do i = len-num_int+1, len
      a2(i+j) = a1(i-1) + ((a1(i) - a1(i-1)) / 2)
      b2(i+j) = b1(i-1) + ((b1(i) - b1(i-1)) / 2)
      j = j + 1
      a2(i+j) = a1(i)
      b2(i+j) = b1(i)
   end do

   write (2, "(2e16.8)") (a2(i), b2(i), i=1,len2)
   close(2)

end program
      
      
   
