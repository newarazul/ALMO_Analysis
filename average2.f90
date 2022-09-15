program average

implicit none

integer                   :: iunit
character(len=1024)       :: filename
character(len=3)          :: dummy
real(8)                   :: energies1(1:3), energies2(1:3), energies3(1:3)
!real(8)                   :: energies4(1:2), energies5(1:2)
real(8)                   :: sumlayer1, sumlayer2, sumlayer3
real(8)                   :: sumlayer4, sumlayer5, sumlayer6
!real(8)                   :: sum4(1:2), sum5(1:2)
!real(8)                   :: sum6(1:8), asymmetry(1:8)
integer                   :: l1count
real(8)                   :: bonds1(1:2)

l1count=0
iunit = 0
do iunit=1, 3500
write (filename, "(A,I6)") "asymmetry_", iunit*100 + 250000-100
print*, filename
open(11,file=filename)
read(11,*) energies1(1:3)
read(11,*) energies2(1:3)
!read(11,*) energies3(1:3)
!print*, energies1(1:3)

!if(energies1(1) .ne. 0) then
!sum1(1:2)=sum1(1:2) + energies1(1:2)
!l1count=l1count+1
!end if
!print*, l1count
sumlayer1=sumlayer1 + energies1(1)
sumlayer2=sumlayer2 + energies1(2)
sumlayer3=sumlayer3 + energies1(3)
sumlayer4=sumlayer4 + energies2(1)
sumlayer5=sumlayer5 + energies2(2)
sumlayer6=sumlayer6 + energies2(3)
!print*, iunit
close(11)
end do
print*, l1count
sumlayer1=sumlayer1/3500
sumlayer2=sumlayer2/3500
sumlayer3=sumlayer3/3500
sumlayer4=sumlayer4/3500
sumlayer5=sumlayer5/3500
sumlayer6=sumlayer6/3500
!sum3(1:2)=sum3(1:2)/3500
!sum6(1:8)=sum6(1:8)/3500
print*, 1, sumlayer1
print*, 2, sumlayer2
print*, 3, sumlayer3
print*, 4, sumlayer4
print*, 5, sumlayer5
print*, 6, sumlayer6
!calculate the interactions.
!bonds1(1)=2*sum2(2)/(2-sum6(2))
!bonds1(2)=2*sum2(2)-bonds1(1)
!print*, bonds1(1:2)



!print*, sum3(1:2)

open(11000, file="averaged.dat")
write(11000,*) sumlayer1, sumlayer2, sumlayer3
write(11000,*) sumlayer4, sumlayer5, sumlayer6
!write(11000,*) 3, sum3(1:2), sum6(3), sum6(7)
!write(11000,*) 4, sum4(1:2), sum6(4)
!write(11000,*) 5, sum5(1:2), sum6(5), sum6(8)
close(11000)
end program
