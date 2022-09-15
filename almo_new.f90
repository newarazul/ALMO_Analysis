program almo_analysis

!!!!use this script to calculate the average donor and acceptor interactions energies for the chose layer distribution. the script
!calls the proximity file to import to divide the molecules into the layers

implicit none
logical                 :: connection

real(8)                 :: L1acceptorsum, L2acceptorsum, L3acceptorsum
real(8)                 :: L1donorsum, L2donorsum, L3donorsum
real(8)                 :: L3donor(1:5,1:384), L3acceptor(1:5,1:384)
real(8)                 :: L1strength, L2strength, L3strength
real(8)                 :: L2donor(1:5,1:384), L2acceptor(1:5,1:384)
real(8)                 :: L1donor(1:5,1:384), L1acceptor(1:5,1:384)
real(8)                 :: acceptor_label(1:5,1:384), donor_label(1:5,1:384)
real(8)                 :: lowest_acceptor(1:5,1:384), lowest_donor(1:5,1:384)
real(8)                 :: sumaverage1(1:2), sumaverage2(1:2), sumaverage3(1:2)
real(8)                 :: sumdonor_acceptor1(1:2), sumdonor_acceptor2(1:2)
real(8)                 :: sumdonor_acceptor3(1:2)
real(8)                 :: donor_acceptor_energy(1:2,1:384)
real(8)                 :: proximity(1:2)
real(8)                 :: L1asymmetry_donor,L2asymmetry_donor,L3asymmetry_donor
real(8)                 :: L1asymmetry_donor_total, L2asymmetry_donor_total, L3asymmetry_donor_total
real(8)                 :: L1asymmetry_acceptor,L2asymmetry_acceptor,L3asymmetry_acceptor
real(8)                 :: L1asymmetry_acceptor_total, L2asymmetry_acceptor_total, L3asymmetry_acceptor_total

integer                 :: layer1_index, layer2_index, layer3_index
integer                 :: layer1(1:384), layer2(1:384), layer3(1:384)


integer                 :: icount, bcount, ccount, lcount, l1count, lxcount


character(len=6)        :: chardum
character(len=6)        :: tmp



!.............................read proximity file.......................!
open(50, FILE ="Proximity.str")
do icount = 1, 13
        read(50,*) 
end do
layer1_index=0
layer2_index=0
layer3_index=0
do icount=1, 384
    read(50,*) chardum, proximity(1:2)
    if (proximity(1) .ge. -2.5 .and. proximity(1) .lt. 0.5) then
      layer1_index=layer1_index +1
      layer1(layer1_index) = icount
      print*, layer1(layer1_index)
    else if (proximity(2) .ge. -2.5 .and. proximity(2) .lt. 0.5) then
      layer1_index=layer1_index +1
      layer1(layer1_index) = icount
      print*, layer1(layer1_index)
    else if (proximity(1) .ge. 3.5 .and. proximity(1) .lt. 6.5) then
      layer2_index=layer2_index +1
      layer2(layer2_index) = icount
    else if (proximity(2) .ge. 3.5 .and. proximity(2).lt. 6.5) then
      layer2_index=layer2_index +1
      layer2(layer2_index) = icount
    else if (proximity(1) .ge. 6.5 .and. proximity(1) .lt. 9.5) then
      layer3_index = layer3_index +1
      layer3(layer3_index) = icount
    else if (proximity(2) .ge. 6.5 .and. proximity(2) .lt. 9.5) then
      layer3_index=layer3_index +1
      layer3(layer3_index) = icount
    end if
end do

print*, "number of molecules in each layer 1-3"
print*, layer1_index, layer2_index, layer3_index
 
close(50)

!.....read in most file data!
!.....read donor and acceptor energy of each molecule!

open(51, FILE="molecules.EDA")

do icount=1, 384
        read(51,*) tmp, donor_acceptor_energy(1:2,icount)
end do
close(51)
!.....read lowest donors/5lowest.......

open(52, FILE="molecules.lowest.donor")

do icount=1, 384
        read(52,*) tmp, lowest_donor(1:5,icount)
end do
close(52)
!.....read lowest acceptor/5lowest.....

open(53, FILE="molecules.lowest.acceptor")
do icount=1, 384
        read(53,*) tmp, lowest_acceptor(1:5,icount)
end do
close(53)

!.....read donor labels..........
open(54, FILE="molecules.lowest.donor.label")
do icount=1, 384
        read(54,*) tmp, donor_label(1:5,icount)
end do
close(54)

!....read acceptor labels........
open(55, FILE="molecules.lowest.acceptor.label")
do icount=1, 384
        read(55,*) tmp, acceptor_label(1:5,icount)
end do
close(55)
!print*, layer1_index
!...energy calculation for layer 1
do icount=1, layer1_index
L1donor(1:5,icount) = lowest_donor(1:5,layer1(icount))
L1acceptor(1:5,icount)=lowest_acceptor(1:5,layer1(icount))
!print*, L1donor(1:5,icount), L1acceptor(1:5,icount)
end do

do icount=1, layer2_index
L2donor(1:5,icount) = lowest_donor(1:5,layer2(icount))
L2acceptor(1:5,icount)=lowest_acceptor(1:5,layer2(icount))
!print*, L2donor(1:5,icount), L2acceptor(1:5,icount)
end do

do icount=1, layer3_index
L3donor(1:5,icount) = lowest_donor(1:5,layer3(icount))
L3acceptor(1:5,icount)=lowest_acceptor(1:5,layer3(icount))
!print*, L3donor(1:5,icount), L3acceptor(1:5,icount)
end do















!..average stregth of the 2 strongest bonds for both donor and acceptor

L1strength = 0
do icount=1, layer1_index
L1donorsum = L1donorsum + L1donor(1,icount) + L1donor(2,icount)
L1acceptorsum = L1acceptorsum + L1acceptor(1,icount) + L1acceptor(2,icount)
end do
L1donorsum=L1donorsum/(layer1_index*2)
L1acceptorsum=L1acceptorsum/(layer1_index*2)
L1strength=(L1acceptorsum+L1donorsum)*0.5

L2strength = 0
do icount=1, layer2_index
L2donorsum = L2donorsum + L2donor(1,icount) + L2donor(2,icount)
L2acceptorsum = L2acceptorsum + L2acceptor(1,icount) + L2acceptor(2,icount)
end do
L2donorsum=L2donorsum/(layer2_index*2)
L2acceptorsum=L2acceptorsum/(layer2_index*2)
L2strength=(L2acceptorsum+L2donorsum)*0.5





L3strength = 0
do icount=1, layer3_index
L3donorsum = L3donorsum + L3donor(1,icount) + L3donor(2,icount)
L3acceptorsum = L3acceptorsum + L3acceptor(1,icount) + L3acceptor(2,icount)
end do
L3donorsum=L3donorsum/(layer3_index*2)
L3acceptorsum=L3acceptorsum/(layer3_index*2)
L3strength=(L3acceptorsum+L3donorsum)*0.5

L1asymmetry_donor=0
L1asymmetry_donor_total=0
do icount=1, layer1_index
L1asymmetry_donor=L1asymmetry_donor+(1-(L1donor(2,icount)/L1donor(1,icount)))
print*, (L1donor(1,icount)/L1donor(2,icount))
end do
L1asymmetry_donor_total=L1asymmetry_donor/layer1_index
print*, L1asymmetry_donor_total

L2asymmetry_donor_total=0
do icount=1, layer2_index
L2asymmetry_donor=L2asymmetry_donor+(1-(L2donor(2,icount)/L2donor(1,icount)))
end do
L2asymmetry_donor_total=L2asymmetry_donor/layer2_index
!print*,L2asymmetry_donor_total

L3asymmetry_donor_total=0
do icount=1, layer3_index
L3asymmetry_donor=L3asymmetry_donor+(1-(L3donor(2,icount)/L3donor(1,icount)))
end do
L3asymmetry_donor_total=L3asymmetry_donor/layer3_index
!print*,L3asymmetry_donor_total

L1asymmetry_acceptor_total=0
do icount=1, layer1_index
L1asymmetry_acceptor=L1asymmetry_acceptor+(1-(L1acceptor(2,icount)/L1acceptor(1,icount)))
!print*, (L1acceptor(1,icount)/L1acceptor(2,icount))
end do
!if(layer1_index .gt. 0) then
L1asymmetry_acceptor_total=L1asymmetry_acceptor/layer1_index
!else
!L1asymmetry_acceptor_total=L1asymmetry_acceptor
!print*, L1asymmetry_acceptor_total

L2asymmetry_acceptor_total=0
do icount=1, layer2_index
L2asymmetry_acceptor=L2asymmetry_acceptor+(1-(L2acceptor(2,icount)/L2acceptor(1,icount)))
end do
L2asymmetry_acceptor_total=L2asymmetry_acceptor/layer2_index
!print*,L2asymmetry_acceptor_total

L3asymmetry_acceptor_total=0
do icount=1, layer3_index
L3asymmetry_acceptor=L3asymmetry_acceptor+(1-(L3acceptor(2,icount)/L3acceptor(1,icount)))
end do
L3asymmetry_acceptor_total=L3asymmetry_acceptor/layer3_index
!print*,L3asymmetry_acceptor_total
!print*, L1L1sumdonor
!!!!!printing

!print*, "donor energies"
!print*, "L1L1", L1L1sumdonor
!print*, "L1L2", L1L2averagedonor
!print*, "acceptor energies"
!print*, "L1L1", L1L1sumacceptor
!print*, "L1L2", L1L2averageacceptor
open(120, FILE="average_energies")
write(120,*) L1strength, L1donorsum, L1acceptorsum
write(120,*) L2strength, L2donorsum, L2acceptorsum
write(120,*) L3strength, L3donorsum, L3acceptorsum 
close(120)

open(123, FILE="asymmetry")
write(123,*) L1asymmetry_donor_total, L2asymmetry_donor_total, L3asymmetry_donor_total
write(123,*) L1asymmetry_acceptor_total, L2asymmetry_acceptor_total, L3asymmetry_acceptor_total
close(123)














!...........average donor/acceptor energy
!do icount=1, layer1_index
!        sumdonor_acceptor1(1:2) = sumdonor_acceptor1(1:2) + donor_acceptor_energy(1:2,(layer1(icount)))
!end do#do icount=1, layer2_index
!        sumdonor_acceptor2(1:2) = sumdonor_acceptor2(1:2) + donor_acceptor_energy(1:2,(layer2(icount)))
!end do
!do icount=1, layer3_index
!        sumdonor_acceptor3(1:2) = sumdonor_acceptor3(1:2) + donor_acceptor_energy(1:2,(layer3(icount)))
!end do
!
!open(52, FILE="molecules.lowest.acceptor")
!open(60, FILE="lowest.acceptor.L1")
!------file with acceptors in layer1
!do icount=1, 384
!        read*,  lowest_acceptor(1:5,icount) 
!        lowest_acceptor1(1:5)=
!end do
!close(52)
!close(60)        











!............average
!sumaverage1(1:2) = sumdonor_acceptor1(1:2)/layer1_index
!sumaverage2(1:2) = sumdonor_acceptor2(1:2)/layer2_index
!sumaverage3(1:2) = sumdonor_acceptor3(1:2)/layer3_index

!print*, "l1:" ,sumaverage1(1:2)
!print*, "l2:" ,sumaverage2(1:2)
!print*, "l3:" ,sumaverage3(1:2)
!----writefile

!open(52, FILE="average_energies")
!write(52,*) "Layer1", sumaverage1(1:2)
!write(52,*) "Layer2", sumaverage2(1:2)
!write(52,*) "Layer3", sumaverage3(1:2)
!close(52)

end program 
 
 


