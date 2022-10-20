program one
implicit none
integer::i,j,k,itot
integer::id,etype
real::dummy
double precision,allocatable::x(:),nh(:)
double precision,allocatable::tgas(:),tdust(:),abun(:,:)
double precision,parameter::pc2cm=3.0856776d+18
double precision,allocatable::N(:)
character(len=100)::filein,filepdr,fileparams
character(len=100)::directory,prefix,filepop
real::Z,step
integer::Nspec,avtot
real::Avobs
real,allocatable::avin(:),pdfin(:)
real::vturb,c,kb,hp,pi,mhp,freq0(1:17),g(1:17,1:2),A(1:17),mh(1:17)
character(len=50)::spec(1:17)
double precision, allocatable::Tex(:,:),Bnu(:,:),pop(:,:,:),tau(:,:)
double precision::phi(1:17),sigma(1:17),tau_incr(1:17),frac(1:17)
double precision::t_r(1:17),t_a_i(1:17),t_a_ip1(1:17),dtau(1:17),Ntr(1:17),Ntot,Nrho
real::Ncol
real:: fuv, cosmicrays,Ntgas, Ntot_nopdf
real::av_bar,lav,s,m, lowerbound
real::tau_cii,tau_ci,tau_co
real::metallicity
!for making prefix
integer::ipref
character(len=2),allocatable::num(:)
character(len=10)::pref(1:1681)


!gfortran -o pdfchem_algorithm pdfchem_algorithm.f90


call constants

!Input parameters-------------------------------
close(1);open(unit=1,file='pdfchem.params',status='old')
read(1,*) av_bar, s, metallicity

if (metallicity == 0.1) then 
   directory = 'Z0p1'
elseif (metallicity == 0.5) then
   directory = 'Z0p5'
elseif (metallicity == 1.0) then 
   directory = 'Z1p0'
elseif (metallicity == 2.0) then 
   directory = 'Z2p0'
endif
directory = trim(adjustl(directory))//"/"
!------------------------------------------------

call makepdf

call makeprefix

open(unit=10,file='output.dat',status='replace')
do ipref = 1,1681 !No. PDR simulations
  prefix = pref(ipref)

  call readfile

  !calculate column densities of species
  N=0;Ntgas=0;Ntot_nopdf=0;Nrho=0
  allocate(Tex(1:17,1:itot),Bnu(1:17,1:itot))
  tau_incr=0
  do i=1,itot-1
   Avobs = 0.06*(nh(i)**0.69) 
   do k=1,avtot
     if (avin(k).gt.Avobs) exit
   enddo
   step = abs(x(i)-x(i+1))*pc2cm
   N(0) = N(0) + 0.5*(nh(i)+nh(i+1))*step*pdfin(k) !total column density
   N(1:Nspec) = N(1:Nspec) + 0.5*(nh(i)*abun(:,i) + nh(i+1)*abun(:,i+1))*step*pdfin(k) !species
                        !for the sequence of species, see the very last comment in this program
   Ntgas = Ntgas + 0.5*(nh(i)*tgas(i) + nh(i+1)*tgas(i+1))*step*pdfin(k) !for <Tgas>
   Nrho = Nrho + 0.5*(nh(i)**2*abun(31,i) + nh(i+1)**2*abun(31,i+1))*step*pdfin(k)

   !Excitation temperatures and Background Radiation for radiative transfer
   do j=1,17
     if (pop(j,2,i).eq.0.or.abs(g(j,2)*pop(j,1,i)/pop(j,2,i)/g(j,1)-1).lt.1e-2) then !handle exploding values
       Tex(j,i)=0
       Bnu(j,i)=0
     else 
       Tex(j,i)=(hp*freq0(j)/kb)/log(g(j,2)*pop(j,1,i)/pop(j,2,i)/g(j,1)) !excitation temperature
       Bnu(j,i)=(2.*hp*freq0(j)**3/c**2)/(exp(hp*freq0(j)/kb/Tex(j,i))-1) !Black body emission
     endif
   enddo
   sigma=(freq0/c)*sqrt(kb*Tgas(i)/mh+vturb**2/2.) !sigma value used in phi
   phi=1./sigma/sqrt(2.*pi) !phi value used in optical depth
   frac=0.5*((pop(:,1,i)+pop(:,1,i+1))*g(:,2)/g(:,1)-(pop(:,2,i)+pop(:,2,i+1)))
   tau_incr=tau_incr+phi*(A*c**2/8./pi/freq0**2)*frac*step !optical depth calculation
   tau(:,i)=tau_incr !record value
  enddo
  
  !solve radiative tranfer equation
  Ntr=0;Ntot=0;t_r=0;Ncol=0
  tau_cii=0;tau_ci=0;tau_co=0
  do i=1,itot-2
    Avobs = 0.06*(nh(i)**0.69) !Av,obs -- nH relation
    do k=1,avtot
      if (avin(k).gt.Avobs) exit
    enddo
    dtau=tau(:,i+1)-tau(:,i)
    do j=1,17 !frequencies explored (see subroutine constants)
      if (dtau(j).gt.1e10) then
        t_r(j)=Bnu(j,i)
      else if (dtau(j).gt.1e-6) then
        t_a_i(j)=Bnu(j,i)*((1-exp(-dtau(j)))/dtau(j)-exp(-dtau(j)))
        t_a_ip1(j)=Bnu(j,i+1)*(1.-(1.-exp(-dtau(j)))/dtau(j))
        t_r(j)=t_r(j)*exp(-dtau(j))+t_a_i(j)+t_a_ip1(j)
      else
        t_r(j)=t_r(j)*(1-dtau(j))+(Bnu(j,i)+Bnu(j,i+1))*dtau(j)/2.
      endif
    enddo
    step = abs(x(i)-x(i+1))*pc2cm
    Ntot = Ntot + pdfin(k)
    Ntr = Ntr + (t_r*c**2/2./kb/freq0**2)*pdfin(k)
    Ncol = Ncol + (0.5*(nh(i)+nh(i+1))*step)
    tau_cii = tau_cii + dtau(1)*pdfin(k)
    tau_ci  = tau_ci  + dtau(2)*pdfin(k)
    tau_co  = tau_co  + dtau(8)*pdfin(k)
  enddo

  !write output.dat
  write(10,'(100ES11.2)') fuv, cosmicrays, Z, & !parameters of PDR simulation
          Ntgas/N(0), N(1:33)/N(0), & !average column densities of species
          Ntr(1), Ntr(2), Ntr(4), Ntr(8:17) !CII, CI(1-0), CI(2-1), CO(1-0...10-9)

  deallocate(tgas); deallocate(tdust); deallocate(abun)
  deallocate(N); deallocate(x); deallocate(nh); deallocate(Tex)
  deallocate(Bnu); deallocate(pop); deallocate(tau)

enddo
write(6,*) 'Finished!'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains


!Make Av-PDF -----------------------------------------------------
subroutine makepdf
m = log(av_bar) - s**2/2.
avtot=500
allocate(avin(1:avtot+1),pdfin(1:avtot+1))
open(unit=111,file='avpdf.dat')
do i=1,avtot+1
  lav = -2 + 4.*real(i)/real(avtot+1)
  avin(i) = 10**lav
  pdfin(i) = pdf(real(avin(i)),s,m)
  if (pdfin(i)<1e-10) pdfin(i)=1e-10
  write(111,*) avin(i),pdfin(i)
enddo
return
end subroutine

!Av,obs - PDF function
real function pdf(x,s,m)
real::x,s,m

pdf=(1./s/sqrt(2.*3.141592653589))*exp(-(log(x)-m)**2/2./s**2) 

return
end function
!-----------------------------------------------------------------


!Make prefix for all inputs --------------------------------------
subroutine makeprefix
allocate(num(0:40))
do i=0,40
    write (num(i),"(I0.2)") i
enddo

k=0
do i=0,40
  do j=0,40
     k=k+1
     pref(k) = "zcr"//trim(adjustl(num(i)))//"_uv"//trim(adjustl(num(j)))
  enddo
enddo
return
end subroutine
!-----------------------------------------------------------------


!Read params, abundances, density, level populations -------------
!If you are using your own PDR code, this is what you need to modify
subroutine readfile
fileparams=trim(adjustl(directory))//trim(adjustl(prefix))//".params" !parameters
open(unit=1,file=fileparams,status='old')
read(1,*) fuv, cosmicrays, Z
close(1)

filein=trim(adjustl(directory))//trim(adjustl(prefix))//".pdr.fin" !abundances
filepdr=trim(adjustl(directory))//trim(adjustl(prefix))
open(unit=1,file=filein,status='old')
itot=0
do 
 read(1,*,end=100) dummy
 itot=itot+1
enddo
100 continue
rewind(1)
allocate(tgas(1:itot),tdust(1:itot),abun(Nspec,1:itot))
allocate(N(0:Nspec),x(1:itot),nh(1:itot))
do i=1,itot
 read(1,*) id,x(i),dummy,tgas(i),tdust(i),etype,nh(i),dummy,abun(1:Nspec,i)
enddo
close(1)
open(unit=1,file='variable_density.dat',status='old') !variable density slab
do i=1,itot
 read(1,*) x(i)
enddo
close(1)
  
filepop=trim(adjustl(directory))//trim(adjustl(prefix))//".spop.fin" !level populations
open(unit=3,file=filepop,status='old')
allocate(pop(1:17,1:2,1:itot),tau(1:17,1:itot))
do i=1,itot
  read(3,*) dummy,dummy,pop(1,1,i),pop(1,2,i),dummy,dummy,dummy,&
       &pop(2,1,i),pop(2,2,i),pop(3,2,i),dummy,dummy,&
       &pop(5,1,i),pop(5,2,i),pop(6,2,i),dummy,dummy,&
       &pop(8,1,i),pop(8,2,i),pop(9,2,i),pop(10,2,i),pop(11,2,i),&
       &pop(12,2,i),pop(13,2,i),pop(14,2,i),pop(15,2,i),pop(16,2,i),pop(17,2,i)
  !record pairs
  pop(3,1,i)=pop(2,1,i);  pop(4,2,i)=pop(3,2,i);  pop(4,1,i)=pop(2,2,i)
  pop(6,1,i)=pop(5,1,i);  pop(7,2,i)=pop(6,2,i);  pop(7,1,i)=pop(5,2,i)
  pop(9,1,i)=pop(8,2,i);  pop(10,1,i)=pop(9,2,i);  pop(11,1,i)=pop(10,2,i)
  pop(12,1,i)=pop(11,2,i);  pop(13,1,i)=pop(12,2,i);  pop(14,1,i)=pop(13,2,i)
  pop(15,1,i)=pop(14,2,i);  pop(16,1,i)=pop(15,2,i);  pop(17,1,i)=pop(16,2,i)
enddo
close(3)

return
end subroutine
!-----------------------------------------------------------------


!Constants for calculations --------------------------------------
subroutine constants
c=2.9979246e10 !cm/s
kb=1.380650e-16 !erg / K  *or*  g cm^2 / K s^2
hp=6.6260696e-27 !erg s  *or*  g cm^2 / s
pi=3.1415927
mhp=1.6726218e-24 !g
Nspec=33

vturb = 1e5 !microturbulent velocity in cm/s
            !NOTE: all given PDR calculations used 1~km/s. 

freq0(1)=1900.5369e9     ;g(1,1)=2.0  ;g(1,2)=4.0   ; spec(1)="CII 158um"
freq0(2)=492.16065e9     ;g(2,1)=1.0  ;g(2,2)=3.0   ; spec(2)="CI (1-0)"
freq0(3)=1301.50262e9    ;g(3,1)=1.0  ;g(3,2)=5.0   ; spec(3)="CI (2-0)"
freq0(4)=809.34197e9     ;g(4,1)=3.0  ;g(4,2)=5.0   ; spec(4)="CI (2-1)"
freq0(5)=4744.77749e9    ;g(5,1)=5.0  ;g(5,2)=3.0   ; spec(5)="OI  1-0 "
freq0(6)=6804.84658e9    ;g(6,1)=5.0  ;g(6,2)=1.0   ; spec(6)="OI  2-0 "
freq0(7)=2060.06909e9    ;g(7,1)=3.0  ;g(7,2)=1.0   ; spec(7)="OI  2-1 "
freq0(8)=115.2712018e9   ;g(8,1)=1.0  ;g(8,2)=3.0   ; spec(8)="CO (1-0)"
freq0(9)=230.538e9       ;g(9,1)=3.0  ;g(9,2)=5.0   ; spec(9)="CO (2-1)"
freq0(10)=345.7959899e9  ;g(10,1)=5.0 ;g(10,2)=7.0  ; spec(10)="CO (3-2)"
freq0(11)=461.040768e9   ;g(11,1)=7.0 ;g(11,2)=9.0  ; spec(11)="CO (4-3)"
freq0(12)=576.2679305e9  ;g(12,1)=9.0 ;g(12,2)=11.0 ; spec(12)="CO (5-4)"
freq0(13)=691.4730763e9  ;g(13,1)=11.0;g(13,2)=13.0 ; spec(13)="CO (6-5)"
freq0(14)=806.6518060e9  ;g(14,1)=13.0;g(14,2)=15.0 ; spec(14)="CO (7-6)"
freq0(15)=921.7997000e9  ;g(15,1)=15.0;g(15,2)=17.0 ; spec(15)="CO (8-7)"
freq0(16)=1036.9123930e9 ;g(16,1)=17.0;g(16,2)=19.0 ; spec(16)="CO (9-8)"
freq0(17)=1151.9854520e9 ;g(17,1)=19.0;g(17,2)=21.0 ; spec(17)="CO (10-9)"

!Einstein A coefficients
A(1)=2.321e-06 ;  mh(1)=12.*mhp  
A(2)=7.880e-08 ;  mh(2)=12.*mhp  
A(3)=1.810e-14 ;  mh(3)=12.*mhp  
A(4)=2.650e-07 ;  mh(4)=12.*mhp  
A(5)=8.910e-05 ;  mh(5)=16.*mhp  
A(6)=1.340e-10 ;  mh(6)=16.*mhp  
A(7)=1.750e-05 ;  mh(7)=16.*mhp  
A(8)=7.203e-08 ;  mh(8)=28.*mhp  
A(9)=6.910e-07 ;  mh(9)=28.*mhp  
A(10)=2.497e-06;  mh(10)=28.*mhp 
A(11)=6.126e-06;  mh(11)=28.*mhp 
A(12)=1.221e-05;  mh(12)=28.*mhp 
A(13)=2.137e-05;  mh(13)=28.*mhp 
A(14)=3.422e-05;  mh(14)=28.*mhp 
A(15)=5.134e-05;  mh(15)=28.*mhp 
A(16)=7.330e-05;  mh(16)=28.*mhp 
A(17)=1.006e-04;  mh(17)=28.*mhp 

return
end subroutine
!-----------------------------------------------------------------


!Sequence of species in the N(1:33) array
! 1,H3+;  2,He+;   3,Mg;     4,H2+;   5,CH5+
! 6,O2;   7,CH4+;  8,O+;     9,OH+;  10,Mg+
!11,C+;  12,CH4;  13,H2O+;  14,H3O+; 15,CO+
!16,O2+; 17,CH2;  18,H2O;   19,H+;   20,CH3+
!21,CH;  22,CH3;  23,HCO+;  24,CH2+; 25,C
!26,He;  27,CH+;  28,CO;    29,OH;   30,O
!31,H2;  32,H;    33,e-; 

end program

