program heat
implicit none
include 'mpif.h'
real*8:: dx,dt,a,b,r,final_t,right,left,sright,sleft,x
integer:: i,k,nt,nx,J,start_node
real*8, allocatable:: u(:,:),sol(:),solt(:,:)
real*8, parameter:: pi=2*asin(1.d0)
integer:: comm,rank,nproc,ierr,request
real*8:: start,finish,truesol,error

comm = MPI_COMM_WORLD
call MPI_INIT(ierr)
call MPI_COMM_RANK(comm, rank, ierr)
call MPI_COMM_SIZE(comm, nproc, ierr)

call cpu_time(start)
nx=5000                    !number of spatial intervals (points from 0 to nx)
nx=nx-mod(nx-1,nproc)      !make nx satisfy: (nx-1)/nproc is integer
final_t=0.1d0              !final integration time     

dx=1.d0/nx                 !spatial step
dt=.5d0*dx**2              !temporal step
nt=int(final_t/dt)         !number of temporal intervals

r=dt/(dx**2)               !parameter for Euler  
a=0.d0                     !left boundary condition
b=0.d0                     !right boudnary condition

open(10,file='heat.dat')
J=(nx-1)/nproc+2
allocate(u(1:J,0:1))
u=0
allocate(solt(1:J,1:nproc),sol(0:nx))

!% % % % % % % % % % % % % % % % % % % % % % 

start_node=(J-2)*rank

!set initial condition
do i=1,J
  x=(start_node+i-1)*dx
  u(i,0)=sin(2*pi*x)+2*sin(5*pi*x)+3*sin(20*pi*x)
enddo


do i=1,nt

!perform Euler method for points in the interior
  do k=2,J-1
      u(k,1)=r*(u(k-1,0)+u(k+1,0))+(1-2*r)*u(k,0)
  enddo

!message swap between adjacent processors
  if (rank==0) then
    
    u(1,1)=a

    sright=u(J-1,1)
    call mpi_send(sright,1,mpi_real8,rank+1,1,comm,ierr)

    call mpi_recv(right,1,mpi_real8,rank+1,1,comm,mpi_status_ignore,ierr)
    u(J,1)=right

  elseif (rank==nproc-1) then

    call mpi_recv(left,1,mpi_real8,rank-1,1,comm,mpi_status_ignore,ierr)
    u(1,1)=left
    
    sleft=u(2,1)
    call mpi_send(sleft,1,mpi_real8,rank-1,1,comm,ierr)
  
    u(J,1)=b

  else 

    call mpi_recv(left,1,mpi_real8,rank-1,1,comm,mpi_status_ignore,ierr)    
    u(1,1)=left    

    sleft=u(2,1)
    call mpi_send(sleft,1,mpi_real8,rank-1,1,comm,ierr)

    sright=u(J-1,1)
    call mpi_send(sright,1,mpi_real8,rank+1,1,comm,ierr)
    
    call mpi_recv(right,1,mpi_real8,rank+1,1,comm,mpi_status_ignore,ierr)    
    u(J,1)=right

  endif

!make current solution the initial one for next iteration
u(:,0)=u(:,1)
enddo

!gather pieces of the solution
call mpi_gather(u(:,1),J,mpi_real8,solt,J,mpi_real8,0,comm,ierr)
if (rank==0) then

!put pieces in order
  do i=1,nproc
    sol((i-1)*(j-2):i*(j-2)-1)=solt(1:j-2,i)
  enddo
  sol(nx-1:nx)=solt(j-1:j,nproc)
  write(10,*) sol
  call cpu_time(finish)
  print *,'Subdivisions         :',nx
  print *,'Processors           :',nproc
  print *,'Points per processor :',J
  print *,'Computational time   :',finish-start

!compute error 2-norm
  error=0.d0
  do i=0,nx
    x=i*dx
    truesol=exp(-4*pi**2*final_t)*sin(2*pi*x)+&
            2*exp(-25*pi**2*final_t)*sin(5*pi*x)+& 
            3*exp(-400*pi**2*final_t)*sin(20*pi*x)
    error=error+(truesol-sol(i))**2
  enddo
  print *,'Error 2-norm         :',sqrt(error)

endif


deallocate(u,solt,sol)
call mpi_finalize(ierr)
end program heat
