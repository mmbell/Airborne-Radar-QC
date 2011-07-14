c
c =====================================
      subroutine chol_inv(a,x,b,n,ierr)
c =====================================

c   solve a.x = b using Choleski decomposition
c
c      real*8 a(n,n),x(n),b(n),p(100)     ! n <= 100
      real*4 a(n,n),x(n),b(n),p(100)     ! n <= 100
c
c  Choleski Decomposition  L. LT = A
c            L is returned in the lower triangle of A, except diagonal
c            elements   stored in p.
c
      ierr=0
      do i=1,n
         do j=1,n
            sum=a(i,j)
            do k=i-1,1,-1
               sum = sum - a(i,k)*a(j,k)
c!!!!	       print *,'SUM aik ajk',sum,a(i,k),a(j,k),i,j,k
            enddo
            if(i.eq.j) then
              if(sum.le.0.0) then
                ierr=-1
                print *,'  Choleski Decomposition Failed'
                go to 1
              else
                p(i)=sqrt(sum)
              endif
            else
              a(j,i) = sum/p(i)
            endif
         enddo
      enddo

   1  if(ierr.eq.0) then     !solve the linear system using forward
                             !and back substitution
c    Solve A x = L (LT x) = b
c                             First solve L y = b
c                             Then  solve LT x = y
c
      do i=1,n         ! solve L y = b, storing y in x
         sum=b(i)
         do k=i-1,1,-1
            sum = sum - a(i,k)*x(k)
         enddo
         x(i) = sum/p(i)
      enddo
c
      do i=n,1,-1      ! solve LT x = y
         sum=x(i)
         do k=i+1,n
            sum = sum - a(k,i)*x(k)
         enddo
         x(i) = sum/p(i)
      enddo
c
      endif

      end
c
