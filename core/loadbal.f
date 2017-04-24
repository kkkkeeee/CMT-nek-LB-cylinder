      subroutine assign_partitions !(gllnid, lelt, nelgt, np)
c     This subroutine is used for update gllnid            
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL' !these include contains gllnid, lelt, and np     

      !integer gllnid(1) !, lelt, nelgt, np
      common /elementload/ gfirst, inoassignd, resetFindpts, pload(lelg)
      integer gfirst, inoassignd, resetFindpts
      integer pload
      integer psum(lelg)
      integer nw, i, part  !the number of particles
c     nw=200 
      
      call izero(pload, lelg)
c      call randGenet(pload, nelgt, nw) !random assign pload
c     call preSum(pload, psum, nelgt) !calculate the prefix sum of pload
c      call ldblce(psum, nelgt, gllnid, np) !equally distribute the load to np processor, assigned to gllnid
      
c     uniformally assign element
      part = nelgt/np
      do i=1, nelgt
         gllnid(i) = (i-1)/part
      enddo
c      call printi(gllnid, nelgt)

      end subroutine

c----------------------------------------------------------------

c This subroutine assign normalized random number to pload of length len; and make the sum of pload equal to total, pload store the actual number of particles in elements
      subroutine randGenet(pload, len, total)
      integer len, total
      integer pload(len)
c     local variable
      real temppload(len)
      integer i, seed, summ2
      real summ, summ1
c     real mu, sigma, r4_normal_ab

      summ=0
      summ1=0
      summ2=0
      !call RANDOM_NUMBER(pload)   !random generate len values and store in pload
      seed = 123456789
      mu = 20.0e+01
      sigma = 170.0e+00

      do i=1, len
c        r4_normal_ab introduced new file normal.f
c        temppload(i) = r4_normal_ab (mu, sigma, seed) !generate normalize number
          if (temppload(i) .lt. 0) then
             temppload(i) = 0
          endif
      enddo
      !print *, 'temppload:'
      !call printr(temppload, len)

      do 10 i=1, len
         summ=summ+temppload(i)
   10 continue

      !print *, 'The sum of pload is: ', summ 
      do 20 i=1, len
        temppload(i)=temppload(i)*total/summ  !make the total of pload equal to total
        summ1=summ1+temppload(i)
   20 continue
      !print *, 'The new sum of pload is: ', summ1

      !print *, 'temppload:'
      !call printr(temppload, len)     

      do i = 1, len-1
         pload(i) = int(temppload(i))! convert real to int, cut off
         summ2 = summ2 + pload(i)
      enddo
      if ( total .le. summ2) then
         pload(len) = 0
      else
         pload(len) = total - summ2
      endif

      return
      end

c------------------------------------------------------------------

c     subroutine to get prefix sum
      subroutine preSum(pload, psum, len)
      integer pload(len)
      integer psum(len)
      integer i
      psum(1)=pload(1)
      do 30 i=2, len
          psum(i)=psum(i-1)+pload(i)
  30  continue

c      do 50 i=1, len
c         pload(i)=psum(i)
c  50  continue

      return
      end

c--------------------------------------------------------
c     assign to corresponding processor
      subroutine ldblce(psum, len, gllnid, np)
         integer psum(len)
         integer np
         integer gllnid(len)

         integer i,j, k, flag
         integer pos(np+1)  !pos(i)-pos(i+1)-1 belong to processor i-1
         real diff, diffn, thresh
         i=1
         flag = 0
         thresh=psum(len)*1.0/np*i
         call izero(pos, np+1)
         pos(1)=1
         do 70 j=2, len
            diff=abs(thresh-psum(j-1))
            diffn=abs(thresh-psum(j))
            if(diff .ge. diffn) then
               pos(i+1)=j+1
               !write(*,*) "i/:", i, "pos(i):", pos(i)
            else
               pos(i+1)=j
               !write(*,*) "i/:", i, "pos(i):", pos(i)
               i=i+1
c              print *, 'thresh', thresh
               thresh=psum(len)*1.0/np*i
            endif
  70      continue
c         print *, 'prefix sum, len: ', len
c         call printi(psum, len)
c         print *, ' i', i
          !call printi(pos, np+1)
          if( i .lt. np) then ! this part is for the partition less than np
              do k = i+2, np+1
                 pos(k) = len + 1
              enddo
          endif 
c         print *, 'printing pos'
c         call printi(pos, np+1)

          do 80 i=1, np
             do 90 j=pos(i), pos(i+1)-1
                gllnid(j)=i-1
  90         continue
  80      continue
c         print *, 'printing gllnid, length: ',len 
c         call printi(gllnid, len)
          
      return
      end

c-----------------------------------------------------------------

c      print array real
       subroutine printr(pload, len)
          real pload(len)
          integer i
          do 40 i=1, len
             print *, pload(i)
   40     continue
       return
       end

c      print array integer
       subroutine printi(pos, len)
          integer pos(len)
          integer i
          do 40 i=1, len
             print *, pos(i)
   40     continue
       return
       end


c------------------------------------------------------------------
c     recompute partitions
       subroutine recompute_partitions
          include 'SIZE'
          include 'INPUT'
          include 'PARALLEL'          
          include 'TSTEP'          
          include 'SOLN'          
 
          !parameter (lr=16*ldim,li=5+6)
          parameter (lr=76+4,li=10)
          common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
          common /elementload/ gfirst, inoassignd, 
     >                 resetFindpts, pload(lelg)
          integer gfirst, inoassignd, resetFindpts, pload

          integer nw
          common /particlenumber/ nw

          integer newgllnid(lelg), trans(3, lelg), trans_n, psum(lelg)
c          integer total
          integer e,eg,eg0,eg1,mdw,ndw
          integer ntuple, i, el, delta, nxyz
      common /ptpointers/ jrc,jpt,je0,jps,jpid1,jpid2,jpid3,jpnn,jpid
     >                   ,jai,nai,    jr,jd,jx,jy,jz,jv0,ju0,jf0,jfusr
     >                   ,jfqs,jfun,jfiu,jtaup,jcd,jdrhodt,jre,jDuDt
     >                   ,jtemp,jrho,jrhop,ja,jvol,jdp,jar,jx1,jx2,jx3
     >                   ,jv1,jv2,jv3,ju1,ju2,ju3,nar,jvol1,jgam
     >                   ,jaa, jab, jac, jad

          common  /cparti/ ipart(li,llpart)
          real   xerange(2,3,lelt)
          common /elementrange/ xerange

          common  /iparti/ n,nr,ni
          integer particleMap(3, lelt)
          integer gaselement(np)
          real ratio
          
c         if(nid .eq. 0) then
c             print *, 'old pload:'
c             call printi(pload, nelgt)
c         endif
          nxyz = nx1*ny1*nz1   !total # of grid point per element
          delta = ceiling(nw*1.0/nelgt)
          ratio = 1.0
          ntuple = nelt
          do i=1,ntuple
             eg = lglel(i)
             particleMap(1,i) = eg
             particleMap(2,i) = 0        ! processor id to send for element eg
             particleMap(3, i) = 0      !  #of particles in each element, reinitialize to 0, otherwise it keeps the previous value
          enddo

          do ip=1,n
             el = ipart(je0, ip) + 1      ! element id in ipart is start from 0, so add 1
             particleMap(3, el) = particleMap(3, el) + 1
          enddo

c         gas_right_boundary = exp(TIME/2.0)



          do i=1,ntuple
c            x_left_boundary = xerange(1,1,i)
c            if (x_left_boundary .lt. gas_right_boundary) then
c            !if (vx(1,1,1,i) .ne. 0) then
                particleMap(3, i) = particleMap(3, i) + delta*ratio
c            else
c               particleMap(3, i) = 0
c               !print *, 'element gas 0'
c            endif
          enddo
          mdw=3
          ndw=nelgt
          key = 2  ! processor id is in wk(2,:)
          call crystal_ituple_transfer(cr_h,particleMap,
     $                                 mdw,ntuple,ndw,key)
          
c          total=lelt*10
          
          if (nid .eq. 0) then
             key=1
             nkey = 1
             call crystal_ituple_sort(cr_h,particleMap,mdw,
     $                                ntuple,key,nkey)
             do i=1,ntuple
                pload(i) = particleMap(3, i)
             enddo

c            print *, 'new pload:'
c            call printi(pload, nelgt)
             !print *, 'new pload/n'
             !call printr(newPload, lelt)
             call izero(psum, lelg)
             call preSum(pload, psum, nelgt)
             print *, 'recompute_partitions: psum(nelgt): ', psum(nelgt)
             print *, 'ratio:', ratio
             !call printr(newPload, lelt)
   
             call izero(newgllnid, lelg) 
             call ldblce(psum, nelgt, newgllnid, np_)

c            print *, 'print new gllnid'
c            call printi(newgllnid, nelgt)
c            call izero(gaselement, np)
c            do i=1, nelgt
c               if (pload(i) .ne. 0) then
c               gaselement(newgllnid(i)+1)=gaselement(newgllnid(i)+1)+1
c               endif
c            enddo 
c            do i=1, np
c                print *, '# gas element on', i-1, 'is: ', gaselement(i)
c            enddo
          endif
          call bcast(newgllnid,4*nelgt)

          call izero(trans, 3*lelg)
          call track_elements(gllnid, newgllnid, nelgt, trans, 
     $                               trans_n, lglel)
c         print *, 'print trans'
c         do 110 i=1, trans_n
c           print *, trans(1, i), trans(2, i), trans(3, i)
c 110  continue

          call track_particles(trans, trans_n)
          call mergePhigArray(newgllnid, trans, trans_n)
          call mergeUArray(newgllnid)

          call copy(gllnid, newgllnid, nelgt)

          return
          end
          
c------------------------------------------------------
c subroutine of track elements to be send and received
       subroutine track_elements(gllnid, newgllnid, len, trans, 
     $                     trans_n, lglel)
       include 'SIZE'
       integer gllnid(len), newgllnid(len), trans(3, len), trans_n
       integer lglel(1)
       !trans: first column stores the source pid, second column stores the target pid, the third column stores the element id
       integer i, j; !local variable

       trans_n=1;
       j=0;
       do i=1, len
          if ((gllnid(i) .eq. nid)) then
             j = j+1;
             if(gllnid(i) .ne. newgllnid(i)) then
c since added gllnid(i) .eq. nid, right now, the trans in processor i only store the elements that he shold send. Not all the processors.            
             trans(1, trans_n) = gllnid(i);
             trans(2, trans_n) = newgllnid(i);
             trans(3, trans_n) = lglel(j)
             trans_n=trans_n+1;
           endif
         endif
       enddo
       trans_n=trans_n-1;  !the length of trans
c      print *, trans_n, 'print again in track elements'
c       do 110 i=1, trans_n
          !do 120 j=1, width
c           print *, i, trans(i,1), trans(i,2), trans(i,3)
c  120     continue
c  110  continue

       return
       end
 
c-----------------------------------------------------
c update the particles that are in the elements to be
c transferred, and set jpt to the destination processor
       subroutine track_particles(trans, trans_n)
       include 'SIZE'
       include 'PARALLEL' 

       integer trans(3, lelg), trans_n
       parameter (lr=76+4,li=10)
       common  /iparti/ n,nr,ni
       common  /cpartr/ rpart(lr,llpart) ! Minimal value of lr = 16*ndim
       common  /cparti/ ipart(li,llpart) ! Minimal value of li = 5
c      common /ptpointers/ jrc,jpt,je0,jps,jpid1,jpid2,jpid3,jpnn,jai
c    >                ,nai,jr,jd,jx,jy,jz,jx1,jx2,jx3,jv0,jv1,jv2,jv3
c    >                ,ju0,ju1,ju2,ju3,jf0,jar,jaa,jab,jac,jad,nar,jpid
      common /ptpointers/ jrc,jpt,je0,jps,jpid1,jpid2,jpid3,jpnn,jpid
     >                   ,jai,nai,    jr,jd,jx,jy,jz,jv0,ju0,jf0,jfusr
     >                   ,jfqs,jfun,jfiu,jtaup,jcd,jdrhodt,jre,jDuDt
     >                   ,jtemp,jrho,jrhop,ja,jvol,jdp,jar,jx1,jx2,jx3
     >                   ,jv1,jv2,jv3,ju1,ju2,ju3,nar,jvol1,jgam
     >                   ,jaa, jab, jac, jad

       common /myparth/ i_fp_hndl, i_cr_hndl
     
       
       integer ip, it, e
       logical partl         ! This is a dummy placeholder, used in cr()
       nl = 0                ! No logicals exchanged

       
c     change ipart(je0,i) to global element id
       ip=0
       do ip = 1, n
           e = ipart(je0, ip) + 1 ! je0 start from 0
           ipart(je0, ip) = lglel(e)
       enddo

       do ip = 1, n
          !e = ipart(je0,ip)
          do it = 1, trans_n
             if(ipart(je0, ip) .eq. trans(3, it)) then  
                ipart(jpt, ip) = trans(2, it) !new processor
                exit
             endif
          enddo
          ipart(jps, ip) = ipart(jpt, ip)
       enddo 

       call crystal_tuple_transfer(i_cr_hndl,n,llpart
     $              , ipart,ni,partl,nl,rpart,nr,jps)


       return 
       end 
c------------------------------------------------------
c subroutine to merge u array
          subroutine mergeUArray(newgllnid)
          include 'SIZE'
          include 'INPUT'
          include 'PARALLEL'          
          include 'CMTDATA'          

          integer newgllnid(lelg), trans(3, lelg)
          real uarray(lx1*ly1*lz1, lelt) 
          integer procarray(3, lelg) !keke changed real to integer to use crystal_ituple_transfer
          real tempu(lx1*ly1*lz1, lelt)
          logical partl
          integer nl, sid, eid
          integer key, nkey, ifirstelement, ilastelement, u_n, trans_n

          index_n=1
          trans_n=0
          nxyz = lx1*ly1*lz1*toteq
          nl=0
          do i=1, nelt
              ieg = lglel(i)
              if ((gllnid(ieg) .eq. nid) .and. 
     $                     (gllnid(ieg) .ne. newgllnid(ieg))) then
                  procarray(1, index_n) = gllnid(ieg)
                  procarray(2, index_n) = newgllnid(ieg)
                  procarray(3, index_n) = ieg
                  trans(1, index_n) = gllnid(ieg)
                  trans(2, index_n) = newgllnid(ieg)
                  trans(3, index_n) = ieg
                  do l=1, toteq
                      do k=1, lz1
                          do n=1, ly1
                             do m=1, lx1
                                ind=m+(n-1)*ly1+(k-1)*lz1*ly1+
     $                                           (l-1)*toteq*lz1*ly1
                                uarray(ind, index_n) = u(m,n,k,l,i)
                             enddo
                          enddo
                      enddo
                  enddo
              index_n=index_n+1
              endif
          enddo
          index_n=index_n-1
          trans_n=index_n
c          print *, index_n, trans_n, lelt, nid
           
          key=2
          call crystal_tuple_transfer(cr_h, trans_n, nelgt, trans,
     $                    3, partl, nl, uarray,nxyz,key)
c         print *, 'nid: ', nid, 'trans_n', trans_n      
          key=3
          nkey=1
          call crystal_tuple_sort(cr_h, trans_n, trans, 3, 
     $               partl, nl, uarray, nxyz, key,nkey)

          !Update u
          ! set sid and eid
          sid=1
          ifirstelement = lglel(sid)
          u_n=0
          do i=1, index_n
              !if (procarray(3,i) .lt. ifirstelement) then
      !            set tempu from uarray
               !    u_n=u_n+1
                !   do j=1, nxyz
                 !     tempu(j,u_n)=uarray(j,i)
                  ! enddo
              if (procarray(3,i) .eq. ifirstelement) then
                       sid=sid+ 1
                       ifirstelement=lglel(sid)
              endif
          enddo
          eid = nelt
          ilastelement = lglel(eid)
          do i=index_n, 1, -1
               if (procarray(3,i) .eq. ilastelement) then
                   eid = eid -1
                   ilastelement = lglel(eid)
               endif
          enddo

          ifirstelement = lglel(sid)
          do i=1, trans_n
              if (trans(3, i) .lt. ifirstelement) then
                  ! set tempu from uarray
                   u_n=u_n+1
                   do j=1, nxyz
                      tempu(j,u_n)=uarray(j,i)
                   enddo
              endif
          enddo

          
          if (sid .le. eid) then
              do i=sid, eid      !set tempu from original u
                  u_n=u_n+1
                  do l=1, toteq
                      do k=1, lz1
                          do n=1, ly1
                             do m=1, lx1
                                ind=m+(n-1)*ly1+(k-1)*lz1*ly1+
     $                                           (l-1)*toteq*lz1*ly1
                                tempu(ind, u_n) = u(m,n,k,l,i)
                             enddo
                          enddo
                      enddo
                  enddo
              enddo 
          endif
          ilastelement = lglel(eid)
          do i=1, trans_n
              if (trans(3, i) .gt. ilastelement) then
                  ! set tempu from uarray
                   u_n=u_n+1
                   do j=1, nxyz
                      tempu(j,u_n)=uarray(j,i)
                   enddo
              endif
          enddo
          !copy tempu to u
              do i=1, u_n
                  do l=1, toteq
                      do k=1, lz1
                          do n=1, ly1
                             do m=1, lx1
                                ind=m+(n-1)*ly1+(k-1)*lz1*ly1+
     $                                           (l-1)*toteq*lz1*ly1
                                uarray(ind, index_n) = u(m,n,k,l,i)
                                u(m,n,k,l,i) = tempu(ind, i)
                             enddo
                          enddo
                      enddo
                  enddo
              enddo
c      print *, "Update u array: u_n", u_n
       end
c----------------------------------------------------------------------
c subroutine to merge phig array
          subroutine mergePhigArray(newgllnid, trans, trans_n)
          include 'SIZE'
          include 'INPUT'
          include 'PARALLEL'          
          include 'TSTEP'          
          include 'CMTDATA'          

          integer newgllnid(lelg), trans(3, lelg), trans_n
          real phigarray(lx1*ly1*lz1, lelt) !keke changed real to integer to use crystal_ituple_transfer
          integer procarray(3, lelg) !keke changed real to integer to use crystal_ituple_transfer
          real tempphig(lx1*ly1*lz1, lelt)
          logical partl
          integer nl, sid, eid
          integer key, nkey, ifirstelement, ilastelement, phig_n

c         Step 1: Build phigarray for the elements to be transferred to neighboring processes.
          index_n=1;
          nl=0
          nxyz = lx1*ly1*lz1
c         for verification !!!!! 
          call copy(phig, vx, nxyz*nelt)
c         for verification !!!!! 
          do i=1, nelt
              ieg = lglel(i)
              if ((gllnid(ieg) .eq. nid) .and. 
     $                     (gllnid(ieg) .ne. newgllnid(ieg))) then
                  procarray(1, index_n) = gllnid(ieg)
                  procarray(2, index_n) = newgllnid(ieg)
                  procarray(3, index_n) = ieg
                  do k=1, lz1
                      do n=1, ly1
                         do m=1, lx1
                            ind=m+(n-1)*ly1+(k-1)*lz1*ly1 
                            phigarray(ind, index_n) = phig(m,n,k,i)
                         enddo
                      enddo
                  enddo
              index_n=index_n+1
              endif
          enddo
          index_n=index_n-1
         
c          print *, index_n, lelg, nid
           
c         Step 2: Send/receive phigarray to/from neighbors
          key=2
          call crystal_tuple_transfer(cr_h, trans_n, nelgt, trans,
     $                    3, partl, nl, phigarray,nxyz,key)
c         print *, 'nid: ', nid, 'received trans_n', trans_n      
c         Step 2: Sort the received phigarray based on global element number
          key=3
          nkey=1
          call crystal_tuple_sort(cr_h, trans_n, trans, 3, 
     $               partl, nl, phigarray, nxyz, key,nkey)

c         Step 4: set start id and end id of phig of existing (not transferred) elements
          !start id is the index of the first element that has not been sent to the left neighbor
          !end id is the index of the last element that has not been sent to the right neighbor
          sid=1
          ifirstelement = lglel(sid) 
          phig_n=0
          do i=1, index_n
              !if (procarray(3,i) .lt. ifirstelement) then
      !            set tempphig from phigarray
                   !phig_n=phig_n+1
                   !do j=1, nxyz
                    !  tempphig(j,phig_n)=phigarray(j,i)
                   !enddo
              if (procarray(3,i) .eq. ifirstelement) then
                       sid=sid+ 1
                       ifirstelement=lglel(sid)
              endif
          enddo
          eid = nelt
          ilastelement = lglel(eid)
          do i=index_n, 1, -1
               if (procarray(3,i) .eq. ilastelement) then
                   eid = eid -1
                   ilastelement = lglel(eid)
               endif
          enddo

c         Step 5: Update local phig based on elements 
c                 a) received from left neighbor
c                 b) existing elements
c                 c) received from right neighbor
          ifirstelement = lglel(sid) 
          do i=1, trans_n
              if (trans(3, i) .lt. ifirstelement) then
                  ! set tempphig from phigarray
                   phig_n=phig_n+1
                   do j=1, nxyz
                      tempphig(j,phig_n)=phigarray(j,i)
                   enddo
              endif
          enddo
          
          if (sid .le. eid) then  !set tempphig from original phig
              do i=sid, eid      
                  phig_n=phig_n+1
                  do k=1, lz1
                      do n=1, ly1
                         do m=1, lx1
                            ind=m+(n-1)*ly1+(k-1)*lz1*ly1 
                            tempphig(ind, phig_n) = phig(m,n,k,i)
                         enddo
                      enddo
                  enddo
              enddo 
          endif
          ilastelement = lglel(eid)
          do i=1, trans_n
              if (trans(3, i) .gt. ilastelement) then
                  ! set tempphig from phigarray
                   phig_n=phig_n+1
                   do j=1, nxyz
                      tempphig(j,phig_n)=phigarray(j,i)
                   enddo
              endif
          enddo
          !copy tempphig to phig
          do i=1, phig_n
             do k=1, lz1
                do n=1, ly1
                    do m=1, lx1
                        ind=m+(n-1)*ly1+(k-1)*lz1*ly1 
                        phig(m,n,k,i) = tempphig(ind, i)
                    enddo
                enddo
             enddo
          enddo
          !if ((istep .eq. 10) .and. (nid .eq. 0)) then
              !OPEN(UNIT=9999,FILE="phig_old_p0.txt",FORM="FORMATTED",
     $         !STATUS="REPLACE",ACTION="WRITE")
          !do i=1, phig_n
             !do k=1, lz1
                !do n=1, ly1
                    !do m=1, lx1
            !WRITE(UNIT=9999, FMT=*) nid, m, n, k, i, phig(m,n,k,i)
                    !enddo
                !enddo
             !enddo
          !enddo
              !CLOSE(UNIT=9999)
          !endif
          !if ((istep .eq. 10) .and. (nid .eq. 1)) then
              !OPEN(UNIT=10000,FILE="phig_old_p1.txt",FORM="FORMATTED",
     $         !STATUS="REPLACE",ACTION="WRITE")
          !do i=1, phig_n
             !do k=1, lz1
                !do n=1, ly1
                    !do m=1, lx1
            !WRITE(UNIT=10000, FMT=*) nid, m, n, k, i, phig(m,n,k,i)
                    !enddo
                !enddo
             !enddo
          !enddo
              !CLOSE(UNIT=10000)
          !endif
c         print *, "Update Phig array: phig_n", phig_n 
          end
c----------------------------------------------------------------------
      subroutine reinitialize
      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
      include 'ZPER'
c
      include 'OPCTR'
      include 'CTIMER'

      logical ifemati,ifsync_
      common /elementload/ gfirst, inoassignd, resetFindpts, pload(lelg)
      integer gfirst, inoassignd, resetFindpts, pload

      inoassignd = 0
      call readat
      etims0 = dnekclock_sync()
      if (nio.eq.0) then
         write(6,12) 'nelgt/nelgv/lelt:',nelgt,nelgv,lelt
         write(6,12) 'lx1  /lx2  /lx3 :',lx1,lx2,lx3
         write(6,'(A,g13.5,A,/)')  ' done :: read .rea file ',
     &                             etims0-etime,' sec'
 12      format(1X,A,4I12,/,/)
      endif

      ifsync_ = ifsync
      ifsync = .true.

      call setvar          ! Initialize most variables !skip 

!#ifdef MOAB
!      if (ifmoab) call nekMOAB_bcs  !   Map BCs
!#endif

      instep=1             ! Check for zero steps
      if (nsteps.eq.0 .and. fintim.eq.0.) instep=0

      igeom = 2
      call setup_topo      ! Setup domain topology  

      call genwz           ! Compute GLL points, weights, etc.

      call io_init         ! Initalize io unit

!      if (ifcvode.and.nsteps.gt.0)
!     $   call cv_setsize(0,nfield) !Set size for CVODE solver

      if(nio.eq.0) write(6,*) 'call usrdat'
      call usrdat
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat'
      call gengeom(igeom)  ! Generate geometry, after usrdat 

      if (ifmvbd) call setup_mesh_dssum ! Set mesh dssum (needs geom)

      if(nio.eq.0) write(6,*) 'call usrdat2'
      call usrdat2
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat2'

      call geom_reset(1)    ! recompute Jacobians, etc.
      call vrdsmsh          ! verify mesh topology

!      call echopar ! echo back the parameter stack
      call setlog  ! Initalize logical flags

      call bcmask  ! Set BC masks for Dirichlet boundaries.

      if (fintim.ne.0.0.or.nsteps.ne.0)
     $   call geneig(igeom) ! eigvals for tolerances

      call vrdsmsh     !     Verify mesh topology

      call dg_setup    !     Setup DG, if dg flag is set.

      !if (ifflow.and.(fintim.ne.0.or.nsteps.ne.0)) then    ! Pressure solver 
         !call estrat                                       ! initialization.
         !if (iftran.and.solver_type.eq.'itr') then         ! Uses SOLN space 
         !   call set_overlap                               ! as scratch!
         !elseif (solver_type.eq.'fdm'.or.solver_type.eq.'pdm')then
         !   ifemati = .true.
         !   kwave2  = 0.0
         !   if (ifsplit) ifemati = .false.
         !   call gfdm_init(nx2,ny2,nz2,ifemati,kwave2)
         !elseif (solver_type.eq.'25D') then
         !   call g25d_init
         !endif
      !endif

!      call init_plugin !     Initialize optional plugin
      if(ifcvode) call cv_setsize

      if(nio.eq.0) write(6,*) 'call usrdat3'
      call usrdat3
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat3'

#ifdef CMTNEK
        call nek_cmt_init
        if (nio.eq.0) write(6,*)'Initialized DG machinery'
#endif

c     call cmt_switch          ! Check if compiled with cmt
c     if (ifcmt) then          ! Initialize CMT branch
c       call nek_cmt_init
c       if (nio.eq.0) write(6,*)'Initialized DG machinery'
c     endif

c         for verification !!!!! 
          call outpost(vx, vy, vz, pr, t, 'xyz')
          call copy(vx, phig, nxyz*nelt)
c         for verification !!!!! 
       end

c------------------------------------------------------------------------
c this function is copyed from David new nek5000 code.
c It is called in the init_stokes_particles
c Sice right now is running David's example 3dbox_back, just rename zufalli and zufall to make them uncalled since they are in the file 3dbox_back
      subroutine zufalli_unused(seed)
      implicit none
c
c  generates initial seed buffer by linear congruential
c  method. Taken from Marsaglia, FSU report FSU-SCRI-87-50
c  variable seed should be 0 < seed <31328
c
      integer seed
      integer ptr
      double precision s,t
      double precision buff(607)
      integer ij,kl,i,ii,j,jj,k,l,m
      common /klotz0/buff,ptr
      data ij/1802/,kl/9373/
c
      if(seed.ne.0) ij = seed
c
      i = mod(ij/177,177) + 2
      j = mod(ij,177) + 2
      k = mod(kl/169,178) + 1
      l = mod(kl,169)
      do 1 ii=1,607
         s = 0.0
         t = 0.5
         do 2 jj=1,24
            m = mod(mod(i*j,179)*k,179)
            i = j
            j = k
            k = m
            l = mod(53*l+1,169)
            if(mod(l*m,64).ge.32) s = s+t
            t = .5*t
2        continue
         buff(ii) = s
1     continue
      
      return 
      end

c------------------------------------------------------------------------


      subroutine zufall_unused(n,a)
      implicit none
c
c portable lagged Fibonacci series uniform random number
c generator with "lags" -273 und -607:
c
c       t    = u(i-273)+buff(i-607)  (floating pt.)
c       u(i) = t - float(int(t))
c
c W.P. Petersen, IPS, ETH Zuerich, 19 Mar. 92
c
      double precision a(*)
      double precision buff(607)
      double precision t
      integer i,k,ptr,VL,k273,k607
      integer buffsz,nn,n,left,q,qq
      integer aptr,aptr0,bptr
c
      common /klotz0/buff,ptr
      data buffsz/607/
c
      aptr = 0
      nn   = n
c
1     continue
c
      if(nn .le. 0) return
c
c factor nn = q*607 + r
c
      q    = (nn-1)/607
      left = buffsz - ptr
c
      if(q .le. 1) then
c
c only one or fewer full segments
c
         if(nn .lt. left) then
            do 2 i=1,nn
               a(i+aptr) = buff(ptr+i)
2           continue
            ptr  = ptr + nn
            return
         else
            do 3 i=1,left
               a(i+aptr) = buff(ptr+i)
3           continue
            ptr  = 0
            aptr = aptr + left
            nn   = nn - left
c  buff -> buff case
            VL   = 273
            k273 = 334
            k607 = 0
            do 4 k=1,3
cdir$ ivdep
*vdir nodep
*VOCL LOOP, TEMP(t), NOVREC(buff)
               do 5 i=1,VL
                  t            = buff(k273+i) + buff(k607+i)
                  buff(k607+i) = t - float(int(t))
5              continue
               k607 = k607 + VL
               k273 = k273 + VL
               VL   = 167
               if(k.eq.1) k273 = 0
4           continue
c
            goto 1
         endif
      else
c
c more than 1 full segment
c
          do 6 i=1,left
             a(i+aptr) = buff(ptr+i)
6         continue
          nn   = nn - left
          ptr  = 0
          aptr = aptr+left
c
c buff -> a(aptr0)
c
          VL   = 273
          k273 = 334
          k607 = 0
          do 7 k=1,3
             if(k.eq.1)then
*VOCL LOOP, TEMP(t)
                do 8 i=1,VL
                   t         = buff(k273+i) + buff(k607+i)
                   a(aptr+i) = t - float(int(t))
8               continue
                k273 = aptr
                k607 = k607 + VL
                aptr = aptr + VL
                VL   = 167
             else
cdir$ ivdep
*vdir nodep
*VOCL LOOP, TEMP(t)
                do 9 i=1,VL
                   t         = a(k273+i) + buff(k607+i)
                   a(aptr+i) = t - float(int(t))
9               continue
                k607 = k607 + VL
                k273 = k273 + VL
                aptr = aptr + VL
             endif
7         continue
          nn = nn - 607
c
c a(aptr-607) -> a(aptr) for last of the q-1 segments
c
          aptr0 = aptr - 607
          VL    = 607
c
*vdir novector
          do 10 qq=1,q-2
             k273 = 334 + aptr0
cdir$ ivdep
*vdir nodep
*VOCL LOOP, TEMP(t), NOVREC(a)
             do 11 i=1,VL
                t         = a(k273+i) + a(aptr0+i)
                a(aptr+i) = t - float(int(t))
11           continue
             nn    = nn - 607
             aptr  = aptr + VL
             aptr0 = aptr0 + VL
10        continue
c
c a(aptr0) -> buff, last segment before residual
c
          VL   = 273
          k273 = 334 + aptr0
          k607 = aptr0
          bptr = 0
          do 12 k=1,3
             if(k.eq.1) then
*VOCL LOOP, TEMP(t)
                do 13 i=1,VL
                   t            = a(k273+i) + a(k607+i)
                   buff(bptr+i) = t - float(int(t))
13              continue
                k273 = 0
                k607 = k607 + VL
                bptr = bptr + VL
                VL   = 167
             else
cdir$ ivdep
*vdir nodep
*VOCL LOOP, TEMP(t), NOVREC(buff)
                do 14 i=1,VL
                   t            = buff(k273+i) + a(k607+i)
                   buff(bptr+i) = t - float(int(t))
14              continue
                k607 = k607 + VL
                k273 = k273 + VL
                bptr = bptr + VL
             endif
12        continue
          goto 1
      endif
      end



c------------------------------------------------------------------------
        subroutine printVerify
            include 'SIZE'
            common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
            
            print *, 'nid: ', nid_, 'nelt: ', nelt
        end




c--------------------------------------------------------------------
c     not called
      subroutine get_vert_map_again(vertex, nlv, nel, suffix, ifgfdm)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      logical ifgfdm
      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
      integer vertex(nlv,1)
      character*4 suffix

      parameter(mdw=2+2**ldim)
      parameter(ndw=7*lx1*ly1*lz1*lelv/mdw)
      common /scrns/ wk(mdw,ndw)   ! room for long ints, if desired
      integer wk,e,eg,eg0,eg1

      character*132 mapfle
      character*1   mapfle1(132)
      equivalence  (mapfle,mapfle1)

      iok = 0
      if (nid.eq.0) then
         lfname = ltrunc(reafle,132) - 4
         call blank (mapfle,132)
         call chcopy(mapfle,reafle,lfname)
         call chcopy(mapfle1(lfname+1),suffix,4)
         open(unit=80,file=mapfle,status='old',err=99)
         read(80,*,err=99) neli,nnzi
         iok = 1
      endif
   99 continue
      iok = iglmax(iok,1)
      if (iok.eq.0) goto 999     ! Mapfile not found

      if (nid.eq.0) then
         neli = iglmax(neli,1)   ! communicate to all procs
      else
         neli = 0
         neli = iglmax(neli,1)   ! communicate neli to all procs
      endif

      npass = 1 + (neli/ndw)
      if (npass.gt.np) then
         if (nid.eq.0) write(6,*) npass,np,neli,ndw,'Error get_vert_map'
         call exitt
      endif

      len = 4*mdw*ndw
      if (nid.gt.0.and.nid.lt.npass) msg_id=irecv(nid,wk,len)
      call nekgsync

      if (nid.eq.0) then
         eg0 = 0
         do ipass=1,npass
            eg1 = min(eg0+ndw,neli)
            m   = 0
            do eg=eg0+1,eg1
               m = m+1
               read(80,*,end=998) (wk(k,m),k=2,mdw)
               if(.not.ifgfdm)  gllnid(eg) = wk(2,m)  !proc map,  must still be divided
               wk(1,m)    = eg
            enddo
            if (ipass.lt.npass) call csend(ipass,wk,len,ipass,0) !send to ipass
            eg0 = eg1
         enddo
         close(80)
         ntuple = m
      elseif (nid.lt.npass) then
         call msgwait(msg_id)
         ntuple = ndw
      else
         ntuple = 0
      endif

      if (.not.ifgfdm) then             ! gllnid is already assigned for gfdm
        lng = isize*neli
        call bcast(gllnid,lng)
        !call assign_gllnid(gllnid,gllel,nelgt,nelgv,np) ! gllel is used as scratch
        !if(nid .eq. 0) then  !keke add
        !call assign_partitions   !(gllnid, lelt, nelgt, np) !keke add, assign gllnid according to the elements load balance
           !endif 
c          if(nid.eq.0) then
c         write(99,*) (gllnid(i),i=1,nelgt)
c       endif
c       call exitt
        call recompute_partitions
      endif

      nelt=0 !     Count number of elements on this processor
      nelv=0
      do eg=1,neli
         if (gllnid(eg).eq.nid) then
            if (eg.le.nelgv) nelv=nelv+1
            if (eg.le.nelgt) nelt=nelt+1
         endif
      enddo
      if (np.le.64) write(6,*) nid,nelv,nelt,nelgv,nelgt,' NELV'

c     NOW: crystal route vertex by processor id

      do i=1,ntuple
         eg=wk(1,i)
         wk(2,i)=gllnid(eg)        ! processor id for element eg
      enddo

      key = 2  ! processor id is in wk(2,:)
      call crystal_ituple_transfer(cr_h,wk,mdw,ntuple,ndw,key)

      if (.not.ifgfdm) then            ! no sorting for gfdm?
         key = 1  ! Sort tuple list by eg := wk(1,:)
         nkey = 1
         call crystal_ituple_sort(cr_h,wk,mdw,nelt,key,nkey)
      endif
      iflag = 0
      if (ntuple.ne.nelt) then
         write(6,*) nid,ntuple,nelv,nelt,nelgt,' NELT FAIL'
         write(6,*) 'Check that .map file and .rea file agree'
         iflag=1
      else
         nv = 2**ndim
         do e=1,nelt
            call icopy(vertex(1,e),wk(3,e),nv)
         enddo
      endif

      iflag = iglmax(iflag,1)
      if (iflag.gt.0) then
         do mid=0,np-1
            call nekgsync
            if (mid.eq.nid)
     $      write(6,*) nid,ntuple,nelv,nelt,nelgt,' NELT FB'
            call nekgsync
         enddo
         call nekgsync
         call exitt
      endif

      return

  999 continue
      if (nid.eq.0) write(6,*) 'ABORT: Could not find map file ',mapfle
      call exitt

  998 continue
      if (nid.eq.0) write(6,*)ipass,npass,eg0,eg1,mdw,m,eg,'get v fail'
      call exitt0  ! Emergency exit

      return
      end

