c     n = length(gdat) i.e. number of markers and log-ratio
c     k is the neighborhood width
c     oSD = how far the outlier is to the nearest observation
c     sSD = how close should it be moved
c     nchr = number of chromosomes
c     cfrq = number of probes in respective chromosomes
      subroutine smoothLR(n, gdat, nchr, cfrq, sgdat, k, oSD, sSD)
      integer n, nchr, cfrq(nchr), k
      double precision gdat(n), sgdat(n), oSD, sSD

      integer i, j, ilo, ihi, k1, k2, j1, ic, cilo, cihi
      double precision mnnbd, mxnbd, distij, xmed

c     temporary array for finding median
      double precision, allocatable :: xnbhd(:)

      k2 = 2*k + 1
      allocate(xnbhd(k2))

c     initial values for start and end of chromosomes
      cilo = 1
      cihi = 0
c     loop over chromomsomes
      do 100 ic = 1, nchr
c     end of the current chromosome
         cihi = cihi + cfrq(ic)
         do 50 i = cilo, cihi
c     range of the neighborhood
            ilo = max(cilo, i-k)
            ihi = min(cihi, i+k)
c     check if ith observation is an outlier
c     initialize the distances to be large
            mxnbd = 100*oSD
            mnnbd = 100*oSD
            do 10 j = ilo, ihi
               if (j .ne. i) then
c     calculate distance from between ith and jth obsn
                  distij = gdat(i) - gdat(j)
c     if distance is less than oSD no smoothing necessary
                  if (abs(distij) .le. oSD) then
                     sgdat(i) = gdat(i)
                     go to 50
c     otherwise calculate distances from above and below
                  else
c     mxnbd is distance from above
                     if (distij .lt. mxnbd) mxnbd = distij
c     mnnbd is distance from below
                     if (-distij .lt. mnnbd) mnnbd = -distij
                  endif
               endif
 10         continue
c     if all the points in the nbhd are above than mxnbd will be negative
c     and mnnbd will be greater than oSD. Vice versa if all points below
c     
c     If both are negative then the ith point is singleton in the middle 
c     but distance oSD away from all points in the nbhd. No smoothing done.
            if ((mxnbd .le. 0) .and. (mnnbd .le. 0)) then
               sgdat(i) = gdat(i)
               go to 50
            else
c     calculate the median of the nbhd
c     number of points in the nbhd
               k1 = ihi - ilo + 1
c     get the data into temporary array
               do 20 j = ilo, ihi
                  xnbhd(j-ilo+1) = gdat(j)
 20            continue
c     sort the data
               call qsort3(xnbhd, 1, k1)
c     median is the midpoint if n is odd and average of the two if even
               j1 = k1/2
               if (k1 .eq. 2*j1) then
                  xmed = (xnbhd(j1) + xnbhd(j1+1))/2
               else
                  xmed = xnbhd(j1+1)
               endif
c     if point is above the nbhd bring it down
               if (mxnbd .gt. 0) sgdat(i) = xmed + sSD
c     if point is below the nbhd bring it up
               if (mnnbd .gt. 0) sgdat(i) = xmed - sSD
            endif
 50      continue
c     beginning of next chromosome
         cilo = cilo + cfrq(ic)
 100  continue
      deallocate(xnbhd)

      return
      end
	  
      subroutine getbdry(eta, m, nperm, mb, ibdry, etastr, tol)
      integer m, nperm, mb, ibdry(mb)
      double precision eta, etastr(m), tol

      double precision eta0, etalo, etahi, plo, phi, pexcd
      integer j, l

      l = 1
      ibdry(1) = nperm-int(dfloat(nperm)*eta)
      etastr(1) = eta
      eta0 = eta
      do 20 j = 2,m
         etahi = eta0*1.1
         call etabdry(nperm, etahi, j, ibdry(l+1))
         call pexceed(nperm, j, ibdry(l+1), phi)
         etalo = eta0*0.25
         call etabdry(nperm, etalo, j, ibdry(l+1))
         call pexceed(nperm, j, ibdry(l+1), plo)
         do 10 while ((etahi-etalo)/etalo .gt. tol)
            eta0 = etalo + (etahi-etalo)*(eta-plo)/(phi-plo)
            call etabdry(nperm, eta0, j, ibdry(l+1))
            call pexceed(nperm, j, ibdry(l+1), pexcd)
            if (pexcd .gt. eta) then
               etahi = eta0
               phi = pexcd
            else
               etalo = eta0
               plo = pexcd
            endif
 10      continue
         etastr(j) = eta0
         l = l+j
 20   continue

      return
      end

      subroutine etabdry(nperm, eta0, n1s, ibdry)
      integer nperm, n1s, ibdry(n1s)
      double precision eta0

      double precision fphypr
      external fphypr

      integer i, k
      double precision di, dn, dn1s, dk, tprob

      dn1s = dfloat(n1s)
      dn = dfloat(nperm-n1s)
      
      k = 0
      dk = 0.0d0
      do 10 i = 1, nperm
         di = dfloat(i)
         tprob = fphypr(dk, dn1s, dn, di)
         if (tprob .le. eta0) then
            k = k+1
            dk = dk + 1.0d0
            ibdry(k) = i
         endif
 10   continue

      return
      end

      subroutine pexceed(nperm, n1s, ibdry, pexcd)
      integer nperm, n1s, ibdry(n1s)
      double precision pexcd

      double precision dn, dk, dn1, dk1, dn2, dk2, dn3, dk3, dlcnk
      integer i

      double precision flchoose
      external flchoose

      dn = dfloat(nperm)
      dk = dfloat(n1s)
      dn1 = dfloat(nperm-ibdry(1))
      dlcnk = flchoose(dn, dk)

      pexcd = exp(flchoose(dn1, dk) - dlcnk)

      if (n1s .ge. 2) then
         dn1 = dfloat(ibdry(1))
         dn = dfloat(nperm-ibdry(2))
         dk = dfloat(n1s-1)
         pexcd = pexcd + exp(log(dn1) + flchoose(dn, dk) - dlcnk)
      endif

      if (n1s .ge. 3) then
         dn1 = dfloat(ibdry(1))
         dn2 = dfloat(ibdry(2))
         dn = dfloat(nperm-ibdry(3))
         dk = dfloat(n1s-2)
         pexcd = pexcd + 
     1        exp(log(dn1) + log(dn1-1.0) - log(2.0) + 
     2                        flchoose(dn, dk) - dlcnk) +
     3        exp(log(dn1) + log(dn2-dn1) + flchoose(dn, dk) - dlcnk)
      endif

      if (n1s .gt. 3) then
         do 10 i = 4, n1s
            dn1 = dfloat(ibdry(i-3))
            dk1 = dfloat(i-1)
            dk2 = dfloat(i-2)
            dk3 = dfloat(i-3)
            dn2 = dfloat(ibdry(i-2))
            dn3 = dfloat(ibdry(i-1))
            dn = dfloat(nperm-ibdry(i))
            dk = dfloat(n1s-i+1)
            pexcd = pexcd + 
     1           exp(flchoose(dn1, dk1) + flchoose(dn, dk) - dlcnk) +
     2           exp(flchoose(dn1, dk2) + log(dn3-dn1) + 
     3                       flchoose(dn, dk) - dlcnk) +
     4           exp(flchoose(dn1, dk3) + log(dn2-dn1) + log(dn3-dn2) +
     3                        flchoose(dn, dk) - dlcnk) +
     5           exp(flchoose(dn1, dk3) + log(dn2-dn1) - log(2.0) + 
     6                    log(dn2-dn1-1.0) + flchoose(dn, dk) - dlcnk)
 10      continue
      endif

      return
      end

c     new approach to maximizing t-statistic
c     dynamic memory allocation using allocatable arrays 
      subroutine tmaxo(n,x,tss,sx,iseg,ostat,al0,ibin)
      integer n,iseg(2),al0
      double precision x(n),tss,sx(n),ostat
      logical ibin
c     
c     look at the partial sums in blocks of size sqrt(n)
c     
      integer ipsmin, ipsmax, ipsmin0, ipsmax0, nb, i, j, k, l, nb1,
     1     nb2, bi, bj, ilo, ihi, jlo, jhi, alenmax, i2j, sxmxi,
     2     alenlo, alenhi, tmaxi, tmaxj, ixlo, ixhi, nal0
      double precision psum, psmin, psmax, psmin0, psmax0, bssmax,
     1     bsslim, rn, rj, rjhi, rjlo, rnjov1, sij1, sij2, sijmx0,
     2     absx, sxmx, bijbss, rnov2, psdiff
c     
c     use local arrays for working within blocks
c     block partial sum max and min 
      double precision, allocatable :: bpsmax(:), bpsmin(:)
c     location of the max and min
      integer, allocatable :: bb(:), ibmin(:), ibmax(:)

c     t statistic corresponding to max for block i,j (and max possible)
      double precision, allocatable :: bssbij(:), bssijmax(:)
c     row, column and order vector for reordering bssbij
      integer, allocatable :: bloci(:), blocj(:), loc(:), alen(:)

c     calculate number of blocks (nb) and block boundaries (vector bb)
      rn = dfloat(n)
      if (n .ge. 50) then
         nb = nint(sqrt(dfloat(n)))
      else
         nb = 1
      endif

c     the number of paiwise block comparison
      nb2 = nb*(nb+1)/2
c     allocate memory
      allocate(bpsmax(nb), bpsmin(nb))
      allocate(bb(nb), ibmin(nb), ibmax(nb))
      allocate(bssbij(nb2), bssijmax(nb2))
      allocate(bloci(nb2), blocj(nb2), loc(nb2), alen(nb2))

c     block boundaries
      do 110 i = 1, nb
         bb(i) = nint(rn*(dfloat(i)/dfloat(nb)))
 110  continue

c     find the max, min of partial sums and their locations within blocks
      ilo = 1
      psum = 0
      psmin0 = 0
      psmax0 = 0
      ipsmin0 = n
      ipsmax0 = n
      do 20 j = 1, nb
         sx(ilo) = psum + x(ilo)
         psmin = sx(ilo)
         ipsmin = ilo
         psmax = sx(ilo)
         ipsmax = ilo
         do 10 i = ilo+1, bb(j)
            sx(i) = sx(i-1) + x(i)
            if (sx(i) .lt. psmin) then 
               psmin = sx(i)
               ipsmin = i
            endif
            if (sx(i) .gt. psmax) then 
               psmax = sx(i)
               ipsmax = i
            endif
 10      continue
c     store the block min, max and locations
         ibmin(j) = ipsmin
         ibmax(j) = ipsmax
         bpsmin(j) = psmin
         bpsmax(j) = psmax
c     adjust global min, max and locations
         if (psmin .lt. psmin0) then
            psmin0 = psmin
            ipsmin0 = ipsmin
         endif
         if (psmax .gt. psmax0) then
            psmax0 = psmax
            ipsmax0 = ipsmax
         endif
c     reset ilo to be the block boundary + 1
         psum = sx(bb(j))
         ilo = bb(j) + 1
 20   continue

c     calculate bss for max s_i - min s_i
      psdiff = psmax0 - psmin0
      rj = dfloat(abs(ipsmax0 - ipsmin0))
      rnjov1 = rn/(rj*(rn-rj))
      if (ibin) then
         bssmax = rnjov1*(psdiff-0.5)**2
      else
         bssmax = rnjov1*psdiff**2
      endif
      tmaxi = min(ipsmax0, ipsmin0)
      tmaxj = max(ipsmax0, ipsmin0)

c     if the segment is all constant then psdiff = 0 and so bssmax = 0
      if (psdiff .le. 0) then
         bssmax = 0
         go to 120
      endif

c     for a pair of blocks (i,j) calculate the max absolute t-statistic
c     at the (min_i, max_j) and (max_i, min_j) locations 
c     for other indices the t-statistic can be bounded using this
c
c     if a block doesn't have the potential to exceed bssmax ignore it
c     calculate the bsslim for each block and include ones >= bssmax

      rnov2 = rn/2
      l = 0
      nal0 = n - al0
      do 40 i = 1, nb
         do 30 j = i, nb
c     calculate bsslim
            if (i .eq. 1) then
               ilo = 1
            else
               ilo = bb(i-1) + 1
            endif
            ihi = bb(i)
            if (j .eq. 1) then
               jlo = 1
            else
               jlo = bb(j-1) + 1
            endif
            jhi = bb(j)
            alenhi = jhi - ilo
            if (alenhi .gt. nal0) alenhi = nal0
            rjhi = dfloat(alenhi)
            if (i .eq. j) then
               alenlo = 1
            else
               alenlo = jlo - ihi
            endif
            if (alenlo .lt. al0) alenlo = al0
c     max S_k over block j - min S_k over block i
            sij1 = abs(bpsmax(j) - bpsmin(i))
c     max S_k over block i - min S_k over block j
            sij2 = abs(bpsmax(i) - bpsmin(j))
c     if i = j then sij1 and sij2 are the same
            sijmx0 = max(sij1, sij2)
            rjlo = dfloat(alenlo)
            rnjov1 = rn/min(rjlo*(rn-rjlo), rjhi*(rn-rjhi))
            if (ibin) then
               bsslim = rnjov1*(sijmx0-0.5)**2
            else
               bsslim = rnjov1*(sijmx0**2)
            endif
c     if its as large as bssmax add block
            if (bssmax .le. bsslim) then
               l = l+1
               loc(l) = l
               bloci(l) = i
               blocj(l) = j
               bssijmax(l) = bsslim
c     max sij in the (i,j) block, t-statistic etc
               if (sij1 .gt. sij2) then
                  alen(l) = abs(ibmax(j) - ibmin(i))
                  rj = dfloat(alen(l))
                  rnjov1 = rn/(rj*(rn-rj))
                  if (ibin) then
                     bssbij(l) = rnjov1*(sij1-0.5)**2
                  else
                     bssbij(l) = rnjov1*(sij1**2)
                  endif
               else
                  alen(l) = abs(ibmin(j) - ibmax(i))
                  rj = dfloat(alen(l))
                  rnjov1 = rn/(rj*(rn-rj))
                  if (ibin) then
                     bssbij(l) = rnjov1*(sij2-0.5)**2
                  else
                     bssbij(l) = rnjov1*(sij2**2)
                  endif
               endif
            endif
 30      continue
 40   continue
      nb1 = l

c     Now sort the t-statistics by their magnitude
      call qsort4(bssbij, loc, 1, nb1)

c     now go through the blocks in reverse order (largest down)
      do 100 l = nb1, 1, -1
         k = loc(l)
c     need to check a block only if it has potential to increase bss
c     rjlo is the smalllest (j-i) in the block and rjhi is the largest
         bsslim = bssijmax(k)
         if (bssmax .le. bsslim) then
c     bi, bj give the block location
            bi = bloci(k)
            bj = blocj(k)
c     max arc length of interest in block
            alenmax = alen(k)
            if (bi .eq. 1) then
               ilo = 1
            else
               ilo = bb(bi-1) + 1
            endif
            ihi = bb(bi)
            if (bj .eq. 1) then
               jlo = 1
            else
               jlo = bb(bj-1) + 1
            endif
            jhi = bb(bj)
            alenhi = jhi - ilo
            if (alenhi .gt. nal0) alenhi = nal0
            rjhi = dfloat(alenhi)
            if (bi .eq. bj) then
               alenlo = 1
            else
               alenlo = jlo - ihi
            endif
            if (alenlo .lt. al0) alenlo = al0
            rjlo = dfloat(alenlo)
c
c     if arc length is larger than n/2 make is n - arc length
c
            if (alenmax .gt. n - alenmax) alenmax = n - alenmax
c
c     if alenlo <= n/2 start from (ihi, jlo) and go up
c     if alenhi >= n/2 start from (ilo, jhi) and go down
c
            if ((rjlo .le. rnov2) .and. (alenlo .le. alenmax)) then
               do 60 i2j = alenlo, alenmax
c     excess calcultaions to set range of i
                  ixlo = max(0, jlo - ilo - i2j)
                  ixhi = max(0, ihi + i2j - jhi)
                  sxmx = 0
                  do 55 i = ilo + ixlo, ihi - ixhi
                     j = i+i2j
                     absx = abs(sx(j) - sx(i))
                     if (sxmx .lt. absx) then
                        sxmx = absx
                        sxmxi = i
                     endif
 55               continue
                  rj = dfloat(i2j)
                  rnjov1 = rn/(rj*(rn-rj))
                  if (ibin) then
                     bijbss = rnjov1*(sxmx-0.5)**2
                  else
                     bijbss = rnjov1*(sxmx**2)
                  endif
                  if (bijbss .gt. bssmax) then
                     bssmax = bijbss
                     tmaxi = sxmxi
                     tmaxj = sxmxi + i2j
                  endif
 60            continue
            endif
c
c     make arclength n - arc length
c
            alenmax = n - alenmax
            if ((rjhi .ge. rnov2) .and. (alenhi .ge. alenmax)) then
               do 70 i2j = alenhi, alenmax, -1
c     excess calcultaions to set range of i
                  ixlo = max(0, jlo - ilo - i2j)
                  ixhi = max(0, ihi + i2j - jhi)
                  sxmx = 0
                  do 65 i = ilo + ixlo, ihi - ixhi
                     j = i + i2j
                     absx = abs(sx(j) - sx(i))
                     if (sxmx .lt. absx) then
                        sxmx = absx
                        sxmxi = i
                     endif
 65               continue
                  rj = dfloat(i2j)
                  rnjov1 = rn/(rj*(rn-rj))
                  if (ibin) then
                     bijbss = rnjov1*(sxmx-0.5)**2
                  else
                     bijbss = rnjov1*(sxmx**2)
                  endif
                  if (bijbss .gt. bssmax) then
                     bssmax = bijbss
                     tmaxi = sxmxi
                     tmaxj = sxmxi + i2j
                  endif
 70            continue
            endif
         endif
 100  continue

 120  if (ibin) then
         if (tss.le.0.0001) tss = 1.0
         bssmax = bssmax/(tss/rn)
      else
         if (tss.le.bssmax+0.0001) tss = bssmax + 1.0
         bssmax = bssmax/((tss-bssmax)/(rn-2.0))
      endif

c     deallocate memory
      deallocate(bpsmax, bpsmin, bb, ibmin, ibmax)
      deallocate(bssbij, bssijmax, bloci, blocj, loc, alen)

      ostat = bssmax
      iseg(1) = tmaxi
      iseg(2) = tmaxj

      return
      end

c     function for calculating the full max t-statistic on permuted data
c     new approach to maximizing t-statistic using allocatable arrays 
      double precision function tmaxp(n,tss,px,sx,al0,ibin)
      integer n,al0
      double precision tss,px(n),sx(n)
      logical ibin
c     
c     look at the partial sums in blocks of size sqrt(n)
c     
      integer ipsmin, ipsmax, ipsmin0, ipsmax0, nb, i, j, k, l, nb1,
     1     nb2, bi, bj, ilo, ihi, jlo, jhi, alenmax, i2j, alenlo,
     2     alenhi, ixlo, ixhi, nal0
      double precision psum, psmin, psmax, psmin0, psmax0, bssmax,
     1     bsslim, rn, rj, rjhi, rjlo, rnjov1, sij1, sij2, sijmx0,
     2     absx, sxmx, bijbss, rnov2, psdiff
c     
c     use local arrays for working within blocks
c     block partial sum max and min 
      double precision, allocatable :: bpsmax(:), bpsmin(:)
c     location of the max and min
      integer, allocatable :: bb(:), ibmin(:), ibmax(:)

c     t statistic corresponding to max for block i,j (and max possible)
      double precision, allocatable :: bssbij(:), bssijmax(:)
c     row, column and order vector for reordering bssbij
      integer, allocatable :: bloci(:), blocj(:), loc(:), alen(:)

c     calculate number of blocks (nb) and block boundaries (vector bb)
      rn = dfloat(n)
      if (n .ge. 50) then
         nb = nint(sqrt(dfloat(n)))
      else
         nb = 1
      endif

c     the number of paiwise block comparison
      nb2 = nb*(nb+1)/2
c     allocate memory
      allocate(bpsmax(nb), bpsmin(nb))
      allocate(bb(nb), ibmin(nb), ibmax(nb))
      allocate(bssbij(nb2), bssijmax(nb2))
      allocate(bloci(nb2), blocj(nb2), loc(nb2), alen(nb2))

c     block boundaries
      do 110 i = 1, nb
         bb(i) = nint(rn*(dfloat(i)/dfloat(nb)))
 110  continue

c     find the max, min of partial sums and their locations within blocks
      ilo = 1
      psum = 0
      psmin0 = 0
      psmax0 = 0
      ipsmin0 = n
      ipsmax0 = n
      do 20 j = 1, nb
         sx(ilo) = psum + px(ilo)
         psmin = sx(ilo)
         ipsmin = ilo
         psmax = sx(ilo)
         ipsmax = ilo
         do 10 i = ilo+1, bb(j)
            sx(i) = sx(i-1) + px(i)
            if (sx(i) .lt. psmin) then 
               psmin = sx(i)
               ipsmin = i
            endif
            if (sx(i) .gt. psmax) then 
               psmax = sx(i)
               ipsmax = i
            endif
 10      continue
c     store the block min, max and locations
         ibmin(j) = ipsmin
         ibmax(j) = ipsmax
         bpsmin(j) = psmin
         bpsmax(j) = psmax
c     adjust global min, max and locations
         if (psmin .lt. psmin0) then
            psmin0 = psmin
            ipsmin0 = ipsmin
         endif
         if (psmax .gt. psmax0) then
            psmax0 = psmax
            ipsmax0 = ipsmax
         endif
c     reset ilo to be the block boundary + 1
         psum = sx(bb(j))
         ilo = bb(j) + 1
 20   continue

c     calculate bss for max s_i - min s_i
      psdiff = psmax0 - psmin0
      rj = dfloat(abs(ipsmax0 - ipsmin0))
      rnjov1 = rn/(rj*(rn-rj))
      if (ibin) then
         bssmax = rnjov1*(psdiff-0.5)**2
      else
         bssmax = rnjov1*psdiff**2
      endif

c     for a pair of blocks (i,j) calculate the max absolute t-statistic
c     at the (min_i, max_j) and (max_i, min_j) locations 
c     for other indices the t-statistic can be bounded using this
c
c     if a block doesn't have the potential to exceed bssmax ignore it
c     calculate the bsslim for each block and include ones >= bssmax

      rnov2 = rn/2
      l = 0
      nal0 = n - al0
      do 40 i = 1, nb
         do 30 j = i, nb
c     calculate bsslim
            if (i .eq. 1) then
               ilo = 1
            else
               ilo = bb(i-1) + 1
            endif
            ihi = bb(i)
            if (j .eq. 1) then
               jlo = 1
            else
               jlo = bb(j-1) + 1
            endif
            jhi = bb(j)
            alenhi = jhi - ilo
            if (alenhi .gt. nal0) alenhi = nal0
            rjhi = dfloat(alenhi)
            if (i .eq. j) then
               alenlo = 1
            else
               alenlo = jlo - ihi
            endif
            if (alenlo .lt. al0) alenlo = al0
c     max S_k over block j - min S_k over block i
            sij1 = abs(bpsmax(j) - bpsmin(i))
c     max S_k over block i - min S_k over block j
            sij2 = abs(bpsmax(i) - bpsmin(j))
c     if i = j then sij1 and sij2 are the same
            sijmx0 = max(sij1, sij2)
            rjlo = dfloat(alenlo)
            rnjov1 = rn/min(rjlo*(rn-rjlo), rjhi*(rn-rjhi))
            if (ibin) then
               bsslim = rnjov1*(sijmx0-0.5)**2
            else
               bsslim = rnjov1*(sijmx0**2)
            endif
c     if its as large as bssmax add block
            if (bssmax .le. bsslim) then
               l = l+1
               loc(l) = l
               bloci(l) = i
               blocj(l) = j
               bssijmax(l) = bsslim
c     max sij in the (i,j) block, t-statistic etc
               if (sij1 .gt. sij2) then
                  alen(l) = abs(ibmax(j) - ibmin(i))
                  rj = dfloat(alen(l))
                  rnjov1 = rn/(rj*(rn-rj))
                  if (ibin) then
                     bssbij(l) = rnjov1*(sij1-0.5)**2
                  else
                     bssbij(l) = rnjov1*(sij1**2)
                  endif
               else
                  alen(l) = abs(ibmin(j) - ibmax(i))
                  rj = dfloat(alen(l))
                  rnjov1 = rn/(rj*(rn-rj))
                  if (ibin) then
                     bssbij(l) = rnjov1*(sij2-0.5)**2
                  else
                     bssbij(l) = rnjov1*(sij2**2)
                  endif
               endif
            endif
 30      continue
 40   continue
      nb1 = l

c     Now sort the t-statistics by their magnitude
      call qsort4(bssbij, loc, 1, nb1)

c     now go through the blocks in reverse order (largest down)
      do 100 l = nb1, 1, -1
         k = loc(l)
c     need to check a block only if it has potential to increase bss
c     rjlo is the smalllest (j-i) in the block and rjhi is the largest
         bsslim = bssijmax(k)
         if (bssmax .le. bsslim) then
c     bi, bj give the block location
            bi = bloci(k)
            bj = blocj(k)
c     max arc length of interest in block
            alenmax = alen(k)
            if (bi .eq. 1) then
               ilo = 1
            else
               ilo = bb(bi-1) + 1
            endif
            ihi = bb(bi)
            if (bj .eq. 1) then
               jlo = 1
            else
               jlo = bb(bj-1) + 1
            endif
            jhi = bb(bj)
            alenhi = jhi - ilo
            if (alenhi .gt. nal0) alenhi = nal0
            rjhi = dfloat(alenhi)
            if (bi .eq. bj) then
               alenlo = 1
            else
               alenlo = jlo - ihi
            endif
            if (alenlo .lt. al0) alenlo = al0
            rjlo = dfloat(alenlo)
c
c     if arc length is larger than n/2 make is n - arc length
c
            if (alenmax .gt. n - alenmax) alenmax = n - alenmax
c
c     if alenlo <= n/2 start from (ihi, jlo) and go up
c     if alenhi >= n/2 start from (ilo, jhi) and go down
c
            if ((rjlo .le. rnov2) .and. (alenlo .le. alenmax)) then
               do 60 i2j = alenlo, alenmax
c     excess calcultaions to set range of i
                  ixlo = max(0, jlo - ilo - i2j)
                  ixhi = max(0, ihi + i2j - jhi)
                  sxmx = 0
                  do 55 i = ilo + ixlo, ihi - ixhi
                     j = i+i2j
                     absx = abs(sx(j) - sx(i))
                     if (sxmx .lt. absx) sxmx = absx
 55               continue
                  rj = dfloat(i2j)
                  rnjov1 = rn/(rj*(rn-rj))
                  if (ibin) then
                     bijbss = rnjov1*(sxmx-0.5)**2
                  else
                     bijbss = rnjov1*(sxmx**2)
                  endif
                  if (bijbss .gt. bssmax) bssmax = bijbss
 60            continue
            endif
c
c     make arclength n - arc length
c
            alenmax = n - alenmax
            if ((rjhi .ge. rnov2) .and. (alenhi .ge. alenmax)) then
               do 70 i2j = alenhi, alenmax, -1
c     excess calcultaions to set range of i
                  ixlo = max(0, jlo - ilo - i2j)
                  ixhi = max(0, ihi + i2j - jhi)
                  sxmx = 0
                  do 65 i = ilo + ixlo, ihi - ixhi
                     j = i + i2j
                     absx = abs(sx(j) - sx(i))
                     if (sxmx .lt. absx) sxmx = absx
 65               continue
                  rj = dfloat(i2j)
                  rnjov1 = rn/(rj*(rn-rj))
                  if (ibin) then
                     bijbss = rnjov1*(sxmx-0.5)**2
                  else
                     bijbss = rnjov1*(sxmx**2)
                  endif
                  if (bijbss .gt. bssmax) bssmax = bijbss
 70            continue
            endif
         endif
 100  continue

      if (ibin) then
         if (tss.le.0.0001) tss = 1.0
         tmaxp = bssmax/(tss/rn)
      else
         if (tss.le.bssmax+0.0001) tss = bssmax + 1.0
         tmaxp = bssmax/((tss-bssmax)/(rn-2.0))
      endif

c     deallocate memory
      deallocate(bpsmax, bpsmin, bb, ibmin, ibmax)
      deallocate(bssbij, bssijmax, bloci, blocj, loc, alen)

      return
      end

c     function for the max (over small arcs) t-statistic on permuted data
c     new code to speed up this part 3/31/2010
      double precision function htmaxp(n,k,tss,px,sx,al0,ibin)
      integer n,k,al0
      double precision tss,px(n),sx(n)
      logical ibin

      integer i, j, nmj
      double precision rn, rj, absx, sxmx, bssmx, psmin, psmax, psdiff,
     1     bsslim, rnjov1

c     create blocks of size k (or k+1) to span 1 thru n
c     block partial sum max and min 
      double precision, allocatable :: bpsmax(:), bpsmin(:)
c     location of the max and min
      integer, allocatable :: bb(:)
c     variables to work on block specific data
      integer nb, ilo, ihi, l
      double precision psum, psdiffsq

      rn = dfloat(n)
c     number of blocks of size k (plus fraction since n/k may not be integer)
      nb = int(rn/dfloat(k))
c     allocate memory
      allocate(bpsmax(nb), bpsmin(nb))
      allocate(bb(nb))
c     block boundaries
      do 110 i = 1, nb
         bb(i) = nint(rn*(dfloat(i)/dfloat(nb)))
 110  continue

c     don't need global min and max
c     find the max, min of partial sums and their locations within blocks
      ilo = 1
      psum = 0
      htmaxp = 0.0d0
      do 20 j = 1, nb
         sx(ilo) = psum + px(ilo)
         psmin = sx(ilo)
         ipsmin = ilo
         psmax = sx(ilo)
         ipsmax = ilo
         do 10 i = ilo+1, bb(j)
            sx(i) = sx(i-1) + px(i)
            if (sx(i) .lt. psmin) then 
               psmin = sx(i)
               ipsmin = i
            endif
            if (sx(i) .gt. psmax) then 
               psmax = sx(i)
               ipsmax = i
            endif
 10      continue
c     store the block min, max and locations
         bpsmin(j) = psmin
         bpsmax(j) = psmax
c     reset ilo to be the block boundary + 1
         psum = sx(bb(j))
         ilo = bb(j) + 1
c     calculate the bss at the block max & min pr
         i = abs(ipsmin - ipsmax)
         if ((i .le. k) .and. (i .ge. al0)) then
            rj = dfloat(i)
            rnjov1 = rn/(rj*(rn-rj))
            if (ibin) then
               bssmx = rnjov1*(bpsmax(j) - bpsmin(j) -0.5)**2
            else
               bssmx = rnjov1*(bpsmax(j) - bpsmin(j))**2
            endif
            if (htmaxp .lt. bssmx) htmaxp = bssmx
         endif
 20   continue

c     check the first block
      ilo = 1
      ihi = bb(1)
      psdiff = bpsmax(1) - bpsmin(1)
      if (ibin) then
         psdiffsq = (psdiff-0.5)**2
      else
         psdiffsq = psdiff**2
      endif
      do 40 j = al0,k
         rj = dfloat(j)
         rnjov1 = rn/(rj*(rn-rj))
         bsslim = rnjov1*psdiffsq
         if (bsslim .lt. htmaxp) go to 50
         sxmx = 0.0d0
         do 30 i = ilo,ihi-j
            absx = abs(sx(i+j) - sx(i))
            if (sxmx.lt.absx) sxmx = absx
 30      continue
         if (ibin) then
            bssmx = rnjov1*(abs(sxmx)-0.5)**2
         else
            bssmx = rnjov1*sxmx**2
         endif
         if (htmaxp.lt.bssmx) htmaxp = bssmx
 40   continue

c     now the minor arcs spanning the end (n)
 50   psdiff = max(abs(bpsmax(1)-bpsmin(nb)), abs(bpsmax(nb)-bpsmin(1)))
      if (ibin) then
         psdiffsq = (psdiff-0.5)**2
      else
         psdiffsq = psdiff**2
      endif
      do 70 j = al0,k
         rj = dfloat(j)
         rnjov1 = rn/(rj*(rn-rj))
         bsslim = rnjov1*psdiffsq
         if (bsslim .lt. htmaxp) go to 100
         sxmx = 0.0d0
         nmj = n-j
         do 60 i = 1,j
            absx = abs(sx(i+nmj) - sx(i))
            if (sxmx.lt.absx) sxmx = absx
 60      continue
         if (ibin) then
            bssmx = rnjov1*(abs(sxmx)-0.5)**2
         else
            bssmx = rnjov1*sxmx**2
         endif
         if (htmaxp.lt.bssmx) htmaxp = bssmx
 70   continue

c     now the other blocks
 100  do 200 l = 2,nb
         ilo = bb(l-1)+1
         ihi = bb(l)
         psdiff = bpsmax(l) - bpsmin(l)
         if (ibin) then
            psdiffsq = (psdiff-0.5)**2
         else
            psdiffsq = psdiff**2
         endif
         do 140 j = al0,k
            rj = dfloat(j)
            rnjov1 = rn/(rj*(rn-rj))
            bsslim = rnjov1*psdiffsq
            if (bsslim .lt. htmaxp) go to 150
            sxmx = 0.0d0
            do 130 i = ilo,ihi-j
               absx = abs(sx(i+j) - sx(i))
               if (sxmx.lt.absx) sxmx = absx
 130        continue
            if (ibin) then
               bssmx = rnjov1*(abs(sxmx)-0.5)**2
            else
               bssmx = rnjov1*sxmx**2
            endif
            if (htmaxp.lt.bssmx) htmaxp = bssmx
 140     continue
 150     psdiff = max(abs(bpsmax(l)-bpsmin(l-1)), 
     1        abs(bpsmax(l-1)-bpsmin(l)))
         if (ibin) then
            psdiffsq = (psdiff-0.5)**2
         else
            psdiffsq = psdiff**2
         endif
         do 170 j = al0,k
            rj = dfloat(j)
            rnjov1 = rn/(rj*(rn-rj))
            bsslim = rnjov1*psdiffsq
            if (bsslim .lt. htmaxp) go to 200
            sxmx = 0.0d0
            nmj = n-j
            do 160 i = ilo-j,ilo-1
               absx = abs(sx(i+j) - sx(i))
               if (sxmx.lt.absx) sxmx = absx
 160        continue
            if (ibin) then
               bssmx = rnjov1*(abs(sxmx)-0.5)**2
            else
               bssmx = rnjov1*sxmx**2
            endif
            if (htmaxp.lt.bssmx) htmaxp = bssmx
 170     continue
 200  continue
      if (ibin) then
         if (tss .le. 0.0001d0) tss = 1.0d0
         htmaxp = htmaxp/(tss/rn)
      else
         if (tss .le. htmaxp+0.0001d0) tss = htmaxp + 1.0d0
         htmaxp = htmaxp/((tss-htmaxp)/(rn-2.0d0))
      endif

c     deallocate memory
      deallocate(bpsmax, bpsmin, bb)

      return
      end
	  
c     these are the subroutines to do the weigthed version of CBS
c     which is useful in order to merge data from multiple platforms
c     --------------------------------------------------------------
c         This is relevant only for log-ratio not binary data
c     --------------------------------------------------------------
c     function for calculating the full max weighted t-statistic
c     new approach to maximizing t-statistic

      subroutine wtmaxo(n,x,wts,tss,sx,cwts,iseg,ostat,al0)
      integer n,iseg(2),al0
      double precision x(n),wts(n),tss,sx(n),cwts(n),ostat
c     
c     look at the partial sums in blocks of size sqrt(n)
c     
      integer ipsmin, ipsmax, ipsmin0, ipsmax0, nb, i, j, k, l, nb1,
     1     nb2, bi, bj, ilo, ihi, jlo, jhi, ihi1, jlo1, jhi1,
     2     tmaxi, tmaxj, nal0
      double precision psum, psmin, psmax, psmin0, psmax0, bssmax,
     1     bsslim, rn, sij1, sij2, sijmx0, bijbss, awtmax, psrnov2,
     2     psdiff, psrj, psrn, psrnj, awtlo, awthi, awt1
c     
c     use local arrays for working within blocks
c     block partial sum max and min 
      double precision, allocatable :: bpsmax(:), bpsmin(:)
c     location of the max and min
      integer, allocatable :: bb(:), ibmin(:), ibmax(:)

c     t statistic corresponding to max for block i,j
      double precision, allocatable :: bssbij(:), bssijmax(:), awt(:)
c     row, column and order vector for reordering bssbij
      integer, allocatable :: bloci(:), blocj(:), loc(:)

c     calculate number of blocks (nb) and block boundaries (vector bb)
      rn = dfloat(n)
      if (n .ge. 50) then
         nb = nint(sqrt(dfloat(n)))
      else
         nb = 1
      endif

c     the number of paiwise block comparison
      nb2 = nb*(nb+1)/2
c     allocate memory
      allocate(bpsmax(nb), bpsmin(nb))
      allocate(bb(nb), ibmin(nb), ibmax(nb))
      allocate(bssbij(nb2), bssijmax(nb2), awt(nb2))
      allocate(bloci(nb2), blocj(nb2), loc(nb2))

c     block boundaries
      do 110 i = 1, nb
         bb(i) = nint(rn*(dfloat(i)/dfloat(nb)))
 110  continue

c     find the max, min of partial sums and their locations within blocks
      ilo = 1
      psum = 0
      psmin0 = 0
      psmax0 = 0
      ipsmin0 = n
      ipsmax0 = n
      do 20 j = 1, nb
         sx(ilo) = psum + x(ilo)*wts(ilo)
         psmin = sx(ilo)
         ipsmin = ilo
         psmax = sx(ilo)
         ipsmax = ilo
         do 10 i = ilo+1, bb(j)
            sx(i) = sx(i-1) + x(i)*wts(i)
            if (sx(i) .lt. psmin) then 
               psmin = sx(i)
               ipsmin = i
            endif
            if (sx(i) .gt. psmax) then 
               psmax = sx(i)
               ipsmax = i
            endif
 10      continue
c     store the block min, max and locations
         ibmin(j) = ipsmin
         ibmax(j) = ipsmax
         bpsmin(j) = psmin
         bpsmax(j) = psmax
c     adjust global min, max and locations
         if (psmin .lt. psmin0) then
            psmin0 = psmin
            ipsmin0 = ipsmin
         endif
         if (psmax .gt. psmax0) then
            psmax0 = psmax
            ipsmax0 = ipsmax
         endif
c     reset ilo to be the block boundary + 1
         psum = sx(bb(j))
         ilo = bb(j) + 1
 20   continue

c     calculate bss for max s_i - min s_i
      psdiff = psmax0 - psmin0
c     if the segment is all constant then psdiff = 0 and so bssmax = 0
      if (psdiff .le. 0) then
         bssmax = 0
         go to 120
      endif
      psrn = cwts(n)
      psrj = abs(cwts(ipsmax0) - cwts(ipsmin0))
      psrnj = psrj*(psrn-psrj)
      bssmax = (psdiff**2)/psrnj
      tmaxi = min(ipsmax0, ipsmin0)
      tmaxj = max(ipsmax0, ipsmin0)

c     for a pair of blocks (i,j) calculate the max absolute t-statistic
c     at the (min_i, max_j) and (max_i, min_j) locations 
c     for other indices the t-statistic can be bounded using this
c
c     if a block doesn't have the potential to exceed bssmax ignore it
c     calculate the bsslim for each block and include ones >= bssmax

      psrnov2 = psrn/2
      l = 0
      nal0 = n - al0
      do 40 i = 1, nb
         do 30 j = i, nb
c     calculate bsslim
            if (i .eq. 1) then
               ilo = 1
            else
               ilo = bb(i-1) + 1
            endif
            ihi = bb(i)
            if (j .eq. 1) then
               jlo = 1
            else
               jlo = bb(j-1) + 1
            endif
            jhi = bb(j)
c     for wCBS calculated hi and lo arc weights instead of lengths
            awthi = cwts(jhi) - cwts(ilo)
            if (jhi - ilo .gt. nal0) then
               awthi = 0
               do 35 k = 1, al0
                  awthi = max(awthi, cwts(nal0+k) - cwts(k))
 35            continue
            endif
            if (i .eq. j) then
               awtlo = cwts(ilo+al0) - cwts(ilo)
               do 36 k = ilo + 1, ihi - al0
                  awtlo = min(awtlo, cwts(k+al0) - cwts(k))
 36            continue
            else if (i+1 .eq. j) then
               awtlo = cwts(jlo) - cwts(jlo-al0)
               do 37 k = jlo - al0 + 1, ihi
                  awtlo = min(awtlo, cwts(k+al0) - cwts(k))
 37            continue
            else
               awtlo = cwts(jlo) - cwts(ihi)
            endif
c     max S_k over block j - min S_k over block i
            sij1 = abs(bpsmax(j) - bpsmin(i))
c     max S_k over block i - min S_k over block j
            sij2 = abs(bpsmax(i) - bpsmin(j))
c     if i = j then sij1 and sij2 are the same
            sijmx0 = max(sij1, sij2)
            psrnj = min(awtlo*(psrn-awtlo), awthi*(psrn-awthi))
            bsslim = (sijmx0**2)/psrnj
c     if its as large as bssmax add block
            if (bssmax .le. bsslim) then
               l = l+1
               loc(l) = l
               bloci(l) = i
               blocj(l) = j
               bssijmax(l) = bsslim
c     max sij in the (i,j) block, t-statistic etc
               if (sij1 .gt. sij2) then
                  awt(l) = abs(cwts(ibmax(j)) - cwts(ibmin(i)))
                  bssbij(l) = (sij1**2)/(awt(l)*(psrn-awt(l)))
               else
                  awt(l) = abs(cwts(ibmin(j)) - cwts(ibmax(i)))
                  bssbij(l) = (sij2**2)/(awt(l)*(psrn-awt(l)))
               endif
            endif
 30      continue
 40   continue

      nb1 = l

c     Now sort the t-statistics by their magnitude
      call qsort4(bssbij, loc, 1, nb1)

c     now go through the blocks in reverse order (largest down)
      do 100 l = nb1, 1, -1
         k = loc(l)
c     need to check a block only if it has potential to increase bss
c     rjlo is the smalllest (j-i) in the block and rjhi is the largest
         bsslim = bssijmax(k)
         if (bssmax .le. bsslim) then
c     bi, bj give the block location
            bi = bloci(k)
            bj = blocj(k)
            awtmax = awt(k)
            if (bi .eq. 1) then
               ilo = 1
            else
               ilo = bb(bi-1) + 1
            endif
            ihi = bb(bi)
            if (bj .eq. 1) then
               jlo = 1
            else
               jlo = bb(bj-1) + 1
            endif
            jhi = bb(bj)
            awthi = cwts(jhi) - cwts(ilo)
            if (bi .eq. bj) then
               awtlo = 0
            else
               awtlo = cwts(jlo) - cwts(ihi)
            endif
c
c     if arc wt is larger than half total wt (psrn/2) make is psrn - arc wt
c
            if (awtmax .gt. psrn - awtmax) awtmax = psrn - awtmax
c
c     if awtlo <= psrn/2 start from (ihi, jlo) and go up
c     if awthi >= psrn/2 start from (ilo, jhi) and go down
c
            if (awtlo .le. psrnov2) then
               if (bi .eq.bj) then 
                  ihi1 = ihi - al0
               else
                  ihi1 = ihi
               endif
               do 60 i = ihi1, ilo, -1
                  jlo1 = max(i + al0, jlo)
                  do 55 j = jlo1, jhi
                     awt1 = cwts(j) - cwts(i)
                     if (awt1 .le. awtmax) then
                        bijbss = (sx(j) - sx(i))**2/(awt1*(psrn-awt1))
                        if (bijbss .gt. bssmax) then
                           bssmax = bijbss
                           tmaxi = i
                           tmaxj = j
                        endif
                     endif
 55               continue
 60            continue
            endif
c
c     make arc wt  psrn - arc wt
c
            awtmax = psrn - awtmax
            if (awthi .ge. psrnov2) then
               do 70 i = ilo, ihi
                  if ((bi .eq. 1) .and. (bj .eq. nb)) then 
                     jhi1 = min(jhi, jhi - al0 + i)
                  else
                     jhi1 = jhi
                  endif
                  do 65 j = jhi1, jlo, -1
                     awt1 = cwts(j) - cwts(i)
                     if (awt1 .ge. awtmax) then
                        bijbss = (sx(j) - sx(i))**2/(awt1*(psrn-awt1))
                        if (bijbss .gt. bssmax) then
                           bssmax = bijbss
                           tmaxi = i
                           tmaxj = j
                        endif
                     endif
 65               continue
 70            continue
            endif
         endif
 100  continue

 120  if (tss.le.bssmax+0.0001) tss = bssmax + 1.0
      bssmax = bssmax/((tss-bssmax)/(rn-2.0))

c     deallocate memory
      deallocate(bpsmax, bpsmin, bb, ibmin, ibmax)
      deallocate(bssbij, bssijmax, bloci, blocj, loc, awt)

      ostat = bssmax
      iseg(1) = tmaxi
      iseg(2) = tmaxj

      return
      end

c     function for calculating the full max wtd t-statistic on permuted data
c     using a new approach to maximizing t-statistic
      double precision function wtmaxp(n,px,wts,sx,cwts,al0)
      integer n,al0
      double precision px(n),wts(n),sx(n),cwts(n)
c     
c     look at the partial sums in blocks of size sqrt(n)
c     
      integer ipsmin, ipsmax, ipsmin0, ipsmax0, nb, i, j, k, l, nb1,
     1     nb2, bi, bj, ilo, ihi, jlo, jhi, ihi1, jlo1, jhi1, nal0
      double precision psum, psmin, psmax, psmin0, psmax0, bssmax,
     1     bsslim, rn, sij1, sij2, sijmx0, bijbss, awtmax, psrnov2,
     2     psdiff, psrj, psrn, psrnj, awtlo, awthi, awt1, ssq, tss
c     
c     use local arrays for working within blocks
c     block partial sum max and min 
      double precision, allocatable :: bpsmax(:), bpsmin(:)
c     location of the max and min
      integer, allocatable :: bb(:), ibmin(:), ibmax(:)

c     t statistic corresponding to max for block i,j
      double precision, allocatable :: bssbij(:), bssijmax(:), awt(:)
c     row, column and order vector for reordering bssbij
      integer, allocatable :: bloci(:), blocj(:), loc(:)

c     calculate number of blocks (nb) and block boundaries (vector bb)
      rn = dfloat(n)
      if (n .ge. 50) then
         nb = nint(sqrt(dfloat(n)))
      else
         nb = 1
      endif

c     the number of paiwise block comparison
      nb2 = nb*(nb+1)/2
c     allocate memory
      allocate(bpsmax(nb), bpsmin(nb))
      allocate(bb(nb), ibmin(nb), ibmax(nb))
      allocate(bssbij(nb2), bssijmax(nb2), awt(nb2))
      allocate(bloci(nb2), blocj(nb2), loc(nb2))

c     block boundaries
      do 110 i = 1, nb
         bb(i) = nint(rn*(dfloat(i)/dfloat(nb)))
 110  continue

c     find the max, min of partial sums and their locations within blocks
      ilo = 1
      psum = 0
      psmin0 = 0
      psmax0 = 0
      ipsmin0 = n
      ipsmax0 = n
      ssq = 0
      do 20 j = 1, nb
         sx(ilo) = psum + px(ilo)*wts(ilo)
         ssq = ssq + (px(ilo)**2)*wts(ilo)
         psmin = sx(ilo)
         ipsmin = ilo
         psmax = sx(ilo)
         ipsmax = ilo
         do 10 i = ilo+1, bb(j)
            sx(i) = sx(i-1) + px(i)*wts(i)
            ssq = ssq + (px(i)**2)*wts(i)
            if (sx(i) .lt. psmin) then 
               psmin = sx(i)
               ipsmin = i
            endif
            if (sx(i) .gt. psmax) then 
               psmax = sx(i)
               ipsmax = i
            endif
 10      continue
c     store the block min, max and locations
         ibmin(j) = ipsmin
         ibmax(j) = ipsmax
         bpsmin(j) = psmin
         bpsmax(j) = psmax
c     adjust global min, max and locations
         if (psmin .lt. psmin0) then
            psmin0 = psmin
            ipsmin0 = ipsmin
         endif
         if (psmax .gt. psmax0) then
            psmax0 = psmax
            ipsmax0 = ipsmax
         endif
c     reset ilo to be the block boundary + 1
         psum = sx(bb(j))
         ilo = bb(j) + 1
 20   continue

c     calculate bss for max s_i - min s_i
      psdiff = psmax0 - psmin0
      psrn = cwts(n)
      tss = ssq - (sx(n)/psrn)**2
      psrj = abs(cwts(ipsmax0) - cwts(ipsmin0))
      psrnj = psrj*(psrn-psrj)
      bssmax = (psdiff**2)/psrnj

c     for a pair of blocks (i,j) calculate the max absolute t-statistic
c     at the (min_i, max_j) and (max_i, min_j) locations 
c     for other indices the t-statistic can be bounded using this
c
c     if a block doesn't have the potential to exceed bssmax ignore it
c     calculate the bsslim for each block and include ones >= bssmax

      psrnov2 = psrn/2
      l = 0
      nal0 = n - al0
      do 40 i = 1, nb
         do 30 j = i, nb
c     calculate bsslim
            if (i .eq. 1) then
               ilo = 1
            else
               ilo = bb(i-1) + 1
            endif
            ihi = bb(i)
            if (j .eq. 1) then
               jlo = 1
            else
               jlo = bb(j-1) + 1
            endif
            jhi = bb(j)
c     for wCBS calculated hi and lo arc weights instead of lengths
            awthi = cwts(jhi) - cwts(ilo)
            if (jhi - ilo .gt. nal0) then
               awthi = 0
               do 35 k = 1, al0
                  awthi = max(awthi, cwts(nal0+k) - cwts(k))
 35            continue
            endif
            if (i .eq. j) then
               awtlo = cwts(ilo+al0) - cwts(ilo)
               do 36 k = ilo + 1, ihi - al0
                  awtlo = min(awtlo, cwts(k+al0) - cwts(k))
 36            continue
            else if (i+1 .eq. j) then
               awtlo = cwts(jlo) - cwts(jlo-al0)
               do 37 k = jlo - al0 + 1, ihi
                  awtlo = min(awtlo, cwts(k+al0) - cwts(k))
 37            continue
            else
               awtlo = cwts(jlo) - cwts(ihi)
            endif
c     max S_k over block j - min S_k over block i
            sij1 = abs(bpsmax(j) - bpsmin(i))
c     max S_k over block i - min S_k over block j
            sij2 = abs(bpsmax(i) - bpsmin(j))
c     if i = j then sij1 and sij2 are the same
            sijmx0 = max(sij1, sij2)
            psrnj = min(awtlo*(psrn-awtlo), awthi*(psrn-awthi))
            bsslim = (sijmx0**2)/psrnj
c     if its as large as bssmax add block
            if (bssmax .le. bsslim) then
               l = l+1
               loc(l) = l
               bloci(l) = i
               blocj(l) = j
               bssijmax(l) = bsslim
c     max sij in the (i,j) block, t-statistic etc
               if (sij1 .gt. sij2) then
                  awt(l) = abs(cwts(ibmax(j)) - cwts(ibmin(i)))
                  bssbij(l) = (sij1**2)/(awt(l)*(psrn-awt(l)))
               else
                  awt(l) = abs(cwts(ibmin(j)) - cwts(ibmax(i)))
                  bssbij(l) = (sij2**2)/(awt(l)*(psrn-awt(l)))
               endif
            endif
 30      continue
 40   continue
      nb1 = l

c     Now sort the t-statistics by their magnitude
      call qsort4(bssbij, loc, 1, nb1)

c     now go through the blocks in reverse order (largest down)
      do 100 l = nb1, 1, -1
         k = loc(l)
c     need to check a block only if it has potential to increase bss
c     rjlo is the smalllest (j-i) in the block and rjhi is the largest
         bsslim = bssijmax(k)
         if (bssmax .le. bsslim) then
c     bi, bj give the block location
            bi = bloci(k)
            bj = blocj(k)
            awtmax = awt(k)
            if (bi .eq. 1) then
               ilo = 1
            else
               ilo = bb(bi-1) + 1
            endif
            ihi = bb(bi)
            if (bj .eq. 1) then
               jlo = 1
            else
               jlo = bb(bj-1) + 1
            endif
            jhi = bb(bj)
            awthi = cwts(jhi) - cwts(ilo)
            if (bi .eq. bj) then
               awtlo = 0
            else
               awtlo = cwts(jlo) - cwts(ihi)
            endif
c
c     if arc wt is larger than half total wt (psrn/2) make is psrn - arc wt
c
            if (awtmax .gt. psrn - awtmax) awtmax = psrn - awtmax
c
c     if awtlo <= psrn/2 start from (ihi, jlo) and go up
c     if awthi >= psrn/2 start from (ilo, jhi) and go down
c
            if (awtlo .le. psrnov2) then
               if (bi .eq.bj) then 
                  ihi1 = ihi - al0
               else
                  ihi1 = ihi
               endif
               do 60 i = ihi1, ilo, -1
                  jlo1 = max(i + al0, jlo)
                  do 55 j = jlo1, jhi
                     awt1 = cwts(j) - cwts(i)
                     if (awt1 .le. awtmax) then
                        bijbss = (sx(j) - sx(i))**2/(awt1*(psrn-awt1))
                        if (bijbss .gt. bssmax) bssmax = bijbss
                     endif
 55               continue
 60            continue
            endif
c
c     make arc wt  psrn - arc wt
c
            awtmax = psrn - awtmax
            if (awthi .ge. psrnov2) then
               do 70 i = ilo, ihi
                  if ((bi .eq. 1) .and. (bj .eq. nb)) then 
                     jhi1 = min(jhi, jhi - al0 + i)
                  else
                     jhi1 = jhi
                  endif
                  do 65 j = jhi1, jlo, -1
                     awt1 = cwts(j) - cwts(i)
                     if (awt1 .ge. awtmax) then
                        bijbss = (sx(j) - sx(i))**2/(awt1*(psrn-awt1))
                        if (bijbss .gt. bssmax) bssmax = bijbss
                     endif
 65               continue
 70            continue
            endif
         endif
 100  continue

      if (tss.le.bssmax+0.0001) tss = bssmax + 1.0
      wtmaxp = bssmax/((tss-bssmax)/(rn-2.0))

c     deallocate memory
      deallocate(bpsmax, bpsmin, bb, ibmin, ibmax)
      deallocate(bssbij, bssijmax, bloci, blocj, loc, awt)

      return
      end

c     function for the max (over small arcs) wtd t-statistic on permuted data
c     new code to speed up this part 4/1/2010
      double precision function hwtmaxp(n,k,px,wts,sx,cwts,mncwt,al0)
      integer n,k,al0
      double precision px(n),wts(n),sx(n),cwts(n),mncwt(k)

      integer i, j, nmj, ipj, ipnmj
      double precision rn, rj, rnj, bssmax, bssij, psmin, psmax, psdiff,
     1     bsslim, ssq, tss

c     create blocks of size k (or k+1) to span 1 thru n
c     block partial sum max and min 
      double precision, allocatable :: bpsmax(:), bpsmin(:)
c     location of the max and min
      integer, allocatable :: bb(:)
c     variables to work on block specific data
      integer nb, ilo, ihi, l
      double precision psum, psdiffsq

      rn = dfloat(n)
c     number of blocks
      nb = int(rn/dfloat(k))
c     allocate memory
      allocate(bpsmax(nb), bpsmin(nb))
      allocate(bb(nb))
c     block boundaries
      do 110 i = 1, nb
         bb(i) = nint(rn*(dfloat(i)/dfloat(nb)))
 110  continue

c     don't need global min and max
c     find the max, min of partial sums and their locations within blocks
      ilo = 1
      psum = 0
      ssq = 0.0d0
      bssmax = 0.0d0
      rn = cwts(n)
      do 20 j = 1, nb
         sx(ilo) = psum + px(ilo)*wts(ilo)
         ssq = ssq + wts(ilo)*px(ilo)**2
         psmin = sx(ilo)
         ipsmin = ilo
         psmax = sx(ilo)
         ipsmax = ilo
         do 10 i = ilo+1, bb(j)
            sx(i) = sx(i-1) + px(i)*wts(i)
            ssq = ssq + wts(i)*px(i)**2
            if (sx(i) .lt. psmin) then 
               psmin = sx(i)
               ipsmin = i
            endif
            if (sx(i) .gt. psmax) then 
               psmax = sx(i)
               ipsmax = i
            endif
 10      continue
c     store the block min, max and locations
         bpsmin(j) = psmin
         bpsmax(j) = psmax
c     reset ilo to be the block boundary + 1
         psum = sx(bb(j))
         ilo = bb(j) + 1
c     calculate the bss at the block max & min pr
         i = abs(ipsmin - ipsmax)
         if ((i .le. k) .and. (i .ge. al0)) then
            rj = abs(cwts(ipsmax) - cwts(ipsmin))
            rnj = rj*(rn-rj)
            bssij = (bpsmax(j) - bpsmin(j))**2/rnj
            if (bssmax .lt. bssij) bssmax = bssij
         endif
 20   continue
      tss = ssq - (sx(n)/rn)**2

c     check the first block
      ilo = 1
      ihi = bb(1)
      psdiff = bpsmax(1) - bpsmin(1)
      psdiffsq = psdiff**2
      do 40 j = al0,k
         rj = mncwt(j)
         bsslim = psdiffsq/(rj*(rn-rj))
         if (bsslim .lt. bssmax) go to 50
         sxmx = 0.0d0
         do 30 i = ilo,ihi-j
            ipj = i+j
            rj = cwts(ipj) - cwts(i)
            bssij = (sx(ipj) - sx(i))**2/(rj*(rn-rj))
            if (bssij .gt. bssmax) bssmax = bssij
 30      continue
 40   continue

c     now the minor arcs spanning the end (n)
 50   psdiff = max(abs(bpsmax(1)-bpsmin(nb)), abs(bpsmax(nb)-bpsmin(1)))
      psdiffsq = psdiff**2
      do 70 j = al0,k
         rj = mncwt(j)
         bsslim = psdiffsq/(rj*(rn-rj))
         if (bsslim .lt. bssmax) go to 100
         nmj = n-j
         do 60 i = 1,j
            ipnmj = i + nmj
            rj = cwts(ipnmj) - cwts(i)
            bssij = (sx(ipnmj) - sx(i))**2/(rj*(rn-rj))
            if (bssij .gt. bssmax) bssmax = bssij
 60      continue
 70   continue

c     now the other blocks
 100  do 200 l = 2,nb
         ilo = bb(l-1)+1
         ihi = bb(l)
         psdiff = bpsmax(l) - bpsmin(l)
         psdiffsq = psdiff**2
         do 140 j = al0,k
            rj = mncwt(j)
            bsslim = psdiffsq/(rj*(rn-rj))
            if (bsslim .lt. bssmax) go to 150
            sxmx = 0.0d0
            do 130 i = ilo,ihi-j
               ipj = i+j
               rj = cwts(ipj) - cwts(i)
               bssij = (sx(ipj) - sx(i))**2/(rj*(rn-rj))
               if (bssij .gt. bssmax) bssmax = bssij
 130        continue
 140     continue
 150     psdiff = max(abs(bpsmax(l)-bpsmin(l-1)), 
     1        abs(bpsmax(l-1)-bpsmin(l)))
         psdiffsq = psdiff**2
         do 170 j = al0,k
            rj = mncwt(j)
            bsslim = psdiffsq/(rj*(rn-rj))
            if (bsslim .lt. bssmax) go to 200
            do 160 i = ilo-j,ilo-1
               ipj = i+j
               rj = cwts(ipj) - cwts(i)
               bssij = (sx(ipj) - sx(i))**2/(rj*(rn-rj))
               if (bssij .gt. bssmax) bssmax = bssij
 160        continue
 170     continue
 200  continue

c      call dblepr("bss max", 7, bssmax, 1)

      if (tss .le. bssmax+0.0001d0) tss = bssmax + 1.0d0
      hwtmaxp = bssmax/((tss-bssmax)/(dfloat(n)-2.0d0))

c     deallocate memory
      deallocate(bpsmax, bpsmin, bb)

      return
      end

c     the new statistic routine doesn't compute mncwt
      subroutine getmncwt(n, cwts, k, mncwt, delta)
      integer n, k
      double precision cwts(n), mncwt(k), delta

      integer i, j, nmj
      double precision rj, rn

      rn = cwts(n)
      do 30 j = 1,k
         mncwt(j) = cwts(j)
         nmj = n-j
         do 10 i = 1,nmj
            rj = cwts(i+j) - cwts(i)
            mncwt(j) = min(mncwt(j), rj)
 10      continue
         do 20 i = 1, j
            rj = cwts(i+nmj) - cwts(i)
            mncwt(j) = min(mncwt(j), rn-rj)
 20      continue
 30   continue

      j = k+1
      nmj = n-j
      delta = cwts(j)
      do 40 i = 1,nmj
         rj = cwts(i+j) - cwts(i)
         delta = min(delta, rj)
 40   continue
      do 50 i = 1, j
         rj = cwts(i+nmj) - cwts(i)
         delta = min(delta, rn-rj)
 50   continue

      delta = delta/cwts(n)

      return
      end	  
	  
      subroutine prune(n,x,nseg,lseg,pcut,sx,ncpt,loc,loc1,pncpt)
      integer n, nseg, lseg(nseg), ncpt, loc(ncpt), loc1(2,ncpt), pncpt
      double precision x(n), pcut, sx(nseg)

      integer i, j, k, kmj
      double precision ssq, wssqk, wssq1, wssqj
      logical jleft

      double precision errssq
      external errssq

      ssq = 0.0
      do 10 i = 1,n
         ssq = ssq + x(i)**2
 10   continue
      k = 0
      do 15 i = 1,nseg
         sx(i) = 0
         do 14 j = 1,lseg(i)
            k = k + 1
            sx(i) = sx(i) + x(k)
 14      continue
 15   continue

      k = nseg - 1
      do 16 i = 1,k
         loc(i) = i
         loc1(2,i) = i
 16   continue
      wssqk = ssq - errssq(nseg,lseg,sx,k,loc)
      do 100 j = k-1, 1, -1
         kmj = k - j
         jleft = .TRUE.
         do 20 i = 1,j
            loc(i) = i
            loc1(1,i) = i
 20      continue
         wssqj = ssq - errssq(nseg,lseg,sx,j,loc)
         do 30 while(jleft) 
            call combn(j, kmj, loc, jleft)
            wssq1 = ssq - errssq(nseg,lseg,sx,j,loc)
            if (wssq1 .le. wssqj) then
               wssqj = wssq1
               do 25 i = 1,j
                  loc1(1,i) = loc(i)
 25            continue
            endif
 30      continue
         if (wssqj/wssqk .gt. 1+pcut) then
            pncpt = j+1
            do 35 i = 1,pncpt
               loc(i) = loc1(2,i)
 35         continue
            return
         else
            do 40 i = 1,j
               loc1(2,i) = loc1(1,i)
 40         continue
         endif
 100  continue
      pncpt = 0
      return
      end

      double precision function errssq(nseg,lseg,sx,k,loc)
      integer nseg, lseg(nseg),k,loc(k)
      double precision sx(nseg)

      double precision segsx
      integer segnx, i, j

      errssq = 0.0
      segsx = 0.0
      segnx = 0
      do 10 i = 1,loc(1)
         segsx = segsx + sx(i)
         segnx = segnx + lseg(i)
 10   continue
      errssq = errssq + segsx**2/dfloat(segnx)
      do 20 j = 2,k
         segsx = 0.0
         segnx = 0
         do 15 i = loc(j-1)+1,loc(j)
            segsx = segsx + sx(i)
            segnx = segnx + lseg(i)
 15      continue
         errssq = errssq + segsx**2/dfloat(segnx)
 20   continue
      segsx = 0.0
      segnx = 0
      do 25 i = loc(k)+1,nseg
         segsx = segsx + sx(i)
         segnx = segnx + lseg(i)
 25   continue
      errssq = errssq + segsx**2/dfloat(segnx)

      return
      end
c
c     This program generates Choose(n,r) combinations one at a time
c     Adapted from Algorithm AS 88  Appl. Statist. (1975) Vol.24, No. 3
c     
      subroutine combn(r, nmr, loc, rleft)
      integer r, nmr, loc(r)
      logical rleft

      integer i,j

      i = r
      do 10 while (loc(i) .eq. nmr+i)
         i = i-1
 10   continue
      loc(i) = loc(i) + 1
      do 20 j = i+1,r
         loc(j) = loc(j-1)+1
 20   continue
      if (loc(1) .eq. nmr+1) rleft = .FALSE.

      return
      end
	  
c     Ternary segmentation with permutation reference distribution
c     probes have weights due to differences in variances

      subroutine wfindcpt(n,x,tss,wts,rwts,cwts,px,sx,nperm,cpval,ncpt,
     1     icpt,hybrid,al0,hk,mncwt,delta,ngrid,sbn,sbdry,tol)
      integer n,nperm,ncpt,icpt(2),al0,hk,ngrid,sbn,sbdry(sbn)
      logical hybrid
      double precision x(n),tss,wts(n),rwts(n),cwts(n),px(n),sx(n),
     1     cpval,mncwt(hk),delta,tol

      integer np,nrej,nrejc,iseg(2),n1,n2,n12,l,k
      double precision ostat,ostat1,pstat,tpval,pval1,pval2

c     new functions to replace tmax and htmax (also tmaxo replaces tmax1)
      double precision tailp, wtmaxp, hwtmaxp, wtpermp
      external tailp, wtmaxp, hwtmaxp, wtpermp

      call rndstart()

      nrej = 0
      ncpt = 0

c     call the observed statistic routine
      call wtmaxo(n,x,wts,tss,sx,cwts,iseg,ostat,al0)
      ostat1 = sqrt(ostat)
      ostat = ostat * 0.99999

c     if maximal t-statistic is too small (for now use 0.1) don't split
      if (ostat1 .le. 0.1) go to 500
c     if maximal t-statistic is too large (for now use 7.0) split
c     also make sure it's not affected by outliers i.e. small seglength
      l = min(iseg(2) - iseg(1), n - iseg(2) + iseg(1))
      if ((ostat1 .ge. 7.0) .and. (l .ge. 10)) go to 200
c     o.w calculate p-value and decide if & how data are segmented
      if (hybrid) then
         call getmncwt(n, cwts, hk, mncwt, delta)
c     delta is a function of arc lengths
         pval1 = tailp(ostat1, delta, n, ngrid, tol)
         if (pval1 .gt. cpval) go to 500
         pval2 = cpval - pval1
         nrejc = int(pval2*dfloat(nperm))
         k=nrejc*(nrejc+1)/2 + 1
         do 50 np = 1,nperm
c     call permutation code for data with weights
            call wxperm(n,x,px,rwts)
c     call the small arc permutation statistic function
            pstat = hwtmaxp(n,hk,px,wts,sx,cwts,mncwt,al0)
            if (ostat.le.pstat) then
               nrej = nrej + 1
               k = k + 1
            endif
            if (nrej.gt.nrejc) go to 500
            if (np .ge. sbdry(k)) go to 200
 50      continue
      else
         nrejc = int(cpval*dfloat(nperm))
         k=nrejc*(nrejc+1)/2 + 1
         do 100 np = 1,nperm
c     call permutation code for data with weights
            call wxperm(n,x,px,rwts)
c     call full data permutation statistic function
            pstat = wtmaxp(n,px,wts,sx,cwts,al0)
            if (ostat.le.pstat) then
               nrej = nrej + 1
               k = k + 1
            endif
            if (nrej.gt.nrejc) go to 500
            if (np .ge. sbdry(k)) go to 200
 100     continue
      endif
 200  if (iseg(2).eq.n) then
         ncpt = 1
         icpt(1) = iseg(1)
      else
         if(iseg(1).eq.0) then
            ncpt = 1
            icpt(1) = iseg(2)
         else
            l = 1
            n1 = iseg(1)
            n12 = iseg(2)
            n2 = n12 - n1
            tpval = wtpermp(n1,n2,n12,x(l),px,wts(l),rwts(l),nperm)
            if (tpval.le.cpval) then
               ncpt = 1
               icpt(1) = iseg(1)
            endif
            l = iseg(1) + 1
            n12 = n - iseg(1)
            n2 = n - iseg(2)
            n1 = n12 - n2
            tpval = wtpermp(n1,n2,n12,x(l),px,wts(l),rwts(l),nperm)
            if (tpval.le.cpval) then
               ncpt = ncpt + 1
               icpt(ncpt) = iseg(2)
            endif
         endif
      endif

 500  call rndend()

      return
      end

c     *******  code to permute the data vector with weights  ********
c     since variance of probe i is inversely proportional to weights 
c     multiply by square root, permute and then divide by square root
      subroutine wxperm(n,x,px,rwts)
      integer n
      double precision x(n),px(n),rwts(n)

      integer i,j
      double precision cc,tmpx

      double precision dunif
      external dunif

      do 10 i = 1,n
         px(i) = x(i)*rwts(i)
 10   continue

      do 20 i = n,1,-1
         cc = dunif()
         j = int(cc*dfloat(i))+1
         tmpx = px(i)
         px(i) = px(j)/rwts(i)
         px(j) = tmpx
 20   continue

      return
      end

c     function for the p-value of t-statistics for removing edge effects
      double precision function wtpermp(n1,n2,n,x,px,wts,rwts,nperm)
      integer n1,n2,n,nperm
      double precision x(n),px(n),wts(n),rwts(n)

      integer np,i,m1,j,nrej
      double precision xsum1,xsum2,xbar,ostat,pstat,rn1,rn2,rm1,
     1     tstat, tss, rn, cc, tmpx

      double precision dunif
      external dunif

      if (n1.eq.1 .or. n2.eq.1) then
         nrej = nperm
         go to 110
      endif
      xsum1 = 0.0
      tss = 0.0
      rn1 = 0.0
      do 10 i=1,n1
         px(i) = x(i)*rwts(i)
         xsum1 = xsum1 + wts(i)*x(i)
         tss = tss + wts(i)*x(i)**2
         rn1 = rn1 + wts(i)
 10   continue
      xsum2 = 0.0
      rn2 = 0.0
      do 20 i=n1+1,n
         px(i) = x(i)
         xsum2 = xsum2 + wts(i)*x(i)
         tss = tss + wts(i)*x(i)**2
         rn2 = rn2 + wts(i)
 20   continue
      rn = rn1 + rn2
      xbar = (xsum1 + xsum2)/rn
      tss = tss - rn*(xbar**2)
      if (n1.le.n2) then
         m1 = n1
         rm1 = rn1
         ostat = 0.99999*abs(xsum1/rn1 - xbar)
         tstat = (ostat**2)*rn1*rn/rn2
      else
         m1 = n2
         rm1 = rn2
         ostat = 0.99999*abs(xsum2/rn2 - xbar)
         tstat = (ostat**2)*rn2*rn/rn1
      endif
      nrej = 0
      tstat = tstat/((tss-tstat)/(dfloat(n)-2.0))
c     if observed t is large (> 5) don't bother with permutation p-value
c     also make sure there are enough observations i.e. m1 >= 10
      if ((tstat .gt. 25) .and. (m1 .ge. 10)) go to 110
      do 100 np = 1,nperm
         xsum1 = 0
         do 30 i = n,n-m1+1,-1
            cc = dunif()
            j = int(cc*dfloat(i))+1
            tmpx = px(i)
            px(i) = px(j)
            px(j) = tmpx
c     the observation should be divided by sqrt(wts(i)) to get the correct 
c     probe variance.  But should be multiplied by wts(i) for statistic
            xsum1 = xsum1 + px(i)*rwts(i)
 30      continue
         pstat = abs(xsum1/rm1 - xbar)
         if (ostat.le.pstat) nrej = nrej + 1
 100  continue
 110  wtpermp = dfloat(nrej)/dfloat(nperm)

      return
      end
	  
c     Ternary segmentation with permutation reference distribution
      subroutine fndcpt(n,x,tss,px,sx,nperm,cpval,ncpt,icpt,ibin,
     1     hybrid,al0,hk,delta,ngrid,sbn,sbdry,tol)
      integer n,nperm,ncpt,icpt(2),al0,hk,ngrid,sbn,sbdry(sbn)
      logical ibin,hybrid
      double precision x(n),tss,px(n),sx(n),cpval,delta,tol

      integer np,nrej,nrejc,iseg(2),n1,n2,n12,l,k
      double precision ostat,ostat1,pstat,tpval,pval1,pval2

c     new functions to replace tmax and htmax (also tmaxo replaces tmax1)
      double precision tailp, tmaxp, htmaxp, tpermp
      external tailp, tmaxp, htmaxp, tpermp

      call rndstart()

      nrej = 0
      ncpt = 0

c      call tmax1(n,twon,x,tss,sx,tx,iseg,ostat,ibin)
      call tmaxo(n,x,tss,sx,iseg,ostat,al0,ibin)
      ostat1 = sqrt(ostat)
      ostat = ostat * 0.99999
c      call dblepr("Max Stat",8,ostat,1)
c      call intpr("Location",8,iseg,2)

c     if maximal t-statistic is too small (for now use 0.1) don't split
      if (ostat1 .le. 0.1) go to 500
c     if maximal t-statistic is too large (for now use 7.0) split
c     also make sure it's not affected by outliers i.e. small seglength
      l = min(iseg(2) - iseg(1), n - iseg(2) + iseg(1))
      if ((ostat1 .ge. 7.0) .and. (l .ge. 10)) go to 200
c     o.w calculate p-value and decide if & how data are segmented
      if (hybrid) then
         pval1 = tailp(ostat1, delta, n, ngrid, tol)
         if (pval1 .gt. cpval) go to 500
         pval2 = cpval - pval1
         nrejc = int(pval2*dfloat(nperm))
         k=nrejc*(nrejc+1)/2 + 1
         do 50 np = 1,nperm
            call xperm(n,x,px)
c            pstat = htmax(n,twon,hk,tss,px,sx,tx,ibin)
            pstat = htmaxp(n,hk,tss,px,sx,al0,ibin)
            if (ostat.le.pstat) then
               nrej = nrej + 1
               k = k + 1
            endif
            if (nrej.gt.nrejc) go to 500
            if (np .ge. sbdry(k)) go to 200
 50      continue
      else
         nrejc = int(cpval*dfloat(nperm))
         k=nrejc*(nrejc+1)/2 + 1
         do 100 np = 1,nperm
            call xperm(n,x,px)
c            pstat = tmax(n,twon,tss,px,sx,tx,ibin)
            pstat = tmaxp(n,tss,px,sx,al0,ibin)
c     call dblepr("Perm Max Stat",13,pstat,1)
            if (ostat.le.pstat) then
               nrej = nrej + 1
               k = k + 1
            endif
c     call intpr("num rej",7,nrej,1)
            if (nrej.gt.nrejc) go to 500
            if (np .ge. sbdry(k)) go to 200
 100     continue
      endif
 200  if (iseg(2).eq.n) then
         ncpt = 1
         icpt(1) = iseg(1)
      else
         if(iseg(1).eq.0) then
            ncpt = 1
            icpt(1) = iseg(2)
         else
            l = 1
            n1 = iseg(1)
            n12 = iseg(2)
            n2 = n12 - n1
            tpval = tpermp(n1,n2,n12,x(l),px,nperm)
c            call dblepr("binseg p-value",14,tpval,1)
            if (tpval.le.cpval) then
               ncpt = 1
               icpt(1) = iseg(1)
            endif
            l = iseg(1) + 1
            n12 = n - iseg(1)
            n2 = n - iseg(2)
            n1 = n12 - n2
            tpval = tpermp(n1,n2,n12,x(l),px,nperm)
c            call dblepr("binseg p-value",14,tpval,1)
            if (tpval.le.cpval) then
               ncpt = ncpt + 1
               icpt(ncpt) = iseg(2)
            endif
         endif
      endif

 500  call rndend()

      return
      end

c     code to permute the data vector
      subroutine xperm(n,x,px)
      integer n
      double precision x(n),px(n)

      integer i,j
      double precision cc,tmpx

      double precision dunif
      external dunif

      do 10 i = 1,n
         px(i) = x(i)
 10   continue

      do 20 i = n,1,-1
         cc = dunif()
         j = int(cc*dfloat(i))+1
         tmpx = px(i)
         px(i) = px(j)
         px(j) = tmpx
 20   continue
      return
      end

c     function for the p-value of t-statistics for removing edge effects
      double precision function tpermp(n1,n2,n,x,px,nperm)
      integer n1,n2,n,nperm
      double precision x(n),px(n)

      integer np,i,m1,j,nrej
      double precision xsum1,xsum2,xbar,ostat,pstat,rn1,rn2,rm1,
     1     tstat, tss, rn, cc, tmpx

      double precision dunif
      external dunif

      rn1 = dfloat(n1)
      rn2 = dfloat(n2)
      rn = rn1 + rn2
      if (n1.eq.1 .or. n2.eq.1) then
         nrej = nperm
         go to 110
      endif
      xsum1 = 0.0
      tss = 0.0
      do 10 i=1,n1
         px(i) = x(i)
         xsum1 = xsum1 + x(i)
         tss = tss + x(i)**2
 10   continue
      xsum2 = 0.0
      do 20 i=n1+1,n
         px(i) = x(i)
         xsum2 = xsum2 + x(i)
         tss = tss + x(i)**2
 20   continue
      xbar = (xsum1 + xsum2)/rn
      tss = tss - rn*(xbar**2)
      if (n1.le.n2) then
         m1 = n1
         rm1 = rn1
         ostat = 0.99999*abs(xsum1/rn1 - xbar)
         tstat = (ostat**2)*rn1*rn/rn2
      else
         m1 = n2
         rm1 = rn2
         ostat = 0.99999*abs(xsum2/rn2 - xbar)
         tstat = (ostat**2)*rn2*rn/rn1
      endif
c      call dblepr("O-Stat",6,ostat,1)
      nrej = 0
      tstat = tstat/((tss-tstat)/(rn-2.0))
c      call dblepr("T-square",8,tstat,1)
c     if observed t is large (> 5) don't bother with permutation p-value
c     also make sure there are enough observations i.e. m1 >= 10
      if ((tstat .gt. 25) .and. (m1 .ge. 10)) go to 110
      do 100 np = 1,nperm
c*******************************************
c     the following is very inefficient
c*******************************************
c         call xperm(n,x,px)
c         xsum1 = 0.0
c         do 30 i=1,m1
c            xsum1 = xsum1 + px(i)
c 30      continue
c*******************************************
c     changed to the following: instead of
c     full permutation sample m1 w.o. repl
c******************************************* 
         xsum1 = 0
         do 30 i = n,n-m1+1,-1
            cc = dunif()
            j = int(cc*dfloat(i))+1
            tmpx = px(i)
            px(i) = px(j)
            px(j) = tmpx
            xsum1 = xsum1 + px(i)
 30      continue
         pstat = abs(xsum1/rm1 - xbar)
c         call dblepr("P-Stat",6,pstat,1)
         if (ostat.le.pstat) nrej = nrej + 1
 100  continue
 110  tpermp = dfloat(nrej)/dfloat(nperm)

      return
      end
	  
      subroutine bsegp(n, gendat, ostat, pval, ng, tol)
      integer n, ng
      double precision gendat(n), ostat, pval, tol

      double precision btmax, btailp
      external btmax, btailp

      ostat = btmax(n, gendat)
c      call dblepr("Max Stat",8,ostat,1)
      pval = btailp(ostat, n, ng, tol)
      if (pval .gt. 1) pval = 1.0d0

      return
      end

      double precision function btmax(n, x)
      integer n
      double precision x(n)

      integer i
      double precision sumxi, btmaxi, dn, di, ostat

      sumxi = x(1)
      ostat = 0.0
      dn = dfloat(n)
      di = 1.0
      do 20 i = 2,n-2
         di = di + 1.0
         sumxi = sumxi + x(i)
         btmaxi = dn*(sumxi**2)/(di*(dn-di))
         if (ostat .lt. btmaxi) then
            ostat = btmaxi
c            ibseg = i
         endif
 20   continue
      btmax = sqrt(ostat)

      return
      end

c     pseudo confidence interval based on permutations
      subroutine bsegci(n, k, sumxk, x, px, sr, vfact, nperm, bsloc)
      integer n, k, sr(2), nperm, bsloc(nperm)
      double precision sumxk, x(n), px(n), vfact(n)

      integer k1, nk, np, ibseg

      call rndstart()
      k1 = k+1
      nk = n-k
      do 10 np = 1, nperm
         call xperm(k,x,px)
         call xperm(nk,x(k1),px(k1))
         call btmxci(n,k,sr,px,vfact,ibseg,sumxk)
         bsloc(np) = ibseg
 10   continue
      call rndend()

      return
      end

      subroutine btmxci(n,k,sr,x,vfact,ibseg,sumxk)
      integer n,k,sr(2),ibseg
      double precision x(n),vfact(n),sumxk

      integer i
      double precision sumxi, ostat, btmaxi

      ostat = vfact(k)*(sumxk**2)
      ibseg = k
      sumxi = sumxk
      do 10 i = k-1,sr(1),-1
         sumxi = sumxi - x(i+1)
         btmaxi = vfact(i)*(sumxi**2)
         if (ostat .lt. btmaxi) then
            ostat = btmaxi
            ibseg = i
         endif
 10   continue

      sumxi = sumxk
      do 20 i = k+1,sr(2)
         sumxi = sumxi + x(i)
         btmaxi = vfact(i)*(sumxi**2)
         if (ostat .lt. btmaxi) then
            ostat = btmaxi
            ibseg = i
         endif
 20   continue
      ostat = sqrt(ostat)

      return
      end
	  
c     tail probability of circular binary segmentation statistic
c     from Siegmund (1988) or Yao (1989) paper
      double precision function tailp(b, delta, m, ngrid, tol)
      double precision b, delta, tol
      integer m, ngrid
c     it1tsq is the integral of 1/(t*(1-t))**2
      double precision nu, it1tsq
      external nu, it1tsq

      double precision t, tl, dincr, bsqrtm, x, nux
      integer i

      dincr = (0.5d0 - delta)/dfloat(ngrid)
      bsqrtm = b/sqrt(dfloat(m))

      tl = 0.5d0 - dincr
      t = 0.5d0 - 0.5d0*dincr
      tailp = 0.0d0
      do 10 i = 1,ngrid
         tl = tl + dincr
         t = t + dincr
         x = bsqrtm/sqrt(t*(1-t))
         nux = nu(x, tol)
         tailp = tailp + (nux**2)*it1tsq(tl, dincr)
 10   continue
      tailp = 9.973557d-2*(b**3)*exp(-b**2/2)*tailp
c     since test is two-sided need to multiply tailp by 2
      tailp = 2.0d0*tailp

      return
      end

c     integral of 1/(t*(1-t))**2 from x to x+a
      double precision function it1tsq(x, a)
      double precision x, a

      double precision y

      y = x + a - 0.5d0
      it1tsq = (8.0d0*y)/(1.0d0 - 4.0d0*y**2) + 
     1     2.0d0*log((1.0d0 + 2.0d0*y)/(1.0d0 - 2.0d0*y))
      y = x - 0.5d0
      it1tsq = it1tsq - (8.0d0*y)/(1.0d0 - 4.0d0*y**2) -
     1     2.0d0*log((1.0d0 + 2.0d0*y)/(1.0d0 - 2.0d0*y))

      return
      end

      double precision function nu(x, tol)
      double precision x, tol

      double precision fpnorm
      external fpnorm

      double precision lnu0, lnu1, dk, xk
      integer i, k

      if (x .gt. 0.01d0) then
         lnu1 = log(2.0d0) - 2*log(x)
         lnu0 = lnu1
         k = 2
         dk = 0
         do 10 i = 1, k
            dk = dk + 1
            xk = -x*sqrt(dk)/2.0d0
            lnu1 = lnu1 - 2.0d0*fpnorm(xk)/dk
 10      continue

         do 50 while (dabs((lnu1-lnu0)/lnu1) .gt. tol)
            lnu0 = lnu1
            do 20 i = 1,k
               dk = dk + 1
               xk = -x*sqrt(dk)/2.0d0
               lnu1 = lnu1 - 2.0d0*fpnorm(xk)/dk
 20         continue
            k = 2*k
 50      enddo
      else
         lnu1 = -0.583d0*x
      endif
      nu = exp(lnu1)

      return
      end

c     tail probability of binary segmentation statistic
c     from page 387 of Siegmund (1986) paper
      double precision function btailp(b, m, ng, tol)
      integer m, ng
      double precision b, tol

      double precision ll, ul, dincr, nulo, nuhi, x, x1, dm
      integer i, k

      double precision fpnorm, nu
      external fpnorm, nu

      dm = dfloat(m)
      k = 2
      ll = b*sqrt(1.0/dfloat(m-k) - 1.0/dfloat(m))
      ul = b*sqrt(1.0/dfloat(k) - 1.0/dfloat(m))
      dincr = (ul - ll)/dfloat(ng)

      btailp = 0.0
      x = ll
      x1 = x + (b**2)/(dm*x)
      nulo = nu(x1, tol)/x
      do 10 i = 1, ng
         x = x + dincr
         x1 = x + (b**2)/(dm*x)
         nuhi = nu(x1, tol)/x
         btailp = btailp + (nuhi + nulo)*dincr
         nulo = nuhi
 10   continue
      btailp = b*exp(-b**2/2)*btailp/2.506628275

      btailp =  btailp + 2*(1.0-fpnorm(b))

      return
      end
