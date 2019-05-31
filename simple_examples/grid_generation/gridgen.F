      PROGRAM GRID_GENERATON
      INTEGER :: i,j,k,nkx

      i=1
      j=1
      x(1)=xl(1)
      y(1)=yl(1)
c
      do 40 k=1,nkx
        if (nxl(k).eq.0) go to 20
        dxml=(xc(k)-xl(k))/nxl(k)
        nt=nxl(k)
        tn=nt
        tn=dmax1(tn,one+em6)
        dxmn1=dmin1(dxmn(k),dxml)
        cmc=(xc(k)-xl(k)-tn*dxmn1)*tn/(tn-1.0)
        if (nt.eq.1) cmc=0.0
        bmc=xc(k)-xl(k)-cmc
        do 10 l=1,nt
          i=i+1
          rln=(float(l)-tn)/tn
   10   x(i)=xc(k)+bmc*rln-cmc*rln*rln
   20   if (nxr(k).eq.0) go to 40
        dxmr=(xl(k+1)-xc(k))/nxr(k)
        nt=nxr(k)
        tn=nt
        tn=dmax1(tn,one+em6)
        dxmn1=dmin1(dxmn(k),dxmr)
        cmc=(xl(k+1)-xc(k)-tn*dxmn1)*tn/(tn-1.0)
        if (nt.eq.1) cmc=0.0
        bmc=xl(k+1)-xc(k)-cmc
        do 30 l=1,nt
          i=i+1
          rln=float(l)/tn
   30   x(i)=xc(k)+bmc*rln+cmc*rln*rln
   40 continue
c
      if (kr.ne.4) go to 50
      i=i+1
      x(i)=x(i-1)+x(2)-x(1)
   50 continue
c
      do 90 k=1,nky
        if (nyl(k).eq.0) go to 70
        dyml=(yc(k)-yl(k))/nyl(k)
        nt=nyl(k)
        tn=nt
        tn=dmax1(tn,one+em6)
        dymn1=dmin1(dymn(k),dyml)
        cmc=(yc(k)-yl(k)-tn*dymn1)*tn/(tn-1.0)
        if (nt.eq.1) cmc=0.0
        bmc=yc(k)-yl(k)-cmc
        do 60 l=1,nt
          j=j+1
          rln=(float(l)-tn)/tn
   60   y(j)=yc(k)+bmc*rln-cmc*rln*rln
   70   if (nyr(k).eq.0) go to 90
        dymr=(yl(k+1)-yc(k))/nyr(k)
        nt=nyr(k)
        tn=nt
        tn=dmax1(tn,one+em6)
        dymn1=dmin1(dymn(k),dymr)
        cmc=(yl(k+1)-yc(k)-tn*dymn1)*tn/(tn-1.0)
        if (nt.eq.1) cmc=0.0
        bmc=yl(k+1)-yc(k)-cmc
        do 80 l=1,nt
          j=j+1
          rln=float(l)/tn
   80   y(j)=yc(k)+bmc*rln+cmc*rln*rln
   90 continue
c
      if (kt.ne.4) go to 100
      j=j+1
      y(j)=y(j-1)+y(2)-y(1)
  100 continue

  END
