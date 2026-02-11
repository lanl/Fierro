!
!************************************************************************
!************************************************************************
!************************************************************************
!      
      subroutine vumat( &
       nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal, &
       stepTime, totalTime, dt, cmname, coordMp, charLength, &
       props, density, strainInc, relSpinInc, &
       tempOld, stretchOld, defgradOld, fieldOld, &
       stressOld, stateOld, enerInternOld, enerInelasOld, &
       tempNew, stretchNew, defgradNew, fieldNew, &
       stressNew, stateNew, enerInternNew, enerInelasNew )

!      include 'vaba_param.inc'
      dimension propsv(nprops), density(nblock), coordMp(nblock,*), &
        charLength(nblock), strainInc(nblock,ndir+nshr), &
        relSpinInc(nblock,nshr), tempOld(nblock), &
        stretchOld(nblock,ndir+nshr), &
        defgradOld(nblock,ndir+nshr+nshr), &
        fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr), &
        stateOld(nblock,nstatev), enerInternOld(nblock), &
        enerInelasOld(nblock), tempNew(nblock), &
        stretchNew(nblock,ndir+nshr), &
        defgradNew(nblock,ndir+nshr+nshr), &
        fieldNew(nblock,nfieldv), &
        stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev), &
        enerInternNew(nblock), enerInelasNew(nblock)

      character*80 cmname


!     loop over points
      do 100 km = 1,nblock      
  100 continue

      return
      end


