FUNCTION Norm2, qin

   return, qin[3]*qin[3] + total(qin[0:2]^2);

END

FUNCTION qtvrotation, vin, qin, invert=invert 

  real = qin[3]
  vector = qin[0:2]

  if keyword_set(invert) then $  ; point rotation = qvq*
  vout = (real^2 - total(vector^2))*vin $
         + total(2*(transpose(vector)#vin))*vector $
         + 2*real*CROSSP(vector,vin) $
  else $
  vout = (2.*real^2 - 1.0)*vin $ ; frame rotation  = q*vq
         + total(2*(transpose(vin)#vector))*vector $
         + 2*real*CROSSP(vin,vector)

  return, vout/norm2(qin)

END
