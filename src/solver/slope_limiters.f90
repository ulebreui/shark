module slope_limiter

contains

function slope_limit(slope_lft,slope_rgt)
  use parameters
  use commons
  use units
  implicit none
  real(dp) :: slope_lft,slope_rgt,slop,dcen,dsgn,slope_limit,dlim

  if(slope_type==1) then

	slope_limit=slope_lft
    if(abs(slope_rgt)<abs(slope_lft)) slope_limit = slope_rgt
    if(slope_rgt*slope_lft<0.0d0) slope_limit     = 0.0d0

   else if(slope_type==2) then

	slope_limit=2.0d0*slope_lft*slope_rgt/(slope_lft+slope_rgt)
  if(slope_rgt*slope_lft<=0.0d0) slope_limit = 0.0d0

  else if(slope_type==3) then

	dcen = half*(slope_lft+slope_rgt)
    dsgn = sign(1.0d0,dcen)
    slop = min(1.5d0*abs(slope_lft),1.5d0*abs(slope_rgt))
    dlim = slop
    if((slope_lft*slope_rgt)<=0.0d0)dlim=0.d0
    slope_limit  = dsgn*min(dlim,abs(dcen))

   else 
	print *, 'unknown slope limiter'
	stop
endif

end	function slope_limit


end	module slope_limiter