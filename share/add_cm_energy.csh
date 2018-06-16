#threshold
foreach x ( `seq  55115 55155`)
  if ( -f $x.cfg ) then
    echo TauTau.CENTER_MASS_ENERGY=3.539482 >> $x.cfg
  endif
end
