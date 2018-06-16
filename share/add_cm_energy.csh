#below threshold P1
foreach x ( `seq  55115 55155`)
  if ( -f $x.cfg ) then
    echo TauTau.CENTER_MASS_ENERGY=3.539482 \; >> $x.cfg
  endif
end
#tune point below threshold
#Points.push_back({"Point1", 55157,55161,3550.872,0.182,0,0});
foreach x ( `seq  55157 55161`)
  if ( -f $x.cfg ) then
    echo TauTau.CENTER_MASS_ENERGY=3.550872 \; >> $x.cfg
  endif
end
#first threshold point P2
#Points.push_back({"Point2", 55162,55199,3552.849,0.093,0,0});
foreach x ( `seq  55162 55199`)
  if ( -f $x.cfg ) then
    echo TauTau.CENTER_MASS_ENERGY=3.552849 \; >> $x.cfg
  endif
end

#Points.push_back({"Point3", 55200,55231,3553.934,0.08,0,0});
foreach x ( `seq  55200 55231`)
  if ( -f $x.cfg ) then
    echo TauTau.CENTER_MASS_ENERGY=3.553934 \; >> $x.cfg
  endif
end

#Points.push_back({"Point4", 55232,55239,3560.356,0.157,0,0});
foreach x ( `seq  55232 55239`)
  if ( -f $x.cfg ) then
    echo TauTau.CENTER_MASS_ENERGY=3.560356 \; >> $x.cfg
  endif
end

#Points.push_back({"Point5", 55240,55257, 3599.572,0.117,0,0});
foreach x ( `seq  55240 55257`)
  if ( -f $x.cfg ) then
    echo TauTau.CENTER_MASS_ENERGY=3.599572 \; >> $x.cfg
  endif
end


