#!MC 1400

# normalize density gradient into tangent vector {t}
$!AlterData 
  Equation = '{Dens Grad Mag} = max(1e-10,sqrt( {X Density Gradient}**2 + {Y Density Gradient}**2 + {Z Density Gradient}**2 ))'
$!AlterData 
  Equation = '{t1} = {X Density Gradient} / {Dens Grad Mag}'
$!AlterData 
  Equation = '{t2} = {Y Density Gradient} / {Dens Grad Mag}'
$!AlterData 
  Equation = '{t3} = {z Density Gradient} / {Dens Grad Mag}'

# calculate normal vector {n}, this won't be normal to {t} but will be in the correct plane
$!AlterData 
  Equation = '{n1} = ddx({t1}) * {t1} + ddy({t1}) * {t2} + ddz({t1}) * {t3}'
$!AlterData 
  Equation = '{n2} = ddx({t2}) * {t1} + ddy({t2}) * {t2} + ddz({t2}) * {t3}'
$!AlterData 
  Equation = '{n3} = ddx({t3}) * {t1} + ddy({t3}) * {t2} + ddz({t3}) * {t3}'

# and normalize it
$!AlterData 
  Equation = '{n_tot} = max( sqrt( {n1}**2 + {n2}**2 + {n3}**2 ), 1.0e-10)'
$!AlterData 
  Equation = '{n1} = {n1} / {n_tot}'
$!AlterData 
  Equation = '{n2} = {n2} / {n_tot}'
$!AlterData 
  Equation = '{n3} = {n3} / {n_tot}'

# calculate binormal vector {b}
$!AlterData 
  Equation = '{b1} = {t2} * {n3} - {t3} * {n2}'
$!AlterData 
  Equation = '{b2} = {t3} * {n1} - {t1} * {n3}'
$!AlterData 
  Equation = '{b3} = {t1} * {n2} - {t2} * {n1}'
$!AlterData 
  Equation = '{b_tot} = max( sqrt ( {b1}**2 + {b2}**2 + {b3}**2 ), 1.0e-6)'

# and normalize it
$!AlterData 
  Equation = '{b1} = {b1} / {b_tot}'
$!AlterData 
  Equation = '{b2} = {b2} / {b_tot}'
$!AlterData 
  Equation = '{b3} = {b3} / {b_tot}'

# recalculate the normal vector {n} now correctly normal to both {b} and {t}
$!AlterData 
  Equation = '{n1} = {b2} * {t3} - {b3} * {t2}'
$!AlterData 
  Equation = '{n2} = {b3} * {t1} - {b1} * {t3}'
$!AlterData 
  Equation = '{n3} = {b1} * {t2} - {b2} * {t1}'

