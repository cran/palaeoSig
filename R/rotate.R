`rotate` <-
function(lats, longs){
 lats<-lats*pi/180
 longs<-longs*pi/180
  x = mean(cos(lats)*cos(longs))
  y = mean(cos(lats)*sin(longs))
  z = mean(sin(lats))
  r = sqrt(x^2 + y^2 + z^2)
  lat_m = asin(z/r)*180/pi
  lon_m = atan2(y,x)*180/pi
  p_lat=lat_m+90
  if(p_lat>90)p_lat<-180-p_lat
  c(p_lat,lon_m+180,180)
}

