# Randomly move and rotate a polygon
## You need to supply a polygon in sf format (poly)
## You can control the maximum extent to which it movesin both the x and y directions using x and y

## The function works as so:
### Define rot function which creates a rotation matrix
### Copy polygon geometry
### Find centroid of polygon
### Move centroid a random amount between -x and +x in the x direction, and then for y direction
### Rotate polygon by a random amount

## The output is also an sf object
## I wrote this quickly so hopefully it works...

random_move <- function(poly, x, y){
  rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
  
  new.poly <- sf::st_geometry(poly)
  
  cntrd <- sf::st_centroid(new.poly)
  
  cntrd <- cntrd + c(runif(1, -x,x), runif(1,-y,y))
  
  ncg2 = (new.poly - cntrd) * rot(pi/runif(1,0,4)) + cntrd
  
  return(ncg2)
}