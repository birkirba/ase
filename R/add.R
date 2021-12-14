#' Adds interpolated values into cells in a boxgrid.
#'
#' Adds interpolated values into cells in a boxgrid. First run boxgrid() and keep the resulting plot active.
#'
#'
#' @param x     Data frame resulting from boxgrid().
#' @param nr    The number(s) of boxes that should get interpolated values.
#' @param ip    The number(s) of boxes used in the avarege of each nr box. If ip argument is missing ip gets all adjacent boxes having values.
#' @param hl    Proportion of the interpolated box(es). Default = 1.
#' @return
#' The data frame in x gets updated with interpolated values in the given(nr) boxes.
#' The boxgrid() plot also gets updated with the newly estimataed boxes and its values.
#' @
#' @examples
#' tmp <- boxgrid(kol,latmin=0.1,lonmin=0.25,dir=getwd())
#' tmp <- add(tmp,nr=c(20,21,31,32,42,43,96,105,125,145,146,148,156,157,160))
#' @

#' @export
"add"<-
   function(x,nr,ip,hl=1)
   {

box1 <- x
n <- nr #The number(s) of the box(es) that get(s) added.
if(length(hl)==1 & length(n)>1) hl <- rep(hl,length(n))
box2 <- box1
#Set the interpolated NASC of the selected boxes.
for(i in n){
   if(missing(ip)) ip1 <- c(box1$box[box1$lat1==box1$lat4[box1$box==i] & box1$lon1==box1$lon4[box1$box==i]],
     box1$box[box1$lat1==box1$lat2[box1$box==i] & box1$lon1==box1$lon2[box1$box==i]],
     box1$box[box1$lat4==box1$lat1[box1$box==i] & box1$lon4==box1$lon1[box1$box==i]],
     box1$box[box1$lat2==box1$lat1[box1$box==i] & box1$lon2==box1$lon1[box1$box==i]]) else ip1 <-ip

   geo::geolines(lat=t(box1[box1$box==i,c("lat1","lat2","lat3","lat4","lat1")]),lon=t(box1[box1$box==i,c("lon1","lon2","lon3","lon4","lon1")]),col="red")
	box2$part[box2$box==i] <- hl[match(i,n)]
   box2$mnasc[box2$box==i] <- mean(box1$mnasc[box1$box %in% ip1],na.rm=T)
   geo::geotext(lat = box2$lat[box2$box==i], lon = box2$lon[box2$box==i], box2$mnasc[box2$box==i], csi = 0.08, angle = 0, col="green")
   cw <- box2$latmin[box2$box==i]/60
	rw <- box2$lonmin[box2$box==i]/60
   geo::geotext(lat = (box2$lat4[box2$box==i] - (0.3 * cw)),
   lon = (box2$lon4[box2$box==i] + (0.22 * rw)), i, csi = 0.08, angle = 0, col ="red")
   #geotext(lat=(box2$lat[box2$box==i]-(0.3*cw)),lon=(box2$lon[box2$box==i]+(0.22*rw)),
	#as.character(unique(box2$part[box2$box==i])),csi=0.08,angle=0,col="green")
   }

return(box2)
}
