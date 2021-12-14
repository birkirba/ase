#' Summarises NASC data (SA) into grid of squares.
#'
#' Summarises NASC data (SA) into grid of squares. This is the main workhorse of the ase package.
#'
#' @param x          data frame with integrated NASC values from Echoview or LSSS.
#' @param cut        logical, if a part of the mapped data should be selected with geoinside().
#' @param cutfile    file containing a polygon used to define data of intrest.

#' @param cutfileout name of created cut-poligon file.
#' @param dir        output directory.
#' @param trackfile  filename for the output of the NASC data used in this analysis.
#' @param gridfile   filename for the output of the grid data used in this analysis.
#' @param latmin     rectangle hight in latitude minutes. Default is 3 min.
#' @param lonmin     rectangle width in longitude minutes. Default is 6 min.
#' @param cn         number of columns in the grid.
#' @param rn         number of rows in the grid.
#' @param latshift   shift of grid rectangles in latitudal direction.
#' @param lonshift   shift of grid rectangles in longitudal direction.
#' @param lat        scalar or vector of latitude range for the grid of boxes. Gives the min value if scalar.
#' @param lon        longitude range for the grid of boxes. Gives the min value if scalar.
#' @param empty      defines the meaning of NASC = -9999 (EchoView).
#'            Need to consider export from EchoView carefully regarding NASC = -9999.
#'            i.e. is it 0-values (empty=0), or should it be excluded (empty="skip").
#'            If -9999 is 0-value then excluded records need to be removed before running
#'            this function or it should be excluded using "cut=T".
#'            If -9999 records are excluded then analysis data (e.g. regions) must exist
#'            for all valid analysis cells.
#' @param cex.val    cex magnification of text for value (NASC) within each box in the boxgrid plot.
#' @param cex.id     cex magnification of text for the id of each box in the boxgrid plot.
#' @param new        geoplot argument. If new is FALSE the plot is added to the current plot otherwise a new plot is made.
#' @param pol        polygon that defines outer border of survey coverage. If pol exists, squares intersecting the poligon are redefined by the intersect.
#' @return
#' Calculates the average value of NASC data within each square box of the grid.
#' Plots the squares geographically and the average value within each square. Also the id
#' of each square is plotted.
#' Returns a dataframe with a row for each calculated square.
#' Can write out the cutfile. Default name is "cut.csv".
#' Can write out the trackfile. Containing the data as used in this analysis
#' Can write out the gridfile. Containing the geographic grid of the rectangles.
#'
#' @examples
#' tmp <- boxgrid(kol,latmin=0.1,lonmin=0.25,dir=getwd())

#' @export
"boxgrid"<-
    function(x, cut = F, cutfile,cutfileout="cut.csv",dir="boxgrid/", trackfile, gridfile, latmin = 3, lonmin = 6, cn, rn, latshift =
               0.5, lonshift = 0, lat, lon, empty,cex.val,cex.id,new=F,pol)
    {

	dat <- x
	#Key column names are: "lat","lon" and "NASC".
	names(dat)[names(dat) == "Lat.M"] <- "lat"
	names(dat)[names(dat) == "Lon.M"] <- "lon"
	names(dat)[names(dat) == "Lat_M"] <- "lat"
	names(dat)[names(dat) == "Lon_M"] <- "lon"
	names(dat)[names(dat) == "PRC.NASC"] <- "NASC"
	names(dat)[names(dat) == "PRC_NASC"] <- "NASC"
	names(dat)[names(dat) == "SA"] <- "NASC"

  if(missing(cex.val)) cex.val <- 0.66
  if(missing(cex.id)) cex.id <- 0.66

	if(missing(empty) & any(dat$NASC == -9999)) stop("NASC with -9999 values. Please define empty argument e.g. skip or 0\n" )

	if(!missing(empty)){
		if(empty == "skip") dat <- dat[dat$NASC != -9999,  ]
		else(dat$NASC[dat$NASC==-9999] <- empty)
	}
  op <- par("mar")
  par(mar=c(0.5,0.5,0,0))
	#####################
	#Cut the analysis area.
	if(cut == T) {
		geoplot(dat,yaxdist=0.9)
		grid <- geodefine()
		write.table(grid,file=file.path(dir,cutfileout),sep=",",row.names=F)
		dat <- geoinside(dat, grid)
	}
	if(!missing(cutfile)){
    cut1 <- read.table(cutfile,sep=",",header=TRUE)
    geoplot(dat,yaxdist=0.9)
		dat <- geoinside(dat, cut1)
	}
  if(!missing(pol)){
    pol1 <- read.csv(file=pol)
  }


	#####################
	#Column interval(lat).
	cbil <- latmin/60
	#Row interval(lon).
	rbil <- lonmin/60

  ####################
  #Start coordinates and shift.
  if(missing(lat)) {
    lat <- min(dat$lat) - (latshift * cbil)
  } else lat <- lat
  if(missing(lon)) {
    lon <- min(dat$lon) - (lonshift * rbil)
  } else lon <- lon

	#Column is added if longitudional shift is used (lonshift).
	#Number of columns (cn).
	if(lonshift != 0) ch<-1 else ch<-0
	if(missing(cn)) {
		cn <- ceiling((max(c(lon,dat$lon), na.rm = T) - min(c(lon,dat$lon), na.rm
			 = T))/rbil) + 1 + ch
	} else cn<-cn
	#Row is added if latitudional shift is used (latshift).
	#Number of rows (rn).
	if(latshift != 0) {
		rh <- 1
	} else (rh <- 0)
	if(missing(rn)) {
		rn <- ceiling((max(c(lat,dat$lat), na.rm = T) - min(c(lat,dat$lat), na.rm
			 = T))/cbil) + 1 + rh
	} else rn <- rn

	#Startpoint of the boxes.
	start <- data.frame(lat=min(lat), lon=min(lon))
	#Vectors based on length of columns and rows.
	c <- 0:(cn - 1)
	r <- 0:(rn - 1)

	#lat coordinates of the points.
	gcol <- rep(cbil, length(r))
	gc <- rep(gcol * r, rep(length(c), length(r)))
	p <- start
	k <- data.frame(lat=rep(NA,length(gc)))
	k$lat <- p$lat + gc
	#lon coordinates of the points.
	grow <- rep(rbil, length(c))
	gr <- rep(grow * c, length(r))
	k$lon <- p$lon + gr
	k$num <- 1:length(k$lat)
	#Finding endpoints to exclude during square creation.
	nu <- k$num
	f <- floor(k$num/length(c))
	h <- k$num/length(c)
	endi <- data.frame(nu, f, h)
	endi$num <- endi$nu
	p2 <- merge(k, endi)
	#Define the point number on first corner of each square.
	n <- p2$num[p2$f != p2$h & p2$num < (length(p2$num) - length(c))]
	for(i in n) {
		p2$r1[p2$num == i] <- i
		p2$r2[p2$num == i] <- i + 1
		p2$r3[p2$num == i] <- i + length(c) + 1
		p2$r4[p2$num == i] <- i + length(c)
		p2$r5[p2$num == i] <- i
	}

	reit <- data.frame(reitn = 1:length(n))
	reit$num <- n
	p2 <- merge(p2, reit, all.x = T)

	#Plot the survey track.
	geoplot(k, grid = F, col = 1, pch = " ",yaxdist=0.9,new=new)
	geopolygon(bisland, col = 16)
	geolines(bisland)
	#geopoints(dat, pch = 16, col = 4)
	#geolines(dat, col = 4)
	geosymbols(lat=dat$lat,lon=dat$lon,z=dat$NASC,perbars=0.5,lwd=2,col="red")
	#geolines(hnit)
	#Make empty dataframe.

	box1 <- data.frame(box = 0, lat = 0, lon = 0, mnasc = 0, area = 0,
	                   count = 0, lat1 = 0, lat2 = 0, lat3 = 0, lat4 = 0, lon1 = 0, lon2 =0,
	                   lon3 = 0, lon4 = 0,area1 = 0)

	#Loop for each square.
	n <- p2$num[p2$f != p2$h & p2$num < (length(p2$num) - length(c))]
	cat(paste0("Working on square: "))
	for(i in n) {
	  cat(paste0(p2$reitn[i],","))

	  box2 <- data.frame(box = 0)
		hnit <- data.frame(num = rep(0, 5))
		hnit$num <- c(p2$r1[p2$num == i], p2$r2[p2$num == i], p2$r3[
			p2$num == i], p2$r4[p2$num == i], p2$r5[p2$num == i])
		hnit$lat <- c(p2$lat[p2$num == hnit$num[1]], p2$lat[p2$num ==
			hnit$num[2]], p2$lat[p2$num == hnit$num[3]], p2$lat[
			p2$num == hnit$num[4]], p2$lat[p2$num == hnit$num[
			5]])
		hnit$lon <- c(p2$lon[p2$num == hnit$num[1]], p2$lon[p2$num ==
			hnit$num[2]], p2$lon[p2$num == hnit$num[3]], p2$lon[
			p2$num == hnit$num[4]], p2$lon[p2$num == hnit$num[
			5]])
		hnit <- hnit[,c("lat","lon")]

		#Coordinates of square corners.
		box2$lat1 <- hnit$lat[1]
		box2$lat2 <- hnit$lat[2]
		box2$lat3 <- hnit$lat[3]
		box2$lat4 <- hnit$lat[4]
		box2$lon1 <- hnit$lon[1]
		box2$lon2 <- hnit$lon[2]
		box2$lon3 <- hnit$lon[3]
		box2$lon4 <- hnit$lon[4]
		box2$area1 <- geoarea(hnit)

		#####################
		# Intersects to polygon.
		# This was changed after boxgrid9c.R
		if(!missing(pol)){
		  test1 <- nrow(geoinside(hnit,pol1))>0
		  try(hnit <- geointersect(pol1,hnit,in.or.out=0),silent=TRUE)
		  try(hnits <- geo.Split.poly(hnit))
		}

		####################################
		# Centroid of each box.

		# The centroid function:
		centroid <- function(x,y)
		{
		  N <- length(x)
		  signedArea = 0.0;
		  x0 = 0.0; # Current vertex X
		  y0 = 0.0; # Current vertex Y
		  x1 = 0.0; # Next vertex X
		  y1 = 0.0; # Next vertex Y
		  a = 0.0;  # Partial signed area
		  cx <- cy <- 0

		  # For all vertices
		  for(i in 1:(N-1)){
		    x0 = x[i]
		    y0 = y[i]
		    x1 = x[(i)%%N+1]
		    y1 = y[(i)%%N+1]
		    a = x0*y1 - x1*y0;
		    signedArea = signedArea + a;
		    cx = cx+ (x0 + x1)*a;
		    cy = cy+ (y0 + y1)*a;
		  }

		  signedArea = signedArea*0.5
		  cx = cx/(6.0*signedArea)
		  cy = cy/(6.0*signedArea)

		  return(list(cx=cx,cy=cy,area=signedArea))
		}

		##############################################
		#Coordinates in centre of square.

		hnit$reitn <- p2$reitn[p2$num == i]

		# Centroid of each box:
		cd <- centroid(y=hnit$lat[!is.na(hnit$lat)],x=hnit$lon[!is.na(hnit$lon)])
		hnit$mlat <- cd$cy
		hnit$mlon <- cd$cx

		#Crop data within the square.
		if(missing(pol)){
		  dat2 <- geoinside(dat, reg = hnit)}else{
		    if(exists("hnits") ){ #Changed after boxgrid9c.R
		      dat2 <- geoinside(data=dat, reg = hnits[[1]])
		      if(length(hnits)>1) for(i in 2:length(hnits)) {
		        dat2b <- geoinside(data=dat,reg=hnits[[i]])
		        dat2 <- rbind(dat2,dat2b)
		        rm(dat2b)
		      }}else{if(test1==TRUE) {dat2 <- geoinside(dat, reg = hnit)}else{dat2 <- dat[0,]}
		      }
		  }
		box2$mnasc <- mean(dat2$NASC)
		box2$box <- hnit$reitn[1]
		box2$lat <- hnit$mlat[1]
		box2$lon <- hnit$mlon[1]

		#Area within each polygon box:
		if(missing(pol) | !exists("hnits")){a1 <- geoarea(hnit)}else{
		  a1 <- geoarea(hnits[[1]])
		  if(length(hnits)>1) for(a in 2:length(hnits)) {
		    a1b <- geoarea(hnits[[a]])
		    a1 <- sum(a1,a1b)
		  }
		}
		if(exists("hnits")) rm(hnits)
		box2$area <- a1
		box2$count <- length(dat2$NASC)
		rm(dat2)
		box1 <- rbind(box1, box2)

		#Plot squares having nasc values and show mean nasc in square centre.
		if(!is.na(box2$mnasc)) {
			geolines(hnit)
			geotext(lat = hnit$mlat, lon = hnit$mlon, unique(box2$
				mnasc), cex = cex.val, angle = 0)
			#geotext(lat = (hnit$mlat + (0.15 * cbil)), lon = (
			#	hnit$mlon - (0.22 * rbil)), unique(box2$
			#	box), cex = cex.id, angle = 0, col ="red")
		}
        if(i==n[1]) hnit1 <- hnit else hnit1 <- rbind(hnit1,hnit)
		    rm(hnit)
	}

	#Convert square area to nmi^2.
	box1$arean <- box1$area/(1.852^2)
	#All sqares get the proportion 1.
	box1$part <- rep(1, length(box1$box))
	box1$part1 <- box1$area/box1$area1
	box1$latmin <- latmin
	box1$lonmin <- lonmin

  par(mar=op)

	########################
  #Writing text output.
  if(!file.exists(dir)) dir.create(dir,recursive=T)
	if(!missing(trackfile)) write.table(dat,file=file.path(dir,trackfile))
  if(!missing(gridfile)) write.table(hnit1,file=file.path(dir,gridfile),sep=",",row.names=F)
  #The returned output.
	return(box1)
}
