#' Length distributions from Hafvog data
#'
#' Extract length distributions from Hafvog data and prepare for stock estimate.
#'
#' @param x data frame with biological measurements (kvarnir in Hafvog)
#' @param stfirst first station number when excluding earlier stations. Default is empty. stfirst overwrites st.
#' @param st vector of station numbers to define stations used.
#' @param teg species number as in Hafro database.
#' @param out output directory.
#' @param save logical, should results be exported.
#' @param lfile character, name for exported lengths csv table. Default: "lengdir.csv".
#' @param lsfile character, name for exorted length distribution csv table. Default: "lendist.sum.csv".
#' @param psfile character, name of image file. Default: "lendist.png".
#' @param lel
#' @param rx
#' @param ytop
#' @param mfrow
#' @author Birkir Bardarson (birkir.bardarson@hafogvatn.is)
#' @return Length distributions
#' @examples
#' ld <- ldist(syni,teg=31,out="syni.out/",mfrow=c(3,5))
#' ldist1 <- ldist(samp1,teg=31,out=bout)
#' @

#' @export
"ldist"<-
function(x, stfirst, st, teg, out = "I:\\tmp\\out\\", save = T, lfile =
	"lengdir.csv", ldfile = "lendist.csv", lsfile = "lendist.sum.csv",psfile = "lendist.png", owrite = F,
	lel, rx, ytop,mfrow)
{
	###################################################
  # Get data from hafvog e.g. length, weight, age, maturity.
	###############################################
  if(!file.exists(out)) dir.create(out,recursive=T)
	# Hafvog data imported from a data frame.
	ldata <- x
	# Stations defined.
	if(missing(st)) {
		st <- unique(ldata$STOD[ldata$TEGUND == teg])
	}
	# Define a selection of stations above given station id (optional).
	if(!missing(stfirst)) st <- unique(ldata$STOD[ldata$STOD >= stfirst &
			ldata$TEGUND == teg])
	# Create lengths data frame, first initial station.
	lengdir <- ldata[ldata$STOD == st[1] & ldata$TEGUND == teg,  ]
	# Then the rest (if exist).
	if(length(st) > 1) {
		sl <- c(2:length(st))
		for(i in sl) {
			lengdir2 <- ldata[ldata$STOD == st[i] & ldata$TEGUND ==
				teg,  ]
			lengdir <- rbind(lengdir, lengdir2)
		}
	}
	if(save == T) {
		write.table(lengdir,file=paste(out,lfile,sep=""),sep=",",row.names=F)
	}
	#####################################################################
	# Summary of measurements. Length groups within each station.
	# Create initial (empty data frame).
	lendist <- data.frame(length = 0, station = 0, count = 0, weight = 0,
		count.i = 0, weight.i = 0, count.m = 0, weight.m = 0, part = 0, part.i = 0, part.m = 0)

	# Plot length distribution: (Should perhaps remove this from the function)
	png(file=paste(out,psfile,sep=""))
	if(!missing(mfrow)) par(mfrow = mfrow)
	for(i in st) {
		# Sample date for given species and station.
		lengdir1 <- ldata[ldata$STOD == i & ldata$TEGUND == teg,  ]

    # Create maturity logic column. Mature = 1, Immature = 0.
    lengdir1$kynth <- ifelse(lengdir1$KYNTHROSKI>2,1,0)
    lengdir.sum <- apply.shrink.dataframe(lengdir1, name.x = c(
  		"LENGD","OSLAEGT"), name.ind = c("LENGD", "STOD"), FUNS = c(length,mean))

    # Maturity.
		if(nrow(lengdir1[lengdir1$kynth==0,])>0){
      lengdir.sum.i <- apply.shrink.dataframe(lengdir1[lengdir1$kynth==0,], name.x = c(
    	"LENGD","OSLAEGT"), name.ind = c("LENGD"), FUNS = c(length,mean))
      distjoin <- join(lengdir.sum,lengdir.sum.i,"LENGD")
		}else
		  {
      distjoin <- lengdir.sum
		  distjoin$LENGD.length.1 <- 0
		  distjoin$OSLAEGT.mean.1 <- 0
		  }

    if(nrow(lengdir1[lengdir1$kynth==1,])>0){
       lengdir.sum.m <- apply.shrink.dataframe(lengdir1[lengdir1$kynth==1,], name.x = c(
    	 "LENGD","OSLAEGT"), name.ind = c("LENGD"), FUNS = c(length,mean))
       distjoin <- join(distjoin,lengdir.sum.m,"LENGD")
    }else
       {
        distjoin$LENGD.length.2 <- 0
        distjoin$OSLAEGT.mean.2 <- 0
       }

    names(distjoin) <- c("length","station","count","weight","count.i","weight.i","count.m","weight.m")
		distjoin$part <- distjoin$count/sum(distjoin$count)
    distjoin$part.i <- distjoin$count.i/sum(distjoin$count)
    distjoin$part.m <- distjoin$count.m/sum(distjoin$count)
		distjoin <- as.data.frame(distjoin)
		# Combine results in one table (lendist).
		lendist <- rbind(lendist, distjoin)

		##################################################################
    # TS ~ length relationship for chosen species.

		# Capelin:
		if(teg == 31) {
			if(missing(lel))
				le1 <- (0:40)/2
		}

		#Herring:
		if(teg == 30) {
			if(missing(lel))
				le1 <- (0:45)
		}

		# Pearlside:
		if(teg == 130) {
			if(missing(lel))
				le1 <- (0:80)/10
		}

		###############################################
		# Length ranges for chosen species.

		# Capelin:
		if(teg == 31) {
			if(missing(rx))
				rx <- c(0, 20)
		}

		#Herring:
		if(teg == 30) {
			if(missing(rx))
				rx <- c(0, 45)
		}

		# Pearlside:
		if(teg == 130) {
			if(missing(rx))
				rx <- c(0, 8)
		}

		f1 <- rep(0, length(le1))
    f1.i <- rep(0, length(le1))
    f1.m <- rep(0, length(le1))
		len <- data.frame(le1, f1, f1.i, f1.m)

		## Match empty lengths with measured lengths.
		len$m1 <- match(len$le1, distjoin$length)
		len$f1[!is.na(len$m1)] <- distjoin$count
    len$m1.i <- match(len$le1, distjoin$length[!is.na(distjoin$count.i)])
    len$f1.i[!is.na(len$m1.i)] <- distjoin$count.i[!is.na(distjoin$count.i)]
    len$m1.m <- match(len$le1, distjoin$length[!is.na(distjoin$count.m)])
    len$f1.m[!is.na(len$m1.m)] <- distjoin$count.m[!is.na(distjoin$count.m)]

		## Plot lengths as proportion of numbers.
		t <- sum(distjoin$count)
		m <- sum(distjoin$count * distjoin$length)/t
    t.i <- sum(distjoin$count.i,na.rm=T)
    m.i <- sum(distjoin$count.i * distjoin$length,na.rm=T)/t
    t.m <- sum(distjoin$count.m,na.rm=T)
    m.m <- sum(distjoin$count.m * distjoin$length,na.rm=T)/t

		if(missing(ytop)) ytop <- max((len$f1/t) * 100) * 1.1
		plot(len$le1, (len$f1/t) * 100, type = "l", xlab = "Length(cm)",
			ylab = "% by numbers", xlim = rx, ylim = c(0, ytop),lwd=4)
    lines(len$le1, (len$f1.i/t)*100,col="blue",lwd=2)
    lines(len$le1, (len$f1.m/t)*100,col="red",lwd=2)

		texti1 <- paste("St ", i)
		mtext(texti1, side = 3, line = -1.5, adj = 0.05, col = 2)
		texti2 <- paste("n = ", t)
		mtext(texti2, side = 3, line = -3, adj = 0.05, col = 2)
		texti3 <- paste("mean = ", round(m, dig = 2))
		mtext(texti3, side = 3, line = -4.4, adj = 0.05, col = 8)
	}
	dev.off()
	#######################################################
	#### Combine in one table.
	lendist$countweight <- lendist$count * lendist$weight
  lendist$countweight.i <- lendist$count.i * lendist$weight.i
  lendist$countweight.m <- lendist$count.m * lendist$weight.m
  lendist[is.na(lendist)] <- 0
	lendist.sum <- apply.shrink.dataframe(lendist, name.x = c("count",
		"countweight","count.i","countweight.i","count.m","countweight.m"), name.ind = c("length"), FUNS = sum)

  lendist.sum$weight <- lendist.sum$countweight.sum/lendist.sum$count.sum
  lendist.sum$weight.i <- lendist.sum$countweight.i.sum/lendist.sum$count.i.sum
  lendist.sum$weight.m <- lendist.sum$countweight.m.sum/lendist.sum$count.m.sum
	lendist.sum$part <- lendist.sum$count.sum/sum(lendist.sum$count.sum)
  lendist.sum$part.i <- lendist.sum$count.i.sum/sum(lendist.sum$count.sum)
  lendist.sum$part.m <- lendist.sum$count.m.sum/sum(lendist.sum$count.sum)
	lendist.sum <- lendist.sum[lendist.sum$length != 0.0,]

	names(lendist.sum)[names(lendist.sum) == "count.sum"] <- "count"
  names(lendist.sum)[names(lendist.sum) == "count.i.sum"] <- "count.i"
  names(lendist.sum)[names(lendist.sum) == "count.m.sum"] <- "count.m"

	if(save == T) {
		write.table(lendist,file=paste(out,ldfile,sep=""),sep=",",row.names=F)
	}
	if(save == T) {
		write.table(lendist.sum,file=paste(out,lsfile,sep=""),sep=",",row.names=F)
	}
return(lendist.sum)
}
