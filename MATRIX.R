

### SOURCE ########################
splots <- read.csv('Saunders.csv')

p.heights <- read.csv('p.heights/p.saund.heights.5.csv')

### SAUNDERS MATRIX ###############
mat <- subset(splots, select = c(X, Y, CCDC, dbh, TrueHt, CLAVG, CRAVG))
mat <- cbind(mat, p.heights[,19])
colnames(mat) <- c("X","Y","CCDC","dbh","TrueHt","CLAVG","CRAVG","P95")


cov(mat[,c("dbh","TrueHt","CLAVG","CRAVG")])
pairs(log(mat[,c("dbh","TrueHt","CLAVG","CRAVG")]),  pch = 21, bg = c("red", "green3", "blue", "cornsilk"))

### VARIOGRAM ####################
vario <- variog(coords = mat[,1:2], data = log(mat$TrueHt),uvec=(seq(0,30, length=15)))
plot(vario)

vario <- variog(coords = mat[,1:2], data = log(mat$CLAVG),uvec=(seq(0,30, length=15)))
plot(vario)

plot(mat$TrueHt, mat$P95, ylim=c(5, 40))
abline(1,1, col="red")

p95 <- mat[,'P95']
ht <- log(mat$TrueHt)
fit <- lm(ht ~ p95)
summary(fit)

vario <- variog(coords = mat[,1:2], data = resid(fit), uvec=(seq(0,30, length=15)))
plot(vario)

p95 <- mat[,'P95']
ht <- log(mat$CRAVG)
fit <- lm(ht ~ p95)
summary(fit)

vario <- variog(coords = mat[,1:2], data = resid(fit), uvec=(seq(0,30, length=20)))
plot(vario)

par(mfrow = c(1,2))
plot(ht, fitted(fit), xlim=c(2,4), ylim=c(2,4), main = 'p.heights at radius 5')
plot(ht, fitted(fit), xlim=c(2,4), ylim=c(2,4), main = 'p.heights at radius 4')
lines(0:5,0:5, col="red")

vario.mod <- variofit(vario, c(300,10), nugget = 1, cov.model = "exponential")
plot(vario.mod)
abline(vario.mod, col = 'red')


## LiDARDATA #############################
dat <- read.lidar("/home/jason/Documents/PEF-LiDAR/PEF_pts") ## takes a few minutes

## GET P95S ###########################
p.hts <- matrix(nrow = 0, ncol = 20)

for(j in 1:nrow(splots)){

    ## GET LiDAR POINTS AROUND TREE ########
    tree.1.lidar <- grab.points(dat, splots[j,1], splots[j,2], 3)

    ## SORT LiDAR POINT HEIGHTS #############
    if(length(tree.1.lidar) > 3){

        ## SUBTRACT GROUND ELEVATION ############
        tree.1.lidar[,3] <- tree.1.lidar[,3] - splots$elev[j]
        sorted.z <- sort(tree.1.lidar[,3])

        ## GET PERCENTILE HEIGHTS FOR 5% INTERVALS ##########
        incs <- seq(.05,1,by=.05)
        p.ht <- c()
        for(i in 1:length(incs)){
            p.ht[i] <- sorted.z[round(length(sorted.z)*incs[i],0)]
        }
        p.hts <- rbind(p.hts,p.ht)
    } else {p.hts <- rbind(p.hts,rep(NA, ncol(p.hts)))
        }
    print(paste(100*(j/nrow(splots)), "% complete.", sep = ""))
}

write.csv(p.hts,"p.saund.heights.3.csv", row.names = FALSE)
