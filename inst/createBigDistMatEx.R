library(SSN2)


## Import data ---------------------------------------------------------
# for examples, copy MiddleFork04.ssn directory to R's temporary directory
copy_lsn_to_temp()

mf04p <- ssn_import(paste0(tempdir(), "/MiddleFork04.ssn"),
                    predpts = c("pred1km", "CapeHorn"),
                    overwrite = TRUE)

## Delete distance folder if it exists
if(dir.exists(paste0(mf04p$path, "/distance"))) {
  closeAllConnections()
  unlink(paste0(mf04p$path, "/distance"), recursive = TRUE)

}


## Calculate old distance matrices ----------------------------------
ssn_create_distmat(mf04p, predpts = c("pred1km", "CapeHorn"),
                    among_predpts = TRUE, overwrite = TRUE)

## Get old distance matrices for CapeHorn
## This will not work if you have old and new distance matrix
## files in /distance
old.preds <- ssn_get_stream_distmat(mf04p, "pred1km")
names(old.preds)

## Extract all old matrices from distance matrix list
net1a.old <- old.preds[[1]]
net1b.old <- old.preds[[2]]
net1.old <- old.preds[[3]]
net2a.old <- old.preds[[4]]
net2b.old <- old.preds[[5]]
net2.old <- old.preds[[6]]
dim(net1a.old)
class(net1a.old)

## Calculate new distance matrices ---------------------------------
## Does not currently work with multiple sets of prediction points
## Ignore the warning con$close() for now...
ssn_create_bigdist(mf04p, predpts = "pred1km", overwrite = TRUE,
                   among_predpts = TRUE, no_cores = 2)

# ssn_create_bigdist(mf04p,
#                    predpts = "CapeHorn",
#                    overwrite = TRUE,
#                    among_predpts = TRUE,
#                    no_cores = 2)

## Compare old and new distance matrices -------------------------------
## Extract new distance matrix for dist.net2.a
library(filematrix)
net2a.fm <- fm.open(paste0(tempdir(),
                           "/MiddleFork04.ssn/distance/pred1km/dist.net2.a.bmat"))
class(net2a.fm)
net2a.fm

# ## Example of how to select by pid
# ind<- rownames(net2a.fm) %in% c("1494", "1495")
# test<- net2a.fm[ind,]

## Convert filematrix to matrix. If object is an 'a' matrix,
## it must be transposed. Not necessary for other matrices.
net2a <- net2a.fm[,]
colnames(net2a) <- colnames(net2a.fm)
rownames(net2a) <- rownames(net2a.fm)
net2a <-t(net2a)

dim(net2a)
class(net2a)

net2a[1:5,1:5]

## Compare old and new distance matrices
sum(net2a != net2a.old)
sum(colnames(net2a) != colnames(net2a.old))
sum(rownames(net2a) != rownames(net2a.old))

## Close the connection
close(net2a.fm)

## More examples ------------------
## Extract net1b. No need to transpose 'b' matrices
net1b.fm <-  fm.open(paste0(tempdir(),
                            "/MiddleFork04.ssn/distance/pred1km/dist.net1.b.bmat"))
net1b <-(net1b.fm[,])
rownames(net1b)<- rownames(net1b.fm[,])
colnames(net1b)<- colnames(net1b.fm[,])

dim(net1b)

sum(net1b != net1b.old)
sum(colnames(net1b) != colnames(net1b.old))
sum(rownames(net1b) != rownames(net1b.old))
close(net1b.fm)

## Extract new net1 matrix (pred1km x pred1km).
## Do not transpose
net1.fm <- fm.open(paste0(tempdir(),
                          "/MiddleFork04.ssn/distance/pred1km/dist.net1.bmat"))
net1 <-(net1.fm[,])
rownames(net1)<- rownames(net1.fm[,])
colnames(net1)<- colnames(net1.fm[,])
dim(net1)

sum(net1 != net1.old)
sum(colnames(net1) != colnames(net1.old))
sum(rownames(net1) != rownames(net1.old))
close(net1.fm)

## Just to make sure...
closeAllConnections()







