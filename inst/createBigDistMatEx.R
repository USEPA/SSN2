library(SSN2)



# NOT RUN
# mf04 <- importSSN(system.file("lsndata/MiddleFork04.ssn",
#	        package = "SSN"), o.write = TRUE)
# # use SpatialStreamNetwork object mf04 that was already created
# data(mf04p)
# # for examples, copy MiddleFork04.ssn directory to R's temporary directory
# copy_lsn_to_temp()

mf04p <- ssn_import("c:/temp/MiddleFork04.ssn",
                    predpts = c("pred1km", "CapeHorn"),
                    overwrite = TRUE)


ssn_create_distmat(mf04p, predpts = c("pred1km", "CapeHorn"),
                   among_predpts = TRUE, overwrite = TRUE)

ssn_create_bigdist(mf04p, predpts = "pred1km", overwrite = TRUE,
                   among_predpts = TRUE, no_cores = 2)

ssn_create_bigdist(mf04p,
                   predpts = c("CapeHorn"),
                   overwrite = TRUE,
                   among_predpts = TRUE,
                   no_cores = 2)

# NOT RUN include prediction to prediction site distances
      # createBigDistMat(mf04p, predpts = "pred1km", o.write = TRUE,
      #     amongpreds = TRUE, no.cores = 4)

###############################
## Check predictions
old.preds <- ssn_get_stream_distmat(mf04p, "CapeHorn")
names(old.preds)
# net1a.old <- old.preds[[1]]
# net1b.old <- old.preds[[2]]
# net1.old <- old.preds[[3]]
net2a.old <- old.preds[[1]]
net2b.old <- old.preds[[2]]
net2.old <- old.preds[[3]]
dim(net1.preds)

## Network 1
library(filematrix)
net1a.fm <- fm.open("c:/temp/MiddleFork04.ssn/distance.fm/CapeHorn/dist.net1.a.bmat")
net1a <-t(net1a.fm[,])
dim(net1a)
net1a

sum(net1a != net1a.old)
max(net1a - net1a.old)
close(net1a.fm)

net1b.fm <- fm.open("c:/temp/MiddleFork04.ssn/distance.fm/CapeHorn/dist.net1.b.bmat")
net1b <-(net1b.fm[,])
dim(net1b)
net1b
sum(net1b != net1b.old)
close(net1b.fm)

net1.fm <- fm.open("c:/temp/MiddleFork04.ssn/distance.fm/CapeHorn/dist.net1.bmat")
net1 <-(net1.fm[,])
dim(net1)
net1
sum(net1 != net1.old)
close(net1.fm)

## Network 2
net2a.fm <- fm.open("c:/temp/MiddleFork04.ssn/distance.fm/CapeHorn/dist.net2.a.bmat")
net2a <-t(net2a.fm[,])
dim(net2a)
net2a

sum(net2a != net2a.old)
max(net2a - net2a.old)
close(net2a.fm)

net2b.fm <- fm.open("c:/temp/MiddleFork04.ssn/distance.fm/CapeHorn/dist.net2.b.bmat")
net2b <-(net2b.fm[,])
dim(net2b)
net2b

sum(net2b != net2b.old)
close(net2b.fm)

net2.fm <- fm.open("c:/temp/MiddleFork04.ssn/distance.fm/CapeHorn/dist.net2.bmat")
net2 <-(net2.fm[,])
dim(net2)
net2
sum(net2 != net2.old)
close(net2.fm)


###########
old.obs <- ssn_get_stream_distmat(mf04p, "obs")
names(old.obs)
net1.old <- old.obs[[1]]
net2.old <- old.obs[[2]]
dim(net1.old)

library(filematrix)
net1.fm <- fm.open("c:/temp/MiddleFork04.ssn/distance.fm/obs/dist.net1.bmat")
net1 <-(net1.fm[,])
dim(net1)
net1
max(net1 - net1.old)
close(net1.fm)

net2.fm <- fm.open("c:/temp/MiddleFork04.ssn/distance.fm/obs/dist.net2.bmat")
net2 <-(net2.fm[,])
dim(net2)
net2
max(net2 - net2.old)
close(net2.fm)




