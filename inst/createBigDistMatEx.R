library(SSN2)



# NOT RUN
# mf04 <- importSSN(system.file("lsndata/MiddleFork04.ssn",
#	        package = "SSN"), o.write = TRUE)
# use SpatialStreamNetwork object mf04 that was already created
data(mf04p)
# for examples, copy MiddleFork04.ssn directory to R's temporary directory
copy_lsn_to_temp()

mf04p <- ssn_import("c:/temp/test.ssn", predpts = "pred1km", overwrite = TRUE)

mf04p <- ssn_update_path(mf04p, paste0(tempdir(),'/MiddleFork04.ssn'))
##createBigDistMat(mf04p, predpts = "pred1km", o.write = TRUE, no.cores = 1)

ssn_create_distmat(mf04p, predpts = "pred1km", among_predpts = TRUE, overwrite = TRUE)

ssn_create_bigdist(mf04p, predpts = "pred1km", overwrite = TRUE,
                   among_predpts = TRUE, no_cores = 2)

# NOT RUN include prediction to prediction site distances
      # createBigDistMat(mf04p, predpts = "pred1km", o.write = TRUE,
      #     amongpreds = TRUE, no.cores = 4)

###############################

library(filematrix)
net1a.fm <- fm.open("c:/temp/test.ssn/distance.fm/pred1km/dist.net1.a.bmat")
net1a <-t(net1a.fm[,])
dim(net1a)
net1a
close(net1a.fm)

sum(net1a != net1a.old)
max(net1a - net1a.old)

net1b.fm <- fm.open("c:/temp/test.ssn/distance.fm/pred1km/dist.net1.bmat")
net1b <-(net1b.fm[,])
dim(net1b)
net1b
close(net1b.fm)

###########
old.obs <- ssn_get_stream_distmat(mf04p, "obs")
names(old.obs)
net1.old <- old.obs[[1]]
dim(net1.old)

library(filematrix)
net1.fm <- fm.open("c:/temp/test.ssn/distance.fm/obs/dist.net1.bmat")
net1 <-(net1.fm[,])
dim(net1)
net1
max(net1 - net1.old)
close(net1.fm)

##
old.preds <- ssn_get_stream_distmat(mf04p, "pred1km")

net1.old <- old.preds[[3]]
dim(net1.old)

library(filematrix)
net1.fm <- fm.open("c:/temp/test.ssn/distance.fm/pred1km/dist.net1.bmat")
net1 <-(net1.fm[,])
dim(net1)
net1
close(net1.fm)

max(net1 - net1.old)

## parLapply
clusterEvalQ(cl, {
  library(SSN2)
  library(filematrix)
})
parLapply(cl, obs.pids, function(x) amongSitesBigDistMat(ssn = ssn, pids = x,
                                                         bin.table = bin.table,
                                                         name = "obs", workspace.name = workspace.name1))


