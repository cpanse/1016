
.xicMassShift <- function(f, mass = c(758.5694, 760.5851)){
  message("reading xic from file ", basename(f), " ...")
  xic_ <- rawrr::readChromatogram(f, mass = mass)
  
  message("reading index from file ", f)
  Idx <- rawrr::readIndex(f)
  Idx.Ms <- Idx[which(Idx$MSOrder=="Ms"), ]
  
  lapply(1:length(mass), function(ii){
    xic <- xic_[[ii]]
    idx <- which(xic$intensities == max(xic$intensities))
    if (length(idx) != 1) {
      warning("more than one hit for ", mass[ii], " in ", f)
      break;
    }
    ## pick APEX and 5 left and right
    idxs <- seq(idx - 5, idx + 5)
    #plot(xic$times[idxs], xic$intensities[idxs],
    #     main = mass[ii],
    #     type='b',
    #     sub=basename(f))
    
    message("merging xic and index data.frames...")
    data.frame(StartTime=xic$times[idxs]) |> merge(Idx.Ms) -> qMs
    if (nrow(qMs) == 0) {
      warning("nothing to merge ", mass[ii], " in ", f)
      break;
    }
    
    message("reading ", nrow(qMs) , " spectra from file ", f)
    
    
    ## find nearest peak
    rawrr::readSpectrum(f, scan = qMs$scan) |>
      sapply(FUN = function(x){
        protViz::findNN(mass[ii], x$centroid.mZ) -> hit;
        x$centroid.mZ[hit]
      }) -> mzDiff
    
    #plot(xic$times[idxs], mzDiff - mass[ii],
    #     sub=basename(f));
    #abline(h = 0, col='grey')
    
    ## compose output in long format
    rbind(
      data.frame(
        times = xic$times[idxs],
        value = xic$intensities[idxs],
        attribute = "intensities",
        file = basename(f),
        mass = mass[ii]
      ),
      data.frame(
        times = xic$times[idxs],
        value = (mzDiff - mass[ii]),
        attribute = "mzDiff",
        file = basename(f),
        mass = mass[ii]
      ))
  }) |>
    Reduce(f = rbind) 
}
