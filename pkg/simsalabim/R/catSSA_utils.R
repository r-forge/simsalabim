.findFreq <- function(e){
### Find Dominant Frequency of an Vector
    nFFT <- length(e)
    hh <- floor(nFFT/2)-1
    f <- (0:(hh-1))/nFFT
    Y <- fft(e-mean(e))
    Pyy <- Re(Y*Conj(Y))/nFFT
    Pyy <- Pyy[1:hh]
    Pyy <- Pyy/sum(Pyy)
    return(f[Pyy==max(Pyy)])
}
