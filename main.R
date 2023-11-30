# .. This program simulates the serial passage of an asexual population in an environment in which only one beneficial
# .. mutation with selection coefficient 's' is available. Mutator individuals (y) are being generated from wild-type
# .. cells (x) at a rate 'mu_m'. These mutators not only exhibit a mutation rate increased by a factor of 'm', but also 
# .. change the effect of the beneficial mutation by a factor of 'sig'. The simulation ends when the beneficial mutation 
# .. is completely fixed in either background. 
# .. Reference:
# .. Couce, A.; Guelfo, J.R. and Blazquez J.(2012). 'Mutational Spectrum Drives the Rise of Mutator Bacteria'. PLoS Genet.

  rm(list=ls())

# .. start the stopwatch
  clock <- Sys.time()

# .. PARAMETERS
  # .. population sizes
  Nmax=1e8
  Btlnck=1e5

  # .. mutation rates 
  m=100
  mu_ben=1e-7
  mu_del=1e-4
  mu_let=1e-5
  mu_m=5e-6

  # .. fitness effects   
  s=0.2
  sig=1.8
  c=-0.04

  # .. genotypic matrices
  n_ben=2
  n_del=3
  rx=matrix(nrow = n_ben, ncol = n_del)
  ry=matrix(nrow = n_ben, ncol = n_del)
  rys=matrix(nrow = n_ben, ncol = n_del)
  
  for (i in 1:n_ben) {
      for (j in 1:n_del) { 
           rx[i,j]=2+(i-1)*s+(j-1)*c
	   rys[i,j]=2+(i-1)*(s*sig)+(j-1)*c
	   ry[i,j]=2+(i-1)*s+(j-1)*c
     			 }
		  }

  # .. VARIABLES
  x=matrix(nrow = n_ben, ncol = n_del)
  y=matrix(nrow = n_ben, ncol = n_del)
  ys=matrix(nrow = n_ben, ncol = n_del)

  ratiox=matrix(nrow = n_ben, ncol = n_del)
  x_t=matrix(nrow = 10000,  ncol = 1)

  ratioy=matrix(nrow = n_ben, ncol = n_del)
  y_t=matrix(nrow = 10000,  ncol = 1) 

  ratioys=matrix(nrow = n_ben, ncol = n_del)
  ys_t=matrix(nrow = 10000,  ncol = 1)

  # .. ALGORITHM
  t=0
  x[]=0
  x[1,1]=1
  y[]=0
  ys[]=0
  ratiox=x/(sum(x)+sum(y)+sum(ys))
  ratioy=y/(sum(x)+sum(y)+sum(ys))
  ratioys=ys/(sum(x)+sum(y)+sum(ys))

  source("serial_culture.R")

  # .. GRAPHIC OUTPUT

  x11(width = 3.5, height = 4, pointsize = 14, canvas = "white")
  plot(c(1:t),x_t[1:t,1],type='b',col='blue',ylim=c(0.0001,1),xlim=c(0,t),log='y',xlab="Time (days)", ylab="Frequency")
  points(c(1:t),ys_t[1:t,1],type='b',col='red')
 	
# .. stop the stopwatch
  print(Sys.time()-clock)


