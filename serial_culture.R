# .. This routine simulates a single serial culture experiment. It is intended to be 
# .. called by the script 'main.R'

while ((sum(ratiox[2,]+ratioy[2,]) < 1) && (sum(ratioys[2,]) < 1)) {
  t=t+1
	
  # .. this loop represents one day in a flask 
  while ( sum(x)+sum(y)+sum(ys)<Nmax ) {
    x[]=rx[]*x
    y[]=ry[]*y
    ys[]=rys[]*ys
    
    for (i in 1:n_ben) {
	for (j in 1:n_del) {

	#.. allocation of mutations
	mbx<- rpois(1,lambda=x[i,j]*mu_ben)
	mdx<- rpois(1,lambda=x[i,j]*mu_del)
	mlx<- rpois(1,lambda=x[i,j]*mu_let)
	mm<- rpois(1,lambda=x[i,j]*mu_m)
	
	mby<- rpois(1,lambda=y[i,j]*mu_ben*m)
	mdy<- rpois(1,lambda=y[i,j]*mu_del*m)
	mly<- rpois(1,lambda=y[i,j]*mu_let*m)

	mbys<- rpois(1,lambda=ys[i,j]*mu_ben*m)
	mdys<- rpois(1,lambda=ys[i,j]*mu_del*m)
	mlys<- rpois(1,lambda=ys[i,j]*mu_let*m)

	#.. this if is for the boundaries of the genotypic matrix
	a=1
	b=1

	if (i==n_ben) { a=0 }
	if (j==n_del) { b=0 }

	x[i,j]=x[i,j]-mbx-mdx-mlx-mm
	x[i+a,j]=x[i+a,j]+mbx
	x[i,j+b]=x[i,j+b]+mdx
	
	y[i,j]=y[i,j]-mby-mdy-mly

	#.. adaptive mutation can only be acquired once!
	if (i<n_ben) {ys[i+a,j]=ys[i+a,j]+mby} else {y[i,j]==y[i,j]+mby}

	y[i,j+b]=y[i,j+b]+mdy
	y[i,j]=y[i,j]+mm        

	ys[i,j]=ys[i,j]-mbys-mdys-mlys
	ys[i+a,j]=ys[i+a,j]+mbys
	ys[i,j+b]=ys[i,j+b]+mdys
	
	x[x<0]=0
	y[y<0]=0
	ys[ys<0]=0
    
	}
    }
  }
        ratiox=x/(sum(x)+sum(y)+sum(ys))
        ratioy=y/(sum(x)+sum(y)+sum(ys))
	ratioys=ys/(sum(x)+sum(y)+sum(ys))
	x_t[t,1]<-(sum(ratiox[2,]))
	y_t[t,1]<-(sum(ratioy[1,]))
	ys_t[t,1]<-(sum(ratioys[2,]))

   # .. the bottleneck
    for (i in 1:n_ben) {
	for (j in 1:n_del) {

	  if ((ratiox[i,j]>0.1)) {
	    x[i,j]=rbinom (1, Btlnck, ratiox[i,j])
	  }
	  else if ((ratiox[i,j]>0)) {

	    x[i,j]=rpois (1, lambda=Btlnck*ratiox[i,j])
	  } 
	  else {
            x[i,j]=0
	  }

	  if ((ratioy[i,j]>0.1)) {
	    y[i,j]=rbinom (1, Btlnck, ratioy[i,j])
	  }
	  else if ((ratioy[i,j]>0)) {

	    y[i,j]=rpois (1, lambda=Btlnck*ratioy[i,j])
	  } 
	  else {
            y[i,j]=0
	  }

	  if ((ratioys[i,j]>0.1)) {
	    ys[i,j]=rbinom (1, Btlnck, ratioys[i,j])
	  }
	  else if ((ratioys[i,j]>0)) {

	    ys[i,j]=rpois (1, lambda=Btlnck*ratioys[i,j])
	  } 
	  else {
            ys[i,j]=0
	  }
	}
     }

}

