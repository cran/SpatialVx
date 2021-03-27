

# calculate NSS value for neigberhood n from enlarged summed field 
calculate_NSS_from_enlarged_summed_fields <- function(fsum1enl, fsum2enl, n, d)
	{

	ffrac1=calculate_fractions_from_enlarged_summed_field(fsum1enl, n, d)
	ffrac2=calculate_fractions_from_enlarged_summed_field(fsum2enl, n, d)

	top=sum(abs(ffrac1-ffrac2))
	bottom=sum(ffrac1+ffrac2)
		
	if (bottom > 0)
		{return(1-top/bottom)}
	else 
		{return(0)}
	}


# calucluate neighborhood size when NSS>0.5 using bisection
calculate_n_NSS05_using_bisection_from_the_summed_fields <- function(fsum1or, fsum2or, d)
	{
	dimx=ncol(fsum1or);
	dimy=nrow(fsum1or);

	rows=nrow(fsum1or)-2*d;
	cols=ncol(fsum1or)-2*d;

	x1=1;
	x2=2*max(rows,cols) - 1;

	if (x2 < 5)
		{
		stop("Domain size needs to be at least 2 grid points")
		}

	NSS1=calculate_NSS_from_enlarged_summed_fields(fsum1or,fsum2or,x1, d);
	NSS2=calculate_NSS_from_enlarged_summed_fields(fsum1or,fsum2or,x2, d);

	#print(fsum1or)
	#print(fsum2or)
	#print(NSS2)


	# special case when NSS>0.5 at n=1
	if (NSS1>0.5)
		return(1);

	# special case when NSS does net ever reach value 0.5
	if (NSS2 <= 0.5)
		{
		stop("NSS does not ever reach value 0.5. There is something wrong.")
		}

	# use bisection
	repeat
		{ 
		
		# select new middle point
		xnew=(x1 + x2)/2
		# if xnew is even add 1
		if ( xnew%%2 == 0) {xnew=xnew+1}
		# calulcate NSS vale at xnew
		NSSnew=calculate_NSS_from_enlarged_summed_fields(fsum1or,fsum2or,xnew,d);
		# move x1 or x2 to the middle point according to NSS value at xnew
		if (NSSnew > 0.5)
			{
			x2=xnew;
			NSS2=NSSnew;
			}
		else
			{
			x1=xnew;
			NSS1=NSSnew;
			}
		
		if ( x2 - x1 <= 2 ) {break}
		}	

	return(x2);		
	}



# caluclate dNSS 
calculate_dNSS <- function(f1, f2)
	{
	# check that the dimensions of the two fields are the same
	if (ncol(f1) != ncol(f2) || nrow(f1) != nrow(f2) )
		{
		stop("The two input fields don't have the same dimensions !!!")
		}

	# check if all values are positive or zero
	if (any(f1 < 0) || any(f2 < 0) )
		{
		stop("At least one of the two fields contains negative values. Only positive or zero values are allowed in the fields.")
		}

	# check that the fields are not empty
	sumf1=sum(f1)
	sumf2=sum(f2)
	if (sumf1 == 0 || sumf2 == 0)
		{
		stop("One or both of the input fields are empty (all values are 0) !!!")
		}

	# remove bias by multipying the second field with bias
	R=sumf1/sumf2
	f2=f2*R

	# test for a special case with identical input fields. In this case dNSS=0
	if (all (f1 == f2))
		{
		dNSS=0
		}
	
	# otherwise calculate dNSS the usual way via decomposition to overlapping/nonoverlapping fields
	else
		{
		# calculate fields with overlap removed 
		foverlaping <- pmin( f1, f2 ) 
		Q=sum(foverlaping)/sum(f1)
		f1or=f1-foverlaping
		f2or=f2-foverlaping
	
		# enlarge
		d=max(ncol(f1or),nrow(f1or));
		f1orenl=enlarge_matrix(f1or,d)
		f2orenl=enlarge_matrix(f2or,d)

		# calculate summed fields from binary fields with overlap removed 
		fsum1or=calculate_summed_field(f1orenl)
		fsum2or=calculate_summed_field(f2orenl)
				
		# use bisection to determine n size when NSS>0.5
		n_NSS05=calculate_n_NSS05_using_bisection_from_the_summed_fields(fsum1or,fsum2or,d)
		
		#print(n_NSS05)
		
		# get the dNSS value
		dNSS=(1-Q)*floor(n_NSS05/2)
		}
	
	return(dNSS)
	}


