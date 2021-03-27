
# enlarge_matrix by d gridpoints in every direction 
enlarge_matrix <- function(f,d)
	{
	fenl=matrix(0, nrow(f)+2*d, ncol(f)+2*d, byrow = FALSE)
	fenl[(1+d):(nrow(f)+d),(1+d):(ncol(f)+d)] = f
	return(fenl)
	}

# calculate summed fields which are needed for the calulation of fractions
calculate_summed_field <- function(f)
	{
	fsum1<-f
	fsum2<-f
	for(k in 1:nrow(f))
		{
		fsum1[k,]<-cumsum(f[k,])
		}
	for(k in 1:ncol(f))
		{
		fsum2[,k]<-cumsum(fsum1[,k])
		}
	return(fsum2)
	}


# calculate fractions from enlarged summed field 
calculate_fractions_from_enlarged_summed_field <- function(fsumenl, n, d)
	{

	if (n %% 2 != 1)
		{
		stop("Neighbourhood size can only be an odd positive integer")
		}
		

	dr=(n-1)/2
	rows=nrow(fsumenl)-2*d;
	cols=ncol(fsumenl)-2*d;
	
	# special case if n > 2 Nx - 1
	if (n > 2*max(rows,cols) - 1) 
		{
		return(calculate_fractions_from_enlarged_summed_field(fsumenl, 2*max(rows,cols) - 1, d))
		}
	
	f_topright = fsumenl[(d+1+dr):(d+rows+dr),(d+1+dr):(d+cols+dr)]
	f_bottomleft = fsumenl[(d+1-dr-1):(d+rows-dr-1),(d+1-dr-1):(d+cols-dr-1)]
	f_topleft =  fsumenl[(d+1+dr):(d+rows+dr),(d+1-dr-1):(d+cols-dr-1)]
	f_bottomright =  fsumenl[(d+1-dr-1):(d+rows-dr-1),(d+1+dr):(d+cols+dr)]
	
	ffrac = f_topright - f_topleft -  f_bottomright + f_bottomleft
	ffrac = ffrac/(n*n)

	return(ffrac)
	}


# calculate FSS value for neigberhood n from enlarged summed field 
calculate_FSS_from_enlarged_summed_fields <- function(fsum1enl, fsum2enl, n, d)
	{

	ffrac1=calculate_fractions_from_enlarged_summed_field(fsum1enl, n, d)
	ffrac2=calculate_fractions_from_enlarged_summed_field(fsum2enl, n, d)

	top=sum((ffrac1-ffrac2)^2)
	bottom=sum(ffrac1^2+ffrac2^2)
		
	if (bottom > 0)
		{return(1-top/bottom)}
	else 
		{return(0)}
	}

# calculate FSSvector for neigberhood vector nvector from enlarged summed field 
calculate_FSSvector_from_enlarged_summed_fields <- function(fsum1enl, fsum2enl, nvector, d)
	{
	FSSvector=nvector
	for (inn in 1:length(nvector))
		FSSvector[[inn]]=calculate_FSS_from_enlarged_summed_fields(fsum1enl, fsum2enl, nvector[inn], d)
	return(FSSvector)
	}


# calculate FSSvector for neigberhood vector from binary fields 
calculate_FSSvector_from_binary_fields <- function(fbin1, fbin2, nvector)
	{
	# check that the dimensions of the two fields are the same
	if (ncol(fbin1) != ncol(fbin2) || nrow(fbin1) != nrow(fbin2) )
		{
		stop("The two input binary fields don't have the same dimensions !!!")
		}

	# check if all values in binary fields are really binary
	if (sum(fbin1==0) + sum(fbin1==1) != length(fbin1) || sum(fbin2==0) + sum(fbin2==1) != length(fbin2))
		{
		stop("At least one of the two fields contains a non-binary value. Only values 0 and 1 are allowed in the binary input fields.")
		}

	# check if nvector contains only odd positive integer values
	if ( identical(nvector,round(nvector)) == FALSE || length(which(nvector %% 2 != 1)) > 0)
		{
		stop("nvector can only contain odd positive integer values representing the neigberhood size!")
		}


	# enlarge fields - in order to avoid the use of embedded C++ code for summed fields and use native R functions
	d=max(ncol(fbin1),nrow(fbin1));
	fbin1enl=enlarge_matrix(fbin1,d)
	fsum1enl=calculate_summed_field(fbin1enl)
	fbin2enl=enlarge_matrix(fbin2,d)
	fsum2enl=calculate_summed_field(fbin2enl)

	return(calculate_FSSvector_from_enlarged_summed_fields(fsum1enl,fsum2enl,nvector,d))
	}


# calucluate neighborhood size when FSS>0.5 using bisection
calculate_n05_using_bisection_from_the_summed_fields <- function(fsum1or, fsum2or, d)
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

	FSS1=calculate_FSS_from_enlarged_summed_fields(fsum1or,fsum2or,x1, d);
	FSS2=calculate_FSS_from_enlarged_summed_fields(fsum1or,fsum2or,x2, d);

	#print(fsum1or)
	#print(fsum2or)
	#print(FSS2)


	# special case when FSS>0.5 at n=1
	if (FSS1>0.5)
		return(1);

	# special case when FSS does never reach value 0.5
	if (FSS2 <= 0.5)
		{
		stop("FSS does never reach value 0.5. There is something wrong.")
		}

	# use bisection
	repeat
		{ 
		
		# select new middle point
		xnew=(x1 + x2)/2
		# if xnew is even add 1
		if ( xnew%%2 == 0) {xnew=xnew+1}
		# calulcate FSS vale at xnew
		FSSnew=calculate_FSS_from_enlarged_summed_fields(fsum1or,fsum2or,xnew,d);
		# move x1 or x2 to the middle point according to FSS value at xnew
		if (FSSnew > 0.5)
			{
			x2=xnew;
			FSS2=FSSnew;
			}
		else
			{
			x1=xnew;
			FSS1=FSSnew;
			}
		
		if ( x2 - x1 <= 2 ) {break}
		}	

	return(x2);		
	}



# caluclate dFSS 
calculate_dFSS <- function(fbin1, fbin2)
	{
	# check that the dimensions of the two fields are the same
	if (ncol(fbin1) != ncol(fbin2) || nrow(fbin1) != nrow(fbin2) )
		{
		stop("The two input binary fields don't have the same dimensions !!!")
		}

	# check if all values in binary fields are really binary
	if (sum(fbin1==0) + sum(fbin1==1) != length(fbin1) || sum(fbin2==0) + sum(fbin2==1) != length(fbin2))
		{
		stop("At least one of the two fields contains a non-binary value. Only values 0 and 1 are allowed in the binary input fields.")
		}


	# enlarge fields - on order to avoid the use of embedded C++ code for summed fields and use native R functions
	d=max(ncol(fbin1),nrow(fbin1));
	fbin1enl=enlarge_matrix(fbin1,d);
	fbin2enl=enlarge_matrix(fbin2,d);

	# calculate summed fields from binary fields
	fsum1 = calculate_summed_field(fbin1enl)
	fsum2 = calculate_summed_field(fbin2enl)


	# check that the freqeuncy bias is not too large
	dimx=ncol(fsum1);
	dimy=nrow(fsum1);
	sum1=fsum1[dimy,dimx]
	sum2=fsum2[dimy,dimx]
	
	if (sum1 == 0 || sum2 == 0)
		{
		stop("One of the two input binary fields has no precipation (all values are 0) !!!")
		}
	
	bias=max(sum1/sum2,sum2/sum1)
	if (bias > 2) 
		{
		stop("The freqeuncy bias for the two input binary field is larger than 2. This is not allowed since it will produce unrealistic results for dFSS. We recomend using freqeuncy threshold to generate the binary input fields with no bias!")
		}
	if (bias > 1.5) 
		{
		warning("The freqeuncy bias for the two input binary field is larger than 1.5. In some cases this might result in unrealistic results for the dFSS. We recomend using freqeuncy threshold to generate the binary input fields with no bias!")
		}

	# calculate FSSn1
	FSSn1=calculate_FSS_from_enlarged_summed_fields(fsum1, fsum2, 1, d)
	
	dFSS=0
	# test for a special case with identical binary input fields when FSSn1=1 -> dFSS=0
	if (FSSn1 == 1)
		{
		dFSS=0
		}
	
	# otherwise calculate dFSS the usual way via overlap removed fields
	else
		{
		# calculate fields with overlap removed 
		funion=fbin1*fbin2;
		fbin1or=fbin1-funion;
		fbin2or=fbin2-funion;
	
		# enlarge
		fbin1orenl=enlarge_matrix(fbin1or,d);
		fbin2orenl=enlarge_matrix(fbin2or,d);

		# calculate summed fields from binary fields with overlap removed 
		fsum1or=calculate_summed_field(fbin1orenl)
		fsum2or=calculate_summed_field(fbin2orenl)

		# calculate the frequeny bias for the two fields with overlap removed 
		sum1or=fsum1or[dimy,dimx]
		sum2or=fsum2or[dimy,dimx]
		if (sum1or == 0 || sum2or == 0)
			{
			stop("One of the two fields with overlap removed has no precipation (all values are 0) !!!. We recomend using freqeuncy threshold to generate the binary input fields with no bias!")
			}

		biasor=max(sum1or/sum2or,sum2or/sum1or)
		if (biasor > 2) 
			{
			stop("The freqeuncy bias for the two fields with overlap removed is larger than 2. This is not allowed since it will produce unrealistic results for dFSS. We recomend using freqeuncy threshold to generate the binary input fields with no bias!")
			}
		if (biasor > 1.5) 
			{
			warning("The freqeuncy bias for the two fields with overlap removed is larger than 1.5. In some cases this might result in unrealistic results for the dFSS. We recomend using freqeuncy threshold to generate the binary input fields with no bias!")
			}
		
				
		# use bisection to setermine n size when FSS>0.5
		n_FSS05=calculate_n05_using_bisection_from_the_summed_fields(fsum1or,fsum2or,d)
		
		#print(n_FSS05)
		
		# use the FSS=0.5 rule to determine dFSS
		dFSS=(1-FSSn1)*floor(n_FSS05/2)
		}
	
	return(dFSS)
	}

# Calculate FSSwind score value for the neighborhood sizes specified in nvector. The matrixes findex1 and findex2 are the two input wind class index fields. The two input fields can only contain positive integer values that represent wind classes. The smallest allowed windclass index is 1 and not 0!
calculate_FSSwind <- function(findex1, findex2, nvector)
	{
	# check that the dimensions of the two input fields are the same
	if (ncol(findex1) != ncol(findex2) || nrow(findex1) != nrow(findex2) )
		{
		stop("The two input wind class index fields don't have the same dimensions !!!")
		}

	# check if all values in index fields are really integers and positive
	max_index=max(max(findex1),max(findex2))
	min_value=min(min(findex1),min(findex2))
	if ( identical(findex1,round(findex1)) == FALSE || identical(findex2,round(findex2)) == FALSE || min_value < 1 )
		{
		stop("At least one of the two index fields contains non-integer or non-positive values. Only positive integer values that represent wind classes are allowed in the input index fields! The smallest allowed windclass index is 1 and not 0.")
		}

	# check if nvector contains only odd positive integer values
	if ( identical(nvector,round(nvector)) == FALSE || length(which(nvector %% 2 != 1)) > 0)
		{
		stop("nvector can only contain odd positive integer values representing the neigberhood size!")
		}

	# decompose index fields into binary fields and calculate summed fields that are used for faster calulation of fractions
	fsum1enl_list = list()
	fsum2enl_list = list()
	d=max(ncol(findex1),nrow(findex1));
	for (iclass in 1:max_index)
		{
		# create binary field for this wind class
		ftemp=findex1
		ftemp[ ftemp != iclass ] <- 0
		ftemp[ ftemp == iclass ] <- 1
		# create enlarged summed field for this wind class
		ftemp2=enlarge_matrix(ftemp,d);
		fsum1enl_list[[iclass]] <- calculate_summed_field(ftemp2)

		# create binary field for this wind class
		ftemp=findex2
		ftemp[ ftemp != iclass ] <- 0
		ftemp[ ftemp == iclass ] <- 1
		# create summed field for this wind class
		ftemp2=enlarge_matrix(ftemp,d);
		fsum2enl_list[[iclass]] <- calculate_summed_field(ftemp2)
		}


	
	# calculate the FSSwind values for all neigberhood sizes specified in nvector
	FSSwind=nvector
	for (inn in 1:length(nvector))
		{
		numerator=0
		denominator=0

		for (iclass in 1:max_index)
			{
			# calculate the fractions for this wind class
			ffrac1 <- calculate_fractions_from_enlarged_summed_field(fsum1enl_list[[iclass]], nvector[inn],d)
			ffrac2 <- calculate_fractions_from_enlarged_summed_field(fsum2enl_list[[iclass]], nvector[inn],d)

			# add the contribution of this wind class to the sums in the numerator and denominator
			numerator = numerator + sum((ffrac1-ffrac2)^2)
			denominator = denominator + sum(ffrac1^2+ffrac2^2)
			}
		
		# calculate the score value
		FSSwind[inn]=1-numerator/denominator	
		}

	return(FSSwind)
	}








