function auto_integrator_trap_romb_hnm(func,a,b,nmax,tol_ae,tol_rae)
// INPUTS
// func=integrand
// a= lower limit of integration
// b= upper limit of integration
// nmax = number of partitions, n=2^nmax
// tol_ae= maximum absolute approximate error acceptable (should be >=0)
// tol_rae=maximum absolute relative approximate error acceptable (should be >=0)


// a= lower limit of integration
	// b= upper limit of integraton
	// nmax= number of partitions, segment is then 2^nmax 
	// tol_ae= tolerance on absolute approximate error 
	// tol_rae=tolerance on percentage absolute relative approximate error 

// OUTPUTS
// integ_value= estimated value of integral

{
//Checking for input errors
	if (typeof a !== 'number') 
		{
		  alert('<lower limit of integration> must be a number');
		  throw new TypeError('<a> must be a number');
		  
		}
    if (typeof b !== 'number') 
		{
		  alert('<upper limit of integration> must be a number');		
		  throw new TypeError('<b> must be a number');
		}
    if ((!Number.isInteger(nmax)) || (nmax<1))
		{
		  alert('<number of partitions, n=2^nmax> must be an integer greater than or equal to one.');	
		  throw new TypeError('<nmax> must be an integer greater than or equal to one.');
		}
	if ((typeof tol_ae !== 'number') || (tol_ae<0)) 
		{
		  alert('<tolerance on absolute approximate error (should be >=0)> must be a number greater than or equal to zero');		 
		  throw new TypeError('<tole_ae> must be a number greater than or equal to zero');
		}
	if ((typeof tol_rae !== 'number') || (tol_rae<=0)) 
		{
		  alert('<tolerance on percentage absolute relative approximate error (should be >=0)> must be a number greater than or equal to zero');			 
		  throw new TypeError('<tole_ae> must be a number greater than or equal to zero');
		}
    
	var h=b-a
	// initialize matrix where the values of integral are stored
	
	var Romb = []; // rows
	for (var i = 0; i < nmax+1; i++) 
	{
		Romb.push([]);
		for (var j = 0; j < nmax+1; j++) 
		{
			Romb[i].push(math.bignumber(0)); 
		}
	}
	
	//calculating the value with 1-segment trapezoidal rule
	Romb[0][0]=0.5*h*(func(a)+func(b))
	var integ_val=Romb[0][0]
	
	for (var i=1; i<=nmax; i++)
	// updating the value with double the number of segments
	// by only using the values where they need to be calculated
	// See https://autarkaw.org/2009/02/28/an-efficient-formula-for-an-automatic-integrator-based-on-trapezoidal-rule/
	{
		h=0.5*h
		var integ=0
		for (var j=1; j<=2**i-1; j+=2)
		{
			var integ=integ+func(a+j*h)
		}
	
		Romb[i][0]=0.5*Romb[i-1][0]+integ*h
		// Using Romberg method to calculate next extrapolatable value
		// See https://young.physics.ucsc.edu/115/romberg.pdf
		for (k=1; k<=i; k++)
		{   
			var addterm=Romb[i][k-1]-Romb[i-1][k-1]
			addterm=addterm/(4**k-1.0)
			Romb[i][k]=Romb[i][k-1]+addterm

			//Calculating absolute approximate error
			var Ea=math.abs(Romb[i][k]-Romb[i][k-1])
			
			//Calculating absolute relative approximate error
			var epsa=math.abs(Ea/Romb[i][k])*100.0
			
			//Assigning most recent value to the return variable
			integ_val=Romb[i][k]
			
			// returning the value if either tolerance is met
			if ((epsa<tol_rae) || (Ea<tol_ae))
			{
				return(integ_val)
			}
		}
	}
	// returning the last calculated value of integral whether tolerance is met or not
	return(integ_val)
}