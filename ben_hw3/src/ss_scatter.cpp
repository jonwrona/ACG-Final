//this is where I will write my subsurface scattering code. Hopefully I cna #include it, but if not I will paste it in
 
 //hit is the hit data from the 'eye ray'. xo is the position. wo is the direction this was hit at.
glm::vec3 RayTracer::ss_scatter(Hit &hit, glm::vec3 xo, glm::vec3 lightPos, glm::vec3 wo, glm::vec3 lightColor){
	//set all the constants

	//for now assume we have a point light source. The first thing we must do is distribute 
	//the incident contributing rays
	glm::vec3 answer(0.0,0.0,0.0);
	glm::vec3 ones(1.0,1.0,1.0);	
	wo=-wo;
	//declare the stuff we do have
	glm::vec3 no = hit.getNormal(); //this is the outgoing surface normal
	float sigma_t=args->sigma_t_prime.y; //this corresponds to the green value sigma_t;
	float nu=args->nu;
	float e=2.718281828459045;
	float sigma_s=args->sigma_s_prime.y;


	//declare the stuff we need to find
	glm::vec3 xInside; //this is the position under the surface of the 'intersection'
	glm::vec3 ni; //this is going to be the incident surface normal
	glm::vec3 xi;
	glm::vec3 wi;


	//first we sample for the diffusion term. This follows
	for (int i=0; i<args->num_ss_samples /* need to add num_ss_samples to the argparser */ ; i++){
		//first compute a point along the plane perpendicular to the normal
		//float disp_length = log(random_length())*ones/args->sigma_t;
		//now calcukate the other terms of the equation
		//glm::vec3 si=

	}

	//next we sample for the single scattering term. This follows
	for (int i=0; i<args->num_ss_samples /* need to add num_ss_samples to the argparser */ ; i++){
		//first send a ray into the surface some distance. For now don't refract. CHANGE THISSSSSS.
		//We are also going to assume that all the colors travel the same distance for now. We can change this
		//by repeating this calc for R, G, and B. We are using the G value for everything right now

		//determing so_prime
		float so_prime=log(GLOBAL_MTRAND())/sigma_t;
		xInside=xo+so_prime*glm::normalize(wo); //need to implement random number and refrect here...
		Ray r_inside(xInside, glm::normalize(lightPos-xInside));
		Hit surface_hit;
		CastRay(r_inside, surface_hit, false);
		//if we did not manage to pass back through the surface the contribution is 0
		if(fabs(surface_hit.getT()-glm::length(lightPos-xInside))<0.0001){
			answer+=0.0;
			continue;
		}


		//here we determine sigma_tc. Once we have the refracted directions this will change
		sigma_tc=sigma_t+glm::length(ni*glm::normalize(lightPos-xi))/glm::length(ni*glm::normalize(lightPos-xo))*sigma_t;

		//if we did pass back through the surface get xi and the normal
		xi=xInside+r_inside.getDirection()*surface_hit*getT();
		ni=surface_hit.getNormal();
		wi=glm::normalize(r_inside.getDirection());

		//si is the percieved distance
		float si=glm::length(xi-xo);
		//si prime is an estimate of the refracted distance through the material
		si_prime=si*(glm::dot(wi, ni)/(float)sqrt(1.0-(float)(1.0/(double)nu)*(float)(1.0/(double)nu)*(float)(1.0-pow((double)glm::dot(wi,ni),2.0)))) ; //not sure what they meant by n(xi) in the paper...
		float F=Fresnel_transmittance(nu, wo, no)*Fresnel_transmittance(nu, wi, ni); //assuming fresnel transmittance is 1-fresnel reflection
		

		//compute the incoming radiance at xi
		float distToLight=glm::length(lightPos-xi);
		glm::vec3 Li=(glm::dot(ni, glm::normalize(lightPos-xi)))/(12.5663706144f*distToLight*distToLight)*lightColor;
		FresnelTerm=Fresnel_transmittance(nu, glm::normalize(lightPos-xi), ni)*Fresnel_transmittance(nu, glm::normalize(lightPos-xo), no);
		Lo=sigma_s*FresnelTerm/sigma_tc*pow(e, -si_prime*sigma_t)*pow(e, -so_prime)*Li;
		answer=Lo;
	}

	//average answer out
	//answer/=float(args->num_ss_samples);
	return answer;

}



