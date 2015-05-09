#include "glCanvas.h"

#include "raytracer.h"
#include "material.h"
#include "argparser.h"
#include "raytree.h"
#include "utils.h"
#include "mesh.h"
#include "face.h"
#include "primitive.h"
#include "photon_mapping.h"

//this function is mostly used for debugging, just prints vecs
void print_vec(const glm::vec3 &a){
  std::cout<<a.x<<" "<<a.y<<" "<<a.z<<std::endl;
}

 //hit is the hit data from the 'eye ray'. xo is the position. wo_ is the direction this was hit a, but we negate it.
glm::vec3 RayTracer::ss_scatter(const Hit &hit, const glm::vec3 xo, const glm::vec3 lightPos_,const glm::vec3 wo_,const glm::vec3 lightColor_, const Face* f) const{
  //set all the constants
  //print_vec(args->sigma_s_prime);

  glm::vec3 answer(0.f,0.f,0.f);
  //we have to take seperate samples for R, G, and B
  for (int j=0; j<3; j++){
    glm::vec3 pseudo_answer(0.f, 0.f, 0.f);
    glm::vec3 pseudo_answer_diff(0.f,0.f,0.f);
    int count=0;
    //these next three terms are color dependent, so we have three seperate cases
    float sigma_t;
    float sigma_s;
    float sigma_a;
    glm::vec3 lightColor;
    if(j==0){
      sigma_t=args->sigma_t_prime.x; //this corresponds to the green value sigma_t;
      sigma_s=args->sigma_s_prime.x;
      sigma_a=args->sigma_a.x;
      lightColor=glm::vec3(lightColor_.x, 0.0f,0.0f);
      //std::cout<<"first lightColor is ";
      //print_vec(lightColor);

    }
    else if (j==1){
      sigma_t=args->sigma_t_prime.y; //this corresponds to the green value sigma_t;
      sigma_s=args->sigma_s_prime.y;
      sigma_a=args->sigma_a.y;
      //lightColor=glm::vec3(lightColor_.x, lightColor_.y, lightColor_.z);
      lightColor=glm::vec3(0.f,lightColor_.y,0.f);
      //std::cout<<"second lightColor is ";
      //print_vec(lightColor);
    }
    else if (j==2){
      sigma_t=args->sigma_t_prime.z; //this corresponds to the green value sigma_t;
      sigma_s=args->sigma_s_prime.z;
      sigma_a=args->sigma_a.z;
      lightColor=glm::vec3(0.f,0.f,lightColor_.z);
      //std::cout<<"third lightColor is ";
      //print_vec(lightColor);

    }
  
    //for now assume we have a point light source. The first thing we must do is distribute 
    //the incident contributing rays
    
    glm::vec3 ones(1.0,1.0,1.0);  
    glm::vec3 wo=-wo_;
    //declare the stuff we do have
    glm::vec3 no = glm::normalize(hit.getNormal()); //this is the outgoing surface normal

    float nu=args->nu;
    float e=2.718281828459045;

    //declare the stuff we need to find
    glm::vec3 xInside; //this is the position under the surface of the 'intersection'
    glm::vec3 ni; //this is going to be the incident surface normal
    glm::vec3 xi;
    glm::vec3 wi;

    //diff_answer is the diffusion term
    int count_diff=0;

    #if 1
    //first we sample for the diffusion term. This follows
    for (int i=0; i<args->num_ss_samples; i++){
      //first compute the coefficents we need
      float A=(1.f-args->Fdr)/(1.f+args->Fdr);
      float D=1.f/(3.f*sigma_t);
      float sigma_tr = glm::sqrt((double)(3.f*sigma_a*sigma_t));  
      
      //calculate xi
      glm::vec3 xi;

      //tang is the unit vector in the direction of the eye ray hit that is tangent to the surface
      glm::vec3 tang=glm::normalize(glm::normalize(wo_)+no*(float)(sin(1.57079632679-acos(glm::dot(wo, no)))));
      float rand_dist=-sigma_tr*(float)log(GLOBAL_MTRAND()/sigma_tr)*.1f;
      //std::cout<<"rand dist is "<<rand_dist<<std::endl;
      //  std::cout<<glm::dot(no, tang)<<std::endl;
      xi=xo+tang*rand_dist; //NEED TO CHANGE THIS
      Ray r_find(xi, -no);
      Hit h_find, h_find2;
      //this is the direction to the surface
      glm::vec3 surface_dir;
      //find the surface
      CastRay(r_find, h_find, false);
      RayTree::AddReflectedSegment (r_find, 0, h_find.getT());
      surface_dir=-no;

      //std::cout<<rand_dist<<std::endl;  
      #if 0
      if (! glm::length(h_find.getMaterial()->getReflectiveColor())>0.1){
        //we went the wrong way
        //std::cout<<"Wrong way"<<std::endl;
        Ray r_find2(xi, no);
        CastRay(r_find, h_find2, false);
        surface_dir=-no;
      }
      if(h_find2.getT()>rand_dist){
        //if we are here we did not find the surface
        std::cout<<"we did not find the surface\n";
        continue;
      }
      #endif
      //std::cout<<"yay\n";
      //set xi and ni
      //std::cout<<h_find.getT()<<std::endl;
      if(h_find.getT()>rand_dist*1.5f){
        continue;
      }
      xi=xi+surface_dir*h_find.getT();
      //print_vec(xi);
      glm::vec3 ni =h_find.getNormal();
      //print_vec(ni);
      //print_vec(xi);
      //calculate zr and zv, the distances above and below the surface we place the 'lights'
      float zr=1.f/sigma_t;      
      float zv=zr+4.f*A*D;
      //std::cout<<"sigma_t is "<<sigma_t<<" and zr is "<<zr<<" and zv is "<<zv<<std::endl;

      //calculate the real and virtual sources
      glm::vec3 real_source = xi-ni*zr;
      glm::vec3 virtual_source = xi+ni*zv;
      //std::cout<<"ni is:"<<std::endl;
      //print_vec(ni);

      //declare the stuff we already know
      float albedo=args->albedo;
      

      float dr=glm::length(xo-real_source);
      float dv=glm::length(xo-virtual_source);
      //std::cout<<"dr is "<<dr<<" and dv is "<<dv<<std::endl;

      //determine the contribution
      float Rd=albedo/(12.5663706144)*((sigma_tr*dr+1.f)* ((float)exp((double)(-sigma_tr*dr))/(sigma_t*dr*dr*dr)) + ( zv*(sigma_tr*dv+1.f) * ( exp((double)-sigma_tr*dv)/(sigma_t*dv*dv*dv) )) );
      //std::cout<<"Term 1 is "<<(sigma_tr*dr+1.f)<<" term 2 is "<<((float)exp((double)(-sigma_tr*dr))/(sigma_t*dr*dr*dr))<<std::endl;
      //std::cout<<"and term 3 is "<<( zv*(sigma_tr*dv+1.f) * ( exp((double)-sigma_tr*dv)/(sigma_t*dv*dv*dv) ))<<std::endl;
      glm::vec3 wi=glm::normalize(lightPos_-xi);
      float Fo=Fresnel_transmittance(nu, wo_, no);
      float Fi=Fresnel_transmittance(1.f/nu, -wi, ni);
      float FresnelTerm=Fo+Fi;
      // /std::cout<<"Fo is "<<Fo<<" and Fi is "<<Fi<<" and nu is "<<nu<<std::endl;
       if( ! isnan(Rd) and ! isinf(Rd)){
        //std::cout<<"Rd seems good"<<std::endl;
        float term=FresnelTerm*Rd/(glm::length(xi-lightPos_ )*glm::length(xi-lightPos_ )*(float)12.5663706144)*((float)glm::dot(xi-lightPos_, -ni));
        //std::cout<<"Term is "<<term<<std::endl;
        if(term>0){ //really should be if term>0
          pseudo_answer_diff+=term*lightColor;
          pseudo_answer_diff*=.1f;
          count_diff+=1;
        }
      }
      else{
        //std::cout<<"We got a nan..."<<std::endl;
      }
    }
    
    count_diff+=1;
    //
    //print_vec(pseudo_answer_diff);
    pseudo_answer_diff/=(float)count_diff;
    #endif
    //std::cout<<count_diff<<std::endl;

    #if 1
    //next we sample for the single scattering term. This follows
    for (int i=0; i<args->num_ss_samples; i++){

      //first send a ray into the surface some distance. For now don't refract. CHANGE THISSSSSS.
      //We are also going to assume that all the colors travel the same distance for now. We can change this
      //by repeating this calc for R, G, and B. We are using the G value for everything right now
      glm::vec3 lightPos=f->RandomPoint();
      //first calculate the refracted direction below the surface
      float c0 = glm::dot(-wo, no);
      float c1 = (float)sqrt(1.0 - pow(nu, 2.0) * (1.0 - pow(c0, 2.0)));
      glm::vec3 wo_prime = -glm::normalize((1.f/nu * -wo) + ((nu * c0 - c1) * no));
      //print_vec(wo_prime);
      //determing so_prime
      float so_prime=log(GLOBAL_MTRAND())/sigma_t;
      //std::cout<<so_prime<<std::endl;
      xInside=xo+so_prime*wo_prime; //need to implement random number and refract here...
      Ray r_inside(xInside, glm::normalize(lightPos-xInside));
      Hit surface_hit;
      CastRay(r_inside, surface_hit, false);
      //if we did not manage to pass back through the surface the contribution is 0
      if(fabs(surface_hit.getT()-glm::length(lightPos-xInside))<0.00001){
        //answer+=glm::vec3(0.0,0.0,0.0);
        count+=1;
        continue;
      }


      //here we determine sigma_tc. Once we have the refracted directions this will change
      //calculate the refracted direction into the surface
      float sigma_tc=sigma_t+glm::length(ni*glm::normalize(lightPos-xi))/glm::length(ni*glm::normalize(wo_prime))*sigma_t;

      //if we did pass back through the surface get xi and the normal
      xi=xInside+r_inside.getDirection()*surface_hit. getT();
      ni=surface_hit.getNormal();
      wi=glm::normalize(r_inside.getDirection());

      //si is the percieved distance
      float si=glm::length(xi-xo);
      //si prime is an estimate of the refracted distance through the material
      float si_prime=si*(glm::dot(wi, ni)/(float)sqrt(1.0-(float)(1.0/(double)nu)*(float)(1.0/(double)nu)*(float)(1.0-pow((double)glm::dot(wi,ni),2.0)))) ; //not sure what they meant by n(xi) in the paper...
      //float F=Fresnel_transmittance(nu, wo, no)*Fresnel_transmittance(nu, wi, ni); //assuming fresnel transmittance is 1-fresnel reflection
      

      //compute the incoming radiance at xi
      float distToLight=glm::length(lightPos-xi);
      glm::vec3 Li=(glm::dot(ni, glm::normalize(lightPos-xi)))/(12.5663706144f*distToLight*distToLight)*lightColor;
      //Li*=5.0;
      //std::cout<<"Li is:";
      //print_vec(Li);
      float FresnelTerm=Fresnel_transmittance(1.f/nu, glm::normalize(lightPos-xi), -ni)*Fresnel_transmittance(nu, -wo_prime, no);
      //std::cout<<FresnelTerm<<std::endl;
      glm::vec3 Lo=sigma_s*FresnelTerm/sigma_tc*(float)pow(e, -si_prime*sigma_t)*(float)pow(e, -so_prime)*Li;
      //std::cout<<"Lo is:";
      //print_vec(Lo);
      if( ! isnan(Lo.x) and ! isinf(Lo.x)){
        count+=1;
        pseudo_answer+=Lo;
      }
      else{
        //std::cout<<j<<std::endl;
        //answer+=answer/((float)i-1.f);
      }
    }
    #endif

    
    if(j==1){
      //std::cout<<"Count was "<<count<<std::endl;
    }
    //count++;
    pseudo_answer/=(float)count;
    answer+=pseudo_answer;
    #if 0
    std::cout<<"pseudo_answer is: "<<std::endl;
    print_vec(pseudo_answer);
    std::cout<<"diff_answer is "<<std::endl;
    print_vec(pseudo_answer_diff);
    #endif
    answer+=pseudo_answer_diff;
  }


  //average answer out
  //std::cout<<"num ss samples is "<<(float)(args->num_ss_samples-1  )<<std::endl;
  //answer/=(float)(args->num_ss_samples-1  );
  //std::cout<<"answer is ";
  //print_vec(answer);

  return answer*10.f;

}

//add my new functions here
glm::vec3 sqrt(const glm::vec3 &x){
  assert (x.x>=0.0 and x.y>=0.0 and x.z>=0.0);
  glm::vec3 ans(sqrt(x.x), sqrt(x.y), sqrt(x.z));
  return ans;
}

float safe_length(glm::vec3 a){

  return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
}

//here w is the direction to the light, n1 is the outward facing normal of the surface, and nu is the index of refraction
float RayTracer::Fresnel_transmittance(const float &nu, const glm::vec3 &w_, const glm::vec3 &n1_) const{

  glm::vec3 n1=glm::normalize(n1_);
  glm::vec3 w= glm::normalize(w_);
  glm::vec3 n2=-n1;
  glm::vec3 ones(1.0,1.0,1.0);

  float theta=glm::acos(glm::dot(n1, w)); //might be wrong...  
  //assert(sin(theta)<1.571f and sin(theta)>0.f);
  if(! 1.0-sin(theta)*sin(theta)>0.0){
    return 0;
  }
  float sqrt_term=(float)sqrt(1.f+sin(theta)*sin(theta));
  //std::cout<<"sqrt_term is "<<sqrt_term<<std::endl;
  //std::cout<< cos(theta) <<std::endl;
  //std::cout<<"we're here\n";
  glm::vec3 RS_VEC=((float)cos(theta)*n1-sqrt_term*n2)/(n1*(float)cos(theta)+n2*sqrt_term);
  //std::cout<<"RS_VEC is "<<RS_VEC.x<<" "<<RS_VEC.y<<" "<<RS_VEC.z<<std::endl;
  float Rs =  safe_length(((float)cos(theta)*n1-sqrt_term*n2)/(n1*(float)cos(theta)+n2*sqrt_term));
  Rs*=Rs;

  float Rp =  safe_length((sqrt_term*n1-(float)cos(theta)*n2)/(sqrt_term*n1+(float)cos(theta)*n2));  
  Rp*=Rp;
  float R=(Rs+Rp)/2.f;
  float T;
  if (R>1.f){
    T=0.f;
  }
  else{
    T=1.f-R;
  }
  //std::cout<<"RS is "<<Rs<<" and Rp is "<<Rp<<std::endl;
  //std::cout<<"we're here\n";
  return T;
}

// ===========================================================================
// casts a single ray through the scene geometry and finds the closest hit
bool RayTracer::CastRay(const Ray &ray, Hit &h, bool use_rasterized_patches) const {
  bool answer = false;

  // intersect each of the quads
  for (int i = 0; i < mesh->numOriginalQuads(); i++) {
    Face *f = mesh->getOriginalQuad(i);
    if (f->intersect(ray,h,args->intersect_backfacing)) answer = true;
  }

  // intersect each of the primitives (either the patches, or the original primitives)
  if (use_rasterized_patches) {
    for (int i = 0; i < mesh->numRasterizedPrimitiveFaces(); i++) {
      Face *f = mesh->getRasterizedPrimitiveFace(i);
      if (f->intersect(ray,h,args->intersect_backfacing)) answer = true;
    }
  } else {
    int num_primitives = mesh->numPrimitives();
    for (int i = 0; i < num_primitives; i++) {
      if (mesh->getPrimitive(i)->intersect(ray,h)) answer = true;
    }
  }
  return answer;
}

// ===========================================================================
// does the recursive (shadow rays & recursive rays) work
// the default index of refraction is set to 1.000277 (air)
glm::vec3 RayTracer::TraceRay(Ray &ray, Hit &hit, int bounce_count, bool inside) const {

  // First cast a ray and see if we hit anything.
  hit = Hit();
  bool intersect = CastRay(ray,hit,false);
    
  // if there is no intersection, simply return the background color
  if (intersect == false) {
    return glm::vec3(srgb_to_linear(mesh->background_color.r),
                     srgb_to_linear(mesh->background_color.g),
                     srgb_to_linear(mesh->background_color.b));
  }
  //try adding the ray...


  // otherwise decide what to do based on the material
  Material *m = hit.getMaterial();
  assert (m != NULL);

  // rays coming from the light source are set to white, don't bother to ray trace further.
  if (glm::length(m->getEmittedColor()) > 0.001) {
    return glm::vec3(1,1,1);
  } 
 
  
  glm::vec3 normal = hit.getNormal();
  glm::vec3 point = ray.pointAtParameter(hit.getT());
  glm::vec3 answer;

  // ----------------------------------------------
  //  start with the indirect light (ambient light)
  glm::vec3 diffuse_color = m->getDiffuseColor(hit.get_s(),hit.get_t());
  if (args->gather_indirect) {
    // photon mapping for more accurate indirect light
    answer = diffuse_color * (photon_mapping->GatherIndirect(point, normal, ray.getDirection()) + args->ambient_light);

  } else {
    // the usual ray tracing hack for indirect light
    if(args->ss_scatter) {
      //answer = ss_scatter(hit);
    }
    else{
      answer = diffuse_color * args->ambient_light;
    }
  }      

  // ----------------------------------------------
  // add contributions from each light that is not in shadow
  int num_lights = mesh->getLights().size();
  for (int i = 0; i < num_lights; i++) {

    Face *f = mesh->getLights()[i];
    glm::vec3 lightColor = f->getMaterial()->getEmittedColor() * f->getArea();
    glm::vec3 myLightColor;
    glm::vec3 lightCentroid = f->computeCentroid();
    glm::vec3 dirToLightCentroid = glm::normalize(lightCentroid-point); 
    //if we are trying to visualize the fresnel transmittance
    //this if statement is a hack for now
    if (m->isSubsurfaceMaterial()) std::cout << "SSS" << std::endl;
    if(args->ss_scatter and glm::length(m->getReflectiveColor())>0.1f){
      #if 0
      float F = this->Fresnel_transmittance(args->nu, -dirToLightCentroid, normal);
      F = F*F;
      //std::cout<<F<<std::endl;
      glm::vec3 tmp_answer(F,F,F);
      return tmp_answer;
      #endif
      #if 1
      glm::vec3 tmp = ss_scatter(hit, point, lightCentroid, ray.getDirection(), lightColor, f);
      //print_vec(tmp);
      return tmp;
      #endif
    }

    for (int j = 0; j < args->num_shadow_samples; j++) {
      glm::vec3 direction = dirToLightCentroid;
      glm::vec3 lightPoint = lightCentroid;
      if (args->num_shadow_samples > 1) {
        // generate random point on light
        lightPoint = f->RandomPoint();
        // get direction towards said point
        direction = glm::normalize(lightPoint - point);
      }
      Ray shadowRay = Ray(point, direction);
      Hit shadowHit = Hit();

      CastRay(shadowRay, shadowHit, false);
      bool blocked = shadowHit.getMaterial() != f->getMaterial();
      bool transparent = shadowHit.getMaterial()->getRefraction() > 0.0f;
      if (!blocked || transparent) {
        float distToLightPoint = glm::length(lightPoint-point);
        myLightColor = lightColor / float(M_PI*distToLightPoint*distToLightPoint);
        if (args->num_shadow_samples > 1) myLightColor /= args->num_shadow_samples;

        // add the lighting contribution from this particular light at this point
        // (fix this to check for blockers between the light & this surface)
        answer += m->Shade(ray,hit,direction,myLightColor,args);
        RayTree::AddShadowSegment(shadowRay, 0, shadowHit.getT());
      }
    }
  }

    // ===========================================
    // ASSIGNMENT:  ADD SHADOW & SOFT SHADOW LOGIC
    // ===========================================


  // ----------------------------------------------
  // add contribution from reflection, if the surface is shiny
  glm::vec3 reflectiveColor = m->getReflectiveColor();


  // // =================================
  // // ASSIGNMENT:  ADD REFLECTIVE LOGIC
  // // =================================
  //std::cout << glm::length(reflectiveColor) << std::endl;
  if (glm::length(reflectiveColor) > 0.001) {
    glm::vec3 ray_dir=glm::normalize(ray.getDirection());
    normal=glm::normalize(normal);
    glm::vec3 reflected_dir=glm::normalize(ray_dir-2*glm::dot(ray_dir, normal)*normal);
    //next we trace this ray and get the color reflected onto it
    Ray reflected_ray(point,reflected_dir);
    Hit reflection_hit;
    answer += glm::normalize(reflectiveColor) * TraceRay(reflected_ray, reflection_hit, bounce_count-1);
    RayTree::AddReflectedSegment(reflected_ray, 0, reflection_hit.getT());
  }

  // =================================
  // ASSIGNMENT:  ADD REFRACTIVE LOGIC
  // =================================
  float refraction = m->getRefraction();
  if (refraction > 0.0 && bounce_count > 0) {
    float air = 1.000277f;
    glm::vec3 incident = glm::normalize(ray.getDirection());

    float n = air / refraction;
    float n1 = air;
    float n2 = refraction;
    if (inside) {
      normal = -1.0f * normal;
      n = refraction / air;
      n1 = refraction;
      n2 = air;
    }

    float cosI = -1.0f * glm::dot(normal, incident);
    float sinT2 = n * n * (1.0f - cosI * cosI);
    float reflectivity;
    if (sinT2 > 1.0f) reflectivity = 1.0f;
    else {
      float cosT = sqrt(1.0f - sinT2);
      float rOrth = (n1 * cosI - n2 * cosT) / (n1 * cosI + n2 * cosT);
      float rPar = (n2 * cosI - n1 * cosT) / (n2 * cosI + n1 * cosT);
      reflectivity = (rOrth * rOrth + rPar + rPar) / 2.0f;
    }
    float transitivity = 1.0f - reflectivity;
    if (transitivity > 0.0) {
      glm::vec3 direction;
      float cosI = -1.0f * glm::dot(normal, incident);
      float sinT2 = n * n * (1.0f - cosI * cosI);

      if (sinT2 < 1.0) {
        float cosT = sqrt(1.0f - sinT2);
        direction = n * incident + (n * cosI - cosT) * normal;
        direction = glm::normalize(direction);
        Ray refractRay = Ray(point, glm::normalize(direction));
        Hit refractHit = Hit();
        answer += transitivity * TraceRay(refractRay, refractHit, bounce_count-1, !inside);
        RayTree::AddTransmittedSegment(refractRay, 0, refractHit.getT());
      }
    }
    if (reflectivity > 0.0 && !inside) {
      float cosI = -1.0f * glm::dot(normal, incident);
      glm::vec3 direction = incident + 2 * cosI * normal;
      Ray reflectRay = Ray(point, glm::normalize(direction));
      Hit reflectHit = Hit();
      answer += reflectivity * TraceRay(reflectRay, reflectHit, bounce_count-1);
      RayTree::AddReflectedSegment(reflectRay, 0, reflectHit.getT());
    }
  }

  
  return answer; 

}



void RayTracer::initializeVBOs() {
  glGenBuffers(1, &pixels_a_VBO);
  glGenBuffers(1, &pixels_b_VBO);
  glGenBuffers(1, &pixels_indices_a_VBO);
  glGenBuffers(1, &pixels_indices_b_VBO);
  render_to_a = true;
}


void RayTracer::resetVBOs() {

  pixels_a.clear();
  pixels_b.clear();

  pixels_indices_a.clear();
  pixels_indices_b.clear();

  render_to_a = true;
}

void RayTracer::setupVBOs() {

  glBindBuffer(GL_ARRAY_BUFFER,pixels_a_VBO); 
  glBufferData(GL_ARRAY_BUFFER,sizeof(VBOPosNormalColor)*pixels_a.size(),&pixels_a[0],GL_STATIC_DRAW); 
  glBindBuffer(GL_ARRAY_BUFFER,pixels_b_VBO); 
  glBufferData(GL_ARRAY_BUFFER,sizeof(VBOPosNormalColor)*pixels_b.size(),&pixels_b[0],GL_STATIC_DRAW); 

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,pixels_indices_a_VBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
	       sizeof(VBOIndexedTri) * pixels_indices_a.size(),
	       &pixels_indices_a[0], GL_STATIC_DRAW);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,pixels_indices_b_VBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
	       sizeof(VBOIndexedTri) * pixels_indices_b.size(),
	       &pixels_indices_b[0], GL_STATIC_DRAW);

}

void RayTracer::drawVBOs() {
  // turn off lighting
  glUniform1i(GLCanvas::colormodeID, 0);
  // turn off depth buffer
  glDisable(GL_DEPTH_TEST);

  if (render_to_a) {
    drawVBOs_b();
    drawVBOs_a();
  } else {
    drawVBOs_a();
    drawVBOs_b();
  }

  glEnable(GL_DEPTH_TEST);
}

void RayTracer::drawVBOs_a() {
  if (pixels_a.size() == 0) return;
  glBindBuffer(GL_ARRAY_BUFFER, pixels_a_VBO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,pixels_indices_a_VBO); 
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)0);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));
  glEnableVertexAttribArray(3);
  glVertexAttribPointer(3, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2 + sizeof(glm::vec4)));
  glDrawElements(GL_TRIANGLES,
                 pixels_indices_a.size()*3,
                 GL_UNSIGNED_INT, 0);
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);
  glDisableVertexAttribArray(3);
}

void RayTracer::drawVBOs_b() {
  if (pixels_b.size() == 0) return;
  glBindBuffer(GL_ARRAY_BUFFER, pixels_b_VBO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,pixels_indices_b_VBO); 
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)0);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));
  glEnableVertexAttribArray(3);
  glVertexAttribPointer(3, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2 + sizeof(glm::vec4)));
  glDrawElements(GL_TRIANGLES,
                 pixels_indices_b.size()*3,
                 GL_UNSIGNED_INT, 0);
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);
  glDisableVertexAttribArray(3);
}


void RayTracer::cleanupVBOs() {
  glDeleteBuffers(1, &pixels_a_VBO);
  glDeleteBuffers(1, &pixels_b_VBO);
}
