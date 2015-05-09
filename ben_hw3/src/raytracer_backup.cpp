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

 //hit is the hit data from the 'eye ray'. xo is the position. wo is the direction this was hit at.
glm::vec3 RayTracer::ss_scatter(const Hit &hit, const glm::vec3 xo, const glm::vec3 lightPos_,const glm::vec3 wo_,const glm::vec3 lightColor_, const Face* f) const{
  //set all the constants
  //print_vec(args->sigma_s_prime);

  glm::vec3 answer(0.f,0.f,0.f);
  //we have to take seperate samples for R, G, and B
  for (int j=0; j<3; j++){
    glm::vec3 pseudo_answer(0.f, 0.f, 0.f);
    int count=0;
    //these next three terms are color dependent, so we have three seperate cases
    float sigma_t;
    float sigma_s;
    glm::vec3 lightColor;
    if(j==0){
      sigma_t=args->sigma_t_prime.x; //this corresponds to the green value sigma_t;
      sigma_s=args->sigma_s_prime.x;
      lightColor=glm::vec3(lightColor_.x, 0.f,0.f);
      //std::cout<<"first lightColor is ";
      //print_vec(lightColor);

    }
    else if (j==1){
      sigma_t=args->sigma_t_prime.y; //this corresponds to the green value sigma_t;
      sigma_s=args->sigma_s_prime.y;
      //lightColor=glm::vec3(lightColor_.x, lightColor_.y, lightColor_.z);
      lightColor=glm::vec3(0.f,lightColor_.y,0.f);
      //std::cout<<"second lightColor is ";
      //print_vec(lightColor);
    }
    else if (j==2){
      sigma_t=args->sigma_t_prime.z; //this corresponds to the green value sigma_t;
      sigma_s=args->sigma_s_prime.z;
      lightColor=glm::vec3(0.f,0.f,lightColor_.z);
      //std::cout<<"third lightColor is ";
      //print_vec(lightColor);

    }
  
    //for now assume we have a point light source. The first thing we must do is distribute 
    //the incident contributing rays
    
    glm::vec3 ones(1.0,1.0,1.0);  
    glm::vec3 wo=-wo_;
    //declare the stuff we do have
    glm::vec3 no = hit.getNormal(); //this is the outgoing surface normal

    float nu=args->nu;
    float e=2.718281828459045;

    //declare the stuff we need to find
    glm::vec3 xInside; //this is the position under the surface of the 'intersection'
    glm::vec3 ni; //this is going to be the incident surface normal
    glm::vec3 xi;
    glm::vec3 wi;

    //first we sample for the diffusion term. This follows
    for (int i=0; i<args->num_ss_samples; i++){
      //first compute a point along the plane perpendicular to the normal
      //float disp_length = log(random_length())*ones/args->sigma_t;
      //now calcukate the other terms of the equation
      //glm::vec3 si=

    }

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
      if(fabs(surface_hit.getT()-glm::length(lightPos-xInside))<0.0001){
        //answer+=glm::vec3(0.0,0.0,0.0);
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
      Li*=5.0;
      //std::cout<<"Li is:";
      //print_vec(Li);
      float FresnelTerm=Fresnel_transmittance(nu, -glm::normalize(lightPos-xi), ni)*Fresnel_transmittance(nu, -glm::normalize(lightPos-xo), no);
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
    
    if(j==1){
      //std::cout<<"Count was "<<count<<std::endl;
    }
    count++;
    pseudo_answer/=(float)count;
    answer+=pseudo_answer;
  }

  //average answer out
  //std::cout<<"num ss samples is "<<(float)(args->num_ss_samples-1  )<<std::endl;
  //answer/=(float)(args->num_ss_samples-1  );
  //std::cout<<"answer is ";
  //print_vec(answer);

  return answer;

}

//add my new functions here
glm::vec3 sqrt(const glm::vec3 &x){
  assert (x.x>=0.0 and x.y>=0.0 and x.z>=0.0);
  glm::vec3 ans(sqrt(x.x), sqrt(x.y), sqrt(x.z));
  return ans;
}

float safe_length(glm::vec3 a){
  /*
  glm::vec3 zeros(0.0,0.0,0.0);
  std::cout<<"in safe length\n";
  float answer;
  if (a==zeros){
    
    answer= 0;
    std::cout<<"in safe length\n";
  }
  else{
    
    answer= glm::length(a);
    std::cout<<"in safe length\n";
  }
  */
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
    T=0;
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
glm::vec3 RayTracer::TraceRay(Ray &ray, Hit &hit, int bounce_count) const {

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
    //this is my hard shadows code
    #if 1
    glm::vec3 lightColor = f->getMaterial()->getEmittedColor() * f->getArea();
    glm::vec3 myLightColor;
    glm::vec3 lightCentroid = f->computeCentroid();
    glm::vec3 dirToLightCentroid = glm::normalize(lightCentroid-point); 
    //if we are trying to visualize the fresnel transmittance
    //this if statement is a hack for now
    if(args->ss_scatter and glm::length(m->getReflectiveColor())>0.01){
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
    // ===========================================
    // ASSIGNMENT:  ADD SHADOW & SOFT SHADOW LOGIC
    // ===========================================
    
    float distToLightCentroid = glm::length(lightCentroid-point);
    myLightColor = lightColor / float(M_PI*distToLightCentroid*distToLightCentroid);
    
    // add the lighting contribution from this particular light at this point
    // (fix this to check for blockers between the light & this surface)
    
    //hard shadow stuff here
    //cast a ray towards the light, if we hit it then light this up
    Hit tmp_hit;
    Ray lightRay(point, dirToLightCentroid);

    CastRay(lightRay,tmp_hit,false);
    if((tmp_hit.getT()<distToLightCentroid+0.001 and tmp_hit.getT()>distToLightCentroid-0.001)){
      answer += m->Shade(ray,hit,dirToLightCentroid,myLightColor,args);
    }

      
    #endif
    //this is my soft shadows code
    #if 0
    for (int i=0; i<4; i++){
      glm::vec3 lightColor = f->getMaterial()->getEmittedColor() * f->getArea();
      glm::vec3 myLightColor;
      glm::vec3 lightCentroid = f->Corner_Num(i);
      //if you want random selection use the code below instead of the above line
      
      //glm::vec3 lightCentroid = f->RandomPoint();
      glm::vec3 dirToLightCentroid = glm::normalize(lightCentroid-point); 

      // ===========================================
      // ASSIGNMENT:  ADD SHADOW & SOFT SHADOW LOGIC
      // ===========================================
      
      float distToLightCentroid = glm::length(lightCentroid-point);
      myLightColor = lightColor / float(M_PI*distToLightCentroid*distToLightCentroid);
      
      // add the lighting contribution from this particular light at this point
      // (fix this to check for blockers between the light & this surface)
      
      //hard shadow stuff here
      //cast a ray towards the light, if we hit it then light this up
      Hit tmp_hit;
      Ray lightRay(point, dirToLightCentroid);

      CastRay(lightRay,tmp_hit,false);
      if((tmp_hit.getT()<distToLightCentroid+0.001 and tmp_hit.getT()>distToLightCentroid-0.001)){
        glm::vec3 tmp=m->Shade(ray,hit,dirToLightCentroid,myLightColor,args);
        tmp/=4.0;
        answer += tmp;
      }
    }
    #endif
  }

  // ----------------------------------------------
  // add contribution from reflection, if the surface is shiny
  glm::vec3 reflectiveColor = m->getReflectiveColor();


  // =================================
  // ASSIGNMENT:  ADD REFLECTIVE LOGIC
  // =================================
  #if 1
  //first we must reflect the ray across the normal
  glm::vec3 ray_dir=ray.getDirection();
  normal=glm::normalize(normal);
  glm::vec3 reflected_dir=glm::normalize(ray_dir-2*glm::dot(ray_dir, normal)*normal);
  //next we trace this ray and get the color reflected onto it
  Ray reflected_ray(point,reflected_dir);
  Hit reflection_hit;

  if (glm::length(reflectiveColor) > 0.001 and bounce_count>0){
    answer += glm::normalize(reflectiveColor) * TraceRay(reflected_ray, reflection_hit, bounce_count-1);

  }
  //next we update answer
  #endif

  
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
