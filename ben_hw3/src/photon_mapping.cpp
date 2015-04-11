#include "glCanvas.h"

#include <iostream>
#include <algorithm>

#include "argparser.h"
#include "photon_mapping.h"
#include "mesh.h"
#include "face.h"
#include "primitive.h"
#include "kdtree.h"
#include "utils.h"
#include "raytracer.h"
#include <glm/gtx/random.hpp>


// ==========
// DESTRUCTOR
PhotonMapping::~PhotonMapping() {
  // cleanup all the photons
  delete kdtree;
}


// ========================================================================
// Recursively trace a single photon

void PhotonMapping::TracePhoton(const glm::vec3 &position, const glm::vec3 &direction, 
				const glm::vec3 &energy, int iter) {




  // ==============================================
  // ASSIGNMENT: IMPLEMENT RECURSIVE PHOTON TRACING
  // ==============================================

  // Trace the photon through the scene.  At each diffuse or
  // reflective bounce, store the photon in the kd tree.

  // One optimization is to *not* store the first bounce, since that
  // direct light can be efficiently computed using classic ray
  // tracing.

  //first check if this photon has enough energy to worry about
  if (glm::length(energy) < 0.0000025 or iter > args->num_bounces){
    return;
  }
  else{
    //set up what we need to cast a ray
    Hit hit;
    Ray ray(position, glm::normalize(direction));
    if (raytracer->CastRay(ray, hit, false)){
      //std::cout<<"we hit something!!!\n";
      //we hit something!
      Material* m=hit.getMaterial();
      glm::vec3 pos=position+hit.getT()*direction;
      glm::vec3 normal=hit.getNormal();
      glm::vec3 ones(1.0,1.0,1.0);

      //we will need these next variables to store the photon hit
      glm::vec3 photon_energy;
      glm::vec3 bounce_dir;
      glm::vec3 bounce_energy;
      //check for diffuse/reflective surface and compute the photon and bounce dir
      if(glm::length(m->getReflectiveColor()) > 0.1){
        //we are on a reflective object
        photon_energy = (ones - m->getReflectiveColor())*energy;
        //std::cout<<m->getReflectiveColor().x<<' '<<m->getReflectiveColor().y<<' '<<m->getReflectiveColor().z<<std::endl;
        bounce_energy = m->getReflectiveColor()*energy;
        //compute the bounce direction
        bounce_dir=glm::normalize(direction-2*glm::dot(direction, normal)*normal);
      }
      else{
        //we are on a diffuse object
        photon_energy = (ones - m->getDiffuseColor())*energy;
        //std::cout<<m->getDiffuseColor().x<<' '<<m->getDiffuseColor().y<<' '<<m->getDiffuseColor().z<<std::endl;
        bounce_energy = m->getDiffuseColor()*energy;
        //compute the random bounce direction
        bounce_dir=RandomDiffuseDirection(normal);
      }
      
      //add the photon to the kd-tree and trace it 
      if (iter>0){
        Photon photon(pos, -direction, photon_energy, iter);
        kdtree->AddPhoton(photon);
      }
      TracePhoton(pos, bounce_dir, bounce_energy, iter+1);
    }

  }



}


// ========================================================================
// Trace the specified number of photons through the scene

void PhotonMapping::TracePhotons() {
  std::cout << "trace photons" << std::endl;

  // first, throw away any existing photons
  delete kdtree;

  // consruct a kdtree to store the photons
  BoundingBox *bb = mesh->getBoundingBox();
  glm::vec3 min = bb->getMin();
  glm::vec3 max = bb->getMax();
  glm::vec3 diff = max-min;
  min -= 0.001f*diff;
  max += 0.001f*diff;
  kdtree = new KDTree(BoundingBox(min,max));

  // photons emanate from the light sources
  const std::vector<Face*>& lights = mesh->getLights();

  // compute the total area of the lights
  float total_lights_area = 0;
  for (unsigned int i = 0; i < lights.size(); i++) {
    total_lights_area += lights[i]->getArea();
  }

  // shoot a constant number of photons per unit area of light source
  // (alternatively, this could be based on the total energy of each light)
  for (unsigned int i = 0; i < lights.size(); i++) {  
    float my_area = lights[i]->getArea();
    int num = args->num_photons_to_shoot * my_area / total_lights_area;
    // the initial energy for this photon
    glm::vec3 energy = my_area/float(num) * lights[i]->getMaterial()->getEmittedColor();
    glm::vec3 normal = lights[i]->computeNormal();
    for (int j = 0; j < num; j++) {
      glm::vec3 start = lights[i]->RandomPoint();
      // the initial direction for this photon (for diffuse light sources)
      glm::vec3 direction = RandomDiffuseDirection(normal);
      TracePhoton(start,direction,energy,0);
    }
  }
}


// ======================================================================

// helper function
bool closest_photon(const std::pair<Photon,float> &a, const std::pair<Photon,float> &b) {
  return (a.second < b.second);
}


// ======================================================================
glm::vec3 PhotonMapping::GatherIndirect(const glm::vec3 &point, const glm::vec3 &normal, const glm::vec3 &direction_from) const {


  if (kdtree == NULL) { 
    std::cout << "WARNING: Photons have not been traced throughout the scene." << std::endl;
    return glm::vec3(0,0,0); 
  }
  glm::vec3 energy(0,0,0);


  // ================================================================
  // ASSIGNMENT: GATHER THE INDIRECT ILLUMINATION FROM THE PHOTON MAP
  // ================================================================

  // collect the closest args->num_photons_to_collect photons
  // determine the radius that was necessary to collect that many photons
  // average the energy of those photons over that radius

  //first setup a bounding box of a certain radius
  float radius = .2;
  int num_photons=args->num_photons_to_collect;

  while(true){
    if (radius>100.0){
      break;
    }
    glm::vec3 bound_min(-radius, -radius, -radius);
    glm::vec3 bound_max(radius, radius, radius);
    bound_min+=point;
    bound_max+=point;
    BoundingBox box(bound_min, bound_max);

    //get all photons potentially in the bounding box
    std::vector<Photon> photons;
    kdtree->CollectPhotonsInBox(box, photons);
    //std::cout<<"we caught "<<photons.size()<<" photons in this box\n";
    
    //trim that down to the photons actually within that radius
    std::vector<std::pair<Photon, float> > p_dists;
    for (unsigned int i=0; i<photons.size(); i++){
      p_dists.push_back(std::make_pair(photons[i], glm::length(photons[i].getPosition()-point)));
    }
    std::sort(p_dists.begin(), p_dists.end(), closest_photon);

    int count_photons=0;
    if (p_dists.size()< num_photons){
      //std::cout<<radius<<" was not big enough\n";
      radius*=2.0;

      continue;
    }
    else if (p_dists[num_photons-1].second>radius){
      //std::cout<<radius<<" was not big enough\n";
      radius*=2.0;
      continue;
    }
    else{
      //in this case we have found enough photons
      radius=p_dists[num_photons-1].second;
      //std::cout<<radius<<std::endl;
      for (unsigned int i=0; i<num_photons; i++){
        energy+=p_dists[i].first.getEnergy();
      }
      energy /= .3*3.14159265359*radius*radius;
      break;
    }
  }
  //while there aren't enough photons double the radius and repeat
  

  // return the color
  return energy;


}


// ======================================================================
// PHOTON VISUALIZATION FOR DEBUGGING
// ======================================================================

void PhotonMapping::initializeVBOs() {
  HandleGLError("enter photonmapping initializevbos()");
  glGenBuffers(1, &photon_direction_verts_VBO);
  glGenBuffers(1, &photon_direction_indices_VBO);
  glGenBuffers(1, &kdtree_verts_VBO);
  glGenBuffers(1, &kdtree_edge_indices_VBO);
  HandleGLError("leave photonmapping initializevbos()");
}

void PhotonMapping::setupVBOs() {
  HandleGLError("enter photonmapping setupvbos()");

  photon_direction_verts.clear();
  photon_direction_indices.clear();
  kdtree_verts.clear();
  kdtree_edge_indices.clear();

  // initialize the data
  BoundingBox *bb = mesh->getBoundingBox();
  float max_dim = bb->maxDim();

  if (kdtree == NULL) return;
  std::vector<const KDTree*> todo;  
  todo.push_back(kdtree);
  while (!todo.empty()) {
    const KDTree *node = todo.back();
    todo.pop_back(); 
    if (node->isLeaf()) {

      // initialize photon direction vbo
      const std::vector<Photon> &photons = node->getPhotons();
      int num_photons = photons.size();
      for (int i = 0; i < num_photons; i++) {
	const Photon &p = photons[i];
	glm::vec3 energy = p.getEnergy()*float(args->num_photons_to_shoot);
        glm::vec4 color(energy.x,energy.y,energy.z,1);
	const glm::vec3 &position = p.getPosition();
	glm::vec3 other = position - p.getDirectionFrom()*0.02f*max_dim;
        addEdgeGeometry(photon_direction_verts,photon_direction_indices,
                        position,other,color,color,max_dim*0.0005f,0);
      }

      // initialize kdtree vbo
      float thickness = 0.001*max_dim;
      glm::vec3 A = node->getMin();
      glm::vec3 B = node->getMax();
      glm::vec4 black(1,0,0,1);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(A.x,A.y,A.z),glm::vec3(A.x,A.y,B.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(A.x,A.y,B.z),glm::vec3(A.x,B.y,B.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(A.x,B.y,B.z),glm::vec3(A.x,B.y,A.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(A.x,B.y,A.z),glm::vec3(A.x,A.y,A.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(B.x,A.y,A.z),glm::vec3(B.x,A.y,B.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(B.x,A.y,B.z),glm::vec3(B.x,B.y,B.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(B.x,B.y,B.z),glm::vec3(B.x,B.y,A.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(B.x,B.y,A.z),glm::vec3(B.x,A.y,A.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(A.x,A.y,A.z),glm::vec3(B.x,A.y,A.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(A.x,A.y,B.z),glm::vec3(B.x,A.y,B.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(A.x,B.y,B.z),glm::vec3(B.x,B.y,B.z),black,black,thickness,thickness);
      addEdgeGeometry(kdtree_verts,kdtree_edge_indices,glm::vec3(A.x,B.y,A.z),glm::vec3(B.x,B.y,A.z),black,black,thickness,thickness);

    } else {
      todo.push_back(node->getChild1());
      todo.push_back(node->getChild2());
    } 
  }



  // copy the data to each VBO
  if (photon_direction_verts.size() > 0) {
    glBindBuffer(GL_ARRAY_BUFFER,photon_direction_verts_VBO); 
    glBufferData(GL_ARRAY_BUFFER,
                 sizeof(VBOPosNormalColor) * photon_direction_verts.size(),
                 &photon_direction_verts[0],
                 GL_STATIC_DRAW); 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,photon_direction_indices_VBO); 
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 sizeof(VBOIndexedTri) * photon_direction_indices.size(),
                 &photon_direction_indices[0], GL_STATIC_DRAW);
  }
  if (kdtree_verts.size() > 0) {
    glBindBuffer(GL_ARRAY_BUFFER,kdtree_verts_VBO); 
    glBufferData(GL_ARRAY_BUFFER,
                 sizeof(VBOPosNormalColor) * kdtree_verts.size(),
                 &kdtree_verts[0],
                 GL_STATIC_DRAW); 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,kdtree_edge_indices_VBO); 
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 sizeof(VBOIndexedTri) * kdtree_edge_indices.size(),
                 &kdtree_edge_indices[0], GL_STATIC_DRAW);
  }

  HandleGLError("leave photonmapping setupvbos()");
}

void PhotonMapping::drawVBOs() {
  HandleGLError("enter photonmapping drawvbos()");

  glUniform1i(GLCanvas::colormodeID, 1);
  if (args->render_photons && photon_direction_verts.size() > 0) {
    glBindBuffer(GL_ARRAY_BUFFER,photon_direction_verts_VBO); 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,photon_direction_indices_VBO); 
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2 + sizeof(glm::vec4)));
    glDrawElements(GL_TRIANGLES,
                   photon_direction_indices.size()*3,
                   GL_UNSIGNED_INT, 0);
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
  }

  if (args->render_kdtree && kdtree_edge_indices.size() > 0) {
    glBindBuffer(GL_ARRAY_BUFFER,kdtree_verts_VBO); 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,kdtree_edge_indices_VBO); 
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2 + sizeof(glm::vec4)));
    glDrawElements(GL_TRIANGLES,
                   kdtree_edge_indices.size()*3,
                   GL_UNSIGNED_INT, 0);
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
  }

  HandleGLError("leave photonmapping drawvbos()");
}

void PhotonMapping::cleanupVBOs() {
  glDeleteBuffers(1, &photon_direction_verts_VBO);
  glDeleteBuffers(1, &photon_direction_indices_VBO);
  glDeleteBuffers(1, &kdtree_verts_VBO);
  glDeleteBuffers(1, &kdtree_edge_indices_VBO);
}

