#include "glCanvas.h"

#include "radiosity.h"
#include "mesh.h"
#include "face.h"
#include "glCanvas.h"
#include "sphere.h"
#include "raytree.h"
#include "raytracer.h"
#include "utils.h"

// ================================================================
// CONSTRUCTOR & DESTRUCTOR
// ================================================================
Radiosity::Radiosity(Mesh *m, ArgParser *a) {
  mesh = m;
  args = a;
  num_faces = -1;  
  formfactors = NULL;
  area = NULL;
  undistributed = NULL;
  absorbed = NULL;
  radiance = NULL;
  max_undistributed_patch = -1;
  total_area = -1;
  Reset();
}

Radiosity::~Radiosity() {
  Cleanup();
  cleanupVBOs();
}

void Radiosity::Cleanup() {
  delete [] formfactors;
  delete [] area;
  delete [] undistributed;
  delete [] absorbed;
  delete [] radiance;
  num_faces = -1;
  formfactors = NULL;
  area = NULL;
  undistributed = NULL;
  absorbed = NULL;
  radiance = NULL;
  max_undistributed_patch = -1;
  total_area = -1;
}

void Radiosity::Reset() {
  delete [] area;
  delete [] undistributed;
  delete [] absorbed;
  delete [] radiance;

  // create and fill the data structures
  num_faces = mesh->numFaces();
  area = new float[num_faces];
  undistributed = new glm::vec3[num_faces];
  absorbed = new glm::vec3[num_faces];
  radiance = new glm::vec3[num_faces];
  for (int i = 0; i < num_faces; i++) {
    Face *f = mesh->getFace(i);
    f->setRadiosityPatchIndex(i);
    setArea(i,f->getArea());
    glm::vec3 emit = f->getMaterial()->getEmittedColor();
    setUndistributed(i,emit);
    setAbsorbed(i,glm::vec3(0,0,0));
    setRadiance(i,emit);
  }

  // find the patch with the most undistributed energy
  findMaxUndistributed();
}


// =======================================================================================
// =======================================================================================

void Radiosity::findMaxUndistributed() {
  // find the patch with the most undistributed energy 
  // don't forget that the patches may have different sizes!
  max_undistributed_patch = -1;
  total_undistributed = 0;
  total_area = 0;
  float max = -1;
  for (int i = 0; i < num_faces; i++) {
    float m = glm::length(getUndistributed(i)) * getArea(i);
    total_undistributed += m;
    total_area += getArea(i);
    if (max < m) {
      max = m;
      max_undistributed_patch = i;
    }
  }
  assert (max_undistributed_patch >= 0 && max_undistributed_patch < num_faces);
}


void Radiosity::ComputeFormFactors() {
  assert (formfactors == NULL);
  assert (num_faces > 0);
  formfactors = new float[num_faces*num_faces];


  // =====================================
  // ASSIGNMENT:  COMPUTE THE FORM FACTORS
  // =====================================

  //cycle through every pair
  for (int i=0; i<num_faces; i++){
    for (int j=0; j<num_faces; j++){
      //a surface does not distribute any light to itself
      float form_factor=0.0;
      if (i==j){
        formfactors[i*num_faces+j]=0.0;
      }
      //check for occlusions...

      
      //here we actually have to calculate the form factor...
      else{
        for (int samples=0; samples<args->num_form_factor_samples; samples++){
          float phi_i, phi_j, Aj, r;
          glm::vec3 disp;
          
          Face* face_i=mesh->getFace(i);
          Face* face_j=mesh->getFace(j);
          glm::vec3 norm_i= glm::normalize(face_i->computeNormal());
          glm::vec3 norm_j= glm::normalize(face_j->computeNormal());
          //first compute the displacement vector between the two faces.
          //points from i to j

          //use below code when you want to trace from center to center

          //disp = face_i->computeCentroid()-face_j->computeCentroid();
          glm::vec3 face_j_point=face_j->RandomPoint();
          disp = face_i->RandomPoint()-face_j_point;
          Hit h;
          r = glm::length(disp);
          
          //if we are only taking one sample we want it from the center for consistency
          if (args->num_form_factor_samples==1){
            disp=face_i->computeCentroid()-face_j->computeCentroid();
            r = glm::length(disp);
            face_j_point=face_j->computeCentroid();
          }
          Ray ray(face_j_point, glm::normalize(disp));
          raytracer->CastRay(ray, h, true);
          //here we test for occlusion
          if ( args->num_shadow_samples>0 and !(h.getT()>r-0.001 and h.getT()<r+0.001) ){
            //setFormFactor(i,j, 0.0);
          }
          else{
          
            //std::cout<<"r is "<<r<<std::endl;
            disp = -glm::normalize(disp);
            //compute the rest of the constants
            phi_i = fabs(glm::dot(norm_i, disp));
            //std::cout<<"phi_i is "<<phi_i<<std::endl;
            disp = -disp;
            phi_j = fabs(glm::dot(norm_j, disp));
            //std::cout<<"phi_j is "<<phi_j<<std::endl;
            Aj = face_j->getArea();
            //std::cout<<"Aj is "<<Aj<<std::endl;
            //compute and enter the form factor
            //setFormFactor(i,j,(phi_i)*(phi_j)*Aj/(3.14159265359*r*r));
            form_factor+=(phi_i)*(phi_j)*Aj/(3.14159265359*r*r);
            // formfactors[i*num_faces+j]=10000000000000000.0*glm::cos(phi_i)*glm::cos(phi_j)*Aj/(3.14159265359*r*r);
          }
        }
      }
    setFormFactor(i,j,form_factor/(float(args->num_form_factor_samples)));

    }
    normalizeFormFactors(i);
  }
  #if 0
  std::cout<<"formfactors matrix is\n";
  for (int i=0; i<num_faces; i++){
    for(int j=0; j<num_faces;j++){
      std::cout<<formfactors[i*num_faces+j]<<' ';
    }
    std::cout<<'\n';
  }
  #endif

}





// ================================================================
// ================================================================

float Radiosity::Iterate() {
  if (formfactors == NULL) 
    ComputeFormFactors();
  assert (formfactors != NULL);




  // ==========================================
  // ASSIGNMENT:  IMPLEMENT RADIOSITY ALGORITHM
  // ==========================================

  //Compute the radiosity matrix. I feel like I shouldn't have to do this 
  //every time. Maybe this is what the form factors matrix was meant for

  //save the radiances
  std::vector<glm::vec3> radiations, absorbeds;
  for (int i=0; i<num_faces; i++){
    
    Face* face=mesh->getFace(i);
    glm::vec3 ro_i=face->getMaterial()->getDiffuseColor();
    glm::vec3 ones(1,1,1);
    glm::vec3 abso_i=ones-ro_i;

    glm::vec3 radiance=face->getMaterial()->getEmittedColor();
    glm::vec3 absorbed_i(0,0,0);//=getAbsorbed(i);
    for (int j=0; j<num_faces; j++){
      //new_row.push_back(-ro_i*getFormFactor(i,j))
      //might need to reverse that j,i
      radiance+=ro_i*getFormFactor(j,i)*getRadiance(j); 
      absorbed_i+=abso_i*getFormFactor(j,i)*getRadiance(j);

    } 
    //the undistributed light is the stuff we got on this iteration
    setUndistributed(i, radiance-getRadiance(i));  
    //std::cout<<glm::length(radiance)<<' ';
    
    radiations.push_back(radiance);
    absorbeds.push_back(absorbed_i);
  }

  //std::cout<<std::endl;
  //set the computed radiances
  for(int i=0; i<num_faces; i++){
    setRadiance(i, radiations[i]);
    setAbsorbed(i, absorbeds[i]);
  }
  //now we must calculate the undistributed light
  float total_undistributed = 0;
  glm::vec3 ambient;
  for (int i=0; i<num_faces;i++){
    total_undistributed+=glm::length(getUndistributed(i));
    ambient+=getUndistributed(i);
  }
  
  //uncomment this for the ambient lighting term. This implementation, while simple
  //breaks my visualization for undistributed
  #if 1
  ambient/=float(num_faces);
  for (int i=0; i<num_faces; i++){
    setRadiance(i, getRadiance(i)+ambient);
  }
  #endif

  // return the total light yet undistributed
  // (so we can decide when the solution has sufficiently converged)
  std::cout<<total_undistributed<<std::endl;
  return total_undistributed;




}


// =======================================================================================
// VBO & DISPLAY FUNCTIONS
// =======================================================================================

// for interpolation
void CollectFacesWithVertex(Vertex *have, Face *f, std::vector<Face*> &faces) {
  for (unsigned int i = 0; i < faces.size(); i++) {
    if (faces[i] == f) return;
  }
  if (have != (*f)[0] && have != (*f)[1] && have != (*f)[2] && have != (*f)[3]) return;
  faces.push_back(f);
  for (int i = 0; i < 4; i++) {
    Edge *ea = f->getEdge()->getOpposite();
    Edge *eb = f->getEdge()->getNext()->getOpposite();
    Edge *ec = f->getEdge()->getNext()->getNext()->getOpposite();
    Edge *ed = f->getEdge()->getNext()->getNext()->getNext()->getOpposite();
    if (ea != NULL) CollectFacesWithVertex(have,ea->getFace(),faces);
    if (eb != NULL) CollectFacesWithVertex(have,eb->getFace(),faces);
    if (ec != NULL) CollectFacesWithVertex(have,ec->getFace(),faces);
    if (ed != NULL) CollectFacesWithVertex(have,ed->getFace(),faces);
  }
}

// different visualization modes
glm::vec3 Radiosity::setupHelperForColor(Face *f, int i, int j) {
  assert (mesh->getFace(i) == f);
  assert (j >= 0 && j < 4);
  if (args->render_mode == RENDER_MATERIALS) {
    return f->getMaterial()->getDiffuseColor();
  } else if (args->render_mode == RENDER_RADIANCE && args->interpolate == true) {
    std::vector<Face*> faces;
    CollectFacesWithVertex((*f)[j],f,faces);
    float total = 0;
    glm::vec3 color = glm::vec3(0,0,0);
    glm::vec3 normal = f->computeNormal();
    for (unsigned int i = 0; i < faces.size(); i++) {
      glm::vec3 normal2 = faces[i]->computeNormal();
      float area = faces[i]->getArea();
      if (glm::dot(normal,normal2) < 0.5) continue;
      assert (area > 0);
      total += area;
      color += float(area) * getRadiance(faces[i]->getRadiosityPatchIndex());
    }
    assert (total > 0);
    color /= total;
    return color;
  } else if (args->render_mode == RENDER_LIGHTS) {
    return f->getMaterial()->getEmittedColor();
  } else if (args->render_mode == RENDER_UNDISTRIBUTED) { 
    return getUndistributed(i);
  } else if (args->render_mode == RENDER_ABSORBED) {
    return getAbsorbed(i);
  } else if (args->render_mode == RENDER_RADIANCE) {
    return getRadiance(i);
  } else if (args->render_mode == RENDER_FORM_FACTORS) {
    if (formfactors == NULL) ComputeFormFactors();
    float scale = 0.2 * total_area/getArea(i);
    float factor = scale * getFormFactor(max_undistributed_patch,i);
    return glm::vec3(factor,factor,factor);
  } else {
    assert(0);
  }
  exit(0);
}



void Radiosity::initializeVBOs() {
  // create a pointer for the vertex & index VBOs
  glGenBuffers(1, &mesh_tri_verts_VBO);
  glGenBuffers(1, &mesh_tri_indices_VBO);
  glGenBuffers(1, &mesh_textured_tri_indices_VBO);
}


void Radiosity::setupVBOs() {
  HandleGLError("enter radiosity setupVBOs()");
  mesh_tri_verts.clear();
  mesh_tri_indices.clear();
  mesh_textured_tri_indices.clear();

  // initialize the data in each vector
  int num_faces = mesh->numFaces();
  assert (num_faces > 0);
  for (int i = 0; i < num_faces; i++) {
    Face *f = mesh->getFace(i);
    Edge *e = f->getEdge();
    glm::vec3 normal = f->computeNormal();

    double avg_s = 0;
    double avg_t = 0;
    glm::vec3 avg_color(0,0,0);

    int start = mesh_tri_verts.size();

    // wireframe is normally black, except when it's the special
    // patch, then the wireframe is red
    glm::vec4 wireframe_color(0,0,0,0.5);
    if (args->render_mode == RENDER_FORM_FACTORS && i == max_undistributed_patch) {
      wireframe_color = glm::vec4(1,0,0,1);
    }

    // add the 4 corner vertices
    for (int j = 0; j < 4; j++) {
      glm::vec3 pos = ((*f)[j])->get();
      double s = (*f)[j]->get_s();
      double t = (*f)[j]->get_t();
      glm::vec3 color = setupHelperForColor(f,i,j);
      color = glm::vec3(linear_to_srgb(color.r),
                        linear_to_srgb(color.g),
                        linear_to_srgb(color.b));
      avg_color += 0.25f * color;
      mesh_tri_verts.push_back(VBOPosNormalColor(pos,normal,
                                                 glm::vec4(color.r,color.g,color.b,1.0),
                                                 wireframe_color,
                                                 s,t));
      avg_s += 0.25 * s;
      avg_t += 0.25 * t;
      e = e->getNext();
    }

    // the centroid (for wireframe rendering)
    glm::vec3 centroid = f->computeCentroid();
    mesh_tri_verts.push_back(VBOPosNormalColor(centroid,normal,
                                               glm::vec4(avg_color.r,avg_color.g,avg_color.b,1),
                                               glm::vec4(1,1,1,1),
                                               avg_s,avg_t));

    if (f->getMaterial()->hasTextureMap()) {
      mesh_textured_tri_indices.push_back(VBOIndexedTri(start+0,start+1,start+4));
      mesh_textured_tri_indices.push_back(VBOIndexedTri(start+1,start+2,start+4));
      mesh_textured_tri_indices.push_back(VBOIndexedTri(start+2,start+3,start+4));
      mesh_textured_tri_indices.push_back(VBOIndexedTri(start+3,start+0,start+4));
    } else {
      mesh_tri_indices.push_back(VBOIndexedTri(start+0,start+1,start+4));
      mesh_tri_indices.push_back(VBOIndexedTri(start+1,start+2,start+4));
      mesh_tri_indices.push_back(VBOIndexedTri(start+2,start+3,start+4));
      mesh_tri_indices.push_back(VBOIndexedTri(start+3,start+0,start+4));
    }
  }
  assert ((int)mesh_tri_verts.size() == num_faces*5);
  assert ((int)mesh_tri_indices.size() + (int)mesh_textured_tri_indices.size() == num_faces*4);
  
  // copy the data to each VBO
  glBindBuffer(GL_ARRAY_BUFFER,mesh_tri_verts_VBO); 
  glBufferData(GL_ARRAY_BUFFER,
	       sizeof(VBOPosNormalColor) * num_faces * 5,
	       &mesh_tri_verts[0],
	       GL_STATIC_DRAW); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_tri_indices_VBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
	       sizeof(VBOIndexedTri) * mesh_tri_indices.size(),
	       &mesh_tri_indices[0], GL_STATIC_DRAW);
  if (mesh_textured_tri_indices.size() > 0) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_textured_tri_indices_VBO); 
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 sizeof(VBOIndexedTri) * mesh_textured_tri_indices.size(),
                 &mesh_textured_tri_indices[0], GL_STATIC_DRAW);
    
  }

  HandleGLError("radiosity setupVBOs() just before texture");

  // WARNING: this naive VBO implementation only allows a single texture
  // FIXME: something still buggy about textures
  int num_textured_materials = 0;
  for (unsigned int mat = 0; mat < mesh->materials.size(); mat++) {
    Material *m = mesh->materials[mat];
    if (m->hasTextureMap()) {
      // FIXME: old gl...
      //glBindTexture(GL_TEXTURE_2D,m->getTextureID());
      num_textured_materials++;
    }
  }
  assert (num_textured_materials <= 1);

  HandleGLError("leave radiosity setupVBOs()");
}


void Radiosity::drawVBOs() {

  // =====================
  // DRAW ALL THE POLYGONS

  assert ((int)mesh_tri_indices.size() + (int)mesh_textured_tri_indices.size() == num_faces*4);

  // render with Phong lighting?
  if (args->render_mode == RENDER_MATERIALS) {
    // yes
    glUniform1i(GLCanvas::colormodeID, 1);
  } else {
    // no
    glUniform1i(GLCanvas::colormodeID, 0);
  }

  // render untextured faces
  glBindBuffer(GL_ARRAY_BUFFER,mesh_tri_verts_VBO); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_tri_indices_VBO); 
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)0);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));
  glEnableVertexAttribArray(3);
  glVertexAttribPointer(3, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2 + sizeof(glm::vec4)));
  glEnableVertexAttribArray(4);
  glVertexAttribPointer(4, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2 + sizeof(glm::vec4)*2));
  glDrawElements(GL_TRIANGLES, mesh_tri_indices.size()*3,GL_UNSIGNED_INT, 0);
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);
  glDisableVertexAttribArray(3);
  glDisableVertexAttribArray(4);


  // render faces with textures
  if (mesh_textured_tri_indices.size() > 0) {

    // FIXME: there is something buggy with textures still
    //glUniform1i(GLCanvas::colormodeID, 2);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, GLCanvas::textureID);
    GLCanvas::mytexture = glGetUniformLocation(GLCanvas::programID, "mytexture");
    glUniform1i(GLCanvas::mytexture, /*GL_TEXTURE*/0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_textured_tri_indices_VBO); 
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2 + sizeof(glm::vec4)));
    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2 + sizeof(glm::vec4)*2));
    glDrawElements(GL_TRIANGLES, mesh_textured_tri_indices.size()*3, GL_UNSIGNED_INT, 0);
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
    glDisableVertexAttribArray(4);

    //glUniform1i(GLCanvas::colormodeID, 1);
  }

  HandleGLError(); 
}


void Radiosity::cleanupVBOs() {
  glDeleteBuffers(1, &mesh_tri_verts_VBO);
  glDeleteBuffers(1, &mesh_tri_indices_VBO);
  glDeleteBuffers(1, &mesh_textured_tri_indices_VBO);

  glDeleteTextures(1, &GLCanvas::textureID);
}

