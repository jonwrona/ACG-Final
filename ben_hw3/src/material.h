#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "glCanvas.h"

#include <cassert>
#include <string>

#include "image.h"

class ArgParser;
class Ray;
class Hit;

// ====================================================================
// ====================================================================
// A simple Phong-like material 

class Material {

public:

  Material(const std::string &texture_file, const glm::vec3 &d_color,
	   const glm::vec3 &r_color, const glm::vec3 &e_color, float ind_of_refract_, bool sss_, float roughness_) {
    textureFile = texture_file;
    if (textureFile != "") {
      image = new Image(textureFile);
      ComputeAverageTextureColor();
    } else {
      diffuseColor = d_color;
      image = NULL;
    }
    reflectiveColor = r_color;
    emittedColor = e_color;
    ind_of_refract = ind_of_refract_;
    sss = sss_;
    roughness = roughness_;
    // need to initialize texture_id after glut has started
    texture_id = 0;
  }
  
  ~Material();

  // ACCESSORS
  const glm::vec3& getDiffuseColor() const { return diffuseColor; }
  const glm::vec3 getDiffuseColor(float s, float t) const;
  const glm::vec3& getReflectiveColor() const { return reflectiveColor; }
  const glm::vec3& getEmittedColor() const { return emittedColor; }
  float getRefraction() const { return ind_of_refract; }
  bool isSubsurfaceMaterial() const { return sss; }
  float getRoughness() const { return roughness; } 
  bool hasTextureMap() const { return (textureFile != ""); } 
  GLuint getTextureID();

  // SHADE
  // compute the contribution to local illumination at this point for
  // a particular light source
  glm::vec3 Shade
  (const Ray &ray, const Hit &hit, const glm::vec3 &dirToLight, 
   const glm::vec3 &lightColor, ArgParser *args) const;
  
protected:

  Material() { exit(0); }
  Material(const Material&) { exit(0); }
  const Material& operator=(const Material&) { exit(0); }

  void ComputeAverageTextureColor();

  // REPRESENTATION
  glm::vec3 diffuseColor;
  glm::vec3 reflectiveColor;
  glm::vec3 emittedColor;
  float ind_of_refract;
  bool sss;
  float roughness;

  std::string textureFile;
  GLuint texture_id;
  Image *image;
};

// ====================================================================
// ====================================================================

#endif

