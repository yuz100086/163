# version 330 core
// Do not use any version older than 330!

/* This is the fragment shader for reading in a scene description, including 
   lighting.  Uniform lights are specified from the main program, and used in 
   the shader.  As well as the material parameters of the object.  */

// Inputs to the fragment shader are the outputs of the same name of the vertex shader.
// Note that the default output, gl_Position, is inaccessible!
in vec3 mynormal; 
in vec4 myvertex; 

// You will certainly need this matrix for your lighting calculations
uniform mat4 modelview;

// This first defined output of type vec4 will be the fragment color
out vec4 fragColor;

uniform vec3 color;

const int numLights = 10; 
uniform bool enablelighting; // are we lighting at all (global).
uniform vec4 lightposn[numLights]; // positions of lights 
uniform vec4 lightcolor[numLights]; // colors of lights
uniform int numused;               // number of lights used

// Now, set the material parameters.
// I use ambient, diffuse, specular, shininess. 
// But, the ambient is just additive and doesn't multiply the lights.  

uniform vec4 ambient; 
uniform vec4 diffuse; 
uniform vec4 specular; 
uniform vec4 emission; 
uniform float shininess; 

// Copied from hw1
vec4 ComputeLight (vec3 direction, vec4 lightcolor, vec3 normal, vec3 halfvec, 
vec4 mydiffuse, vec4 myspecular, float myshininess) {

        float nDotL = dot(normal, direction)  ;         
        vec4 lambert = mydiffuse * lightcolor * max (nDotL, 0.0) ;  

        float nDotH = dot(normal, halfvec) ; 
        vec4 phong = myspecular * lightcolor * pow (max(nDotH, 0.0), myshininess) ; 

        vec4 retval = lambert + phong ; 
        return retval ;
}       


void main (void) 
{       
    if (enablelighting) {       
		//ambient is additive
        vec4 finalcolor = ambient; 
		
		const vec3 eyepos = vec3(0,0,0);
		vec3 mypos = (modelview * myvertex).xyz / (modelview * myvertex).w;
		vec3 eyedirn = normalize(eyepos - mypos);
		vec3 normal = normalize((transpose(inverse(modelview)) * vec4(mynormal,0)).xyz);

        // YOUR CODE FOR HW 2 HERE
        // A key part is implementation of the fragment shader
		
		vec3 direction;
		for (int i = 0; i < numused; i++) {
			// Directional light
			if(lightposn[i].w == 0){
				direction = normalize(lightposn[i].xyz);
			}else{
			// Point light
				vec3 position = lightposn[i].xyz / lightposn[i].w;
				direction = normalize(position - mypos);	
			}
			vec3 half = normalize(direction + eyedirn);
			vec4 col = ComputeLight (direction, lightcolor[i], normal, half, diffuse, specular, shininess);
			finalcolor = finalcolor + col;
		}

        fragColor = finalcolor; 
    } else {
        fragColor = vec4(color, 1.0f); 
    }
}
