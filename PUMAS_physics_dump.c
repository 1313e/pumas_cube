// Standard libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// PUMAS
#include "pumas.h"

// Main function for dumping physics to a binary file
int main(){
    // Obtain name of XML file
    const char *XML_file = "./materials/materials.xml";
    const char *dump_file = "./materials/materials_dump";

    // Declare physics object
    struct pumas_physics *physics;

    // Create physics
    pumas_physics_create(&physics, PUMAS_PARTICLE_MUON, XML_file,
                         "materials", NULL);

    // Dump to binary file
    FILE *file = fopen(dump_file, "wb+");
    pumas_physics_dump(physics, file);
    fclose(file);

    // Release memory
    pumas_physics_destroy(&physics);
}